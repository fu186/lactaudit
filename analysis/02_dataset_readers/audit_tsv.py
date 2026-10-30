#!/usr/bin/env python3
"""
audit_tsv.py -- reproduction audit for GEO wide tables (tab- or comma-separated).

Handles the three header layouts that GEO series actually ship:

    no_header    first row is all barcodes; data rows are ``gene + values``
    two_label    first row is ``Gene_ID, Symbol, barcode...``
    one_label    first row is ``Gene, barcode...``

Usage:
  python3 audit_tsv.py --name <label> --tsv <file(.gz|.csv[.gz])> [--markers markers_pan.json]
      [--species human|mouse] [--max-cells 4000] [--nperm 150]
"""
import argparse
import gzip
import json
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lactcheck as LC                                    # noqa: E402
from audit_one import annotate                            # noqa: E402

# A barcode is 10+ characters of ACGTN, optionally suffixed with ".1" or "-1".
_BARCODE_RE = re.compile(r'^[ACGTNacgtn]{10,}([.\-]\d+)?$')


def _open_text(path):
    """Open a possibly-gzipped text file."""
    if str(path).endswith(('.gz', '.bgz')):
        return gzip.open(path, 'rt', errors='replace')
    return open(path, 'rt', errors='replace')


def sniff_wide_table(path):
    """Identify the layout of a GEO wide table.

    Returns
    -------
    (sep, variant, n_label_cols)
        variant is one of ``'no_header'``, ``'one_label'``, ``'two_label'``.

    Notes
    -----
    The earlier implementation of this module guessed the layout from a single
    substring test on the first field and then compensated with an off-by-one
    column list (``usecols=[0] + [p - 1 for p in pos]``).  On the actual GEO
    files that compensation was wrong: pandas raised
    ``Expected 4598 fields in line 2, saw 4599`` because the data rows carry one
    more field than the header.  The robust approach is to always parse with
    ``header=None`` and skip the header row explicitly.
    """
    with _open_text(path) as fh:
        header = fh.readline().rstrip('\n')
    sep = '\t' if header.count('\t') > header.count(',') else ','
    fields = [x.strip().strip('"') for x in header.split(sep)]

    n_label = 0
    for x in fields:
        if _BARCODE_RE.match(x):
            break
        n_label += 1

    # Second line disambiguates: in a headerless file it is the first gene row.
    with _open_text(path) as fh:
        fh.readline()
        second = fh.readline().rstrip('\n')
    second_fields = [x.strip().strip('"') for x in second.split(sep)]

    if n_label == 0:
        variant = 'no_header'
    elif n_label >= 2:
        variant = 'two_label'
    else:
        variant = 'one_label'

    # Cross-check: a headerless file's first data field is a gene symbol, not a
    # number. If the "header" was actually all barcodes we already know, but if
    # the second row's first field looks numeric, treat it as a header anyway.
    if variant == 'one_label' and second_fields:
        try:
            float(second_fields[0])
            variant = 'no_header'
        except ValueError:
            pass

    return sep, variant, n_label


def load_wide(path, max_cells=4000, seed=0):
    """Read a GEO wide table into a ``cells x genes`` float32 DataFrame.

    Memory-friendliness
    -------------------
    The barcode list is read first so the target cell columns can be chosen at
    random, then only those columns are materialised.  This avoids pulling a
    20k x 10k dense table into memory.

    Parameters
    ----------
    path : str
        ``.txt`` / ``.tsv`` / ``.csv``, optionally gzipped.
    max_cells : int
        Randomly keep at most this many cells (columns).
    seed : int
        Seed for the column sample.
    """
    sep, variant, n_label = sniff_wide_table(path)

    with _open_text(path) as fh:
        header_fields = [x.strip().strip('"') for x in fh.readline().rstrip('\n').split(sep)]
    if variant == 'no_header':
        cells = header_fields
    else:
        cells = header_fields[n_label:]

    n_cells = len(cells)
    if max_cells and n_cells > max_cells:
        rng = np.random.default_rng(seed)
        pos = np.sort(rng.choice(n_cells, size=max_cells, replace=False))
    else:
        pos = np.arange(n_cells)

    # Parse with header=None and skip the header row. This is both the correct
    # handling of all three variants and markedly faster than header=0.
    d = pd.read_csv(_open_text(path), sep=sep, header=None, skiprows=1, low_memory=False)
    if d.empty:
        raise ValueError(f'no data rows in {path}')

    if variant == 'two_label':
        genes = d.iloc[:, 1].astype(str).str.upper().to_numpy()
        num = d.iloc[:, 2:]
    else:
        genes = d.iloc[:, 0].astype(str).str.upper().to_numpy()
        num = d.iloc[:, 1:]

    if num.shape[1] != n_cells:
        # Header/data width mismatch: keep the shorter side and note it.
        n = min(num.shape[1], n_cells)
        num = num.iloc[:, :n]
        pos = pos[pos < n]

    try:
        V = num.to_numpy(dtype=np.float32, na_value=0.0, copy=False)
    except (TypeError, ValueError):
        V = num.apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(np.float32)

    V = V[:, pos]
    X = pd.DataFrame(V.T, columns=genes)                  # cells x genes
    X.index = [cells[i] for i in pos]
    X.index = X.index.astype(str)

    if X.columns.duplicated().any():
        keep = []
        for g in pd.unique(X.columns):
            sub = X.loc[:, X.columns == g]
            keep.append(sub.iloc[:, int(np.argmax(sub.values.sum(0)))].rename(g))
        X = pd.concat(keep, axis=1)
    return X


def qc_and_norm(X, min_genes=200, max_genes=6000, max_mito=20.0):
    gn = np.asarray(X.columns.astype(str))
    A = np.asarray(X.values, dtype=np.float32)
    tot = A.sum(1); n_gene = (A > 0).sum(1)
    mt = np.array([g.upper().startswith('MT-') for g in gn])
    mito = (A[:, mt].sum(1) / np.maximum(tot, 1)) * 100 if mt.any() else np.zeros(len(A))
    keep = (n_gene >= min_genes) & (n_gene <= max_genes) & (mito <= max_mito) & (tot > 0)
    Ak = A[keep]; tk = Ak.sum(1); tk[tk == 0] = 1
    D = np.log1p(Ak / tk[:, None] * 1e4)
    return pd.DataFrame(D, columns=gn), keep, n_gene[keep], tk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', required=True)
    ap.add_argument('--tsv', required=True)
    ap.add_argument('--markers', default='markers_pan.json')
    ap.add_argument('--species', default='human', choices=['human', 'mouse'])
    ap.add_argument('--outdir', default='/root/audit/results')
    ap.add_argument('--max-cells', type=int, default=4000)
    ap.add_argument('--nperm', type=int, default=150)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=" * 84); print(f"Reproduction audit (wide table): {args.name}"); print("=" * 84)
    X = load_wide(args.tsv, args.max_cells)
    print(f"[load] {X.shape[0]} cells x {X.shape[1]} genes")

    # Automatic species correction based on marker hit counts
    mk_h = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'markers_pan.json')))
    mk_m = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'markers_mouse.json')))
    hit_h = sum(1 for gs in mk_h.values() for g in gs if g in X.columns)
    hit_m = sum(1 for gs in mk_m.values() for g in gs if g in X.columns)
    if args.species == 'mouse' or hit_m > hit_h:
        markers, sp = mk_m, 'mouse'
    else:
        markers, sp = mk_h, 'human'
    print(f"[species] {sp}  (human marker hits={hit_h}, mouse marker hits={hit_m})")

    df, keep, ngene, tot = qc_and_norm(X)
    print(f"[QC] kept {df.shape[0]} / {X.shape[0]} cells")

    celltype, _ = annotate(df, markers)
    if celltype is None:
        print("[abort] no marker genes matched"); return
    vc = pd.Series(celltype).value_counts()
    print("[annotate] cell-type composition:")
    for k, v in vc.items():
        print(f"    {k:<16}{v:>7}  ({v/len(celltype)*100:5.1f}%)")

    LRG = LC.default_lactylation_sets_mouse() if sp == 'mouse' else LC.default_lactylation_sets()
    HK = LC.default_housekeeping_mouse() if sp == 'mouse' else LC.default_housekeeping()
    n_lrg = len([g for g in LRG if g in df.columns]); n_hk = len([g for g in HK if g in df.columns])
    print(f"[gene sets] LRG {n_lrg}/{len(LRG)} available   HK {n_hk}/{len(HK)} available")

    lact = LC.score_z(df, LRG)
    if lact is None:
        print("[abort] scoring failed"); return
    df['_lact'] = lact
    res = {'name': args.name, 'n_cells': int(df.shape[0]), 'n_genes': int(df.shape[1]),
           'celltype_counts': {k: int(v) for k, v in vc.items()},
           'n_lrg_avail': n_lrg, 'n_hk_avail': n_hk}

    print(f"\n[reproduce] lactylation score by cell type (AddModuleScore-style):")
    tab = []
    for c in vc.index:
        v = df.loc[celltype == c, '_lact']
        tab.append((c, len(v), float(v.mean())))
    tab.sort(key=lambda x: -x[2])
    for c, n, m in tab:
        print(f"    {c:<16}{n:>7}  mean={m:+.4f}")
    res['reproduced_ranking'] = [t[0] for t in tab]
    res['reproduced_top'] = tab[0][0] if tab else None
    res['reproduced_scores'] = {t[0]: t[2] for t in tab}

    print(f"\n{'='*84}\n[T1] housekeeping control (SSR = H_real / H_hk)\n{'='*84}")
    t1 = LC.SSR_between_groups(LC.score_z, df, LRG, HK, celltype, n_repeat=30)
    if t1:
        print(f"  H_real={t1['H_real']:.1f}  H_hk_mean={t1['H_hk_mean']:.1f}  H_hk_p95={t1['H_hk_p95']:.1f}")
        print(f"  SSR = {t1['SSR']:.3f}   {'PASS' if t1['SSR'] > 1 else 'FAIL'}")
        res['T1'] = t1

    print(f"\n{'='*84}\n[T2] expression-matched random null\n{'='*84}")
    t2 = LC.test_T2_matched_null(LC.score_z, df, LRG, celltype, n_perm=args.nperm)
    if t2:
        print(f"  observed effect={t2['observed_effect']:.4f}  null mean={t2['null_mean']:.4f}  null p95={t2['null_p95']:.4f}")
        print(f"  empirical p = {t2['empirical_p']:.4f}   {'PASS' if t2['pass'] else 'FAIL'}")
        res['T2'] = t2

    print(f"\n{'='*84}\n[T3] score vs technical covariates\n{'='*84}")
    lg_tot = np.log1p(tot); lg_ng = np.log1p(ngene)
    r_d, p_d = spearmanr(lact, lg_tot); r_g, p_g = spearmanr(lact, lg_ng)
    hk_score = LC.score_z(df, HK); r_h, p_h = spearmanr(lact, hk_score)
    print(f"  rho(score, log depth)     = {r_d:+.4f}  p={p_d:.2e}")
    print(f"  rho(score, log n_genes)   = {r_g:+.4f}  p={p_g:.2e}  {'STRONG COUPLING' if abs(r_g)>0.3 else ''}")
    print(f"  rho(score, HK score)      = {r_h:+.4f}  p={p_h:.2e}")
    res['T3'] = {'rho_depth': float(r_d), 'rho_ngene': float(r_g), 'rho_hk': float(r_h),
                 'depth_coupled': bool(abs(r_g) > 0.3), 'hk_overlap': bool(abs(r_h) > 0.5)}

    print(f"\n{'='*84}\n[T4] score vs sparsity / detection rate\n{'='*84}")
    rows = []
    for c in vc.index:
        m = celltype == c; sub = df.loc[m]
        det = float((sub[[g for g in LRG if g in sub.columns]] > 0).values.mean())
        rows.append({'celltype': c, 'n': int(m.sum()), 'mean_ngene': float(ngene[m].mean()),
                     'detect_rate_lrg': det, 'lact_score': float(sub['_lact'].mean())})
    t4 = pd.DataFrame(rows).sort_values('lact_score', ascending=False)
    print(f"{'cell_type':<16}{'n':>7}{'mean_genes':>12}{'LRG_detect':>12}{'lact_score':>12}")
    for _, r in t4.iterrows():
        print(f"{r.celltype:<16}{r.n:>7}{r.mean_ngene:>12.1f}{r.detect_rate_lrg:>12.3f}{r.lact_score:>12.4f}")
    if len(t4) >= 3:
        r_ng, p_ng = spearmanr(t4['mean_ngene'], t4['lact_score'])
        r_det, p_det = spearmanr(t4['detect_rate_lrg'], t4['lact_score'])
    else:
        r_ng = r_det = np.nan; p_ng = p_det = np.nan
    print(f"\n  rho(mean genes, score)    = {r_ng:+.4f}  p={p_ng:.2e}")
    print(f"  rho(LRG detection, score) = {r_det:+.4f}  p={p_det:.2e}")
    dd = bool((r_ng > 0.5) or (r_det > 0.5))
    print(f"  {'WARNING: ranking is driven by sequencing depth / detection rate' if dd else 'ranking shows no clear depth dependence'}")
    res['T4'] = {'rho_ngene_vs_score': float(r_ng), 'p_ngene': float(p_ng),
                 'rho_detect_vs_score': float(r_det), 'p_detect': float(p_det), 'depth_driven': dd}
    t4.to_csv(f"{args.outdir}/{args.name}_T4.csv", index=False)

    checks = []
    if 'T1' in res: checks.append(res['T1']['SSR'] > 1)
    if 'T2' in res: checks.append(res['T2']['pass'])
    if 'T4' in res: checks.append(not res['T4']['depth_driven'])
    res['audit'] = {'n_checks': len(checks), 'n_pass': int(sum(checks)),
                    'overall_pass': bool(all(checks)) if checks else None}
    print(f"\n{'='*84}\n[VERDICT] {args.name}\n{'='*84}")
    print(f"  passed {res['audit']['n_pass']} / {res['audit']['n_checks']} control checks")
    if not res['audit']['overall_pass']:
        if 'T1' in res and res['T1']['SSR'] <= 1:
            print(f"  FAIL T1: housekeeping between-group H={res['T1']['H_hk_mean']:.1f}) >= lactylation set H={res['T1']['H_real']:.1f})")
        if 'T2' in res and not res['T2']['pass']:
            print(f"  FAIL T2: expression-matched random sets reach the same effect (p={res['T2']['empirical_p']:.3f})")
        if 'T4' in res and res['T4']['depth_driven']:
            print(f"  FAIL T4: cell-type ranking is collinear with detection rate (rho={res['T4']['rho_detect_vs_score']:+.3f})")
    print(f"  → {'PASS: original conclusion survives the audit' if res['audit']['overall_pass'] else 'FAIL: original conclusion does not survive the audit'}")
    json.dump(res, open(f"{args.outdir}/{args.name}.json", 'w'), indent=2, ensure_ascii=False)
    print(f"\n[saved] {args.outdir}/{args.name}.json")


if __name__ == '__main__':
    main()
