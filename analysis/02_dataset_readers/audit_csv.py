#!/usr/bin/env python3
"""
audit_csv.py -- reproduction audit for the wide-table CSV format "Gene_ID / Symbol / barcode..."
(a common shape when 10x CellRanger output has been exported to CSV by the authors, e.g. GSE156632)

Usage:
  python3 audit_csv.py --name <label> --csv-dir <dir> --markers markers_pan.json \
      [--pattern '*.csv.gz'] [--max-cells 3000] [--nperm 150]
"""
import numpy as np, pandas as pd, json, argparse, os, sys, glob, gzip, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lactcheck as LC
from scipy.stats import spearmanr


def load_one(path, max_cells=3000, seed=0):
    """Read a single CSV -> returns genes(list) and a cells x genes DataFrame"""
    op = gzip.open if path.endswith('.gz') else open
    df = pd.read_csv(op(path), index_col=0)
    sym = None
    if 'Symbol' in df.columns:
        sym = df['Symbol'].astype(str).values
        df = df.drop(columns=[c for c in ('Gene_ID', 'Symbol') if c in df.columns])
    else:
        for c in ('Gene_Name', 'gene', 'symbol'):
            if c in df.columns:
                sym = df[c].astype(str).values
                df = df.drop(columns=[c]); break
    if sym is None:
        sym = df.index.astype(str).values
        # If the index is an Ensembl ID and a symbol column exists, capture it from the column names
    X = df.T                      # cells x genes
    X.columns = sym
    # The same Symbol may appear more than once (multiple Ensembl IDs mapping to it); keep only the highest-expressed column
    if X.columns.duplicated().any():
        keep_cols = []
        for g in pd.unique(X.columns):
            sub = X.loc[:, X.columns == g]
            keep_cols.append(sub.iloc[:, int(np.argmax(sub.values.sum(0)))].rename(g))
        X = pd.concat(keep_cols, axis=1)
    if max_cells and X.shape[0] > max_cells:
        rng = np.random.default_rng(seed)
        idx = rng.choice(X.shape[0], size=max_cells, replace=False)
        X = X.iloc[np.sort(idx)]
    return X


def load_counts_csv(csv_dir, pattern='*.csv.gz', max_cells=3000):
    files = sorted(glob.glob(os.path.join(csv_dir, pattern)))
    if not files:
        files = sorted(glob.glob(os.path.join(csv_dir, '*.csv')))
    frames = []
    for f in files:
        try:
            x = load_one(f, max_cells=max_cells)
            frames.append(x)
            print(f"    [{os.path.basename(f)}] {x.shape[0]} cells x {x.shape[1]} genes")
        except Exception as e:
            print(f"    [skipped] {os.path.basename(f)}: {e}")
    if not frames:
        return None
    # Take the intersection of genes so that columns line up
    common = set(frames[0].columns)
    for f in frames[1:]:
        common &= set(f.columns)
    common = [g for g in pd.unique(frames[0].columns) if g in common]
    frames = [f.loc[:, common] for f in frames]
    # Force a unique index before merging, to avoid clashing cell barcodes with the same name
    for i, f in enumerate(frames):
        f.index = [f"{i}_{b}" for b in f.index]
    return pd.concat(frames, axis=0)


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
    ap.add_argument('--csv-dir', required=True)
    ap.add_argument('--pattern', default='*.csv.gz')
    ap.add_argument('--markers', default='markers_pan.json')
    ap.add_argument('--outdir', default='/root/audit/results')
    ap.add_argument('--max-cells', type=int, default=3000)
    ap.add_argument('--nperm', type=int, default=150)
    args = ap.parse_args()

    from audit_one import annotate
    os.makedirs(args.outdir, exist_ok=True)
    print("=" * 84); print(f"Reproduction audit (CSV): {args.name}"); print("=" * 84)

    X = load_counts_csv(args.csv_dir, args.pattern, args.max_cells)
    if X is None:
        print("[abort] no readable CSV")
        return
    print(f"[loaded] {X.shape[0]} cells x {X.shape[1]} genes")
    df, keep, ngene, tot = qc_and_norm(X)
    print(f"[QC] kept {df.shape[0]} / {X.shape[0]} cells")

    markers = json.load(open(args.markers))
    celltype, _ = annotate(df, markers)
    if celltype is None:
        print("[abort] no marker hit"); return
    vc = pd.Series(celltype).value_counts()
    print("[annotation] cell type distribution:")
    for k, v in vc.items():
        print(f"    {k:<16}{v:>7}  ({v/len(celltype)*100:5.1f}%)")

    LRG = LC.default_lactylation_sets(); HK = LC.default_housekeeping()
    n_lrg = len([g for g in LRG if g in df.columns]); n_hk = len([g for g in HK if g in df.columns])
    print(f"[gene set] LRG {n_lrg}/{len(LRG)} available   HK {n_hk}/{len(HK)} available")

    lact = LC.score_z(df, LRG)
    if lact is None:
        print("[abort] scoring failed"); return
    df['_lact'] = lact
    res = {'name': args.name, 'n_cells': int(df.shape[0]), 'n_genes': int(df.shape[1]),
           'celltype_counts': {k: int(v) for k, v in vc.items()},
           'n_lrg_avail': n_lrg, 'n_hk_avail': n_hk}

    print(f"\n[reproduction] lactylation score per cell type (AddModuleScore style):")
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

    print(f"\n{'='*84}\n[T1] housekeeping-gene control (SSR = H_real / H_hk)\n{'='*84}")
    t1 = LC.SSR_between_groups(LC.score_z, df, LRG, HK, celltype, n_repeat=30)
    if t1:
        print(f"  H_real={t1['H_real']:.1f}  H_hk_mean={t1['H_hk_mean']:.1f}  H_hk_p95={t1['H_hk_p95']:.1f}")
        print(f"  SSR = {t1['SSR']:.3f}   {'✅ PASS' if t1['SSR'] > 1 else '❌ FAIL'}")
        res['T1'] = t1

    print(f"\n{'='*84}\n[T2] expression-matched random null distribution\n{'='*84}")
    t2 = LC.test_T2_matched_null(LC.score_z, df, LRG, celltype, n_perm=args.nperm)
    if t2:
        print(f"  observed effect={t2['observed_effect']:.4f}  null mean={t2['null_mean']:.4f}  95th pctile={t2['null_p95']:.4f}")
        print(f"  empirical p = {t2['empirical_p']:.4f}   {'✅ PASS' if t2['pass'] else '❌ FAIL'}")
        res['T2'] = t2

    print(f"\n{'='*84}\n[T3] score vs technical covariates\n{'='*84}")
    lg_tot = np.log1p(tot); lg_ng = np.log1p(ngene)
    r_d, p_d = spearmanr(lact, lg_tot); r_g, p_g = spearmanr(lact, lg_ng)
    hk_score = LC.score_z(df, HK); r_h, p_h = spearmanr(lact, hk_score)
    print(f"  rho(score, log depth)      = {r_d:+.4f}  p={p_d:.2e}")
    print(f"  rho(score, log #genes)     = {r_g:+.4f}  p={p_g:.2e}  {'⚠️ strongly coupled' if abs(r_g)>0.3 else ''}")
    print(f"  rho(score, HK score)       = {r_h:+.4f}  p={p_h:.2e}  {'⚠️ high overlap' if abs(r_h)>0.5 else ''}")
    res['T3'] = {'rho_depth': float(r_d), 'rho_ngene': float(r_g), 'rho_hk': float(r_h),
                 'depth_coupled': bool(abs(r_g) > 0.3), 'hk_overlap': bool(abs(r_h) > 0.5)}

    print(f"\n{'='*84}\n[T4] relationship between the score and sparsity / detection rate\n{'='*84}")
    rows = []
    for c in vc.index:
        m = celltype == c; sub = df.loc[m]
        det = float((sub[[g for g in LRG if g in sub.columns]] > 0).values.mean())
        rows.append({'celltype': c, 'n': int(m.sum()), 'mean_ngene': float(ngene[m].mean()),
                     'detect_rate_lrg': det, 'lact_score': float(sub['_lact'].mean())})
    t4 = pd.DataFrame(rows).sort_values('lact_score', ascending=False)
    print(f"{'cell type':<16}{'n':>7}{'mean genes':>12}{'LRG det.':>12}{'lact. score':>12}")
    for _, r in t4.iterrows():
        print(f"{r.celltype:<16}{r.n:>7}{r.mean_ngene:>12.1f}{r.detect_rate_lrg:>12.3f}{r.lact_score:>12.4f}")
    if len(t4) >= 3:
        r_ng, p_ng = spearmanr(t4['mean_ngene'], t4['lact_score'])
        r_det, p_det = spearmanr(t4['detect_rate_lrg'], t4['lact_score'])
    else:
        r_ng = r_det = np.nan; p_ng = p_det = np.nan
    print(f"\n  rho(mean #genes, score)   = {r_ng:+.4f}  p={p_ng:.2e}")
    print(f"  rho(LRG det. rate, score) = {r_det:+.4f}  p={p_det:.2e}")
    dd = bool((r_ng > 0.5) or (r_det > 0.5))
    res['T4'] = {'rho_ngene_vs_score': float(r_ng), 'p_ngene': float(p_ng),
                 'rho_detect_vs_score': float(r_det), 'p_detect': float(p_det), 'depth_driven': dd}
    t4.to_csv(f"{args.outdir}/{args.name}_T4.csv", index=False)

    checks = []
    if 'T1' in res: checks.append(res['T1']['SSR'] > 1)
    if 'T2' in res: checks.append(res['T2']['pass'])
    if 'T4' in res: checks.append(not res['T4']['depth_driven'])
    res['audit'] = {'n_checks': len(checks), 'n_pass': int(sum(checks)),
                    'overall_pass': bool(all(checks)) if checks else None}
    print(f"\n{'='*84}\n[AUDIT VERDICT] {args.name}\n{'='*84}")
    print(f"  {res['audit']['n_pass']} / {res['audit']['n_checks']} control tests passed")
    if not res['audit']['overall_pass']:
        if 'T1' in res and res['T1']['SSR'] <= 1:
            print(f"  ❌ T1: between-group difference of housekeeping genes (H={res['T1']['H_hk_mean']:.1f}) ≥ that of the lactylation gene set (H={res['T1']['H_real']:.1f})")
        if 'T2' in res and not res['T2']['pass']:
            print(f"  ❌ T2: expression-matched random gene sets reach an equal effect (p={res['T2']['empirical_p']:.3f})")
        if 'T4' in res and res['T4']['depth_driven']:
            print(f"  ❌ T4: cell type ranking collinear with detection rate (rho={res['T4']['rho_detect_vs_score']:+.3f})")
    print(f"  → {'✅ original conclusion passes the audit' if res['audit']['overall_pass'] else '❌ original conclusion fails the audit'}")
    json.dump(res, open(f"{args.outdir}/{args.name}.json", 'w'), indent=2, ensure_ascii=False)
    print(f"\n[saved] {args.outdir}/{args.name}.json")


if __name__ == '__main__':
    main()
