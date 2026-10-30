#!/usr/bin/env python3
"""
audit_one.py -- run the "reproduction + added controls" audit on a single dataset

Usage:
  python3 audit_one.py --name <label> --matrix <mtx> --genes <genes.tsv> --barcodes <bc.tsv> \
      --markers <marker_json> [--ct-key ...]

Workflow:
  1. Load raw counts -> QC -> log1p CPM normalisation
  2. Annotate cell types with literature markers (reproducing their classification)
  3. Compute the lactylation score (AddModuleScore style, i.e. the mainstream practice under audit)
  4. Reproduce their core conclusion: which cell type has the highest lactylation score
  5. Add the T1/T2/T4 controls
  6. Export JSON + CSV
"""
import numpy as np, pandas as pd, json, argparse, os, sys, warnings
from scipy.io import mmread
from scipy.sparse import csr_matrix
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lactcheck as LC

def load_counts(matrix, genes, barcodes):
    M = mmread(matrix)
    M = csr_matrix(M)
    g = pd.read_csv(genes, sep='\t', header=None)
    b = pd.read_csv(barcodes, sep='\t', header=None)
    # Determine orientation: rows=genes or rows=cells?
    if M.shape[0] == len(g):
        X = M.T.tocsr()          # -> cells x genes
    elif M.shape[1] == len(g):
        X = M.tocsr()
    else:
        # Common case: treat the smaller dimension as the gene count
        X = M.T.tocsr() if M.shape[0] < M.shape[1] else M.tocsr()
    return X, np.array(g[1].values if g.shape[1] > 1 else g[0].values), np.array(b[0].values)

def qc_and_norm(X, gene_names, min_genes=200, max_genes=6000, max_mito=20.0):
    n_gene = np.asarray((X > 0).sum(1)).ravel()
    tot = np.asarray(X.sum(1)).ravel()
    mt = np.array([g.upper().startswith('MT-') for g in gene_names])
    mito = (np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(tot, 1)) * 100
    keep = (n_gene >= min_genes) & (n_gene <= max_genes) & (mito <= max_mito) & (tot > 0)
    Xk = X[keep]
    totk = np.asarray(Xk.sum(1)).ravel(); totk[totk == 0] = 1
    D = np.log1p(Xk.toarray().astype(np.float32) / totk[:, None] * 1e4)
    df = pd.DataFrame(D, columns=gene_names)
    return df, keep, n_gene[keep], totk

def annotate(df, markers):
    """Annotate cell types by marker scoring"""
    scores = {}
    for ct, mk in markers.items():
        v = LC.score_z(df, mk)
        if v is not None:
            scores[ct] = v
    if not scores:
        return None, None
    sdf = pd.DataFrame(scores)
    return sdf.idxmax(axis=1).values, sdf

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', required=True)
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--genes', required=True)
    ap.add_argument('--barcodes', required=True)
    ap.add_argument('--markers', required=True, help='JSON: {celltype: [genes]}')
    ap.add_argument('--outdir', default='/root/audit/results')
    ap.add_argument('--nperm', type=int, default=150)
    ap.add_argument('--species', default='human', choices=['human', 'mouse'])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=" * 84); print(f"Reproduction audit: {args.name}"); print("=" * 84)

    X, gnames, bcs = load_counts(args.matrix, args.genes, args.barcodes)
    print(f"[loaded] {X.shape[0]} cells x {X.shape[1]} genes")
    df, keep, ngene, tot = qc_and_norm(X, gnames)
    print(f"[QC] kept {df.shape[0]} / {X.shape[0]} cells")

    markers = json.load(open(args.markers))
    # Automatic species correction: if all human markers miss but the gene names look like
    # mouse, switch to the mouse marker panel
    mdir = os.path.dirname(os.path.abspath(args.markers))
    mfiles = {'human': os.path.join(mdir, 'markers_pan.json'),
              'mouse': os.path.join(mdir, 'markers_mouse.json')}
    if args.species == 'mouse' or args.markers.endswith('markers_pan.json'):
        probe = json.load(open(args.markers if args.species != 'mouse' else mfiles['mouse']))
        hit = sum(1 for gs in probe.values() for g in gs if g in df.columns)
        if args.species == 'mouse' and os.path.exists(mfiles['mouse']):
            markers = json.load(open(mfiles['mouse']))
            print(f"[species] using the mouse marker panel ({hit} probe hits)")
    celltype, sdf = annotate(df, markers)
    if celltype is None:
        print("[abort] none of the marker genes were detected")
        return
    print(f"[annotation] cell type distribution:")
    vc = pd.Series(celltype).value_counts()
    for k, v in vc.items():
        print(f"    {k:<16}{v:>7}  ({v/len(celltype)*100:5.1f}%)")

    if args.species == 'mouse':
        LRG = LC.default_lactylation_sets_mouse()
        HK = LC.default_housekeeping_mouse()
    else:
        LRG = LC.default_lactylation_sets()
        HK = LC.default_housekeeping()
    n_lrg = len([g for g in LRG if g in df.columns])
    n_hk = len([g for g in HK if g in df.columns])
    print(f"[gene set] LRG {n_lrg}/{len(LRG)} available   HK {n_hk}/{len(HK)} available")

    # ---- Reproduce the core conclusion ----
    lact = LC.score_z(df, LRG)
    res = {'name': args.name, 'n_cells': int(df.shape[0]), 'n_genes': int(df.shape[1]),
           'celltype_counts': {k: int(v) for k, v in vc.items()},
           'n_lrg_avail': n_lrg, 'n_hk_avail': n_hk}
    if lact is None:
        print("[abort] lactylation score could not be computed")
        return
    df['_lact'] = lact
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

    # ---- T1: SSR ----
    print(f"\n{'='*84}\n[T1] housekeeping-gene control (SSR = H_real / H_hk)\n{'='*84}")
    t1 = LC.SSR_between_groups(LC.score_z, df, LRG, HK, celltype, n_repeat=30)
    if t1:
        print(f"  H_real={t1['H_real']:.1f}  H_hk_mean={t1['H_hk_mean']:.1f}  H_hk_p95={t1['H_hk_p95']:.1f}")
        print(f"  SSR = {t1['SSR']:.3f}   {'✅ PASS' if t1['SSR'] > 1 else '❌ FAIL'}")
        res['T1'] = t1

    # ---- T2: expression-matched null distribution ----
    print(f"\n{'='*84}\n[T2] expression-matched random null distribution\n{'='*84}")
    t2 = LC.test_T2_matched_null(LC.score_z, df, LRG, celltype, n_perm=args.nperm)
    if t2:
        print(f"  observed effect={t2['observed_effect']:.4f}  null mean={t2['null_mean']:.4f}  95th pctile={t2['null_p95']:.4f}")
        print(f"  empirical p = {t2['empirical_p']:.4f}   {'✅ PASS' if t2['pass'] else '❌ FAIL'}")
        res['T2'] = t2

    # ---- T3: depth / gene-count confounding ----
    print(f"\n{'='*84}\n[T3] lactylation score vs technical covariates (is the score merely a proxy for sequencing depth?)\n{'='*84}")
    from scipy.stats import spearmanr
    lg_tot = np.log1p(tot); lg_ng = np.log1p(ngene)
    r_d, p_d = spearmanr(lact, lg_tot)
    r_g, p_g = spearmanr(lact, lg_ng)
    hk_score = LC.score_z(df, HK)
    r_h, p_h = spearmanr(lact, hk_score)
    print(f"  rho(lactylation score, log depth)     = {r_d:+.4f}  p={p_d:.2e}  {'⚠️ strongly coupled' if abs(r_d)>0.3 else ''}")
    print(f"  rho(lactylation score, log #genes)    = {r_g:+.4f}  p={p_g:.2e}  {'⚠️ strongly coupled' if abs(r_g)>0.3 else ''}")
    print(f"  rho(lactylation score, HK score)      = {r_h:+.4f}  p={p_h:.2e}  {'⚠️ high overlap' if abs(r_h)>0.5 else ''}")
    res['T3'] = {'rho_depth': float(r_d), 'rho_ngene': float(r_g), 'rho_hk': float(r_h),
                 'depth_coupled': bool(abs(r_d) > 0.3), 'hk_overlap': bool(abs(r_h) > 0.5)}

    # ---- T4: global diagnostic ----
    print(f"\n{'='*84}\n[T4] relationship of the score with cell-type sparsity / expression level\n{'='*84}")
    rows = []
    for c in vc.index:
        m = celltype == c
        sub = df.loc[m]
        det = float((sub[LRG[:0] if False else [g for g in LRG if g in sub.columns]] > 0).values.mean())
        rows.append({'celltype': c, 'n': int(m.sum()),
                     'mean_ngene': float(ngene[m].mean()),
                     'detect_rate_lrg': det,
                     'lact_score': float(sub['_lact'].mean())})
    t4 = pd.DataFrame(rows).sort_values('lact_score', ascending=False)
    print(f"{'cell type':<16}{'n':>7}{'mean genes':>12}{'LRG det.':>12}{'lact. score':>12}")
    for _, r in t4.iterrows():
        print(f"{r.celltype:<16}{r.n:>7}{r.mean_ngene:>12.1f}{r.detect_rate_lrg:>12.3f}{r.lact_score:>12.4f}")
    r_ng, p_ng = spearmanr(t4['mean_ngene'], t4['lact_score'])
    r_det, p_det = spearmanr(t4['detect_rate_lrg'], t4['lact_score'])
    print(f"\n  cross-cell-type rho(mean #genes, lactylation score) = {r_ng:+.4f}  p={p_ng:.2e}")
    print(f"  cross-cell-type rho(LRG det. rate, lactylation score) = {r_det:+.4f}  p={p_det:.2e}")
    # Verdict: a positive rank correlation with |rho|>=0.5 is taken as depth/detection-rate
    # driven (with n=9, p<0.1 already indicates the direction of the problem)
    depth_driven = bool((r_ng > 0.5) or (r_det > 0.5))
    print(f"  {'⚠️ score ranking at the cell-type level is mainly driven by sequencing depth / detection rate' if depth_driven else 'score ranking shows no clear relationship with depth'}")
    res['T4'] = {'rho_ngene_vs_score': float(r_ng), 'p_ngene': float(p_ng),
                 'rho_detect_vs_score': float(r_det), 'p_detect': float(p_det),
                 'depth_driven': depth_driven}
    t4.to_csv(f"{args.outdir}/{args.name}_T4.csv", index=False)

    # ---- Verdict ----
    flags = {}
    checks = []
    if 'T1' in res:
        checks.append(res['T1']['SSR'] > 1)
        flags['T1_SSR:real_le_hk' if not checks[-1] else 'T1_SSR_ok'] = True
    if 'T2' in res:
        checks.append(res['T2']['pass'])
    if 'T4' in res:
        checks.append(not res['T4']['depth_driven'])
    res['audit'] = {'n_checks': len(checks), 'n_pass': int(sum(checks)),
                    'overall_pass': bool(all(checks)) if checks else None,
                    'fails': [k for k, v in flags.items() if 'le_hk' in k]}
    print(f"\n{'='*84}\n[AUDIT VERDICT] {args.name}\n{'='*84}")
    print(f"  {res['audit']['n_pass']} / {res['audit']['n_checks']} control tests passed")
    if not res['audit']['overall_pass']:
        if 'T1' in res and res['T1']['SSR'] <= 1:
            print(f"  ❌ T1: between-group difference of housekeeping genes (H={res['T1']['H_hk_mean']:.1f}) ≥ that of the lactylation gene set (H={res['T1']['H_real']:.1f})")
        if 'T2' in res and not res['T2']['pass']:
            print(f"  ❌ T2: expression-matched random gene sets reach an equal effect (p={res['T2']['empirical_p']:.3f})")
        if 'T4' in res and res['T4']['depth_driven']:
            print(f"  ❌ T4: cell type ranking collinear with detection rate / gene count (rho={res['T4']['rho_detect_vs_score']:+.3f})")
    print(f"  → {'✅ original conclusion passes the audit' if res['audit']['overall_pass'] else '❌ original conclusion fails the audit'}")

    json.dump(res, open(f"{args.outdir}/{args.name}.json", 'w'), indent=2, ensure_ascii=False)
    print(f"\n[saved] {args.outdir}/{args.name}.json")

if __name__ == '__main__':
    main()
