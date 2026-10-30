#!/usr/bin/env python3
"""
audit_prescored.py -- runs the same T1/T2/T4 audit on an already-scored cell table (e.g. our own processing of GSE159677).
Used to place an in-house dataset on the same yardstick as the external papers.

Usage:
  python3 audit_prescored.py --name AS_GSE159677_ours \
      --csv /root/as_test/as_cells_with_region.csv \
      --score-col META --hk-col CTRL_HK --ct-col celltype \
      --gene-matrix <optional: cells x genes csv> --lrg-json <optional>
"""
import numpy as np, pandas as pd, json, argparse, os, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lactcheck as LC
from scipy.stats import spearmanr, kruskal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', required=True)
    ap.add_argument('--csv', required=True)
    ap.add_argument('--score-col', default='META')
    ap.add_argument('--hk-col', default='CTRL_HK')
    ap.add_argument('--ct-col', default='celltype')
    ap.add_argument('--outdir', default='/root/audit/results')
    ap.add_argument('--nperm', type=int, default=200)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=" * 84); print(f"Reproduction audit (pre-scored table): {args.name}"); print("=" * 84)

    d = pd.read_csv(args.csv)
    print(f"[loaded] {d.shape[0]} cells, columns: {list(d.columns)[:14]} ...")
    ct = d[args.ct_col].values
    y = d[args.score_col].values.astype(float)
    hk = d[args.hk_col].values.astype(float)
    vc = pd.Series(ct).value_counts()

    print(f"\n[reproduction] {args.score_col} score per cell type:")
    tab = [(c, int((ct == c).sum()), float(y[ct == c].mean())) for c in vc.index]
    tab.sort(key=lambda x: -x[2])
    for c, n, m in tab:
        print(f"    {c:<16}{n:>7}  mean={m:+.4f}")

    res = {'name': args.name, 'n_cells': int(len(d)), 'score_col': args.score_col,
           'celltype_counts': {k: int(v) for k, v in vc.items()},
           'reproduced_ranking': [t[0] for t in tab],
           'reproduced_scores': {t[0]: t[2] for t in tab}}

    # ---- T1: H(real) vs H(hk) ----
    print(f"\n{'='*84}\n[T1] housekeeping control (SSR = H_real / H_hk)\n{'='*84}")
    groups = [g for g in vc.index if (ct == g).sum() >= 10]
    H_real = kruskal(*[y[ct == g] for g in groups]).statistic
    # Use HK as the stand-in for a "real gene set" and compute its between-group H directly
    H_hk = kruskal(*[hk[ct == g] for g in groups]).statistic
    SSR = H_real / H_hk if H_hk > 0 else np.nan
    print(f"  H_real={H_real:.1f}   H_hk={H_hk:.1f}")
    print(f"  SSR = {SSR:.3f}   {'✅ PASS (real gene set differs more between groups than housekeeping)' if SSR > 1 else '❌ FAIL (between-group difference of housekeeping genes is not lower than that of the lactylation gene set)'}")
    res['T1'] = {'H_real': float(H_real), 'H_hk': float(H_hk), 'SSR': float(SSR), 'pass': bool(SSR > 1)}

    # ---- T2: expression-matched null (here an empirical null built from HK and several control pathways) ----
    print(f"\n{'='*84}\n[T2] effect-size comparison against built-in control pathways\n{'='*84}")
    ctrl_cols = [c for c in d.columns if c.startswith('CTRL_')]
    if ctrl_cols:
        def eff(v):
            m = [v[ct == g].mean() for g in groups]
            return max(m) - min(m)
        obs = eff(y)
        null = np.array([eff(d[c].values.astype(float)) for c in ctrl_cols])
        p = (np.sum(null >= obs) + 1) / (len(null) + 1)
        print(f"  observed effect (max between-group mean difference) = {obs:.4f}")
        for c, v in zip(ctrl_cols, null):
            print(f"    control {c:<14} effect = {v:.4f}   {'⚠️ ≥ observed' if v >= obs else ''}")
        print(f"  empirical p = {p:.4f}   {'✅ PASS' if p < 0.05 else '❌ FAIL (a control set has an effect not lower than the target set)'}")
        res['T2'] = {'observed': float(obs), 'null': {c: float(v) for c, v in zip(ctrl_cols, null)},
                     'empirical_p': float(p), 'pass': bool(p < 0.05)}

    # ---- T3: is the score merely a depth proxy? ----
    print(f"\n{'='*84}\n[T3] score vs technical covariates\n{'='*84}")
    for cov in [c for c in ('nCount', 'nFeature', 'ngene', 'depth') if c in d.columns]:
        r, pv = spearmanr(y, d[cov].values)
        print(f"  rho(score, {cov}) = {r:+.4f}  p={pv:.2e}  {'⚠️ strongly coupled' if abs(r) > 0.3 else ''}")
    r_h, p_h = spearmanr(y, hk)
    print(f"  rho(score, HK score) = {r_h:+.4f}  p={p_h:.2e}  {'⚠️ high overlap' if abs(r_h) > 0.5 else ''}")
    res['T3'] = {'rho_hk': float(r_h), 'p_hk': float(p_h), 'hk_overlap': bool(abs(r_h) > 0.5)}

    json.dump(res, open(f"{args.outdir}/{args.name}.json", 'w'), indent=2, ensure_ascii=False)
    print(f"\n[saved] {args.outdir}/{args.name}.json")


if __name__ == '__main__':
    main()
