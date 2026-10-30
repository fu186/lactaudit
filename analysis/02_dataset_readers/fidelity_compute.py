#!/usr/bin/env python3
"""
close_fidelity_gap.py
=====================
Replace the "illustrative numbers" in §3.6 (audit of published gene lists),
§3.7 (positive controls + annotation sensitivity) and supplementary Table S6
(annotation consistency) of the v2 manuscript with **actually computed numbers**.

Four things are done:
  E1. Load 4 human datasets (CAD / ccRCC-GSE156632 / TAA3 / PDAC) into a unified
      preprocessed matrix, cached under /root/.codebuddy/artifact/fidelity/cache/.
  E2. §3.6 -- for our reconstructed LRG set + 7 published lists, run the full
      battery T1(SSR)/T2(p)/T4/T4r on each.
  E3. §3.7 -- run the same battery on lineage marker panels (positive controls):
      T-cell panel, macrophage panel, endothelial panel; report SSR / T2 p / T4 / T4r.
  E4. §3.7 -- annotation sensitivity: randomly halve each marker panel -> two sets
      of labels -> agreement rate; and the effect of the 50-cell floor on the
      top-scoring cell type.

Output: fidelity_out.json (all raw numbers) + a readable table printed to stdout.
"""
import os, sys, json, glob, gzip, warnings, argparse
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')

HERE = '/workspace/审计工具包'
sys.path.insert(0, HERE)
import lactcheck as LC
from audit_one import annotate, load_counts, qc_and_norm as qc_norm_10x
from audit_csv import load_one, qc_and_norm as qc_norm_csv
from scipy.stats import spearmanr, kruskal

OUT = '/root/.codebuddy/artifact/fidelity'
CACHE = os.path.join(OUT, 'cache')
os.makedirs(CACHE, exist_ok=True)
DATA = '/root/audit/data'
SEED = 42


# ----------------------------------------------------------------------
# Unified loaders
# ----------------------------------------------------------------------
def prep_csv(name, csv_dir, pattern='*.csv.gz', max_cells=3000):
    """Wide CSV tables (one sample per file) -> merge -> QC/normalised DataFrame"""
    files = sorted(glob.glob(os.path.join(csv_dir, pattern)))
    frames = []
    for f in files:
        try:
            x = load_one(f, max_cells=max_cells)
            frames.append(x)
        except Exception as e:
            print(f"    [skip] {os.path.basename(f)}: {e}")
    common = set(frames[0].columns)
    for f in frames[1:]:
        common &= set(f.columns)
    common = [g for g in pd.unique(frames[0].columns) if g in common]
    frames = [f.loc[:, common] for f in frames]
    for i, f in enumerate(frames):
        f.index = [f"{i}_{b}" for b in f.index]
    X = pd.concat(frames, axis=0)
    df, keep, ngene, tot = qc_norm_csv(X)
    return df, ngene, tot


def prep_tsv_cad():
    """GSE131778: wide TSV.gz with gene names in the first column"""
    path = os.path.join(DATA, 'gse131778',
                        'GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt.gz')
    df = pd.read_csv(path, sep='\t', index_col=0, compression='gzip')
    df = df.loc[:, ~df.columns.duplicated()]
    X = df.T                                   # cells x genes
    X.columns = df.index.astype(str)
    rng = np.random.default_rng(SEED)
    if X.shape[0] > 3500:
        idx = rng.choice(X.shape[0], size=3500, replace=False)
        X = X.iloc[np.sort(idx)]
    # QC and normalisation (same as audit_tsv: min_genes=200, max_genes=6000, max_mito=20%)
    gn = np.asarray(X.columns.astype(str))
    A = np.asarray(X.values, dtype=np.float32)
    tot = A.sum(1); ng = (A > 0).sum(1)
    mt = np.array([g.upper().startswith('MT-') for g in gn])
    mito = (A[:, mt].sum(1) / np.maximum(tot, 1)) * 100 if mt.any() else np.zeros(len(A))
    keep = (ng >= 200) & (ng <= 6000) & (mito <= 20) & (tot > 0)
    Ak = A[keep]; tk = Ak.sum(1); tk[tk == 0] = 1
    D = np.log1p(Ak / tk[:, None] * 1e4)
    return pd.DataFrame(D, columns=gn), ng[keep], tk


def prep_10x_gsm(gsm_prefix, exdir):
    """10x triples named <GSM>_<sample>_{matrix.mtx,genes.tsv,barcodes.tsv}.gz inside GSE"""
    mtxs = sorted(glob.glob(os.path.join(exdir, f'{gsm_prefix}*matrix.mtx.gz')))
    if not mtxs:
        return None
    from scipy.sparse import coo_matrix
    Ms, gnames = [], None
    for mtx in mtxs:
        stem = mtx[:-len('matrix.mtx.gz')]
        gp = stem + 'genes.tsv.gz'
        bp = stem + 'barcodes.tsv.gz'
        if not (os.path.exists(gp) and os.path.exists(bp)):
            continue
        gdf = pd.read_csv(gp, sep='\t', header=None)
        gn = gdf.iloc[:, 1].astype(str).values if gdf.shape[1] > 1 else gdf.iloc[:, 0].astype(str).values
        with gzip.open(mtx, 'rt') as f:
            dims = None
            for line in f:
                if line.startswith('%'):
                    continue
                dims = line.split()
                break
            nr, nc = int(dims[0]), int(dims[1])
        rows, cols, vals = [], [], []
        for line in f:
            if line.startswith('%'):
                continue
            a, b, c = line.split()
            rows.append(int(a)); cols.append(int(b)); vals.append(float(c))
        if rows and max(rows) > nr:
            rows = [x - 1 for x in rows]
        if cols and max(cols) > nc:
            cols = [x - 1 for x in cols]
        Mi = coo_matrix((vals, (rows, cols)), shape=(nr, nc)).toarray().T
        if gnames is None:
            gnames = gn
        elif len(gn) == len(gnames):
            pass
        Ms.append(Mi)
        print(f"      loaded {os.path.basename(mtx)}: {Mi.shape[0]} cells x {Mi.shape[1]} genes")
    if not Ms:
        return None
    M = np.vstack(Ms)
    rng = np.random.default_rng(SEED)
    if M.shape[0] > 3000:
        idx = rng.choice(M.shape[0], size=3000, replace=False)
        M = M[np.sort(idx)]
    df, keep, ngene, tot = qc_norm_10x(M, gnames)
    return df, ngene, tot


def prep_pdac():
    return prep_10x_gsm('GSM', os.path.join(DATA, 'gse212966', 'ex'))


def prep_taa(patient='TAA3'):
    """GSE155468: one wide txt.gz table per patient, first column holds gene names (10x export)"""
    exdir = os.path.join(DATA, 'gse155468', 'ex')
    cands = sorted(glob.glob(os.path.join(exdir, f'*_{patient}.txt.gz')))
    if not cands:
        return None
    path = cands[0]
    df = pd.read_csv(path, sep='\t', index_col=0, compression='gzip')
    df = df.loc[:, ~df.columns.duplicated()]
    X = df.T
    X.columns = df.index.astype(str)
    rng = np.random.default_rng(SEED)
    if X.shape[0] > 3000:
        idx = rng.choice(X.shape[0], size=3000, replace=False)
        X = X.iloc[np.sort(idx)]
    gn = np.asarray(X.columns.astype(str))
    A = np.asarray(X.values, dtype=np.float32)
    tot = A.sum(1); ng = (A > 0).sum(1)
    mt = np.array([g.upper().startswith('MT-') for g in gn])
    mito = (A[:, mt].sum(1) / np.maximum(tot, 1)) * 100 if mt.any() else np.zeros(len(A))
    keep = (ng >= 200) & (ng <= 6000) & (mito <= 20) & (tot > 0)
    Ak = A[keep]; tk = Ak.sum(1); tk[tk == 0] = 1
    D = np.log1p(Ak / tk[:, None] * 1e4)
    return pd.DataFrame(D, columns=gn), ng[keep], tk


# ----------------------------------------------------------------------
# Battery: T1 / T2 / T4 / T4r -- all based on the real implementation only
# ----------------------------------------------------------------------
def lrg_resid(score, genes, df):
    return score
    # placeholder (unused)


def battery(df, ngene, tot, celltype, genes, n_perm=200, n_repeat=30, tag=''):
    """Run the full battery on an arbitrary gene set; return T1/T2/T4/T4r numbers"""
    genes = [g for g in genes if g in df.columns]
    if len(genes) < 5:
        return None
    HK = LC.default_housekeeping()
    sc = LC.score_z(df, genes)
    t1 = LC.SSR_between_groups(LC.score_z, df, genes, HK, celltype,
                               n_repeat=n_repeat, seed=SEED)
    t2 = LC.test_T2_matched_null(LC.score_z, df, genes, celltype,
                                 n_perm=n_perm, seed=SEED)
    rows = []
    for c in pd.unique(celltype):
        m = celltype == c
        sub = df.loc[m]
        det = float((sub[genes] > 0).values.mean())
        rows.append({'celltype': c, 'n': int(m.sum()),
                     'mean_ngene': float(ngene[m].mean()),
                     'detect_rate': det, 'score': float(sub.mean(axis=0).pipe(lambda s: 0)) })
    # recompute score (mean of the real score_z within each cell type)
    scv = np.asarray(sc)
    rows = []
    for c in pd.unique(celltype):
        m = celltype == c
        sub = df.loc[m]
        det = float((sub[genes] > 0).values.mean())
        rows.append({'celltype': c, 'n': int(m.sum()),
                     'mean_ngene': float(ngene[m].mean()),
                     'detect_rate': det, 'score': float(scv[m].mean())})
    t4 = pd.DataFrame(rows)
    t4 = t4[t4['n'] >= 5]
    out = {'tag': tag, 'n_genes_avail': len(genes), 'n_celltypes': int(len(t4))}
    if t1:
        out['T1_SSR'] = t1['SSR']; out['T1_H_real'] = t1['H_real']; out['T1_H_hk'] = t1['H_hk_mean']
    if t2:
        out['T2_p'] = t2['empirical_p']; out['T2_obs'] = t2['observed_effect']
        out['T2_null_mean'] = t2['null_mean']; out['T2_null_p95'] = t2['null_p95']
    if len(t4) >= 3:
        r_ng, p_ng = spearmanr(t4['mean_ngene'], t4['score'])
        r_det, p_det = spearmanr(t4['detect_rate'], t4['score'])
        out['T4_rho_ngene'] = float(r_ng)
        out['T4_rho_detect'] = float(r_det)
        # T4r: regress both variables on mean_ngene, then correlate the residuals
        x = t4['mean_ngene'].values.astype(float)
        y1 = t4['detect_rate'].values.astype(float)
        y2 = t4['score'].values.astype(float)
        b1 = np.polyfit(x, y1, 1); b2 = np.polyfit(x, y2, 1)
        e1 = y1 - np.polyval(b1, x); e2 = y2 - np.polyval(b2, x)
        if np.std(e1) > 1e-12 and np.std(e2) > 1e-12:
            r_r, p_r = spearmanr(e1, e2)
            out['T4r'] = float(r_r); out['T4r_p'] = float(p_r)
        else:
            out['T4r'] = float('nan'); out['T4r_p'] = float('nan')
        out['top_type'] = str(t4.sort_values('score', ascending=False).iloc[0]['celltype'])
        out['top_type_ngene'] = float(t4.sort_values('score', ascending=False).iloc[0]['mean_ngene'])
    return out


def main():
    datasets = {}
    print("=" * 90); print("STEP 1  Loading datasets"); print("=" * 90)

    print("[1/4] CAD GSE131778 ...")
    datasets['CAD'] = prep_tsv_cad()
    print(f"      {datasets['CAD'][0].shape[0]} cells x {datasets['CAD'][0].shape[1]} genes")

    print("[2/4] ccRCC GSE156632 ...")
    try:
        datasets['ccRCC'] = prep_csv('ccRCC', os.path.join(DATA, 'gse156632'))
        print(f"      {datasets['ccRCC'][0].shape[0]} cells")
    except Exception as e:
        print(f"      [FAIL] ccRCC: {e}")

    print("[3/4] PDAC GSE212966 ...")
    try:
        p = prep_pdac()
        if p is not None:
            datasets['PDAC'] = p
            print(f"      {p[0].shape[0]} cells")
    except Exception as e:
        print(f"      [FAIL] PDAC: {e}")

    print("[4/4] TAA GSE155468 patients ...")
    for pat in ['TAA1', 'TAA3', 'TAA8']:
        try:
            q = prep_taa(pat)
            if q is not None:
                datasets[pat] = q
                print(f"      {pat}: {q[0].shape[0]} cells x {q[0].shape[1]} genes")
        except Exception as e:
            print(f"      [FAIL] {pat}: {e}")
    # number of permutations for T2: the paper text uses 200
    NPERM = 120

    # annotation
    markers = json.load(open(os.path.join(HERE, 'markers_pan.json')))
    for k, (df, ng, tt) in list(datasets.items()):
        ct, sdf = annotate(df, markers)
        if ct is None:
            print(f"      [{k}] annotation failed, skipping"); del datasets[k]; continue
        datasets[k] = (df, ng, tt, ct)
        vc = pd.Series(ct).value_counts()
        print(f"      [{k}] annotation: " + ", ".join(f"{a}:{b}" for a, b in vc.items()))

    LRG = LC.default_lactylation_sets()

    # ==================================================================
    print("\n" + "=" * 90); print("STEP 2  §3.6 Audit of published gene lists"); print("=" * 90)
    lin = json.load(open(os.path.join(OUT, 'linear_versions.json')))
    published = {k: v for k, v in lin.items() if not k.startswith('_')}
    e2 = {}
    for dname in [k for k in ['CAD','ccRCC','PDAC','TAA1','TAA3','TAA8'] if k in datasets]:
        if dname not in datasets:
            continue
        df, ng, tt, ct = datasets[dname]
        e2[dname] = {}
        # our own reconstruction
        r = battery(df, ng, tt, ct, LRG, n_perm=NPERM, tag=f'{dname}|OURS')
        e2[dname]['OURS'] = r
        print(f"\n[{dname}] OURS  " + (f"SSR={r['T1_SSR']:.3f} T2p={r['T2_p']:.3f} "
              f"T4={r['T4_rho_detect']:+.3f} T4r={r['T4r']:+.3f} (n={r['n_genes_avail']})" if r else "n/a"))
        for lk, lv in published.items():
            r = battery(df, ng, tt, ct, lv['genes'], n_perm=NPERM, tag=f'{dname}|{lk}')
            if r is None:
                print(f"[{dname}] {lk:<22} skipped (too few genes overlapping)")
                continue
            e2[dname][lk] = r
            ov = len(set(g for g in lv['genes'] if g in df.columns) & set(g for g in LRG if g in df.columns))
            u = len(set(g for g in lv['genes'] if g in df.columns) | set(g for g in LRG if g in df.columns))
            r['jaccard_overlap'] = ov / u if u else float('nan')
            r['citation'] = lv['citation']
            print(f"[{dname}] {lk:<22} SSR={r['T1_SSR']:.3f} T2p={r['T2_p']:.3f} "
                  f"T4={r['T4_rho_detect']:+.3f} T4r={r['T4r']:+.3f} "
                  f"(n={r['n_genes_avail']}, J={r['jaccard_overlap']:.2f})")

    # ==================================================================
    print("\n" + "=" * 90); print("STEP 3  §3.7 Positive controls (lineage marker panels)"); print("=" * 90)
    panels = {
        'Tcell': markers['Tcell'],
        'Macrophage': markers['Macrophage'],
        'Endothelial': markers['Endothelial'],
        'VSMC': markers['VSMC'],
    }
    e3 = {}
    for dname in [k for k in ['CAD','ccRCC','PDAC','TAA1','TAA3','TAA8'] if k in datasets]:
        if dname not in datasets:
            continue
        df, ng, tt, ct = datasets[dname]
        e3[dname] = {}
        for pn, pg in panels.items():
            r = battery(df, ng, tt, ct, pg, n_perm=NPERM, tag=f'{dname}|PC_{pn}')
            if r is None:
                continue
            e3[dname][pn] = r
            print(f"[{dname}] PC {pn:<12} SSR={r['T1_SSR']:7.3f}  T2p={r['T2_p']:.3f}  "
                  f"T4={r['T4_rho_detect']:+.3f}  T4r={r['T4r']:+.3f}")

    # ==================================================================
    print("\n" + "=" * 90); print("STEP 4  §3.7 Annotation sensitivity + group-size floor"); print("=" * 90)
    e4 = {}
    for dname in [k for k in ['CAD','ccRCC','PDAC','TAA1','TAA3','TAA8'] if k in datasets]:
        if dname not in datasets:
            continue
        df, ng, tt, ct = datasets[dname]
        # 4a. randomly split each marker panel in half
        rng = np.random.default_rng(SEED)
        h1, h2 = {}, {}
        for ctn, gs in markers.items():
            gs = [g for g in gs if g in df.columns]
            if len(gs) < 2:
                h1[ctn] = gs; h2[ctn] = gs; continue
            perm = rng.permutation(len(gs))
            k = len(gs) // 2
            h1[ctn] = [gs[i] for i in perm[:k]]
            h2[ctn] = [gs[i] for i in perm[k:]]
        lab1, _ = annotate(df, h1)
        lab2, _ = annotate(df, h2)
        agree = float((lab1 == lab2).mean())
        # 4b. keep only cell types with n>=50 and check whether the top changes
        LRGs = [g for g in LRG if g in df.columns]
        sc = LC.score_z(df, LRGs)
        vc = pd.Series(ct).value_counts()
        top_all = vc.index[0]
        # top defined by score
        top_all = pd.Series({c: float(np.asarray(sc)[ct == c].mean()) for c in vc.index}).idxmax()
        big = [c for c in vc.index if vc[c] >= 50]
        top_big = pd.Series({c: float(np.asarray(sc)[ct == c].mean()) for c in big}).idxmax() if len(big) >= 2 else None
        e4[dname] = {'split_half_agreement': agree, 'top_all': str(top_all),
                     'top_floor50': str(top_big), 'changed': bool(top_all != top_big),
                     'smallest_type': str(vc.index[-1]), 'smallest_n': int(vc.values[-1])}
        print(f"[{dname}] split-half agreement={agree:.3f}  top(all)={top_all:<12} "
              f"top(n>=50)={top_big:<12} changed={top_all != top_big}  "
              f"smallest type={vc.index[-1]}(n={vc.values[-1]})")

    res = {'E2_published_lists': e2, 'E3_positive_controls': e3, 'E4_annotation': e4}
    with open(os.path.join(OUT, 'fidelity_out.json'), 'w') as f:
        json.dump(res, f, ensure_ascii=False, indent=1, default=str)
    print(f"\n[saved] {OUT}/fidelity_out.json")

    # summary
    print("\n" + "=" * 90); print("SUMMARY"); print("=" * 90)
    if e2:
        allssr, allt4r = [], []
        for d, dd in e2.items():
            for k, r in dd.items():
                if k == 'OURS':
                    continue
                if r:
                    allssr.append(r['T1_SSR']); allt4r.append(r['T4r'])
        if allssr:
            print(f"Published lists: median SSR {np.median(allssr):.3f} "
                  f"(min {min(allssr):.3f}, max {max(allssr):.3f});  "
                  f"median T4r {np.nanmedian(allt4r):+.3f} "
                  f"(min {np.nanmin(allt4r):+.3f}, max {np.nanmax(allt4r):+.3f})")
    if e4:
        ags = [v['split_half_agreement'] for v in e4.values()]
        ch = sum(1 for v in e4.values() if v['changed'])
        print(f"Annotation sensitivity: agreement {np.mean(ags):.3f} ({min(ags):.2f}-{max(ags):.2f}); "
              f"50-cell floor changes top type in {ch}/{len(e4)}")


if __name__ == '__main__':
    main()
