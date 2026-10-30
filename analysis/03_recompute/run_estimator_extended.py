#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extended estimator comparison: scale the 3-dataset pilot up to the full set of
datasets that can be re-loaded, and add the depth-spread diagnostic.

For each dataset x estimator we report:
  T3   Spearman rho between the cell's score and its log detected-gene count
  T4   cross-cell-type Spearman rho between LRG detection rate and mean score
  T4r  same after residualising both on the cell type's mean detected-gene count
  top  top-ranked cell type
  det_span  detection-rate span across cell types (T4r interpretability gate)

Estimators:
  E1 score_z          manuscript's estimator (mean of per-gene z-scores)
  E2 scanpy.score_genes  Seurat-AddModuleScore style, bin-matched control draw
  E3 rank (UCell)     Mann-Whitney U over within-cell top-1500 ranks

Checkpoints after every dataset so a crash loses nothing.
"""
import gzip, os, sys, json, time
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import coo_matrix
from scipy.stats import spearmanr, mannwhitneyu, pearsonr

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '模拟'))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '06_verification'))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / '模拟'))
import audit_paths as AP            # noqa: E402

HERE = Path(__file__).resolve().parent
LIB = AP.TOOL
DATA_ROOT = AP.datadir()

sys.path.insert(0, str(LIB))
import lactcheck as LC

DATA = str(DATA_ROOT)
OUT  = str(HERE)
os.makedirs(OUT, exist_ok=True)
sc.settings.verbosity = 0

LRG = LC.default_lactylation_sets(336)


# ---------------------------------------------------------------- readers
def _isnum(x):
    try:
        float(x); return True
    except ValueError:
        return False


def read_wide_txt_gz(path, max_cells=3000, seed=0):
    """Headerless genes x cells, first column = gene symbol."""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        first = f.readline().rstrip('\n').split('\t')
        body0 = first[1:]
        headerless = bool(body0) and all(_isnum(x) for x in body0[:50])
        syms, rows = [], []
        if headerless:
            syms.append(first[0]); rows.append(body0)
        for line in f:
            p = line.rstrip('\n').split('\t')
            syms.append(p[0]); rows.append(p[1:])
    width = min(len(r) for r in rows)
    rows = [r[:width] for r in rows]
    df = pd.DataFrame(np.array(rows, dtype=np.float32), index=syms,
                      columns=[f'c{i}' for i in range(width)])
    df = df.loc[~df.index.duplicated()].T
    if df.shape[0] > max_cells:
        rng = np.random.default_rng(seed)
        df = df.loc[rng.choice(df.index, max_cells, replace=False)]
    return df


def read_10x_gsm(exdir, max_cells=3000, seed=0):
    mtxs = sorted(f for f in os.listdir(exdir) if f.endswith('matrix.mtx.gz'))
    if not mtxs:
        raise FileNotFoundError(exdir)
    mtx = os.path.join(exdir, mtxs[0])
    pre = mtx[:-len('matrix.mtx.gz')]
    genes = os.path.join(exdir, pre + 'features.tsv.gz')
    if not os.path.exists(genes):
        genes = os.path.join(exdir, pre + 'genes.tsv.gz')
    with gzip.open(mtx, 'rt') as f:
        dims = None
        for line in f:
            if line.startswith('%'):
                continue
            dims = line.split(); break
        nr, nc = int(dims[0]), int(dims[1])
        rows, cols, vals = [], [], []
        for line in f:
            if line.startswith('%'):
                continue
            a, b, c = line.split()
            rows.append(int(a)); cols.append(int(b)); vals.append(float(c))
    mn = min(min(rows), min(cols))
    if mn >= 1 and (max(rows) >= nr or max(cols) >= nc):
        rows = [x - 1 for x in rows]; cols = [x - 1 for x in cols]
    if max(rows) >= nr or max(cols) >= nc:
        raise ValueError('index out of range')
    Mi = coo_matrix((vals, (rows, cols)), shape=(nr, nc)).tocsr()   # genes x cells
    with gzip.open(genes, 'rt') as f:
        g = [l.split('\t')[1].strip() if '\t' in l else l.strip() for l in f]
    g = (g + [f'g{i}' for i in range(Mi.shape[0])])[:Mi.shape[0]]
    if Mi.shape[1] > max_cells:
        rng = np.random.default_rng(seed)
        Mi = Mi[:, np.sort(rng.choice(Mi.shape[1], max_cells, replace=False))]
    df = pd.DataFrame(Mi.T.toarray(), columns=g)
    return df.loc[:, ~df.columns.duplicated()]


# ---------------------------------------------------------------- estimators
def est_scorez(df, genes):
    return LC.score_z(df, genes)


def est_scanpy(df, genes, n_bins=25, ctrl_size=50, seed=0):
    common = [g for g in genes if g in df.columns]
    ad = sc.AnnData(df.values.astype(np.float32))
    ad.var_names = list(df.columns); ad.obs_names = list(df.index)
    sc.pp.normalize_total(ad); sc.pp.log1p(ad)
    sc.tl.score_genes(ad, common, ctrl_size=min(ctrl_size, max(1, len(common))),
                      n_bins=n_bins, random_state=seed, score_name='sg')
    return ad.obs['sg'].values


def est_rank(df, genes, n_top=1500):
    tot = df.sum(axis=1).values
    tot[tot == 0] = 1
    X = np.log1p(df.values / tot[:, None] * 1e4)
    present = [g for g in genes if g in df.columns]
    gi = [df.columns.get_loc(g) for g in present]
    out = np.full(X.shape[0], np.nan)
    for i in range(X.shape[0]):
        row = X[i]
        nz = np.nonzero(row)[0]
        if len(nz) <= len(gi):
            continue
        top = nz[np.argsort(row[nz])[-min(n_top, len(nz)):]]
        ins = np.isin(top, gi)
        a, b = top[ins], top[~ins]
        if len(a) == 0 or len(b) == 0:
            continue
        out[i] = mannwhitneyu(row[a], row[b], alternative='greater').statistic / (len(a) * len(b))
    return out


# ---------------------------------------------------------------- battery
def type_table(df, score, celltype):
    ng = (df > 0).sum(axis=1).values
    avail = [g for g in LRG if g in df.columns]
    rec = []
    for c in pd.unique(celltype):
        m = celltype == c
        if m.sum() < 5:
            continue
        sub = df.loc[m]
        det = float((sub[[g for g in avail if g in sub.columns]] > 0).values.mean())
        rec.append(dict(celltype=str(c), n=int(m.sum()),
                        mean_ngene=float(ng[m].mean()),
                        detect_rate_lrg=det,
                        lact_score=float(np.nanmean(score[m]))))
    return pd.DataFrame(rec)


def run(name, df, celltype):
    logng = np.log1p((df > 0).sum(axis=1).values)
    recs = []
    for est, fn in (('score_z', est_scorez),
                    ('scanpy.score_genes', est_scanpy),
                    ('rank_UCell', est_rank)):
        try:
            s = fn(df, LRG)
        except Exception as e:
            print(f'    [skip] {est}: {type(e).__name__}: {e}', flush=True)
            continue
        if s is None or np.all(np.isnan(s)):
            print(f'    [skip] {est}: all NaN', flush=True)
            continue
        r3 = spearmanr(s, logng, nan_policy='omit')[0]
        t = type_table(df, s, celltype)
        if len(t) < 3:
            print(f'    [skip] {est}: only {len(t)} types', flush=True)
            continue
        r4 = LC.test_T4r_depth_residualised(t)
        rec = dict(dataset=name, estimator=est,
                   n_cells=int(df.shape[0]),
                   T3=float(r3),
                   T4_rho_ngene=float(spearmanr(t['mean_ngene'], t['lact_score'])[0]),
                   T4_rho_detect=float(spearmanr(t['detect_rate_lrg'], t['lact_score'])[0]),
                   T4r=r4['T4r'], T4r_span=r4['det_rate_span'],
                   T4r_interpretable=LC.T4r_verdict(r4)[0],
                   n_types=int(len(t)),
                   top=str(t.sort_values('lact_score', ascending=False).iloc[0]['celltype']))
        recs.append(rec)
        print(f"    {est:20s} T3={r3:+.3f}  T4det={rec['T4_rho_detect']:+.3f}  "
              f"T4r={rec['T4r'] if rec['T4r'] is None else format(rec['T4r'],'+.3f')}  "
              f"span={rec['T4r_span']:.3f}  top={rec['top']}", flush=True)
    return recs


def main():
    jobs = []
    cad = os.path.join(DATA, 'gse131778', 'GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt.gz')
    if os.path.exists(cad):
        jobs.append(('GSE131778_CAD', 'wide', cad))
    for ex, nm in (('gse212966', 'GSE212966_PDAC'), ('gse166676', 'GSE166676_AAA')):
        d = os.path.join(DATA, ex, 'ex')
        if os.path.isdir(d):
            jobs.append((nm, '10x', d))
    taa = os.path.join(DATA, 'gse155468', 'ex')
    if os.path.isdir(taa):
        for gsm in sorted(os.listdir(taa)):
            if gsm.endswith('.txt.gz'):
                label = gsm.split('_')[-1].replace('.txt.gz', '')
                jobs.append((f'GSE155468_{label}', 'wide', os.path.join(taa, gsm)))

    print(f'{len(jobs)} datasets queued', flush=True)
    allrecs, done = [], set()
    for name, kind, path in jobs:
        print(f'\n=== {name} ===', flush=True)
        try:
            df = read_wide_txt_gz(path) if kind == 'wide' else read_10x_gsm(path)
        except Exception as e:
            print(f'  [FAIL] load: {type(e).__name__}: {e}', flush=True)
            continue
        print(f'  loaded {df.shape[0]} cells x {df.shape[1]} genes', flush=True)
        try:
            markers = json.load(open(str(LIB / 'markers_pan.json')))
            sdf = pd.DataFrame({k: LC.score_z(df, v) for k, v in markers.items()})
            celltype = sdf.idxmax(axis=1).values
        except Exception as e:
            print(f'  [FAIL] annotate: {type(e).__name__}: {e}', flush=True)
            continue
        try:
            recs = run(name, df, celltype)
        except Exception as e:
            print(f'  [FAIL] battery: {type(e).__name__}: {e}', flush=True)
            continue
        allrecs += recs
        done.add(name)
        with open(os.path.join(OUT, 'estimator_comparison_extended.json'), 'w') as f:
            json.dump(allrecs, f, indent=1)
        print(f'  checkpoint: {len(allrecs)} records / {len(done)} datasets', flush=True)

    print(f'\nDONE: {len(allrecs)} records from {len(done)} datasets', flush=True)
    if allrecs:
        d = pd.DataFrame(allrecs)
        print('\n=== median by estimator ===')
        print(d.groupby('estimator')[['T3', 'T4_rho_detect', 'T4r']].median().round(3))


if __name__ == '__main__':
    t0 = time.time(); main(); print(f'elapsed {time.time()-t0:.0f}s')
