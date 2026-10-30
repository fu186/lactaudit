#!/usr/bin/env python3
"""
lactcheck.py -- Specificity testing toolkit for epigenetic-modification scores (LactCheck v0.1)

Four classes of tests:
  T1 housekeeping-gene control      -> SSR (Signal Specificity Ratio)
  T2 expression-matched null        -> empirical p-value
  T3 confounder residualisation     -> residual partial correlation
  T4 cross-cell-type negative control -> specificity verdict + sparse-artifact detection

Design principles:
  - Input is always an adata-like (cells x genes) dense/sparse matrix + cell type labels + group labels
  - All scores use a mean of z-scores with "baseline = control group" (consistent with the AddModuleScore family, i.e. the common practice under audit)
  - A rank-based score (UCell style) is also provided as a benchmark, to separate "method problems" from "gene set problems"
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.linear_model import LinearRegression

# ------------------------------------------------------------------
# Scoring functions
# ------------------------------------------------------------------
def score_z(df, genes):
    """AddModuleScore style: per-gene z-scores, then averaged (the mainstream practice under audit)"""
    f = [g for g in genes if g in df.columns]
    if len(f) < 3:
        return None
    s = df[f]
    return ((s - s.mean()) / (s.std() + 1e-12)).mean(axis=1).values

def score_rank(df, genes, rmax=1500):
    """UCell style: Mann-Whitney U on within-cell ranks (insensitive to depth)"""
    f = [g for g in genes if g in df.columns]
    if len(f) < 3:
        return None
    R = df.rank(axis=1, ascending=False)
    R = R.where(R <= rmax, rmax + 1)
    sub = R[f]
    n = len(f)
    # Normalise the U statistic
    u = sub.sum(axis=1) - n * (n + 1) / 2
    return 1 - u / (n * rmax)

def score_aucell(df, genes):
    """AUCell style: AUC of the gene set within per-cell rankings"""
    return score_rank(df, genes)


# ------------------------------------------------------------------
# T1 housekeeping-gene control -> SSR
# ------------------------------------------------------------------
def test_T1_housekeeping(score_func, df, gene_set, hk_set, n_repeat=20, seed=42):
    """
    SSR = effect size of the real gene set / effect size of housekeeping genes
    Effect size definition: |Cohen's d| or |between-group mean difference| / pooled SD
    Housekeeping genes are drawn to match the size of the gene set, averaged over
    n_repeat random draws (to reduce the chance effect of a single draw)
    """
    rng = np.random.default_rng(seed)
    real = score_func(df, gene_set)
    if real is None:
        return None
    hk_pool = [g for g in hk_set if g in df.columns]
    if len(hk_pool) < 3:
        return None
    n = len([g for g in gene_set if g in df.columns])
    effs = []
    for _ in range(n_repeat):
        gs = rng.choice(hk_pool, size=min(n, len(hk_pool)), replace=False)
        v = score_func(df, list(gs))
        if v is None:
            continue
        effs.append(np.std(v))
    return {
        'real_sd': float(np.std(real)),
        'hk_sd_mean': float(np.mean(effs)) if effs else np.nan,
    }

def SSR_between_groups(score_func, df, gene_set, hk_set, group_labels,
                       n_repeat=20, seed=42):
    """
    SSR: under the given grouping (e.g. each cell type), compare "how well the real
    score separates groups" vs "how well housekeeping genes separate groups".
    Uses the Kruskal-Wallis H statistic (dimensions differ between methods, so H keeps
    them comparable), then normalises.
    Returns ssr = H_real / H_hk_mean
    """
    from scipy.stats import kruskal
    rng = np.random.default_rng(seed)
    groups = pd.Series(group_labels)
    n = len([g for g in gene_set if g in df.columns])

    def H_of(genes):
        v = score_func(df, genes)
        if v is None:
            return np.nan
        samples = [v[groups.values == g] for g in groups.unique() if (groups.values == g).sum() >= 5]
        if len(samples) < 2:
            return np.nan
        try:
            return kruskal(*samples).statistic
        except Exception:
            return np.nan

    h_real = H_of(gene_set)
    hk_pool = [g for g in hk_set if g in df.columns]
    hs = []
    for _ in range(n_repeat):
        gs = list(rng.choice(hk_pool, size=min(n, len(hk_pool)), replace=False))
        h = H_of(gs)
        if not np.isnan(h):
            hs.append(h)
    h_hk = np.mean(hs) if hs else np.nan
    return {
        'H_real': float(h_real) if not np.isnan(h_real) else np.nan,
        'H_hk_mean': float(h_hk),
        'H_hk_p95': float(np.quantile(hs, 0.95)) if hs else np.nan,
        'SSR': float(h_real / h_hk) if h_hk and h_hk > 0 and not np.isnan(h_real) else np.nan,
        'n_genes': n,
    }


# ------------------------------------------------------------------
# T2 expression-matched random null distribution
# ------------------------------------------------------------------
def test_T2_matched_null(score_func, df, gene_set, group_labels,
                         n_perm=200, seed=42, use_effect=None):
    """
    Randomly draw gene sets of equal size from an expression-matched gene pool, compute
    their scores, and test whether the effect size of the real gene set exceeds the null.
    use_effect: a function (score_vector) -> scalar; defaults to the |largest between-group mean difference|
    """
    rng = np.random.default_rng(seed)
    groups = pd.Series(group_labels)
    gs_in = [g for g in gene_set if g in df.columns]
    if len(gs_in) < 3:
        return None

    def default_effect(v):
        m = [v[groups.values == g].mean() for g in groups.unique() if (groups.values == g).sum() >= 5]
        return (max(m) - min(m)) if len(m) >= 2 else np.nan

    eff = use_effect or default_effect

    real_v = score_func(df, gs_in)
    obs = eff(real_v)

    # Expression-matched pool: means within the 5%-95% range of the real gene set
    # Note: df may contain duplicated gene column names, so column de-duplication must
    # come first before taking means
    dfe = df.loc[:, ~df.columns.duplicated()]
    means = dfe.mean(axis=0)
    gs_in_u = [g for g in gs_in if g in means.index]
    lo, hi = np.quantile(means[gs_in_u], 0.05), np.quantile(means[gs_in_u], 0.95)
    mvals = means.values
    pool = [g for g, mv in zip(means.index, mvals)
            if lo <= mv <= hi and g not in set(gs_in_u)]
    if len(pool) < len(gs_in_u) * 2:
        pool = list(means.index)

    null = np.full(n_perm, np.nan)
    for t in range(n_perm):
        gs = list(rng.choice(pool, size=len(gs_in), replace=False))
        v = score_func(df, gs)
        if v is not None:
            null[t] = eff(v)
    p_emp = (np.nansum(null >= obs) + 1) / (n_perm + 1)
    return {
        'observed_effect': float(obs),
        'null_mean': float(np.nanmean(null)),
        'null_p95': float(np.nanquantile(null, 0.95)),
        'null_max': float(np.nanmax(null)),
        'empirical_p': float(p_emp),
        'pass': bool(p_emp < 0.05),
    }


# ------------------------------------------------------------------
# T3 confounder residualisation
# ------------------------------------------------------------------
def test_T3_residual(df, x_col, y_genes, hk_set, depth=None, ngene=None):
    """
    Test whether the association between x and y (gene-set score) is independent of
    the HK score + sequencing depth + gene count

    Returns
    -------
    resid_rho / resid_p : result controlling for all supplied covariates simultaneously (primary result)
    per_covariate       : dict of rho/p when each covariate is controlled **individually**.
                          Section 2.5 of the paper states that "the number of detected genes
                          (ngene) is the pre-specified primary variable"; that claim can only
                          be verified if per_covariate['ngene'] can be read off on its own and
                          agrees in direction with the primary result. The key name is the
                          covariate name as supplied.
    primary             : name of the pre-specified primary variable ('ngene') together with
                          its individually-controlled result, so reviewers can compare
                          directly without re-running.
    """
    y = score_z(df, y_genes)
    if y is None:
        return None
    x = df[x_col].values
    hk = score_z(df, hk_set)
    rho0, p0 = spearmanr(x, y)

    covariates = {'HK': hk}
    if depth is not None:
        covariates['depth'] = depth
    if ngene is not None:
        covariates['ngene'] = ngene

    def residualised_rho(extra):
        """Rank correlation after controlling for const + HK + the specified extra covariates (extra is a list of names)."""
        cols = [np.ones(len(x)), hk]
        for nm in extra:
            if nm != 'HK':
                cols.append(np.asarray(covariates[nm]))
        Z = np.column_stack(cols)
        rx = x - LinearRegression().fit(Z, x).predict(Z)
        ry = y - LinearRegression().fit(Z, y).predict(Z)
        return spearmanr(rx, ry)

    # Primary result: control for all covariates
    r1, p1 = residualised_rho([nm for nm in covariates if nm != 'HK'])

    # Each covariate controlled individually (HK is always in the design, so the 'HK'
    # entry means controlling for HK only)
    per_cov = {}
    for nm in covariates:
        extra = [] if nm == 'HK' else [nm]
        rr, pp = residualised_rho(extra)
        per_cov[nm] = {'rho': float(rr), 'p': float(pp), 'controls': ['HK'] + extra}

    out = {
        'raw_rho': float(rho0), 'raw_p': float(p0),
        'resid_rho': float(r1), 'resid_p': float(p1),
        'attenuation': float(1 - abs(r1) / abs(rho0)) if rho0 != 0 else np.nan,
        'pass': bool(p1 < 0.05),
        'covariates': ['const'] + list(covariates.keys()),
        'per_covariate': per_cov,
    }
    if 'ngene' in per_cov:
        out['primary'] = {'name': 'ngene', **per_cov['ngene']}
    return out


# ------------------------------------------------------------------
# T4 cross-cell-type negative control + sparse-artifact detection
# ------------------------------------------------------------------
def test_T4_celltype_specificity(df, gene_set, cells_of_interest, hk_set,
                                 celltype_labels, depth=None, ngene=None,
                                 min_cells=150, sparsity_warn=0.01):
    """
    Score each cell type separately and test whether "the association within that type"
    is a sparse artifact.
    cells_of_interest: gene(s) used to define the association (e.g. TREM2)
    Returns, per cell type: positive rate of the target gene, score effect, verdict
    """
    out = []
    if cells_of_interest not in df.columns:
        return None
    ct = pd.Series(celltype_labels)
    for c in ct.unique():
        idx = np.where(ct.values == c)[0]
        if len(idx) < min_cells:
            continue
        sub = df.iloc[idx].reset_index(drop=True)
        t = sub[cells_of_interest].values
        pos_frac = float((t > 0).mean())
        if t.std() == 0:
            out.append({'celltype': c, 'n': len(idx), 'pos_frac': pos_frac,
                        'rho': np.nan, 'p': np.nan, 'resid_rho': np.nan, 'artifact': True})
            continue
        y = score_z(sub, gene_set)
        rho, p = spearmanr(t, y)
        # Residualise
        hk = score_z(sub, hk_set)
        cols = [np.ones(len(t)), hk]
        if depth is not None:
            cols.append(np.asarray(depth)[idx])
        if ngene is not None:
            cols.append(np.asarray(ngene)[idx])
        Z = np.column_stack(cols)
        rt = t - LinearRegression().fit(Z, t).predict(Z)
        ry = y - LinearRegression().fit(Z, y).predict(Z)
        rr, pr = spearmanr(rt, ry)
        out.append({
            'celltype': c, 'n': len(idx), 'pos_frac': pos_frac,
            'rho': float(rho), 'p': float(p),
            'resid_rho': float(rr), 'resid_p': float(pr),
            # Sparse-artifact verdict: when the positive rate is too low, any "strong association" should be suspected
            'artifact': bool(pos_frac < sparsity_warn and abs(rr) > 0.2),
        })
    res = pd.DataFrame(out).sort_values('resid_rho', key=lambda s: -s.abs())
    return res


# ------------------------------------------------------------------
# Overall verdict
# ------------------------------------------------------------------
def verdict(t1=None, t2=None, t3=None, t4=None):
    """
    Overall verdict. SSR>1 and T2 passes and T3 passes => PASS
    """
    flags = {}
    if t1 is not None and not np.isnan(t1.get('SSR', np.nan)):
        flags['T1_SSR_gt1'] = bool(t1['SSR'] > 1.0)
    if t2 is not None:
        flags['T2_pass'] = bool(t2.get('pass', False))
    if t3 is not None:
        flags['T3_pass'] = bool(t3.get('pass', False))
    if t4 is not None and isinstance(t4, pd.DataFrame):
        flags['T4_all_artifact'] = bool(t4['artifact'].all())
    core = [flags.get('T1_SSR_gt1'), flags.get('T2_pass'), flags.get('T3_pass')]
    core = [c for c in core if c is not None]
    return {'flags': flags, 'n_checks': len(core),
            'n_pass': sum(core), 'overall_pass': bool(all(core)) if core else None}


# ------------------------------------------------------------------
# Default gene sets
# ------------------------------------------------------------------
def default_lactylation_sets(size=336):
    """
    Reproduces the "LRG gene set" commonly used in the literature. Several papers use
    332-336 genes, sourced from MsigDB and literature curation.
    The set given here is a reproducible approximation: lactylation machinery + lactate
    metabolism + commonly reported LRG genes.
    """
    machine = ['EP300','CREBBP','KAT2A','KAT2B','KAT5','KAT6A','KAT6B','KAT7','KAT8',
               'HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HDAC8','HDAC9','HDAC10','HDAC11',
               'SIRT1','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7',
               'AARS1','AARS2','TRIM33','DPF2','NEDD4']
    meta = ['LDHA','LDHB','LDHC','SLC16A1','SLC16A2','SLC16A3','SLC16A4','SLC16A5','SLC16A6','SLC16A7',
            'SLC16A8','SLC16A9','SLC16A10','SLC16A11','SLC16A12','SLC16A13','SLC16A14',
            'PKM','PKLR','PFKP','PFKL','PFKM','HK1','HK2','HK3','GPI','ENO1','ENO2','ENO3',
            'PDK1','PDK2','PDK3','PDK4','HIF1A','HIF1AN','MPC1','MPC2','PDHA1','PDHB','PDHX','DLAT','DLD']
    targets = ['IL10','ARG1','PDHA1','IDH2','MRC1','CD163','TGFB1','SDHA','FH','ACLY','OGDH','DLST',
               'SUCLG1','MDH2','CDKN1A','CDKN2A','SERPINE1','MMP9','VEGFA','LRG1','YTHDF2','METTL3',
               'SLC7A11','NLRP3','GSDMD','CASP1','HMGB1','MECP2','RUNX1','SOX10','NR4A3','TRAP1',
               'CHI3L1','SNAI1','SCD','VCAM1','IL6','IL1B','TNF','NFKB1','RELA','STAT3','PKM2']
    return sorted(set(machine + meta + targets))


def default_housekeeping():
    return ['ACTB','GAPDH','TUBB','PPIA','RPL13A','RPS18','B2M','HPRT1','TBP','YWHAZ',
            'PGK1','EEF1A1','RPLP0','UBC','RPS3','RPL8','RPSA','RPL10','RPL11','RPL12',
            'RPS5','RPS6','RPS7','RPS9','RPL4','RPL6','RPS2','RPS4X','RPL3','RPL5',
            'ACTG1','TUBB4B','TUBA1B','PPIB','CCT3','CCT5','PSMA1','PSMB1','RAN','NAA10']


# ------------------------------------------------------------------
# Species conversion: human Symbol -> mouse Symbol (first letter upper, rest lower; the most common rule)
# ------------------------------------------------------------------
def to_mouse(symbols):
    """
    Crude human->mouse Symbol conversion: first letter upper-case + rest lower-case.
    Works for the vast majority of genes (e.g. LDHA->Ldha, SLC16A3->Slc16a3, EP300->Ep300).
    """
    out = []
    for s in symbols:
        if not s:
            continue
        out.append(s[0].upper() + s[1:].lower())
    return out


def default_lactylation_sets_mouse():
    return to_mouse(default_lactylation_sets())


def default_housekeeping_mouse():
    return to_mouse(default_housekeeping())


# ------------------------------------------------------------------
def test_T4r_depth_residualised(celltype_table, method='spearman'):
    """
    T4r -- cell-type-level depth-residualised correlation (definition in Section 2.5 of the paper).

    For each cell type, regress two quantities separately on "the mean number of detected
    genes of that type" by OLS, take the residuals of the two regressions, then compute the
    rank correlation between the residuals.

    Parameters
    ----------
    celltype_table : DataFrame that must contain the columns
        mean_ngene       mean number of detected genes per cell type
        detect_rate_lrg  LRG detection rate of that type
        lact_score       mean lactylation score of that type
    method : 'spearman' (default, consistent with the main text) or 'pearson'

    Returns
    ----
    dict(n_types, T4r, p, method, det_rate_span, residual_sd_ratio)
    When residual variance collapses (detection rate saturated) T4r is None, with an
    explanation in note.
    """
    need = {'mean_ngene', 'detect_rate_lrg', 'lact_score'}
    if not need.issubset(set(celltype_table.columns)):
        raise ValueError(f'Missing columns: {need - set(celltype_table.columns)}')
    t = celltype_table.dropna(subset=list(need))
    if len(t) < 3:
        return dict(n_types=len(t), T4r=None, p=None, method=method,
                    det_rate_span=None, residual_sd_ratio=None,
                    note='fewer than 3 cell types')

    x  = t['mean_ngene'].values.astype(float)
    y1 = t['detect_rate_lrg'].values.astype(float)
    y2 = t['lact_score'].values.astype(float)

    b1 = np.polyfit(x, y1, 1)
    b2 = np.polyfit(x, y2, 1)
    e1 = y1 - np.polyval(b1, x)
    e2 = y2 - np.polyval(b2, x)

    sd_collapse = (np.std(e1) < 1e-12) or (np.std(e2) < 1e-12)
    if sd_collapse:
        return dict(n_types=len(x), T4r=None, p=None, method=method,
                    det_rate_span=float(y1.max() - y1.min()),
                    residual_sd_ratio=None,
                    note='residual variance collapsed; detection rate saturated')

    if method == 'pearson':
        from scipy.stats import pearsonr
        r, p = pearsonr(e1, e2)
    else:
        r, p = spearmanr(e1, e2)

    return dict(n_types=len(x), T4r=float(r), p=float(p), method=method,
                det_rate_span=float(y1.max() - y1.min()),
                residual_sd_ratio=float(np.std(e1) / (np.std(y1) + 1e-12)),
                note='')


def T4r_verdict(r, min_det_span=0.15, min_types=6):
    """
    Interpretability criterion for T4r (domain of applicability declared in Sections 2.5 / 3.6 of the paper).

    T4r is interpretable only when the detection rate of the gene set has sufficient dynamic
    range across cell types; when the detection rate is near saturation (narrow span) or
    there are too few cell types, residuals are unstable and T4r inflates.
    Returns (interpretable: bool, reason: str).
    """
    if r is None or r.get('T4r') is None:
        return False, 'not computable'
    if r['n_types'] < min_types:
        return False, f"only {r['n_types']} cell types"
    if r['det_rate_span'] is not None and r['det_rate_span'] < min_det_span:
        return False, (f"detection-rate span {r['det_rate_span']:.3f} "
                       f"< {min_det_span} (saturation)")
    return True, ''
