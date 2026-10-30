"""
The audit itself.

``audit_gene_set`` is the one function most users need.  It takes an expression
matrix plus per-cell labels and returns every statistic reported in the paper:
the size-matched null, the specificity ratio under three panel-size regimes,
depth coupling, and the within-cell-type check that distinguishes a genuine
signal from a depth artefact.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Sequence

import numpy as np
import pandas as pd

from .genesets import HOUSEKEEPING_SET, HOUSEKEEPING_POOL
from .scoring import (
    score_add_module,
    score_pearson_residual,
    score_rank,
    SCORERS,
)
from .stats import kruskal_h, signal_specificity_ratio, cohen_d, spearman_depth

__all__ = ["size_matched_null", "audit_gene_set", "AuditResult", "cluster_mat"]


def size_matched_null(matrix, gene_index, labels, h_observed,
                      n_control: int | None = None, n_repeat: int = 200,
                      seed: int = 42):
    """Build the null distribution of H by repeatedly scoring random gene sets.

    Parameters
    ----------
    gene_index : sequence of int
        Columns of the *control pool* to draw from.  This should already be
        restricted to housekeeping genes if you want a housekeeping null, or to
        expression-matched genes if you want an expression-matched null.
    n_control : int, optional
        Size of each random draw.  Defaults to ``len(gene_index)``, i.e. draws
        the same number of genes as the pool.  Pass the target set's size to get
        a size-matched null.
    n_repeat : int
        Number of random draws.

    Returns
    -------
    dict with ``h_mean``, ``h_p95``, ``ssr``, ``percentile_of_observed``.
    """
    pool = np.asarray(gene_index, dtype=int).ravel()
    if pool.size == 0:
        return {"h_mean": np.nan, "h_p95": np.nan, "ssr": np.nan,
                "percentile_of_observed": np.nan}

    k = int(n_control if n_control is not None else pool.size)
    k = min(k, pool.size)
    rng = np.random.default_rng(seed)

    hs = []
    for _ in range(n_repeat):
        pick = rng.choice(pool, k, replace=False)
        hs.append(kruskal_h(score_add_module(matrix, pick), labels))
    hs = np.asarray(hs, dtype=float)
    ok = hs[np.isfinite(hs)]

    h_mean = float(ok.mean()) if ok.size else np.nan
    return {
        "n_control_genes": k,
        "h_mean": h_mean,
        "h_p95": float(np.percentile(ok, 95)) if ok.size else np.nan,
        "ssr": signal_specificity_ratio(h_observed, h_mean),
        "percentile_of_observed": (
            float((ok < h_observed).mean() * 100) if ok.size else np.nan
        ),
    }


@dataclass
class AuditResult:
    """Container for one gene set's audit outcome."""

    name: str
    n_genes: int
    n_cells: int

    # T1: cell-type specificity vs size-matched housekeeping
    h_target: float = np.nan
    h_hk40: float = np.nan
    ssr_original: float = np.nan
    h_target_down40: float = np.nan
    ssr_shrink: float = np.nan
    h_hk_upsized: float = np.nan
    ssr_grow: float = np.nan

    # T3: depth coupling per scorer
    depth_coupling: dict = field(default_factory=dict)

    # within-cell-type check
    within_type: list = field(default_factory=list)

    # null distribution
    null_hk40: dict = field(default_factory=dict)

    verdict: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        lines = [
            f"{self.name}: {self.n_genes} genes, {self.n_cells} cells",
            f"  H(target)={self.h_target:.1f}   H(HK40)={self.h_hk40:.1f}",
            f"  SSR (original)     = {self.ssr_original:.3f}",
            f"  SSR (shrink target)= {self.ssr_shrink:.3f}",
            f"  SSR (grow control) = {self.ssr_grow:.3f}",
            "  depth coupling (Spearman vs genes detected):",
        ]
        for k, v in self.depth_coupling.items():
            lines.append(f"    {k:18s} {v:+.3f}")
        lines.append(f"  verdict: {self.verdict}")
        return "\n".join(lines)


def _verdict(ssr: float) -> str:
    if not np.isfinite(ssr):
        return "undetermined"
    if ssr >= 1.0:
        return "passes: target set separates cell types at least as well as housekeeping"
    if ssr >= 0.8:
        return "marginal: barely worse than housekeeping; interpret with caution"
    return "fails: target set is less cell-type specific than housekeeping"


def audit_gene_set(matrix, genes: Sequence[str], labels: Sequence[str],
                   gene_set: Sequence[str],
                   housekeeping: Sequence[str] = HOUSEKEEPING_SET,
                   housekeeping_pool: Sequence[str] = HOUSEKEEPING_POOL,
                   scorers: dict | None = None,
                   n_repeat: int = 200, seed: int = 42,
                   within_type: bool = True) -> AuditResult:
    """Run the full depth-bias audit on one gene set.

    Parameters
    ----------
    matrix : ndarray, shape (n_cells, n_genes)
        **Log-normalised** expression (see :func:`lactaudit.log_normalize`).
        The Pearson-residual scorer will convert this back to counts internally.
    genes : sequence of str
        Column names of ``matrix``, uppercase.
    labels : sequence of str
        Per-cell cell-type labels, length ``n_cells``.
    gene_set : sequence of str
        The set under test.
    housekeeping / housekeeping_pool : sequence of str
        Control panels.  The 40-gene panel is the original comparator; the
        139-gene pool enables the "grow the control" arm.
    scorers : dict, optional
        Subset of :data:`lactaudit.scoring.SCORERS` to evaluate.  Defaults to
        all four.
    within_type : bool
        Also run the within-cell-type depth check, which splits each cell type
        at its median genes-detected and compares scores across halves.  This
        is the cheapest way to see whether a "cell-type signal" survives once
        depth is held roughly constant.

    Returns
    -------
    AuditResult
    """
    genes_up = np.asarray([str(g).upper() for g in genes])
    index = {g: i for i, g in enumerate(genes_up)}
    labels = np.asarray(labels).astype(str)
    matrix = np.asarray(matrix)
    rng = np.random.default_rng(seed)

    tgt = [index[g] for g in gene_set if g in index]
    hk = [index[g] for g in housekeeping if g in index]
    pool = sorted({index[g] for g in housekeeping_pool if g in index} - set(tgt))

    if len(tgt) < 20:
        raise ValueError(
            f"only {len(tgt)} of {len(gene_set)} gene-set members found in the matrix; "
            "check that gene symbols are uppercase and match the species"
        )
    if len(hk) < 10:
        raise ValueError(f"only {len(hk)} housekeeping genes found; the control panel is unusable")

    res = AuditResult(name="custom", n_genes=len(tgt), n_cells=int(matrix.shape[0]))

    # ---- T1: specificity vs housekeeping, three panel-size regimes ----
    score_t = score_add_module(matrix, tgt)
    score_h = score_add_module(matrix, hk)

    res.h_target = kruskal_h(score_t, labels)
    res.h_hk40 = kruskal_h(score_h, labels)
    res.ssr_original = signal_specificity_ratio(res.h_target, res.h_hk40)

    if len(tgt) > len(hk):
        hs = [kruskal_h(score_add_module(matrix, rng.choice(tgt, len(hk), replace=False)), labels)
              for _ in range(min(n_repeat, 100))]
        res.h_target_down40 = float(np.nanmedian(hs))
        res.ssr_shrink = signal_specificity_ratio(res.h_target_down40, res.h_hk40)

    if len(pool) >= len(tgt):
        hs = [kruskal_h(score_add_module(matrix, rng.choice(pool, len(tgt), replace=False)), labels)
              for _ in range(min(n_repeat, 100))]
        res.h_hk_upsized = float(np.nanmedian(hs))
        res.ssr_grow = signal_specificity_ratio(res.h_target, res.h_hk_upsized)

    res.null_hk40 = size_matched_null(
        matrix, pool, labels, res.h_target,
        n_control=len(tgt), n_repeat=n_repeat, seed=seed,
    )

    # ---- T3: depth coupling per scorer ----
    n_detected = (matrix > 0).sum(axis=1)
    for name, fn in (scorers or SCORERS).items():
        try:
            s = fn(matrix, tgt)
            res.depth_coupling[name] = spearman_depth(s, n_detected)
        except Exception as exc:                       # noqa: BLE001
            res.depth_coupling[name] = np.nan
            res.depth_coupling[f"{name}_error"] = str(exc)[:100]

    # ---- within-cell-type check ----
    if within_type:
        for ct in pd.unique(labels):
            m = labels == ct
            if m.sum() < 100:
                continue
            d = n_detected[m]
            med = np.median(d)
            hi = np.where(m)[0][d > med]
            lo = np.where(m)[0][d <= med]
            if len(hi) < 30 or len(lo) < 30:
                continue
            res.within_type.append({
                "cell_type": str(ct),
                "n": int(m.sum()),
                "d_score": cohen_d(score_t[hi], score_t[lo]),
                "d_depth": cohen_d(n_detected[hi], n_detected[lo]),
            })

    res.verdict = _verdict(res.ssr_original)
    return res


def cluster_mat(matrix, k: int = 10, n_hvg: int = 500, seed: int = 0):
    """Cheap unsupervised cell labels via k-means on the most variable genes.

    Provided so an audit can run end-to-end without any annotation.  If you have
    real cell-type labels from the paper, use those instead -- the audit is only
    as meaningful as the groups you compare.

    Parameters
    ----------
    matrix : ndarray, (n_cells, n_genes)
        Log-normalised expression.
    k : int
        Number of clusters.  Usually set to the number of cell types expected.
    n_hvg : int
        Number of most variable genes used to define the clusters.

    Returns
    -------
    np.ndarray of str, length n_cells
    """
    from scipy.cluster.vq import kmeans2

    X = np.asarray(matrix)
    var = X.var(axis=0)
    idx = np.argsort(var)[-min(n_hvg, X.shape[1]):]
    sub = X[:, idx]
    sub = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-9)
    _, lab = kmeans2(sub, k, minit="++", seed=seed)
    return np.array([f"c{i}" for i in lab])
