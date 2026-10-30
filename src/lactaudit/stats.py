"""
Statistics used by the depth-bias audit.
"""

from __future__ import annotations

import numpy as np
from scipy import stats as sps

__all__ = [
    "kruskal_h",
    "signal_specificity_ratio",
    "cohen_d",
    "max_gap",
    "spearman_depth",
]


def _groups_by_label(score, labels, min_cells: int = 5):
    score = np.asarray(score, dtype=float).ravel()
    labels = np.asarray(labels).astype(str).ravel()
    uniq, counts = np.unique(labels, return_counts=True)
    keep = uniq[counts >= min_cells]
    return [score[labels == u] for u in keep]


def kruskal_h(score, labels, min_cells: int = 5) -> float:
    """Kruskal-Wallis H statistic of ``score`` across cell types.

    H answers the question this audit actually cares about: *does this score
    separate cell types?*  A gene set that is truly cell-type specific should
    give a large H; a panel of housekeeping genes should give a small one.

    Cell types with fewer than ``min_cells`` cells are dropped, because H is
    unreliable for tiny groups.

    Returns ``nan`` when fewer than two cell types remain.
    """
    groups = _groups_by_label(score, labels, min_cells)
    if len(groups) < 2:
        return float("nan")
    try:
        return float(sps.kruskal(*groups).statistic)
    except ValueError:
        return float("nan")


def signal_specificity_ratio(h_target: float, h_control: float) -> float:
    """Signal Specificity Ratio (SSR).

        SSR = H(target gene set) / H(size-matched housekeeping)

    Interpretation
    --------------
    SSR > 1  the target set separates cell types better than housekeeping does
             -> it carries genuine cell-type signal
    SSR = 1  indistinguishable from housekeeping
             -> the apparent signal is fully explained by depth
    SSR < 1  worse than housekeeping
             -> the set is *not* specific; using it to rank cell types is invalid

    A value of 1 is the decision boundary, not a target.
    """
    if not np.isfinite(h_target) or not np.isfinite(h_control) or h_control <= 0:
        return float("nan")
    return float(h_target / h_control)


def cohen_d(a, b) -> float:
    """Standardised mean difference between two groups (pooled SD)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 2 or b.size < 2:
        return float("nan")
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2.0)
    if pooled <= 0:
        return float("nan")
    return float((a.mean() - b.mean()) / pooled)


def max_gap(score, labels, min_cells: int = 5) -> float:
    """Difference between the highest and lowest cell-type mean score.

    A more interpretable effect size than H when reporting to readers who are
    not statisticians.
    """
    groups = _groups_by_label(score, labels, min_cells)
    if len(groups) < 2:
        return float("nan")
    means = [g.mean() for g in groups]
    return float(max(means) - min(means))


def spearman_depth(score, n_genes_detected) -> float:
    """Spearman rho between a per-cell score and genes detected per cell.

    This is the T3 statistic of the audit: a direct measure of how much of the
    score is explained by sequencing depth.
    """
    score = np.asarray(score, dtype=float).ravel()
    ng = np.asarray(n_genes_detected, dtype=float).ravel()
    ok = np.isfinite(score) & np.isfinite(ng)
    if ok.sum() < 10:
        return float("nan")
    return float(sps.spearmanr(score[ok], ng[ok]).statistic)
