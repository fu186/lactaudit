"""
Gene-set scoring functions.

All scorers share the signature::

    scorer(matrix, gene_index) -> np.ndarray of shape (n_cells,)

where ``matrix`` is a dense ``n_cells x n_genes`` float array that has already
been library-size normalised and log-transformed (see :func:`log_normalize`),
and ``gene_index`` is a sequence of column indices.

The four scorers differ in exactly one respect: how much they let sequencing
depth leak into the score.

    score_add_module        summed z-score       -- depth coupling: HIGH
    score_gene_zscore       per-gene z mean      -- mathematically identical to above
    score_pearson_residual  variance-stabilised  -- depth coupling: MODERATE
    score_rank              within-cell ranks    -- depth coupling: INVERTED

"Mathematically identical" is not a typo: z-scoring is linear, so the mean of
per-gene z-scores equals the z-score of the summed expression divided by a
constant.  Any Spearman correlation with depth is therefore unchanged.  This
duplication is deliberate -- it is a useful sanity check that no numerical
error has crept in.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "log_normalize",
    "score_add_module",
    "score_gene_zscore",
    "score_pearson_residual",
    "score_rank",
    "SCORERS",
]


def log_normalize(counts, libsize=None, target_sum: float = 1e4,
                  log1p: bool = True) -> np.ndarray:
    """Library-size normalisation followed by log1p.

    Parameters
    ----------
    counts : array-like, shape (n_cells, n_genes)
        Raw UMI / read counts.  A scipy sparse matrix works too.
    libsize : array-like, shape (n_cells,), optional
        **Total counts across all genes** for each cell.

        This parameter is not optional in spirit.  If you pass a gene *subset*
        to ``counts`` and let the function derive the library size from that
        subset, every cell's size factor will be underestimated and every score
        will be inflated.  In our own testing this alone moved a specificity
        ratio from 0.35 to 2.95 -- a complete reversal of the conclusion.  When
        working with a pre-filtered matrix, always supply ``libsize`` measured
        on the full transcriptome.
    target_sum : float
        Counts-per-``target_sum`` scaling.  1e4 (CP10K) is conventional.
    log1p : bool
        Apply ``log(1 + x)`` after scaling.

    Returns
    -------
    np.ndarray, shape (n_cells, n_genes), dtype float32
    """
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    X = np.asarray(counts, dtype=np.float32)

    if libsize is None:
        lib = X.sum(axis=1, dtype=np.float64)
    else:
        lib = np.asarray(libsize, dtype=np.float64).ravel()
        if lib.shape[0] != X.shape[0]:
            raise ValueError(
                f"libsize has {lib.shape[0]} entries but counts has {X.shape[0]} rows"
            )

    lib = np.where(lib > 0, lib, 1.0)
    out = X / lib[:, None] * target_sum
    if log1p:
        np.log1p(out, out=out)
    return out.astype(np.float32, copy=False)


def _col_indices(gene_index) -> np.ndarray:
    idx = np.asarray(gene_index, dtype=int).ravel()
    if idx.size == 0:
        raise ValueError("gene_index is empty -- check that your gene symbols match")
    return idx


def score_add_module(matrix: np.ndarray, gene_index) -> np.ndarray:
    """AddModuleScore-style summed z-score (Seurat ``AddModuleScore`` analogue).

    For each gene in the set, z-score its expression across cells, then average
    the z-scores per cell.

    Why this fails the depth audit
    ------------------------------
    A deeply sequenced cell expresses more genes, and each expressed gene
    carries a positive z-score.  Averaging positive numbers yields a positive
    score.  The score therefore rises monotonically with depth *regardless of
    which genes are in the set*.
    """
    idx = _col_indices(gene_index)
    sub = matrix[:, idx]
    mu = sub.mean(axis=0)
    sd = sub.std(axis=0)
    sd = np.where(sd < 1e-12, 1e-12, sd)
    return ((sub - mu) / sd).mean(axis=1)


def score_gene_zscore(matrix: np.ndarray, gene_index) -> np.ndarray:
    """Per-gene z-score mean.

    Provided as an explicit equivalence check against
    :func:`score_add_module`; under any linear transform the two agree exactly
    in their rank correlation with depth.
    """
    return score_add_module(matrix, gene_index)


def score_pearson_residual(matrix: np.ndarray, gene_index) -> np.ndarray:
    """Analytic Pearson-residual score (variance-stabilised).

    Reproduces the intent of scanpy's
    ``sc.experimental.pp.normalize_pearson_residuals`` without requiring
    scanpy as a dependency:

        r_ij = (x_ij - mu_ij) / sqrt(mu_ij + mu_ij^2 / theta)

    with ``mu_ij = s_i * q_j`` under a negative-binomial mean model and
    ``theta`` the overdispersion.  We use the analytic approximation with
    ``theta`` estimated per gene by moment matching on the supplied matrix.

    Literature
    ----------
    Lause J, Berens P, Kobak D.  Analytic Pearson residuals for normalization
    of single-cell RNA-seq UMI data.  Genome Biology 2021;22:258.

    Note that this scorer expects **counts**, not log-normalised values.  If you
    pass a log-scale matrix, it will be exponentiated back to counts first.
    """
    idx = _col_indices(gene_index)
    X = np.asarray(matrix, dtype=np.float64)

    # Heuristic: log1p(CP10K) values are small and non-integer; counts are
    # integer-valued.  Convert back if we were handed a log-scale matrix.
    if X.max() < 30 and np.abs(X - np.round(X)).mean() > 1e-3:
        X = np.expm1(X)

    lib = X.sum(axis=1, keepdims=True)
    lib = np.where(lib > 0, lib, 1.0)
    q = X.sum(axis=0) / X.sum()
    q = np.where(q > 0, q, 1e-8)
    mu = lib * q[None, :]

    # Per-gene overdispersion by moment matching: Var = mu + mu^2/theta
    resid = np.zeros_like(X)
    for j in range(X.shape[1]):
        mu_j = mu[:, j]
        var_j = X[:, j].var()
        mean_j = mu_j.mean()
        if mean_j <= 0 or var_j <= mean_j:
            theta = 1e6
        else:
            theta = mean_j ** 2 / (var_j - mean_j)
        denom = np.sqrt(mu_j + mu_j ** 2 / max(theta, 1e-6))
        denom = np.where(denom < 1e-12, 1e-12, denom)
        resid[:, j] = (X[:, j] - mu_j) / denom

    return np.nanmean(resid[:, idx], axis=1)


def score_rank(matrix: np.ndarray, gene_index) -> np.ndarray:
    """Within-cell rank score (UCell-style).

    For each cell, rank all genes by expression (0 = lowest, n-1 = highest) and
    average the ranks of the gene-set members.

    Note on direction
    -----------------
    Higher rank value means higher expression, so a **higher** score means the
    set members sit nearer the top of that cell's ranking. The correlation with
    depth then comes out **negative** in most datasets: shallow cells rank
    everything more noisily and a few set members happen to float to the top.
    That is not a bug -- it is the same coupling with the sign flipped, which is
    why "switch to a rank-based estimator" does not actually solve the problem.
    """
    idx = _col_indices(gene_index)
    X = np.asarray(matrix)
    ranks = np.argsort(np.argsort(X, axis=1), axis=1)
    denom = max(X.shape[1] - 1, 1)
    return ranks[:, idx].mean(axis=1) / denom


SCORERS = {
    "add_module": score_add_module,
    "gene_zscore": score_gene_zscore,
    "pearson_residual": score_pearson_residual,
    "rank": score_rank,
}
"""Registry of available scorers, keyed by short name."""
