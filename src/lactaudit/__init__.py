"""
lactaudit — Depth-bias audit toolkit for single-cell gene-set scoring.

Public API
----------
Gene sets
    LACTYLATION_SET          the 116-gene lactylation-related set (LRG)
    HOUSEKEEPING_SET         the 40-gene housekeeping panel (HK)
    HOUSEKEEPING_POOL        a 139-gene expanded housekeeping pool
    CELL_TYPE_MARKERS        marker dict for ten major cell types
    REFERENCE_SIGNATURES     non-metabolic reference signatures

Scoring
    score_add_module         AddModuleScore-style summed z-score
    score_gene_zscore        per-gene z-score mean (mathematically equivalent)
    score_pearson_residual   analytic Pearson residual (variance-stabilised)
    score_rank               rank-based score (UCell-style)

Statistics
    kruskal_h                Kruskal-Wallis H across cell types
    signal_specificity_ratio SSR = H(target) / H(size-matched housekeeping)

Audit
    size_matched_null        build the size-matched housekeeping null
    audit_gene_set           full T1/T3 audit of one expression matrix

I/O
    read_10x                 read a 10x mtx directory triplet
    read_wide_table          read a GEO-style wide table (three header variants)
    load_cached              load a matrix cached by ``lactaudit cache``
"""

from .genesets import (
    LACTYLATION_SET,
    HOUSEKEEPING_SET,
    HOUSEKEEPING_POOL,
    CELL_TYPE_MARKERS,
    REFERENCE_SIGNATURES,
)
from .scoring import (
    log_normalize,
    score_add_module,
    score_gene_zscore,
    score_pearson_residual,
    score_rank,
    SCORERS,
)
from .stats import (
    kruskal_h,
    signal_specificity_ratio,
    cohen_d,
    max_gap,
    spearman_depth,
)
from .audit import (
    size_matched_null,
    audit_gene_set,
    AuditResult,
    cluster_mat,
)
from .io import (
    read_10x,
    read_wide_table,
    sniff_wide_table,
    qc_filter,
    save_cached,
    load_cached,
)

__version__ = "1.0.0"

__all__ = [
    # gene sets
    "LACTYLATION_SET", "HOUSEKEEPING_SET", "HOUSEKEEPING_POOL",
    "CELL_TYPE_MARKERS", "REFERENCE_SIGNATURES",
    # scoring
    "log_normalize", "score_add_module", "score_gene_zscore",
    "score_pearson_residual", "score_rank", "SCORERS",
    # stats
    "kruskal_h", "signal_specificity_ratio", "cohen_d", "max_gap",
    "spearman_depth",
    # audit
    "size_matched_null", "audit_gene_set", "AuditResult", "cluster_mat",
    # io
    "read_10x", "read_wide_table", "sniff_wide_table", "qc_filter",
    "save_cached", "load_cached",
]
