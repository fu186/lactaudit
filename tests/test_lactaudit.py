"""
Tests for lactaudit.

Run with:  pytest -q
"""

import numpy as np
import pytest

import lactaudit as la


@pytest.fixture(scope="module")
def synthetic():
    """Matrix where depth is the only real structure, plus 30 planted markers."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 600, 800
    labels = rng.choice(["a", "b", "c"], n_cells)
    depth = rng.lognormal(8.0, 0.8, n_cells)

    base = rng.gamma(0.6, 1.0, n_genes)
    counts = rng.poisson(np.outer(depth, base / base.sum()) * 1e3)

    marker_idx = rng.choice(n_genes, 30, replace=False)
    for i, lab in enumerate(labels):
        if lab == "a":
            counts[i, marker_idx] += 60

    genes = np.array([f"G{j}" for j in range(n_genes)], dtype=object)
    # Assign LRG to the *high-expression* genes and HK to mid-expression genes,
    # mirroring reality: lactylation genes are broadly expressed, housekeeping
    # genes are abundant but far fewer in number.
    order = np.argsort(-base)
    lrg_slots = order[30:30 + len(la.LACTYLATION_SET)]
    hk_slots = order[200:200 + len(la.HOUSEKEEPING_SET)]
    pool_slots = order[200:200 + len(la.HOUSEKEEPING_POOL)]
    for j, g in zip(lrg_slots, la.LACTYLATION_SET):
        genes[j] = g
    for j, g in zip(hk_slots, la.HOUSEKEEPING_SET):
        genes[j] = g
    for j, g in zip(pool_slots, la.HOUSEKEEPING_POOL):
        genes[j] = g

    mat = la.log_normalize(counts)
    return mat, genes, labels, [genes[j] for j in marker_idx]


def test_gene_sets_are_uppercase_and_unique():
    for name in ("LACTYLATION_SET", "HOUSEKEEPING_SET", "HOUSEKEEPING_POOL"):
        gs = getattr(la, name)
        assert len(gs) == len(set(gs)), f"{name} has duplicates"
        assert all(g == g.upper() for g in gs), f"{name} has non-uppercase symbols"


def test_housekeeping_is_subset_of_pool():
    assert set(la.HOUSEKEEPING_SET) <= set(la.HOUSEKEEPING_POOL)


def test_lactylation_set_size():
    assert len(la.LACTYLATION_SET) >= 100


def test_log_normalize_uses_supplied_libsize():
    """A wrong library size must change the result -- guards against the bug
    where a gene subset is normalised by its own sum."""
    counts = np.array([[10, 20, 30], [1, 2, 3]], dtype=np.float32)
    a = la.log_normalize(counts, libsize=np.array([60, 6]))
    b = la.log_normalize(counts, libsize=np.array([600, 60]))
    assert not np.allclose(a, b)


def test_log_normalize_rejects_mismatched_libsize():
    counts = np.zeros((3, 4), dtype=np.float32)
    with pytest.raises(ValueError):
        la.log_normalize(counts, libsize=np.array([1, 2]))


def test_add_module_and_gene_zscore_agree(synthetic):
    mat, genes, labels, _ = synthetic
    idx = {g: i for i, g in enumerate(genes)}
    gi = [idx[g] for g in la.LACTYLATION_SET if g in idx]
    a = la.score_add_module(mat, gi)
    b = la.score_gene_zscore(mat, gi)
    assert np.allclose(a, b, equal_nan=True)


def test_scores_track_depth(synthetic):
    """The core claim: summed z-scores rise with genes detected."""
    mat, genes, labels, _ = synthetic
    idx = {g: i for i, g in enumerate(genes)}
    gi = [idx[g] for g in la.LACTYLATION_SET if g in idx]
    n_det = (mat > 0).sum(axis=1)
    rho = la.spearman_depth(la.score_add_module(mat, gi), n_det)
    assert rho > 0.1, f"expected positive depth coupling, got {rho:.3f}"


def test_rank_score_orientation():
    """Rank score must rise with expression, not fall.

    This is a unit test of the implementation, not of any biological claim: a
    cell whose gene-set members sit at the top of its own expression ranking
    must score higher than one where they sit at the bottom.
    """
    mat = np.array([
        [9.0, 5.0, 1.0, 0.0],   # set members (cols 0,1) are highest
        [0.0, 1.0, 5.0, 9.0],   # set members are lowest
        [2.0, 2.0, 2.0, 2.0],   # tie
    ])
    s = la.score_rank(mat, [0, 1])
    assert s[0] > s[1], "cell with high set expression must score higher"
    assert 0.0 <= s.min() <= s.max() <= 1.0


def test_rank_score_has_no_negative_values():
    """Rank scores are bounded in [0, 1] by construction."""
    rng = np.random.default_rng(0)
    mat = rng.random((100, 50)).astype(np.float32)
    s = la.score_rank(mat, list(range(10)))
    assert s.min() >= 0.0 and s.max() <= 1.0


def test_pearson_residual_weakens_but_does_not_remove(synthetic):
    mat, genes, labels, _ = synthetic
    idx = {g: i for i, g in enumerate(genes)}
    gi = [idx[g] for g in la.LACTYLATION_SET if g in idx]
    n_det = (mat > 0).sum(axis=1)
    base = la.spearman_depth(la.score_add_module(mat, gi), n_det)
    pr = la.spearman_depth(la.score_pearson_residual(mat, gi), n_det)
    assert pr < base, "variance stabilisation should reduce coupling"


def test_planted_markers_pass_and_noise_fails(synthetic):
    mat, genes, labels, markers = synthetic
    good = la.audit_gene_set(mat, genes, labels, markers, n_repeat=30)
    bad = la.audit_gene_set(mat, genes, labels, la.LACTYLATION_SET, n_repeat=30)
    assert good.ssr_original > bad.ssr_original
    assert "passes" in good.verdict
    assert good.ssr_original > 2.0 * bad.ssr_original


def test_kruskal_h_needs_two_groups():
    score = np.arange(20, dtype=float)
    labels = np.array(["x"] * 20)
    assert np.isnan(la.kruskal_h(score, labels))


def test_ssr_boundary():
    assert la.signal_specificity_ratio(100, 100) == pytest.approx(1.0)
    assert la.signal_specificity_ratio(50, 100) == pytest.approx(0.5)
    assert np.isnan(la.signal_specificity_ratio(100, 0))


def test_sniff_wide_table(tmp_path):
    """All three GEO header variants must be identified correctly."""
    import gzip

    # (a) no_header: first row is all barcodes
    p1 = tmp_path / "a.txt.gz"
    with gzip.open(p1, "wt") as f:
        f.write("AAACCTGAGTAGATGT.1\tAAACCTGAGTTAGCGG.1\n")
        f.write("ACTB\t5\t3\n")
        f.write("GAPDH\t2\t7\n")
    sep, variant, nlabel = la.sniff_wide_table(p1)
    assert sep == "\t" and variant == "no_header" and nlabel == 0

    # (b) two_label: Gene_ID, Symbol, barcodes
    p2 = tmp_path / "b.csv.gz"
    with gzip.open(p2, "wt") as f:
        f.write("Gene_ID,Symbol,AAACCTGCAATGAATG.1,AAACCTGGTGAAAGAG.1\n")
        f.write("ENSG00000000003,TSPAN6,0,1\n")
    sep, variant, nlabel = la.sniff_wide_table(p2)
    assert sep == "," and variant == "two_label" and nlabel == 2

    # (c) one_label: Gene, barcodes
    p3 = tmp_path / "c.tsv.gz"
    with gzip.open(p3, "wt") as f:
        f.write("Gene\tAAACCTGCAATGAATG.1\tAAACCTGGTGAAAGAG.1\n")
        f.write("ACTB\t0\t1\n")
    sep, variant, nlabel = la.sniff_wide_table(p3)
    assert variant == "one_label" and nlabel == 1


def test_read_wide_table_no_header(tmp_path):
    import gzip

    p = tmp_path / "x.txt.gz"
    with gzip.open(p, "wt") as f:
        f.write("AAACCTGAGTAGATGT.1\tAAACCTGAGTTAGCGG.1\n")
        f.write("ACTB\t5\t3\n")
        f.write("GAPDH\t2\t7\n")
    X, genes, cells = la.read_wide_table(p)
    assert X.shape == (2, 2)
    assert list(genes) == ["ACTB", "GAPDH"]
    assert list(cells) == ["AAACCTGAGTAGATGT.1", "AAACCTGAGTTAGCGG.1"]
    assert X.toarray()[0].tolist() == [5, 2]      # cell 1: ACTB=5, GAPDH=2


def test_read_wide_table_two_label_uses_symbol(tmp_path):
    import gzip

    p = tmp_path / "y.csv.gz"
    with gzip.open(p, "wt") as f:
        f.write("Gene_ID,Symbol,AAACCTGCAATGAATG.1,AAACCTGGTGAAAGAG.1\n")
        f.write("ENSG00000000003,TSPAN6,0,1\n")
        f.write("ENSG00000000005,TNMD,4,0\n")
    X, genes, cells = la.read_wide_table(p)
    assert list(genes) == ["TSPAN6", "TNMD"], "must use Symbol, not Gene_ID"
    assert X.shape == (2, 2)


def test_qc_filter_downsamples_without_bias():
    from scipy.sparse import csr_matrix
    rng = np.random.default_rng(0)
    counts = rng.poisson(3.0, (500, 300)).astype(np.float32)
    X = csr_matrix(counts)
    genes = np.array([f"G{i}" for i in range(300)])
    Xk, mask, stats = la.qc_filter(X, genes, max_cells=100)
    assert Xk.shape[0] == 100
    assert stats["n_used"] == 100


def test_save_and_load_cache(tmp_path):
    from scipy.sparse import csr_matrix

    counts = csr_matrix(np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32))
    genes = ["A", "B", "C"]
    lib = np.array([6.0, 15.0])
    la.io.save_cached(tmp_path / "m.npz", counts, genes, lib, {"note": "test"})
    X, g, l = la.load_cached(tmp_path / "m.npz")
    assert X.shape == counts.shape
    # genes come back as an ndarray, so compare element-wise rather than with ==,
    # which would return an array and make the assertion ambiguous.
    assert np.array_equal(np.asarray(g), np.asarray(genes))
    assert np.allclose(l, lib)


def test_audit_gene_set_rejects_tiny_overlap():
    rng = np.random.default_rng(0)
    mat = rng.random((100, 50)).astype(np.float32)
    genes = np.array([f"G{i}" for i in range(50)])
    labels = rng.choice(["a", "b"], 100)
    with pytest.raises(ValueError):
        la.audit_gene_set(mat, genes, labels, la.LACTYLATION_SET)
