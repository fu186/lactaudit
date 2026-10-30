# lactaudit

**Depth-bias audit for single-cell gene-set scoring.**

[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-19%20passed-brightgreen)]()

---

## The problem this solves

Single-cell gene-set scores are usually computed as the mean expression of a
gene list. That score is **not** normalised against sequencing depth. A deeply
sequenced cell expresses more genes, so it accumulates more positive z-scores
and therefore scores higher — *regardless of which genes are in the list.*

The practical consequence is that any **broadly expressed** gene set (metabolism,
hypoxia, housekeeping, ribosomal) will appear to "define" whichever cell type
happens to be captured most deeply in each dataset. The cell type that wins
changes from dataset to dataset, and the change looks like biology.

`lactaudit` quantifies this. It compares the cell-type specificity of your gene
set against a **size-matched housekeeping panel** and reports a single ratio:

```
SSR = H(target gene set) / H(size-matched housekeeping)
```

where `H` is the Kruskal–Wallis statistic of the per-cell score across cell
types. SSR > 1 means the set separates cell types better than housekeeping does,
so it carries real signal. SSR ≤ 1 means the apparent signal is fully explained
by depth, and ranking cell types with it is invalid.

---

## Install

```bash
git clone https://github.com/fu186/lactaudit.git
cd lactaudit
pip install -e .
```

> From the source archive supplied with the paper, `pip install .` inside the unpacked
> directory works without any clone. The repository is the same code under version control.

Dependencies are deliberately minimal: `numpy`, `scipy`, `pandas`. No scanpy, no
AnnData, no PyTorch. The Pearson-residual scorer is implemented analytically so
you do not inherit scanpy's version pinning.

Verify the installation in ten seconds without downloading anything:

```bash
lactaudit demo
```

```
Synthetic matrix: 1200 cells x 2000 genes, 4 cell types, library size log-normal(8, 0.8)
  - housekeeping block: expressed equally in all cells (H ~ 0)
  - lactylation set:    broadly expressed, NO cell-type structure
  - planted markers:    15 genes per cell type, boosted 600x

planted markers        SSR = 14.959   passes
LRG (depth-driven)     SSR =  0.460   fails
housekeeping (HK40)    SSR =  1.000   passes
```

The planted markers recover strongly, the broad expression set fails, and the
housekeeping panel sits exactly at 1.0 — as it must, by definition.

---

## Quick start

### Python API

```python
import lactaudit as la
from scipy.sparse import csr_matrix

# 1. Load counts. Any (cells x genes) matrix works.
X, genes, cells = la.read_wide_table("GSE155468_TAA2.txt.gz", max_cells=5000)

# 2. Quality control. IMPORTANT: run this on the FULL transcriptome.
Xk, mask, stats = la.qc_filter(X, genes, max_cells=3000)

# 3. Normalise. Supply the FULL library size, not the subset's row sums.
lib = Xk.sum(axis=1)
mat = la.log_normalize(Xk, libsize=lib)

# 4. Have cell-type labels from anywhere -- clustering, marker scoring, or the paper.
labels = my_annotation_function(mat)

# 5. Audit your gene set.
res = la.audit_gene_set(mat, genes, labels, la.LACTYLATION_SET)
print(res.summary())
```

```
LRG: 109 genes, 3000 cells
  H(target)=159.5   H(HK40)=424.8
  SSR (original)     = 0.376
  SSR (shrink target)= 0.253
  SSR (grow control) = 0.404
  depth coupling (Spearman vs genes detected):
    add_module          +0.588
    gene_zscore         +0.588
    pearson_residual    +0.010
    rank                -0.331
  verdict: fails: target set is less cell-type specific than housekeeping
```

### Command line

```bash
# Audit the built-in lactylation set
lactaudit audit --matrix cache/TAA2/counts.npz --labels labels.tsv

# Audit your own gene list (whitespace-separated file)
lactaudit audit --matrix data.npz --labels labels.tsv --gene-set my_genes.txt --json out.json

# Compare LRG against hypoxia / senescence / stemness / apoptosis
lactaudit compare --matrix data.npz --labels labels.tsv
```

---

## What the three panel-size regimes mean

The obvious objection to comparing a 110-gene set against a 40-gene panel is that
the larger set has more chances to find a difference. `lactaudit` answers this
directly by reporting SSR under three regimes:

| Regime | Target set | Control set | Interpretation |
|---|---|---|---|
| `SSR (original)` | all available | 40-gene HK | the naive comparison |
| `SSR (shrink target)` | down-sampled to 40 | 40-gene HK | equal size, target loses information |
| `SSR (grow control)` | all available | pool expanded to match | equal size, target keeps information |

If all three agree, the conclusion does not depend on panel size. In our own
validation on 16 GEO datasets, shrinking the target **lowered** SSR (0.825 →
0.596) and the direction of the verdict was unchanged in **16/16** datasets. The
size advantage was never the source of the result.

---

## The four scorers

| Name | What it does | Depth coupling |
|---|---|---|
| `add_module` | mean of per-gene z-scores (Seurat `AddModuleScore`) | strong positive |
| `gene_zscore` | mathematically identical to the above | identical — provided as a check |
| `pearson_residual` | analytic Pearson residuals (Lause et al. 2021) | weakened, **not** removed |
| `rank` | within-cell ranks (UCell-style) | sign often flips |

The `gene_zscore` duplication is intentional. Z-scoring is linear, so averaging
per-gene z-scores equals the z-score of the summed expression divided by a
constant. If the two ever disagree, there is a numerical bug.

The `rank` scorer is worth reading about in the source: it commonly yields a
**negative** depth correlation, which people misread as "the depth problem is
solved". It is not solved; the sign has flipped.

Variance stabilisation is the most common proposed fix. It works — we measured a
56% reduction in depth coupling — but in our 16-dataset validation the
correlation remained positive in **16/16** datasets. Coupling is a property of
summing over a gene list, not of the normalisation method.

---

## The three GEO wide-table layouts

If you have ever had a re-analysis silently produce wrong numbers, this is
probably why. GEO series ship matrices in at least three header layouts, and two
of the three are mishandled by the obvious `pd.read_csv(..., index_col=0)` call.

| Layout | First row | Correct read | Consequence of getting it wrong |
|---|---|---|---|
| `no_header` | all barcodes | `header=None, skiprows=1, index_col=0` | pandas raises `Expected N fields, saw N+1` |
| `two_label` | `Gene_ID, Symbol, barcode...` | `index_col=1`, drop column 0 | **silently uses Ensembl IDs as symbols** |
| `one_label` | `Gene, barcode...` | `header=0, index_col=0` | works by luck |

Always check:

```python
print(la.sniff_wide_table("matrix.txt.gz"))   # ('\t', 'no_header', 0)
```

`read_wide_table` dispatches automatically and returns uppercase symbols.

---

## Two traps that will silently corrupt your numbers

Both of these produced wrong results during this package's own development, and
both are easy to miss because the output looks plausible.

### 1. Normalise with the full library size

```python
# WRONG: X here is already a 396-gene subset
mat = la.log_normalize(X_subset)

# RIGHT
mat = la.log_normalize(X_subset, libsize=full_library_size)
```

Using a subset's row sums underestimates every size factor, inflating every
score. In our testing this alone moved a specificity ratio from 0.35 to 2.95 — a
complete reversal of the conclusion. `save_cached` / `load_cached` store
`libsize` for exactly this reason.

### 2. Quality-control on the full transcriptome

```python
# WRONG: a 396-gene subset can never have 200 genes detected
qc_filter(X_subset, genes, min_genes=200)     # discards ~99% of cells

# RIGHT
qc_filter(X_full, all_genes, min_genes=200)
```

Filtering on a gene subset makes the `min_genes` threshold meaningless. In our
testing this reduced one dataset from 4,598 cells to 2 — and analyses of 2 cells
still return numbers.

---

## API reference

### Gene sets

| Object | Size | Description |
|---|---|---|
| `LACTYLATION_SET` | 116 | lactylation-related genes (chromatin writers/erasers, lactate production/transport, reported targets) |
| `HOUSEKEEPING_SET` | 40 | ribosomal, cytoskeletal, universal metabolic genes |
| `HOUSEKEEPING_POOL` | 139 | expanded panel for the "grow the control" arm |
| `CELL_TYPE_MARKERS` | 10 types | canonical markers (5–10 per type) |
| `REFERENCE_SIGNATURES` | 4 | hypoxia, senescence, stemness, apoptosis (negative control) |

### Functions

```python
la.log_normalize(counts, libsize=None, target_sum=1e4)      # -> ndarray
la.score_add_module(matrix, gene_index)                      # -> ndarray
la.score_gene_zscore(matrix, gene_index)                     # -> ndarray
la.score_pearson_residual(matrix, gene_index)                # -> ndarray
la.score_rank(matrix, gene_index)                            # -> ndarray

la.kruskal_h(score, labels, min_cells=5)                     # -> float
la.signal_specificity_ratio(h_target, h_control)             # -> float
la.spearman_depth(score, n_genes_detected)                   # -> float
la.cohen_d(a, b)                                             # -> float
la.max_gap(score, labels)                                    # -> float

la.size_matched_null(matrix, pool_index, labels, h_observed) # -> dict
la.audit_gene_set(matrix, genes, labels, gene_set, ...)      # -> AuditResult

la.read_10x(matrix, genes, barcodes)                         # -> (X, genes, cells)
la.read_wide_table(path, max_cells=None)                     # -> (X, genes, cells)
la.sniff_wide_table(path)                                    # -> (sep, variant, n_label)
la.qc_filter(X, gene_names, min_genes=200, max_cells=3000)   # -> (X, mask, stats)
la.save_cached(path, X, genes, libsize, meta=None)
la.load_cached(path)                                         # -> (X, genes, libsize)
```

---

## Repository layout

```
lactaudit/
├── src/lactaudit/
│   ├── __init__.py       public API
│   ├── genesets.py       gene set definitions + provenance
│   ├── scoring.py        the four scorers
│   ├── stats.py          H, SSR, Cohen's d, Spearman
│   ├── audit.py          the audit pipeline
│   ├── io.py             readers, QC, caching
│   └── cli.py            command line interface
├── tests/
│   └── test_lactaudit.py 19 tests
├── examples/
│   ├── basic_audit.py
│   └── compare_signatures.py
├── pyproject.toml
├── LICENSE
└── README.md
```

---

## Testing

```bash
pip install -e ".[dev]"
pytest -q            # 19 passed
lactaudit demo       # end-to-end sanity check
```

The test suite covers gene-set integrity, the linear-algebra equivalence of the
two z-score scorers, the direction and bounds of the rank scorer, header sniffing
for all three GEO layouts, normalisation edge cases (mismatched library sizes,
all-zero cells), QC down-sampling, cache round-trips, and an end-to-end check
that a planted marker set passes while a broad expression set fails.

---

## When SSR ≤ 1 is the honest answer

An SSR below 1 does not mean your biology is wrong. It means **this scoring
approach cannot attribute the apparent signal to cell identity.** Three responses
are reasonable:

1. Report the set as an **unresolved composite** and drop cell-type rankings from
   the claims (this is what we did).
2. Switch to a method that models depth explicitly — differential expression,
   decoupling, or a generative model — rather than to a different ad-hoc score.
   Changing the scorer only moves the coupling around; we measured this.
3. Report SSR alongside any cell-type claim, so readers can judge.

What does not work: choosing a scorer because it produces the ranking you
expected. That is how the original problem arose.

---

## Citation

If you use this package, please cite the companion paper (details on acceptance)
and the method it implements:

> Lause J, Berens P, Kobak D. Analytic Pearson residuals for normalization of
> single-cell RNA-seq UMI data. *Genome Biology* 2021;22:258.

## License

MIT — see [LICENSE](LICENSE). Use it, fork it, ignore it; just check your SSR.
