# Reproduction Guide

This directory contains all audit code, the per-dataset result records, and the
three self-check scripts used before submission. **Every number in the paper can
be recomputed from the contents of this directory alone; no external repository
is required.**

---

## Directory layout

| Subdirectory | Contents |
|---|---|
| `01_audit_library/` | `lactcheck.py` — the core audit library. Implements the four control tests (T1–T4), the T4r type-level partialling, the SSR metric and the verdict logic. Also holds the marker panels and the seven published lactylation gene lists under audit. |
| `02_dataset_readers/` | Readers for four input formats, plus aggregation and normalisation scripts. |
| `03_recompute/` | Three-estimator comparison (13 datasets x 3 implementations) and the T4r type-level bootstrap. |
| `04_simulation/` | The T4/T4r simulation script and its output JSON (the sole source of the Figure 5 and main-text simulation values). |
| `05_figures/` | Journal-specification figure scripts (300 dpi). |
| `06_verification/` | The three self-check scripts (see below). |
| `07_result_records/` | One audit record per dataset: `*.json` holds every T1–T4 intermediate quantity, `*_T4.csv` the per-cell-type detail. |
| `08_raw_data_download/` | Instructions for obtaining the raw data. |

---

## Environment

Python 3.11 with:

```
numpy  pandas  scipy  scikit-learn  scanpy  matplotlib
```

The original environment pinned exact versions so that the audit reproduces
inside a single container.

---

## Raw data

All raw data are public GEO deposits and are **not distributed with this
package**. The full set of 9 accessions is about 5.8 GB; the 7 accessions needed
for the main analysis are about 1.2 GB. See
`08_raw_data_download/README_data_acquisition.md`.

The audit library and all three verification scripts **do not need the raw
data** — they read `07_result_records/` directly. Only
`03_recompute/run_estimator_extended.py` requires the expression matrices.

---

## Reproduction steps

### A. Without raw data (run this first to confirm package self-consistency)

```bash
cd <this directory>

# 1) Full numeric re-verification.
#    Independently recomputes every Table 1 field, all T4r columns, the BH
#    correction and each aggregate, and scans the text for stale values.
#    Expected: "No findings. All checked numbers reconcile."
python3 06_verification/verify_all_numbers.py

# 2) Cross-manuscript consistency check.
#    21 shared facts must appear in all 3 manuscripts simultaneously; also
#    covers text-to-artifact provenance, figure-script sourcing, Table 1
#    internal consistency, abstract-vs-main-text, and supplementary inventory.
#    Expected: all checks pass.
python3 06_verification/crosscheck_all.py

# 3) Figure-versus-text simulation consistency.
#    Triangulates from the simulation JSON, the plotting script, and the three
#    manuscript versions.
#    Expected: simulation values agree across all three sources.
python3 06_verification/check_fig_text_consistency.py

# 4) Re-run the T4r bootstrap (about 10 s).
#    Reproduces Table S5b byte for byte; expected "median width = 0.920,
#    intervals containing zero = 7/15".
python3 03_recompute/rerun_t4r_bootstrap.py

# 5) Re-run the simulation (about 2 min), producing sim_t4_t4r.json.
python3 04_simulation/sim_t4_t4r.py

# 6) Redraw all figures (requires matplotlib).
python3 05_figures/make_bib_figs.py
```

### B. With raw data

```bash
# Three-estimator comparison (13 datasets x 3 implementations, about 150 s).
# Data must sit in <this directory>/../data/ or be named by AUDIT_DATA.
python3 03_recompute/run_estimator_extended.py

# Audit a single new dataset (10x MTX directory shown).
python3 02_dataset_readers/audit_one.py \
    --indir <mtx_dir> --name MyDataset --species human --outdir ./results
# For mouse data add --species mouse to switch to the mouse marker panel.
```

---

## Path resolution

Every script locates the result records and manuscript directories inside the
package automatically, so **no path editing is needed**. Both the workspace
layout and the submission layout are supported.

If the directory has been moved such that resolution fails, set:

```bash
export AUDIT_ROOT=<directory containing 07_result_records/>
export AUDIT_DATA=<raw data directory>   # only 03_recompute/run_estimator_extended.py needs this
```

---

## Key aspects of the audit

### The four control tests

| Code | What it does | Pass condition |
|---|---|---|
| **T1** | Housekeeping substitution: SSR = H(LRG) / mean H(40 housekeeping genes) | SSR > 1 |
| **T2** | Expression-matched permutation null: 150 equally sized random sets drawn from the 5th–95th percentile expression band of the LRG set itself | p < 0.05 |
| **T3** | Within-cell residualisation against technical covariates. **Genes detected is the pre-specified primary covariate** | \|rho\| > 0.3 flagged as strong coupling (diagnostic only) |
| **T4** | Cross-cell-type collinearity with detection rate | **Not used as a test** — the simulation reaches +0.991 even with no biological signal |
| **T4r** | Type-level partialling: both variables are residualised on per-type mean genes detected by OLS, then correlated | Interpretable only when **span >= 0.15 and n_types >= 6** |

### Positive control

The same battery run against lineage marker panels **known to be cell-type
specific**: SSR 1.32–62.1 (median 4.38), T2 p = 0.008. The positive control
**passes**, showing the tests can detect real signal.

### Conventions that matter

- **T2 uses 150 permutations**, so p values fall on a **1/151 grid** with minimum
  non-zero value **0.0066**.
- **The BH correction is computed over 16 p values**, giving minimum q =
  **0.106**. (Computing it over 15 gives 0.099, which is wrong.)
- **The pre-scored row `AS_GSE159677_ours` is not comparable**: it was produced by
  a different pipeline and the text excludes it from **all aggregate statistics**,
  so the comparable corpus is **15 datasets**. The `comparable()` predicate in
  `crosscheck_all.py` implements that exclusion.
- **The T4r simulation null band sits below zero** (−0.750 to +0.141), so the
  observed +0.770 **cannot** be read as biological evidence. The paper states this
  limitation; the null band is not a usable reference.

---

## Coverage of the three verification scripts

| Script | Coverage |
|---|---|
| `verify_all_numbers.py` | 9 check sections: field-by-field Table 1 recomputation, independent T4r recomputation, BH correction, bootstrap interval consistency, and the round-four additions (T2 matching quality, real depth span, the Poisson claim, selection-bias wording) |
| `crosscheck_all.py` | 6 check sections: (A) shared facts across the three manuscripts, (B) text-to-artifact provenance, (C) figure-script sourcing, (D) Table 1 internal consistency, (E) abstract vs main text, (F) supplementary inventory vs actual files |
| `check_fig_text_consistency.py` | 3-source triangulation: simulation JSON, plotting script (which must not hard-code simulation values), and the simulation context in each manuscript |

Any non-zero exit code from these scripts means an inconsistency was found.

---

## Verification record

After packaging, the scripts were re-run from a freshly unpacked copy:

| Check | Result |
|---|---|
| `verify_all_numbers.py` | No findings. All checked numbers reconcile. |
| `crosscheck_all.py` | All checks pass: 3 manuscripts consistent, values traceable, figures sourced, table self-consistent. |
| `check_fig_text_consistency.py` | Figure and text simulation values agree; plotting script is sourced from the JSON. |
| `rerun_t4r_bootstrap.py` | Median width = 0.920 (range 0.340–1.565), intervals containing zero = 7/15 — **byte-identical** to the packaged `t4r_bootstrap_corrected.csv` |
