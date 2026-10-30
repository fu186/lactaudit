# Raw Data Acquisition

The 16 datasets used in this audit all come from **NCBI Gene Expression Omnibus
(GEO)** across **9 accessions**. No new data were generated.

Total download size is approximately **5.8 GB for all 9 accessions**, or
**1.2 GB for the 7 accessions needed to reproduce the main analysis** (the two
remaining accessions are described under [Skipping GSE159677 and
GSE214611](#skipping-gse159677-and-gse214611)). Because of the size, the raw
matrices are **not distributed with the submission package**; reviewers can
download them directly from GEO.

---

## Accessions

| Accession | Disease | Tissue / species | Label in this audit | Cells |
|---|---|---|---|---|
| GSE166676 | Abdominal aortic aneurysm (AAA) | Human aorta | AAA | 3,000 |
| GSE159677 | Atherosclerosis (AS) | Human carotid | AS (pre-scored row, not comparable) | 3,000 |
| GSE214611 | Atherosclerosis (AS) | Mouse aorta | AS_mouse | 3,500 |
| GSE131778 | Coronary artery disease (CAD) | Human coronary artery | CAD | 3,500 |
| GSE165816 | Diabetic foot ulcer (DFU) | Human skin | DFU | 4,173 |
| GSE212966 | Pancreatic ductal adenocarcinoma (PDAC) | Human pancreas | PDAC | 2,999 |
| GSE155468 | Thoracic aortic aneurysm (TAA) | Human aorta | TAA1-TAA8 (8 patients) | 8 x 2,999 |
| GSE156632 | Clear cell renal cell carcinoma (ccRCC) | Human kidney | ccRCC | 4,749 |
| GSE171306 | Clear cell renal cell carcinoma (ccRCC) | Human kidney | ccRCC2 | 3,000 |

> Exact cell counts, tissue provenance and genes detected for each dataset are
> recorded in the JSON files under `../07_result_records/`; the summary table is
> **Table S1** in the supplementary material.

---

## Input formats

Data enter the pipeline in three formats. Readers live in `../02_dataset_readers/`:

| Format | Reader | Datasets in this audit |
|---|---|---|
| 10x Genomics MTX triplet | `audit_one.py` | GSE214611 |
| Wide CSV (`Gene_ID, Symbol, barcode...`) | `audit_csv.py` | GSE156632 |
| Wide TSV (header variant auto-detected) | `audit_tsv.py` | GSE131778 / GSE165816 / GSE155468 |
| Pre-scored cell table | `audit_prescored.py` | GSE159677 |

### The three wide-table layouts

`audit_tsv.py` sniffs the header before parsing. The three variants that occur in
GEO series are:

| Variant | First row | Correct read |
|---|---|---|
| `no_header` | all barcodes | `header=None, skiprows=1, index_col=0` |
| `two_label` | `Gene_ID, Symbol, barcode...` | use the **Symbol** column, drop `Gene_ID` |
| `one_label` | `Gene, barcode...` | `header=0, index_col=0` |

Two of the three are mishandled by a plain `pd.read_csv(..., index_col=0)` call.
The `no_header` case raises `Expected N fields in line 2, saw N+1`; the
`two_label` case fails **silently**, using Ensembl gene IDs as if they were gene
symbols. `audit_tsv.sniff_wide_table(path)` reports which variant a file uses so
the choice can be verified rather than assumed.

---

## Downloading

### Option 1: GEO web interface

Open `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>` and use the
"Supplementary file" links near the bottom of the page.

### Option 2: command line (recommended)

Downloading a GEO series requires **two** URL patterns. The series-level
directory does not contain the per-sample 10x matrices:

```
series  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166676/suppl/
samples https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5077nnn/GSM5077729/suppl/
```

A `filelist.txt` in the series `suppl/` directory names every file; the `File`
and `Archive` columns give the name and size, and files whose names begin with
`GSM` must be fetched from the **sample** directory, not the series directory.

The helper script below handles both patterns, skips the redundant `RAW.tar`
archives, and is resumable:

```bash
mkdir -p data && cd data

for acc in GSE166676 GSE131778 GSE165816 GSE212966 \
           GSE155468 GSE156632 GSE171306; do
    stub=$(echo "$acc" | sed 's/...$//')
    series="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}nnn/${acc}/suppl"

    # 1. Per-sample files live under samples/GSMxxxnnn/GSMxxx/suppl/
    if curl -sf "${series}/filelist.txt" -o /tmp/fl.txt; then
        awk -F'\t' '($1=="File") && $2 ~ /^GSM/ {
            n=$2; split(n,a,"_"); g=a[1];
            gs=substr(g,1,length(g)-3);
            printf "https://ftp.ncbi.nlm.nih.gov/geo/samples/%snnn/%s/suppl/%s\n", gs, g, n
        }' /tmp/fl.txt | while read -r url; do
            wget -c -q --show-progress -P "../data/${acc}/" "$url"
        done

        # 2. Series-level files that are not per-sample
        awk -F'\t' '($1=="File") && $2 !~ /^GSM/ {
            printf "%s/%s\n", "'"$series"'", $2
        }' /tmp/fl.txt | while read -r url; do
            wget -c -q --show-progress -P "../data/${acc}/" "$url"
        done
    else
        # 3. Fallback for accessions without filelist.txt (e.g. GSE131778):
        #    parse the HTML directory listing instead.
        curl -s "${series}/" \
          | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"$//' \
          | grep -vE '^\?|^/|^https?:|^\.\.' \
          | while read -r n; do
              wget -c -q --show-progress -P "../data/${acc}/" "${series}/${n}"
            done
    fi
done
```

> **Do not use `wget -r` on the series directory alone.** It downloads only the
> series `suppl/` folder and silently omits the 10x matrices for GSE166676,
> GSE212966, GSE171306 and GSE214611, because those files are served from the
> `samples/` tree. This is the single most common reason a re-analysis of these
> datasets ends up with missing matrices.

> **Do not download `RAW.tar` alongside the individual files.** For GSE159677
> (2.12 GB) and GSE214611 (2.90 GB) the archive is the concatenation of the
> individual files; fetching both doubles the disk footprint for identical data.

### Skipping GSE159677 and GSE214611

Both are excluded from this reproduction, which is why the main analysis needs
1.2 GB rather than 5.8 GB:

- **GSE159677** ships a Cell Ranger **pre-scored** matrix (`molecule_info.h5`,
  1.98 GB) rather than a count matrix. The manuscript already marks this row with
  a dagger symbol as "different pipeline, not comparable"; it enters no summary
  statistic.
- **GSE214611** is the mouse atherosclerosis dataset (2.90 GB, 10x). It
  contributes exactly one row to Table 1 and would require switching to the mouse
  marker panel.

Downloading them is only necessary if a reviewer asks for those two rows to be
reproduced independently.

---

## Where the data go

`03_recompute/run_estimator_extended.py` needs the raw matrices. It searches, in
order:

1. the directory named by the `AUDIT_DATA` environment variable
2. `<this package>/data/` (that is, `04_code_and_data/../data/`)
3. `/root/audit/data` (the original working directory)

Place the per-accession directories in any one of these locations:

```bash
export AUDIT_DATA=/path/to/your/data
python3 03_recompute/run_estimator_extended.py
```

> **No other script needs the raw data.** The audit library and the three
> verification scripts read the per-dataset JSON records under
> `07_result_records/`. A reviewer can therefore check every number in the paper
> without downloading any matrices.

---

## Quality control

Processing applied before auditing:

- **Filtering**: 200-6,000 genes per cell, mitochondrial reads <= 20%
- **Normalisation**: library-size scaling followed by log transformation
- **Down-sampling**: 11 of the 16 datasets are down-sampled to a fixed cell count.
  Sampling is uniform and **blind to cell type**; the cap is the smallest sample
  in the cohort (3,000 for the eight GSE155468 patients; 2,999 / 3,500 / 4,173 /
  4,749 elsewhere)
- **Cell-type annotation**: each cell is scored against a marker panel covering
  10 human (11 mouse) populations and assigned its highest-scoring label. **The
  same panel is used for every dataset** so that annotation behaviour is
  comparable across studies

> Note: the round cell counts in Table 1 are a consequence of the down-sampling
> rule, **not** a property of the raw data. Section 2.1 of the manuscript states
> this explicitly.

---

## Citation

If you use these data, cite the original dataset publications (see the manuscript
reference list) and note GEO as the source. This audit generated no new
experimental data.
