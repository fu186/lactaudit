# Analysis code (Supplementary Software S1)

The scripts that produced every number, table and figure in the paper.

**Start with [`README_reproduction.md`](README_reproduction.md)** — it is the full
reproduction guide: directory layout, environment, run order and the three self-check
scripts. This file only records the three things that guide does not cover.

## 1. Relationship to `lactaudit` in the repository root

They are the same controls, packaged for different readers.

| | Repository root `lactaudit` | This directory |
|---|---|---|
| What it is | installable package, generalised | the analysis exactly as it was run |
| Audience | someone who wants to audit their own gene set | someone reproducing this paper |
| Scope | T1–T4, SSR, null generator, CLI | those plus dataset readers, result records, figures |
| Language | English throughout | English, with historical path literals (see §3) |
| Tests | 19 unit tests, CI on three Python versions | three self-check scripts over the deposited results |

If you only want to apply the audit to your own data, `pip install -e .` and use
`lactaudit` — you do not need anything in this directory.

## 2. The input data is not here

The nine GEO accessions are public but not redistributable with the manuscript;
roughly 5 GB compressed is needed for a full re-run. `08_raw_data_download/` holds the
download instructions and Methods describes the QC filter and the downsampling rule.

**Reading the deposited results needs no data at all** — `07_result_records/` contains the
completed run, and the self-check scripts verify the manuscript against it.

## 3. Some path literals are historical

Several scripts carry absolute paths from the machine the analysis originally ran on
(`/workspace/...`), and some directory names on that machine were in Chinese
(`审计工具包`, `论文数据表`, `模拟`, `补充材料`). Those literals are kept **verbatim**: rewriting
them would make the scripts differ from what actually produced the deposited results.

Comments beside most of them read `directory name kept as-is`. To re-run, replace those
constants with paths on your own machine; no logic depends on them.
