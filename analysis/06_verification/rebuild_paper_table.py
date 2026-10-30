#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild paper_table.csv from the authoritative audit results (the data source for the plotting script).

Data sources (all are raw artefacts that package scripts can recompute):
  - /workspace/审计工具包/results/*.json      -> SSR / T1 / T2 p / T3 rho_ngene / T4 rho_detect
  - /workspace/BIB论文_评审回应_第三轮/复算/t4r_recomputed.csv -> T4r / span / cluster
  (both paths above are on-disk paths with Chinese directory names, kept verbatim)

Output columns are kept identical to what make_bib_figs.py expects:
  dataset, T4_rho, SSR, T2_p, T3_rho_ngene, T4r, span
"""
import csv
import json
import os
from pathlib import Path
import sys
import glob

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_paths as AP            # noqa: E402

RESULTS = str(AP.RES)
REC = str(AP.find('t4r_recomputed.csv'))
OUT = str(AP.RES / 'paper_table.csv')

# JSON filename -> dataset label used in the manuscript
LABEL = {
    'AAA_GSE166676': 'AAA_GSE166676',
    'AS_GSE159677_ours': 'AS_GSE159677_ours',
    'AS_mouse_GSE214611': 'AS_mouse_GSE214611',
    'CAD_GSE131778': 'CAD_GSE131778',
    'DFU_GSE165816': 'DFU_GSE165816',
    'PDAC_GSE212966': 'PDAC_GSE212966',
    'ccRCC_GSE156632': 'ccRCC_GSE156632',
    'ccRCC2_GSE171306': 'ccRCC2_GSE171306',
    'TAA_GSE155468': 'TAA_TAA1',
}
for i in range(2, 9):
    LABEL[f'TAA_TAA{i}_GSE155468'] = f'TAA_TAA{i}'

# GEO accession (cannot be inferred from the source column of the T4 CSV, so parse it from the filename)
def geo_of(stem):
    for g in ('GSE166676', 'GSE159677', 'GSE214611', 'GSE131778', 'GSE165816',
              'GSE212966', 'GSE156632', 'GSE171306', 'GSE155468'):
        if g in stem:
            return g
    return ''


def main():
    # ---- T4r recomputation table (includes span) ----
    t4r = {}
    with open(REC, encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            t4r[r['Dataset'].strip()] = r

    rows = []
    for jf in sorted(glob.glob(f'{RESULTS}/*.json')):
        stem = os.path.basename(jf)[:-5]
        if stem in ('audit_summary', 'fidelity_out'):
            continue
        d = json.load(open(jf, encoding='utf-8'))

        t1 = d.get('T1') or {}
        t2 = d.get('T2') or {}
        t3 = d.get('T3') or {}
        t4 = d.get('T4') or {}

        # GSE accession -> label (de-duplicated)
        geo = geo_of(stem)
        # match against the T4r recomputation table: its Dataset column is the GSE accession;
        # the multiple TAA patients must be matched in order
        label = LABEL.get(stem, stem)

        # T4r / span: preferably matched by the GSE extracted from label + patient index
        t4r_v = span_v = cluster = ''
        key = geo
        if stem.startswith('TAA_') and stem != 'TAA_GSE155468':
            # each TAA patient has its own T4 CSV; the recomputation table lists them in
            # order as repeated records of GSE155468
            pass

        if geo in t4r:
            pass  # with multiple patients, overwritten in order below

        rows.append({
            'dataset': label,
            'geo': geo,
            'stem': stem,
            'T1_SSR': t1.get('SSR'),
            'T2_p': t2.get('empirical_p'),
            'T3_rho_ngene': t3.get('rho_ngene'),
            'T4_rho': t4.get('rho_ngene_vs_score'),
            'patient': label.split('_')[-1] if label.startswith('TAA_TAA') else '',
            'top': d.get('reproduced_top', ''),
        })

    # ---- merge the T4r recomputation table in order ----
    # the Dataset column of t4r_recomputed.csv holds the GSE accession; TAA patients
    # appear in order as TAA1..TAA8
    rec = []
    with open(REC, encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            rec.append(r)

    # matching order: TAA patients are matched uniquely by (TAA1..TAA8) suffix;
    # all others are matched uniquely by GSE accession
    non_taa = {r['Dataset'].strip(): r for r in rec
               if 'TAA' not in r['Dataset']}
    taa_recs = {}
    for r in rec:
        ds = r['Dataset'].strip()
        if 'TAA' in ds:
            # of the form 'GSE155468 (TAA1)'
            tag = ds.split('(')[-1].rstrip(')').strip()
            taa_recs[tag] = r

    for row in rows:
        geo = row.pop('geo')
        stem = row.pop('stem')
        if geo == 'GSE155468':
            tag = row['dataset'].split('_')[-1]   # TAA1 .. TAA8
            r = taa_recs.get(tag)
            if r:
                row['T4r'] = r['T4r_Spearman']
                row['span'] = r['detection_rate_span']
                row['cluster'] = r['Disease']
            else:
                row['T4r'] = row['span'] = row['cluster'] = ''
        elif geo in non_taa:
            r = non_taa[geo]
            row['T4r'] = r['T4r_Spearman']
            row['span'] = r['detection_rate_span']
            row['cluster'] = r['Disease']
        else:
            row['T4r'] = row['span'] = row['cluster'] = ''

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cols = ['dataset', 'cluster', 'T1_SSR', 'T2_p', 'T3_rho_ngene',
            'T4_rho', 'T4r', 'span', 'patient', 'top']
    with open(OUT, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f'wrote {OUT}  ({len(rows)} rows)')
    for row in rows:
        print('  ', row['dataset'], '| T4r=', row['T4r'], '| span=', row['span'])


if __name__ == '__main__':
    main()
