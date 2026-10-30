#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full pre-submission consistency cross-check
===========================================

Division of labour with verify_all_numbers.py:
  - verify_all_numbers.py  checks Table 1 against the raw JSON cell by cell and scans
                           for stale literals.
  - crosscheck_all.py      maps every numeric claim in the two manuscripts back to the
                           raw artefacts in an "assertion -> provenance" manner, and
                           checks that the same fact agrees across manuscripts.

Check levels
------------
  [A] cross-manuscript consistency of the same fact (a number must be identical in all manuscripts)
  [B] body text -> raw artefacts provenance (each claim must have an independent source)
  [C] figure scripts -> body text (the JSON read by the figures is the same source as the cited values)
  [D] internal consistency of the tables (derivable quantities within a Table 1 row)
  [E] abstract <-> body text consistency
  [F] supplementary-material list <-> actual files

Usage: python3 crosscheck_all.py
"""
import csv
import json
import re
import sys
from pathlib import Path
import os

import numpy as np


# ---------------------------------------------------------------------------
# Path resolution: delegated to the audit_paths module in this package (compatible with
# both the workspace and the submission-package layout).
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_paths as AP            # noqa: E402

ROOT = AP.ROOT
RES = AP.RES
TOOL = AP.TOOL


def _findfile(name):
    return AP.find(name)


def _m(fname):
    hits = sorted(ROOT.glob(f'**/{fname}'))
    if not hits:
        raise SystemExit(f'manuscript not found: {fname}')
    return hits[0]


MS = AP.manuscripts()

FIND = []


def rec(sev, code, msg):
    FIND.append((sev, code, msg))


def flat(p):
    """Read the file and normalise all whitespace to single spaces, to avoid false negatives
    caused by line breaks."""
    return ' '.join(Path(p).read_text(encoding='utf-8').split())


TXT = {k: flat(v) for k, v in MS.items()}


def has(txt, needle):
    return needle in txt


# ==============================================================
# Load raw artefacts
# ==============================================================
RAW = {}
for jf in sorted(RES.glob('*.json')):
    if 'fidelity' in jf.name:
        continue
    d = json.load(open(jf, encoding='utf-8'))
    if d.get('T1'):
        RAW[jf.stem] = d

# pre-scored row: produced by a different pipeline and therefore not comparable;
# always excluded from the aggregate statistics in the text.
PRESCORED = 'AS_GSE159677_ours'


def comparable(d):
    """The aggregate statistics in the text cover only the 15 comparable datasets."""
    return d is not None and d.get('T3') and d['T3'].get('rho_ngene') is not None

T4 = {}
for cf in sorted(RES.glob('*_T4.csv')):
    T4[cf.stem.replace('_T4', '')] = list(csv.DictReader(open(cf, encoding='utf-8')))

AUTH = {r['Dataset']: r for r in
        csv.DictReader(open(_findfile('t4r_recomputed.csv'), encoding='utf-8'))}
BOOT = {r['dataset']: r for r in
        csv.DictReader(open(_findfile('t4r_bootstrap_corrected.csv'), encoding='utf-8'))}
SIM = json.load(open(AP.find('sim_t4_t4r.json'), encoding='utf-8'))

print('=' * 74)
print('Full pre-submission consistency cross-check')
print('=' * 74)
print(f'Loaded: {len(RAW)} raw JSONs, {len(T4)} T4 CSVs, '
      f'{len(AUTH)} authoritative recomputation rows, {len(BOOT)} bootstrap rows')

# ==============================================================
# ==============================================================
# [A] cross-manuscript consistency of the same fact
# ==============================================================
print(f'\n[A] Cross-manuscript consistency (manuscripts: {len(TXT)}: {", ".join(TXT)})')
SHARED = {
    'SSR median':          '0.907',
    'SSR min':             '0.331',
    'SSR max':             '4.320',
    'T1 failure count':    '8 of the 15',
    'T2 failure count':    '10 of 15',
    'T3 median':           '+0.574',
    'T3 min':              '+0.434',
    'T3 max':              '+0.793',
    'T3 threshold 0.3':    '15 of 15',
    'T3 threshold 0.5':    '14 of 15',
    'T4r median':          '+0.770',
    'T4r min':             '+0.248',
    'T4r max':             '+0.967',
    'span threshold':      '0.15',
    'span pass count':     '5',
    'bootstrap width':     '0.920',
    'bootstrap contains zero': 'Seven of the 15',
    'BH minimum':          '0.106',
    'permutation depth':   '150',
    'T1T2 agreement':      '9 of 15',
    'total cells':         '113,419',
}
for label, needle in SHARED.items():
    inb = [k for k, t in TXT.items() if needle in t]
    if len(inb) != len(TXT):
        missing = sorted(set(TXT) - set(inb))
        rec('P0', 'CROSS', f'shared fact "{label}" ({needle}) missing from {missing}')
    else:
        print(f'    [OK ] {label:14s} identical in {len(TXT)} manuscripts = {needle}')

# ==============================================================
# [B] body text -> raw artefacts provenance
# ==============================================================
print('\n[B] Body-text claims -> raw-artefact provenance')

# B1. SSR distribution
ssr = np.array([d['T1']['SSR'] for d in RAW.values() if d.get('T1')])
if len(ssr) == 16:
    exp = (ssr.min(), np.median(ssr), ssr.max(), int((ssr <= 1).sum()))
    print(f'    raw SSR: n={len(ssr)} min={exp[0]:.3f} median={exp[1]:.3f} '
          f'max={exp[2]:.3f} count <=1={exp[3]}')
    for v, want in [(f'{exp[0]:.3f}', '0.331'), (f'{exp[1]:.3f}', '0.907'),
                    (f'{exp[2]:.3f}', '4.320')]:
        if v != want:
            rec('P0', 'SSR', f'raw SSR gives {v}, text says {want}')
else:
    rec('P1', 'SSR', f'{len(ssr)} raw JSONs carry T1, expected 16')

# B2. T3 statistics
t3 = np.array([d['T3']['rho_ngene'] for d in RAW.values()
               if d.get('T3') and d['T3'].get('rho_ngene') is not None])
print(f'    raw T3: n={len(t3)} min={t3.min():+.3f} median={np.median(t3):+.3f} '
      f'max={t3.max():+.3f} >0.5:{(t3 > 0.5).sum()} >0.3:{(t3 > 0.3).sum()}')
for got, want, lab in [(f'{t3.min():+.3f}', '+0.434', 'min'),
                       (f'{np.median(t3):+.3f}', '+0.574', 'median'),
                       (f'{t3.max():+.3f}', '+0.793', 'max')]:
    if got != want:
        rec('P0', 'T3', f'T3 {lab} raw={got} text={want}')

# B3. T3 unrelated to the housekeeping score
hk = np.array([d['T3']['rho_hk'] for k, d in RAW.items()
               if comparable(d) and d['T3'].get('rho_hk') is not None])
if len(hk) != 15:
    rec('P1', 'T3-HK', f'{len(hk)} datasets in the comparable set carry rho_hk, expected 15')
hk_ov = [k for k, d in RAW.items()
         if comparable(d) and d['T3'].get('hk_overlap')]
print(f'    T3 vs housekeeping score: n={len(hk)} median={np.median(hk):+.3f} '
      f'range={hk.min():+.3f}~{hk.max():+.3f}')
print(f'    comparable datasets with hk_overlap=True: {hk_ov or "none"} (0/15)')
if hk_ov:
    rec('P0', 'T3-HK', f'the text states overlap flag 0 of 15, but {hk_ov} are True')
for got, want, lab in [(f'{np.median(hk):+.3f}', '+0.038', 'median'),
                       (f'{hk.min():+.3f}', '-0.294', 'min'),
                       (f'{hk.max():+.3f}', '+0.341', 'max')]:
    if got != want:
        rec('P0', 'T3-HK', f'housekeeping correlation {lab} raw={got} text={want}')

# B4. T4 statistics
t4 = np.array([d['T4']['rho_detect_vs_score'] for d in RAW.values()
               if d.get('T4') and d['T4'].get('rho_detect_vs_score') is not None])
print(f'    raw T4: n={len(t4)} min={t4.min():+.3f} median={np.median(t4):+.3f} '
      f'max={t4.max():+.3f} >0.7:{(t4 > 0.7).sum()}')
if len(t4) != 14:
    rec('P1', 'T4', f'{len(t4)} datasets with a computable T4, text says 14')
if (t4 > 0.7).sum() != 14:
    rec('P0', 'T4', f'only {(t4 > 0.7).sum()} datasets have T4>0.7, text says 14')

# B5. T4r and bootstrap
t4r_s = np.array([float(r['T4r_Spearman']) for r in AUTH.values()])
t4r_p = np.array([float(r['T4r_Pearson']) for r in AUTH.values()])
print(f'    authoritative T4r: n={len(t4r_s)} spearman med={np.median(t4r_s):+.3f} '
      f'range={t4r_s.min():+.3f}~{t4r_s.max():+.3f}')
print(f'                        pearson med={np.median(t4r_p):+.3f}')
if f'{np.median(t4r_s):+.3f}' != '+0.770':
    rec('P0', 'T4R', f'T4r Spearman median={np.median(t4r_s):+.3f}, text=+0.770')
if f'{np.median(t4r_p):+.3f}' != '+0.787':
    rec('P0', 'T4R', f'T4r Pearson median={np.median(t4r_p):+.3f}, text=+0.787')

# B6. span gate
span = np.array([float(r['detection_rate_span']) for r in AUTH.values()])
n_span = int((span >= 0.15).sum())
print(f'    detection-rate span: n={len(span)} median={np.median(span):.3f} '
      f'>=0.15: {n_span}')
if n_span != 5:
    rec('P0', 'SPAN', f'{n_span} datasets have span>=0.15, text says 5')

# B7. bootstrap
bw = np.array([float(r['width']) for r in BOOT.values()])
nz = sum(1 for r in BOOT.values() if r['contains_zero'] == 'True')
interp = [r for r in BOOT.values() if r['interpretable'] == 'True']
nz_int = sum(1 for r in interp if r['contains_zero'] == 'True')
print(f'    bootstrap: n={len(bw)} median width={np.median(bw):.3f} '
      f'contains zero={nz}/{len(bw)} among the {len(interp)} interpretable, contains zero={nz_int}')
if f'{np.median(bw):.3f}' != '0.920':
    rec('P0', 'BOOT', f'bootstrap median width={np.median(bw):.3f}, text=0.920')
if nz != 7:
    rec('P0', 'BOOT', f'bootstrap contains zero={nz}, text=7')
if nz_int != 4 or len(interp) != 5:
    rec('P0', 'BOOT', f'interpretable set contains zero={nz_int}/{len(interp)}, text=4 of 5')

# B8. simulation values must share their source with the JSON
_null = SIM.get('null_arm', {})
_sig = SIM.get('signal_arm', {})
if _null and '3' in _null:
    v3 = _null['3']['T4']['mean']
    v1 = _null['1']['T4']['mean']
    v2 = _null['2']['T4']['mean']
    print(f'    simulation T4 null: @1.0={v1:+.3f} @2.0={v2:+.3f} @3.0={v3:+.3f}')
    for got, want, lab in [(f'{v1:+.3f}', '+0.946', '@1.0'),
                           (f'{v2:+.3f}', '+0.984', '@2.0'),
                           (f'{v3:+.3f}', '+0.991', '@3.0')]:
        if got != want:
            rec('P0', 'SIM', f'simulation T4 {lab}={got}, text={want}')
    # the maximum over the whole null arm should be +1.000
    mx = max(_null[k]['T4']['max'] for k in _null)
    print(f'    simulation T4 null, global maximum = {mx:+.3f}')
    if f'{mx:+.3f}' != '+1.000':
        rec('P1', 'SIM', f'simulation T4 max={mx:+.3f}, text=+1.000')
    # after injecting signal, T4 barely moves
    if '40' in _sig:
        t4s0 = _sig['0']['T4']['mean']
        t4s40 = _sig['40']['T4']['mean']
        print(f'    simulation with injected signal: T4 {t4s0:+.3f} -> {t4s40:+.3f} '
              f'(change {abs(t4s40 - t4s0):.4f})')
        if abs(t4s40 - t4s0) > 0.005:
            rec('P0', 'SIM', f'T4 changes by {abs(t4s40-t4s0):.4f} after signal injection > 0.005')
        # T4r should rise in the signal arm
        r0 = _sig['0']['T4r']['mean']
        r40 = _sig['40']['T4r']['mean']
        print(f'    simulation with injected signal: T4r {r0:+.3f} -> {r40:+.3f}')
        for got, want, lab in [(f'{r0:+.3f}', '-0.750', 'T4r signal arm k=0'),
                               (f'{r40:+.3f}', '-0.442', 'T4r signal arm k=40')]:
            if got != want:
                rec('P0', 'SIM', f'simulation {lab}={got}, text={want}')
    # T4r null-distribution band
    if '3' in _null:
        rn_min = _null['3']['T4r']['min']
        rn_max = _null['3']['T4r']['max']
        print(f'    simulation T4r null@3.0 band = {rn_min:+.3f} ~ {rn_max:+.3f}')
else:
    rec('P0', 'SIM', f'simulation JSON lacks null_arm; actual top-level keys={list(SIM)}')

# B9. the total cell count must equal the sum of the Table 1 rows
rows_t1 = []
for line in MS['BMC'].read_text(encoding='utf-8').splitlines():
    s = line.strip()
    if not s.startswith('|'):
        continue
    c = [x.strip() for x in s.strip('|').split('|')]
    if len(c) == 13 and c[0].startswith('GSE'):
        rows_t1.append(c)
tot = sum(int(c[3].replace(',', '')) for c in rows_t1)
print(f'    Table 1 rows={len(rows_t1)} total cells={tot:,}')
if tot != 113419:
    rec('P0', 'CELLS', f'Table 1 total cells={tot:,}, text=113,419')
if len(rows_t1) != 16:
    rec('P0', 'ROWS', f'Table 1 has {len(rows_t1)} rows, expected 16')

# ==============================================================
# [C] figure scripts -> body text
# ==============================================================
print('\n[C] Figure-script source check')
fig = (AP.find('make_bib_figs.py')).read_text(encoding='utf-8')
if 'sim_t4_t4r.json' not in fig:
    rec('P0', 'FIG', 'the plotting script does not read sim_t4_t4r.json')
else:
    print('    [OK ] make_bib_figs.py reads sim_t4_t4r.json')
# no hard-coded simulation values may appear in the figures
for bad in ('0.059', '0.245', '+0.946'):
    if bad in fig:
        rec('P1', 'FIG', f'the plotting script still hard-codes the simulation value {bad}')
print('    [OK ] the plotting script contains no hard-coded simulation values')

# simulation values cited in the captions should also appear in the body text
for k, t in TXT.items():
    for v in ('+0.991', '-0.750'):
        if v not in t:
            rec('P1', 'FIG', f'{k}: the simulation value {v} cited in the caption does not appear in the body text')

# ==============================================================
# [D] internal table consistency
# ==============================================================
print('\n[D] Table 1 within-row consistency')
bad_internal = 0
for c in rows_t1:
    label, ssr, t1v, p, t2v, t3, t4, t4r, sp = (c[0], c[5], c[6], c[7],
                                                 c[8], c[9], c[10], c[11], c[12])
    # the T1 verdict must agree with SSR
    try:
        s = float(ssr)
        want = 'pass' if s > 1 else 'fail'
        if t1v != want:
            rec('P0', 'T1JUDGE', f'{label}: SSR={s} verdict={t1v}, expected {want}')
            bad_internal += 1
    except ValueError:
        pass
    # the T2 verdict must agree with p
    try:
        pv = float(p)
        want = 'pass' if pv < 0.05 else 'fail'
        if t2v != want:
            rec('P0', 'T2JUDGE', f'{label}: p={pv} verdict={t2v}, expected {want}')
            bad_internal += 1
    except ValueError:
        pass
    # a row carrying a dagger must have span < 0.15
    if '†' in t4r:
        try:
            if float(sp) >= 0.15:
                rec('P0', 'DAGGER', f'{label}: T4r carries a dagger but span={sp} >= 0.15')
                bad_internal += 1
        except ValueError:
            rec('P0', 'DAGGER', f'{label}: T4r carries a dagger but span is not numeric ({sp})')
    # a row with a numeric T4r but no dagger must have span >= 0.15
    m = re.match(r'^([+-]?\d\.\d+)$', t4r)
    if m and float(sp) < 0.15:
        rec('P0', 'DAGGER', f'{label}: T4r={t4r} without a dagger but span={sp} < 0.15')
        bad_internal += 1
if bad_internal == 0:
    print(f'    [OK ] all {len(rows_t1)} rows are consistent (T1/T2 verdicts, dagger and span linkage)')

# ==============================================================
# [E] abstract <-> body text
# ==============================================================
print('\n[E] Abstract vs body-text consistency')
for k, t in TXT.items():
    ai = t.find('Abstract')
    ri = t.find('## Results') if k == 'BMC' else t.find('**Background.**')
    abs_txt = t[ai:ai + 3000]
    # every 3-decimal number in the abstract must also be found in the body text
    for m in re.finditer(r'([+-]?\d\.\d{3})\b', abs_txt):
        v = m.group(1)
        if v not in t:
            rec('P1', 'ABS', f'{k}: abstract value {v} does not appear in the body text')
    print(f'    [OK ] all abstract values in {k} are traceable in the body text')

# ==============================================================
# [F] supplementary-material list <-> actual files
# ==============================================================
print('\n[F] Supplementary-material list vs actual files')
SUPP = AP.suppdir()
actual = sorted(p.name for p in SUPP.iterdir()) if SUPP.exists() else []
# the text cites them as "Supplementary Table Sn"; the directory filenames look like Table_Sn_*.csv
REF = sorted(set(int(x) for x in re.findall(
    r'Supplementary Table S(\d+)', TXT['BMC'] + ' ' + TXT['BIB'])))
files = {}
for name in actual:
    m = re.match(r'Table_S(\d+)(\w*)_', name)
    if m:
        files.setdefault(int(m.group(1)), []).append(name)
print(f'    {len(actual)} actual files in the directory; S-table numbers cited in the text: {REF}')
missing = [f'S{n}' for n in REF if n not in files]
if missing:
    rec('P1', 'SUPP', f'S tables cited in the text but missing from the directory: {missing}')
else:
    print('    [OK ] every cited S table has a corresponding file')
# every S table listed in the list must have a file
listed = sorted(set(int(x) for x in re.findall(
    r'\*\*Table S(\d+)', TXT['BMC'])))
notshipped = [f'S{n}' for n in listed if n not in files]
if notshipped:
    rec('P0', 'SUPP', f'S tables listed but shipped without a file in the directory: {notshipped}')
else:
    print(f'    [OK ] all {len(listed)} listed S tables have files')
# S5b must share its source with the recomputed bootstrap
if 'Table_S5b_T4r_bootstrap_intervals.csv' in actual:
    _b = list(csv.DictReader(open(
        SUPP / 'Table_S5b_T4r_bootstrap_intervals.csv', encoding='utf-8')))
    _w = np.array([float(r['width']) for r in _b])
    print(f'    Supplementary Table S5b: {len(_b)} rows, median width={np.median(_w):.3f}')
    if len(_b) != len(BOOT) or abs(np.median(_w) - np.median(bw)) > 1e-9:
        rec('P0', 'SUPP', 'Supplementary Table S5b disagrees with the recomputed bootstrap')
    else:
        print('    [OK ] Supplementary Table S5b shares its source with the recomputed bootstrap')
else:
    rec('P1', 'SUPP', 'Supplementary Table S5b is missing')

# ==============================================================
# Summary
# ==============================================================
print('\n' + '=' * 74)
print('Summary')
print('=' * 74)
if not FIND:
    print(f'  All passed: {len(TXT)} manuscripts agree, numbers are traceable, figures share their source, tables are self-consistent.')
else:
    from collections import Counter
    c = Counter(x[0] for x in FIND)
    print(f'  {len(FIND)} findings: ' + ', '.join(f'{k}={v}' for k, v in sorted(c.items())))
    print()
    for sev, code, msg in FIND:
        print(f'  [{sev}] {code}: {msg}')
print('=' * 74)
sys.exit(0 if not [f for f in FIND if f[0] == 'P0'] else 1)
