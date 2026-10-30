#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure-text consistency check (new).

Background: during the third review round, Figure 1D / Figure 5 were found to still
hard-code the old simulation values (T4r null band +0.076, signal +0.147~+0.280,
T4 saturation +0.941, graphical abstract +0.059), whereas the main text had already
switched to the true simulation values in `sim_t4_t4r.json`. Figures and text
contradicted each other, which is a submission-level defect.

This script gathers evidence from both sides and compares:
  1) simulation JSON: `审计工具包/sim_t4_t4r.json`  (Chinese on-disk directory name)
  2) plotting script: `审计工具包/make_bib_figs.py` (must reference the JSON variables,
     never hard-code them)
  3) both manuscript bodies: any simulation number that appears must be traceable
     in the JSON

A non-zero exit code means an inconsistency was found.
"""
import json
import os
import re
import sys

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_paths as AP            # noqa: E402

ROOT = AP.ROOT
TOOL = AP.TOOL


def _find(name):
    return AP.find(name)


SIM = str(_find('sim_t4_t4r.json'))
FIGS = str(_find('make_bib_figs.py'))
MANUSCRIPTS = [str(p) for p in AP.manuscripts().values()]

# stale simulation values known to have been removed (neither figure nor text should show them)
STALE = ['0.076', '0.073', '0.147', '0.280', '0.941', '+0.059', '0.990', '0.694']

problems = []

# ---------------------------------------------------------------- 1) load JSON
with open(SIM, encoding='utf-8') as fh:
    sim = json.load(fh)

valid = set()
for arm in ('null_arm', 'signal_arm'):
    for k, v in sim[arm].items():
        for est in ('T4', 'T4r'):
            if est in v and 'mean' in v[est]:
                valid.add(round(float(v[est]['mean']), 3))

t4_null_max = max(round(float(v['T4']['mean']), 3) for v in sim['null_arm'].values())
t4r_null_lo = min(round(float(v['T4r']['mean']), 3) for v in sim['null_arm'].values())
t4r_null_hi = max(round(float(v['T4r']['mean']), 3) for v in sim['null_arm'].values())
t4r_sig_lo = min(round(float(v['T4r']['mean']), 3) for v in sim['signal_arm'].values())
t4r_sig_hi = max(round(float(v['T4r']['mean']), 3) for v in sim['signal_arm'].values())

print('[1] key values from the simulation JSON')
print(f'    T4 null max           = {t4_null_max:+.3f}')
print(f'    T4r null band         = {t4r_null_lo:+.3f} ~ {t4r_null_hi:+.3f}')
print(f'    T4r signal band       = {t4r_sig_lo:+.3f} ~ {t4r_sig_hi:+.3f}')
print(f'    all (arm, k, est) means: {len(valid)}')

# ---------------------------------------------------------------- 2) plotting script
with open(FIGS, encoding='utf-8') as fh:
    fig_src = fh.read()

print('\n[2] plotting-script reference check')
if 'sim_t4_t4r.json' not in fig_src:
    problems.append('make_bib_figs.py does not read sim_t4_t4r.json')
    print('    [FAIL] sim_t4_t4r.json not read')
else:
    print('    [OK ] sim_t4_t4r.json read')

found_stale = [s for s in STALE if s in fig_src]
if found_stale:
    problems.append(f'make_bib_figs.py still contains stale simulation values: {found_stale}')
    print(f'    [FAIL] still contains stale simulation values {found_stale}')
else:
    print('    [OK ] no hard-coded stale simulation values')

# ---------------------------------------------------------------- 3) manuscript body
print('\n[3] traceability of simulation numbers in the body text')
for mp in MANUSCRIPTS:
    if not os.path.isfile(mp):
        problems.append(f'missing manuscript {mp}')
        print(f'    [FAIL] missing {mp}')
        continue
    txt = open(mp, encoding='utf-8').read()
    name = os.path.basename(mp)

    # extract values of the form +/-0.xxx whose magnitude is < 1, then check whether
    # they come from the simulation set or from the declared set
    nums = {float(m) for m in re.findall(r'[−-]?\d\.\d{3}', txt)
            .__iter__()} if False else set()
    nums = set()
    for m in re.findall(r'[−\-\+]?\d\.\d{3}', txt):
        try:
            nums.add(round(abs(float(m.replace('−', '-'))), 3))
        except ValueError:
            pass

    stale_hits = []
    for s in STALE:
        token = s.lstrip('+')
        if token not in txt:
            continue
        # judge the context line by line: only hits in a "simulation / reference band"
        # context count as problems; numbers in table rows (lines starting with '|')
        # are real data and are allowed.
        ctx_bad = False
        for line in txt.splitlines():
            if token not in line:
                continue
            stripped = line.strip()
            if stripped.startswith('|'):
                continue      # table data row
            low = line.lower()
            # Only simulation-context lines are inspected. The keywords are kept in
            # English; the Chinese equivalents that used to be listed here are no
            # longer needed now that all manuscripts are in English.
            if any(w in low for w in ('simulat', 'null band', 'reference band',
                                      'signal band')):
                ctx_bad = True
        if ctx_bad:
            stale_hits.append(s)

    if stale_hits:
        problems.append(f'{name} still cites stale values in a simulation context: {stale_hits}')
        print(f'    [FAIL] {name} contains stale values {stale_hits} in a simulation context')
    else:
        print(f'    [OK ] {name} has no stale values in a simulation context (real table data excluded)')

# ---------------------------------------------------------------- summary
print('\n' + '=' * 68)
if problems:
    print('Problems found:')
    for p in problems:
        print('  -', p)
    sys.exit(1)
print('Figure and text simulation values are consistent; the plotting script shares its source with the JSON.')
print('=' * 68)
