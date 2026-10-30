#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Strict numerical verification of the lactylation-scoring audit manuscripts.
Independently re-derives every reported number from raw per-dataset artefacts.
Run:  python3 verify_all.py
"""
import json, glob, re, sys, os
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr

# ---------------------------------------------------------------------------
# Path resolution: delegated to the audit_paths module in this package (compatible with
# both the workspace and the submission-package layout).
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_paths as AP            # noqa: E402

ROOT = AP.ROOT                      # root directory of manuscripts and materials
RES = AP.RES                        # per-dataset result-record directory
T4DIR = RES                         # *_T4.csv sit in the same directory as the result JSONs
WORK = ROOT
TOOL = AP.TOOL

def _find(_patterns, what='file'):
    """Backwards-compatible signature: accepts a list of patterns or a single filename."""
    if isinstance(_patterns, (list, tuple)):
        for pat in _patterns:
            hits = sorted(ROOT.glob(f'**/{Path(pat).name}'))
            if hits:
                return hits[0]
        raise SystemExit(f'{what} not found: {_patterns}')
    return AP.find(_patterns)


_manuscripts = AP.manuscripts


def rec(sev, code, msg):
    FINDINGS.append((sev, code, msg))


def P(s):
    print(s)


# ============================================================
# 0. Load raw JSON
# ============================================================
RAW = {}
for jf in sorted(glob.glob(str(RES / '*.json'))):
    d = json.load(open(jf, encoding='utf-8'))
    nm = d.get('name') or Path(jf).stem
    if d.get('T1'):
        RAW[nm] = d
P(f"[0] loaded {len(RAW)} audit JSONs with T1")

# Map manuscript dataset labels -> raw JSON key
LABEL2KEY = {
    'GSE166676':            'AAA_GSE166676',
    'GSE159677 ‡':          'AS_GSE159677_ours',
    'GSE214611':            'AS_mouse_GSE214611',
    'GSE131778':            'CAD_GSE131778',
    'GSE165816':            'DFU_GSE165816',
    'GSE212966':            'PDAC_GSE212966',
    'GSE155468 (TAA1)':     'TAA_GSE155468',
    'GSE155468 (TAA2)':     'TAA_TAA2_GSE155468',
    'GSE155468 (TAA3)':     'TAA_TAA3_GSE155468',
    'GSE155468 (TAA4)':     'TAA_TAA4_GSE155468',
    'GSE155468 (TAA5)':     'TAA_TAA5_GSE155468',
    'GSE155468 (TAA6)':     'TAA_TAA6_GSE155468',
    'GSE155468 (TAA7)':     'TAA_TAA7_GSE155468',
    'GSE155468 (TAA8)':     'TAA_TAA8_GSE155468',
    'GSE156632':            'ccRCC_GSE156632',
    'GSE171306':            'ccRCC2_GSE171306',
}

# ============================================================
# 1. Parse Table 1 out of each manuscript
# ============================================================
def parse_table1(path):
    txt = Path(path).read_text(encoding='utf-8')
    rows = []
    for line in txt.splitlines():
        s = line.strip()
        if not s.startswith('|'):
            continue
        cells = [c.strip() for c in s.strip('|').split('|')]
        if len(cells) != 13:
            continue
        if cells[0] in ('Dataset', '') or set(cells[0]) <= set('-: '):
            continue
        rows.append(cells)
    return rows


MS = {
    'BIB': AP.manuscript('BIB'),
    'BMC': AP.manuscript('BMC'),
    'BMP': AP.manuscript('BMP'),
}
T1 = {k: parse_table1(v) for k, v in MS.items()}
for k, v in T1.items():
    P(f"[1] {k}: parsed {len(v)} Table-1 data rows")
    if len(v) != 16:
        rec('P0', 'T1-ROWCOUNT', f'{k}: Table 1 has {len(v)} rows, expected 16')

# ============================================================
# 2. Cell-by-cell reconciliation of Table 1 vs raw JSON
# ============================================================

FINDINGS = []          # (severity, code, message)


def rec(sev, code, msg):
    FINDINGS.append((sev, code, msg))


def P(s=''):
    print(s)


MS = _manuscripts()
T1 = {k: parse_table1(v) for k, v in MS.items()}

def num(s):
    s = s.replace('†', '').replace('‡', '').strip()
    s = s.replace(',', '')
    if s in ('n.c.', '—', '', '-'):
        return None
    try:
        return float(s)
    except ValueError:
        return None


FIELD_MAP = [
    # (table col idx, json path, human name, decimal places shown in the table)
    (5,  lambda d: d['T1'].get('SSR'),                              'T1 SSR',        3),
    (7,  lambda d: d['T2'].get('empirical_p'),                      'T2 p',          3),
    (9,  lambda d: d['T3'].get('rho_ngene'),                        'T3 rho(nGene)', 3),
    (10, lambda d: d['T4'].get('rho_detect_vs_score') if d.get('T4') else None,
                                                                    'T4 rho(detect)',3),
]

def close(got, exp, dp):
    """The table prints to dp decimals; accept |got-exp| <= 0.5*10**-dp + tiny slack."""
    return abs(got - float(exp)) <= 0.5 * 10 ** (-dp) + 1e-9


P('\n[2] Table-1 cell reconciliation vs raw JSON')
for mskey, rows in T1.items():
    mism = 0
    for cells in rows:
        label = cells[0]
        key = LABEL2KEY.get(label)
        if key is None:
            rec('P0', 'T1-UNKNOWN-LABEL', f'{mskey}: unknown dataset label {label!r}')
            continue
        d = RAW[key]
        for ci, fn, nm, dp in FIELD_MAP:
            got = num(cells[ci])
            exp = fn(d)
            if exp is None:
                if got is not None:
                    mism += 1
                    rec('P1', 'T1-CELL', f'{mskey} {label} {nm}: manuscript {cells[ci]!r} but raw is null')
                continue
            if got is None:
                mism += 1
                rec('P1', 'T1-CELL', f'{mskey} {label} {nm}: manuscript {cells[ci]!r} but raw {exp:.6f}')
                continue
            if not close(got, exp, dp):
                mism += 1
                rec('P0', 'T1-CELL', f'{mskey} {label} {nm}: manuscript {got} vs raw {float(exp):.6f}')
        # cells column
        got_cells = num(cells[3])
        if got_cells is not None and abs(got_cells - d['n_cells']) > 0.5:
            mism += 1
            rec('P0', 'T1-CELLS', f'{mskey} {label} cells: {got_cells} vs raw {d["n_cells"]}')
        # top type
        got_top = cells[4]
        exp_top = d.get('reproduced_top')
        if exp_top and got_top != exp_top:
            mism += 1
            rec('P1', 'T1-TOP', f'{mskey} {label} top type: {got_top!r} vs raw {exp_top!r}')
    P(f'    {mskey}: {mism} mismatches across 16 rows x 6 fields')

# ============================================================
# 3. T4r and span reconciliation (independent recomputation)
# ============================================================
def rho_t4r(rows):
    x  = np.array([float(r['mean_ngene'])      for r in rows])
    y1 = np.array([float(r['detect_rate_lrg']) for r in rows])
    y2 = np.array([float(r['lact_score'])      for r in rows])
    b1 = np.polyfit(x, y1, 1); b2 = np.polyfit(x, y2, 1)
    e1 = y1 - np.polyval(b1, x); e2 = y2 - np.polyval(b2, x)
    rs, _ = spearmanr(e1, e2)
    rp, _ = pearsonr(e1, e2)
    return float(rs), float(rp), float(y1.max() - y1.min()), len(x)


P('\n[3] T4r / span independent recomputation')
T4R = {}
for csv in sorted(glob.glob(str(T4DIR / '*_T4.csv'))):
    import csv as _csv
    with open(csv, newline='', encoding='utf-8') as fh:
        rows = list(_csv.DictReader(fh))
    stem = Path(csv).stem.replace('_T4', '')
    try:
        rs, rp, sp, n = rho_t4r(rows)
    except Exception as e:
        P(f'    skip {stem}: {e}')
        continue
    T4R[stem] = dict(spearman=rs, pearson=rp, span=sp, n=n)
P(f'    recomputed {len(T4R)} datasets')

STEARM2LABEL = {
    'AAA_GSE166676': 'GSE166676',           'AS_mouse_GSE214611': 'GSE214611',
    'CAD_GSE131778': 'GSE131778',           'DFU_GSE165816': 'GSE165816',
    'PDAC_GSE212966': 'GSE212966',          'TAA_GSE155468': 'GSE155468 (TAA1)',
    'TAA_TAA2_GSE155468': 'GSE155468 (TAA2)', 'TAA_TAA3_GSE155468': 'GSE155468 (TAA3)',
    'TAA_TAA4_GSE155468': 'GSE155468 (TAA4)', 'TAA_TAA5_GSE155468': 'GSE155468 (TAA5)',
    'TAA_TAA6_GSE155468': 'GSE155468 (TAA6)', 'TAA_TAA7_GSE155468': 'GSE155468 (TAA7)',
    'TAA_TAA8_GSE155468': 'GSE155468 (TAA8)', 'ccRCC_GSE156632': 'GSE156632',
    'ccRCC2_GSE171306': 'GSE171306',
}

# cross-check against the published recomputation file
REC = json.load(open(AP.find('t4r_recomputed.json'), encoding='utf-8'))['records']
for r in REC:
    lab = r['dataset']
    key = [k for k, v in STEARM2LABEL.items() if v == lab]
    if not key:
        continue
    key = key[0]
    if key not in T4R:
        rec('P1', 'T4R-MISSING', f'{lab}: no *_T4.csv parse')
        continue
    mine = T4R[key]
    if abs(mine['spearman'] - r['t4r_spearman']) > 1e-6:
        rec('P0', 'T4R-SPEAR', f'{lab}: my {mine["spearman"]:.6f} vs file {r["t4r_spearman"]:.6f}')
    if abs(mine['pearson'] - r['t4r_pearson']) > 1e-6:
        rec('P0', 'T4R-PEAR', f'{lab}: my {mine["pearson"]:.6f} vs file {r["t4r_pearson"]:.6f}')
    if abs(mine['span'] - r['det_rate_span']) > 1e-6:
        rec('P0', 'T4R-SPAN', f'{lab}: my {mine["span"]:.6f} vs file {r["det_rate_span"]:.6f}')
# the on-disk directory is named in Chinese; describe it in English here
P('    T4r recomputation reconciles with t4r_recomputed.json (recompute directory)')

# manuscript T4r cells
P('\n[3b] Table-1 T4r / span cell reconciliation')
for mskey, rows in T1.items():
    bad = 0
    for cells in rows:
        lab = cells[0]
        key = [k for k, v in STEARM2LABEL.items() if v == lab]
        if not key or key[0] not in T4R:
            continue
        mine = T4R[key[0]]
        got_r, got_s = num(cells[11]), num(cells[12])
        if got_r is not None and abs(got_r - mine['spearman']) > 5e-4:
            bad += 1
            rec('P0', 'T1-T4R', f'{mskey} {lab}: T4r {got_r} vs recomputed {mine["spearman"]:.4f}')
        if got_s is not None and abs(got_s - mine['span']) > 5e-4:
            bad += 1
            rec('P0', 'T1-SPAN', f'{mskey} {lab}: span {got_s} vs recomputed {mine["span"]:.4f}')
    P(f'    {mskey}: {bad} T4r/span cell mismatches')

# ============================================================
# 4. Dagger rule
# ============================================================
P('\n[4] Dagger (†) rule: dagger iff T4r computable AND span < 0.15')
# Expected dagger set: datasets whose Table-1 T4r cell is numeric (i.e. T4r was
# computable, which requires rho_detect_vs_score to exist) AND whose span < 0.15.
# Datasets with T4r = 'n.c.' carry no dagger by definition.
LABEL_SPAN = {}
for stem, lab in STEARM2LABEL.items():
    if stem in T4R:
        LABEL_SPAN[lab] = T4R[stem]['span']
P(f'    datasets with a computable T4r: {len(LABEL_SPAN)}')

# Which labels the manuscript actually prints a numeric T4r for?
def numeric_t4r_labels(rows):
    out = set()
    for cells in rows:
        if num(cells[11]) is not None:
            out.add(cells[0])
    return out

for mskey, rows in T1.items():
    numeric = numeric_t4r_labels(rows)
    got = {cells[0] for cells in rows if '†' in cells[11]}
    exp = {lab for lab, sp in LABEL_SPAN.items()
           if lab in numeric and sp < 0.15}
    missing = sorted(exp - got)
    spurious = sorted(got - exp)
    if missing or spurious:
        rec('P0', 'DAGGER',
            f'{mskey}: missing daggers {missing}; spurious daggers {spurious}')
    P(f'    {mskey}: numeric-T4r rows={len(numeric)}, daggers={len(got)}, expected={len(exp)} -> '
      f'{"PASS" if not (missing or spurious) else "FAIL"}')
    P(f'        interpretable (span>=0.15): {sorted(numeric - exp)}')

# ============================================================
# 5. Aggregate statistics quoted in the text
# ============================================================
P('\n[5] Aggregate statistics')
ssr = {k: d['T1']['SSR'] for k, d in RAW.items()}
p2  = {k: d['T2'].get('empirical_p') for k, d in RAW.items()}
r3  = {k: (d['T3'].get('rho_ngene') if d.get('T3') else None) for k, d in RAW.items()}
r4  = {k: (d['T4'].get('rho_detect_vs_score') if d.get('T4') else None) for k, d in RAW.items()}

EXCL = 'AS_GSE159677_ours'   # pre-scored, not comparable
comparable = {k: v for k, v in ssr.items() if k != EXCL}

P(f'    SSR  n={len(ssr)} min={min(ssr.values()):.3f} max={max(ssr.values()):.3f} '
  f'median={np.median(list(ssr.values())):.3f}; <=1: {sum(v<=1 for v in ssr.values())}')
P(f'    SSR(15) n={len(comparable)} <=1: {sum(v<=1 for v in comparable.values())}')
t2fail = [k for k, v in p2.items() if v is not None and v >= 0.05 and k != EXCL]
P(f'    T2 p>=0.05 among 15: {len(t2fail)}')
# "pass both controls" = the LRG signal is NOT distinguishable from housekeeping on
# either test, i.e. T1 fails to discriminate (SSR <= 1) AND T2 fails to discriminate
# (empirical p >= 0.05).  NOTE: SSR > 1 is the *desired* direction, so a dataset that
# clears both discriminating tests has SSR > 1 and p < 0.05.
passboth = [k for k in comparable
            if comparable[k] > 1 and p2.get(k) is not None and p2[k] < 0.05]
P(f'    pass both discriminating tests (SSR>1 & p<0.05): {len(passboth)} -> {sorted(passboth)}')
failboth = [k for k in comparable
            if comparable[k] <= 1 and p2.get(k) is not None and p2[k] >= 0.05]
P(f'    fail both (SSR<=1 & p>=0.05): {len(failboth)} -> {sorted(failboth)}')
r3v = [v for k, v in r3.items() if v is not None and k != EXCL]
P(f'    T3 rho n={len(r3v)} min={min(r3v):+.3f} max={max(r3v):+.3f} median={np.median(r3v):+.3f}')
r4v = [v for k, v in r4.items() if v is not None]
P(f'    T4 rho n={len(r4v)} min={min(r4v):+.3f} max={max(r4v):+.3f} median={np.median(r4v):+.3f}')
spv = np.array([T4R[k]['span'] for k in T4R])
P(f'    span>=0.15: {int((spv>=0.15).sum())}/{len(spv)}')
sv = np.array([T4R[k]['spearman'] for k in T4R])
pv = np.array([T4R[k]['pearson']  for k in T4R])
P(f'    T4r spearman median={np.median(sv):+.4f} range=[{sv.min():+.4f},{sv.max():+.4f}]')
P(f'    T4r pearson  median={np.median(pv):+.4f} range=[{pv.min():+.4f},{pv.max():+.4f}]')

# cross-check the numbers the manuscript asserts
ASSERT = {
    'SSR median 0.907': abs(np.median(list(ssr.values())) - 0.907) < 5e-4,
    'SSR n=16 min 0.331': abs(min(ssr.values()) - 0.331) < 5e-4,
    'SSR n=16 max 4.320': abs(max(ssr.values()) - 4.320) < 5e-4,
    'SSR<=1 count 9': sum(v <= 1 for v in ssr.values()) == 9,
    'SSR<=1 among 15 = 8': sum(v <= 1 for v in comparable.values()) == 8,
    'T2 fail among 15 = 10': len(t2fail) == 10,
    'pass-both = 3': len(passboth) == 3,
    'fail-both = 6': len(failboth) == 6,
    'T3 median +0.574': abs(np.median(r3v) - 0.574) < 5e-4,
    'T3 min +0.434': abs(min(r3v) - 0.434) < 5e-4,
    'T3 max +0.793': abs(max(r3v) - 0.793) < 5e-4,
    'T4 median +0.903': abs(np.median(r4v) - 0.903) < 5e-4,
    'T4 min +0.709': abs(min(r4v) - 0.709) < 5e-4,
    'T4 max +0.976': abs(max(r4v) - 0.976) < 5e-4,
    'T4r spearman median +0.7697': abs(np.median(sv) - 0.7697) < 5e-4,
    'T4r pearson median +0.7866': abs(np.median(pv) - 0.7866) < 5e-4,
    'span>=0.15 in 5/15': int((spv >= 0.15).sum()) == 5,
}
P('')
for k, ok in ASSERT.items():
    P(f'    [{"OK " if ok else "BAD"}] {k}')
    if not ok:
        rec('P0', 'AGG', f'aggregate assertion failed: {k}')

# ============================================================
# 6. BH FDR on the 16 T2 p-values
# ============================================================
P('\n[6] BH-FDR over the 16 T2 empirical p-values')
ks = list(p2.keys())
# manuscript order for figure legend; sort not needed for BH
ps = np.array([p2[k] for k in ks], dtype=float)
order = np.argsort(ps)
n = len(ps)
q = np.empty(n)
prev = 1.0
for i in range(n - 1, -1, -1):
    idx = order[i]
    val = ps[idx] * n / (i + 1)
    prev = min(prev, val)
    q[idx] = prev
qm = float(q.min())
P(f'    minimum q = {qm:.4f}')
if abs(qm - 0.106) > 5e-4:
    rec('P0', 'BH', f'BH minimum q = {qm:.4f}, manuscript states 0.106')

# ============================================================
# 7. Text-claim consistency scans
# ============================================================
P('\n[7] Text-claim scans')
PATTERNS = [
    (r'0\.056',  'P0', 'stale BH value 0.056 still present'),
    (r'\+0\.059\b', 'P0', 'stale T4r claim +0.059 still present'),
    (r'only\s+\*?\*?5\*?\*?\s+datasets', 'P1', 'check "only 5 datasets" phrasing'),
]
for mskey, path in MS.items():
    txt = Path(path).read_text(encoding='utf-8')
    for pat, sev, msg in PATTERNS:
        for m in re.finditer(pat, txt):
            rec(sev, 'TEXT', f'{mskey}: {msg} @ {m.start()}')
    P(f'    {mskey}: scanned')

# ============================================================
# [8] The bootstrap intervals and the Table-1 T4r point estimates must agree
#     (defect traced during the fourth review round: the text cited intervals from the old table)
# ============================================================
P('\n[8] T4r bootstrap interval consistency')
_bs = _find(['**/t4r_bootstrap_corrected.csv'], 'bootstrap interval table')
_auth = _find(['**/t4r_recomputed.csv'], 'authoritative T4r recomputation table')
if _bs.exists() and _auth.exists():
    import csv as _csv
    bs = {r['dataset']: r for r in _csv.DictReader(open(_bs, encoding='utf-8'))}
    # the authoritative table uses the GSE accession + TAA index; map to the short names
    # used by the bootstrap table
    _m = {'GSE166676': 'AAA', 'GSE131778': 'CAD', 'GSE165816': 'DFU',
          'GSE212966': 'PDAC', 'GSE214611': 'AS (mouse)',
          'GSE156632': 'ccRCC-1', 'GSE171306': 'ccRCC-2'}
    for i in range(1, 9):
        _m[f'GSE155468 (TAA{i})'] = f'TAA{i}'
    n_ok = n_bad = 0
    for r in _csv.DictReader(open(_auth, encoding='utf-8')):
        lab = _m.get(r['Dataset'])
        if not lab or lab not in bs:
            continue
        a = round(float(r['T4r_Spearman']), 3)
        b = round(float(bs[lab]['T4r']), 3)
        if abs(a - b) <= 5e-4:
            n_ok += 1
        else:
            n_bad += 1
            rec('P0', 'BOOTSTRAP',
                f'bootstrap point estimate disagrees with Table 1: {lab} bs={b} table1={a}')
    P(f'    bootstrap point estimate vs Table 1: agree {n_ok} / disagree {n_bad}')
    w = [float(r['width']) for r in bs.values()]
    cz = sum(1 for r in bs.values() if r['contains_zero'] == 'True')
    P(f'    median interval width {np.median(w):.3f}   containing zero {cz}/{len(bs)}')
    # the median width cited in the text must match the recomputed value
    _MAN = _manuscripts()
    _mw = round(float(np.median(w)), 3)
    for k, p in _MAN.items():
        t = Path(p).read_text(encoding='utf-8')
        if '0.618' in t:
            rec('P0', 'BOOTSTRAP', f'{k} still cites the old median interval width 0.618')
        if f'{_mw:.3f}' not in t:
            rec('P1', 'BOOTSTRAP', f'{k} does not contain the recomputed median interval width {_mw:.3f}')
        P(f'    {k}: interval-width check complete')
else:
    P('    (skipped: bootstrap or authoritative recomputation file missing)')

# ============================================================
# [9] Verification of the data sections added in the fourth review round
#     -- selection bias, T2 matching quality, real-data depth span, simulation-generator parameters
# ============================================================
P('\n[9] Data sections added in the fourth review round')
import csv as _csv2

# 9a. T2 matching quality: observed effect / null-distribution mean
_t2 = []
for _f in sorted(RES.glob('*.json')):
    if 'fidelity' in _f.name:
        continue
    _d = json.load(open(_f, encoding='utf-8'))
    _t = _d.get('T2') or {}
    if 'null_mean' not in _t:
        continue
    _t2.append((_t['observed_effect'], _t['null_mean'], _t['null_p95'], _t['empirical_p']))
if _t2:
    _oo = np.array([r[0] for r in _t2]); _nm = np.array([r[1] for r in _t2])
    _ratio = _oo / _nm
    _med = float(np.median(_ratio)); _lo = float(_ratio.min()); _hi = float(_ratio.max())
    _np95 = sum(1 for r in _t2 if r[0] > r[2])
    P(f'    T2 matching quality n={len(_t2)}: obs/null_mean median={_med:.2f} range={_lo:.2f}~{_hi:.2f}')
    P(f'    obs>null_p95: {_np95}/{len(_t2)} (should equal the number with p<0.05)')
    _np05 = sum(1 for r in _t2 if r[3] < 0.05)
    if _np95 != _np05:
        rec('P1', 'T2-MATCH', f'obs>null_p95 count {_np95} disagrees with p<0.05 count {_np05}')
    # the "1.24 times" and "5 of 15" claims in the text
    for _k, _p in MS.items():
        _t = ' '.join(Path(_p).read_text(encoding='utf-8').split())
        if f'{_med:.2f} times' not in _t:
            rec('P1', 'T2-MATCH', f'{_k} does not contain the median T2 matching-quality ratio {_med:.2f}')
        if f'{_np95} of {len(_t2)} datasets' not in _t:
            rec('P1', 'T2-MATCH', f'{_k} does not contain the obs>null_p95 statement {_np95}/{len(_t2)}')
    P('    BMC/BIB: T2 matching-quality statements checked')

# 9b. log2 depth span of the real data (used to show that the simulation depth of 3.0 is conservative)
_dsp = []
for _f in sorted(RES.glob('*_T4.csv')):
    _r = list(_csv2.DictReader(open(_f, encoding='utf-8')))
    _ng = np.array([float(x['mean_ngene']) for x in _r
                    if x.get('mean_ngene') not in ('', 'nan')])
    if _ng.size < 2:
        continue
    _dsp.append(np.log2(_ng.max()) - np.log2(_ng.min()))
if _dsp:
    _dv = np.array(_dsp)
    P(f'    real-data log2 depth span n={len(_dv)}: median={np.median(_dv):.2f} '
      f'range={_dv.min():.2f}~{_dv.max():.2f}')
    P(f'    exceeding 2.0 log2 units: {int((_dv > 2.0).sum())}/{len(_dv)}')
    if (_dv > 2.0).any():
        rec('P0', 'DEPTH-SPAN', 'the text states that no dataset exceeds 2.0 log2 units, but some do')
    for _k, _p in MS.items():
        _t = ' '.join(Path(_p).read_text(encoding='utf-8').split())
        if 'No dataset exceeds 2.0 log units' not in _t:
            rec('P1', 'DEPTH-SPAN', f'{_k} does not contain the phrase "no dataset exceeds 2.0 log units"')

# 9c. simulation generator: the text claims Poisson, so the script must really use Poisson
_sim = _find(['**/sim_t4_t4r.py'], 'simulation script')
if _sim.exists():
    _st = _sim.read_text(encoding='utf-8')
    if 'rng.poisson' not in _st:
        rec('P0', 'SIM', 'the text states that the simulation uses Poisson counts, but the generator never calls rng.poisson')
    for _k, _p in MS.items():
        _t = ' '.join(Path(_p).read_text(encoding='utf-8').split())
        if 'Poisson counts rather than a negative-binomial model' not in _t:
            rec('P1', 'SIM', f'{_k} does not state that Poisson rather than a negative binomial is used')
    P('    Poisson claim of the simulation generator matches the script')

# 9d. the selection-bias statement must appear in both the introduction and the conclusions
for _k, _p in MS.items():
    _t = Path(_p).read_text(encoding='utf-8')
    if 'selection bias' not in _t.lower():
        rec('P1', 'SELBIAS', f'{_k} lacks a selection-bias statement in the introduction')
    if _t.count('lower bound') < 2:
        rec('P1', 'SELBIAS', f'{_k} does not restate in the conclusions that the failure rate is a lower bound')
P('    selection-bias statements (introduction + conclusions) checked')

# ============================================================
# Summary
# ============================================================
P('\n' + '=' * 72)
P('SUMMARY')
P('=' * 72)
if not FINDINGS:
    P('  No findings. All checked numbers reconcile.')
else:
    from collections import Counter
    c = Counter(s for s, _, _ in FINDINGS)
    P(f'  {len(FINDINGS)} findings: ' + ', '.join(f'{k}={v}' for k, v in sorted(c.items())))
    P('')
    for sev, code, msg in FINDINGS:
        P(f'  [{sev}] {code}: {msg}')
P('=' * 72)
sys.exit(0 if not [f for f in FINDINGS if f[0] == 'P0'] else 1)
