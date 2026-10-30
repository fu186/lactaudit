#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Corrected bootstrap: recompute the intervals under the authoritative T4r
definition of §2.5 of the manuscript.

Background (found while tracing back during the fourth review round):
  The archived `Table_S5b_T4r_bootstrap_intervals.csv` was produced under the
  **old T4r definition / old sample**. Its point estimates do not match the
  current authoritative T4r column:
      CAD   bootstrap 0.915  vs  authoritative 0.8667   x
      ccRCC bootstrap 0.418  vs  authoritative 0.4182   ok
      TAA1  bootstrap 0.212  vs  authoritative 0.8424   x
      TAA3  bootstrap 0.806  vs  authoritative 0.7212   x
      TAA8  bootstrap 0.309  vs  authoritative 0.2485   x
  Yet the manuscript text still cites this stale set of intervals
  ("[-0.152, +0.843] for TAA1", etc.), which contradicts its own Table 1.

  Note: the point estimates in the text (CAD 0.867, TAA1 0.842, TAA8 0.249)
  agree with the authoritative column; only the intervals come from the old
  table -- a case of "new point estimates paired with old intervals".

What this script does:
  Cell-type-level bootstrap (resample the 10-11 cell-type rows with replacement
  B times); after each resample, recompute T4r under the §2.5 definition
  (OLS residuals on mean_ngene -> Spearman).
  This is the **only bootstrap executable with the available data** (there is no
  raw expression matrix, so cell-level resampling is impossible), and it targets
  exactly the dominant error source of T4r -- the use of only 10 aggregated points.

Honest caveat: cell-type-level resampling captures "cell-type sampling variability";
at B=2000 the percentile intervals are fairly stable, but this is still not the full
uncertainty (annotation uncertainty cannot be estimated from it).
"""
import csv
import glob
import json
import math
import os
import sys
import statistics as st

import numpy as np

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '模拟'))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '06_verification'))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / '模拟'))
import audit_paths as AP            # noqa: E402

HERE = Path(__file__).resolve().parent
RESULTS = str(AP.RES)
OUT_JSON = str(HERE / 't4r_bootstrap_corrected.json')
OUT_CSV = str(HERE / 't4r_bootstrap_corrected.csv')

B = 2000
SEED = 20250916

# dataset label mapping
LABELS = {
    'AAA_GSE166676': 'AAA', 'CAD_GSE131778': 'CAD', 'DFU_GSE165816': 'DFU',
    'PDAC_GSE212966': 'PDAC', 'AS_mouse_GSE214611': 'AS (mouse)',
    'ccRCC_GSE156632': 'ccRCC-1', 'ccRCC2_GSE171306': 'ccRCC-2',
    'TAA_GSE155468': 'TAA1', 'TAA_TAA2_GSE155468': 'TAA2',
    'TAA_TAA3_GSE155468': 'TAA3', 'TAA_TAA4_GSE155468': 'TAA4',
    'TAA_TAA5_GSE155468': 'TAA5', 'TAA_TAA6_GSE155468': 'TAA6',
    'TAA_TAA7_GSE155468': 'TAA7', 'TAA_TAA8_GSE155468': 'TAA8',
}


def _rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = _rank(list(x)), _rank(list(y))
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    if den == 0:
        return None
    return num / den


def t4r_of(ng, dr, sc):
    """Per §2.5: residualise on mean_ngene with OLS, then compute Spearman."""
    ng, dr, sc = np.asarray(ng, float), np.asarray(dr, float), np.asarray(sc, float)
    if len(ng) < 3:
        return None
    A = np.vstack([np.ones(len(ng)), ng]).T
    rd = dr - A @ np.linalg.lstsq(A, dr, rcond=None)[0]
    rs = sc - A @ np.linalg.lstsq(A, sc, rcond=None)[0]
    if np.std(rd) == 0 or np.std(rs) == 0:
        return None
    return spearman(rd, rs)


def main():
    rng = np.random.default_rng(SEED)
    out = []
    for f in sorted(glob.glob(f'{RESULTS}/*_T4.csv')):
        stem = os.path.basename(f)[:-7]
        label = LABELS.get(stem, stem)
        rows = list(csv.DictReader(open(f, encoding='utf-8')))
        ng = [float(r['mean_ngene']) for r in rows]
        dr = [float(r['detect_rate_lrg']) for r in rows]
        sc = [float(r['lact_score']) for r in rows]
        n = len(ng)
        if n < 3:
            continue

        point = t4r_of(ng, dr, sc)
        span = max(dr) - min(dr)

        boots = []
        idx = np.arange(n)
        for _ in range(B):
            take = rng.choice(idx, size=n, replace=True)
            if len(set(take.tolist())) < 3:
                continue
            v = t4r_of([ng[i] for i in take], [dr[i] for i in take], [sc[i] for i in take])
            if v is not None and not math.isnan(v):
                boots.append(v)

        if not boots:
            continue
        lo = float(np.percentile(boots, 2.5))
        hi = float(np.percentile(boots, 97.5))
        out.append({
            'dataset': label,
            'n_types': n,
            'T4r': round(point, 4) if point is not None else None,
            'ci_lo': round(lo, 4),
            'ci_hi': round(hi, 4),
            'width': round(hi - lo, 4),
            'contains_zero': bool(lo <= 0 <= hi),
            'detect_span': round(span, 4),
            'interpretable': bool(span >= 0.15 and n >= 6),
        })

    out.sort(key=lambda r: (r['T4r'] is None, r['T4r'] or 0))
    widths = [r['width'] for r in out]
    cz = [r for r in out if r['contains_zero']]
    interp = [r for r in out if r['interpretable']]

    with open(OUT_JSON, 'w', encoding='utf-8') as fh:
        json.dump({'records': out,
                   'summary': {
                       'n': len(out),
                       'B': B,
                       'seed': SEED,
                       'median_width': round(st.median(widths), 4),
                       'range_width': [round(min(widths), 4), round(max(widths), 4)],
                       'n_contains_zero': len(cz),
                       'n_interpretable': len(interp),
                       'n_interpretable_contains_zero': sum(1 for r in interp if r['contains_zero']),
                   }}, fh, ensure_ascii=False, indent=1)

    with open(OUT_CSV, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    print(f'B={B}  seed={SEED}  datasets n={len(out)}')
    hdr = ('dataset', 'T4r', '95% CI', 'width', 'zero', 'span', 'interp')
    print(f'{hdr[0]:<12}{hdr[1]:>8}{hdr[2]:>20}{hdr[3]:>8}{hdr[4]:>6}{hdr[5]:>8}{hdr[6]:>8}')
    for r in out:
        ci = '[{:+.3f}, {:+.3f}]'.format(r['ci_lo'], r['ci_hi'])
        print(f"{r['dataset']:<12}{r['T4r']:>+8.3f}{ci:>20}{r['width']:>8.3f}"
              f"{'YES' if r['contains_zero'] else 'no':>6}"
              f"{r['detect_span']:>8.3f}"
              f"{'YES' if r['interpretable'] else 'no':>8}")
    print()
    print(f'Median width = {st.median(widths):.3f} (range {min(widths):.3f}-{max(widths):.3f})')
    print(f'Intervals containing zero = {len(cz)}/{len(out)}')
    print(f'Interpretable (span>=0.15) = {len(interp)}, of which {sum(1 for r in interp if r["contains_zero"])} contain zero')


if __name__ == '__main__':
    main()
