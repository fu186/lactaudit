#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
T4 / T4r simulation, with fixed seeds and a machine-readable output record.

WHY THIS SCRIPT EXISTS
----------------------
The manuscript reports a simulation (Figure 5) with specific numbers -- T4 reaching
+1.000 under a pure null, T4 being unchanged by injected signal, and T4r taking
values of -0.073 / +0.076 under a null and +0.147 / +0.280 with signal.  No script
producing those numbers existed anywhere in the repository: the figure hard-coded
them as literal Python lists.  This script supplies the missing generator so the
simulation can be reproduced, checked, and cited.

WHAT IT DOES
------------
Builds synthetic single-cell count matrices in which cell types differ ONLY in
sequencing depth (the null arm), or additionally carry a cell-type-specific signal
in k genes (the signal arm), then computes T4 and T4r exactly as the manuscript
defines them:

    T4  = Spearman( per-type LRG detection rate , per-type mean score )
    T4r = Spearman of the OLS residuals after regressing BOTH quantities on the
          cell type's mean detected-gene count

Scoring uses score_z, the mean of per-gene z-scores -- the estimator the manuscript
audits and the one its simulation claims to cover.

USAGE
-----
    python3 sim_t4_t4r.py              # run and write sim_t4_t4r.json
    python3 sim_t4_t4r.py --seeds 50   # more seeds

OUTPUT
------
sim_t4_t4r.json with, for every (span, k_genes) cell of the design, the per-seed
T4 and T4r values plus their mean/min/max, so every number in Figure 5 can be
traced back to this run.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

# ----------------------------------------------------------------------------
# Fixed design parameters (all recorded in the JSON output)
# ----------------------------------------------------------------------------
SEED_BASE   = 20250915   # every replicate seed is SEED_BASE + replicate index
N_TYPES     = 10         # cell types per synthetic dataset
N_CELLS     = 300        # cells per type
N_BACKGROUND= 2000       # background genes
N_LRG       = 100        # lactylation gene set size
LRG_SCALE   = 0.10       # sets per-cell detection of an LRG gene at ~10%
DEPTH_CV    = 0.25       # log-normal spread of per-cell depth within a type
COUNTS_PER_CELL = 3000   # nominal total counts per cell
SIGNAL_FOLD = 6.0        # over-expression fold in the signal arm

SPANS   = [0.0, 0.5, 1.0, 2.0, 3.0]        # imposed depth span, log units
K_GENES = [0, 10, 20, 30, 40]              # genes over-expressed in one type


def simulate(depth_span_log, n_signal, seed, n_types=N_TYPES, n_cells=N_CELLS):
    """One synthetic dataset.  Returns (cell_type_index, counts, n_lrg_offset)."""
    rng = np.random.default_rng(SEED_BASE + seed)

    # cell types laid out evenly across [0, depth_span_log] in log space
    type_depth = np.linspace(0.0, depth_span_log, n_types)
    depth_mult = np.exp(type_depth)

    ct = np.repeat(np.arange(n_types), n_cells)
    lib = np.maximum(
        (depth_mult[ct] * np.exp(rng.normal(0, DEPTH_CV, ct.size)) * COUNTS_PER_CELL).astype(int),
        200,
    )

    # gene baselines: background broad, LRG sparse so detection varies with depth
    base_bg  = rng.gamma(0.4, 1.0, size=N_BACKGROUND)
    base_lrg = rng.gamma(0.35, 1.1, size=N_LRG)

    fold = np.ones(N_LRG)
    if n_signal:
        fold[:n_signal] = SIGNAL_FOLD

    n_tot = ct.size
    X_bg = rng.poisson(np.outer(lib / lib.mean(), base_bg) * 0.35).astype(np.float32)

    is_type0 = (ct == 0)[:, None].astype(float)
    mod = 1.0 + is_type0 * (fold[None, :].repeat(n_tot, axis=0) - 1.0)
    X_lrg = rng.poisson(np.outer(lib / lib.mean(), base_lrg) * LRG_SCALE * mod).astype(np.float32)

    return ct, np.concatenate([X_bg, X_lrg], axis=1), N_BACKGROUND


def t4_and_t4r(counts, ct, lrg_cols):
    """Return (T4, T4r, detection-rate span) for one dataset."""
    sub = counts[:, lrg_cols]
    score = ((sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-12)).mean(axis=1)

    detect, ngene, mean_score = [], [], []
    for t in np.unique(ct):
        m = ct == t
        detect.append((counts[m][:, lrg_cols] > 0).any(axis=1).mean())
        ngene.append((counts[m] > 0).sum(axis=1).mean())
        mean_score.append(score[m].mean())
    detect, ngene, mean_score = map(np.asarray, (detect, ngene, mean_score))

    r4 = float(spearmanr(detect, mean_score)[0])

    e1 = detect     - np.polyval(np.polyfit(ngene, detect,     1), ngene)
    e2 = mean_score - np.polyval(np.polyfit(ngene, mean_score, 1), ngene)
    r4r = float(spearmanr(e1, e2)[0]) if (e1.std() > 0 and e2.std() > 0) else None

    return r4, r4r, float(detect.max() - detect.min())


def run(n_seeds):
    out = {
        'design': {
            'seed_base': SEED_BASE, 'n_seeds': n_seeds, 'n_types': N_TYPES,
            'n_cells_per_type': N_CELLS, 'n_background_genes': N_BACKGROUND,
            'n_lrg': N_LRG, 'lrg_scale': LRG_SCALE, 'depth_cv': DEPTH_CV,
            'counts_per_cell': COUNTS_PER_CELL, 'signal_fold': SIGNAL_FOLD,
            'spans_log_units': SPANS, 'k_genes': K_GENES,
            'estimator': 'score_z (mean of per-gene z-scores)',
        },
        'null_arm': {},      # depth varied, no biology
        'signal_arm': {},    # span fixed at 3.0, k genes over-expressed
    }

    for span in SPANS:
        t4s, t4rs, sps = [], [], []
        for seed in range(n_seeds):
            ct, X, off = simulate(span, 0, seed)
            r4, r4r, sp = t4_and_t4r(X, ct, list(range(off, off + N_LRG)))
            t4s.append(r4); t4rs.append(r4r); sps.append(sp)
        t4s = np.array(t4s, dtype=float)
        t4rs = np.array([v for v in t4rs if v is not None], dtype=float)
        out['null_arm'][f'{span:g}'] = {
            'T4':  {'values': t4s.tolist(),  'mean': float(t4s.mean()),
                    'min': float(t4s.min()), 'max': float(t4s.max())},
            'T4r': {'values': t4rs.tolist(), 'mean': float(t4rs.mean()),
                    'min': float(t4rs.min()),'max': float(t4rs.max()),
                    'n_valid': int(t4rs.size)},
            'detect_span_mean': float(np.mean(sps)),
        }

    for k in K_GENES:
        t4s, t4rs = [], []
        for seed in range(n_seeds):
            ct, X, off = simulate(3.0, k, seed)
            r4, r4r, sp = t4_and_t4r(X, ct, list(range(off, off + N_LRG)))
            t4s.append(r4); t4rs.append(r4r)
        t4s = np.array(t4s, dtype=float)
        t4rs = np.array([v for v in t4rs if v is not None], dtype=float)
        out['signal_arm'][f'{k}'] = {
            'T4':  {'values': t4s.tolist(),  'mean': float(t4s.mean()),
                    'min': float(t4s.min()), 'max': float(t4s.max())},
            'T4r': {'values': t4rs.tolist(), 'mean': float(t4rs.mean()),
                    'min': float(t4rs.min()),'max': float(t4rs.max()),
                    'n_valid': int(t4rs.size)},
        }

    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--seeds', type=int, default=10, help='replicates per design cell')
    ap.add_argument('--out', default='sim_t4_t4r.json')
    args = ap.parse_args()

    res = run(args.seeds)

    print('=' * 70)
    print('NULL ARM -- depth varied across cell types, no biological difference')
    print('=' * 70)
    print(f"{'span':>5s} | {'T4 mean':>8s} {'T4 range':>18s} | "
          f"{'T4r mean':>9s} {'T4r range':>18s}")
    for span, d in res['null_arm'].items():
        t4, t4r = d['T4'], d['T4r']
        print(f"{span:>5s} | {t4['mean']:8.3f} [{t4['min']:+.3f},{t4['max']:+.3f}] | "
              f"{t4r['mean']:9.3f} [{t4r['min']:+.3f},{t4r['max']:+.3f}]")

    print()
    print('=' * 70)
    print('SIGNAL ARM -- k LRG genes over-expressed 6-fold in one type, span 3.0')
    print('=' * 70)
    print(f"{'k':>4s} | {'T4 mean':>8s} {'T4 range':>18s} | {'T4r mean':>9s}")
    for k, d in res['signal_arm'].items():
        t4, t4r = d['T4'], d['T4r']
        print(f"{k:>4s} | {t4['mean']:8.3f} [{t4['min']:+.3f},{t4['max']:+.3f}] | "
              f"{t4r['mean']:9.3f}")

    Path(args.out).write_text(json.dumps(res, indent=1), encoding='utf-8')
    print()
    print(f'written: {args.out}')


if __name__ == '__main__':
    main()
