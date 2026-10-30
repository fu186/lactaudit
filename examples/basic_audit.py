#!/usr/bin/env python3
"""
Basic usage example: audit the lactylation set on a real GEO dataset.

Run from the repository root::

    python examples/basic_audit.py --data /path/to/GSE155468

``--data`` should point at a directory containing the GEO supplementary files
for one sample (e.g. ``GSM4704935_TAA2.txt.gz``).  If omitted, the script falls
back to a cached matrix at ``--cache``.

What it prints
--------------
    H(target) and H(HK40)   the Kruskal-Wallis statistics
    SSR                     the ratio; 1.0 is the decision boundary
    depth coupling          Spearman rho between score and genes detected,
                            for each of the four scorers

Expected pattern on real data: SSR below 1, and depth coupling positive for
the summed scorers, reduced but not eliminated by variance stabilisation.
"""

import argparse
from pathlib import Path

import numpy as np

import lactaudit as la


def load(data_dir: Path):
    """Read the first GEO supplementary matrix found in ``data_dir``."""
    patterns = ("*.txt.gz", "*.tsv.gz", "*.csv.gz", "*.txt", "*.tsv", "*.csv", "*.mtx.gz")
    for pat in patterns:
        hits = sorted(data_dir.glob(pat))
        if hits:
            print(f"reading {hits[0].name}")
            return la.read_wide_table(hits[0])
    raise SystemExit(f"no matrix file found in {data_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=None,
                    help="directory with GEO supplementary files")
    ap.add_argument("--cache", default=None,
                    help="alternative: an npz written by lactaudit.io.save_cached")
    ap.add_argument("--max-cells", type=int, default=3000)
    args = ap.parse_args()

    if args.cache:
        X, genes, lib = la.load_cached(Path(args.cache))
    elif args.data:
        X, genes, cells = load(Path(args.data))
        X, mask, stats = la.qc_filter(X, genes, max_cells=args.max_cells)
        lib = np.asarray(X.sum(axis=1)).ravel()
        print(f"QC: {stats['n_total']} -> {stats['n_used']} cells")
    else:
        raise SystemExit("supply --data or --cache")

    mat = la.log_normalize(X, libsize=lib)
    n_det = (mat > 0).sum(axis=1)

    # Cell types: here we cluster, because we do not have the paper's labels.
    labels = la.cluster_mat(mat, k=10)
    print(f"matrix: {mat.shape[0]} cells x {mat.shape[1]} genes")
    print(f"genes detected per cell: median {np.median(n_det):.0f}")
    print(f"cell types (k=10 clusters): {sorted(set(labels))}\n")

    res = la.audit_gene_set(mat, genes, labels, la.LACTYLATION_SET, n_repeat=200)
    print(res.summary())

    print("\n--- size-matched regimes")
    print(f"  original (all LRG vs HK40)  SSR = {res.ssr_original:.3f}")
    print(f"  shrink target to 40         SSR = {res.ssr_shrink:.3f}")
    print(f"  grow control to LRG size    SSR = {res.ssr_grow:.3f}")


if __name__ == "__main__":
    main()
