"""
Command-line interface.

    lactaudit audit    --matrix X.npz --labels labels.tsv --gene-set genes.txt
    lactaudit compare  --matrix X.npz --labels labels.tsv
    lactaudit demo

``demo`` runs the whole pipeline on a synthetic matrix with a known answer, so
you can verify the installation without downloading anything.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _load_matrix(args):
    """Load a matrix from npz cache or a wide table."""
    import numpy as np

    from . import io as lio
    from .scoring import log_normalize

    path = Path(args.matrix)
    if path.suffix == ".npz" or path.is_dir():
        X, genes, lib = lio.load_cached(path)
        return log_normalize(X, lib), genes, lib
    X, genes, cells = lio.read_wide_table(path, max_cells=args.max_cells)
    Xf, _, stats = lio.qc_filter(X, genes, max_cells=args.max_cells)
    lib = np.asarray(Xf.sum(axis=1)).ravel()
    return log_normalize(Xf, lib), genes, stats


def _load_labels(path, n_cells):
    import numpy as np
    import pandas as pd

    p = Path(path)
    if p.suffix in (".tsv", ".txt", ".csv"):
        sep = "\t" if p.suffix != ".csv" else ","
        df = pd.read_csv(p, sep=sep, header=None)
        labels = df.iloc[:, -1].astype(str).to_numpy()
    else:
        labels = np.array(json.loads(p.read_text()), dtype=str)
    if len(labels) != n_cells:
        sys.exit(f"labels has {len(labels)} entries but matrix has {n_cells} cells")
    return labels


def _cluster(mat, k=10, seed=0):
    """Cheap unsupervised labels so the demo needs no annotation."""
    import numpy as np
    from scipy.cluster.vq import kmeans2

    X = np.asarray(mat)
    # use the most variable genes to define the clusters
    var = X.var(axis=0)
    idx = np.argsort(var)[-min(400, X.shape[1]):]
    sub = X[:, idx]
    sub = (sub - sub.mean(0)) / (sub.std(0) + 1e-9)
    _, lab = kmeans2(sub, k, minit="++", seed=seed)
    return np.array([f"c{i}" for i in lab])


def cmd_audit(args):
    from .audit import audit_gene_set

    mat, genes, _ = _load_matrix(args)
    labels = _load_labels(args.labels, mat.shape[0])

    if args.gene_set:
        gene_set = [g.strip().upper() for g in Path(args.gene_set).read_text().split()
                    if g.strip()]
    else:
        from .genesets import LACTYLATION_SET
        gene_set = LACTYLATION_SET

    res = audit_gene_set(mat, genes, labels, gene_set,
                         n_repeat=args.repeats, seed=args.seed)
    print(res.summary())
    if args.json:
        Path(args.json).write_text(
            json.dumps(res.to_dict(), indent=1, default=float), encoding="utf-8"
        )
        print(f"\nwritten: {args.json}", file=sys.stderr)


def cmd_compare(args):
    """Compare all reference signatures plus the LRG set on the same matrix."""
    import numpy as np

    from .audit import size_matched_null
    from .genesets import (
        LACTYLATION_SET, HOUSEKEEPING_SET, HOUSEKEEPING_POOL, REFERENCE_SIGNATURES,
    )
    from .scoring import score_add_module
    from .stats import kruskal_h, signal_specificity_ratio

    mat, genes, _ = _load_matrix(args)
    labels = _load_labels(args.labels, mat.shape[0])
    index = {g: i for i, g in enumerate(np.asarray(genes))}
    hk = [index[g] for g in HOUSEKEEPING_SET if g in index]
    pool = sorted({index[g] for g in HOUSEKEEPING_POOL if g in index} - set(hk))

    h_hk = kruskal_h(score_add_module(mat, hk), labels)
    rng = np.random.default_rng(args.seed)

    sets = {"LRG": LACTYLATION_SET}
    sets.update(REFERENCE_SIGNATURES)

    print(f"{'set':20s} {'n':>4s} {'H':>9s} {'SSR':>7s}")
    print("-" * 44)
    rows = []
    for name, gs in sets.items():
        idx = [index[g] for g in gs if g in index]
        if len(idx) < 8:
            continue
        h = kruskal_h(score_add_module(mat, idx), labels)
        hs = [kruskal_h(score_add_module(mat, rng.choice(pool, min(len(idx), len(pool)),
                                                          replace=False)), labels)
              for _ in range(20)] if len(pool) >= 8 else []
        h_null = float(np.nanmean(hs)) if hs else float("nan")
        ssr = signal_specificity_ratio(h, h_null) if h_null == h_null else signal_specificity_ratio(h, h_hk)
        rows.append({"set": name, "n": len(idx), "H": h, "SSR": ssr})
        print(f"{name:20s} {len(idx):4d} {h:9.1f} {ssr:7.3f}")

    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=1), encoding="utf-8")


def cmd_demo(args):
    """Synthetic sanity check: LRG must fail, a real marker set must pass."""
    import numpy as np

    from .audit import audit_gene_set
    from .genesets import (
        LACTYLATION_SET, HOUSEKEEPING_SET, HOUSEKEEPING_POOL,
    )

    rng = np.random.default_rng(args.seed)
    n_cells, n_genes = 1200, 2000
    labels = rng.choice([f"ct{i}" for i in range(4)], n_cells)

    # Every cell gets a random library size.  This is the confounder.
    depth = rng.lognormal(mean=8.0, sigma=0.8, size=n_cells)

    genes = np.array([f"G{j}" for j in range(n_genes)], dtype=object)
    counts = np.zeros((n_cells, n_genes), dtype=np.float64)

    # --- 1) Housekeeping block: expressed in EVERY cell at equal level. ---
    # Equal level => no cell-type structure => H should be ~0.  This is the
    # yardstick the audit measures everything else against.
    hk_slots = rng.choice(n_genes, len(HOUSEKEEPING_SET), replace=False)
    genes[hk_slots] = HOUSEKEEPING_SET
    counts[:, hk_slots] = 400

    # --- 2) Background: depth-driven noise for all remaining genes. ---
    taken = set(hk_slots.tolist())
    rest = np.array([j for j in range(n_genes) if j not in taken])
    base = rng.gamma(shape=0.6, scale=1.0, size=rest.size)
    counts[:, rest] = rng.poisson(
        np.clip(np.outer(depth, base / base.sum()) * 1e3, 0, None)
    )

    # --- 3) Lactylation set: broadly expressed (high baseline), no specificity. ---
    lrg_slots = rng.choice(rest, len(LACTYLATION_SET), replace=False)
    genes[lrg_slots] = LACTYLATION_SET
    counts[:, lrg_slots] += rng.poisson(60, (n_cells, len(lrg_slots)))

    # --- 4) Planted markers: each cell type gets its own 15-gene block. ---
    left = np.array([j for j in rest if j not in set(lrg_slots.tolist())])
    pool_slots = left[:len(HOUSEKEEPING_POOL)]
    genes[pool_slots] = HOUSEKEEPING_POOL
    left2 = left[len(HOUSEKEEPING_POOL):]
    marker_genes = []
    for c in range(4):
        block = left2[c * 15:(c + 1) * 15]
        names = [f"MARK{c}_{j}" for j in range(15)]
        genes[block] = names
        marker_genes.extend(names)
        counts[labels == f"ct{c}"][:, block] += 0        # placeholder
        idx = np.where(labels == f"ct{c}")[0]
        counts[np.ix_(idx, block)] += 600

    from .scoring import log_normalize
    mat = log_normalize(counts)

    print("Synthetic matrix: %d cells x %d genes, 4 cell types, "
          "library size log-normal(8, 0.8)" % mat.shape)
    print("  - housekeeping block: expressed equally in all cells (H ~ 0)")
    print("  - lactylation set:    broadly expressed, NO cell-type structure")
    print("  - planted markers:    15 genes per cell type, boosted 600x")
    print("Any apparent cell-type signal outside the planted block is artefact.")
    print()

    marker_set = [g for g in marker_genes]
    for name, gs in [("planted markers", marker_set),
                     ("LRG (depth-driven)", LACTYLATION_SET),
                     ("housekeeping (HK40)", HOUSEKEEPING_SET)]:
        res = audit_gene_set(mat, genes, labels, gs, n_repeat=50, seed=args.seed)
        print(f"{name:22s} SSR = {res.ssr_original:6.3f}   {res.verdict}")
    print("\nExpected: planted >> 1, LRG ~ 1 or below.")


def main(argv=None):
    from . import __version__

    p = argparse.ArgumentParser(
        prog="lactaudit",
        description="Depth-bias audit for single-cell gene-set scoring.",
    )
    p.add_argument("--version", action="version", version=f"lactaudit {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--matrix", required=True,
                        help="npz cache (lactaudit.io.save_cached) or a GEO wide table")
    common.add_argument("--labels", required=True,
                        help="per-cell labels: .tsv/.csv or a JSON array")
    common.add_argument("--max-cells", type=int, default=3000,
                        help="down-sample before analysis (default 3000)")
    common.add_argument("--repeats", type=int, default=200,
                        help="resampling repeats for the null (default 200)")
    common.add_argument("--seed", type=int, default=42)
    common.add_argument("--json", default=None, help="write results here as JSON")

    for name, fn, help_text in [
        ("audit", cmd_audit, "audit one gene set"),
        ("compare", cmd_compare, "compare LRG against reference signatures"),
    ]:
        sp = sub.add_parser(name, parents=[common], help=help_text)
        if name == "audit":
            sp.add_argument("--gene-set", default=None,
                            help="whitespace-separated gene symbols; defaults to the built-in LRG")
        sp.set_defaults(func=fn)

    sp = sub.add_parser("demo", help="run a synthetic sanity check (no data needed)")
    sp.add_argument("--seed", type=int, default=42)
    sp.set_defaults(func=cmd_demo)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
