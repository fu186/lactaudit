#!/usr/bin/env python3
"""
Compare the lactylation set against non-metabolic reference signatures.

This is the analysis that answers "is the depth failure specific to metabolism?"
If a hypoxia or stemness signature fails just as badly, the answer is no -- the
problem lies with the estimator, not with the biology.

Run::

    python examples/compare_signatures.py
"""

import numpy as np

import lactaudit as la


def build_synthetic(n_cells=900, n_genes=1500, seed=1):
    rng = np.random.default_rng(seed)
    labels = rng.choice([f"type{i}" for i in range(5)], n_cells)
    depth = rng.lognormal(8.2, 1.0, n_cells)

    genes = np.array([f"G{j}" for j in range(n_genes)], dtype=object)
    counts = np.zeros((n_cells, n_genes), dtype=np.float64)

    hk_slots = rng.choice(n_genes, len(la.HOUSEKEEPING_SET), replace=False)
    genes[hk_slots] = la.HOUSEKEEPING_SET
    counts[:, hk_slots] = rng.poisson(400, (n_cells, len(hk_slots)))

    rest = np.array([j for j in range(n_genes) if j not in set(hk_slots.tolist())])
    base = rng.gamma(0.6, 1.0, rest.size)
    counts[:, rest] = rng.poisson(
        np.clip(np.outer(depth, base / base.sum()) * 1e3, 0, None)
    )

    # Give every reference signature a broad, non-specific expression profile.
    used = set()
    all_sets = {"LRG": la.LACTYLATION_SET, **la.REFERENCE_SIGNATURES,
                "HK": la.HOUSEKEEPING_POOL}
    for name, gs in all_sets.items():
        free = np.array([j for j in rest if j not in used])
        slots = free[:len(gs)]
        used.update(slots.tolist())
        for j, g in zip(slots, gs):
            genes[j] = g
        counts[:, slots] += rng.poisson(45, (n_cells, len(slots)))

    return la.log_normalize(counts), genes, labels


def main():
    mat, genes, labels = build_synthetic()
    index = {g: i for i, g in enumerate(genes)}
    print(f"matrix: {mat.shape[0]} cells x {mat.shape[1]} genes")
    print(f"cell types: {sorted(set(labels))}\n")

    hk = [index[g] for g in la.HOUSEKEEPING_SET if g in index]
    h_hk = la.kruskal_h(la.score_add_module(mat, hk), labels)
    print(f"housekeeping reference: H = {h_hk:.2f}\n")

    print(f"{'signature':22s} {'n':>4s} {'H':>9s} {'SSR':>8s}  verdict")
    print("-" * 78)
    for name, gs in [("Lactylation (LRG)", la.LACTYLATION_SET), *la.REFERENCE_SIGNATURES.items()]:
        gi = [index[g] for g in gs if g in index]
        if len(gi) < 8:
            print(f"{name:22s} {len(gi):4d}  (too few genes present)")
            continue
        h = la.kruskal_h(la.score_add_module(mat, gi), labels)
        ssr = la.signal_specificity_ratio(h, h_hk)
        flag = "PASS" if ssr >= 1 else ("marginal" if ssr >= 0.8 else "FAIL")
        print(f"{name:22s} {len(gi):4d} {h:9.1f} {ssr:8.3f}  {flag}")

    print("\nIf every broadly expressed signature lands below 1,")
    print("the failure belongs to the estimator, not to the gene set.")


if __name__ == "__main__":
    main()
