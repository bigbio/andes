# Benchmarks

Reproducible engine comparisons for andes against the open-source field. Each
benchmark documents the datasets, the exact per-engine parameters, and the FDR
methodology so the numbers can be regenerated.

| Benchmark | Date | Engines | Datasets |
|---|---|---|---|
| [Public benchmark](2026-06-15-public-benchmark.md) | 2026-06-15 | andes (top-1 + `--chimeric`), Java MS-GF+, Sage, Comet, ProSE | Astral (HCD high-res), TMT a05058 (CID low-res), UPS1/PXD001819 (CID low-res) |

Every engine is re-scored through one uniform Percolator (3.7.1, `--seed 42 -Y`),
and FDR honesty is verified independently with a 1:1 entrapment search (true FDP
≈ 1% at the nominal 1% q-value). Per-engine configuration files live under
[`configs/`](configs/); reproducibility scripts under [`scripts/`](scripts/).
