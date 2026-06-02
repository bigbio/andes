# Benchmarks

Reproducible engine comparisons for msgf-rust. Each benchmark documents the
datasets, the exact per-engine parameters, and the FDR methodology so the
numbers can be regenerated.

| Benchmark | Date | Engines | Datasets |
|---|---|---|---|
| [4-engine native-format](2026-06-01-4engine-native-format.md) | 2026-06-01 | Java MS-GF+, Sage, MSFragger, msgf-rust | Orbitrap Astral (`.raw`), Bruker timsTOF (`.d`) |

Configuration files used by each benchmark live under [`configs/`](configs/).
