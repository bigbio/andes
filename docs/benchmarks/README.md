# Benchmarks

Reproducible engine comparisons for andes. Each benchmark documents the
datasets, the exact per-engine parameters, and the FDR methodology so the
numbers can be regenerated.

| Benchmark | Date | Engines | Datasets |
|---|---|---|---|
| [4-engine native-format](2026-06-01-4engine-native-format.md) | 2026-06-01 | Java MS-GF+, Sage, MSFragger, andes | Orbitrap Astral (`.raw`), Bruker timsTOF (`.d`) |
| [PXD016999 TMT tissue (4-engine)](2026-06-03-pxd016999-tmt-4engine.md) | 2026-06-03 | Java MS-GF+, MSFragger, Sage, andes | PXD016999 human tissue TMT, ion-trap CID-MS2 (`.raw`) |

Configuration files used by each benchmark live under [`configs/`](configs/).
