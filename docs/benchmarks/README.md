# Benchmarks

Reproducible engine comparisons for andes. Each benchmark documents the
datasets, the exact per-engine parameters, and the FDR methodology so the
numbers can be regenerated.

| Benchmark | Date | Engines | Datasets |
|---|---|---|---|
| [4-engine native-format](2026-06-01-4engine-native-format.md) | 2026-06-01 | Java MS-GF+, a comparison search engine, the reference engine, andes | Orbitrap Astral (`.raw`), Bruker timsTOF (`.d`) |
| [PXD016999 TMT tissue (4-engine)](2026-06-03-pxd016999-tmt-4engine.md) | 2026-06-03 | Java MS-GF+, the reference engine, a comparison search engine, andes | PXD016999 human tissue TMT, ion-trap CID-MS2 (`.raw`) |
| [Astral (7-engine, uniform Percolator)](2026-06-04-astral-7engine.md) | 2026-06-04 | andes (top-1 + `--chimeric`), the reference engine, a comparison search engine, a comparison engine, Java MS-GF+, ProSE | Orbitrap Astral (`.raw` / mzML) |
| [TMT a05058 (6-engine, low-res CID)](2026-06-04-tmt-a05058-6engine.md) | 2026-06-04 | andes (top-1 + `--chimeric`), the reference engine, Java MS-GF+, a comparison engine, a comparison search engine, ProSE | a05058 low-res CID/TMT (mzML) |

Configuration files used by each benchmark live under [`configs/`](configs/).
