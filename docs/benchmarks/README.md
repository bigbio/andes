# Benchmarks

How andes is measured, and what the current numbers are. Start with
**[HOWTO.md](HOWTO.md)** — it has the exact commands for all four datasets, the shared
Percolator protocol, how to compute an entrapment FDP correctly, and the mistakes that
have produced wrong published numbers here before.

## Current

| Document | Date | What it establishes |
|---|---|---|
| **[HOWTO.md](HOWTO.md)** | 2026-09-04 | How to run every benchmark: the three standard datasets, the glyco datasets, the FDR protocol, the entrapment arithmetic, and the pitfalls. |
| [Tree-count default + entrapment audit](2026-09-04-tree-count-default-and-entrapment-audit.md) | 2026-09-04 | Astral **1.64x faster** at unchanged identifications (400 s → 244 s, −8 PSMs) on merged `main`. Audits the entrapment metric: Astral has no entrapment component, UPS1's is not 1:1. |
| [Configuration matrix](2026-08-28-config-matrix.md) | 2026-08-28 | Per-configuration counts, and the corrected UPS1 entrapment FDP (~2.5% at a nominal 1%). |
| [andes vs Comet refresh](2026-08-23-andes-vs-comet-refresh.md) | 2026-08-23 | Dataset, file and database provenance for the three standard sets. |

## Glyco

| Document | Date | What it establishes |
|---|---|---|
| [Glyco benchmark summary](2026-08-27-benchmark-summary.md) | 2026-08-27 | Where andes stands on intact N-glycopeptides, and where it loses. |
| [Glyco algorithm conclusions](glyco-algorithm-conclusions.md) | 2026-09-03 | What has been tried on the glyco path and what the measurements refuted. |

The harness that produces these is [`../../benchmarks/glyco/`](../../benchmarks/glyco/);
its README documents each script and the two rules they encode (pool before Percolator;
never ship on yield alone).

## Historical

| Document | Date | Note |
|---|---|---|
| [Public benchmark](2026-06-15-public-benchmark.md) | 2026-06-15 | The three-engine comparison as it stood in June. Predates the 2026-09-04 speedup and the entrapment-metric correction — read those first. |
| [Own-geometry A/B](2026-06-26-owngeometry-ab.md) | 2026-06-26 | Sole record of the geometry retrain that shipped. |

Four early multi-engine validation reports (2026-06-01 to 2026-06-04) were removed on
2026-09-04: they were already banner-marked superseded, used anonymised engine names that
made them uninformative outside the project, and are recoverable from git history.

Per-engine configuration files are in [`configs/`](configs/); driver scripts in
[`scripts/`](scripts/).
