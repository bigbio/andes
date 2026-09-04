# Benchmarks

Reproducible engine comparisons for andes against the open-source field. Each
benchmark documents the datasets, the exact per-engine parameters, and the FDR
methodology so the numbers can be regenerated.

**Canonical public benchmark — andes vs Java MS-GF+ vs a comparison engine:**

| Benchmark | Date | Engines | Datasets |
|---|---|---|---|
| [Tree-count default + entrapment audit](2026-09-04-tree-count-default-and-entrapment-audit.md) | 2026-09-04 | andes (merged `main`) | Astral, UPS1 — measured 1.64x speedup at unchanged IDs, plus an audit of the entrapment metric |
| [Public benchmark](2026-06-15-public-benchmark.md) | 2026-06-15 | andes (top-1 + `--chimeric`), Java MS-GF+, a comparison engine | Astral (HCD high-res), TMT a05058 (CID low-res), UPS1/PXD001819 (CID low-res) |

Every engine is re-scored through one uniform Percolator (3.7.1, `--seed 42 -Y`).
FDR honesty is checked against entrapment databases **where one exists** — but not
uniformly, and not always at the ratio once assumed. A 2026-09-04 audit
([`2026-09-04-tree-count-default-and-entrapment-audit.md`](2026-09-04-tree-count-default-and-entrapment-audit.md))
found the Astral benchmark database carries no entrapment component at all, and the UPS1
one is 6,733:4,531 rather than 1:1 — so its true FDP at a nominal 1% is 2.6–3.6%, not ~1%.
Size `T/E` per database; never assume it is 1. Per-engine configuration files live under
[`configs/`](configs/); reproducibility scripts under [`scripts/`](scripts/).

<details><summary>Internal / superseded development reports</summary>

Early multi-engine validation runs (include the reference engine / a comparison search engine / ProSE for development
context), kept for reproducibility but superseded by the public benchmark above:

- [2026-06-04 · Astral 7-engine](2026-06-04-astral-7engine.md)
- [2026-06-04 · TMT a05058 6-engine](2026-06-04-tmt-a05058-6engine.md)
- [2026-06-03 · PXD016999 TMT 4-engine](2026-06-03-pxd016999-tmt-4engine.md)
- [2026-06-01 · 4-engine native-format](2026-06-01-4engine-native-format.md)

</details>
