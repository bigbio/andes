# MS-GF+ Shipped Work — Short Retrospective

Condensed history of recent iterations. For long-form, see `docs/changelog.md` (user-facing) or `~/.claude/plans/<topic>/` (archived).

## Current state (dev-tip @ `2216bbb`)

| Dataset | Wall (s) | RSS | 1 % FDR PSMs |
|---|---:|---:|---:|
| PXD001819 (Velos, 4 MB) | 105 | 2.2 GB | 15 157 |
| Astral (ProteoBench, 32 MB) | ~620 | 7.6 GB | 35 627 |
| TMT PXD007683 (Lumos, 17 MB) | 321 | 3.7 GB | 10 176 |

Output is `.pin` only (mzIdentML removed). Sensitivity leads Sage at 1 % FDR on every dataset; **speed/RAM gap on Astral (~7.9× behind Sage on wall) is the open frontier.**

## Iteration log

**PR #15-#20 + PR #22 — primitives optimization (Achievements A + B).** GF inner loop ported to primitive arrays. Pin feature additions (longest_b/y). Two-pass precursor mass calibration. `Hashtable`→`HashMap` in `NewRankScorer` killed ~43 % of CPU previously lost to synchronized lookup contention. **Impact:** +254 / +913 / +1 375 PSMs at 1 % FDR (PXD001819 / Astral / TMT).

**PR #23 — speed-v2 cleanup + output consolidation** (`feat/msgfplus-speed-v2`). mzIdentML reader/writer removed; `.pin` is default and only modern format. Pin ion-series run-length features (`longest_b`, `longest_y`, `longest_y_pct`). Tighter `CandidatePeptideGrid` allocation, `Partition.hashCode` cache.

**PR #24 — Astral OOM fix + BuildSA scaling** (`feature/improve-mzid-suffix-big-fasta`). mzML parser MS-level preload filter (cache MS2 only by default) + bounded cache: solves Astral OOM at 8 GB Xmx. BuildSA parallel per-thread bucket sort + merge, no `Suffix[]` boxing, `.cseq` `readFully`. Defer per-task `ScoredSpectraMap` construction to worker thread. Finished removing `jmzidml` dep. *Caveat:* the MS-level filter excludes MS1 — future MS1-aware work must widen filter or add an MS1 accessor.

**PR #25 — search-sync-cleanup + parameter-modernization** (`perf/search-sync-cleanup`). Per-task wall stats + tail-imbalance summary; per-task result buffers (drops shared `synchronizedList`); opt-in ForkJoinPool path. Dropped redundant `synchronized` wrappers in `DBScanner` and `ScoredSpectraMap`. CLI rewritten on picocli (`MSGFPlusOptions`); typed converters/enums for tolerance, int-ranges, `-outputFormat`, `-precursorCal`; `edu.ucsd.msjava.params` hierarchy deleted; `ParamManager` retired from the hot path. Audit pass dropped ~2 074 LOC.

## Abandoned

**Fragment-index (abandoned 2026-04-20).** Sage-style inverted index as Tier-1 candidate generator. Failed all three gates: 1.78× *slower* on PXD001819, OOM on Astral, recall 95.3 % vs ≥ 99.5 % target. Five follow-up speed ideas distilled (graph-skeleton caching, adaptive precursor tolerance, Vector API, parallelism ceiling, SpecEValue caching) — current `feat/astral-speed-improvements` draws from these. Post-mortem: `~/.claude/plans/msgfplus-fragment-index/ABANDONED-2026-04-20.md`.

## Active

- [`astral-speed-improvements.md`](astral-speed-improvements.md) — gate B (1.3-1.5× Astral wall, no PSM regression). TMT-as-inner-loop, Astral-as-phase-gate.
