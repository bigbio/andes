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

**Fragment-index (abandoned 2026-04-20).** Sage-style inverted index as Tier-1 candidate generator. Failed all three gates: 1.78× *slower* on PXD001819, OOM on Astral, recall 95.3 % vs ≥ 99.5 % target. Five follow-up speed ideas distilled (graph-skeleton caching, adaptive precursor tolerance, Vector API, parallelism ceiling, SpecEValue caching). Post-mortem: `~/.claude/plans/msgfplus-fragment-index/ABANDONED-2026-04-20.md`.

**Phase A — deisotoping + peak cap + GF candidate cap + scorer hot-path opt (attempted, reverted 2026-04-28).** Three independent optimization angles tried on `feat/astral-speed-improvements`. None moved Astral wall above run-to-run noise (six measured variants vs OFF baseline 690 s; best Phase A variant was 693 s). TMT showed 1.41× wall but with −0.25 % target / −4.6 % decoy drift — not a clean win. JFR-identified `HashMap.getNode` hot spot did not translate to wall improvement after elimination (JIT already optimizes the path). Branch reset to `eee9fa6`. Retrospective with measurements + lessons + what's untried: [`astral-phase-a-retrospective.md`](astral-phase-a-retrospective.md). Reverted code recoverable via `git show 5cdd21e` (walks back through 11 commits).

**Phase E — parallelism / ForkJoin smart-default (attempted, reverted 2026-04-28; final disproof 2026-04-29).** Initial measurements suggested default `ThreadPoolExecutor` anti-scaled past 6 threads on Astral (4t=690 s, 8t=884 s, +28 %), and the opt-in ForkJoin path (`-Dmsgfplus.useForkJoin=true`) gave 521 s at 8t (1.32×). Implemented auto-default `numThreads >= 8 → ForkJoin`; reverted same day when confirmation runs showed ~30 % wall variance on the same JAR. Multi-run replication on quieter machine (2026-04-29) proved both initial findings were noise: 4t=963 s, 8t=918 s, 8t-FJ=979 s — all within 6.5 % of each other, with 8t-default *faster* than 4t-default. **The yesterday-morning 690 s baseline and 521 s ForkJoin were outliers, not signal.** No Phase E shippable change exists. Retrospective has the full corrected Phase E section.

## Active

**Phase B (calibrated precursor-window tightening) — shipped on `feat/astral-speed-improvements` 2026-04-29.** Four enabling commits:

- `781738e` opt-in `PhaseBTelemetry` counter (pairing fan-out verification via `-Dmsgfplus.phaseBTelemetry=true`)
- `05ec066` calibrator pre-pass uses iso=[0,0] (rejects isotope-error contamination); +50 ppm outlier filter
- `7c027f8` Phase B formula constants exposed as system properties (`-Dmsgfplus.tighteningSigmaMultiplier=<float>` etc.)
- `aac389c` stratify residuals by spec_eValue, keep top MIN_CONFIDENT_PSMS — drops Astral sigma 4× (3.99 → 0.99 ppm)

Astral measurements on `pride-linux-vm.ebi.ac.uk` (5 OFF + 3 AUTO replicates):

| Workload | Window | Sigma | Tightened | Wall Δ | Targets Δ | T/D Δ |
|---|---:|---:|---:|---:|---:|---:|
| **Astral** (ProteoBench Module 8) | 10 ppm | 0.99 ppm | → 3.48 ppm | **−10.4 %** | +0.11 % | +3.6 % ✓ |
| **TMT** (PXD007683, Lumos) | 20 ppm | 2.05 ppm | → 6.67 ppm | **−18.0 %** | −2.05 % ⚠ | +1.3 % ✓ |
| **PXD001819** (Velos) | 5 ppm | 2.15 ppm | safely no-tighten | ~0 % | +0.17 % | +0.5 % ✓ |

Pattern: Phase B wins when calibrated sigma is materially smaller than the user's precursor window; safely no-ops otherwise. TMT's −2.05 % target drift is a known yellow flag — Lumos's wider residual tails are not fully covered by 3-σ. Mitigations for Phase B's broader rollout: instrument-aware k (e.g., k=4 for Lumos) or stricter stratification (top-100 by spec_eValue). T/D ratio still favors target on all three workloads.

OFF-mode (`-precursorCal off`) is bit-identical to dev-tip. Tunable per-workload via `-Dmsgfplus.tighteningSigmaMultiplier=<float>` (default 3.0; k=2 was tested as falsification before stratification fix).

- [`astral-next-experiments.md`](astral-next-experiments.md) — Phase B status notes; Experiment 2 (mass-interval pruning) still untried.
- [`astral-speed-5x-roadmap.md`](astral-speed-5x-roadmap.md) — long-horizon roadmap; Phase B now shipped.
