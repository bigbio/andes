# Phase A — Retrospective (attempted, reverted 2026-04-28)

**Attempt date:** 2026-04-27 to 2026-04-28
**Branch:** `feat/astral-speed-improvements` (reset to `eee9fa6` = consolidated 5× roadmap; Phase A code reverted)
**Decision:** Reverted. None of three independent optimization angles moved Astral wall above noise. TMT/Lumos win was real but not clean enough to justify shipping the surface area.

This retrospective is the artifact future agents should read before re-attempting Astral speed work.

## What was attempted

Three independent angles, all with bit-identical OFF-mode behaviour, gated by Astral measurement:

### Angle 1 — Phase A: in-engine MS2 deisotoping + dense-peak retention cap
- New classes: `Deisotoper`, `Spectrum.deisotope(ppm, maxCharge)`, `Spectrum.capByIntensity(topN)`.
- New CLI: `-deisotopeMS2 on|off`, `-maxPeaksPerSpectrum N`.
- Wired into `ScoredSpectraMap.preProcessSpectra` (main pass only, NOT `MassCalibrator` pre-pass — defended by `Spectrum.isDeisotoped()` idempotence guard).
- Hardcoded 20 ppm spacing tolerance, max charge 6.

### Angle 2 — Iteration 0.5: Tier-1.5 GF candidate cap
- Static field `DBScanner.NUM_CANDIDATES_FOR_GF`, set via `-Dmsgfplus.numCandidatesForGF=N` system property (default 0 = unlimited).
- After cheap-score collection, sort `matchQueue` by score descending, truncate to top-N, then proceed to GF.
- Idea: tighter `minScore` → tighter `setUpScoreThreshold` → smaller GF DP table.

### Angle 3 — NewRankScorer hot-path optimization
- Profile-driven: JFR showed `NewRankScorer.getIonExistenceScore` dispatching `HashMap.get` was ~14 % of Astral CPU.
- Fix: pre-resolve `Float[] ionExistenceProb` per spectrum in `DBScanScorer` and `NewScoredSpectrum` constructors. New overload `getIonExistenceScore(Float[], int, float)` skips the per-edge HashMap lookup.

### Also added (and retained in the abandoned attempt)
- `SearchTelemetry` thread-safe counter class with `-Dmsgfplus.telemetry=true` toggle and `<output>.telemetry.tsv` emission. Used to measure per-spectrum candidates and cheap-score calls. Built into the iteration but never made it past the reset since it was useful only for the killed measurement campaign.

## Astral measurements (clean idle box, 4 threads, 8 GB Xmx, dev-tip @ `2216bbb`)

All runs used the same JAR build per angle, same machine state, same FASTA, same mzML.

| Run | Wall (s) | Peak RSS (MB) | Native targets | Native decoys | Δ wall vs OFF |
|---|---:|---:|---:|---:|---:|
| **OFF (baseline)** | **690.1** | **7 789** | **89 360** | **46 913** | — |
| Phase A (deisotope + cap=200) | 693.4 | 7 088 | 86 134 | 48 497 | +0.5 % |
| Deisotope only (no cap) | 741.3 | 6 832 | 88 941 | 50 819 | +7.4 % |
| GF candidate cap=10 | 714.5 | 6 924 | 89 360 | 46 913 | +3.5 % |
| GF candidate cap=5 | 733.7 | 7 408 | 89 338 | 46 913 | +6.3 % |
| Scorer-opt (cache `ionExistenceProb`) | 719.3 | 6 312 | 89 360 | 46 913 | +4.2 % |

**No variant beats OFF on wall by more than run-to-run noise (~3-5 %).** Three variants (GF cap=10, GF cap=5, scorer-opt) preserve native target/decoy counts bit-identically; Phase A and deisotope-only drift on counts.

JFR profile of Astral OFF (600 s run, 116 K samples) is at `~/work/msgfplus-workspace/benchmark/results/phaseA/astral_off.jfr`.

## TMT measurements (PXD007683, same machine state)

| Run | Wall (s) | Peak RSS (MB) | Native targets | Native decoys | Δ wall vs OFF |
|---|---:|---:|---:|---:|---:|
| OFF | 330.7 | 2 762 | 28 790 | 14 768 | — |
| Phase A (deisotope + cap=200) | 234.5 | 2 820 | 28 719 | 14 081 | **−29 %** |

TMT did show a 1.41× wall reduction, but with **−0.25 % targets and −4.6 % decoys**. The decoy-pool contraction is the bigger concern: it changes Percolator's FDR-calibration shape. A "1.41× faster" claim that comes with non-trivial recall drift is not a clean win.

## Why each angle failed Astral

### Phase A flags
- Astral spectra are already cleaner than TMT's at the resolution where deisotoping is meaningful. Most apparent isotope clusters at TMT's CID resolution are partially merged at the instrument on Astral. Less to deisotope → less benefit.
- Cap=200 too aggressive for Astral. Astral peptides extend to high m/z; mid-intensity diagnostic peaks above the top-200 cutoff drop, hence the −3.6 % target count.
- Net: deisotoping adds per-spectrum overhead that exceeds the cheap-score savings on Astral. Cap throws away signal.

### GF candidate cap
- Astral match queues are typically ≤5–10 entries (10 ppm precursor + small isotope window + 32 MB FASTA). The cap=10 didn't bite (`size > cap` guard skipped the cap path on most spectra).
- cap=5 did bite a small fraction of spectra. The sort+truncate overhead exceeded the GF DP-table savings; Astral wall went up, not down.
- Conclusion: capping is a workload optimization for cases with large per-spectrum candidate sets. Astral's tight precursor window doesn't have that shape.

### Scorer optimization
- JFR showed `NewRankScorer.getIonExistenceScore` → `HashMap.getNode` was ~14 % of Astral CPU samples.
- Fix correctly eliminated those calls (verified via post-fix profile not run, but field cached and used at the call sites). Native counts bit-identical.
- Wall did **not** improve. Likely the JIT was already inlining/escape-analyzing the HashMap lookup; the "fix" replaced a JIT-optimized call with a field load, equivalent cost in real terms.
- This is the post-mortem-fragment-index lesson #3 hitting again: *"three session-worth of micro-opts each measured NEGATIVE impact despite looking sensible on paper. The JVM's JIT optimizer is sophisticated; we reach for machine-level tuning too early."*
- A real fix would need to eliminate the HashMap *invocation overhead* not just the lookup — e.g., split the per-Partition tables into a `PartitionScoringContext` value object created once and held by reference. But the JIT may already handle that for us; need to instrument before betting.

## Lessons learned

1. **TMT is not a reliable Astral proxy on per-spectrum optimizations.** TMT's 20 ppm precursor window + lower MS2 resolution + Lumos peak density gave us a 1.41× win on Phase A that did not transfer. This is the post-mortem-fragment-index lesson #4 again: *"small-FASTA benchmark is NOT a proxy for large-FASTA"* — restated as "high-precursor-tolerance ≠ low-precursor-tolerance for per-spectrum work." The TMT-as-inner-loop strategy from the 5× roadmap §3.1 is unsafe for any optimization whose leverage depends on candidate-density dynamics.
2. **Astral wall on dev-tip is at or near the JIT-optimized floor for the current SA-walk + GF architecture.** Six measurement variants, none beat baseline by more than noise. Phase B (calibrated tolerance), Phase C (branch-and-bound), Phase E (parallelism) — all from the 5× roadmap — remain candidates, but each requires architectural change, not micro-optimization.
3. **The post-mortem-fragment-index's lessons #3 and #4 are the dominant risks** for any future Astral attempt. JIT already compiles aggressively; profile-sample counts overstate optimization headroom; small-FASTA-or-different-instrument benchmarks lie.
4. **Profile before betting on a hot-spot fix.** The JFR profile correctly identified the dominant hot spot, but eliminating it didn't translate to wall improvement. Future profile-driven attempts should run a *post-fix profile* before trusting the JFR delta.
5. **Native target/decoy drift is a leading indicator.** Phase A's −0.25 % targets / −4.6 % decoys on TMT is the same shape, in miniature, as the recall regression that would have killed the experiment in production. If counts drift more than 0.5 % vs OFF on a measurement run, the optimization is not bit-identical-correctness and needs deeper recall validation before shipping.

## Phase E parallelism investigation (added 2026-04-28, also reverted)

After the Phase A retrospective above was committed, a follow-up Phase E
attempt was made: thread-scaling sweep + ForkJoin-pool default selection.
Findings recorded here for completeness; the code change was reverted
because measurement variance was too high to confidently ship.

**Thread-scaling sweep (default `ThreadPoolExecutor`, no flag overrides):**

| Threads | Wall (s) | Note |
|---:|---:|---|
| 4 | 690.1 | morning baseline |
| 6 | 675.0 | within noise of 4t |
| 8 | 884.0 | **+28 % vs 4t — anti-scaling** |

**ForkJoin opt-in (`-Dmsgfplus.useForkJoin=true`):**

| Threads | Wall (s) | Note |
|---:|---:|---|
| 4 | 872.3 | +26 % vs default 4t — ForkJoin loses badly here |
| 6 | (killed) | run was at >1500 s wall when stopped; either hung or extreme regression |
| 8 | 520.9 | 1.32× vs default 4t baseline — only variant that cleared the 1.15× gate |

**Smart-default attempt:** modified `MSGFPlus.runMSGFPlus` to auto-pick ForkJoin
when `numThreads >= 8` (preserving 4t default-executor behaviour, activating
ForkJoin only at the measured-win threshold). Code compiled, scoped tests
passed (9/9 incl. concurrent + telemetry + precursor-cal scaffolding).

**Confirmation runs (same JAR, smart-default change in flight):**

| Run | Threads | Wall (s) | Expected | Δ |
|---|---:|---:|---:|---:|
| auto-FJ | 8 | 861.5 | ~520 | **+65 %** vs morning explicit-FJ |
| auto-default | 4 | 904.3 | ~690 | **+31 %** vs morning measurement |

Same JAR semantically (verified via `unzip -p ... | strings` finding the new
`useForkJoinProp` symbol in the bytecode), same `-thread N` args, same
spectrum/FASTA/mods. **Both metrics regressed ~30 % vs morning.** The
machine state degraded across the day's benchmarking — likely thermal,
accumulated process state, or background macOS work.

**Conclusion:** the morning's ForkJoin-8t = 521 s measurement may have been
real or may have been an outlier. With 30+ % run-to-run variance on the
same JAR across hours, point measurements cannot distinguish a genuine
1.3× ForkJoin win from a 30 % machine-state fluctuation. Reverted the
smart-default change; the underlying `-Dmsgfplus.useForkJoin=true` opt-in
remains in dev unchanged.

**GC-pressure follow-up (2026-04-28, end of iteration):**

After the smart-default revert, JFR analysis of the morning 4t profile
showed **zero `JavaMonitorEnter` contention events and 100 %
RUNNABLE samples** — confirming the 8t regression is not synchronized-lock
contention. But 588 K `GCPhaseParallel` events suggested GC could be the
cause. Tested by re-running 8t and 4t with `-Xmx16g` (double the heap):

| Run | Wall (s) | RSS (MB) | GC count |
|---|---:|---:|---:|
| 8t + Xmx16g | 776.1 | 5 067 | 182 |
| 4t + Xmx16g | 870.0 | 6 120 | 184 |
| (compare) 8t + Xmx8g afternoon | 861.5 | 6 083 | (n/a) |
| (compare) 4t + Xmx8g afternoon | 904.3 | 5 953 | (n/a) |

GC-pressure hypothesis is *partially* confirmed: bigger heap helped 8t by
~12 % wall (and dropped peak RSS by ~17 % because G1GC ran fewer
collections) but only ~4 % at 4t. So GC contributes to the 8t regression
but is not the entire story. Even with -Xmx16g, 8t is slower than the
morning's 4t-Xmx8g baseline (776 vs 690 s). **No actionable recommendation
to ship: heap-tuning helps 8t, but 8t still isn't competitive with 4t at
default heap.**

The afternoon's 4t-Xmx8g (904 s) vs morning's 4t-Xmx8g (690 s) is a
+31 % gap on the same JAR / same args / same machine — confirming the
day's accumulated machine-state degradation dwarfs any code-level signal.
Six hours of benchmarking has hit the noise floor.

**Replication batch (2026-04-29 morning, quieter machine, the iteration's final shot):**

To bound how much of the apparent ForkJoin win was machine-state vs real,
ran three Astral variants in tight back-to-back sequence on a less-loaded
machine:

| Run | Wall (s) | RSS (MB) |
|---|---:|---:|
| 4t default | 963.1 | 5 519 |
| 8t default | 918.3 | 5 740 |
| **8t ForkJoin** | **978.8** | 5 204 |

All three within 6.5 % of each other (within noise). 8t-default is now
*faster* than 4t-default by 4.7 % — directly opposite to yesterday's
"anti-scaling" finding. **The yesterday-morning 4t=690 s baseline was an
outlier**, not the truth — the 921 s machine reality was masked by a
single fortunate quiet-machine measurement that morning. **The 521 s
ForkJoin-8t was likewise an outlier**, not a real 1.32× win — three
independent re-measurements (afternoon 861 s, today's 978 s) put it
solidly above 850 s.

**Corrected conclusion:** there is no Phase E win to ship. The "default
executor anti-scales past 6 threads" claim earlier in this retrospective
was *wrong*; it was a one-day correlation between morning-quiet-machine +
4t and afternoon-noisy-machine + 8t, not a real algorithmic relationship.
The ForkJoin path doesn't outperform default executor on Astral when
measured in clean within-batch conditions. The single 521 s ForkJoin
data point was unreplicable noise.

**What future agents need to do this safely:**

1. **Stable benchmark environment.** A reserved CI runner, an idle box with
   thermal headroom, or a cloud VM with fixed CPU allocation. Not a
   developer workstation that's been running benchmarks for hours.
2. **Multi-run statistics, not point measurements.** Each variant run 3-5
   times; report median + IQR. A single 521 s measurement that doesn't
   replicate is a noise artefact, not a discovery.
3. **Same-day sweep with fixed ordering.** Run all variants back-to-back
   in the same machine state so cross-variant comparisons are valid.
4. **Anti-scaling at 8t default-executor IS reproducible** (884 s and 861 s
   in two measurements at different machine states; the relative slowdown
   vs 4t survives the variance). That finding is real and worth digging
   into — what's the contention point in `ThreadPoolExecutorWithExceptions`
   that causes 8t to lose to 4t? `jfr print --events jdk.JavaMonitorWait`
   on the 8t default-executor profile would identify the lock.
5. **The post-mortem-fragment-index lesson #3 strikes again:** *"the JVM's
   JIT optimizer is sophisticated; we reach for machine-level tuning too
   early."* Wall-time deltas at the 30 % level are below the noise floor
   for a single-machine benchmark of this size. Don't claim a win from
   one measurement.

## What's still untried (for future agents)

The 5× roadmap (`astral-speed-5x-roadmap.md`) specified five phases. Only Phase A was attempted. Remaining:

- **Phase B — calibrated precursor-window tightening.** Use Achievement B's calibration σ to shrink the effective precursor window post-calibration. Reduces candidate fan-out at the `pepMassSpecKeyMap.subMap(...)` site, which IS measurable in the current JFR profile (TreeMap operations ~4 % of CPU). Recall-risky; needs an integration test that asserts no FDR-1 % PSM survives outside the tightened window.
- **Phase C — branch-and-bound during peptide extension.** The roadmap's centerpiece (1.5–2.5× projected). My review of the roadmap (in the git history before the reset, see commit `eee9fa6`'s plan) flagged three concrete sub-problems: dynamic threshold rises late in the SA walk, admissible-yet-selective upper bound is hard to define for a rank-based scorer, per-spectrum bookkeeping cost may exceed savings. Research-grade; should be planned as a multi-iteration investigation with a kill-by-exactness-audit clause.
- **Phase D — GF threshold tightening via `setUpScoreThreshold`.** The current code already passes `minScore` to GF; tightening this further requires raising minScore by capping candidates (Angle 2 in this retrospective), which we showed doesn't bite on Astral. Phase D is unlikely to be useful as a standalone lever on Astral.
- **Phase E — parallelism ceiling investigation.** Attempted 2026-04-28, multi-run replicated 2026-04-29 (see "Phase E parallelism investigation" + "Replication batch" above). **Initial "anti-scaling" finding was disproved by the replication batch** — when measured back-to-back in the same machine state, 8t-default is actually *faster* than 4t-default. The ForkJoin path also did not show any advantage in within-batch comparison. Both initial findings (anti-scaling + ForkJoin win) were noise artefacts. Future agents wanting a parallelism win must build a stable benchmark environment first; the conclusion changes between runs done at different times of day on this machine.
- **Workload retargeting** — the original branch-name framing ("feat/big-fasta-peptide-candidate") was about metaproteomics / proteogenomics big-FASTA workloads, not Astral. Astral was a redirect during brainstorming. The big-FASTA framing has different bottlenecks (peptide redundancy across organisms, candidate dedup) that may be more amenable to per-spectrum optimization. Worth profiling on a metaproteomics dataset before assuming any per-spectrum lever is dead.
- **HashMap-elimination in NewRankScorer (deeper version).** Angle 3 in this retrospective tried the shallow version (cache the array). A deeper version would refactor all 10 per-Partition `HashMap`s in `NewRankScorer` into a `PartitionScoringContext` record, looked up *once per spectrum* and held by reference for the duration of scoring. The shallow fix didn't move wall, but the deeper refactor *might* — JIT optimization of the lookup vs an entire object indirection chain is the open question. Should not be attempted without a post-fix profile to confirm the win.

## Files and artifacts

- This retrospective: `.claude/plans/astral-phase-a-retrospective.md`
- Original Phase A implementation plan (now reverted; recoverable): `git show 6510f08:.claude/plans/astral-speed-phase-a-plan.md`
- Active 5× roadmap (still authoritative for future iterations): `.claude/plans/astral-speed-5x-roadmap.md`
- Earlier shipped retrospective: `.claude/plans/SHIPPED.md`
- JFR Astral profile: `~/work/msgfplus-workspace/benchmark/results/phaseA/astral_off.jfr`
- All measurement summary TSV: `~/work/msgfplus-workspace/benchmark/results/phaseA/summary.tsv`
- Reverted Phase A code recoverable from: `git show 5cdd21e` and walking back through `b78e275..5cdd21e` (11 commits: SearchTelemetry, telemetry CLI/refactor/wiring, Deisotoper, Spectrum.deisotope/capByIntensity, deisotope CLI flag, ScoredSpectraMap wiring).
