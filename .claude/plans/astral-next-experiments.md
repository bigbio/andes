# Astral Next Experiments — Post-Retrospective Action Plan

**Status:** Historical staging plan — Phase B shipped; Experiment 2 retrospective completed
**Date:** 2026-04-30
**Purpose:** define the next implementation after shipping Phase B and retiring the sub-threshold Experiment 2 runtime scaffolding

> **Update (2026-04-30): Phase B shipped; Experiment 2 did not graduate.** After the calibrator iso=0 fix (`05ec066`), the configurable formula constants (`7c027f8`), and the spec_eValue stratification (`aac389c`), the AUTO-mode stratified calibrator delivers **−10.4 % Astral wall** and is the durable improvement from this branch. Experiment 2 later produced a real but smaller **−2.27 %** add-on in a 5-trial bench, but stayed below the 5 % default-on gate and was removed from the cleaned shipping runtime path. The next implementation should therefore aim at a larger algorithmic reduction, not more branch-local cleanup.

## 1. What changed

Two earlier ideas have now been materially de-risked in the wrong direction:

- **Phase A is disproven on Astral.**
  Deisotoping + peak cap, GF candidate cap, and the shallow scorer hot-path tweak all failed the Astral wall gate and were reverted.

- **Phase E is not a current win.**
  The later replication batch showed the initial executor/ForkJoin signal was noise. We should not spend another immediate iteration on pool-default tuning on this workstation.

The practical implication is:

- **do not start with spectrum cleanup**
- **do not start with executor tuning**
- **do not start with another shallow hotspot fix**

The next experiments should attack the real multiplicative fan-out in [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:189) with exact or near-exact levers.

## 2. Updated priority order

### Priority 1 — Phase B: calibrated precursor-window tightening

This is the best next coding experiment.

Why it survives the latest comments:

- still untried
- already has clean seams in `MassCalibrator` and `ScoredSpectraMap`
- reduces fan-out at the real pairing site (`pepMassSpecKeyMap.subMap(...)`)
- exact OFF-mode path can be preserved

### Priority 2 — Exact prefix mass-interval pruning

This is the safest version of the earlier Phase C thinking.

Important correction from the retrospective:

- **do not** open with full score-based branch-and-bound
- **do** start with exact mass reachability pruning on partial peptide prefixes

Why:

- avoids the hardest “admissible score upper bound” problem
- exact by construction
- still attacks the multiplicative fan-out before cheap scoring

### Priority 3 — Persistent mass-indexed peptide DB design spike

This is the strongest “crazy but plausible” architectural alternative still on the table.

But it should start as a design/prototype exercise, not a full implementation branch.

## 3. Experiments we should not do next

- re-attempt Phase A spectrum cleanup on Astral
- another GF candidate cap variant
- another shallow scorer-map optimization
- executor/ForkJoin default changes on this machine
- full score-bound branch-and-bound as the opening pruning experiment

## 4. Experiment 1 — Phase B implementation

> **Status (2026-04-29): core implementation already shipped in dev; telemetry added in this iteration.**
>
> Inspecting `MSGFPlus.runMSGFPlus` lines 396–423 shows Phase B's tightening logic is already in place: when `MassCalibrator.CalibrationStats.hasReliableStats()` is true and the precursor tolerance is ppm-based, `MassCalibrator.tightenedTolerancePpm(...)` is computed for left and right tolerances using the canonical formula `min(userPpm, max(floorPpm, k·robustSigma + marginPpm))` with the documented constants (`floor=2 ppm`, `margin=0.5 ppm`, `k=3`). The `effectiveLeftPrecursorMassTolerance` / `effectiveRightPrecursorMassTolerance` finals are then captured by the per-task `ScoredSpectraMap` Supplier lambda (line 510-511) so the main pass uses the tightened window. OFF mode is bit-identical (early-return at line 362 when `precursorCalMode == OFF`).
>
> Missing piece — **telemetry to verify Phase B's effect on pairing fan-out** — added in commit on this branch via `PhaseBTelemetry`. Enable with `-Dmsgfplus.phaseBTelemetry=true`; emits `pairing_calls`, `matched_speckeys`, and `mean_per_call` summary at end of search. Hooked at `DBScanner.dbSearch:489` (the `pepMassSpecKeyMap.subMap(...)` site). 5 unit tests + the existing `TestPrecursorCalScaffolding` integration confirm OFF-mode bit-identical.
>
> Original Experiment 1 spec preserved below for context; the success/kill gates still apply (the next agent runs the bench with telemetry on, then verifies the gate).

## 4.1 Goal

Shrink the effective precursor tolerance after calibration so the engine does less work at:

1. peptide↔spectrum pairing
2. precursor-mass-index GF expansion

## 4.2 Files to touch

- [MassCalibrator.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/MassCalibrator.java:37)
- [ScoredSpectraMap.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/ScoredSpectraMap.java:14)
- [SearchParams.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/SearchParams.java:18)
- [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:471)

## 4.3 Implementation sketch

1. Extend calibration output from:
   - `shiftPpm`

   to:
   - `shiftPpm`
   - robust spread estimate (`mad`, `robustSigma`)

2. Compute tightened ppm window only when:
   - user tolerance is ppm-based
   - calibration produced enough confident PSMs
   - tightened window is smaller than the user window

3. Suggested initial formula:
   - `tightenedPpm = min(userPpm, max(floorPpm, k * robustSigma + marginPpm))`

4. Preserve the exact no-op path for:
   - `-precursorCal off`
   - insufficient calibration evidence

## 4.4 Telemetry

Add behind a debug flag:

- original precursor window ppm
- tightened precursor window ppm
- matched `SpecKey` count per candidate peptide
- GF precursor-mass-index span per spectrum
- count of spectra where no tightening occurred

## 4.5 Success gate

- Astral median window width shrinks materially
- matched `SpecKey` count drops materially
- GF mass-index span drops materially
- Astral wall improves by at least ~10 %
- no meaningful native target/decoy drift
- no regression below the Astral 1 % FDR gate

## 4.6 Kill gate

- window shrinks but pairing count barely changes
- pairing count drops but wall barely changes
- or recall drifts beyond gate

## 5. Experiment 2 — Exact prefix mass-interval pruning

## 5.1 Goal

Kill peptide-extension branches early when the current prefix cannot possibly end in a mass that overlaps any surviving spectrum window.

## 5.2 Why this is the right Phase C opening

The retrospective correctly flagged that full score-bound pruning has three hard problems:

- dynamic thresholds rise late
- admissible-yet-selective score bounds are hard for a rank-based scorer
- per-spectrum bookkeeping may exceed savings

Exact mass-interval pruning avoids those first two problems entirely.

## 5.3 Files to touch

- [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:189)
- [CandidatePeptideGrid.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGrid.java:152)
- [CandidatePeptideGridConsideringMetCleavage.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGridConsideringMetCleavage.java:6)

## 5.4 Implementation sketch

For each partial peptide prefix:

1. compute the minimum reachable final peptide mass
2. compute the maximum reachable final peptide mass
3. account for:
   - remaining peptide length range
   - modification budget
   - enzyme / terminal constraints
   - Met-cleavage branch if active

If the reachable interval cannot intersect any spectrum mass window, stop extending that branch.

This should be implemented before the inner cheap-score fan-out loop, not after.

## 5.5 Telemetry

- prefixes considered
- prefixes killed by exact mass-interval test
- cheap-score calls avoided
- branch kill ratio by peptide length
- runtime overhead of interval bookkeeping

## 5.6 Success gate

- substantial prefix-prune ratio on Astral
- substantial cheap-score call reduction
- Astral wall improves by at least ~15 %
- zero correctness drift by construction

## 5.7 Kill gate

- pruning ratio too small to matter
- interval bookkeeping overhead cancels the gain

## 6. Experiment 3 — Persistent mass-indexed peptide DB design spike

## 6.1 Goal

Test whether there is a viable middle ground between:

- current live SA walk
- abandoned fragment index

The target concept is:

- store a persistent peptide catalog keyed by precursor mass slabs
- query only relevant slabs at search time
- avoid rebuilding digestion state every run
- avoid storing fragment-index-style heavy Tier-1 structures

## 6.2 Scope of the spike

Do **not** build the full system in this experiment.

Instead produce:

1. file-format sketch
2. build-time complexity estimate
3. query-time complexity estimate
4. memory model
5. variable-mod handling strategy

## 6.3 Constraints

- do not pre-expand all modified variants if that recreates fragment-index memory blow-up
- prefer storing unique peptide backbones plus cleavage/source metadata
- treat variable modifications lazily inside selected precursor-mass slabs

## 6.4 Success gate

- design shows a plausible path to lower repeated runtime work
- memory model looks much safer than fragment index
- mod strategy does not immediately collapse into full runtime expansion

## 6.5 Kill gate

- design complexity explodes immediately
- or lazy-mod generation just recreates current runtime cost

## 7. Recommended implementation order

1. **Phase B implementation**
2. **Exact prefix mass-interval pruning prototype**
3. **Persistent peptide-DB design spike**

This order reflects the latest retrospective comments:

- start with the cleanest still-untried exact lever
- then try the safest pruning form of Phase C
- only then invest in a larger architectural alternative

## 8. Benchmark rules for these experiments

The latest comments changed the benchmark protocol:

1. **Astral is the primary truth dataset.**
   Do not accept TMT as a transfer proxy for these optimizations.

2. **Use TMT only as auxiliary signal** if the optimization is clearly not per-spectrum-shape-sensitive.

3. **Measure variants back-to-back in the same machine state** when possible.

4. **Do not trust single point measurements** for threading or wall claims on this workstation.

5. **Native target/decoy drift is an early warning signal.**

## 9. What I recommend we do now

If we are spending one serious coding week, I would use it on:

- **Phase B implementation plus telemetry**

If that shows the expected drop in pairing fan-out, then the next week goes to:

- **exact prefix mass-interval pruning**

If Phase B does **not** move the pairing counts enough, then I would pause before any more Astral coding and do the peptide-DB design spike instead of forcing Phase C.

## 10. Reference

- Phase A retrospective: [astral-phase-a-retrospective.md](/Users/yperez/work/msgfplus-workspace/astral-speed/.claude/plans/astral-phase-a-retrospective.md:1)
- Long-horizon roadmap: [astral-speed-5x-roadmap.md](/Users/yperez/work/msgfplus-workspace/astral-speed/.claude/plans/astral-speed-5x-roadmap.md:1)
- Short retrospective: [SHIPPED.md](/Users/yperez/work/msgfplus-workspace/astral-speed/.claude/plans/SHIPPED.md:1)
