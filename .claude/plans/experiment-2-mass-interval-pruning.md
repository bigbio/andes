# Experiment 2 — Exact Prefix Mass-Interval Pruning

**Status:** Design + Checkpoint 1 + Checkpoint 2 shipped 2026-04-29 (commit `4241fbb`); off by default (opt-in via system property). **Wall gate FAILED** — bookkeeping cost exceeds savings. Checkpoint 3 (overhead optimization) is the open follow-on.

> **Result summary (Astral, remote pride-linux-vm.ebi.ac.uk):** native counts bit-identical to baseline in all variants (exact-by-construction validated ✓); 12.22 % prune rate at Checkpoint 1, 1.84 % with actual break; **but Phase B + Experiment 2 pruning regresses wall +11 %** vs Phase B alone (549 s vs 494 s). The bound test runs ~1.4 B times at ~40 ns each = ~55 s of pure overhead, dwarfing the cheap-score savings (the pruned pairings had ~0 matches anyway). Phase B remains the iteration's shippable Astral wall lever; Experiment 2 is parked as opt-in scaffolding for Checkpoint 3 follow-on work.
**Date:** 2026-04-29
**Context:** Phase B (commits `aac389c` and earlier) shipped −10.4 % Astral wall via calibrated precursor-window tightening. Plan §5 names this as the natural next attack — exact-by-construction pruning that attacks SA-walk fan-out *before* Phase B's pairing fan-out reduction kicks in. The two compose: Phase B reduces matched_speckeys per pairing call; Experiment 2 reduces the number of pairing calls.

## 1. Goal

For a partial peptide prefix of length `L` (currently being extended by `DBScanner.dbSearch`), compute the interval `[minMass, maxMass]` of all final-peptide masses reachable by extending this prefix. If the interval cannot intersect any spectrum's precursor-mass window, the entire branch is dead — stop extending.

Exact by construction: the bound is the actual reachable interval, not a heuristic upper score bound. No recall risk. Skips peptide variants that would produce zero matches.

## 2. Where the code lives

The SA walk happens inside `DBScanner.dbSearch(...)` ([src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:189](../../src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java)). The relevant inner loop is around lines 370–490:

```java
// Loop iterates over residues in the SA walk
for (...) {
    // 1. Extend prefix by one residue:
    candidatePepGrid.addResidue(peptideLengthIndex, residue);   // line 389

    if (peptideLengthIndex < minPeptideLength) continue;        // line 412

    // 2. For each variant in the grid, look up matching SpecKeys:
    for (int j = 0; j < candidatePepGrid.size(); j++) {
        float theoPeptideMass = candidatePepGrid.getPeptideMass(j); // line 466
        // ... compute tolerance window, subMap query, cheap-score loop
        // (PhaseBTelemetry.recordPairing(matchedSpecKeyList.size()) hook here)
    }
}
```

The pruning hook goes **between extending the prefix and entering the variant loop** — i.e., right after `addResidue` succeeds and before line 412's `continue` / line 466's variant loop.

## 3. Bound construction

For a prefix of length `L` with current variant masses `{m_1, ..., m_k}` (one per modification variant in the grid):

```
prefixMinMass = min(m_i)
prefixMaxMass = max(m_i)
```

Remaining residues can be at most `R_max = maxPeptideLength - L` and at least `R_min = max(0, minPeptideLength - L)`. Each remaining residue adds an amino-acid mass; with modifications, the maximum addition per residue is `maxAaMass + maxModMass` and the minimum is `minAaMass`.

```
reachableMin = prefixMinMass + R_min * minAaMass
reachableMax = prefixMaxMass + R_max * (maxAaMass + maxResidueModMass) + maxFixedTermModMass
```

Two simplifications keep the bound construction cheap:

1. Cache `minAaMass`, `maxAaMass`, `maxResidueModMass`, `maxFixedTermModMass` as fields of `DBScanner` at construction time (once per task).
2. If the grid maintains `getMinPeptideMass()` / `getMaxPeptideMass()` accessors that scan the variants array, that's `O(numVariants)` per call (~tens of variants). Pre-cached if hot.

## 4. Intersection test with spectrum windows

`specScanner.getPepMassSpecKeyMap()` is a `TreeMap<Double, SpecKey>` keyed on peptide mass. Each spectrum has tolerance windows `[leftThr, rightThr]` around its precursor peptide mass.

For the pruning test we need: *"does any spectrum's window touch the reachable interval `[reachableMin, reachableMax]`?"*

Two equivalent formulations:
- **Per-spectrum view**: for each SpecKey with peptide mass `p`, its window is `[p - tolDaLeft(p), p + tolDaRight(p)]`. Branch is alive iff `[reachableMin, reachableMax] ∩ [p - tolDaLeft(p), p + tolDaRight(p)] ≠ ∅` for some SpecKey.
- **Aggregate view**: precompute the *expanded* TreeMap key = `p` (unchanged) but query with widened bounds: `pepMassSpecKeyMap.subMap(reachableMin - maxToleranceDa, reachableMax + maxToleranceDa)`. If empty, branch is dead.

The aggregate view is `O(log N)` in TreeMap size; the per-spectrum view would be `O(N)`. Use aggregate.

`maxToleranceDa` can be precomputed at task start using the post-Phase-B effective tolerance and the largest peptide mass we'd query at: `effectiveLeftPrecursorMassTolerance.getToleranceAsDa(maxPeptideMass)` plus the right-tolerance equivalent.

## 5. Where the bound is most effective

The pruning saves work proportional to how often it fires. Heuristic estimate:

- Long-peptide branches: when `prefixMass` is already large and the remaining-residue reach can't bring it down enough to touch any spectrum. Bound is loose for short prefixes (lots of headroom) but tight for prefixes near `maxPeptideLength` where there's little room to add mass.
- Off-mass branches: when the prefix's accumulated mass is in a "gap" of the spectrum mass distribution. With Astral's ~50 K spectra spanning ~4 kDa, the spectrum mass distribution is dense; gaps narrow.

**Decision:** instrument the prune rate via a counter (similar to `PhaseBTelemetry`) before optimizing. If pruning fires < 1 % of pairing-call sites, the bookkeeping cost wins. If it fires > 5 %, we have a real lever.

## 6. Implementation checkpoints

Bounded scope, in order:

### Checkpoint 1 — instrument first

Add `Experiment2Telemetry` (mirrors `PhaseBTelemetry`):
- `prefixesEvaluated` — how many prefix-extension steps reach the pruning hook
- `prefixesPruned` — how many were eliminated by the mass-interval test
- `pruneRatio` printed at end of search

Implement WITHOUT actually pruning (just compute the bound, count would-be prunes). Run once on Astral with Phase B AUTO. Decide whether to proceed based on the rate.

### Checkpoint 2 — minimal pruning

If Checkpoint 1 shows ≥ 5 % prune rate, add the actual `break` statement in the SA walk when the bound test fails. Re-measure on Astral OFF + AUTO; verify no recall regression (target/decoy counts bit-identical to Phase B baseline).

### Checkpoint 3 — sharpening

Tighten the bound by:
- Per-residue mod-mass cap (some residues admit specific mods; the global `maxResidueModMass` overestimates)
- Cleavage-site constraints (if the next residue isn't cleavable for the enzyme, `R_min` floor rises)

Only pursue if Checkpoint 2 shows wall improvement but the prune ratio is below the theoretical maximum.

## 7. Acceptance / kill gates (from plan §5.6 / §5.7)

**Acceptance:**
- Astral prune rate ≥ 5 % (Checkpoint 1 telemetry)
- Astral wall improves ≥ 5 % vs Phase B baseline (Checkpoint 2 wall)
- Native target counts bit-identical (exact-by-construction)

**Kill:**
- Prune rate < 1 % (bookkeeping > savings)
- Or prune rate adequate but wall doesn't move (downstream still bottleneck)
- Or correctness drift (target/decoy counts differ from Phase B baseline)

## 8. Files to touch

- `src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java` — pruning hook in dbSearch loop; cached aa-mass bounds
- `src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGrid.java` — `getMinPeptideMass()` / `getMaxPeptideMass()` if not already exposed
- `src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGridConsideringMetCleavage.java` — same accessor in the Met-cleavage variant
- `src/main/java/edu/ucsd/msjava/msdbsearch/Experiment2Telemetry.java` (new) — `LongAdder` counters
- Tests: scoped unit + integration verifying OFF-mode bit-identical

## 9. Why this is safe to ship as-designed

The bound is **exact-by-construction**: a peptide whose final mass falls outside `[reachableMin, reachableMax]` cannot be the result of extending this prefix. This is mathematically certain, not a probabilistic argument. So the only failure mode is "bound is correct but bookkeeping cost > savings," which the Checkpoint 1 telemetry catches before any production code path changes.

This is the property that makes Experiment 2 distinct from Phase A's deisotoping (which trades correctness for speed) and Phase B's tightening (which trades a small recall risk via 3-σ envelope for speed). Experiment 2 is purely a work-elimination optimization.

## 10. Reference

- Plan: [`astral-next-experiments.md`](astral-next-experiments.md) §5
- Phase B (the lever this composes with): [`SHIPPED.md`](SHIPPED.md)
- Long-horizon roadmap: [`astral-speed-5x-roadmap.md`](astral-speed-5x-roadmap.md)
