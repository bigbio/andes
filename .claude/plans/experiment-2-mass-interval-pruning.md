# Experiment 2 — Exact Prefix Mass-Interval Pruning

**Status:** Design + Checkpoints 1, 2, 3, 4 shipped 2026-04-29; off by default (opt-in via system property). **Wall gate FAILED across all four checkpoints** — bookkeeping cost approaches but never beats savings. Phase B (commit `aac389c`) remains the durable Astral wall lever; Experiment 2 stays as opt-in scaffolding for future work.

> **Result summary (Astral, remote pride-linux-vm.ebi.ac.uk):** native counts bit-identical to baseline in all variants (exact-by-construction validated ✓); 12.22 % prune rate at Checkpoint 1, 1.84 % with actual break.
>
> **Checkpoint 2** (TreeMap.subMap bound test): Phase B + E2 pruning = 549 s vs Phase B alone 494 s (**+11 % wall regression**). Bound test ~150 ns × 1.4 B = ~210 s of overhead.
>
> **Checkpoint 3** (commit `0c697dd`, binary-search via `ScoredSpectraMap.hasSpecMassInRange`): bound test ~30 ns × 1.4 B = ~42 s overhead. Phase B + E2 pruning = 511 s vs Phase B alone 494 s (**+3.4 % wall regression** — still narrowly negative but ~75 % of the gap closed). OFF + E2 pruning = 559 s vs OFF baseline 551 s = +1.5 % (break-even within noise).
>
> **Checkpoint 4** (commit `8478651`, gate bound test on `peptideLengthIndex >= minPeptideLength`): short prefixes (length 1 to minPeptideLength-1) have reachable intervals many kDa wide and almost never prune; bound test there is dead weight. Code change skips ~3.7 % of evaluations (1.61 B → 1.55 B) without recall risk (the moved test is still a sound necessary condition). Tests pass (37/37 scoped suite).
>
> **Checkpoint 4 confirmation (5-trial interleaved bench, 2026-04-30):** with run-to-run variance properly accounted for, the effect is real but small.
>
> | config              | trials (s)              | n | mean (s) | σ (s) |
> |---------------------|-------------------------|---|---------:|------:|
> | `phaseB_only`       | 522, 518, 517, 529, 513 | 5 |    519.8 |  6.06 |
> | `phaseB_plus_e2`    | 504, 507, 514, 503, 512 | 5 |    508.0 |  4.85 |
>
> **Δ = 11.8 s = −2.27 % vs Phase B alone**, 5/5 trials phaseB+E2 < phaseB_only, Welch's t ≈ 3.4 (p ≈ 0.01). Native target/decoy counts **bit-identical 89580 / 45292 across all 10 runs** — exact-by-construction validated at scale.
>
> **Verdict:** four checkpoints of optimization. The pruning is a real, statistically significant ~2.3 % improvement, but doesn't clear the plan's ≥5 % gate for default-on. Phase B remains the durable Astral wall lever (−10.4 % vs OFF). Experiment 2 stays as opt-in via `-Dmsgfplus.experiment2Pruning=true` — costs nothing in OFF mode, and is available for users who want to stack a small additional gain on top of Phase B. Future paths if revisited: incremental prefix-mass cache (avoid the per-extension grid scan), or coarse-grained per-LCP-block bound (amortize across many SA traversals).
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
