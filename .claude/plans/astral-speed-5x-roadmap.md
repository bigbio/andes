# Astral 5X Roadmap — Search-Space Reduction Fast Path

**Status:** Design / exploratory roadmap
**Date:** 2026-04-28
**Scope:** first credible path toward a **5× Astral wall-time reduction** without giving back MS-GF+'s sensitivity lead

## 0. Shipping model

This iteration ships as **milestone commits** on `feat/astral-speed-improvements`, with **one closing PR** opened at the end of the iteration. Phases do not become individual PRs.

Each phase milestone uses a commit message of the form:

```
feat(astral-speed): MILESTONE Phase <id> — <one-line achievement>

<2–4 lines of measurement detail>
- TMT inner-loop wall delta
- Astral phase-gate result (if run)
- Recall delta on Astral 1 % FDR
- Any new memory or RSS constraint observed
```

Strategy: try the highest-EV phase first; fall back to smaller wins inside the same branch if a phase fails its kill gate.

- **Attempt order:** Phase A → (success: Phase B or C; failure: Iteration 0.5 fallback below) → ...
- **Iteration 0.5 fallback** (used only when a "big-win" phase fails its kill gate): graph-skeleton memoization in `PrimitiveAminoAcidGraph` (~10–15% Astral, recall-neutral) + Tier-1.5 GF candidate cap in `DBScanner.computeSpecEValue` (15–30%, recall-gated). Both are single-site changes and ship as their own milestone commits before this branch's closing PR.
- **Closing PR** is opened only after measured Astral wall improvement on the branch passes the whole-roadmap gate (§8) or after the fallback path delivers a defensible improvement.

Throughout the iteration the branch is visible to reviewers via its commit log; no per-phase PR review.

## 1. Executive view

A real 5× Astral gain means moving from roughly **620 s** to **124 s** on the clean 4-thread baseline.

That is **not** a "next hotspot fix" target.

The current architecture spends most of its time doing legitimate work:

1. walking the suffix-array-derived peptide space
2. matching many peptide masses to many spectra
3. cheap-scoring the matched peptide/spectrum pairs
4. computing GF over the retained precursor-mass window

Even perfect implementation-level tuning will not get us to 124 s. The only credible path is to do **much less work**.

This roadmap proposes an **Astral fast path** that keeps the current SA-walk engine, but adds three major forms of search-space reduction:

1. **cleaner spectra** before scoring
2. **tighter precursor windows** before peptide↔spectrum pairing
3. **branch-and-bound pruning inside the peptide-extension walk** before cheap scoring / GF

The key decision is architectural:

- **Do not** revive the standalone fragment index
- **Do** insert pruning logic *inside* the current `DBScanner` + `CandidatePeptideGrid` path

## 2. Why 5× is hard in the current shape

The benchmark and profiling history give us two hard constraints:

1. **Parallelism alone is not enough.**
   Astral's clean baseline is about 620 s wall. Earlier measurements showed about 2366 CPU-seconds of real work on 4 threads. Even if we reached perfect 8-core scaling with no other improvements, wall would still be roughly 296 s.

2. **Micro-optimizations are no longer enough.**
   The old big bottleneck (`Hashtable` contention in `NewRankScorer`) has already been addressed on `dev`. The remaining work is spread across candidate generation, cheap scoring, and GF. That means further 5-15% wins are still worth doing, but they will not compound to 5× by themselves.

Conclusion:

- **5× requires both**
  - materially lower CPU work
  - materially better parallel efficiency after that work is reduced

## 3. Working thesis

The best shot at 5× is an **Astral-specific fast path** with this sequence:

1. **MS2 deisotoping + dense-peak retention cap**
2. **calibrated precursor-window tightening**
3. **spectrum-aware branch-and-bound during peptide extension**
4. **score-threshold tightening into GF**
5. **follow-up parallel scaling after the search space is smaller**

The core idea is not to replace the current engine. It is to stop feeding it so many hopeless candidates.

## 4. Where the current code multiplies work

The hottest multiplicative loop in the current search path is in [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:189):

1. extend peptide prefixes along the suffix-array walk
2. materialize peptide variants in `CandidatePeptideGrid`
3. for each candidate peptide variant:
   - compute theoretical peptide mass
   - lookup matched `SpecKey`s via `pepMassSpecKeyMap.subMap(...)`
   - cheap-score each matched spectrum with `scorer.getScore(...)`
   - keep top scoring matches per spectrum

The key inner fan-out is here:

- peptide extension and variant materialization: [CandidatePeptideGrid.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGrid.java:152)
- spectrum matching and cheap scoring: [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:488)
- GF pass over surviving precursor-mass indices: [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:563)

This is the choke point we need to change.

## 5. The proposed fast path

## 5.1 Phase A — Spectrum cleanup before search

### A1. In-engine MS2 deisotoping

Goal:

- collapse isotope clusters so Astral spectra look closer to the effective evidence Sage scores

Why it matters:

- reduces peak density
- reduces noisy evidence in cheap scoring
- should close part of the candidate-generation mismatch seen in the benchmark notes

Expected effect:

- lower cheap-score cost
- stronger score separation for real matches
- modest recall upside on Astral

Likely classes to touch:

- [Spectrum.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msutil/Spectrum.java:18)
- [ScoredSpectraMap.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/ScoredSpectraMap.java:203)
- scorer construction path in `NewScoredSpectrum` / `NewRankScorer`

### A2. Dense-peak retention cap

Goal:

- after deisotoping, keep only the most informative peaks for dense Astral MS2 scans

Suggested initial policy:

- configurable top-N by intensity, with optional windowed cap
- start conservative, e.g. 200-300 peaks

This should be treated as a measured extension of deisotoping, not a separate headline feature.

## 5.2 Phase B — Shrink precursor pairing earlier

### B1. Calibrated precursor-window tightening

Use the existing calibration seam to reduce the peptide↔spectrum pairing fan-out before cheap scoring.

This should be applied in two places:

1. when building `pepMassSpecKeyMap`
2. when choosing the precursor-mass index window for GF

Likely classes to touch:

- [MassCalibrator.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/MassCalibrator.java:37)
- [ScoredSpectraMap.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/ScoredSpectraMap.java:14)
- [SearchParams.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/SearchParams.java:18)
- [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:471)

This is the cleanest already-supported lever for reducing search-space width.

## 5.3 Phase C — Branch-and-bound inside the SA walk

This is the centerpiece of the 5× roadmap.

### C1. The idea

Today we extend peptide prefixes largely on enzyme/modification feasibility, then only later cheap-score the full candidate against all matched spectra.

Instead, we should prune branches *during* extension when they cannot possibly beat the current per-spectrum threshold.

That means we need to attach an **optimistic upper bound** to a partial peptide prefix.

### C2. Bounding model

For a peptide prefix of length `L`, define:

- `partialScore(prefix, specKey)` = cheap score already explained by the prefix
- `upperBoundRemaining(prefix, specKey)` = optimistic best-case contribution from residues not yet appended
- `bound(prefix, specKey)` = `partialScore + upperBoundRemaining + cleavage bonuses`

If:

- `bound(prefix, specKey) < currentWorstTopN(specKey)`

then that prefix cannot produce a retained match for that spectrum, so we stop extending it.

### C3. The practical challenge

We cannot afford to track detailed state for every spectrum on every branch.

So the fast path needs a staged pruning model:

1. **Mass gate**
   Keep only spectra whose tightened precursor window still overlaps the reachable peptide-mass interval from this prefix.

2. **Lightweight evidence gate**
   Maintain a coarse prefix evidence score from the current PRM grid against the spectrum scorer.

3. **Top-N bound gate**
   Prune only when the optimistic bound is safely below the current per-spectrum threshold.

This must be done with compact data structures and aggressive reuse.

### C4. Implementation shape

Introduce a small, explicit pruning helper owned by `DBScanner`, for example:

- `SpectrumPruningState`
- `PrefixBoundCalculator`
- `PrefixCandidateWindow`

Likely responsibilities:

- map prefix mass ranges to candidate `SpecKey` subsets
- maintain current worst top-N threshold per `SpecKey`
- compute an optimistic completion bound
- return `KEEP`, `PRUNE_FOR_SPEC`, or `PRUNE_BRANCH`

Likely classes to touch:

- [DBScanner.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java:189)
- [CandidatePeptideGrid.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGrid.java:11)
- [CandidatePeptideGridConsideringMetCleavage.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGridConsideringMetCleavage.java:6)
- scorer interfaces:
  - [SimpleDBSearchScorer.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msscorer/SimpleDBSearchScorer.java:1)
  - [FastScorer.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msscorer/FastScorer.java:11)
  - [DBScanScorer.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msscorer/DBScanScorer.java:1)

### C5. Important constraint

The first branch-and-bound version must be **conservative**:

- never prune a branch unless the bound is mathematically safe
- if a safe bound proves too weak to save real work, stop and reassess

It is better to discover that a bound is too loose than to ship a fast but recall-damaging heuristic disguised as an exact optimization.

## 5.4 Phase D — Tighten the GF stage after pruning

Once the prefix pruning has already removed much of the cheap-score fan-out, then score-threshold tightening into GF becomes more realistic.

This version of the idea is code-accurate:

- use the retained candidate set to raise the minimum score threshold
- pass that threshold to `PrimitiveGeneratingFunction.setUpScoreThreshold`
- verify that the DP state shrinks materially

This is a **Phase D** optimization, not the centerpiece.

## 5.5 Phase E — Recover parallel scaling after search-space shrinkage

Only after A-D have reduced the amount of work per thread should we chase higher scaling.

Why:

- otherwise we parallelize waste
- contention and overhead are harder to reason about while candidate fan-out is still large

Phase E scope:

- measure 1/2/4/8 thread scaling after branch pruning
- identify any remaining serialization in orchestration or scorer access
- only then tune task scheduling / minimum spectra per thread / map ownership

## 6. Expected payoff by phase

These are directional planning numbers, not commitments:

| Phase | Astral wall impact | Recall risk | Notes |
|---|---:|---|---|
| A: deisotope + peak cap | 1.15-1.35× | low-medium | likely helps sensitivity if deisotoping is correct |
| B: calibrated window tightening | 1.15-1.30× | medium | must be heavily recall-gated |
| C: branch-and-bound SA walk | 1.5-2.5× | medium-high | only if the bound is both safe and meaningfully selective |
| D: GF threshold tightening | 1.05-1.15× | low-medium | follow-on effect after C |
| E: better scaling | 1.2-1.8× | low | depends on new post-pruning profile |

Compounded, this is the first roadmap that can plausibly reach **3.5× to 6×**.

The dominant uncertainty is Phase C.

## 7. Telemetry we must add before betting on this

Before major coding, add instrumentation that can run on TMT and Astral:

### Search-space telemetry

- candidate peptide variants considered per SA index
- matched `SpecKey` count per candidate peptide
- cheap-score calls per spectrum
- top-N threshold evolution per spectrum
- precursor-mass index span per spectrum in GF

### Pruning telemetry

- branches considered
- branches pruned by mass gate
- branches pruned by bound gate
- retained branches that produce at least one final top-N match
- false-alarm audit on debug runs:
  - prefixes that would have been pruned
  - whether any descendant became a final retained match

### Spectrum-shape telemetry

- peaks before and after deisotoping
- peaks before and after dense-peak cap
- calibrated precursor-window widths

This telemetry should be written behind a debug flag, not always-on.

## 8. Acceptance and kill gates

This roadmap needs hard stop conditions.

### Phase A gates

- Astral wall improves measurably
- Astral 1% FDR PSMs do not regress below 35 600
- PXD001819 remains within existing gate

Kill:

- if deisotoping reduces Astral recall materially without compensating wall win

### Phase B gates

- precursor-window median width shrinks materially on Astral
- candidate pairing count drops materially
- recall stays within gate

Kill:

- if tightened windows do not meaningfully reduce pairing fan-out

### Phase C gates

- branch pruning removes a large fraction of cheap-score calls
- debug audit shows no exact-bound violations
- Astral wall improves by at least 1.5× over the pre-Phase-C branch baseline

Kill:

- if the safe bound is too weak to prune enough work
- if the bound becomes heuristic and starts threatening recall
- if implementation state balloons memory beyond the 8 GB target

### Whole-roadmap gate

Proceed only while the compounded measured gain is tracking toward at least **3×** by the time Phase C is working. If A+B+C together cannot plausibly clear 3×, stop and reassess instead of polishing a dead branch.

## 9. Proposed implementation order

### Iteration 0 — telemetry-only branch

Goal:

- quantify where Astral fan-out really happens on `dev`

Touches:

- `DBScanner`
- `ScoredSpectraMap`
- optional debug output helpers

### Iteration 1 — deisotoping + peak-cap scaffold

Goal:

- validate that spectrum cleanup helps candidate density and cheap-score separation

Touches:

- `Spectrum`
- scorer preprocessing path
- tests with synthetic isotope clusters

### Iteration 2 — calibrated window tightening

Goal:

- reduce precursor pairing width and GF mass-index span

Touches:

- `MassCalibrator`
- `ScoredSpectraMap`
- `SearchParams`
- `DBScanner`

### Iteration 3 — branch-and-bound prototype

Goal:

- prove that a conservative bound can prune real Astral work

Touches:

- `CandidatePeptideGrid`
- `DBScanner`
- scorer helpers
- new pruning-state classes

Deliverable:

- prototype guarded by an OFF-by-default flag

### Iteration 4 — exactness audit + optimization

Goal:

- prove correctness and reduce overhead of the pruning machinery itself

This is where we decide whether the branch becomes the main path or gets abandoned.

### Iteration 5 — GF tightening and scaling follow-up

Goal:

- exploit the smaller retained candidate set

Touches:

- `PrimitiveGeneratingFunction`
- `DBScanner.computeSpecEValue`
- orchestration / task sizing if needed

## 10. What I would not do next

- **Do not re-open the fragment-index branch.**
  The post-mortem is still right: too much Tier-1 cost, too much memory, too much architectural risk.

- **Do not start with another GF-local optimization.**
  Useful later, but it does not solve the multiplicative fan-out earlier in the search.

- **Do not start with a concurrency rewrite.**
  That risks parallelizing waste before we have shrunk the search space.

## 11. My recommendation

> **Update 2026-04-28: Phase A was attempted and reverted.** See [`astral-phase-a-retrospective.md`](astral-phase-a-retrospective.md) for measurements, lessons, and what's still untried. Three independent angles (deisotope+cap, GF candidate cap, scorer hot-path) all failed the Astral wall gate. TMT-as-inner-loop turned out unsafe — TMT's 1.41× win did not transfer to Astral. The 5× roadmap below is preserved for future agents but the strategy of "start with Phase A" is now disproven; future attempts should pick Phase B, C, or E and re-profile before betting on micro-optimizations.

Phase B (calibrated precursor-window tightening) and Phase E (parallelism ceiling) are the remaining lower-risk shots. Phase C (branch-and-bound) is the highest-variance / highest-upside option but needs the upfront design work the retrospective flags. Phase D is unlikely to be useful as a standalone lever on Astral given the GF candidate cap measurements.

Original recommendation, preserved for context:

> Try **Phase A first** as the opening big-win attempt:
> 1. telemetry milestone commit (Iteration 0)
> 2. spectrum cleanup milestone commit (Iteration 1, Phase A)
>
> If Phase A delivers, continue with Phase B then Phase C as further milestone commits on the same branch. If Phase A fails its kill gate (no measurable wall win and no recall upside), drop to Iteration 0.5 fallback (memoization + GF candidate cap; see §0) and ship those as the iteration's deliverable.
>
> Phase C is the centerpiece of 5× but the highest-variance phase; do not attempt it before Phase A is in place because cleaner spectra make C's upper bounds tighter.

## 12. Reference

- **Phase A attempt retrospective (read first):** [`astral-phase-a-retrospective.md`](astral-phase-a-retrospective.md)
- Iteration retrospective: [SHIPPED.md](/Users/yperez/work/msgfplus-workspace/astral-speed/.claude/plans/SHIPPED.md:1)
- Benchmark summary: `~/.claude/plans/benchmarks/3engine-results.md`
- Fragment-index post-mortem: `~/.claude/plans/msgfplus-fragment-index/ABANDONED-2026-04-20.md`
- Historical Astral profile: `~/.claude/plans/msgfplus-primitives-optimization/profile-astral.md`
- Earlier short-horizon plan (superseded; consolidated into §0 fallback): recoverable via `git show 878b0cb:.claude/plans/astral-speed-improvements.md`
