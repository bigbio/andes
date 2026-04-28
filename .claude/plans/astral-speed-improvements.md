# Astral Speed Improvements — Design Doc

**Status:** Design / awaiting approach selection
**Branch:** `feat/astral-speed-improvements` (off `dev` @ `2216bbb`, post-PR-#25)
**Date:** 2026-04-27

## 1. Why this exists

We just merged PR #23 (`feat/msgfplus-speed-v2`), PR #24 (`feature/improve-mzid-suffix-big-fasta`), and PR #25 (`perf/search-sync-cleanup`) into `dev`. Those landed Achievements A+B (pin features + precursor calibration), parallel BuildSA bucket sort, mzML parser MS-level preload filter, scorer Hashtable→HashMap, per-task search infrastructure, and a pile of dead-code cleanup. MS-GF+ is in a clean state.

The next iteration targets **Astral wall-time and memory**. Astral (ProteoBench Module 8: Orbitrap Astral, 32 MB FASTA, 50 K spectra, 10 ppm precursor / 20 ppm fragment) is where MS-GF+ trails Sage most visibly:

- **Wall:** MS-GF+ ~620 s vs Sage 78 s (**7.9× gap**)
- **Memory:** MS-GF+ 7.6 GB peak RSS vs Sage 3.4 GB (**2.2× gap**)
- **Sensitivity:** MS-GF+ 35 627 PSMs vs Sage 32 074 PSMs at 1 % FDR (**MS-GF+ wins +11.1 %**)

Sensitivity is our moat. Speed/memory is the gap we need to narrow without sacrificing it.

## 2. What is *not* in scope

- **Fragment-index** (Sage-style inverted index). Abandoned 2026-04-20 after failing speed/recall/memory gates on Astral; see `~/.claude/plans/msgfplus-fragment-index/ABANDONED-2026-04-20.md`. Do not revisit without new evidence.
- **Big-FASTA scalability for metaproteomics / proteogenomics**. Different problem (peptide redundancy, variant DBs). May share infrastructure later but is not this iteration.
- **PXD001819 / TMT speed work** as a primary goal. Wins on those datasets are welcome side-effects; regressions there are not blockers unless they exceed gate constraints.
- **Cross-engine parity work** (deisotoping comparisons, candidate-gap analyses against Sage). Future iteration.

## 3. Success gate

Adopted from brainstorming session, gate B (moderate, single mergeable PR):

| Metric | Target | Measurement |
|---|---|---|
| Astral wall (clean idle box, 4 threads, 8 GB Xmx) | **≤ 460 s** (≥ 1.35× speedup vs 620 s baseline) | `/usr/bin/time -l` on dev-tip vs branch head |
| Astral peak RSS | **≤ 7.6 GB** | same |
| Astral 1 % FDR PSMs (Percolator-rescored) | **≥ 35 600** (no regression vs 35 627) | `compare_metrics.py` integration test |
| PXD001819 1 % FDR PSMs | **≥ 15 100** (no regression vs 15 157) | CI benchmark |
| TMT PXD007683 1 % FDR PSMs | **≥ 10 100** (no regression vs 10 176) | manual run |
| Bit-identical OFF-mode behaviour | **required** for any new flag | unit + integration tests |

Stretch (not gating, but tracked): RSS reduction toward Sage's 3.4 GB; ScoreDist allocation rate.

### 3.1 Two-tier benchmark cadence

Astral runs are too slow (~10 min/run) for fast iteration. We split the workload:

- **Inner loop — TMT PXD007683** (~321 s current-dev wall, ~50 K spectra, 17 MB FASTA, Lumos high-res MS2). Used during day-to-day development for measuring wall-time deltas, RSS deltas, and engine-internal target/decoy counts. Closest available analog to Astral candidate-density dynamics; faster turnaround.
- **Phase gate — Astral ProteoBench Module 8** (~620 s baseline). Run only at end-of-phase (Phase 1 acceptance, Phase 2 acceptance, Phase 3 final) to confirm wall/RSS/PSM gates §3 hold on the actual target dataset.
- **Smoke baseline — PXD001819** (~96 s baseline). Run alongside TMT on every iteration for "no regression on small-FASTA" sanity (CI benchmark already automates this).

**Caveat on TMT signal reliability:**

| TMT signal | Reliable? | Notes |
|---|---|---|
| Wall-time delta | ✓ | engine-internal; no decoy-pool dependency |
| Peak RSS delta | ✓ | engine-internal |
| Native target / decoy counts | ✓ | engine-internal |
| Percolator 1 % FDR PSM count | ⚠️ | the documented Sage decoy-pool artefact (§Astral conclusions in `3engine-results.md`) means TMT decoy pools are less calibrated than Astral's. Use as a directional indicator, not a hard recall gate |

Recall regression decisions are made on Astral, not TMT.

## 4. Anchor data

Numbers used to size the approaches below come from:

- `~/.claude/plans/benchmarks/3engine-results.md` — clean Astral baseline, 3-engine table
- `~/.claude/plans/msgfplus-primitives-optimization/profile-astral.md` — pre-Hashtable-fix profile (2026-04-17). Granular CPU/alloc breakdown is **stale on dev-tip** (the 43 % Hashtable contention it identified is gone after `8442f2c`). Macro wall/RSS still anchor; per-method % needs a fresh run before deep tuning.
- `~/.claude/plans/msgfplus-fragment-index/ABANDONED-2026-04-20.md` — five follow-up speed ideas with expected ROI/risk

A re-profile on dev-tip is **Phase 0** of any combo we pick (see §6).

## 5. Catalog of approaches considered

Six well-formed approaches plus smaller levers were considered. Each is sized for gate B (single mergeable PR) and assessed for Astral wall/RSS/recall and PXD001819 side-effects.

### Approach 1 — Graph-skeleton memoization

Cache the read-only arrays of `PrimitiveAminoAcidGraph` (`reachable`, `inEdgeCountByMass`, `activeNodes`, `massToNodeIdx`) keyed on `(peptideMassIndex, enzyme, aaSetHash, useProtNTerm, useProtCTerm)`. Per-spectrum scoring fields stay per-graph. Cache lives on `DBScanner` (per-task). Single construction site at `DBScanner.java:620` inside the per-mass-index inner loop.

| Dimension | Value |
|---|---|
| Astral wall | ~10–15 % |
| RSS | +50–150 MB per task |
| Recall risk | None (pure memoization) |
| PXD001819 | +5–7 % |
| Blast radius | `PrimitiveAminoAcidGraph` + new cache class + tests |
| Standalone hits gate B? | No (1.15× < 1.3×) |

### Approach 2 — Adaptive precursor-tolerance tightening

After `MassCalibrator` finishes its second pass, use the learned ppm shift mean/σ to narrow the effective search window: `effective_tol = min(user_tol, calibrated_shift + k·σ)` with `k=3`, default-on for AUTO mode, opt-out via `-precursorCal aggressive|conservative|off`. Astral typical: 10 ppm window collapses to ~3 ppm post-calibration → ~2–3× fewer candidates per spectrum.

| Dimension | Value |
|---|---|
| Astral wall | ~20–30 % |
| RSS | slight reduction |
| Recall risk | Real but bounded — k=3 keeps a 3-σ envelope. Mitigation: integration test enforces ≥ 99.5 % PSM recall vs OFF on Astral and PXD001819 |
| PXD001819 | smaller % win (Velos σ wider) |
| Blast radius | `MassCalibrator` + `SearchParams.getEffective*Tolerance()` + `DBScanner` + tests |
| Standalone hits gate B? | Yes, with margin (1.25–1.4×) |

### Approach 3 — Parallelism-ceiling investigation

The 2026-04-17 profile ran 4 threads on 11 cores; 7 cores idle. Winkelhardt-paper parity suggests MS-GF+ caps at 4–6 effective cores. Phase 0 measures dev-tip 1→4→8-thread scaling on Astral. If linear, drop this approach. If it caps, root-cause via per-task wall stats (already plumbed) and remove the bottleneck.

| Dimension | Value |
|---|---|
| Astral wall | 0–50 % (high variance; depends on what we find) |
| RSS | +20–30 % (more in-flight tasks) — could push past gate |
| Recall risk | Concurrency bugs only |
| PXD001819 | similar story, smaller absolute |
| Blast radius | `ConcurrentMSGFPlus` + possibly `MSGFPlus.runMSGFPlus` + thread-safety audit |
| Standalone hits gate B? | Unknown — research-shaped, doesn't fit "single mergeable PR" if rewrite is needed |

### Approach 4 — In-engine MS2 deisotoping

Collapse `(M, M+1, M+2…)` isotope clusters in MS2 into the monoisotope before scoring. Sage / MaxQuant / Comet all do this; MS-GF+ trusts the mzML peak list. The 3-engine analysis already identified this as the dominant cause of the Astral candidate-generation gap with Sage — so this is recall-positive AND speed-positive.

| Dimension | Value |
|---|---|
| Astral wall | ~15–25 % (peak count drop ~3×; cheap-score sublinear) |
| RSS | slight reduction |
| Recall risk | None expected (established prior art); risk is implementation bugs |
| PXD001819 | smaller win (Velos lower-res, partial natural merging) |
| FDR sensitivity | **+ on Astral** (closes candidate-gen gap) |
| Blast radius | New `Deisotoper` + `Spectrum.deisotope()` hook + `-deisotopeMS2 on\|off` flag + tests |
| Standalone hits gate B? | Borderline on speed; clearly hits if "+PSMs at same wall" is also a win |

### Approach 5 — Tier-1.5 candidate cap before GF

Today, every match surviving cheap-score top-K reaches `PrimitiveGeneratingFunction`. Tighten the cap two ways:

1. **Hard cap**: `numCandidatesForGF` (default e.g. 10) — only top-N by cheap score reach GF.
2. **Score-gap pruning**: skip GF for candidates whose cheap score is more than Δ below the per-spectrum top score (Δ tunable, gate-tested).

Most aligned with the original "small preprocessing of candidates" framing.

| Dimension | Value |
|---|---|
| Astral wall | ~15–30 % (GF was ~60 % of CPU; cutting 30 % of GF inputs ≈ 18 % wall) |
| RSS | slight reduction |
| Recall risk | Real and quantifiable. Mitigation: integration test asserts no PSMs at 1 % FDR rank below the new cap; if any, raise cap |
| PXD001819 | smaller win (smaller pool to begin with) |
| Blast radius | `DBScanner.computeSpecEValue` (~25 lines around line 600) + sort + `SearchParams.numCandidatesForGF` + tests |
| Standalone hits gate B? | Possibly; comfortable when paired with Approach 1 |

### Approach 6 — Astral-tuned NewRankScorer parameter file

`NewRankScorer`'s rank-distribution / ion-existence tables are trained on Velos-era data. Astral's peak quality, b/y ratios, and fragment-error distributions differ. Retrain on a clean Astral PSM corpus (use current MS-GF+ 1 % FDR PSMs as labels; existing training pipeline supports this) and ship `Astral_*.param` with auto-detect via mzML instrument metadata.

| Dimension | Value |
|---|---|
| Astral wall | ~5–15 % + likely +1–3 % FDR sensitivity |
| RSS | none (data file swap) |
| Recall risk | Minimal — auto-detect + Velos fallback |
| PXD001819 | none (different param file selected) |
| Blast radius | Training script (offline; majority of the work) + auto-detect logic + new .param data file |
| Standalone hits gate B? | No — force multiplier for Approaches 1, 2, 5 |

### Smaller levers

Folded into the chosen approach as nice-to-haves, or kept as Phase 2 follow-ups:

- **GF reuse across same-mass candidates within a single spectrum.** Same nominal mass + same spectrum = identical score distribution; currently recomputed. Tiny code change, ~3–5 % wall.
- **Top-N peak retention for dense MS2.** Cap peaks per spectrum at e.g. 200 highest-intensity. Distinct from deisotoping. ~5 % wall on Astral; needs recall test.
- **`PrimitiveGeneratingFunction` early termination.** Abort when partial score distribution proves SpecEValue is far above the rank-1 threshold. Algorithmic; needs correctness proof. ~5–10 % wall.
- **Vector API in `NewRankScorer.getScore` peak-intersection loop.** Hotter than `ScoreDist.addProbDist`. High variance, JVM-version-sensitive.
- **Charge-state pre-filter on Astral.** Astral reports charge cleanly; trust it more aggressively. Tiny win, near-zero risk.

## 6. Recommended combination

The combinations evaluated for gate B:

| Combo | Astral wall projection | RSS | Sensitivity | PR size |
|---|---|---|---|---|
| **1 + 5** (memo + GF cap) | 1.3–1.5× | ≤ baseline | flat (recall-gated) | Medium |
| **1 + 4** (memo + deisotoping) | 1.25–1.4× | ≤ baseline | **+** (positive) | Medium |
| **2 + 4** (tolerance + deisotoping) | 1.5–1.8× | ≤ baseline | **+** | Larger |
| **4 + 5 + 6** ("Astral pack") | 1.4–1.7× | ≤ baseline | **+** | Larger |

**Primary recommendation: Approach 1 + Approach 5** (graph-skeleton memoization + Tier-1.5 candidate cap before GF).

Rationale:

1. **Two well-known, well-bounded levers.** Each has a single hot site in `DBScanner`, a clear test surface, and zero ambiguity about the cache/cap mechanism.
2. **Layered correctness.** Memoization is provably equivalent (same arrays, same content). The cap has a recall test that fails CI if it would drop a current 1 %-FDR PSM.
3. **Independent commits.** If Approach 5 fails its recall test at any cap value, ship Approach 1 alone — still a measurable Astral improvement, no rework.
4. **Smallest-blast-radius combo that hits gate B.** Touches `DBScanner` + `PrimitiveAminoAcidGraph` + a new cache class + a new `numCandidatesForGF` knob. Reviewable as one PR.
5. **Clear next-iteration runway.** Approach 4 (deisotoping), Approach 2 (adaptive tolerance), and Approach 6 (Astral-tuned scorer) are all natural follow-ups that compose cleanly with this PR's work.

Alternative if sensitivity is a higher priority than raw speed: **Approach 1 + Approach 4** — ships +PSMs alongside ~25 % wall improvement.

## 7. Implementation phases (for the recommended 1+5 combo)

Phases here are scoped at the design level, not as commits. Detailed plan goes to `superpowers:writing-plans` after this design is approved.

### Phase 0 — Re-measure on dev-tip (1 commit, no production code change)

Profile on **both** TMT (inner-loop reference) and Astral (phase-gate reference) so subsequent phases can compare TMT wins against Astral wins and detect divergence early.

- Run async-profiler + JFR on `dev` HEAD with **TMT** (4 threads, 8 GB Xmx, full run, 120 s steady-state CPU + 120 s alloc).
- Run the same profile on **Astral** (same threads/Xmx, 180 s windows).
- Record top-30 self-time methods, top-20 alloc sites, GC summary on each.
- Confirm `PrimitiveAminoAcidGraph.<init>` is still a measurable line item on at least Astral (post-mortem said 7.4 % on PXD001819 at the time; expected higher on Astral).
- Confirm `PrimitiveGeneratingFunction.computeGeneratingFunction` + `ScoreDist.addProbDist` are still in the top 5 on both.
- Compute and record the wall-time ratio `TMT_wall / Astral_wall` on dev-tip — used in Phase 1 / 2 to sanity-check that TMT speedup translates roughly to Astral speedup.
- Save artifacts under `~/.claude/plans/astral-speed-improvements/profile-2026-04-XX/`.

**Gate to proceed:** if the profile shows a different dominant cost (e.g. a new bottleneck introduced by PR #23/#24/#25) **or** TMT's hot-spot ranking diverges materially from Astral's (a sign that TMT-as-inner-loop will mislead), pause and either re-rank approaches or pick a different inner-loop dataset before coding.

### Phase 1 — Graph-skeleton memoization (Approach 1)

1. Add `GraphSkeletonCache` keyed on `(peptideMassIndex, enzymeId, aaSetVersion, ntermFlag, ctermFlag)`.
   - Cached value: the four immutable arrays (`reachable`, `inEdgeCountByMass`, `activeNodes`, `massToNodeIdx`) packaged as a small record.
   - Per-task instance, no cross-thread sharing (preserves the post-PR-#25 lock-free hot path).
   - Bounded by an LRU with a generous default (e.g., 4 096 entries — covers Astral's ~3 000 distinct nominal masses with headroom).
2. Refactor `PrimitiveAminoAcidGraph` constructor to:
   - Accept a pre-built skeleton (new ctor) **or** build one from scratch (existing ctor — kept for tests + as fallback).
   - Apply per-spectrum scoring fields after attachment.
3. Update the `DBScanner.java:620` site to consult the cache, falling back to direct construction on cache miss.
4. Tests:
   - **Unit:** cache hit returns object equal to from-scratch build (deep array equality).
   - **Unit:** cache miss populates correctly.
   - **Integration:** PXD001819 CI benchmark — bit-identical native target counts vs baseline.

**Iteration cadence:** measure each tuning change on TMT (~5 min/run). Run Astral exactly once at end of phase, before merging the phase commit.

**Acceptance:**
- Inner-loop (TMT): wall ↓ ≥ 5 %; native target/decoy counts bit-identical to dev-tip OFF-mode.
- Phase gate (Astral): wall ↓ ≥ 8 % vs Phase 0 measured baseline; PXD001819 native-T count bit-identical (CI benchmark).

### Phase 2 — Tier-1.5 candidate cap before GF (Approach 5)

1. Introduce `SearchParams.numCandidatesForGF` (default `Integer.MAX_VALUE` = current behaviour) and `SearchParams.gfScoreGapPrune` (default disabled).
2. In `DBScanner.computeSpecEValue` (around line 600), before the `for (DatabaseMatch match : matchQueue)` loop:
   - Sort `matchQueue` by cheap score descending.
   - Truncate to `numCandidatesForGF`.
   - If `gfScoreGapPrune` is set, drop entries whose cheap score is below `topScore - gap`.
3. Set defaults conservatively for the released config (e.g. `numCandidatesForGF=20`, gap pruning off) — values must clear the recall test.
4. Tests:
   - **Unit:** `computeSpecEValue` with cap=2 produces SpecEValues for exactly the top 2 cheap-scored matches; remainder marked as filtered.
   - **Integration:** Astral 1 % FDR PSM count ≥ 35 600 at the chosen default cap. PXD001819 ≥ 15 100. TMT ≥ 10 100.
   - **Recall regression:** scan dev-tip OFF-mode pin and verify every 1 %-FDR PSM survives the cap on the same data.

**Iteration cadence:** sweep cap values (e.g. 5, 10, 20, 50) on TMT to find the wall-vs-recall knee; pick the cap value, then run Astral once to confirm it holds.

**Cap-tuning loop (TMT inner-loop):**
1. Run with `numCandidatesForGF=5`, record TMT wall, native targets, native decoys.
2. Repeat with cap = 10, 20, 50.
3. Pick the smallest cap whose native-target count is within 0.2 % of OFF-mode on TMT.
4. Run Astral once at the chosen cap; confirm Astral 1 % FDR PSMs ≥ 35 600.
5. If Astral fails, increase cap one tier (e.g. 10 → 20) and re-run Astral.

**Acceptance:** all three benchmark FDR counts within gate (§3); Astral wall ≤ 460 s on the Phase 0 reference machine, combined with Phase 1.

### Phase 3 — Final benchmark + docs

1. Run full 3-engine matrix (PXD001819, Astral, TMT) on branch HEAD; commit results to `docs/benchmarks/`.
2. Update `docs/changelog.md` with the gate-B numbers.
3. Document the two new flags in `docs/msgfplus.md`.

## 8. Verification strategy

- **Bit-identical OFF-mode.** Both new behaviours behind their flags; default cap value uses `Integer.MAX_VALUE` for the truly-OFF mode tested by integration. (We will only switch the *production* default cap after the recall test demonstrates safety on all three benchmarks.)
- **TMT inner-loop benchmark:** local script that runs TMT with feature ON and OFF, records wall + RSS + native target/decoy counts. Run on every meaningful code change. Not in CI (TMT data is not staged for CI runners).
- **PXD001819 CI benchmark:** existing `benchmark/ci/PXD001819/run_ci.sh` extended with a "feature ON" run; comparator gates 1 % / 5 % FDR counts. Runs on every push.
- **Astral phase-gate:** scripted end-to-end on the existing Astral dataset; results fed to `compare_metrics.py` against the baseline.tsv updated with Astral baseline. Run at end of Phase 1, end of Phase 2, end of Phase 3 — not on every code change.
- **Unit tests:** per-class for cache and cap logic. Run on every push.
- **Profile re-confirm at end:** async-profiler shows `PrimitiveAminoAcidGraph.<init>` and GF self-time both reduced on TMT and Astral relative to Phase 0 baseline.

## 9. Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Phase 0 reveals a different dominant cost (e.g. PR #25 introduced a new hot spot) | Medium | Re-rank approaches before coding; don't proceed on stale assumptions |
| TMT-inner-loop wins fail to translate to Astral (different precursor tolerance: 20 ppm vs 10 ppm; different mod density) | Medium | Phase 0 records the TMT/Astral wall ratio; phase gates check Astral explicitly. If divergence emerges, fall back to running Astral at higher frequency for that specific change |
| Approach 5 cap drops real PSMs | Medium | Recall integration test; conservative default; fallback knob |
| Memoization correctness — graph skeleton silently differs from from-scratch build | Low | Unit-level deep equality test; OFF-mode bit-identical integration test |
| Astral wins on a feature flag but PXD001819 regresses | Low (memoization is dataset-agnostic) | CI benchmark gates regression on PXD001819 |
| Memory bloat from cache | Low | LRU bound; size monitored in test |
| Sensitivity drops below MS-GF+'s lead over Sage | Low (Approach 5 is recall-gated; Approach 1 is recall-neutral) | Same gates as above |

## 10. Open questions / decisions for ypriverol

1. **Approach selection.** Confirm Approach 1 + Approach 5, or pick a different combo from §6. If Approach 4 (deisotoping) appeals more for the sensitivity boost, we can swap.
2. **Datasets staged.** Do we have TMT PXD007683 mzML + FASTA staged for fast iteration on the dev box? Do we have the dev-tip Astral mzML + FASTA staged on a CI-equivalent box for the Phase 0 re-profile and phase-gate benchmarks? If either is missing, that's a prerequisite step.
3. **Default cap value.** I've sketched `numCandidatesForGF=20` as a safe-feeling default. This needs to be picked from the TMT cap-sweep + Astral confirmation, not chosen up-front. Approval to leave it TBD until Phase 2 measurement?
4. **Approach 6 (scorer retraining)** as a follow-up iteration — should it be tracked here as a "next-after-this-PR" item, or kept entirely separate? It composes cleanly with Approach 5 if we do it later.
5. **Worktree path.** I created the worktree at `~/work/msgfplus-workspace/astral-speed`. Confirm or move.

## 11. Reference

- Abandoned fragment-index post-mortem: `~/.claude/plans/msgfplus-fragment-index/ABANDONED-2026-04-20.md`
- Stale Astral profile (pre-Hashtable-fix; macro numbers still valid): `~/.claude/plans/msgfplus-primitives-optimization/profile-astral.md`
- 3-engine benchmark: `~/.claude/plans/benchmarks/3engine-results.md`
- Existing perf-PR plans for reference style: `.claude/plans/parameter-modernization.md`, `.claude/plans/search-sync-cleanup.md`
