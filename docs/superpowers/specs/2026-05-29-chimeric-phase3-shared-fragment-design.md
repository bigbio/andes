# Chimeric Phase 3 — greedy shared-fragment competition + unique-evidence filter

**Date:** 2026-05-29
**Branch:** `feat/chimeric-dda-plus`
**Type:** Feature design (single bench-gated sub-project)
**Supersedes:** Phase 3 sketch in `2026-05-28-chimeric-dda-plus-integration-design.md` §"Phase 3"

## Motivation

The chimeric front-end (Phase 1) emits multiple distinct peptides per scan by
searching the full isolation window. PSM-level target/decoy FDR over ~5 PSMs/scan
is **structurally inflated** — the Astral chimeric count rises +97% (~36.7k →
~72k) even though Astral's narrow isolation windows make real co-isolation rare.
Phase 2 proved a *soft* MS1 isotope feature cannot deflate this (hard-filter test:
72,457 → 69,250, still +89%). The Phase-2 post-mortem identified the true
discriminator as **fragment-level**: spurious co-identified peptides are inflated
because they are scored against MS2 peaks they did not generate.

The **2026-05-29 Astral overlap diagnostic** (commit `a8a7489f`, note
`2026-05-28-chimeric-fragment-overlap-diagnostic.md`) measured this directly on
n=121,423 co-emitting Astral scans: mean fraction-of-smaller 0.367, **38% of
scans ≥0.5 overlap**, with a **bimodal** distribution — a high-overlap "fragment
theft" mode (peaking near-total at [0.9-1.0)=11.4%) **plus** a ~28% coincidental
low-overlap mode. The fragment-theft premise behind Phase 3 is **confirmed** for a
substantial fraction of scans (BSA's low-overlap preview was an artifact of a
single-protein fixture with no co-isolation).

## Goal & definition of done

**Goal:** a *trustworthy* chimeric search — one whose Astral control returns to a
plausible value — that **beats Java on PXD001819** @1% FDR.

**Definition of done (all must hold):**
1. **Canary (trustworthiness):** at the chosen filter threshold, the Astral
   chimeric PSM count returns toward its plausible baseline (~36.7k), not the
   inflated ~72k. Target/decoy balance restored toward ~baseline.
2. **Win:** at that same threshold, **PXD001819 @1% FDR > Java's 14,974**.
3. **Invariant:** `--chimeric off` PIN is bit-identical to current on all 3 datasets.
4. **Speed:** wall within ~3% of current chimeric-on.

**Out of gate-scope (acknowledged):** per `merge-gate-beat-java`, chimeric does
not by itself clear the merge gate (it is a no-op on narrow-window Astral and was
net-negative on TMT). This sub-project is **research toward a trustworthy chimeric
search**, not a merge. Nothing merges to `dev` under this work.

## Design (Approach A — filter + additive columns; no score modification)

Chosen after weighing three approaches:
- **A (chosen):** peak-claiming drives a hard pre-Percolator filter + new additive
  PIN columns; existing scores untouched. Cheapest (reuses `matched_peak_keys`, no
  GF re-run), audit-safe (no Rule-2 distribution change on surviving rows).
- **B (rejected):** full re-score on reduced peaks, emit reduced RawScore + GF
  SpecE. Literature-faithful (Yu 2023) but Rule-2 regression-prone and most wall.
- **C (fallback):** rescore drives the filter but emit additive only. Use only if
  A's count-based filter proves too coarse on the bench.

### §1 Architecture & data flow

New module `crates/search/src/shared_fragment.rs`, called from
`run_chunk_inner` at the point where the per-spectrum merged queue exists and
`fill_post_topn` has populated features (`match_engine.rs:~684`, where the
`chim_overlap` diagnostic already computes the needed primitive). Guarded on
`params.chimeric`; the off path never enters it.

Per scan (chimeric on):
1. Existing path yields emitted PSMs in **confidence order** (rank-1 = best
   SpecEValue — the current merge order).
2. Walk PSMs rank-ascending, maintaining a per-scan `claimed: FxHashSet<PeakKey>`.
3. For each PSM compute its matched peak set via `matched_peak_keys` → split into
   **unique** (∉ claimed) vs **shared** (∈ claimed) → then insert all its peaks
   into `claimed`.
4. **Rank-1 is never filtered and is claimed first** (its peaks seed `claimed`);
   the filter applies only to rank ≥2. This guarantees the single-best-peptide
   result per scan is identical to today.

### §2 Competition algorithm

Greedy, most-confident-first (Yu 2023). For each rank-k PSM (k ≥ 2), after
ranks 1..k-1 have claimed:
- `UniqueMatchedIons` = |matched peaks ∉ claimed|
- `SharedFracClaimed` = |matched ∩ claimed| / |matched|   (0 if |matched|=0)
- `UniqueExplainedFraction` = Σ intensity(unique peaks) / Σ intensity(all matched
  peaks of this PSM)   (0 if denominator 0)

Peak keys are the charge-1 b/y keys `matched_peak_keys` already produces
(identical to the validated diagnostic, so the design is grounded in the measured
overlap distribution).

### §3 Filter + additive columns + decoy symmetry

- **Hard filter (pre-Percolator):** drop any rank ≥2 PSM with
  `UniqueMatchedIons < T`. `T` is a swept CLI knob `--chimeric-min-unique-ions`
  (default chosen by the Astral canary sweep). Theft runner-up → ~0 unique →
  dropped; coincidental weak runner-up → few unique → dropped; real co-isolated
  peptide → substantial unique → kept.
- **Decoy symmetry (critical):** the filter is applied identically to target and
  decoy rank ≥2 rows. Restoring the inflated decoy fraction toward baseline is the
  mechanism that makes the FDR credible; an asymmetric filter would bias it.
- **Additive PIN columns (chimeric-on schema only):** `UniqueMatchedIons`,
  `UniqueExplainedFraction`, `SharedFracClaimed`, emitted for surviving rows.
  Existing column values are **unchanged** (Rule-2-safe). These join the existing
  chimeric-on extra columns (PrecursorIsotopeKL); the `--chimeric off` schema is
  untouched.

### §4 Invariants, gating, success criteria

- **Invariant:** `--chimeric off` → Phase 3 unreachable → PIN bit-identical. Rank-1
  always survives → single-best-peptide result per scan never changes (even
  chimeric-on).
- **Primary gate (canary first):** tune `T` against the Astral count returning
  toward ~36.7k *before* trusting any PXD delta. Track the target/decoy ratio as
  the leading trustworthiness indicator.
- **Win condition:** at the canary-validated `T`, PXD > 14,974 @1% FDR with
  restored T/D balance; wall within ~3%.
- **Speed:** set ops over already-computed matched peaks, no GF re-run → expected
  negligible; measure and hold the budget.

### §5 Testing (TDD)

1. Off-path bit-identity (extend the existing parity/off-mode test).
2. Subset theft: runner-up B's peaks ⊂ A's → after A claims, B has 0 unique → B
   dropped.
3. Disjoint survivor: runner-up C with disjoint strong peaks → survives,
   `SharedFracClaimed ≈ 0`, `UniqueExplainedFraction` high.
4. Threshold boundary: exactly `T` unique → kept; `T−1` → dropped.
5. Decoy symmetry: a decoy runner-up and an identical-evidence target runner-up
   are filtered identically.
6. Rank-1 protection: rank-1 is never dropped (it claims first, is never
   filtered) and its emitted feature values are unchanged from the non-Phase-3 path.

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Count-based filter too coarse (drops real co-IDs / keeps spurious) | Fall back to Approach C (rescore-driven filter); tune `T` on the canary; `UniqueExplainedFraction` available as a secondary criterion |
| Filter asymmetry biases FDR | Apply identically to target+decoy; assert via test 5 + bench T/D ratio |
| Wall regression | No GF re-run; set ops only; measure, hold ~3% |
| Rule-2 regression on surviving rows | Additive-only emission; existing scores untouched |
| `--chimeric off` drift | Path guarded on `params.chimeric`; bit-identity test |

## Out of scope

- Per-scan/peptide-level FDR model change (the unique-evidence hard filter is the
  chosen tail-control mechanism; a structural FDR change is a separate effort).
- TMT (the chimeric path was net-negative on TMT; not a target here).
- Full GF re-score (Approach B); fragment-index speed enabler; DL rescoring; DIA.

## References

- `2026-05-28-chimeric-fragment-overlap-diagnostic.md` — Astral overlap result
  (theft confirmed, bimodal).
- `2026-05-28-chimeric-phase2-bench.md` — why soft features fail; the three
  requirements for trustworthy chimeric.
- `2026-05-28-chimeric-dda-plus-integration-design.md` — Phase 1/2/3 overview.
- MSFragger-DDA+ (Nat Commun 2025); Yu et al. 2023 (greedy shared-fragment method).
