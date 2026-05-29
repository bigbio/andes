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

## Design (two-layer discrimination: residual SpecEValue in-engine + Percolator on all features)

**No hand-tuned filter, no parameter.** A magic-number cutoff
(`--chimeric-min-unique-ions`) was rejected — it won't generalize across
instruments/datasets. Instead the discrimination is the *score itself*, at two
complementary layers:

1. **In-engine — residual SpecEValue.** Greedy peak-claiming (most-confident
   first) re-scores each rank≥2 peptide's RawScore + GF **SpecEValue** on the peaks
   it can *uniquely* claim (those not already taken by a more confident peptide). A
   theft/coincidental peptide, stripped of stolen/coincidental peaks, gets a *bad*
   residual SpecEValue → poor q-value → falls out of the 1% set on its own. A real
   co-isolated peptide retains genuine signal → stays significant. SpecEValue is
   MS-GF+'s native length/quality-robust chance model — exactly the "good vs bad"
   logic, applied symmetrically to targets and decoys.

2. **Percolator — all features.** Percolator receives the full vector: the
   re-scored RawScore + lnSpecEValue *and* the new additive unique-evidence columns
   (`UniqueMatchedIons`, `UniqueExplainedFraction`, `SharedFracClaimed`), and learns
   the final boundary / computes the q-value over the (now-credible) set.

This is the synthesis of the earlier "Approach A (additive columns)" and "Approach
B (re-score SpecE)": we re-score SpecEValue (B — the piece that actually deflates,
which Phase 2's *soft-feature-only* path lacked) **and** emit the unique-evidence
columns (A — orthogonal Percolator signal). The earlier "Rule-2 regression" caveat
on re-scoring does **not** apply: `--chimeric off` is bit-identical and rank-1
claims first (its peaks are never reduced), so the *only* scores that change are the
chimeric rank≥2 extra rows — which exist in no shipped baseline and which we
explicitly intend to deflate.

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
ranks 1..k-1 have claimed, compute the unique-evidence metrics:
- `UniqueMatchedIons` = |matched peaks ∉ claimed|
- `SharedFracClaimed` = |matched ∩ claimed| / |matched|   (0 if |matched|=0)
- `UniqueExplainedFraction` = Σ intensity(unique peaks) / Σ intensity(all matched
  peaks of this PSM)   (0 if denominator 0)

Peak keys are the charge-1 b/y keys `matched_peak_keys` already produces
(identical to the validated diagnostic, so the design is grounded in the measured
overlap distribution). After computing its metrics, the peptide inserts ALL its
matched peaks into the per-scan `claimed` set (so the next peptide sees them as
taken).

### §3 Residual re-score + additive columns + decoy symmetry

- **Residual SpecEValue re-score (the discriminator):** for each rank≥2 PSM, build
  a peak-reduced view of the spectrum with the `claimed` peaks removed, recompute
  the peptide's RawScore (`score_psm` + cleavage + `psm_edge_score`) and its GF
  SpecEValue (`compute_spec_e_values_for_spectrum` on a one-PSM queue) against that
  residual spectrum, and overwrite the PSM's `score`/`rank_score`/`spec_e_value`.
  No hard filter, no threshold — a bad residual SpecEValue simply fails to clear
  the FDR cutoff downstream.
- **Rank-1 untouched:** rank-1 claims first; its residual = full spectrum → its
  scores are unchanged. (Implementation may skip the re-score for rank-1 entirely.)
- **Decoy symmetry (critical):** the re-score is applied identically to target and
  decoy rank≥2 rows. Restoring the inflated decoy fraction toward baseline is the
  mechanism that makes the FDR credible; asymmetric treatment would bias it.
- **Additive PIN columns (chimeric-on schema only):** `UniqueMatchedIons`,
  `UniqueExplainedFraction`, `SharedFracClaimed`, emitted for all surviving rows as
  extra Percolator signal. These join the existing chimeric-on extra columns
  (PrecursorIsotopeKL); the `--chimeric off` schema is untouched.
- **Back-end-window caveat (known limitation):** a genuinely off-precursor
  co-isolated peptide whose nominal mass falls outside the GF mass window already
  gets `spec_e_value = 1.0` (single-precursor-centered back-end, per the Phase-2
  addendum). The residual re-score does not change that; fixing the GF window is out
  of scope here.

### §4 Invariants, gating, success criteria

- **Invariant:** `--chimeric off` → Phase 3 unreachable → PIN bit-identical. Rank-1
  scores/features never change (even chimeric-on).
- **Primary gate (canary first):** the Astral count must return toward ~36.7k
  *before* trusting any PXD delta. Track the target/decoy ratio as the leading
  trustworthiness indicator (the residual SpecEValue should deflate spurious
  near-precursor co-emissions).
- **Win condition:** PXD > 14,974 @1% FDR with restored T/D balance; wall within
  ~3% of current chimeric-on.
- **Speed:** residual re-score runs only for rank≥2 PSMs on scans with ≥2 distinct
  emitted peptides (a minority); each is one ScoredSpectrum rebuild + one GF group.
  Spend the Astral speed surplus; measure; if it blows the ~3% budget, fall back to
  re-scoring only the residual node-score sum (skip the full GF rebuild).

### §5 Testing (TDD)

1. Off-path bit-identity (extend the existing parity/off-mode test).
2. Subset theft (unique-metrics core): runner-up B's peaks ⊂ A's → after A claims,
   B has 0 unique, `SharedFracClaimed = 1.0`, `UniqueExplainedFraction = 0`.
3. Disjoint survivor (unique-metrics core): runner-up C with disjoint peaks →
   `SharedFracClaimed ≈ 0`, `UniqueMatchedIons` = its full count,
   `UniqueExplainedFraction` high.
4. Partial overlap: mixed shared/unique → fractions computed correctly.
5. Rank-1 protection: rank-1 metrics = full-spectrum (claims first), scores
   unchanged.
6. Integration smoke (under `--chimeric`): a synthetic 2-peptide scan where B
   steals A's peaks → B's residual SpecEValue is worse than its full-spectrum
   SpecEValue (the deflation actually happens). Decoy and target runner-ups with
   identical evidence get identical residual scores (symmetry).

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Residual SpecE re-score over-deflates real co-IDs | Canary tuning measures net effect; `UniqueExplainedFraction` gives Percolator a counter-signal; back-end-window caveat noted |
| Re-score asymmetry biases FDR | Apply identically to target+decoy; assert via test 6 + bench T/D ratio |
| Wall regression (GF rebuild per rank≥2 PSM) | Scoped to ≥2-distinct-peptide scans; spend Astral surplus; measure; fallback to node-score-only residual |
| Rule-2 regression on shipped path | `--chimeric off` bit-identical; rank-1 untouched; only chimeric extra rows change |
| `--chimeric off` drift | Path guarded on `params.chimeric`; bit-identity test |

## Out of scope

- Per-scan/peptide-level FDR model change (the residual SpecEValue re-score + the
  existing PSM-level FDR is the chosen mechanism; a structural FDR change is a
  separate effort).
- TMT (the chimeric path was net-negative on TMT; not a target here).
- Fixing the single-precursor-centered GF mass window (back-end-window caveat);
  fragment-index speed enabler; DL rescoring; DIA.

## References

- `2026-05-28-chimeric-fragment-overlap-diagnostic.md` — Astral overlap result
  (theft confirmed, bimodal).
- `2026-05-28-chimeric-phase2-bench.md` — why soft features fail; the three
  requirements for trustworthy chimeric.
- `2026-05-28-chimeric-dda-plus-integration-design.md` — Phase 1/2/3 overview.
- MSFragger-DDA+ (Nat Commun 2025); Yu et al. 2023 (greedy shared-fragment method).
