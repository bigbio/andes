# Chimeric search (MSFragger-DDA+ style) integration design

**Date:** 2026-05-28
**Branch:** `feat/chimeric-dda-plus`
**Type:** Phased feature design (large; each phase a separate bench-gated PR)

## Goal

Recover co-fragmented ("chimeric") peptide identifications in msgf-rust by adopting
the MSFragger-DDA+ strategy: search the **full isolation window** of each MS2
(not just ±tol around the selected precursor), emit multiple distinct-peptide
PSMs per spectrum, then refine/rescore them. Target the wide-isolation-window
datasets where co-isolation is common — **PXD001819** (LTQ-Orbitrap Velos, ~2-3 Da
windows) and **TMT** (Lumos) — while remaining a provable no-op on narrow-window
data (Astral) and when the feature is disabled.

Why now: the parity investigation (PR #40, `2026-05-28-phase2-peak-rank-parity.md`)
proved msgf-rust's per-peptide scoring is at/near parity with Java — so further ID
gains require a *new capability*, not parity fixes. Chimeric search is the
highest-evidence algorithmic lever for DDA (+57% peptides reported for DDA+, and
it composes with DL rescoring).

## What MSFragger-DDA+ does (and what we adopt)

DDA+ **reverses** the conventional "detect precursors from MS1, then search" order:

1. **Full-isolation-window database search FIRST** — score the MS2 against every
   peptide whose mass falls in `[selected_mz − lower_offset, selected_mz + upper_offset]`
   (converted to neutral mass per candidate charge). No MS1 needed here.
2. **MS1 targeted-XIC refinement AFTER** — for each candidate PSM, extract the
   isotope XIC at its *theoretical* m/z from the MS1 survey scans and compare the
   observed isotope envelope to the theoretical one via Kullback–Leibler
   divergence; downweight/discard low-quality. (Targeted → avoids untargeted
   deconvolution ambiguity.)
3. **Greedy shared-fragment rescoring** — remove fragment peaks shared by
   co-identified peptides (most-confident first), re-score each (Yu et al. 2023,
   the MSFragger-DIA method).
4. **Percolator (+ optional DL rescoring) → FDR.**

Speed comes from a fragment-ion index (the wide window explodes candidates). We
adopt steps 1–3; DL rescoring is complementary and out of scope here.

## Current msgf-rust shape (what we build on)

- `crates/input/src/mzml.rs`: streaming parser, defaults `ms_level 2..=2`
  (`with_ms_level_range` can widen). Captures the selected precursor m/z but
  **not** the isolation-window width (no `<isolationWindow>` handling).
- `crates/model/src/spectrum.rs`: `Spectrum { precursor_mz, precursor_charge,
  precursor_intensity, ... }` — single precursor, no window bounds.
- `crates/search/src/match_engine.rs::run_chunk_inner`: per spectrum, derives a
  candidate **mass window** (`window_cand_indices` via `bucket_index.range(min..=max)`)
  around the selected precursor ± precursor tolerance ± isotope errors; scores
  candidates into **`per_charge_queues: FxHashMap<u8, TopNQueue>`**; computes
  per-charge GF SpecEValue; merges to one spectrum queue; `fill_post_topn`.
- `crates/search/src/psm.rs`: `TopNQueue`, `PsmMatch`, `PsmFeatures`.
- `crates/output/src/pin.rs`: one PIN row per emitted PSM; `SpecId` per scan.

Two pieces of existing machinery make this tractable:
- The candidate window is already a **mass-range scan** over `bucket_index` — widening
  it to the isolation window is a range-bound change, not new infrastructure.
- Emission is already **multi-queue** (per charge) and Percolator is **PSM-level**,
  so multiple PSMs per scan flow through `.pin` without FDR-model changes.

## Phases

Each phase is a separate, independently shippable, bench-gated PR. The feature is
behind a `--chimeric` flag (default off) until Phase 1 proves out.

### Phase 1 — Full-window search + multi-PSM emission (the core lever)

**1a. Parse isolation-window width.** In `mzml.rs`, add an `IsolationWindow`
parse state and capture cvParams: `MS:1000827` (isolation window target m/z),
`MS:1000828` (lower offset), `MS:1000829` (upper offset). Store on the spectrum
as `isolation_lower_offset: Option<f64>`, `isolation_upper_offset: Option<f64>`
(Da). Fall back to a configurable default width when absent (`--isolation-width`).

**1b. Spectrum model.** Add the two offset fields to `Spectrum`. Default `None`
→ existing behavior.

**1c. Widen candidate enumeration (gated).** In `run_chunk_inner`, when
`params.chimeric` is set, derive the candidate mass window from the **isolation
window** (`selected_mz − lower .. selected_mz + upper`, per charge) instead of
`selected_mz ± precursor_tol`. Reuse the existing `bucket_index.range(...)` scan.
When off → unchanged (bit-identical).

**1d. Multi-distinct-peptide emission.** Today the per-spectrum result collapses
to the top PSMs; the R-2 pepSeq dedup keeps one row per (charge,peptide). For
chimeric, retain the **top-N distinct peptides** across the window (not just the
single best), each with its own SpecEValue. `TopNQueue` already keeps N; the
change is (i) ensure N>1 is honored end-to-end for chimeric, (ii) keep distinct
peptides rather than collapsing to the single best precursor. Emit one PIN row
per retained distinct peptide with `SpecId = <scan>_<rank>`.

**1e. Output.** `pin.rs` already writes one row per PSM; ensure `SpecId`
uniqueness per (scan, rank). No FDR-model change (Percolator is PSM-level).

**Gate:** `--chimeric` off → sorted-row PIN bit-identical to current on all 3
datasets. `--chimeric` on → measure PXD001819 + TMT PSM gain @1% FDR and the wall
cost; Astral expected ~flat (narrow windows). Ship if PXD/TMT gain without
Astral/​wall regression beyond the agreed envelope.

**Risk:** candidate explosion → wall. Mitigations: window-width cap; the
fragment-index enabler (see Cross-cutting). Measure first; if wall is
unacceptable, gate Phase 1 ship on the fragment-index landing.

### Phase 2 — MS1 targeted-XIC isotope refinement (additive)

**2a. Load + link MS1.** Widen the parser to `ms_level 1..=2` under `--chimeric`;
link each MS2 to its preceding MS1 (scan order / RT). Stream MS1 and retain only
a sliding window of recent MS1 scans to bound memory (MS1 scans are large; the
memory note `ms1-precursor-refinement-caveat` applies).

**2b. Targeted isotope XIC + KL divergence.** For each PSM, at its theoretical
neutral mass + charge, extract the observed isotope-envelope intensities from the
linked MS1 (and ±a few neighboring MS1s for the XIC apex), compute the
**KL-divergence** between observed and theoretical isotope distributions.

**2c. Additive PIN feature.** Emit `PrecursorIsotopeKL` (and optionally
`PrecursorXICApexCorr`) as **new** PIN columns — never modifying existing column
values (the audit-safe additive class). Percolator learns to downweight chimeric
false co-IDs. Optionally hard-filter PSMs above a KL threshold.

**Lean-Rust option (per user preference):** Phase 2 can instead be an **external
handoff** — msgf-rust emits the per-PSM theoretical (m/z, charge, RT) list; an
external tool computes the XIC/KL feature and appends the PIN column. The Phase-1
PIN is the interface; this keeps MS1 handling out of the Rust hot path.

**Gate:** additive column → PIN search output unchanged pre-Percolator except the
new column; keep if it gains @1% FDR, revert if it regresses (should be ≥ flat).

### Phase 3 — Greedy shared-fragment rescoring

For each spectrum's co-identified peptides, iterate **most-confident first**; mark
the fragment peaks it matched as "claimed"; exclude claimed peaks when scoring the
next peptide; re-score. Mirrors the Yu et al. 2023 DDA approach.

**Risk:** this **modifies the score distribution** (audit class that historically
regresses Percolator). Bench-gate hard; prefer emitting the shared-fragment-adjusted
score as an *additional* column first (additive) before letting it replace
RawScore. Revert if Percolator regresses.

### Cross-cutting — fragment-index candidate generator (speed enabler)

Widening to the full window multiplies candidates per spectrum; the protein-walk +
`bucket_index` scan will slow down. MSFragger pays for this with a fragment-ion
index. This repo already lists "fragment-index as candidate generator" as planned
speed work. **Decision needed:** (A) implement Phase 1 on the existing bucket scan,
accept the wall cost, add the index later; or (B) land the fragment-index first as
a prerequisite. Recommendation: measure Phase 1 wall on the existing scan first
(cheap), then decide — don't pre-build the index speculatively.

## Success criteria

- **Per phase:** `--chimeric` off is bit-identical; `--chimeric` on gains
  PXD001819 and/or TMT @1% FDR without regressing Astral or blowing the wall
  budget (target: within a multiple agreed with the user, since wide-window search
  is inherently more work).
- **Headline target:** a meaningful chimeric gain on PXD001819/TMT (DDA+ reports
  large gains on wide-window data; we set a concrete numeric gate after the Phase-1
  wall measurement).

## Risks & mitigations (summary)

| Risk | Mitigation |
|---|---|
| Candidate explosion → wall | window-width cap; fragment-index enabler; measure-first |
| Multi-PSM/scan FDR | Percolator is PSM-level; validate target/decoy balance on the wider PIN |
| Phase 3 modifies score dist (audit regress) | additive-first; bench-gate; revert path |
| MS1 memory | stream + sliding-window retention; or external-handoff for Phase 2 |
| False co-IDs | Phase-2 KL filter + Phase-3 shared-fragment removal + Percolator |

## Out of scope

- DL rescoring (MSBooster/Prosit-style) — complementary, separate effort.
- A full fragment-index rewrite — referenced as the speed enabler; its own roadmap item.
- DIA.

## References

- MSFragger-DDA+ (Nat Commun 2025): full-isolation-window search + MS1 XIC
  refinement + greedy shared-fragment rescoring. https://www.nature.com/articles/s41467-025-58728-z
- Yu et al. 2023 (MSFragger-DIA): the greedy shared-fragment rescoring method.
- `docs/parity-analysis/notes/2026-05-28-phase2-peak-rank-parity.md`: why parity
  fixes are exhausted (motivates this feature).
