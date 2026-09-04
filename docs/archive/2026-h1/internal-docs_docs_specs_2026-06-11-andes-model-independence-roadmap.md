# Andes model independence roadmap — full MS-GF+ replacement for Apache-2.0 relicensing

Date: 2026-06-11
Status: approved program plan (decisions locked) → Phase 1 to writing-plans
Owner: ypriverol

## Why
Andes is a clean-room Rust reauthoring of MS-GF+ whose **bundled scoring models are still
derived from MS-GF+'s released parameter files**, so the project remains a derivative work under
the **UC Regents non-commercial license** (`LICENSE`). That license explicitly permits
copy/modify/distribute for *educational, research, and non-profit* use, but **reserves commercial
use**. To relicense Andes under **Apache-2.0** (commercial-friendly, patent grant), every shipped
model's MS-GF+-derived content must be replaced with independently-trained content.

## Definition of done ("in stone")
A `models.parquet` with **zero MS-GF+-derived bytes**:
- **Learned tables** (rank / intensity / noise / existence distributions) trained from our own
  public-data high-confidence PSMs.
- **Structural skeleton** (parent-mass partition boundaries, segment count, `frag_off` offset
  tables, `mme` fragment tolerance) **re-derived from our own data** — NOT inherited from the
  MS-GF+ `--seed-model` template. *(Decision locked: full structural independence.)*
- Every shipped model **matches or beats the MS-GF+ seed within tolerance** on a held-out
  benchmark: own ≥ 0.99 × seed target-PSMs at 1% FDR, with FDP ≤ seed. *(Decision locked.)*
- License flips to **Apache-2.0** only when all shipped models pass and no MS-GF+ bytes remain.
  *(Decision locked.)*

## Key lesson driving the ordering
This session's benchmark showed the current generation method overfits thin/narrow corpora
(retrained `hcd_qexactive_tryp` regressed −8.3% on Astral; trained on a single HeLa Q Exactive
file). Three agents converged: the estimator arithmetic is correct, but the method **fully
replaces the seed with corpus-only estimates, with no independent prior, no sparse-aware
smoothing, no min-spectra floor, and no acceptance gate**. **Therefore we FIX THE METHOD BEFORE
WE MASS-TRAIN** — otherwise we retrain ~40 models twice.

## Phases

### Phase 0 — Acceptance harness (the gate tool)
- Uniform held-out benchmark: own vs MS-GF+ seed at 1% FDR + FDP, per slug, via Percolator.
- Reuse the VM benchmark harness (3 datasets) + extend to a held-out set per slug.
- Deliverable: a one-command `slug → {own_psms, seed_psms, pass?}` gate.

### Phase 1 — Fix the generation method (code, `crates/model-train`, TDD) — GATES EVERYTHING
1. **Independent pooled prior** — build a broad "global own-model" from aggregated own
   high-confidence PSMs across many datasets/instruments; thin per-slug models shrink toward
   **it** (relicense-safe replacement for "regularize toward the MS-GF+ seed"). The blend
   machinery already exists (`estimate.rs:420-429`); add the pooled-own prior as the parent.
2. **Re-derive structure from own data** — adaptive parent-mass binning + min-spectra-per-
   partition floor, computed from the training corpus rather than read from the seed template
   (`estimate.rs:144-160`). This is *both* the structural-independence cut and the robustness fix.
3. **Sparse-aware (rank-window) smoothing** replacing flat Laplace (`estimate.rs:402`) — smooth
   the sparse high-rank tail, leave discriminative ranks 1–3 untouched.
4. **Fix noise-model over-softening** (the dilution issue — noise rank model flattened in the
   empirical diff; it sits in the LLR denominator).
5. **Acceptance gate wired into `train-from-msnet`** — evaluate on held-out, refuse/flag a model
   that loses to the seed. (`gate.rs:75-93` exists; wire into the training path, not just
   `--update`.)
Each item is TDD with a failing test first. Validation: retrain `hcd_qexactive_tryp` with the
fixed method and recover the Astral regression (own ≥ seed within tolerance) while holding the
UPS1/TMT gains.

### Phase 2 — Corpus program (breadth, not just volume)
- Per slug, curate a **diverse multi-instrument** public corpus (the single-file overfit is the
  cautionary tale). Priority order: high-res HCD/Astral, low-res CID, TMT → then ETD/UVPD/iTRAQ/
  phospho → then enzyme variants.
- **Enzyme models are gated on the enzyme feature**: Sub-project B (fix per-file harvest enzyme
  routing — the PXD000900 mis-routing bug, 20 flats quarantined) and Sub-project A (multi-enzyme
  search, spec `docs/specs/2026-06-11-andes-multi-enzyme-search-design.md`). Fold both in here.
- Harvest on EBI Codon (`/hps/nobackup/juan/pride/reanalysis/andes-training`, the existing
  pipeline) with the corrected routing.

### Phase 3 — Train → gate → iterate
- Retrain every slug with the Phase-1 method on its Phase-2 corpus. **The 10 models trained this
  session used the un-fixed method and are redone here.**
- Each passes the Phase-0 gate vs seed. Ledger: `slug → independent? → passes-gate?`. Iterate
  corpus/params on failures.

### Phase 4 — Assemble & full benchmark
- Assemble the final independent `models.parquet` (verified zero MS-GF+ bytes). Full 3-dataset,
  multi-engine benchmark vs seed + Java/Sage/MSFragger/Comet.

### Phase 5 — Relicense to Apache-2.0
- Replace `LICENSE`/`NOTICE`, update README badges, flip to Apache-2.0. Legal/IP sign-off.

## Branching
- Method fixes (Phase 1) land on `feat/enzyme-support` (already holds the enzyme spec +
  train-data-type-override) or a dedicated `feat/model-independence` branch — decide at plan time.
- Corpus/harvest work is on Codon (operational), tracked in the campaign memory.

## Open items folded in from this session
- Fresh-Java benchmark numbers (msgf2pin 3.7.1 broken → use `build_pins.py`) — finish for Phase 4.
- `feat/train-data-type-override` must be pushed + andes-bin rebuilt on Codon so protocol/enzyme
  labels come from code, not the manifest patch used for `cid_lowres_tryp_tmt`.
