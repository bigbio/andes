# SP-A: Hybrid Glyco Candidate-Gen + Backbone — Implementation Plan

**Goal:** Add a DB-constrained backbone branch (`backbone = precursor − known glycan`) that, unioned with the existing de-novo Y-ladder fallback, lifts hybrid searchable-backbone from 60% (de-novo-alone) to ≥90% on PXD025455.

**Architecture:** New clean-room glycan-composition enumerator + a DB-branch backbone enumerator, unioned with the existing `solve_backbone` de-novo solver, measured by the existing `glyco_probe` harness. All in `crates/andes-glyco/`. Reuse `oxonium.rs`, `backbone.rs`, `glycan_mass.rs`.

## Global Constraints
- Clean-room: glycan masses are combinatorial sums of monosaccharide monoisotopic masses (facts). No copied vendor list.
- Deterministic (total-order sorts); pure-Rust; no deps outside the crate beyond `input` (for mzML in the harness).
- Symmetric ±20 ppm searchable window: `|cand − truth| ≤ max(truth·20e-6, 0.01)`.
- Touch nothing outside `crates/andes-glyco/`.

---

### Task 1: Clean-room glycan composition list (`src/glycan_db.rs`)
**Files:** Create `crates/andes-glyco/src/glycan_db.rs`; register in `lib.rs`.
- `pub struct GlycanComp { pub hexnac:u8, pub hex:u8, pub fuc:u8, pub neuac:u8, pub neugc:u8, pub mass:f64 }`
- `pub fn n_glycan_list() -> Vec<GlycanComp>`: enumerate HexNAc 2..=8, Hex 3..=12, Fuc 0..=3, NeuAc 0..=5, NeuGc 0..=2 with N-glycan plausibility (`fuc ≤ hexnac`, `neuac+neugc ≤ max(0, hexnac−2)`), mass 500..=6000, sorted by mass asc. Masses from `glycan_mass.rs` constants.
- Tests: list non-empty (~1500–1700 comps); sorted; a known composition (HexNAc2Hex3 core = 892.317... wait compute) present; determinism (two calls equal).

### Task 2: DB-branch backbone enumeration (`src/hybrid.rs`)
**Files:** Create `crates/andes-glyco/src/hybrid.rs`.
- `pub struct BackboneHit { pub backbone_mass:f64, pub glycan:Option<GlycanComp>, pub source:Source }` where `pub enum Source { Db, DeNovo }`.
- `pub fn db_branch(precursor_neutral:f64, glycans:&[GlycanComp], min_backbone:f64) -> Vec<BackboneHit>`: for each glycan, `bb = precursor_neutral − glycan.mass`; keep if `bb ≥ min_backbone` (e.g. 500.0, smallest tryptic peptide); emit `BackboneHit{ source: Db }`. Sorted by backbone_mass.
- Tests: a synthetic precursor minus a known glycan yields that backbone with `Source::Db`; below-min filtered.

### Task 3: Hybrid union (`src/hybrid.rs`)
- `pub fn hybrid_candidates(peaks:&[(f64,f32)], precursor_neutral:f64, precursor_z:u8, glycans:&[GlycanComp], tol_ppm:f64, top_k:usize) -> Vec<BackboneHit>`: oxonium-gate (reuse `oxonium_gate`); run `db_branch`; ALSO run de-novo `solve_backbone` and append its candidates as `Source::DeNovo` (glycan None). Dedup backbones within 0.02 Da (prefer Db source). Return unioned set.
- Tests: union contains both sources; dedup keeps Db when a de-novo cluster coincides.

### Task 4: Hybrid gate in `glyco_probe` harness (`src/bin/glyco_probe.rs`)
- Extend the harness: for each truth scan, call `hybrid_candidates`; `searchable = any candidate within symmetric ±20ppm of truth backbone`. Report: hybrid searchable-backbone OVERALL %, the DB-hit vs de-novo-only split, and the fraction of searchable scans whose hit came from `Source::Db`. Keep the prior de-novo-only number for comparison.

### Task 5: Run the gate + record
- `cargo build -p andes-glyco`; run `glyco_probe` on the staged `HCC_pool_Late_Fc3_r1.mzML` + `truth.tsv`.
- Update `crates/andes-glyco/tests/data/PHASE1_RESULT.md` with hybrid number + source split. **Gate: hybrid searchable-backbone ≥ 90%.**
- Commit all SP-A code + result.
