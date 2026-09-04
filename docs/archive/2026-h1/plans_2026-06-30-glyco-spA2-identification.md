# SP-A2: Backbone Search + Glyco-PSM Identification — Implementation Plan

**Branch:** `glyco-phase1` (continue from SP-A, head 3aa93a6c)
**Predecessor gate:** SP-A hybrid searchable-backbone 90.4% (PASS).
**Goal:** Turn SP-A backbone candidates into scored glyco-PSMs, emit a Percolator PIN, run Percolator at 1% FDR, and compare glyco-PSM + glyco-peptide counts vs MSFragger-Glyco on PXD025455 `HCC_pool_Late_Fc3_r1`.

## Global Constraints
- New glyco logic lives ONLY in `crates/andes-glyco/` (candidate gen, sequon, glycan meta) and `crates/andes/src/glyco_search.rs` + `crates/andes/src/glyco_pin.rs` (driver, PIN). **Do NOT modify** any file in `crates/search/src/`, `crates/scoring/src/`, `crates/output/src/`, `crates/model/src/`.
- Fragment scoring: **bare backbone b/y only** — the Asn carries NO glycan mod in the peptide passed to `score_psm`. Correct for stepped-HCD (Q Exactive HF): the high-energy step strips glycan from b/y ions; glycan appears only in oxonium + Y-ladder ions (handled as additive PIN features). This matches MSFragger-Glyco's strategy and needs zero scoring changes.
- Glycan drives ONLY the precursor-mass lookup (`backbone = precursor − glycan`). Do NOT register the glycan as a fixed Asn mod (would force rebuilding `PreparedSearch` per composition — impractical).
- FDR: standard Percolator on the glyco PIN. (In-house GBDT 2D-FDR is SP-D.)
- Backbone candidate cap: **top-20 backbone hits per spectrum** by `core_y_hits DESC, backbone_mass ASC` before the FASTA loop (runtime guard against 200-glycan × 200-candidate explosion). Add a per-spectrum warn+truncate if >5,000 scoring events.
- Decoy: peptide-axis only (existing reversed-protein `PreparedSearch` decoys).
- Deterministic; pure-Rust; no new deps on `crates/andes-glyco/Cargo.toml`.

---

### Task 1: Sequon filter (`crates/andes-glyco/src/sequon.rs`)
**Create** `sequon.rs`; **modify** `lib.rs` (`pub mod sequon;`).
```rust
/// True iff residues (raw one-letter AA bytes) contain an N-X-S/T sequon, X≠P.
pub fn has_nxst_sequon(residues: &[u8]) -> bool;
```
Impl: `for i in 0..len.saturating_sub(2)`: `residues[i]==b'N' && residues[i+1]!=b'P' && (residues[i+2]==b'S'||residues[i+2]==b'T')`.
Tests: `SVNLTK`→true, `SVNPLTK`→false, `PEPTIDE`→false, `NST`→true, `NPT`→false, `NN`→false, ``→false.

### Task 2: Glycan-metadata carrier (`crates/andes-glyco/src/glyco_psm.rs`)
**Create**; **modify** `lib.rs`. Carries glycan provenance + glyco features only (NO `search` dep).
```rust
pub struct GlycoPsmKey {
    pub spectrum_idx: usize,
    pub glycan: Option<GlycanComp>,      // None for de-novo
    pub glycan_source: Source,           // Db | DeNovo
    pub oxonium_summed_frac: f32,
    pub n_core_oxonium_ions: u8,
    pub y_ladder_intensity_score: f32,
    pub core_y_hits: u8,
    pub glycan_mass: f64,                // glycan.map(|g|g.mass).unwrap_or(0.0)
    pub backbone_mass: f64,
}
```
Tests: construction; `glycan_mass`==0.0 for `None`.

### Task 3: Scoring driver (`crates/andes/src/glyco_search.rs`) — the integration seam
**Create**; **modify** `crates/andes/Cargo.toml` (add `andes-glyco = { path = "../andes-glyco" }`); **modify** `crates/andes/src/bin/andes.rs` (`mod glyco_search;`).
```rust
pub struct FullGlycoPsm { pub glycan_key: GlycoPsmKey, pub candidate_idx: u32, pub charge_used: u8,
    pub backbone_mass_error_ppm: f64, pub backbone_score: f32, pub rank_score: f32, pub edge_score: i32,
    pub is_decoy: bool, pub features: search::psm::PsmFeatures }
pub struct GlycoSpectrumResult { pub spectrum_idx: usize, pub hits: Vec<FullGlycoPsm> }
pub fn glyco_search_run(spectra:&[Spectrum], glycan_list:&[GlycanComp], prepared:&PreparedSearch,
    scorer:&RankScorer, fragment_tolerance_da:f64, tol_ppm:f64, backbone_top_k:usize) -> Vec<GlycoSpectrumResult>;
```
Per spectrum (rayon parallel): oxonium_gate once → for each charge `precursor_neutral=(mz−PROTON)*z−H2O`, `hybrid_candidates(...,50)` → union+dedup → cap top-20 by core_y_hits → for each backbone hit: `bucket_index.range(nominal±1)` → filter `matches_backbone_mass(cand, bb, 20ppm)` (|pep.mass()−bb| ≤ max(bb·20e-6,0.01)) → filter `has_nxst_sequon` → `score_psm` (bare peptide) + `psm_edge_score` + `compute_psm_features` → emit `FullGlycoPsm`. Dedup same `(candidate_idx, glycan)` keep max rank_score.
Note: `compute_psm_features` is `pub(crate)` in `match_engine.rs` — callable from the andes binary; if not, promote to `pub` (one line, allowed since it's a visibility-only change in search — EXCEPTION to the no-touch rule, get controller sign-off first).
Tests: synthetic spectrum+known glycan finds the N-X-S/T backbone; non-sequon excluded; wrong-mass excluded; top-20 cap honored; all-decoy → no targets.

### Task 4: Glyco PIN writer (`crates/andes/src/glyco_pin.rs`)
**Create**; **modify** `andes.rs` (`mod glyco_pin;`). Superset of standard PIN + appended columns: `OxoniumScore`(f64), `NCoreOxoniumIons`(int), `YLadderScore`(f64), `CoreYHits`(int), `GlycanMass`(f64), `IsGlycanDb`(int). `Peptide` col = `pre.SEQ[HexNAc{n}Hex{m}Fuc{f}NeuAc{a}NeuGc{g}].post` (omit bracket for de-novo). `CalcMass` = `pep.mass()+glycan_mass`.
Tests: header has glyco cols; CalcMass correct for known backbone+glycan; `[HexNAc..]` suffix for DB hit; `IsGlycanDb` 1/0; `Label` −1 decoy / +1 target.

### Task 5: CLI (`crates/andes/src/bin/andes.rs`)
Add `--glyco` (bool) + `--glyco-backbone-top-k` (default 20, hidden). When set: build `n_glycan_list()`, call `glyco_search_run`, `write_glyco_pin` to `<output>.glyco.pin`, skip standard PIN, return. Smoke test `crates/andes/tests/glyco_cli.rs`: `--glyco` on the staged mzML produces a `.glyco.pin` with `OxoniumScore` header + >0 PSM rows.

### Task 6: Decoy (no code) — reversed-protein decoys serve the peptide axis; document the known sequon-density TDC bias (deferred to SP-D 2D decoy).

### Task 7: Comparative experiment (Codon)
Provenance: record andes commit, MSFragger jar SHA, mzML SHA, FASTA SHA. Run andes `--glyco` on `HCC_pool_Late_Fc3_r1.mzML` + human-serum+cRAP FASTA → Percolator (record mode). Run MSFragger-Glyco (jar at `/hps/nobackup/juan/pride/reanalysis/andes-training/bin/`) with `glyco_search=true, N_glyco=true, precursor_mass_units=1 (ppm), ±20, trypsin, 2 missed, Cam-C fixed, Ox-M var` → Percolator (record mode). Count glyco-PSMs + distinct backbone peptides at q≤0.01 for both.
**Gotchas:** both Percolator runs MUST be Concatenated mode (grep logs); use ONE canonical target+decoy FASTA with matching prefix for both engines (verify decoy rows >0); both engines 2 missed cleavages; record both glycan-list sizes (coverage affects the gap).

### Task 8: Gate + record
**Gate:** andes glyco-PSMs@1%FDR ≥ **50%** of MSFragger-Glyco (baseline floor; full parity is SP-D after learned scoring + 2D-FDR). Below 50% ⇒ structural bug. Write `crates/andes-glyco/tests/data/SPA2_RESULT.md` (both counts, ratio, de-novo-sourced count, list sizes, Percolator mode, provenance). Commit code + result.

## Build sequence
Phase 1 (sequon+glyco_psm) → `cargo test -p andes-glyco`. Phase 2 (glyco_search) → `cargo test -p andes`. Phase 3 (glyco_pin) → test. Phase 4 (CLI + smoke on staged mzML) → build + local glyco PIN + andes glyco-PSM count. Phase 5 (Codon comparative) → gate.

## Top risks
1. **Combinatorial explosion** — top-20 cap + >5k-event truncation (highest runtime risk).
2. **MSFragger FASTA decoy-prefix mismatch** — one canonical target+decoy FASTA, verify decoy rows.
3. **`compute_psm_features` visibility** — `pub(crate)`, callable from andes binary; promote to `pub` only if needed (controller sign-off).
