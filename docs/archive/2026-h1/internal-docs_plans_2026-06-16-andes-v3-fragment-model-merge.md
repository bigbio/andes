# Andes v3 Fragment-Model Merge (MVP) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. Apply `internal-docs/experiment-protocol.md` (provenance banners, deploy-by-commit, milestone commits) for every train/benchmark step.

**Goal:** Replace the coarse `IntensityModel` lookup table inside `RawScore`'s `intensity_signal` cosine with ONE peptide-conditioned GBDT that predicts per-fragment expected relative intensity (regression), so a single trained fragment model powers RawScore (the v1 peptide-agnostic peak classifier is subsumed).

**Architecture:** New peptide-conditioned per-fragment feature vector → a regression GBDT (squared-error objective added to the existing trainer) trained on (fragment-features, observed-rel-intensity) pairs from gold PSMs → at inference, a per-fragment `predict_log_rel`-shaped call inside `intensity_signal` replaces the table lookup. The cosine, the additive `IntensitySignal` PIN column, and `RankScore` (rank-LLR) all stay byte-identical; only what feeds the cosine changes. Pure-Rust, rescoring (top-N), no NN.

**Tech Stack:** Rust (crates `scoring`, `model-train`, `andes`); existing zero-dep GBDT SoA walker (`gbdt_eval.rs`) + histogram regression trees (`gbdt/tree.rs`); Parquet model store; Percolator 3.7.1 for FDR; Codon (train) + VM (bench).

**Baseline to beat:** Astral `--score strong` RawScore = **38,909 PSMs@1%** (RE-VERIFY in Task 0 after the 64710d17 clobber fix). Gate: the v3 model must add on top of that at honest 1% entrapment-FDP on Astral, and not regress TMT/UPS.

> **⚠️ CORRECTION (2026-06-16, empirically confirmed):** `--score strong` is **HIGH-RES-ONLY**. On low-res it TANKS: TMT rank 11,428 → strong 9,624 (−16%); UPS rank 16,844 → strong 9,307 (−45%). The strong top-1 reranking is an Astral-specific win. **Do NOT A/B the frag model via `--score strong` on TMT/UPS** (Step 3 below says to — that is WRONG and would sandbag low-res). Instead deploy the frag features as the **additive Tier-2 LLR battery (commit 446c08a3) with `--score rank`** on low-res (ranking stays on the rank-LLR; the FragPred* columns ride along for Percolator). With that correct deployment the frag battery gave UPS +133 PSMs @ LOWER FDP (1.18→1.05%) and TMT flat. See [[andes-algorithm-audit-2026-06-16]].

**Out of MVP scope (later phases):** relational complement/neighbor features (box 2), explicit noise/chance-match LLR denominator (box 3), extra fragment vocabulary a-ions/losses/multi-charge (box 4), per-regime model proliferation beyond the 3 benchmark slugs.

---

## File structure

| File | Responsibility | Change |
|---|---|---|
| `crates/scoring/src/frag_features.rs` | peptide-CONDITIONED per-fragment feature vector (shared train+infer) | **CREATE** |
| `crates/scoring/src/lib.rs` | export `frag_features` | modify |
| `crates/model-train/src/gbdt/train.rs` | add `Objective::{Logistic,Regression}`; squared-error gradient; held-out R²/Pearson log | modify |
| `crates/model-train/src/gbdt/frag_dataset.rs` | build regression dataset (frag-features, observed rel-intensity) from gold PSMs | **CREATE** |
| `crates/model-train/src/gbdt/mod.rs` | export `frag_dataset` | modify |
| `crates/scoring/src/gbdt_eval.rs` | `predict_value` (regression: raw tree sum, no sigmoid/isotonic) | modify |
| `crates/scoring/src/scoring/strong_score.rs` | `intensity_signal`: use the frag-GBDT (if present) instead of `IntensityModel.predict_log_rel`; fallback-safe | modify |
| `crates/andes/src/bin/andes.rs` | `train-intensity-gbdt` subcommand (build frag dataset → regression train → store blob) + load the blob into `Param` | modify |
| model store (`scoring` store I/O) | carry the frag-GBDT blob (reuse the existing nullable `gbdt_model_bytes` column with an objective flag) | modify |

Design boundary: `frag_features.rs` is the single source of the feature vector (imported by BOTH the trainer's dataset builder and `intensity_signal`), guaranteeing train/infer parity (the same discipline as `PeakFeatureCtx::for_spectrum`).

---

## Task 0: Re-verify the RawScore strong baseline (hygiene prerequisite)

No code; establishes the honest baseline AFTER the 64710d17 clobber fix changed strong-mode PIN content.

- [ ] **Step 1:** Wait for the Java-Astral run (`b6qav62q3`) to free the VM. Then re-sync the VM repo to committed HEAD **whole-tree** (kill drift): from Mac, `git -C msgf-rust archive HEAD | ssh pride-linux-vm 'tar -x -C /srv/data/msgf-bench/repo/msgf-rust'` (or `git fetch+checkout` if the VM repo is a clone). Record `git rev-parse --short HEAD`.
- [ ] **Step 2:** Rebuild: `ssh pride-linux-vm 'cd <repo> && cargo +1.95.0 build --release -p andes --features thermo'`.
- [ ] **Step 3:** Print a provenance banner (`internal-docs/scripts/prov.sh`: `prov_bin` the andes binary, `prov_file model` the corpus43 store + the HCD intensity model, `prov_file data` the Astral mzML).
- [ ] **Step 4:** Re-run Astral `--score rank` and `--score strong` (store `store_corpus43_off.parquet`, `--intensity-model intensity_model.parquet`, entrapment DB) → Percolator → PSMs@1% + ENT-FDP. Record the NEW strong baseline (expected ≈ 38,909; the RankScore column is now the distinct rank-LLR, so it may shift slightly).
- [ ] **Step 5:** Log the result + provenance in `internal-docs/MILESTONES.md`. This number is the v3 gate target.

---

## Task 1: Per-fragment feature vector (`frag_features.rs`)

**Files:**
- Create: `crates/scoring/src/frag_features.rs`
- Modify: `crates/scoring/src/lib.rs` (add `pub mod frag_features;`)
- Test: inline `#[cfg(test)] mod tests` in `frag_features.rs`

Peptide-conditioned features for ONE annotated b/y ion. Residues are integer-encoded (0..25 from `aa.residue - b'A'`); the GBDT splits on them. All `f32`. NCE passed as a parsed `f32` (use `0.0` sentinel when unknown — matches the inference path that has no NCE yet).

- [ ] **Step 1: Write the failing test**
```rust
#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use crate::scoring::fragment_ions::IonKind;

    fn pep(seq: &str) -> Peptide {
        Peptide::new(seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect(), b'K', b'R')
    }

    #[test]
    fn frag_features_stable_shape_and_context() {
        let p = pep("PEPTIDE"); // n=7
        // b2: cleavage after residue 2 (E|P); flanks (E, P), proline on C-side.
        let f = extract_frag_features(&p, IonKind::B, 2, 1, 0.0);
        assert_eq!(f.len(), N_FRAG_FEATURES);
        // ion_type b = 0.0
        assert_eq!(f[FEAT_ION_TYPE], 0.0);
        // position fraction = 2/7
        assert!((f[FEAT_POS_FRAC] - 2.0 / 7.0).abs() < 1e-6);
        // proline flag set (P is the C-flank of b2)
        assert_eq!(f[FEAT_PROLINE_FLANK], 1.0);
        // length 7
        assert_eq!(f[FEAT_PEP_LEN], 7.0);
    }

    #[test]
    fn frag_features_reflect_modification() {
        // A Cam-C (+57) on the N-flank residue must change the mod-delta feature.
        let plain = extract_frag_features(&pep("PEACDEK"), IonKind::B, 3, 1, 0.0);
        // (build a modified peptide via the same helper used in gbdt::labels tests)
        let modded = extract_frag_features(&pep_modC("PEACDEK", 3, 57.02146), IonKind::B, 3, 1, 0.0);
        assert_ne!(plain[FEAT_NFLANK_MOD], modded[FEAT_NFLANK_MOD]);
    }
}
```
(Provide `pep_modC` test helper analogous to `gbdt/labels.rs::pep_with_mod`.)

- [ ] **Step 2: Run to verify it fails** — `cargo test -p scoring --lib frag_features` → FAIL (module/symbols not defined).

- [ ] **Step 3: Implement**
```rust
//! Peptide-CONDITIONED per-fragment features for the v3 intensity regressor.
//! The SINGLE source of the feature vector — imported by both the trainer's
//! frag-dataset builder and `intensity_signal` — so train/infer features match
//! by construction (mirrors `PeakFeatureCtx::for_spectrum`).

use model::peptide::Peptide;
use crate::scoring::fragment_ions::IonKind;

pub const FEAT_ION_TYPE: usize = 0;     // 0=b, 1=y
pub const FEAT_CHARGE: usize = 1;       // fragment charge
pub const FEAT_NFLANK: usize = 2;       // N-side residue index (0..25)
pub const FEAT_CFLANK: usize = 3;       // C-side residue index
pub const FEAT_PROLINE_FLANK: usize = 4;// 1.0 if P is N- or C-flank, else 0
pub const FEAT_POS_FRAC: usize = 5;     // cleavage position / length
pub const FEAT_PEP_LEN: usize = 6;      // peptide length
pub const FEAT_NFLANK_MOD: usize = 7;   // mod mass-delta on the N-flank residue
pub const FEAT_CFLANK_MOD: usize = 8;   // mod mass-delta on the C-flank residue
pub const FEAT_NCE: usize = 9;          // normalized collision energy (0.0 = unknown)
pub const N_FRAG_FEATURES: usize = 10;

/// Residues N- and C-side of the cleavage for `kind`/`position` (1-based),
/// mirroring `strong_score::flank_residues`. Returns indices into `residues`.
fn flank_indices(n: usize, kind: IonKind, position: u32) -> Option<(usize, usize)> {
    let i = position as usize;
    if i < 1 || i >= n { return None; }
    match kind {
        IonKind::B => Some((i - 1, i)),         // cleavage after residue i
        IonKind::Y => { let left = n - i; Some((left - 1, left)) }
    }
}

fn res_idx(b: u8) -> f32 { (b.wrapping_sub(b'A')) as f32 }
fn mod_delta(p: &Peptide, idx: usize) -> f32 {
    p.residues[idx].mod_.as_ref().map_or(0.0, |m| m.mass_delta) as f32
}

/// Feature vector for one annotated b/y ion. `position` is 1-based (b_i/y_i),
/// `charge` the fragment charge, `nce` the parsed NCE (0.0 when unknown).
pub fn extract_frag_features(p: &Peptide, kind: IonKind, position: u32, charge: u8, nce: f32) -> [f32; N_FRAG_FEATURES] {
    let n = p.residues.len();
    let mut f = [0.0f32; N_FRAG_FEATURES];
    let (ni, ci) = match flank_indices(n, kind, position) { Some(v) => v, None => return f };
    f[FEAT_ION_TYPE] = match kind { IonKind::B => 0.0, IonKind::Y => 1.0 };
    f[FEAT_CHARGE] = charge as f32;
    f[FEAT_NFLANK] = res_idx(p.residues[ni].residue);
    f[FEAT_CFLANK] = res_idx(p.residues[ci].residue);
    f[FEAT_PROLINE_FLANK] = if p.residues[ni].residue == b'P' || p.residues[ci].residue == b'P' { 1.0 } else { 0.0 };
    f[FEAT_POS_FRAC] = position as f32 / n as f32;
    f[FEAT_PEP_LEN] = n as f32;
    f[FEAT_NFLANK_MOD] = mod_delta(p, ni);
    f[FEAT_CFLANK_MOD] = mod_delta(p, ci);
    f[FEAT_NCE] = nce;
    f
}
```

- [ ] **Step 4: Run to verify pass** — `cargo test -p scoring --lib frag_features` → PASS.
- [ ] **Step 5: Commit** — `git commit -m "feat(scoring): peptide-conditioned per-fragment feature vector (v3 frag model)"`

---

## Task 2: Regression objective in the GBDT trainer

**Files:**
- Modify: `crates/model-train/src/gbdt/train.rs`
- Test: inline tests in `train.rs`

The current trainer is logistic (gradient `p − y`, output sigmoid+isotonic). Add `Objective::Regression` (squared-error: gradient `pred − y`, output = raw tree sum, NO sigmoid, NO isotonic) and a held-out **Pearson r** + **R²** log mirroring the existing AUC line.

- [ ] **Step 1: Write the failing test**
```rust
#[test]
fn regression_fits_linear_target_high_r2() {
    // y = 2*x0 - 1 with light noise; regression GBDT should reach high R² on val.
    let mut st = 0x1234u64;
    let n = 4000; let nf = 2;
    let mut x = Vec::new(); let mut y = Vec::new(); let mut g = Vec::new();
    for i in 0..n {
        let x0 = lcg(&mut st); let x1 = lcg(&mut st);
        x.push(x0); x.push(x1);
        y.push((2.0*x0 - 1.0) as f32); // continuous target (regression Dataset stores y as f32)
        g.push(i as u32); // each row its own group (val split is per-group)
    }
    let ds = RegressionDataset { x, y, groups: g, n_features: nf };
    let p = TrainParams { objective: Objective::Regression, ..Default::default() };
    let model = train_gbdt_regression(&ds, &p, 42);
    // predict ~ 2*x0-1 ⇒ high correlation; spot-check two points
    let lo = model.predict_value(&[0.0, 0.5]);
    let hi = model.predict_value(&[1.0, 0.5]);
    assert!(hi > lo + 0.5, "monotone in x0: lo={lo} hi={hi}");
}
```
(Decision: the regression target is continuous, so either (a) add a `RegressionDataset { y: Vec<f32> }` + `train_gbdt_regression`, or (b) generalize `Dataset` with `y: Vec<f32>` + an `objective` field. Prefer (b) if it stays clean; the test above assumes a regression entry point returning a model with `predict_value`. Keep the classification path byte-identical.)

- [ ] **Step 2: Run to verify it fails** — `cargo test -p model-train regression_fits_linear` → FAIL.

- [ ] **Step 3: Implement** — In `train.rs`:
  - Add `pub enum Objective { Logistic, Regression }` and `pub objective: Objective` to `TrainParams` (`Default` → `Logistic`, preserving existing behavior).
  - In the boosting loop, branch the gradient: logistic `sigmoid(raw) − y`; regression `raw − y` (raw = current prediction). The histogram regression tree (`tree.rs`) already fits residuals, so only the gradient + the final-prediction transform differ.
  - Final model: regression sets `apply_sigmoid = false`, leaves isotonic empty (`iso_x/iso_y` empty → `predict_value` returns the raw tree sum; see Task 4).
  - Validation metric: add `fn pearson_r2(pred: &[f32], y: &[f32]) -> (f64, f64)` and log `train-gbdt: regression val Pearson r = {r:.4}  R² = {r2:.4}  (n_val=…)` — the regression analogue of the AUC gate line.
  - Test helper for Pearson/R² (perfect fit → r=1, R²=1; constant pred → R²≤0).

- [ ] **Step 4: Run to verify pass** — `cargo test -p model-train regression` and the existing `auc_*` / determinism tests (classification path unchanged).
- [ ] **Step 5: Commit** — `git commit -m "feat(model-train): squared-error regression objective + Pearson/R² gate for the GBDT trainer"`

---

## Task 3: Frag-regression training-data builder (`frag_dataset.rs`)

**Files:**
- Create: `crates/model-train/src/gbdt/frag_dataset.rs`
- Modify: `crates/model-train/src/gbdt/mod.rs` (`pub mod frag_dataset;`)
- Test: inline tests

For each gold PSM: enumerate annotated b/y ions (`predict_by_ions(peptide, 1..=1)`), find the matched observed peak (`ScoredSpectrum::nearest_peak_full` within `param.mme`), and emit `(extract_frag_features(peptide, kind, position, charge, nce), target)` where `target = ln(observed_intensity / base_peak)` (the same log-rel-intensity space `IntensityModel` uses). Unmatched ions are skipped (no observed intensity). Group id = peptide sequence + charge (leakage-free split, reuse `dataset::group_id`).

- [ ] **Step 1: Write the failing test** — synthetic spectrum with peaks exactly at two b/y ion m/z of a known peptide; assert the builder emits 2 rows, each with `N_FRAG_FEATURES` features and a finite target ≈ `ln(obs/base)`.
- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement** — `pub fn build_frag_dataset(rows: &[PsmRow<'_>], scorer: &RankScorer) -> RegressionDataset`. Reuse `PsmRow` (already carries `&Peptide`). NCE = `0.0` (matches inference until NCE threading lands). Base peak = max active-peak intensity (as in `intensity_signal`).
- [ ] **Step 4: Run to verify pass.**
- [ ] **Step 5: Commit** — `git commit -m "feat(model-train): frag-intensity regression dataset builder from gold PSMs"`

---

## Task 4: Regression inference in `gbdt_eval.rs` + store/load the frag model

**Files:**
- Modify: `crates/scoring/src/gbdt_eval.rs`
- Modify: model-store I/O (the function that loads `gbdt_model_bytes` into `Param`)
- Test: inline tests

- [ ] **Step 1: Write the failing test** — build a `GbdtPeakModel` with `apply_sigmoid=false`, empty isotonic, and a one-leaf constant tree (value 0.7); assert `predict_value(x) == 0.7` (raw sum, no sigmoid/isotonic), and `predict_value` of a 2-tree model = sum of leaves.
- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement**
```rust
impl GbdtPeakModel {
    /// Raw regression prediction: sum of tree outputs, NO sigmoid, NO isotonic.
    /// Used by the v3 fragment-intensity regressor (predicts log-rel-intensity).
    pub fn predict_value(&self, x: &[f32]) -> f32 {
        self.trees.iter().map(|t| t.eval(x)).sum()
    }
}
```
  - Blob: reuse the existing `gbdt_model_bytes` column + the `AGBD` blob format; set the `apply_sigmoid=false`/empty-isotonic flags already round-tripped by `from_bytes`. Add an objective marker if needed so loaders know this is the frag regressor (a 1-byte flag in the blob header, default 0=classifier for back-compat). Load it into `Param` as `frag_intensity_model: Option<Arc<GbdtPeakModel>>` (separate from the existing `gbdt_peak_model`).
- [ ] **Step 4: Run to verify pass** — `cargo test -p scoring --lib gbdt_eval`.
- [ ] **Step 5: Commit** — `git commit -m "feat(scoring): regression predict_value + frag-intensity model store/load"`

---

## Task 5: Wire the frag regressor into `intensity_signal` (fallback-safe)

**Files:**
- Modify: `crates/scoring/src/scoring/strong_score.rs`
- Test: inline tests in `strong_score.rs`

In `intensity_signal`, when a frag-GBDT model is present, predict per-fragment log-rel-intensity from `extract_frag_features` instead of `IntensityModel.predict_log_rel`; everything else (cosine, observed matching, `IntensitySignal` PIN emission) unchanged. `RankScore` stays byte-identical (this is the strong-score numerator only).

- [ ] **Step 1: Write the failing test** — with a frag-GBDT that returns a higher predicted intensity for a specific (flank,pos) than the coarse table, assert the cosine for a peptide whose observed peaks match that prediction is ≥ the table-based cosine; assert `intensity_signal` returns 0.0 when neither model is present; assert byte-identical fallback to the table when the frag model is `None`.
- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement** — extend the `intensity_signal` signature to accept `frag_model: Option<&GbdtPeakModel>`; per predicted ion, `pred = if let Some(g) = frag_model { g.predict_value(&extract_frag_features(peptide, ion.kind, ion.position, precursor_charge, nce)).exp() } else { table.predict_log_rel(...).exp() }`. Thread `frag_model` from `match_engine::compute_psm_features` (it already loads optional models from `Param`). NCE = `"unknown"`/`0.0` as today.
- [ ] **Step 4: Run to verify pass** — `cargo test -p scoring --lib strong_score`.
- [ ] **Step 5: Commit** — `git commit -m "feat(scoring): intensity_signal uses the v3 frag-intensity regressor when present (fallback to table)"`

---

## Task 6: `andes train-intensity-gbdt` + load into search

**Files:**
- Modify: `crates/andes/src/bin/andes.rs`
- Test: `crates/andes/tests/` integration (small) or a CLI smoke

- [ ] **Step 1: Write the failing test** — a small integration test: build a `RegressionDataset` from 2 synthetic PSMs via `build_frag_dataset`, train, write the blob into a temp store, reload, assert a `frag_intensity_model` is present and `predict_value` is finite.
- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement** — add `train-intensity-gbdt` subcommand: read gold-PSM flats (reuse `read_msnet_parquet` → `PsmRow{ peptide:&psm.peptide, .. }`), `build_frag_dataset` → `train_gbdt_regression` → log Pearson/R² → write the frag-GBDT blob into `--out-store` (preserving existing models, like `train --gbdt`). Mirror `run_train_intensity`'s store-write discipline.
- [ ] **Step 4: Run to verify pass.**
- [ ] **Step 5: Commit** — `git commit -m "feat(andes): train-intensity-gbdt subcommand (v3 frag-intensity regressor)"`

---

## Task 7: Train per-slug on Codon + A/B gate on the VM

No new code; the experiment. Apply `internal-docs/experiment-protocol.md` throughout: provenance banner at the top of every script; deploy the andes binary by committed SHA; record a milestone.

- [ ] **Step 1:** On Codon (rebuild andes-bin from the pushed/rsynced commit), train the frag-GBDT per slug on the gold-PSM flats: `andes train-intensity-gbdt --in <slug flats> --out-store <store> ...`. Capture the held-out Pearson/R² per slug; sanity-gate (R² must be clearly > 0, e.g. ≥ 0.3, else the features/target are wrong — STOP and fix before benchmarking).
- [ ] **Step 2:** Name + register the models: `<slug>__fraggbdt__<commit8>__<date>.parquet` in `internal-docs/model-registry.tsv`.
- [ ] **Step 3:** scp the stores to the VM. A/B on Astral / TMT / UPS: `--score strong` WITH the frag-GBDT store vs the Task-0 RawScore baseline (table-based). Uniform Percolator; PSMs@1% at honest entrapment-FDP (Astral) + sane FDR (TMT/UPS); equal wall-time check (`strong ≤ ~110%` of rank wall).
- [ ] **Step 4:** Decision gate: the frag-GBDT RawScore must **exceed the Task-0 strong baseline on Astral at honest FDP** and **not regress TMT/UPS**. Record the table + provenance + verdict in `internal-docs/MILESTONES.md`.
- [ ] **Step 5:** If it passes → proceed to the box-2+ enrichments (relational/noise/fusion) as a follow-on plan, then the Phase-B release PR (gated). If it's flat/negative → the honest finding is that peptide-conditioned intensity regression doesn't beat the coarse table on this corpus; record it and stop (do not ship).

---

## Self-review

**Spec coverage:** Box 1 (sequence-conditioned fragment expectation) → Tasks 1–6. Regression-target GBDT risk → Task 2. Train/infer feature parity → Task 1 (shared module) + Task 3/5 (both import it). Gate (honest FDP A/B vs baseline) → Tasks 0 + 7. Prereqs (mod-aware labels, mme) → already DONE (5002df10), referenced not re-planned. Boxes 2–7 explicitly deferred (stated in header + Task 7 Step 5). ✓

**Placeholder scan:** Task 2 leaves one open design choice (RegressionDataset vs generalize Dataset) — flagged explicitly with a preference + a concrete test contract (`train_gbdt_regression` returns a model with `predict_value`), not a vague TODO. Task 4 blob objective-marker is specified (1-byte flag, default 0). No bare "TBD"/"handle edge cases". ✓

**Type consistency:** `extract_frag_features(&Peptide, IonKind, u32, u8, f32) -> [f32; N_FRAG_FEATURES]` used identically in Tasks 1/3/5. `predict_value(&[f32]) -> f32` defined in Task 4, used in Tasks 2-test/5. `PsmRow{ peptide: &Peptide }` consistent with the committed dataset.rs. `Objective`/`TrainParams.objective` introduced in Task 2, used in 3/6. ✓
