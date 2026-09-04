# andes Phase A — Stronger Own Models (the public-release gate) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship andes's own, independently-trained models — a cleanly retrained rank core plus a peptide-agnostic GBDT per-peak signal/noise model that sharpens `RawScore`/`DeltaRawScore` — and prove they beat the open-source field at honest 1% FDR, which is the gate for the public release (Phase B, separate plan).

**Architecture:** The rank-LLR core stays (the natural low-res advantage). A LightGBM binary classifier scores each observed peak `P(signal | peptide-agnostic features)`, is isotonic-calibrated offline on Codon, transcoded to a flat struct-of-arrays tree table, and evaluated by a hand-rolled zero-native-dep Rust walker. Its calibrated probability becomes a **new additive LLR term** `log(s/(1−s))` per matched fragment, computed once per peak at spectrum-prep and cached by intensity rank, so the GBDT never runs in the inner candidate loop. When a model carries no GBDT blob the scoring path is byte-identical to today. Training/data happen on Codon; benchmarking on the VM.

**Tech Stack:** Rust (workspace crates `scoring`, `model-train`, `andes`; Arrow/Parquet 53 already present — **zero new native deps**), Python 3 + LightGBM + scikit-learn `IsotonicRegression` (Codon anaconda only, offline), SLURM (Codon), Percolator 3.7.1 + entrapment-FDP scripts (VM).

**Source of truth:** spec at `internal-docs/specs/2026-06-15-andes-stronger-models-public-release-design.md`. Repo: `/Users/yperez/work/msgfplus-workspace/msgf-rust`. This plan lives in `internal-docs/` (Claude planning docs never ship in the public repo). Training **code** (Python pipeline) DOES live in the repo — it is part of the reproducible tool, documented by `TRAIN.md`.

**Two-host split:** Codon = generate data + train models. VM = benchmark only. See `reference_andes_infra_layout` and the `codon-cluster` skill for access patterns (`become pst_prd`, base64-ship scripts, writes only on compute nodes, tag `--comment=andes:gbdt`).

**Critical parity invariant (read before any A2 task):** the Rust inference feature extractor (Task A2.1) and the Python training feature extractor (Task A2.6) MUST compute byte-for-byte the same feature vector in the same order. The ordered list in `crates/scoring/src/peak_features.rs::FEATURE_NAMES` is the single source of truth; Task A2.6 hardcodes the identical list and Task A2.10 is a cross-language parity gate. If you change a feature, you change it in both and re-run A2.10.

---

## File Structure

**Phase 0 (hygiene) — modify:**
- `crates/andes/src/bin/andes.rs` — `bundled_store_path()` (~line 3387) drops the `ionstat/` segment.
- `crates/model-train/src/store/read.rs:920`, `crates/model-train/src/store/migrate.rs` (param name + comments), and ~15 test/example files holding `../../resources/ionstat/models.parquet`.
- `README.md`, `DOCS.md`, `TRAIN.md`, `docs/benchmarks/2026-06-15-public-benchmark.md` — path + naming scrub.
- Filesystem: `resources/ionstat/models.parquet` → `resources/models.parquet`.

**A2 Rust (new files):**
- `crates/scoring/src/peak_features.rs` — peptide-agnostic per-peak feature contract + extractor.
- `crates/scoring/src/gbdt_eval.rs` — `GbdtPeakModel` (SoA trees + isotonic), `from_bytes`/`to_bytes`, `predict_logit`.

**A2 Rust (modify):**
- `crates/scoring/src/lib.rs` — declare the two new modules.
- `crates/scoring/src/param_model.rs` — add `Param.gbdt_peak_model: Option<GbdtPeakModel>`.
- `crates/scoring/src/scoring/scored_spectrum.rs` — `gbdt_logit_by_rank` cache + additive term in `directional_node_score_inner`.
- `crates/model-train/src/store/schema.rs` — nullable `gbdt_model_bytes` Binary column.
- `crates/model-train/src/store/write.rs` — write the blob.
- `crates/model-train/src/store/read.rs` — read the blob → `Param.gbdt_peak_model`.

**A2 Python / A1 / A3 training pipeline (new, in repo under `training/gbdt/`):**
- `training/gbdt/feature_spec.py` — the FEATURE_NAMES list (mirror of Rust) + shared constants.
- `training/gbdt/extract_features.py` — per-peak feature + label extraction from gold PSMs.
- `training/gbdt/train_gbdt.py` — LightGBM + isotonic, group-split.
- `training/gbdt/transcode.py` — LightGBM dump + isotonic → SoA blob; emits Rust parity fixtures.
- `training/gbdt/README.md` — the reproducible flow (referenced from `TRAIN.md`).

**A4 (VM, reuse existing):**
- `scripts/bench_astral_competitors.sh`, `scripts/bench_tmt_ups_competitors.sh`, `scripts/astral_entrapment_experiment.sh`, `scripts/entrap_fdp.py` (already on the VM; add a LysC variant).

---

## Phase 0 — Hygiene

### Task 0.1: Archive the MS-GF+-derived seed store outside the repo

**Files:**
- Read: `resources/ionstat/models.parquet` (the 10.8 MB seed store)
- Create: `/Users/yperez/work/msgfplus-workspace/internal-docs/model-archive/msgf-derived-seed.parquet`

- [ ] **Step 1: Create the archive directory and copy the seed store**

```bash
mkdir -p /Users/yperez/work/msgfplus-workspace/internal-docs/model-archive
cp /Users/yperez/work/msgfplus-workspace/msgf-rust/resources/ionstat/models.parquet \
   /Users/yperez/work/msgfplus-workspace/internal-docs/model-archive/msgf-derived-seed.parquet
```

- [ ] **Step 2: Verify the copy is byte-identical and record provenance**

```bash
cmp /Users/yperez/work/msgfplus-workspace/msgf-rust/resources/ionstat/models.parquet \
    /Users/yperez/work/msgfplus-workspace/internal-docs/model-archive/msgf-derived-seed.parquet \
  && echo "ARCHIVE OK"
```
Expected: `ARCHIVE OK`. Then create `internal-docs/model-archive/README.md`:

```markdown
# Model archive (internal, never shipped)

- `msgf-derived-seed.parquet` — the original MS-GF+-derived bundled store
  (copied from `msgf-rust/resources/ionstat/models.parquet`, 2026-06-15).
  Retained ONLY as the training **seed/prior** for `andes train-from-msnet`
  (`--seed-model <path>`). It is MS-GF+-derived IP and MUST NOT be shipped in
  the public release. The public `resources/models.parquet` is rebuilt from
  own-trained slugs in Phase A (A3) and replaces it entirely.
```

- [ ] **Step 3: Commit (the archive is outside the repo; commit only the pointer note in internal-docs which is its own area — no repo change here)**

This task makes no change to the `msgf-rust` repo, so there is nothing to `git add` in the repo. The archive + README live under `internal-docs/`. Proceed.

---

### Task 0.2: Rename `resources/ionstat/` → `resources/models.parquet`

**Files:**
- Modify: `crates/andes/src/bin/andes.rs` (`bundled_store_path()`, ~3387–3403)
- Modify: `crates/model-train/src/store/read.rs:920` (fallback default path)
- Modify (test/example paths, ~15 files): `crates/scoring/tests/param_loads_all_bundled.rs`, `crates/scoring/examples/dump_prefix_cache.rs`, `crates/model-train/tests/roundtrip.rs:173`, `crates/model-train/examples/{train_dump.rs:76,local_yield.rs:43,gen_bundled_store.rs}`, `crates/model-train/tests/yield_nonregression.rs`, `crates/model-train/tests/migration_parity.rs`, `crates/andes/tests/store_selection_equivalence.rs`
- Filesystem: move the store file

- [ ] **Step 1: Move the store file and remove the empty dir**

```bash
cd /Users/yperez/work/msgfplus-workspace/msgf-rust
git mv resources/ionstat/models.parquet resources/models.parquet 2>/dev/null \
  || mv resources/ionstat/models.parquet resources/models.parquet
rmdir resources/ionstat 2>/dev/null || true
ls resources/
```
Expected: `models.parquet  unimod.obo` (no `ionstat/`).

- [ ] **Step 2: Update `bundled_store_path()` in andes.rs**

Replace the two `resources/ionstat/models.parquet` literals (the `next_to_binary` join and the `CARGO_MANIFEST_DIR` fallback) with `resources/models.parquet`:

```rust
fn bundled_store_path() -> PathBuf {
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let next_to_binary = dir.join("resources/models.parquet");
            if next_to_binary.exists() {
                return next_to_binary;
            }
        }
    }
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models.parquet")
}
```

- [ ] **Step 3: Update the read.rs fallback default path and every test/example literal**

Run a workspace-wide replace of the path segment (it is unambiguous):

```bash
cd /Users/yperez/work/msgfplus-workspace/msgf-rust
grep -rl 'resources/ionstat/models.parquet' crates/ \
  | xargs sed -i '' 's#resources/ionstat/models.parquet#resources/models.parquet#g'
```
(On Linux/Codon drop the empty `''` after `-i`.)

- [ ] **Step 4: Verify no code path still references the old store location**

```bash
grep -rn 'resources/ionstat' crates/ ; echo "exit=$?"
```
Expected: only matches (if any) are in comments mentioning history; `exit=1` if none. Remaining matches must be comment-only and are handled in Task 0.3.

- [ ] **Step 5: Build + run the store-loading tests**

```bash
cargo test -p scoring --test param_loads_all_bundled
cargo test -p model-train --test roundtrip
```
Expected: PASS — the bundled store loads from its new path.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "chore: move bundled store to resources/models.parquet (drop ionstat/ segment)"
```

---

### Task 0.3: Scrub `ionstat` / MS-GF+ naming from code and tool docs

**Files:**
- Modify: `crates/model-train/src/store/migrate.rs` (rename the `ionstat: &Path` parameter), and any comment referencing `ionstat`/`*.param` migration
- Modify: `README.md`, `DOCS.md`, `TRAIN.md`, `docs/benchmarks/2026-06-15-public-benchmark.md` (`resources/ionstat/` → `resources/`; remove "ionstat" wording)

- [ ] **Step 1: Rename the `migrate_dir` parameter for clarity**

In `crates/model-train/src/store/migrate.rs`, rename the parameter `ionstat: &Path` to `param_dir: &Path` (and its uses inside the function body). Update the doc comment to say "Migrate a directory of legacy `.param` files into a combined Parquet store" with no "ionstat" reference.

- [ ] **Step 2: Scrub the four doc files**

In `README.md`, `DOCS.md`, `TRAIN.md`, `docs/benchmarks/2026-06-15-public-benchmark.md` replace `resources/ionstat/models.parquet` → `resources/models.parquet` and `resources/ionstat/` → `resources/`. In `docs/benchmarks/2026-06-15-public-benchmark.md` the "Models:" line (currently `resources/ionstat/models.parquet, which is MS-GF+-derived and in transition`) is rewritten in Phase B; for now only fix the path.

- [ ] **Step 3: Verify the scrub**

```bash
cd /Users/yperez/work/msgfplus-workspace/msgf-rust
grep -rn 'ionstat' crates/ README.md DOCS.md TRAIN.md docs/ ; echo "exit=$?"
```
Expected: `exit=1` (no matches). If a comment still mentions it for historical reasons, rewrite it to not use the word.

- [ ] **Step 4: Full build + clippy + test to confirm nothing broke**

```bash
cargo build --workspace
cargo clippy --workspace -- -D warnings
cargo test -p scoring -p model-train
```
Expected: build clean, clippy clean, tests pass (the pre-existing `integer_mass_scaler_matches_residue_table_mean` failure is the only known exception — see memory; do not "fix" it here).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "chore: scrub ionstat naming from code + tool docs"
```

---

## Phase A2 (Rust) — the GBDT engine

### Task A2.1: Peptide-agnostic per-peak feature contract + extractor

**Files:**
- Create: `crates/scoring/src/peak_features.rs`
- Modify: `crates/scoring/src/lib.rs` (add `pub mod peak_features;`)
- Test: inline `#[cfg(test)]` in `peak_features.rs`

This is the **single source of truth** for the feature vector. The extractor takes the active peak list (post precursor-filter / post-deconv, ascending m/z), the precursor m/z, charge, and the ranks vector, and returns one `[f32; N_FEATURES]` per peak, indexed parallel to the active peak list. Every feature is peptide-agnostic (no candidate peptide is consulted).

- [ ] **Step 1: Write the failing test (feature count + a hand-checked vector)**

Add to a new file `crates/scoring/src/peak_features.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_names_count_matches() {
        assert_eq!(FEATURE_NAMES.len(), N_FEATURES);
    }

    #[test]
    fn extracts_expected_features_on_tiny_scan() {
        // Three peaks; precursor m/z 500.0, charge 2 → neutral 997.985.
        // peaks ascending by m/z: (100.0, 10.0), (200.05, 100.0), (300.0, 50.0)
        let peaks = vec![(100.0_f64, 10.0_f32), (200.05, 100.0), (300.0, 50.0)];
        // ranks by intensity desc: idx1 (100)→1, idx2 (50)→2, idx0 (10)→3
        let ranks = vec![3u32, 1, 2];
        let ctx = PeakFeatureCtx {
            precursor_mz: 500.0,
            charge: 2,
            parent_neutral_mass: (500.0 - PROTON) * 2.0,
            total_intensity: 160.0,
            base_peak_intensity: 100.0,
            window_da: 50.0,
            match_tol_da: 0.5,
        };
        let f = extract_peak_features(&peaks, &ranks, &ctx);
        assert_eq!(f.len(), 3);
        // peak 1 (mz 200.05, intensity 100, global rank 1):
        let p1 = &f[1];
        // intensity_over_basepeak = 100/100 = 1.0
        assert!((p1[idx("intensity_over_basepeak")] - 1.0).abs() < 1e-6);
        // intensity_over_tic = 100/160
        assert!((p1[idx("intensity_over_tic")] - (100.0 / 160.0)).abs() < 1e-6);
        // global_rank_frac = (1-1)/3 = 0.0
        assert!((p1[idx("global_rank_frac")] - 0.0).abs() < 1e-6);
        // mass_defect = 200.05 - 200.0 = 0.05
        assert!((p1[idx("mass_defect")] - 0.05).abs() < 1e-4);
        // mz_frac_of_precursor = 200.05/500.0
        assert!((p1[idx("mz_frac_of_precursor")] - (200.05 / 500.0)).abs() < 1e-6);
        // is_top1_in_window: within ±50 Da of 200.05 the peaks are {200.05}
        // (100.0 is >50 away, 300.0 is >50 away) → it is the top → 1.0
        assert!((p1[idx("is_top1_in_window")] - 1.0).abs() < 1e-6);
    }

    // Test-only helper to look up a feature column by name.
    fn idx(name: &str) -> usize {
        FEATURE_NAMES.iter().position(|&n| n == name).expect("feature exists")
    }
}
```

- [ ] **Step 2: Run it to confirm it fails to compile (symbols undefined)**

Run: `cargo test -p scoring --lib peak_features`
Expected: FAIL — `cannot find function extract_peak_features`, `FEATURE_NAMES`, etc.

- [ ] **Step 3: Implement the contract + extractor**

Write the rest of `crates/scoring/src/peak_features.rs` (above the test module). `PROTON` is re-exported from the spectrum constants; import the same constant the rest of `scoring` uses (`use crate::scoring::scored_spectrum::PROTON;` if public, else define `const PROTON: f64 = 1.007276466;` locally with a comment that it must equal the engine's proton mass).

```rust
//! Peptide-AGNOSTIC per-peak signal/noise features (the GBDT input contract).
//!
//! THIS LIST IS THE SINGLE SOURCE OF TRUTH. The Python training extractor
//! (`training/gbdt/feature_spec.py`) MUST mirror `FEATURE_NAMES` in the same
//! order; the cross-language parity test (Task A2.10) enforces it. Every
//! feature is computable once per spectrum (no candidate peptide is consulted),
//! so the GBDT is evaluated once per peak at spectrum-prep, never in the inner
//! candidate loop.

/// Ordered feature names. Index in this array == feature index used by the
/// GBDT tree splits. DO NOT reorder without retraining + retranscoding.
pub const FEATURE_NAMES: [&str; 18] = [
    "log_intensity",          // 0  ln(intensity)
    "intensity_over_basepeak",// 1  intensity / max intensity in scan
    "intensity_over_tic",     // 2  intensity / summed kept intensity
    "global_rank_frac",       // 3  (rank-1) / kept_count
    "local_rank_frac",        // 4  (rank within ±window) / count in window
    "is_top1_in_window",      // 5  1.0 if most intense in ±window else 0.0
    "is_top3_in_window",      // 6  1.0 if among top-3 in ±window else 0.0
    "mz",                     // 7  observed m/z
    "mz_frac_of_precursor",   // 8  mz / precursor_mz
    "local_peak_density",     // 9  peaks per Da in ±window
    "spacing_left",           // 10 mz - previous peak mz (SENTINEL if none)
    "spacing_right",          // 11 next peak mz - mz (SENTINEL if none)
    "mass_defect",            // 12 mz - floor(mz)
    "has_isotope_plus1",      // 13 peak at mz + 1.00235/charge within tol
    "has_isotope_minus1",     // 14 peak at mz - 1.00235/charge within tol
    "has_complement",         // 15 peak at (M + 2*PROTON - mz) within tol
    "has_h2o_loss_partner",   // 16 peak at mz - 18.010565/charge within tol
    "has_nh3_loss_partner",   // 17 peak at mz - 17.026549/charge within tol
];

pub const N_FEATURES: usize = FEATURE_NAMES.len();

/// Sentinel for spacing when there is no neighbor on that side (Da). Large so
/// trees can isolate edge peaks.
pub const SPACING_SENTINEL: f32 = 1000.0;

/// Half-width (Da) of the local window for local-rank / density / top-K
/// features. MUST equal `training/gbdt/feature_spec.py::FEATURE_WINDOW_DA` and
/// the value passed by `build_dataset.py` at training time — diverging would
/// silently break feature parity (A2.10 guards it).
pub const FEATURE_WINDOW_DA: f64 = 50.0;

const ISOTOPE_SPACING: f64 = 1.00235; // average mass diff between isotopes (Da)
const H2O: f64 = 18.010565;
const NH3: f64 = 17.026549;
const PROTON: f64 = 1.007276466; // must equal the engine proton mass

/// Per-spectrum constants needed to compute the features.
pub struct PeakFeatureCtx {
    pub precursor_mz: f64,
    pub charge: u8,
    /// (precursor_mz - PROTON) * charge — the observed neutral parent mass.
    pub parent_neutral_mass: f64,
    /// Summed intensity of the active (kept) peaks.
    pub total_intensity: f64,
    /// Maximum intensity among the active peaks.
    pub base_peak_intensity: f32,
    /// Half-width (Da) of the local window for local-rank/density/top-K feats.
    pub window_da: f64,
    /// Tolerance (Da) for partner-peak existence flags (isotope/complement/loss).
    pub match_tol_da: f64,
}

/// Compute one feature vector per peak. `peaks` MUST be ascending by m/z and
/// aligned 1:1 with `ranks` (rank 1 = most intense; `u32::MAX` = filtered out).
/// Filtered-out peaks still get a row (so indices stay aligned) but their
/// rank-based features use the kept-count denominator.
pub fn extract_peak_features(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    ctx: &PeakFeatureCtx,
) -> Vec<[f32; N_FEATURES]> {
    let n = peaks.len();
    let kept_count = ranks.iter().filter(|&&r| r != u32::MAX).count().max(1);
    let tol = ctx.match_tol_da;
    let mut out = Vec::with_capacity(n);

    for i in 0..n {
        let (mz, intensity) = peaks[i];
        let rank = ranks[i];
        let mut f = [0.0_f32; N_FEATURES];

        f[0] = (intensity.max(1e-6)).ln();
        f[1] = if ctx.base_peak_intensity > 0.0 { intensity / ctx.base_peak_intensity } else { 0.0 };
        f[2] = if ctx.total_intensity > 0.0 { (intensity as f64 / ctx.total_intensity) as f32 } else { 0.0 };
        f[3] = if rank == u32::MAX { 1.0 } else { (rank.saturating_sub(1)) as f32 / kept_count as f32 };

        // Local window: peaks with m/z in [mz - window, mz + window].
        let lo = peaks.partition_point(|&(m, _)| m < mz - ctx.window_da);
        let hi = peaks.partition_point(|&(m, _)| m <= mz + ctx.window_da);
        let win = &peaks[lo..hi];
        let win_count = win.len().max(1);
        // local rank = #peaks in window strictly more intense than this one, +1.
        let more_intense = win.iter().filter(|&&(_, pint)| pint > intensity).count();
        f[4] = more_intense as f32 / win_count as f32;
        f[5] = if more_intense == 0 { 1.0 } else { 0.0 };
        f[6] = if more_intense < 3 { 1.0 } else { 0.0 };

        f[7] = mz as f32;
        f[8] = if ctx.precursor_mz > 0.0 { (mz / ctx.precursor_mz) as f32 } else { 0.0 };
        f[9] = (win_count as f64 / (2.0 * ctx.window_da)) as f32;

        f[10] = if i > 0 { (mz - peaks[i - 1].0) as f32 } else { SPACING_SENTINEL };
        f[11] = if i + 1 < n { (peaks[i + 1].0 - mz) as f32 } else { SPACING_SENTINEL };
        f[12] = (mz - mz.floor()) as f32;

        let z = ctx.charge.max(1) as f64;
        f[13] = has_peak(peaks, mz + ISOTOPE_SPACING / z, tol);
        f[14] = has_peak(peaks, mz - ISOTOPE_SPACING / z, tol);
        // Complement of a singly-charged fragment: b_i + y_(n-i) = M + 2*PROTON.
        f[15] = has_peak(peaks, ctx.parent_neutral_mass + 2.0 * PROTON - mz, tol);
        f[16] = has_peak(peaks, mz - H2O / z, tol);
        f[17] = has_peak(peaks, mz - NH3 / z, tol);

        out.push(f);
    }
    out
}

/// 1.0 if any peak lies within `tol` Da of `target`, else 0.0. `peaks` ascending.
fn has_peak(peaks: &[(f64, f32)], target: f64, tol: f64) -> f32 {
    if target <= 0.0 {
        return 0.0;
    }
    let lo = peaks.partition_point(|&(m, _)| m < target - tol);
    if lo < peaks.len() && peaks[lo].0 <= target + tol {
        1.0
    } else {
        0.0
    }
}
```

Add `pub mod peak_features;` to `crates/scoring/src/lib.rs` (next to the other `pub mod` lines).

- [ ] **Step 4: Run the tests**

Run: `cargo test -p scoring --lib peak_features`
Expected: PASS (both tests).

- [ ] **Step 5: Commit**

```bash
git add crates/scoring/src/peak_features.rs crates/scoring/src/lib.rs
git commit -m "feat(scoring): peptide-agnostic per-peak feature contract for GBDT"
```

---

### Task A2.2: SoA GBDT tree walker (`GbdtPeakModel`)

**Files:**
- Create: `crates/scoring/src/gbdt_eval.rs`
- Modify: `crates/scoring/src/lib.rs` (add `pub mod gbdt_eval;`)
- Test: inline `#[cfg(test)]` in `gbdt_eval.rs`

Defines the binary blob format (v1), the SoA in-memory model, `from_bytes`/`to_bytes` round-trip, and `predict_logit`. The Python transcoder (A2.9) writes the same format; the parity fixtures it emits feed Step 1's exact-match test in A2.10. Zero new deps — only `std` + `byteorder` (already a scoring dep).

**Blob format v1 (little-endian):**
```
magic        [u8;4] = b"AGBD"
version      u32     = 1
n_features   u32
flags        u32     bit0 = apply sigmoid to raw sum (LightGBM binary objective)
n_trees      u32
repeat n_trees times:
  n_nodes    u32
  feature    i32 * n_nodes   (-1 => leaf node)
  threshold  f32 * n_nodes   (numeric split: go left if x <= threshold)
  left       i32 * n_nodes   (child node index; leaf => -1)
  right      i32 * n_nodes   (leaf => -1)
  value      f32 * n_nodes   (leaf output; internal nodes => 0.0)
  default_left u8 * n_nodes   (1 => NaN feature goes left)
n_iso        u32             (isotonic calibration breakpoints)
iso_x        f32 * n_iso     (ascending uncalibrated prob breakpoints)
iso_y        f32 * n_iso     (monotone nondecreasing calibrated probs)
```

- [ ] **Step 1: Write the failing tests (round-trip + known prediction)**

In `crates/scoring/src/gbdt_eval.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    /// One tree, one split on feature 0 at threshold 0.5:
    ///   x0 <= 0.5 -> leaf value -1.0 ; else leaf value +2.0
    /// flags: sigmoid ON. isotonic: identity over [0,1] (two breakpoints).
    fn toy_model() -> GbdtPeakModel {
        GbdtPeakModel {
            n_features: 1,
            apply_sigmoid: true,
            trees: vec![Tree {
                feature: vec![0, -1, -1],
                threshold: vec![0.5, 0.0, 0.0],
                left: vec![1, -1, -1],
                right: vec![2, -1, -1],
                value: vec![0.0, -1.0, 2.0],
                default_left: vec![1, 1, 1],
            }],
            iso_x: vec![0.0, 1.0],
            iso_y: vec![0.0, 1.0],
        }
    }

    #[test]
    fn roundtrip_bytes() {
        let m = toy_model();
        let bytes = m.to_bytes();
        let back = GbdtPeakModel::from_bytes(&bytes).expect("decode");
        assert_eq!(back.n_features, m.n_features);
        assert_eq!(back.trees.len(), 1);
        assert_eq!(back.trees[0].value, m.trees[0].value);
        assert_eq!(back.iso_x, m.iso_x);
    }

    #[test]
    fn predict_matches_manual() {
        let m = toy_model();
        // x0 = 0.0 -> leaf -1.0 -> raw=-1.0 -> sigmoid(-1)=0.26894 -> iso identity
        let s_lo = m.predict_proba(&[0.0]);
        assert!((s_lo - 0.2689414).abs() < 1e-5, "got {s_lo}");
        // x0 = 1.0 -> leaf +2.0 -> sigmoid(2)=0.880797
        let s_hi = m.predict_proba(&[1.0]);
        assert!((s_hi - 0.8807971).abs() < 1e-5, "got {s_hi}");
        // logit(s) recovers the raw sum for the identity isotonic map.
        let lg = m.predict_logit(&[1.0]);
        assert!((lg - 2.0).abs() < 1e-4, "logit got {lg}");
    }

    #[test]
    fn empty_iso_is_identity() {
        // No isotonic breakpoints -> calibrated prob == sigmoid(raw).
        let mut m = toy_model();
        m.iso_x.clear();
        m.iso_y.clear();
        let s = m.predict_proba(&[1.0]);
        assert!((s - 0.8807971).abs() < 1e-5);
    }
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p scoring --lib gbdt_eval`
Expected: FAIL — types/methods undefined.

- [ ] **Step 3: Implement the model + walker**

Above the test module in `crates/scoring/src/gbdt_eval.rs`:

```rust
//! Hand-rolled struct-of-arrays GBDT evaluator (zero native deps).
//!
//! Decodes the `AGBD` v1 blob (produced offline by `training/gbdt/transcode.py`
//! from a LightGBM binary classifier + scikit-learn IsotonicRegression) and
//! evaluates it on a peptide-AGNOSTIC per-peak feature vector
//! (`crate::peak_features`). Output is `log(s/(1-s))`, the additive LLR term the
//! rank scorer folds in (`scored_spectrum`).

use std::io::{Cursor, Read};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use thiserror::Error;

const MAGIC: &[u8; 4] = b"AGBD";
const FORMAT_VERSION: u32 = 1;
const FLAG_SIGMOID: u32 = 1;
/// Clamp for the calibrated probability before taking the logit, so the LLR
/// term stays finite (a degenerate 0 or 1 would give ±inf).
const PROB_EPS: f32 = 1e-6;

#[derive(Debug, Error)]
pub enum GbdtError {
    #[error("bad magic (not an AGBD blob)")]
    BadMagic,
    #[error("unsupported AGBD version {0}")]
    BadVersion(u32),
    #[error("truncated AGBD blob: {0}")]
    Io(#[from] std::io::Error),
}

/// One regression tree in struct-of-arrays layout. All vecs have length
/// `n_nodes`; node 0 is the root.
#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    pub feature: Vec<i32>,     // -1 => leaf
    pub threshold: Vec<f32>,
    pub left: Vec<i32>,
    pub right: Vec<i32>,
    pub value: Vec<f32>,       // leaf output
    pub default_left: Vec<u8>, // 1 => NaN feature descends left
}

impl Tree {
    /// Sum of the leaf value reached for `x`.
    fn eval(&self, x: &[f32]) -> f32 {
        let mut node = 0usize;
        loop {
            let feat = self.feature[node];
            if feat < 0 {
                return self.value[node];
            }
            let v = x.get(feat as usize).copied().unwrap_or(f32::NAN);
            let go_left = if v.is_nan() {
                self.default_left[node] == 1
            } else {
                v <= self.threshold[node]
            };
            node = if go_left { self.left[node] } else { self.right[node] } as usize;
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct GbdtPeakModel {
    pub n_features: u32,
    pub apply_sigmoid: bool,
    pub trees: Vec<Tree>,
    pub iso_x: Vec<f32>,
    pub iso_y: Vec<f32>,
}

impl GbdtPeakModel {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, GbdtError> {
        let mut c = Cursor::new(bytes);
        let mut magic = [0u8; 4];
        c.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(GbdtError::BadMagic);
        }
        let version = c.read_u32::<LittleEndian>()?;
        if version != FORMAT_VERSION {
            return Err(GbdtError::BadVersion(version));
        }
        let n_features = c.read_u32::<LittleEndian>()?;
        let flags = c.read_u32::<LittleEndian>()?;
        let apply_sigmoid = flags & FLAG_SIGMOID != 0;
        let n_trees = c.read_u32::<LittleEndian>()?;

        let mut trees = Vec::with_capacity(n_trees as usize);
        for _ in 0..n_trees {
            let n = c.read_u32::<LittleEndian>()? as usize;
            let feature = read_i32_vec(&mut c, n)?;
            let threshold = read_f32_vec(&mut c, n)?;
            let left = read_i32_vec(&mut c, n)?;
            let right = read_i32_vec(&mut c, n)?;
            let value = read_f32_vec(&mut c, n)?;
            let mut default_left = vec![0u8; n];
            c.read_exact(&mut default_left)?;
            trees.push(Tree { feature, threshold, left, right, value, default_left });
        }
        let n_iso = c.read_u32::<LittleEndian>()? as usize;
        let iso_x = read_f32_vec(&mut c, n_iso)?;
        let iso_y = read_f32_vec(&mut c, n_iso)?;
        Ok(Self { n_features, apply_sigmoid, trees, iso_x, iso_y })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(MAGIC);
        b.write_u32::<LittleEndian>(FORMAT_VERSION).unwrap();
        b.write_u32::<LittleEndian>(self.n_features).unwrap();
        b.write_u32::<LittleEndian>(if self.apply_sigmoid { FLAG_SIGMOID } else { 0 }).unwrap();
        b.write_u32::<LittleEndian>(self.trees.len() as u32).unwrap();
        for t in &self.trees {
            let n = t.feature.len();
            b.write_u32::<LittleEndian>(n as u32).unwrap();
            write_i32_vec(&mut b, &t.feature);
            write_f32_vec(&mut b, &t.threshold);
            write_i32_vec(&mut b, &t.left);
            write_i32_vec(&mut b, &t.right);
            write_f32_vec(&mut b, &t.value);
            b.extend_from_slice(&t.default_left);
        }
        b.write_u32::<LittleEndian>(self.iso_x.len() as u32).unwrap();
        write_f32_vec(&mut b, &self.iso_x);
        write_f32_vec(&mut b, &self.iso_y);
        b
    }

    /// Calibrated P(signal) in [0,1].
    pub fn predict_proba(&self, x: &[f32]) -> f32 {
        let raw: f32 = self.trees.iter().map(|t| t.eval(x)).sum();
        let p = if self.apply_sigmoid { 1.0 / (1.0 + (-raw).exp()) } else { raw };
        self.isotonic(p)
    }

    /// The additive LLR term `log(s/(1-s))` with s clamped to (eps, 1-eps).
    pub fn predict_logit(&self, x: &[f32]) -> f32 {
        let s = self.predict_proba(x).clamp(PROB_EPS, 1.0 - PROB_EPS);
        (s / (1.0 - s)).ln()
    }

    /// Piecewise-linear interpolation of the isotonic map. Empty map => identity.
    fn isotonic(&self, p: f32) -> f32 {
        let n = self.iso_x.len();
        if n == 0 {
            return p;
        }
        if p <= self.iso_x[0] {
            return self.iso_y[0];
        }
        if p >= self.iso_x[n - 1] {
            return self.iso_y[n - 1];
        }
        // binary search for the segment containing p
        let mut lo = 0usize;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if self.iso_x[mid] <= p { lo = mid; } else { hi = mid; }
        }
        let (x0, x1) = (self.iso_x[lo], self.iso_x[hi]);
        let (y0, y1) = (self.iso_y[lo], self.iso_y[hi]);
        if (x1 - x0).abs() < f32::EPSILON {
            return y0;
        }
        y0 + (y1 - y0) * (p - x0) / (x1 - x0)
    }
}

fn read_i32_vec(c: &mut Cursor<&[u8]>, n: usize) -> Result<Vec<i32>, GbdtError> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n { v.push(c.read_i32::<LittleEndian>()?); }
    Ok(v)
}
fn read_f32_vec(c: &mut Cursor<&[u8]>, n: usize) -> Result<Vec<f32>, GbdtError> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n { v.push(c.read_f32::<LittleEndian>()?); }
    Ok(v)
}
fn write_i32_vec(b: &mut Vec<u8>, v: &[i32]) {
    for &x in v { b.write_i32::<LittleEndian>(x).unwrap(); }
}
fn write_f32_vec(b: &mut Vec<u8>, v: &[f32]) {
    for &x in v { b.write_f32::<LittleEndian>(x).unwrap(); }
}
```

Add `pub mod gbdt_eval;` to `crates/scoring/src/lib.rs`. Confirm `byteorder` already includes `WriteBytesExt` (it does — same crate as the existing `ReadBytesExt` import in `param_model.rs`).

- [ ] **Step 4: Run the tests**

Run: `cargo test -p scoring --lib gbdt_eval`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add crates/scoring/src/gbdt_eval.rs crates/scoring/src/lib.rs
git commit -m "feat(scoring): SoA GBDT tree walker + AGBD v1 blob format"
```

---

### Task A2.3: Carry the GBDT model on `Param`

**Files:**
- Modify: `crates/scoring/src/param_model.rs` (struct field + every constructor/`rebuild_cache` site)
- Test: extend the inline tests in `param_model.rs`

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)] mod tests` in `param_model.rs` (find the existing test module; if `tiny_param()` lives in `testutil`, this still compiles there):

```rust
#[test]
fn param_defaults_gbdt_model_to_none() {
    let p = crate::testutil::tiny_param();
    assert!(p.gbdt_peak_model.is_none(), "fresh param must carry no GBDT model");
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p scoring --lib param_defaults_gbdt_model_to_none`
Expected: FAIL — `no field gbdt_peak_model on type Param`.

- [ ] **Step 3: Add the field and default it everywhere**

In `crates/scoring/src/param_model.rs`:
1. Add the import near the top: `use crate::gbdt_eval::GbdtPeakModel;`
2. Add the field to the `Param` struct (after `partition_ion_types_cache`):

```rust
    /// Optional peptide-agnostic GBDT per-peak signal/noise model. Populated by
    /// the store reader from the manifest row's `gbdt_model_bytes` blob; `None`
    /// for legacy stores and for any slug without a trained GBDT (scoring is
    /// then byte-identical to the pre-GBDT engine).
    pub gbdt_peak_model: Option<GbdtPeakModel>,
```

3. Every place that builds a `Param { .. }` literal must set `gbdt_peak_model: None`. Find them:

```bash
grep -rn 'Param {' crates/scoring/src crates/scoring/tests crates/model-train/src | grep -v 'fn '
```
Add `gbdt_peak_model: None,` to each struct literal (the binary `.param` loader `load_from_bytes`, `tiny_param()` in `testutil`, and any test builders). The Parquet reader in `model-train` (`reconstruct_param`, Task A2.4) will set it to the decoded blob instead.

4. `Param` derives `PartialEq`; `GbdtPeakModel`/`Tree` already derive `PartialEq` (Task A2.2), so this still compiles.

- [ ] **Step 4: Run the test + the existing param tests**

Run: `cargo test -p scoring --lib param`
Expected: PASS (new test + existing param tests unaffected).

- [ ] **Step 5: Commit**

```bash
git add crates/scoring/src/param_model.rs
git commit -m "feat(scoring): Param carries optional GbdtPeakModel (defaults None)"
```

---

### Task A2.4: Store the GBDT blob in `models.parquet` (schema + write + read)

**Files:**
- Modify: `crates/model-train/src/store/schema.rs` (`combined_schema()`, ~157–230)
- Modify: `crates/model-train/src/store/write.rs` (`build_manifest_batch`, ~161–312; `write_models`, ~73–113)
- Modify: `crates/model-train/src/store/read.rs` (`reconstruct_param`, ~221–519; `parse_manifest_row`, ~521–588)
- Test: `crates/model-train/tests/roundtrip.rs` (add a GBDT-blob round-trip test)

- [ ] **Step 1: Write the failing round-trip test**

Add to `crates/model-train/tests/roundtrip.rs`:

```rust
#[test]
fn gbdt_blob_roundtrips_through_store() {
    use scoring_crate::gbdt_eval::{GbdtPeakModel, Tree};
    // A minimal model identical in shape to gbdt_eval's toy model.
    let model = GbdtPeakModel {
        n_features: 1,
        apply_sigmoid: true,
        trees: vec![Tree {
            feature: vec![0, -1, -1],
            threshold: vec![0.5, 0.0, 0.0],
            left: vec![1, -1, -1],
            right: vec![2, -1, -1],
            value: vec![0.0, -1.0, 2.0],
            default_left: vec![1, 1, 1],
        }],
        iso_x: vec![0.0, 1.0],
        iso_y: vec![0.0, 1.0],
    };
    let blob = model.to_bytes();

    let tmp = tempfile::tempdir().unwrap();
    let store = tmp.path().join("with_gbdt.parquet");

    // Build a tiny Param + attach the blob, write, reopen, assert the model
    // came back and predicts identically.
    let mut param = scoring_crate::testutil::tiny_param();
    param.gbdt_peak_model = None; // the blob is supplied to the writer separately
    model_train::store::write_models_with_gbdt(
        &store,
        &[("toy", &param, Some(blob.clone()))],
    )
    .expect("write");

    let loaded = model_train::store::ModelStore::open(&store)
        .unwrap()
        .load_param("toy")
        .expect("load toy");
    let gm = loaded.gbdt_peak_model.expect("gbdt model present after roundtrip");
    assert!((gm.predict_logit(&[1.0]) - 2.0).abs() < 1e-4);
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p model-train --test roundtrip gbdt_blob_roundtrips_through_store`
Expected: FAIL — `write_models_with_gbdt` undefined / `gbdt_peak_model` not populated.

- [ ] **Step 3: Add the schema column**

In `crates/model-train/src/store/schema.rs`, inside the `combined_schema()` field list (next to the other nullable manifest columns like `ion_loss_class`), add:

```rust
        nf("gbdt_model_bytes", DataType::Binary),
```

- [ ] **Step 4: Write the blob**

In `crates/model-train/src/store/write.rs`:
1. Add a `BinaryBuilder` column to `build_manifest_batch`. Change the function signature so manifest building accepts an optional blob per model. Concretely, add a parallel slice param `gbdt_blobs: &[Option<Vec<u8>>]` and build the column:

```rust
    use arrow::array::BinaryBuilder;
    let mut gbdt_b = BinaryBuilder::new();
    for blob in gbdt_blobs {
        match blob {
            Some(bytes) => gbdt_b.append_value(bytes),
            None => gbdt_b.append_null(),
        }
    }
    // ...push Arc::new(gbdt_b.finish()) as the column aligned with "gbdt_model_bytes"
```
Ensure the column is pushed in the SAME positional order as the field appears in `combined_schema()` (the manifest batch builds columns to match the schema order — append it at the matching position; the manifest-only columns are followed by null-filled table columns, so place `gbdt_model_bytes` exactly where the schema lists it).

2. Add the public entry point used by the test and the assembler:

```rust
/// Like `write_models`, but each model may carry a transcoded GBDT blob
/// (`AGBD` bytes) stored on its manifest row. `None` => null column (legacy).
pub fn write_models_with_gbdt(
    path: &Path,
    models: &[(&str, &Param, Option<Vec<u8>>)],
) -> Result<Vec<String>, TrainError> {
    let plain: Vec<(&str, &Param)> = models.iter().map(|(id, p, _)| (*id, *p)).collect();
    let blobs: Vec<Option<Vec<u8>>> = models.iter().map(|(_, _, b)| b.clone()).collect();
    write_models_inner(path, &plain, &blobs)
}
```
Refactor the existing `write_models` to delegate: `write_models(path, models)` calls `write_models_inner(path, models, &vec![None; models.len()])`. `write_models_inner` is the old body with `build_manifest_batch(schema, models, blobs)` now taking the blob slice. Re-export `write_models_with_gbdt` from `crates/model-train/src/store/mod.rs` (next to the existing `write_models` re-export).

- [ ] **Step 5: Read the blob back onto `Param`**

In `crates/model-train/src/store/read.rs`:
1. In `parse_manifest_row` (or wherever the manifest row is decoded into `ManifestRow`), read the optional binary column:

```rust
    use arrow::array::BinaryArray;
    let gbdt_bytes: Option<Vec<u8>> = batch
        .column_by_name("gbdt_model_bytes")
        .and_then(|c| c.as_any().downcast_ref::<BinaryArray>())
        .filter(|arr| !arr.is_null(i))
        .map(|arr| arr.value(i).to_vec());
```
Carry `gbdt_bytes` on `ManifestRow` (add a field `gbdt_bytes: Option<Vec<u8>>`).

2. In `reconstruct_param`, after the `Param { .. }` literal is built and `rebuild_cache()` runs, decode + attach:

```rust
    if let Some(bytes) = manifest.gbdt_bytes.as_ref() {
        param.gbdt_peak_model = Some(
            scoring_crate::gbdt_eval::GbdtPeakModel::from_bytes(bytes)
                .map_err(|e| TrainError::Other(format!("decode gbdt_model_bytes: {e}")))?,
        );
    }
```
(`reconstruct_param` builds `param` with `gbdt_peak_model: None` from Task A2.3; this overrides it when a blob exists. A legacy store has a null column → stays `None`.)

- [ ] **Step 6: Run the round-trip test + the legacy-store test**

Run:
```bash
cargo test -p model-train --test roundtrip
cargo test -p scoring --test param_loads_all_bundled
```
Expected: PASS — new blob round-trip passes AND the existing bundled store (which has no `gbdt_model_bytes`, or a null one) still loads with `gbdt_peak_model == None`.

- [ ] **Step 7: Commit**

```bash
git add crates/model-train/src/store
git commit -m "feat(store): nullable gbdt_model_bytes column with write/read round-trip"
```

---

### Task A2.5: Fold the GBDT term into scoring (additive, parity-safe)

**Files:**
- Modify: `crates/scoring/src/scoring/scored_spectrum.rs` (`ScoredSpectrum` struct + `new` + `directional_node_score_inner` + `directional_node_score` + `rank_kept`/`new_without_filtering`)
- Test: inline `#[cfg(test)]` in `scored_spectrum.rs` + run the parity golden

The per-peak GBDT logit is computed once in `new()` over the active peak list and stored in a **rank-indexed** vector `gbdt_logit_by_rank` (index = intensity rank, since rank↔active-peak is a bijection). The closure in `directional_node_score_inner` already receives the matched peak's `rank`; it adds `gbdt_logit_by_rank[rank]` to each matched ion's score. Empty vector (no model) ⇒ adds nothing ⇒ byte-identical.

- [ ] **Step 1: Write the failing test (signal peak gains score; absent model unchanged)**

Add to the `#[cfg(test)] mod tests` in `scored_spectrum.rs`:

```rust
#[test]
fn gbdt_term_raises_node_score_for_signal_peak() {
    use crate::gbdt_eval::{GbdtPeakModel, Tree};
    use crate::param_model::IonType;
    // Build a param whose prefix ion at rank 1 scores some base LLR, then attach
    // a GBDT model that always returns a large positive logit, and assert the
    // matched-node score increases by ~that logit.
    let mut param = tiny_param();
    // Constant-output GBDT: single leaf value 0 -> sigmoid(0)=0.5; instead use a
    // leaf that yields a big logit: leaf value 4.0, sigmoid(4)=0.982, iso identity
    // -> logit ~= 4.0.
    let gbdt = GbdtPeakModel {
        n_features: 1,
        apply_sigmoid: true,
        trees: vec![Tree {
            feature: vec![-1],
            threshold: vec![0.0],
            left: vec![-1],
            right: vec![-1],
            value: vec![4.0],
            default_left: vec![1],
        }],
        iso_x: vec![],
        iso_y: vec![],
    };

    // A spectrum with one peak placed at the theoretical prefix m/z of a node.
    let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
    let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
    let theo_mz = ion.mz(300.0); // nominal mass 300
    let spec = Spectrum {
        // minimal: one strong peak at theo_mz so it ranks 1
        peaks: vec![(theo_mz, 1000.0_f32)],
        precursor_mz: 751.0,
        precursor_charge: Some(2),
        ..Spectrum::default()
    };

    let scorer_plain = RankScorer::new(&param);
    let s_plain = ScoredSpectrum::new(&spec, &scorer_plain, 2);
    let base = s_plain.directional_node_score(
        300.0, true, &scorer_plain, 2, s_plain.parent_mass, 0.5,
    );

    param.gbdt_peak_model = Some(gbdt);
    let scorer_gbdt = RankScorer::new(&param);
    let s_gbdt = ScoredSpectrum::new(&spec, &scorer_gbdt, 2);
    let with = s_gbdt.directional_node_score(
        300.0, true, &scorer_gbdt, 2, s_gbdt.parent_mass, 0.5,
    );

    // The matched prefix ion at rank 1 now carries an extra ~+4.0 logit.
    assert!(with > base + 3.0, "expected GBDT bump: base={base}, with={with}");
}
```
(If `Spectrum`/`Partition`/`IonType` import paths differ in the test module, mirror the imports already used by the neighbouring `node_score_*` tests in this file.)

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p scoring --lib gbdt_term_raises_node_score_for_signal_peak`
Expected: FAIL — currently the GBDT model is ignored, so `with == base`.

- [ ] **Step 3: Add the cache field**

In the `ScoredSpectrum` struct (around line 181, after `observed_mass_cache`), add:

```rust
    /// Per-rank GBDT signal/noise LLR term `log(s/(1-s))`. Index = intensity
    /// rank (1-based); `gbdt_logit_by_rank[r]` is the term for the active peak
    /// whose rank is `r`. Empty when the model carries no GBDT (then the term
    /// is never added and scoring is byte-identical). Index 0 is unused.
    gbdt_logit_by_rank: Vec<f32>,
```

- [ ] **Step 4: Populate it in `new()` (over the active peak list)**

In `ScoredSpectrum::new`, AFTER `cache_peaks`/`cache_ranks` are chosen (line ~409–413) and BEFORE the prefix/suffix cache-fill loop (line 414), compute the term:

```rust
        // Per-peak GBDT signal/noise logit, computed ONCE over the active peak
        // list and indexed by intensity rank so the cache-fill loop and the
        // fallback path can add it per matched ion. Empty (no-op) when the
        // model has no GBDT — preserving byte-identical scoring.
        let gbdt_logit_by_rank: Vec<f32> = match scorer.param().gbdt_peak_model.as_ref() {
            None => Vec::new(),
            Some(model) => {
                use crate::peak_features::{extract_peak_features, PeakFeatureCtx};
                let max_r = cache_ranks.iter().copied().filter(|&r| r != u32::MAX).max().unwrap_or(0);
                let mut by_rank = vec![0.0_f32; (max_r as usize) + 1];
                let base_peak_intensity = cache_peaks.iter().map(|&(_, i)| i).fold(0.0_f32, f32::max);
                let ctx = PeakFeatureCtx {
                    precursor_mz: spec.precursor_mz,
                    charge,
                    parent_neutral_mass: parent_mass,
                    total_intensity,
                    base_peak_intensity,
                    window_da: crate::peak_features::FEATURE_WINDOW_DA,
                    // Per-spectrum partner-match tolerance; `build_dataset.py`
                    // MUST pass the identical value at training time (A2.10).
                    match_tol_da: param.mme.as_da(parent_mass.max(1.0)),
                };
                let feats = extract_peak_features(cache_peaks, cache_ranks, &ctx);
                for (i, &r) in cache_ranks.iter().enumerate() {
                    if r != u32::MAX && (r as usize) < by_rank.len() {
                        by_rank[r as usize] = model.predict_logit(&feats[i]);
                    }
                }
                by_rank
            }
        };
```

Add `gbdt_logit_by_rank,` to the `Self { .. }` constructor literal at the end of `new()` (line ~447–462). Add `gbdt_logit_by_rank: Vec::new(),` to the `rank_kept` `Self { .. }` literal (line ~528–545) so `new_without_filtering` compiles.

- [ ] **Step 5: Thread the slice into the scoring inner loop**

`directional_node_score_inner` is a static fn; add a `gbdt_logit_by_rank: &[f32]` parameter and add the term in the closure. Change the signature and the closure body (around lines 840–877):

```rust
    fn directional_node_score_inner(
        peaks: &[(f64, f32)],
        ranks: &[u32],
        segment_partition_cache: SegmentPartitionSlice<'_>,
        scorer: &RankScorer,
        nominal_mass: f64,
        is_prefix: bool,
        charge: u8,
        parent_mass: f64,
        gbdt_logit_by_rank: &[f32],
    ) -> f32 {
        let max_rank = scorer.max_rank();
        let max_rank_idx = max_rank as usize;
        let mut total = 0.0_f32;
        visit_directional_node_ion_matches(
            peaks, ranks, segment_partition_cache, scorer,
            nominal_mass, is_prefix, charge, parent_mass, false,
            |_, _, rank, logs, _, _| {
                let score = match rank {
                    Some(rank) => {
                        let idx = rank.min(max_rank).max(1) as usize - 1;
                        let base = if idx < logs.len() { logs[idx] } else { 0.0 };
                        // Additive per-peak GBDT term for the MATCHED peak.
                        let gbdt = gbdt_logit_by_rank.get(rank as usize).copied().unwrap_or(0.0);
                        base + gbdt
                    }
                    None => {
                        if max_rank_idx < logs.len() { logs[max_rank_idx] } else { 0.0 }
                    }
                };
                total += score;
            },
        );
        total
    }
```
Note: the GBDT term is added ONLY on a match (`Some(rank)`), never on the absent-ion slot — a missing ion has no observed peak, so no signal probability applies.

- [ ] **Step 6: Pass the slice at both call sites**

First find EVERY caller of the inner fn so none are missed:

```bash
grep -rn 'directional_node_score_inner' crates/scoring/src
```
Expected callers: the two `new()` cache-fill calls (lines ~416, ~426) and `directional_node_score` (line ~827). In `new()`'s cache-fill loop pass `&gbdt_logit_by_rank` as the new last argument to both calls. In `directional_node_score` pass `&self.gbdt_logit_by_rank`:

```rust
        Self::directional_node_score_inner(
            peaks, ranks, &self.segment_partition_cache, scorer,
            nominal_mass, is_prefix, charge, parent_mass,
            &self.gbdt_logit_by_rank,
        )
```

If the grep turns up any other caller, pass `&[]` there if it has no spectrum-rank context, or the appropriate `gbdt_logit_by_rank` slice if it does — the empty slice is the no-op (byte-identical) path.

- [ ] **Step 7: Run the new test + the parity golden**

Run:
```bash
cargo test -p scoring --lib gbdt_term_raises_node_score_for_signal_peak
cargo test -p andes --test precursor_cal_bit_identical
```
Expected: the GBDT test PASSES; the parity golden `precursor_cal_off_pin_tsv_match_golden_after_sort` STILL PASSES (the bundled store carries no GBDT blob, so `gbdt_logit_by_rank` is empty and scoring is byte-identical).

- [ ] **Step 8: Full scoring + workspace test sweep**

Run:
```bash
cargo test -p scoring
cargo build --workspace && cargo clippy --workspace -- -D warnings
```
Expected: PASS (modulo the pre-existing `integer_mass_scaler` failure noted in memory). Clippy clean.

- [ ] **Step 9: Commit**

```bash
git add crates/scoring/src/scoring/scored_spectrum.rs
git commit -m "feat(scoring): additive per-peak GBDT LLR term, byte-identical when absent"
```

---

## Phase A2 (Rust) — the GBDT producer (training)

**Rust-only (decision 6/7/8 in the spec, 2026-06-15).** No Python. The GBDT is
trained in `crates/model-train` and emitted directly as a `GbdtPeakModel` (the
SoA walker type from A2.2) — no LightGBM, no `model.txt`, no transcode, no
cross-language parity gate. There is ONE feature extractor
(`scoring::peak_features`) and ONE ion model (the engine's), shared by training
and inference: the dataset builder constructs its `PeakFeatureCtx` via
`PeakFeatureCtx::for_spectrum` (A2.5 follow-up), so train-time and inference-time
features are identical **by construction**.

New module tree: `crates/model-train/src/gbdt/{mod.rs, labels.rs, isotonic.rs,
tree.rs, train.rs}`. `model-train` already depends on `scoring` (as
`scoring_crate`) and `model`, so `GbdtPeakModel`/`Tree`, `extract_peak_features`,
`PeakFeatureCtx`, and residue masses are all reachable.

### Task A2.6′: Rust theoretical-ion labeling

**Files:**
- Create: `crates/model-train/src/gbdt/mod.rs` (declares the submodules: `pub mod labels; pub mod isotonic; pub mod tree; pub mod train;`)
- Create: `crates/model-train/src/gbdt/labels.rs`
- Modify: `crates/model-train/src/lib.rs` (add `pub mod gbdt;`)
- Test: inline `#[cfg(test)]` in `labels.rs`

Label each active peak signal(1)/noise(0): signal = the peak matches ANY
generously-enumerated theoretical ion (b/y, +a, −H2O, −NH3) of the confident
peptide(s) within the model's fragment tolerance, at charges 1..=max(1,z−1).
For chimeric IDs, label against the UNION of confident peptides (so a peak
explained by either is signal). Co-isolation exclusion happens upstream (the
dataset builder decides which peptides feed this).

Use residue masses from the `model` crate (single source of truth) — find the
residue-mass accessor (`model::amino_acid` / `model::mass`); do NOT hardcode a
second AA table. Use `model::mass::{PROTON, H2O}` and the local `NH3 = 17.026549`,
`CO = 27.994915` (no model equivalents). I/L are isobaric (same mass) — fine.

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn b_ion_peak_labeled_signal() {
        // "PEPTIDE" b2 (P+E) singly charged ≈ 227.1026. Place a peak there and a
        // junk peak; expect [signal, noise].
        let theo = theoretical_ion_mzs("PEPTIDE", 2);
        assert!(theo.iter().any(|m| (m - 227.1026).abs() < 0.02), "b2 not enumerated: {theo:?}");
        let peaks = [(227.1026_f64, 500.0_f32), (999.999, 10.0)];
        let labels = label_peaks(&peaks, &["PEPTIDE"], 2, 0.02);
        assert_eq!(labels, vec![1u8, 0]);
    }

    #[test]
    fn chimeric_union_labels_either_peptide() {
        // A peak matching peptide B but not A is still signal when both are confident.
        let a = "PEPTIDE";
        let b = "SAMPLER";
        let b_theo = theoretical_ion_mzs(b, 2);
        let only_b = b_theo[0]; // some b/y ion of B
        let peaks = [(only_b, 100.0_f32)];
        let labels = label_peaks(&peaks, &[a, b], 2, 0.02);
        assert_eq!(labels, vec![1u8]);
    }
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p model-train --lib gbdt::labels`
Expected: FAIL — `theoretical_ion_mzs`/`label_peaks` undefined.

- [ ] **Step 3: Implement `labels.rs`**

Signature contract (fill the body using the `model` crate's residue masses):
```rust
//! Theoretical-ion oracle for GBDT signal/noise labels. Reuses the engine's
//! residue-mass table (the `model` crate) so labels are defined by the same ion
//! chemistry the engine scores. Enumerates b/y (+a, −H2O, −NH3) at charges
//! 1..=max(1, z-1).
use model::mass::{PROTON, H2O};

const NH3: f64 = 17.026549;
const CO: f64 = 27.994915;

/// Sorted theoretical ion m/z for `peptide` at precursor charge `z`. Peptide is
/// an uppercase AA string (mods ignored — labeling is structural). Unknown chars
/// are skipped. Returns ascending m/z (b/y + a + water/ammonia losses).
pub fn theoretical_ion_mzs(peptide: &str, z: u8) -> Vec<f64> {
    // residue masses via the model crate (e.g. model::amino_acid::residue_mass(c)
    // or model::mass table) — use whatever the crate exposes; do not hardcode.
    // total = sum(res) + H2O; for i in 1..n: prefix += res[i-1];
    //   b = prefix + PROTON; y = (total - prefix) + PROTON;
    //   for base in [b, y], for zz in 1..=max(1,z-1):
    //     push (base + (zz-1)*PROTON)/zz, and the −H2O and −NH3 variants;
    //   a = b - CO at each charge.
    // sort ascending, return.
    todo!("enumerate as specified")
}

/// 1 if a peak matches any theoretical ion of any confident peptide within
/// `tol_da`, else 0. `peptides` is the confident set (union for chimeric).
pub fn label_peaks(peaks: &[(f64, f32)], peptides: &[&str], z: u8, tol_da: f64) -> Vec<u8> {
    let mut theo: Vec<f64> = peptides.iter().flat_map(|p| theoretical_ion_mzs(p, z)).collect();
    theo.sort_by(|a, b| a.partial_cmp(b).unwrap());
    peaks.iter().map(|&(mz, _)| {
        let lo = theo.partition_point(|&t| t < mz - tol_da);
        u8::from(lo < theo.len() && theo[lo] <= mz + tol_da)
    }).collect()
}
```
(Replace the two `todo!`/comment bodies with the real enumeration. The `label_peaks` body is complete.)

- [ ] **Step 4: Run the tests** — `cargo test -p model-train --lib gbdt::labels` → PASS.

- [ ] **Step 5: Commit**
```bash
git add crates/model-train/src/gbdt/mod.rs crates/model-train/src/gbdt/labels.rs crates/model-train/src/lib.rs
git commit -m "feat(model-train): theoretical-ion labeling for GBDT signal/noise"
```

---

### Task A2.7′a: PAVA isotonic calibration

**Files:**
- Create: `crates/model-train/src/gbdt/isotonic.rs`
- Test: inline `#[cfg(test)]`

Pool-Adjacent-Violators producing the `(iso_x, iso_y)` breakpoints
`GbdtPeakModel` consumes (ascending x, monotone-nondecreasing y). Calibrates raw
predicted probability → empirical signal fraction.

- [ ] **Step 1: Failing test**
```rust
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn pava_is_monotone_and_fits_steps() {
        // x sorted; y has a violation that PAVA must pool.
        let x = vec![0.1_f32, 0.2, 0.3, 0.4, 0.5];
        let y = vec![0.0_f32, 1.0, 0.0, 1.0, 1.0]; // non-monotone
        let (bx, by) = pava(&x, &y);
        assert_eq!(bx.len(), by.len());
        // output must be nondecreasing
        for w in by.windows(2) { assert!(w[1] >= w[0] - 1e-6, "not monotone: {by:?}"); }
        // breakpoints span the input x-range
        assert!((bx[0] - 0.1).abs() < 1e-6 && (*bx.last().unwrap() - 0.5).abs() < 1e-6);
    }
}
```

- [ ] **Step 2: Run → FAIL** (`pava` undefined). `cargo test -p model-train --lib gbdt::isotonic`

- [ ] **Step 3: Implement PAVA**
```rust
//! Pool-Adjacent-Violators isotonic regression → (x, y) breakpoints for
//! `GbdtPeakModel` calibration. Inputs MUST be sorted ascending by x.
/// Returns (xs, ys): ascending xs, nondecreasing ys, one point per pooled block
/// boundary. Empty input → empty (the walker treats an empty iso map as identity).
pub fn pava(x: &[f32], y: &[f32]) -> (Vec<f32>, Vec<f32>) {
    assert_eq!(x.len(), y.len());
    if x.is_empty() { return (Vec::new(), Vec::new()); }
    // Blocks of (sum_y, weight, right_x). Merge while previous mean > current.
    let mut means: Vec<f32> = Vec::new();
    let mut weights: Vec<f32> = Vec::new();
    let mut right_x: Vec<f32> = Vec::new();
    for i in 0..x.len() {
        means.push(y[i]);
        weights.push(1.0);
        right_x.push(x[i]);
        while means.len() > 1 && means[means.len()-2] > means[means.len()-1] {
            let n = means.len();
            let w = weights[n-2] + weights[n-1];
            let m = (means[n-2]*weights[n-2] + means[n-1]*weights[n-1]) / w;
            means[n-2] = m; weights[n-2] = w; right_x[n-2] = right_x[n-1];
            means.pop(); weights.pop(); right_x.pop();
        }
    }
    // Emit one (x, y) per block at the block's right edge; prepend the first x
    // so the map covers the full input domain.
    let mut bx = vec![x[0]];
    let mut by = vec![means[0]];
    for b in 0..means.len() {
        bx.push(right_x[b]);
        by.push(means[b]);
    }
    (bx, by)
}
```

- [ ] **Step 4: Run → PASS.**

- [ ] **Step 5: Commit**
```bash
git add crates/model-train/src/gbdt/isotonic.rs
git commit -m "feat(model-train): PAVA isotonic calibration for GBDT"
```

---

### Task A2.7′b: histogram regression tree (one boosting weak learner)

**Files:**
- Create: `crates/model-train/src/gbdt/tree.rs`
- Test: inline `#[cfg(test)]`

Fits ONE depth-limited regression tree to per-row gradients/hessians using
pre-binned features, emitting a `gbdt_eval::Tree` (SoA, pre-order node indices so
A2.2's `validate` passes: each internal node's children have a strictly greater
index). Newton split gain: `G²/(H+λ)`; leaf weight `−G/(H+λ)`.

- [ ] **Step 1: Failing test** — a tree fit to gradients that separate on one binned feature splits on that feature and assigns opposite-sign leaf weights.
```rust
#[cfg(test)]
mod tests {
    use super::*;
    use scoring_crate::gbdt_eval::Tree;
    #[test]
    fn tree_splits_on_discriminative_feature() {
        // 100 rows, 2 features. Feature 0 (binned) separates two groups whose
        // gradients have opposite sign; feature 1 is noise.
        let n = 100;
        let n_features = 2;
        let n_bins = 4;
        let mut binned = vec![0u8; n * n_features]; // row-major [row*nf + f]
        let mut grad = vec![0.0_f32; n];
        let hess = vec![1.0_f32; n];
        for r in 0..n {
            let hi = r % 2 == 0;
            binned[r*n_features + 0] = if hi { 3 } else { 0 };  // feature 0 separates
            binned[r*n_features + 1] = (r % n_bins) as u8;      // feature 1 noise
            grad[r] = if hi { -1.0 } else { 1.0 };              // opposite gradients
        }
        let params = TreeParams { max_depth: 2, n_bins, lambda: 1.0, min_hessian: 1.0, n_features };
        // bin_uppers maps (feature, bin) -> upper threshold; here just bin index as f32.
        let bin_uppers = vec![vec![0.5_f32, 1.5, 2.5, 3.5]; n_features];
        let tree: Tree = fit_tree(&binned, &grad, &hess, &params, &bin_uppers);
        // root must split on feature 0
        assert_eq!(tree.feature[0], 0, "expected split on feature 0, got {:?}", tree.feature);
        // a tree with at least one split has >=3 nodes
        assert!(tree.feature.len() >= 3);
    }
}
```

- [ ] **Step 2: Run → FAIL.** `cargo test -p model-train --lib gbdt::tree`

- [ ] **Step 3: Implement the tree**

Define:
```rust
use scoring_crate::gbdt_eval::Tree;

pub struct TreeParams {
    pub max_depth: u32,
    pub n_bins: usize,
    pub lambda: f32,      // L2 on leaf weights
    pub min_hessian: f32, // min summed hessian per child (≈ min_child_weight)
    pub n_features: usize,
}
```
Algorithm (depth-first; assign node indices in pre-order so children index > parent — required by `Tree::validate`):
- For the current node's row set, accumulate per-(feature,bin) gradient/hessian histograms.
- For each feature, scan bins left→right accumulating (GL, HL); right side is (G−GL, H−HL). Split gain at the bin boundary = `GL²/(HL+λ) + GR²/(HR+λ) − G²/(H+λ)`; require `HL≥min_hessian` and `HR≥min_hessian`. Track the best (feature, bin, gain>0).
- If no positive-gain split or depth==max_depth: make a LEAF with weight `−G/(H+λ)`.
- Else: internal node at this index; recurse left (rows with `binned[f] ≤ bin`) then right, appending children after the parent (pre-order). `threshold` = `bin_uppers[feature][bin]` (the upper edge of the left bin), so inference's `x <= threshold` reproduces the binned split. `default_left = 1`.

Emit into the SoA `Tree { feature, threshold, left, right, value, default_left }`.
Build with a recursive helper that pushes the parent, reserves child slots, fills
recursively, and back-patches `left`/`right` to the recorded child indices (mirror
A2.9 transcode's pre-order assignment from the deleted Python — but native here).
Keep helper fns small and tested.

- [ ] **Step 4: Run → PASS.**

- [ ] **Step 5: Commit**
```bash
git add crates/model-train/src/gbdt/tree.rs
git commit -m "feat(model-train): histogram regression tree (GBDT weak learner)"
```

---

### Task A2.7′c: gradient boosting + assembly into `GbdtPeakModel`

**Files:**
- Create: `crates/model-train/src/gbdt/train.rs`
- Test: inline `#[cfg(test)]`

Boost `T` trees on the logistic loss (Newton steps), with negative
undersampling, a group-disjoint validation split, early stopping, and PAVA
calibration on the validation set. Output a `GbdtPeakModel { n_features,
apply_sigmoid: true, trees, iso_x, iso_y }`.

- [ ] **Step 1: Failing test (separable synthetic → high AUC + signal>noise logit)**
```rust
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn separable_data_trains_and_separates() {
        // feature 0 ~ class; rest noise. Deterministic LCG (no rand dep).
        let n = 4000; let nf = 6usize;
        let mut x = vec![0.0_f32; n*nf];
        let mut y = vec![0u8; n];
        let mut s: u64 = 0x1234_5678;
        let mut next = || { s = s.wrapping_mul(6364136223846793005).wrapping_add(1); ((s>>33) as f32)/(u32::MAX as f32) };
        for r in 0..n {
            for f in 0..nf { x[r*nf+f] = next(); }
            let lab = if x[r*nf] + 0.2*next() > 0.6 { 1u8 } else { 0 };
            y[r] = lab;
        }
        let groups: Vec<u32> = (0..n as u32).map(|r| r % 50).collect();
        let params = TrainParams::default();
        let ds = Dataset { x, y, groups, n_features: nf };
        let model = train_gbdt(&ds, &params, 0xC0FFEE);
        assert!(model.apply_sigmoid);
        assert_eq!(model.n_features as usize, nf);
        // a clearly-positive row scores higher logit than a clearly-negative row
        let mut hi = vec![0.0_f32; nf]; hi[0] = 0.95;
        let mut lo = vec![0.0_f32; nf]; lo[0] = 0.05;
        assert!(model.predict_logit(&hi) > model.predict_logit(&lo) + 1.0,
            "hi={} lo={}", model.predict_logit(&hi), model.predict_logit(&lo));
    }
}
```

- [ ] **Step 2: Run → FAIL.** `cargo test -p model-train --lib gbdt::train`

- [ ] **Step 3: Implement boosting**

Interfaces:
```rust
use scoring_crate::gbdt_eval::GbdtPeakModel;
use super::tree::{fit_tree, TreeParams};
use super::isotonic::pava;

/// Row-major feature matrix `x` (`n × n_features`), binary labels `y`, and a
/// `groups` id per row (run+peptide) for leakage-free validation splitting.
pub struct Dataset { pub x: Vec<f32>, pub y: Vec<u8>, pub groups: Vec<u32>, pub n_features: usize }

pub struct TrainParams {
    pub n_rounds: u32,        // max boosting rounds (e.g. 300)
    pub learning_rate: f32,   // shrinkage (e.g. 0.05)
    pub max_depth: u32,       // 4..=7
    pub n_bins: usize,        // ≤256 (e.g. 64)
    pub lambda: f32,          // L2 (e.g. 1.0)
    pub min_hessian: f32,     // min child weight in hessian (e.g. 100*p(1-p)-ish, use ~5..50)
    pub neg_pos_ratio: f32,   // undersample negatives to this ratio (e.g. 4.0)
    pub val_fraction: f32,    // group fraction held out (e.g. 0.2)
    pub early_stop_rounds: u32,// stop if val logloss doesn't improve for N (e.g. 30)
}
impl Default for TrainParams { /* the values above */ }

pub fn train_gbdt(ds: &Dataset, p: &TrainParams, seed: u64) -> GbdtPeakModel { /* below */ }
```
Steps inside `train_gbdt`:
1. **Quantile-bin** each feature into ≤`n_bins` bins once: compute per-feature bin upper edges (`bin_uppers[f]`) from sorted sampled values; map every row to `binned[row*nf+f] = bin index`. (Uniform-quantile is fine; ties go to the lower bin via `<= upper`.)
2. **Group split**: deterministically (seeded) assign whole `groups` to train/val by `val_fraction` (hash the group id with the seed; no leakage — a group is wholly in one side).
3. **Negative undersample** the TRAIN rows to `neg_pos_ratio` (seeded), keeping all positives. Validation is NOT undersampled.
4. **Boost**: `raw[r] = 0`; per round: `p_r = sigmoid(raw[r])`, `grad = p_r − y`, `hess = p_r*(1−p_r)`; `fit_tree(train_binned, grad, hess, TreeParams{...}, bin_uppers)`; scale all leaf `value` by `learning_rate`; add tree; update `raw[r] += learning_rate * tree.eval(row_features_continuous)` — NOTE: evaluate the tree on the SAME continuous features the walker will see (use `bin_uppers` thresholds, so `tree.eval` on continuous x reproduces the binned decision). Track validation logloss each round; keep the best round; stop after `early_stop_rounds` without improvement; truncate `trees` to the best round.
5. **Calibrate**: on the VALIDATION set compute `p = sigmoid(raw_val)`; sort `(p, y)` by `p`; `pava(&sorted_p, &sorted_y)` → `(iso_x, iso_y)`.
6. Return `GbdtPeakModel { n_features: nf as u32, apply_sigmoid: true, trees, iso_x, iso_y }`.

Determinism: seed all sampling/splitting with the `seed` arg (LCG or a tiny seeded PRNG; do NOT use `rand` thread RNG — the workspace forbids nondeterministic `Date.now`/`Math.random`-style sources and we want reproducible models). No external crates needed.

- [ ] **Step 4: Run → PASS** (AUC-ish separation assertion).

- [ ] **Step 5: Commit**
```bash
git add crates/model-train/src/gbdt/train.rs
git commit -m "feat(model-train): logistic gradient-boosting trainer + calibration -> GbdtPeakModel"
```

---

### Task A2.8′: trainer quality + end-to-end-to-scoring integration test

**Files:**
- Test: `crates/model-train/tests/gbdt_trainer.rs`

Replaces the deleted cross-language parity gate. Proves: (a) determinism (same
seed → identical model bytes), (b) round-trip through the store, (c) the trained
model, attached to a `Param`, raises a matched peak's node score (ties A2.7′c →
A2.4 → A2.5 together).

- [ ] **Step 1: Write the test**
```rust
use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};

fn toy_dataset() -> Dataset { /* same separable LCG construction as A2.7'c, nf=18 */ }

#[test]
fn trainer_is_deterministic() {
    let ds = toy_dataset();
    let m1 = train_gbdt(&ds, &TrainParams::default(), 42);
    let m2 = train_gbdt(&ds, &TrainParams::default(), 42);
    assert_eq!(m1.to_bytes(), m2.to_bytes(), "same seed must yield identical model");
}

#[test]
fn trained_model_roundtrips_and_separates() {
    use scoring_crate::gbdt_eval::GbdtPeakModel;
    let ds = toy_dataset();
    let model = train_gbdt(&ds, &TrainParams::default(), 7);
    let back = GbdtPeakModel::from_bytes(&model.to_bytes()).expect("decode");
    assert_eq!(back, model);
    let mut hi = vec![0.0_f32; 18]; hi[0] = 0.95;
    let mut lo = vec![0.0_f32; 18]; lo[0] = 0.05;
    assert!(back.predict_logit(&hi) > back.predict_logit(&lo));
}
```
(The Param→scoring integration is already covered by A2.5's
`gbdt_term_raises_node_score_for_signal_peak`; this test focuses on the trainer's
own contract. Optionally add a test that writes the trained model via
`write_models_with_gbdt` and reloads a `Param` with it — reuse A2.4's pattern.)

- [ ] **Step 2: Run → PASS.** `cargo test -p model-train --test gbdt_trainer`

- [ ] **Step 3: Commit**
```bash
git add crates/model-train/tests/gbdt_trainer.rs
git commit -m "test(model-train): GBDT trainer determinism + round-trip + separation"
```

---

### Task A2.9′: unify `andes train` (source-agnostic, rank-core + GBDT in one pass)

**Files:**
- Modify: `crates/andes/src/bin/andes.rs` (the `Command` enum, `run_train`, `run_train_from_msnet`, `TrainFromMsnetArgs`)
- Modify: `TRAIN.md` (document `andes train` + `--source` + `train-from-search`)
- Test: `crates/andes/tests/cli_smoke.rs` or a new `crates/andes/tests/train_cli.rs`

Resolve the CLI per spec decision 7:
1. **Rename the existing `train`** (`Command::Train` / `run_train`, "train from spectra + FASTA, self-labeling") to subcommand name **`train-from-search`** (keep the function, just rename the clap `#[command(name = "...")]` and the enum variant for clarity, e.g. `TrainFromSearch`/`run_train_from_search`).
2. **Promote `train-from-msnet`** (`Command::TrainFromMsnet` / `run_train_from_msnet`) to subcommand name **`train`** (rename variant → `Train`, fn → `run_train`, args → `TrainArgs`). Add a `--source <STR>` arg (default `msnet` for back-compat) and use it for the source-ledger tag that is currently hardcoded `"source 'msnet'"` (andes.rs:~2829) — so PRIDE/MSnet/quantms harvests are recorded with their real provenance.
3. **Add GBDT training to the promoted `train`**: after the rank-table estimation + before/at store write, build a per-peak dataset from the SAME input PSMs + spectra and train the GBDT, attaching its blob to the written model. Concretely:
   - For each confident PSM's spectrum, build the active peak list with the SAME prep the engine uses (reuse `ScoredSpectrum`/the windowed filter + ranks), construct `PeakFeatureCtx::for_spectrum(...)`, `extract_peak_features(...)`, and `label_peaks(...)` (A2.6′) against the PSM's peptide(s); accumulate rows into a `Dataset` (x, y, group = hash(run+peptide)).
   - `train_gbdt(&dataset, &TrainParams::default(), seed)` → `GbdtPeakModel`.
   - Write the store via `write_models_with_gbdt(out, &[(model_id, &param, Some(model.to_bytes()))])` (A2.4) instead of the plain `write_models`.
   - Gate behind a `--gbdt {auto,on,off}` flag (default `on` for `train`; `off` reproduces the rank-core-only store). When `off`, write the plain store (no blob) — keeps a fast path + an A/B lever for A4.
4. The `train-from-msnet` name MAY be kept as a hidden deprecated alias of `train` (clap `visible_alias`/hidden) so existing scripts don't break — optional; note it in TRAIN.md.

- [ ] **Step 1: Failing test (CLI shape)**
```rust
// crates/andes/tests/train_cli.rs
use std::process::Command;
#[test]
fn train_help_shows_source_and_gbdt() {
    let out = Command::new(env!("CARGO_BIN_EXE_andes")).args(["train","--help"]).output().unwrap();
    let h = String::from_utf8_lossy(&out.stdout);
    assert!(h.contains("--source"), "train --help must document --source:\n{h}");
    assert!(h.contains("--gbdt"), "train --help must document --gbdt:\n{h}");
}
#[test]
fn train_from_search_exists() {
    let out = Command::new(env!("CARGO_BIN_EXE_andes")).args(["train-from-search","--help"]).output().unwrap();
    assert!(out.status.success(), "train-from-search subcommand must exist");
}
```

- [ ] **Step 2: Run → FAIL.** `cargo test -p andes --test train_cli`

- [ ] **Step 3: Implement the rename + `--source`/`--gbdt` + the one-pass GBDT training.** Read `run_train_from_msnet` (the PSM-reading + estimation + store-write flow) and thread the dataset build + `train_gbdt` + `write_models_with_gbdt` in. Keep the rank-core path identical when `--gbdt off`.

- [ ] **Step 4: Verify**
```bash
cargo test -p andes --test train_cli
cargo build --workspace && cargo clippy --workspace -- -D warnings
```
Also run any existing train tests (`crates/andes/tests/train_end_to_end.rs`) — adapt them to the renamed subcommands; they must pass (or be updated to the new names).

- [ ] **Step 5: Commit**
```bash
git add crates/andes/src/bin/andes.rs crates/andes/tests/train_cli.rs TRAIN.md
git commit -m "feat(andes): unify `train` (source-agnostic, rank-core + GBDT one-pass); rename old train -> train-from-search"
```

---

## Phase A3 — source-agnostic harvest + reproducible training (Codon, REMOTE)

**Surface to the user — runs on Codon with real data (per the codon-cluster skill);
not auto-run here.** The dataset/harvest wiring is now just feeding `andes train`.

### Task A3.1: Normalized harvested-PSM input contract + per-source harvest
Define ONE input contract that all provenances normalize to (parquet/TSV):
`run_id, scan_id, peptide, charge, precursor_mz, coisolation, q_value`, plus the
spectra (mzML/MGF/peaks) for each run. Harvests:
- **PRIDE reanalysis** (quantms/SDRF pipelines) → export confident PSMs (q≤0.005) in the contract.
- **MSnet** → the existing `train-from-msnet` PSM parquets already fit (source `msnet`).
- **quantms** → export its PSM table to the contract (source `quantms`).
Co-isolation exclusion (drop high-co-isolation scans) and chimeric-union selection
happen here, before `andes train` consumes the table.
Acceptance: a smoke harvest for ONE slug produces a valid contract file; `andes
train --source <x> --gbdt on` on it writes a `models.parquet` whose loaded `Param`
has `gbdt_peak_model.is_some()`.

### Task A3.2: Per-slug clean retrain on Codon
Build `andes` on Codon (rsync `crates/`+`Cargo.{toml,lock}`; the VM/Codon repos are
source copies, not git — see memory). For each public own slug
(hcd_qexactive_tryp, cid_lowres_tryp, cid_lowres_tryp_tmt, hcd_highres_lysc, + the
GluC/LysC/ETD variants), run `andes train --inputs <harvested>.parquet
--seed-model <archived prior or slug> --model-id <slug> --source <provenance>
--gbdt on --out-store <slug>.parquet`. One pass fits rank tables + GBDT.

### Task A3.3: Assemble the own-only public store
Combine the per-slug stores into one `resources/models.parquet` (own slugs only).
A small `andes`/`model-train` utility or `--update` flow merges per-slug models
(each carrying its `gbdt_model_bytes`) into the bundled store. Acceptance: the
combined store opens locally, every shipped slug has `gbdt_peak_model.is_some()`,
and the parity golden / no-GBDT behavior is unaffected for any slug trained with
`--gbdt off`.

---

## Phase A4 — competitive-advantage gate (VM)

### Task A4.1: Benchmark own models vs the field at 1% true entrapment-FDP

**Files:**
- Uses (VM): `scripts/bench_astral_competitors.sh`, `scripts/bench_tmt_ups_competitors.sh`, `scripts/astral_entrapment_experiment.sh`, `scripts/entrap_fdp.py`
- Create (VM): `scripts/bench_lysc_competitors.sh` (LysC variant of the TMT/UPS competitor script, `--enzyme lysc`)

This is THE GATE. No Phase B PR until it passes. The benchmark must use the own-only `resources/models.parquet` assembled in A3.3 (Step 4) — NOT the MS-GF+-derived store.

- [ ] **Step 1: Deploy the own store + andes binary to the VM**

Copy the own-only `resources/models.parquet` to the VM and build/deploy the andes binary there (the VM repo is a source copy — rsync `crates/`+`Cargo.{toml,lock}`, `cargo build --release -p andes --features thermo`, per memory's GOTCHA). Verify `andes` loads the own store and that at least one slug reports `gbdt_peak_model` present (add a one-line `andes --print-model <slug>` debug or confirm via a search log line).

- [ ] **Step 2: Run all four benchmarks through the uniform Percolator**

```bash
# on the VM (/srv/data/msgf-bench):
bash scripts/bench_astral_competitors.sh        # Astral high-res HCD
bash scripts/bench_tmt_ups_competitors.sh       # TMT a05058 + UPS1 low-res CID
bash scripts/bench_lysc_competitors.sh          # LysC non-tryptic showcase
```
Each emits andes (top-1 + `--chimeric`) + Java MS-GF+ + Sage + Comet + ProSE PINs → one Percolator 3.7.1 (`--seed 42 -Y`), counts at q ≤ 0.01.

- [ ] **Step 3: Run the entrapment-FDP honesty check with the own models**

```bash
bash scripts/astral_entrapment_experiment.sh    # 1:1 ENT_ database
python3 scripts/entrap_fdp.py --qcuts 0.005 0.01 0.02 0.05
# repeat the entrapment build for the LysC run to confirm ~1% true FDP non-tryptic
```
Acceptance for honesty: true FDP at nominal q≤1% stays ≤ ~1.2% on every dataset (Astral/TMT/UPS/LysC), and the chimeric gain holds at that FDP (per the spec's existing 1.14% result).

- [ ] **Step 4: Evaluate the gate**

Record PSMs@1% for andes (own models) vs the best open-source engine on each dataset. The gate PASSES iff ALL hold:
- andes ≥ best-of-field on PSMs on **Astral AND TMT AND UPS1** (all three), AND
- andes **strictly > field on Astral**, AND
- andes **beats field on LysC**, AND
- the entrapment-FDP honesty check above holds (no FDR violation).

Baseline from the spec (existing own models): already > field on Astral (+4.6%) and TMT (+3.7%); UPS1 top-1 is −1.5% (chimeric +6%). So A2's GBDT term must close the UPS1 top-1 gap and capture Astral headroom. The **honest A/B** that proves the GBDT term earned its place: run each dataset with the GBDT blob present vs a store with the blob stripped (same rank core), at matched entrapment-FDP — the GBDT must be net-positive on PSMs, gated on Astral (not TMT-only, per memory's astral-gates rule).

- [ ] **Step 5: Iterate (research loop) until the gate passes**

If the gate fails on a dataset, the lever is the model, not the engine: revisit label quality (co-isolation threshold, chimeric union, theoretical-ion enumeration breadth), the feature set (add a discriminator only after adding it to BOTH extractors + re-running A2.10), LightGBM depth/leaves/undersampling, or per-slug specialization (TMT, low-res CID are the largest levers per the spec). Re-train on Codon, re-transcode, re-assemble, redeploy, re-benchmark. Record each round's numbers. The spec explicitly flags this as multi-round research.

- [ ] **Step 6: Gate decision (no commit; this is a measurement + GO/NO-GO)**

When the gate passes, record the final table (own-models, all four datasets, + entrapment-FDP) — that table is the evidence Phase B's public benchmark doc and README are rewritten from. **Only then** does Phase B (the public-release PR) proceed.

---

## Testing & validation summary

- **Rust unit/integration (run after every Rust task):**
  `cargo test -p scoring && cargo test -p model-train && cargo test -p andes`
- **The correctness invariants (all Rust — no Python, so no cross-language gate):**
  1. **No-GBDT byte-identity.** The GBDT term must be byte-identical to the
     pre-GBDT engine when a model carries no blob (empty `gbdt_logit_by_rank` →
     `+0`). The committed golden `precursor_cal_off_pin_tsv_match_golden_after_sort`
     guards this — BUT it needs the 40 MB `test-fixtures/test.mgf`, ABSENT in this
     workspace (memory `reference-parity-fixture-gap`); a fixture-missing panic is
     NOT a pass. Verify instead with a before/after sorted PIN/TSV diff on
     `iprg-2013/F13.mgf` + `ecoli.fasta` (9808 rows) across a `git worktree` at the
     baseline commit — done for A2.5; re-run if a later task touches the no-GBDT path.
  2. **Walker round-trip + validate-at-decode** (`gbdt_eval` tests) and the
     **store round-trip** of `gbdt_model_bytes` (`model-train` roundtrip test).
  3. **Trainer contract** (`gbdt_trainer` test): determinism given a seed, model
     round-trip, signal>noise separation on synthetic data.
  4. **Once-per-peak caching:** the GBDT runs at spectrum-prep, never in the inner
     candidate loop (precomputed into `gbdt_logit_by_rank`).
- **Benchmark gate (VM):** A4.1 — the GO/NO-GO for the public release.
- **Known exception:** the pre-existing `integer_mass_scaler_matches_residue_table_mean`
  failure (memory) and the `bootstrap_*` model-train tests (missing
  `test-fixtures/test.mgf`) are NOT introduced by this work; do not "fix" them here.

## Risks (from the spec)

- A2 is genuine research: feature design, label quality (co-isolation), calibration, and beating the field will iterate across sessions on Codon+VM. The A4 gate may take several rounds.
- Label noise (real-but-unmatched ions labeled "noise") is the main modeling risk — mitigated by generous theoretical enumeration + co-isolation exclusion + robust tree settings.
- Own-slug coverage: the clean own-only store omits combos without an own model → loud WARN + `--model` (already implemented). Acceptable; expand in later releases.
