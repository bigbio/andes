# MS-GF+-style rank-window smoothing — Implementation Plan (Phase 1, increment 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make own-trained scoring models *generalize* like the MS-GF+ seed by re-introducing MS-GF+'s widening rank-window smoothing of the signal rank distributions — a clean-room re-implementation of a published statistical method (PepNovo/MS-GF+), trained on our own data (Apache-safe). This addresses the diagnosed defect: our dense empirical models are *over-sharp* on b/y signal (peakier than the seed), which over-fits training spectra and transfers poorly to other high-res data (Q Exactive→Astral).

**Architecture:** Add an opt-in `rank_smoothing` flag to `EstimatorConfig`. After each **signal** ion's rank distribution is built (post-blend, post-normalize) in `build_rank_dist_table`, apply a moving-average whose half-width *widens with rank* — leaving the discriminative head (ranks 1–3) and the missing-ion sentinel untouched, smoothing only the noisy tail — then renormalize. The **noise** distribution is NOT touched here (its over-softening is a separate increment). Default `false` keeps every existing model byte-identical.

**Tech Stack:** Rust (`crates/model-train`), `cargo test`. TDD.

**Why this is Apache-safe:** rank-window smoothing is a published statistical technique (Frank & Pevzner 2005; Kim & Pevzner 2014), not the patented generating function (already removed) and not MS-GF+'s code. We implement it from the method description and train on own data. Same clean-room basis as the existing Phase-2 scoring reauthoring.

---

### Task 1: `smooth_rank_window` helper + unit test

**Files:**
- Modify: `crates/model-train/src/estimate.rs` (add a private helper)
- Test: `crates/model-train/tests/estimate.rs` (add one test)

- [ ] **Step 1: Write the failing test** — append to `crates/model-train/tests/estimate.rs`. It calls the new helper (re-export it for the test via `pub(crate)` + a thin test shim, OR test it through `estimate_with_prior` once Task 2 wires it; here we test the function directly by making it `pub` in the module and importing it). Add:

```rust
use model_train::estimate::smooth_rank_window; // make pub in Task 1 Step 3

/// Widening rank-window smoothing must (a) leave ranks 1-3 (indices 0..3) and the
/// missing-ion sentinel (last index) untouched, (b) smooth the tail (reduce a lone
/// spike at a high rank by averaging neighbors), and (c) renormalize to sum 1.
#[test]
fn rank_window_smoothing_preserves_head_smooths_tail() {
    let max_rank = 150usize;
    let n = max_rank + 1; // observed ranks 0..150 + missing slot at index 150
    // Distribution: sharp mass at index 0 (rank 1), a lone spike at index 40, missing slot.
    let mut d = vec![0.0f32; n];
    d[0] = 0.50;     // discriminative head — must be preserved
    d[40] = 0.40;    // tail spike — must be spread to neighbors
    d[150] = 0.10;   // missing-ion sentinel — must be preserved
    let out = smooth_rank_window(&d, max_rank);

    // sums to 1
    let s: f32 = out.iter().sum();
    assert!((s - 1.0).abs() < 1e-4, "must renormalize, sum={s}");
    // head preserved (relative to its own pre-norm value; index 0 has window 0)
    assert!(out[0] > 0.45, "rank-1 head must stay sharp, got {}", out[0]);
    // missing-ion sentinel preserved
    assert!(out[150] > 0.08, "missing slot must be preserved, got {}", out[150]);
    // tail spike spread: index 40 reduced, neighbors raised
    assert!(out[40] < 0.40, "tail spike must be smoothed down, got {}", out[40]);
    assert!(out[39] > 0.0 && out[41] > 0.0, "tail neighbors must receive mass");
}
```

- [ ] **Step 2: Run, confirm it FAILS** (function missing): `cargo test -p model-train --test estimate rank_window_smoothing_preserves_head_smooths_tail`

- [ ] **Step 3: Implement** the helper in `crates/model-train/src/estimate.rs` (place it near `normalize_with_pseudo`, make it `pub` so the integration test can import it):

```rust
/// MS-GF+-style widening rank-window smoothing of a rank-distribution vector.
///
/// `dist` has length `max_rank + 1`: indices `0..max_rank` are observed ranks
/// 1..max_rank, index `max_rank` is the missing-ion sentinel. The half-width of
/// the moving average widens with rank — the discriminative head (ranks 1-3) and
/// the missing-ion sentinel are left untouched; the noisy tail is averaged with
/// neighbors. The result is renormalized to sum 1.0.
///
/// Window schedule (from the published method, smoothingRanks {3,5,10,20,50}):
/// rank<3→hw 0 (no smoothing), <5→1, <10→2, <20→3, <50→4, else→5.
pub fn smooth_rank_window(dist: &[f32], max_rank: usize) -> Vec<f32> {
    let n = dist.len();
    let last = max_rank.min(n.saturating_sub(1)); // missing-ion sentinel index
    let mut out = dist.to_vec();
    let halfwidth = |r: usize| -> usize {
        if r < 3 { 0 } else if r < 5 { 1 } else if r < 10 { 2 }
        else if r < 20 { 3 } else if r < 50 { 4 } else { 5 }
    };
    // Smooth only observed-rank slots [0..last); never the missing-ion sentinel.
    for i in 0..last {
        let hw = halfwidth(i);
        if hw == 0 { continue; }
        let lo = i.saturating_sub(hw);
        let hi = (i + hw + 1).min(last); // exclusive; excludes sentinel
        let mut s = 0.0f32;
        let mut c = 0usize;
        for v in dist.iter().take(hi).skip(lo) { s += *v; c += 1; }
        if c > 0 { out[i] = s / c as f32; }
    }
    let tot: f32 = out.iter().sum();
    if tot > 0.0 {
        for x in &mut out { *x /= tot; }
    }
    out
}
```

- [ ] **Step 4: Run, confirm PASS**: `cargo test -p model-train --test estimate rank_window_smoothing_preserves_head_smooths_tail`

- [ ] **Step 5: Commit**
```bash
git add crates/model-train/src/estimate.rs crates/model-train/tests/estimate.rs
git commit -m "feat(model-train): add MS-GF+-style widening rank-window smoothing helper"
```

---

### Task 2: wire `rank_smoothing` into the estimator (signal ions only, opt-in)

**Files:**
- Modify: `crates/model-train/src/estimate.rs` (`EstimatorConfig`, `build_rank_dist_table`)
- Test: `crates/model-train/tests/estimate.rs` (add one test)

- [ ] **Step 1: Write the failing test** — append:

```rust
/// With rank_smoothing enabled, a signal ion's trained rank distribution is
/// smoother (lower peak) than with it disabled; the NOISE distribution is
/// unchanged (smoothing must not touch noise). Uses a dense single partition.
#[test]
fn rank_smoothing_softens_signal_not_noise() {
    let max_rank = 150;
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let template = one_partition_template(max_rank);

    // Dense counts: signal concentrated at a mid rank (so smoothing has an effect),
    // noise concentrated at slot 0.
    let mut counts = CountStats::new();
    for _ in 0..2000 {
        counts.bump_rank(part, prefix1, 25); // mid-tail signal
        counts.bump_rank(part, IonType::Noise, 0);
    }

    let mut cfg_on = EstimatorConfig::default();
    cfg_on.rank_smoothing = true;
    let on = Estimator::new(cfg_on).estimate(&counts, &template);
    let off = Estimator::new(EstimatorConfig::default()).estimate(&counts, &template);

    let sig_on = on.rank_dist_table[&part][&prefix1][25];
    let sig_off = off.rank_dist_table[&part][&prefix1][25];
    assert!(sig_on < sig_off, "smoothing must lower the signal peak: on={sig_on} off={sig_off}");

    // Noise distribution must be identical (smoothing is signal-only).
    let noise_on = &on.rank_dist_table[&part][&IonType::Noise];
    let noise_off = &off.rank_dist_table[&part][&IonType::Noise];
    assert_eq!(noise_on, noise_off, "noise dist must be unchanged by rank smoothing");
}
```

- [ ] **Step 2: Run, confirm it FAILS** (no `rank_smoothing` field): `cargo test -p model-train --test estimate rank_smoothing_softens_signal_not_noise`

- [ ] **Step 3: Implement**
  (a) Add to `EstimatorConfig` (with a doc comment) and to its `Default`:
```rust
    /// Apply MS-GF+-style widening rank-window smoothing to SIGNAL rank
    /// distributions (not noise). Default: `false` (byte-identical to before).
    pub rank_smoothing: bool,
```
```rust
            rank_smoothing: false,
```
  (b) In `build_rank_dist_table`, capture the flag (`let rank_smoothing = self.cfg.rank_smoothing;`) and, where each SIGNAL ion's `blended` vector is inserted into `ion_table`, apply smoothing first:
```rust
                let final_dist = if rank_smoothing {
                    smooth_rank_window(&blended, n_slots - 1)
                } else {
                    blended
                };
                ion_table.insert(ion, final_dist);
```
  Leave the `IonType::Noise` insertion path UNCHANGED (do not smooth noise). (`n_slots - 1 == max_rank`, the missing-ion sentinel index.)

- [ ] **Step 4: Run, confirm PASS**: `cargo test -p model-train --test estimate rank_smoothing_softens_signal_not_noise`

- [ ] **Step 5: No-regression**: `cargo test -p model-train` (the `rank_smoothing=false` default keeps all existing tests green).

- [ ] **Step 6: Commit**
```bash
git add crates/model-train/src/estimate.rs crates/model-train/tests/estimate.rs
git commit -m "feat(model-train): opt-in rank_smoothing wires widening smoothing into signal rank dists"
```

---

### Task 3: CLI `--rank-smoothing` for `train-from-msnet`

**Files:**
- Modify: `crates/andes/src/bin/andes.rs` (`TrainFromMsnetArgs`, `run_train_from_msnet` EstimatorConfig construction)
- Test: `crates/andes/tests/train_from_msnet.rs` (add one test)

- [ ] **Step 1: Write a failing integration test** that runs `train-from-msnet --rank-smoothing` on the existing fixtures and asserts success + a model is written (reuse the file's existing command/fixture helpers — open the test file first and reuse its constants like the sibling tests do).

- [ ] **Step 2: Run, confirm FAIL** (unknown arg): `cargo test -p andes --test train_from_msnet <name>`

- [ ] **Step 3: Implement** — add the flag to `TrainFromMsnetArgs`:
```rust
    /// Apply MS-GF+-style rank-window smoothing to signal rank distributions
    /// (improves generalization of own-trained high-res models).
    #[arg(long)]
    rank_smoothing: bool,
```
and in `run_train_from_msnet`, set it on the `EstimatorConfig` that is built before `Estimator::new(cfg)` (the struct literal around the existing `EstimatorConfig { ... }`): add `rank_smoothing: args.rank_smoothing,`.

- [ ] **Step 4: Run, confirm PASS**; **Step 5:** `cargo build -p andes && cargo clippy -p model-train -p andes -- -D warnings`.

- [ ] **Step 6: Commit**
```bash
git add crates/andes/src/bin/andes.rs crates/andes/tests/train_from_msnet.rs
git commit -m "feat(andes): --rank-smoothing flag for train-from-msnet"
```

---

### Task 4: Empirical re-measure (operational, Codon + VM)

After Tasks 1–3 land and `andes-bin-prior` is rebuilt on Codon from this branch:
- [ ] Retrain `hcd_qexactive_tryp` on the broadened corpus **with `--rank-smoothing`** into `stores/hcd_qexactive_tryp_smooth/`.
- [ ] Ship to VM; benchmark BOTH the **held-out QE dataset** (from the finder agent) and **Astral** vs the seed and vs the un-smoothed own model.
- [ ] Acceptance signal: smoothed own ≥ un-smoothed own on QE AND on Astral, and ideally closes toward the seed (36,243 Astral). Record numbers + verdict in `memory/project_independence_license_status.md`. If smoothing helps Astral *without* Astral data, it confirms the defect was model-construction, not data.

---

## Self-review
- **Back-compat:** `rank_smoothing` defaults `false`; the `estimate`/`estimate_with_prior` paths are byte-identical when off. Noise distribution never smoothed.
- **Apache-safe:** published method, clean-room, own-data training (see header).
- **Scope:** signal-ion rank smoothing only. The noise-model over-softening fix and any adaptive-binning/min-spectra-floor work are separate increments (the min-spectra floor is low-value for our *dense* models — 0% sparse cells — so it is deliberately NOT in this plan).
- **Type consistency:** `smooth_rank_window(dist: &[f32], max_rank: usize) -> Vec<f32>`; `EstimatorConfig.rank_smoothing: bool`; smoothing applied to signal ions in `build_rank_dist_table`, never to `IonType::Noise`.
