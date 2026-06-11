# Independent pooled-own prior — Implementation Plan (Phase 1, increment 1)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `train-from-msnet` shrink sparse-partition estimates toward an *independent, own-data* broad prior model instead of only the corpus-internal global pool, so thin/narrow corpora stop over-specializing (the Astral −8.3% regression) — without re-introducing any MS-GF+-derived content.

**Architecture:** Add an optional `prior: Option<&Param>` to the estimator. The prior is a broad "global own-model" (a previously-trained own Param sharing the template's partition layout). In the existing backoff, a new **Level 0** consults the prior's already-normalized per-`(partition, ion)` distribution as the blend parent before falling back to the current segment-collapse (Level 1) and global-pool (Level 2). The Bayesian blend `(n·emp + w·prior)/(n+w)` and all other machinery are reused unchanged. When `prior` is `None`, behavior is byte-identical to today.

**Tech Stack:** Rust (`crates/model-train`, `crates/andes`), `cargo test`. TDD.

---

### Task 1: Thread an optional prior into the rank-distribution backoff

**Files:**
- Modify: `crates/model-train/src/estimate.rs` (`estimate`, `build_rank_dist_table`)
- Test: `crates/model-train/tests/estimate.rs` (add one test)

- [ ] **Step 1: Write the failing test**

Append to `crates/model-train/tests/estimate.rs` (reuses the existing `one_partition_template` fixture and `CountStats` API already imported in that file):

```rust
/// A sparse partition (n < min_count) must shrink toward the INDEPENDENT PRIOR's
/// distribution, not the corpus-internal pool. Here the corpus empirical mass is
/// all on slot 0, but the prior is peaked on slot 5; the blended result must
/// carry materially more mass on slot 5 than the no-prior estimate does.
#[test]
fn sparse_partition_shrinks_toward_independent_prior() {
    let max_rank = 150;
    let n_slots = (max_rank + 1) as usize;
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let template = one_partition_template(max_rank);

    // Prior model: same layout, but its rank dist for (part, prefix1) peaks at slot 5.
    let mut prior = one_partition_template(max_rank);
    let mut prior_dist = vec![0.0_f32; n_slots];
    prior_dist[5] = 1.0;
    let mut ion_map = FxHashMap::default();
    ion_map.insert(prefix1, prior_dist);
    // Noise entry required so the prior is self-consistent for lookups.
    let mut noise_dist = vec![0.0_f32; n_slots];
    noise_dist[0] = 1.0;
    ion_map.insert(IonType::Noise, noise_dist);
    prior.rank_dist_table.insert(part, ion_map);

    // Sparse corpus: 10 observations (< default min_count 50), all on slot 0.
    let mut counts = CountStats::new();
    for _ in 0..10 {
        counts.bump_rank(part, prefix1, 0);
        counts.bump_rank(part, IonType::Noise, 0);
    }

    let est = Estimator::new(EstimatorConfig::default());
    let with_prior = est.estimate_with_prior(&counts, &template, Some(&prior));
    let no_prior = est.estimate_with_prior(&counts, &template, None);

    let p5_with = with_prior.rank_dist_table[&part][&prefix1][5];
    let p5_without = no_prior.rank_dist_table[&part][&prefix1][5];
    assert!(
        p5_with > p5_without + 0.05,
        "prior must pull mass toward slot 5: with={p5_with} without={p5_without}"
    );
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test -p model-train --test estimate sparse_partition_shrinks_toward_independent_prior`
Expected: FAIL to compile — `estimate_with_prior` does not exist yet.

- [ ] **Step 3: Add `estimate_with_prior` and thread `prior` into the rank table**

In `crates/model-train/src/estimate.rs`, replace the `estimate` method (currently starting at line 113) so the existing body moves into a prior-aware variant and `estimate` delegates with `None`:

```rust
    /// Build a [`Param`] from accumulated counts (no independent prior).
    pub fn estimate(&self, counts: &CountStats, template: &Param) -> Param {
        self.estimate_with_prior(counts, template, None)
    }

    /// Build a [`Param`] from counts, optionally shrinking sparse partitions
    /// toward an independent `prior` model (a broad own-data Param sharing the
    /// template's partition layout). When `prior` is `None` the result is
    /// identical to the prior-free estimate.
    pub fn estimate_with_prior(
        &self,
        counts: &CountStats,
        template: &Param,
        prior: Option<&Param>,
    ) -> Param {
```

(Keep the entire existing body of the old `estimate` below this new signature, unchanged except the four internal builder calls, which now forward `prior`:)

```rust
        let rank_dist_table = self.build_rank_dist_table(counts, template, max_rank, prior);
        let (ion_err_dist_table, noise_err_dist_table) =
            self.build_error_tables(counts, template, esf, prior);
        let ion_existence_table = self.build_existence_table(counts, template, prior);
```

Update `build_rank_dist_table`'s signature and its `parent_vec` closure. Change the signature (line 175) to add `prior: Option<&Param>`:

```rust
    fn build_rank_dist_table(
        &self,
        counts: &CountStats,
        template: &Param,
        max_rank: i32,
        prior: Option<&Param>,
    ) -> FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> {
```

Inside the `for &part in &all_partitions` loop, prepend a **Level 0** lookup to the `parent_vec` closure (the prior's distributions are already normalized probability vectors of length `n_slots`):

```rust
            let parent_vec = |ion: IonType, ps: f32| -> Vec<f32> {
                // Level 0: independent prior model (own-data broad prior).
                if let Some(pr) = prior {
                    if let Some(dist) = pr.rank_dist_table.get(&part).and_then(|m| m.get(&ion)) {
                        if dist.len() == n_slots {
                            return dist.clone();
                        }
                    }
                }
                // Level 1: segment-collapse.
                if let Some(seg_map) = seg_parent {
                    if let Some(raw) = seg_map.get(&ion) {
                        let n: u64 = raw.iter().sum();
                        if n >= min_count {
                            return normalize_with_pseudo(raw, n_slots, ps);
                        }
                    }
                }
                // Level 2: global pool.
                let graw = global_pool.get(&ion).map(|v| v.as_slice()).unwrap_or(&[]);
                normalize_with_pseudo(graw, n_slots, ps)
            };
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test -p model-train --test estimate sparse_partition_shrinks_toward_independent_prior`
Expected: PASS.

- [ ] **Step 5: Run the full estimate + incremental suites for no-regression (prior=None path)**

Run: `cargo test -p model-train --test estimate --test incremental --test yield_nonregression`
Expected: all PASS (the `estimate(counts, template)` delegation keeps existing behavior identical).

- [ ] **Step 6: Commit**

```bash
git add crates/model-train/src/estimate.rs crates/model-train/tests/estimate.rs
git commit -m "feat(model-train): shrink sparse rank partitions toward an independent prior model"
```

---

### Task 2: Extend the prior to error and existence tables

**Files:**
- Modify: `crates/model-train/src/estimate.rs` (`build_error_tables`, `build_existence_table`)
- Test: `crates/model-train/tests/estimate.rs` (add one test)

- [ ] **Step 1: Write the failing test**

Append to `crates/model-train/tests/estimate.rs`:

```rust
/// The independent prior must also drive the existence-table backoff for a
/// sparse partition. The prior puts existence mass on index 3; the (empty)
/// corpus would otherwise back off to a flat/global existence shape.
#[test]
fn sparse_existence_shrinks_toward_independent_prior() {
    let template = one_partition_template(150);
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };

    let mut prior = one_partition_template(150);
    prior.ion_existence_table.insert(part, vec![0.0, 0.0, 0.0, 1.0]); // mass on idx 3

    // No existence counts at all → n=0 < min_count → must use the prior.
    let counts = CountStats::new();
    let est = Estimator::new(EstimatorConfig::default());
    let with_prior = est.estimate_with_prior(&counts, &template, Some(&prior));

    let ex = &with_prior.ion_existence_table[&part];
    assert!(ex[3] > 0.5, "existence should follow the prior's idx-3 peak, got {ex:?}");
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test -p model-train --test estimate sparse_existence_shrinks_toward_independent_prior`
Expected: FAIL — `build_existence_table` does not accept a prior yet (compile error on the new call) OR assertion fails.

- [ ] **Step 3: Thread `prior` into error and existence tables**

Change `build_error_tables` (line 290) signature to add `prior: Option<&Param>`, and in its per-partition loop replace the parent arguments to `blend` with a prior-aware parent. After computing `global_ion_norm` / `global_noise_norm`, define helpers inside the `for &part in &template.partitions` loop:

```rust
            let ion_parent: Vec<f32> = prior
                .and_then(|p| p.ion_err_dist_table.get(&part))
                .filter(|d| d.len() == dist_len)
                .cloned()
                .unwrap_or_else(|| global_ion_norm.clone());
            let ion_dist = if ion_n < min_count {
                blend(&ion_emp, &ion_parent, ion_n as f32, w)
            } else {
                ion_emp
            };
```

and symmetrically for noise:

```rust
            let noise_parent: Vec<f32> = prior
                .and_then(|p| p.noise_err_dist_table.get(&part))
                .filter(|d| d.len() == dist_len)
                .cloned()
                .unwrap_or_else(|| global_noise_norm.clone());
            let noise_dist = if noise_n < min_count {
                blend(&noise_emp, &noise_parent, noise_n as f32, w)
            } else {
                noise_emp
            };
```

Change `build_existence_table` (line 356) signature to add `prior: Option<&Param>`, and in its per-partition loop replace the parent:

```rust
            let parent: Vec<f32> = prior
                .and_then(|p| p.ion_existence_table.get(&part))
                .filter(|d| d.len() == N_EX)
                .cloned()
                .unwrap_or_else(|| global_norm.clone());
            let dist = if n < min_count {
                blend(&emp, &parent, n as f32, w)
            } else {
                emp
            };
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test -p model-train --test estimate sparse_existence_shrinks_toward_independent_prior`
Expected: PASS.

- [ ] **Step 5: Full no-regression run**

Run: `cargo test -p model-train`
Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/model-train/src/estimate.rs crates/model-train/tests/estimate.rs
git commit -m "feat(model-train): extend independent-prior backoff to error and existence tables"
```

---

### Task 3: CLI — supply the prior to `train-from-msnet`

**Files:**
- Modify: `crates/andes/src/bin/andes.rs` (`TrainFromMsnetArgs` struct, `run_train_from_msnet`)
- Test: `crates/andes/tests/train_from_msnet.rs` (add one test)

- [ ] **Step 1: Write the failing test**

Add to `crates/andes/tests/train_from_msnet.rs` a test that invokes the binary with `--prior-model-store <store> --prior-model <id>` on the existing fixtures and asserts the run succeeds and writes a model. (Follow the existing invocation pattern already in that file — reuse its fixture store and `assert_cmd`/output-parsing helpers; the new flags are additive and default to no prior.)

```rust
#[test]
fn train_from_msnet_accepts_prior_model_flags() {
    // Arrange: reuse the same fixtures the other tests in this file use for
    // --in / --out-store / --seed-model; pass the seed store ALSO as the prior
    // store and the same slug as --prior-model (mechanism smoke test: a valid
    // prior is loaded and the run still produces a model).
    // (Use the existing helper that builds the base command in this test file.)
    let out = run_train_from_msnet_cmd(&[
        "--prior-model-store", FIXTURE_STORE,
        "--prior-model", FIXTURE_SLUG,
    ]);
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
    assert!(model_written_to_out_store(), "expected a trained model in the out store");
}
```

> Note for the implementer: open `crates/andes/tests/train_from_msnet.rs` first and reuse its existing constants/helpers (`FIXTURE_STORE`, `FIXTURE_SLUG`, the command builder, and the out-store assertion) rather than inventing new ones — they already exist for the sibling tests.

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test -p andes --test train_from_msnet train_from_msnet_accepts_prior_model_flags`
Expected: FAIL — the `--prior-model-store` / `--prior-model` args are unknown.

- [ ] **Step 3: Add the flags and load + pass the prior**

In `TrainFromMsnetArgs` (the clap args struct for the subcommand), add:

```rust
    /// Optional path to an independent prior model store. Sparse partitions in
    /// the trained model shrink toward the matching prior model instead of the
    /// corpus-internal pool. Must be own-data (NOT the MS-GF+ seed) to stay
    /// relicense-safe.
    #[arg(long)]
    prior_model_store: Option<std::path::PathBuf>,

    /// Model id to load from `--prior-model-store` (defaults to the trained
    /// model id when omitted).
    #[arg(long)]
    prior_model: Option<String>,
```

In `run_train_from_msnet`, just before the `estimator.estimate(&stats, &seed_param)` call (line 2511), load the optional prior and switch to `estimate_with_prior`:

```rust
    let prior_param: Option<Param> = match &args.prior_model_store {
        Some(store_path) => {
            let prior_id = args.prior_model.clone().unwrap_or_else(|| model_id.clone());
            let (_pid, p) = load_param_from_store(store_path, &prior_id)
                .map_err(|e| format!("loading --prior-model '{prior_id}': {e}"))?;
            Some(p)
        }
        None => None,
    };
    let mut trained_param =
        estimator.estimate_with_prior(&stats, &seed_param, prior_param.as_ref());
```

(Where `model_id` is the trained model's id already in scope in this function; confirm its binding name when editing and reuse it.)

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test -p andes --test train_from_msnet train_from_msnet_accepts_prior_model_flags`
Expected: PASS.

- [ ] **Step 5: Build + clippy**

Run: `cargo build -p andes && cargo clippy -p andes -p model-train -- -D warnings`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add crates/andes/src/bin/andes.rs crates/andes/tests/train_from_msnet.rs
git commit -m "feat(andes): --prior-model-store/--prior-model to shrink train-from-msnet toward an own-data prior"
```

---

### Task 4: Empirical validation — build a pooled-own prior and recover the Astral regression (operational, on Codon + VM)

This is not a unit test; it is the acceptance check for the whole increment. Run it after Tasks 1–3 land and a fresh `andes-bin` is built on Codon from this branch.

- [ ] **Step 1: Build a broad pooled-own prior on Codon.** Aggregate the harvested own flats across many high-res slugs into one broad corpus and train a single "global own" model on the `hcd_qexactive_tryp` structural template (reuse `CountStats::add` semantics / the `--update` multi-source path, or pass many `--in` flats). Write it to `stores/_priors/global_highres.parquet` with model id `hcd_qexactive_tryp`.

- [ ] **Step 2: Retrain `hcd_qexactive_tryp` with the prior.** Re-run the slug's `train-from-msnet` with `--prior-model-store stores/_priors/global_highres.parquet --prior-model hcd_qexactive_tryp` into a fresh store. Keep all other flags identical to the un-prior run.

- [ ] **Step 3: Re-benchmark Astral (and UPS1/TMT as guards).** Assemble the prior-trained `hcd_qexactive_tryp` into a combined own store, ship to the VM, run the andes-own Astral search, percolate at 1% FDR via `run_percolator_docker.sh`.

- [ ] **Step 4: Check the acceptance criterion.** PASS if Astral own ≥ 0.99 × seed (recovers the −8.3% toward parity) AND UPS1/TMT hold their gains. Record the numbers in the campaign ledger. If it under-recovers, that is the signal to proceed to Phase-1 increment 2 (re-derived structure + sparse-aware smoothing), not to widen the prior weight blindly.

---

## Self-review notes
- **Spec coverage:** implements roadmap Phase-1 item 1 (independent pooled-own prior). Items 2–5 (re-derive structure, sparse-aware smoothing, noise fix, wire acceptance gate) are separate plans.
- **Back-compat:** `estimate(counts, template)` is preserved verbatim via delegation to `estimate_with_prior(..., None)`; all existing call sites and tests are unaffected. The `prior=None` path is byte-identical.
- **Relicense safety:** the prior MUST be own-data; the CLI doc string and Task 4 Step 1 enforce building it from own flats, never the MS-GF+ seed. (Blending toward the seed would defeat independence — that is exactly why this increment uses a pooled-own prior, not `--seed-model`.)
- **Type consistency:** `estimate_with_prior` / `build_rank_dist_table` / `build_error_tables` / `build_existence_table` all take `prior: Option<&Param>`; `Param.rank_dist_table` / `ion_err_dist_table` / `noise_err_dist_table` / `ion_existence_table` are the exact field names from `scoring_crate::param_model::Param`.
