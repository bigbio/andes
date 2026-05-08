# Rayon-Parallel match_spectra Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Parallelize the per-spectrum search loop in `match_spectra` over a Rayon worker pool, mirroring Java MS-GF+'s `-thread N` execution model. Targets the dominant single-thread bottleneck identified in profiling (PXD001819 200-spectrum slice: GF 50.7%, full-dataset projection ~14× per-thread vs Java).

**Architecture:** Convert the outer `for (spec_idx, spec) in spectra.iter().enumerate()` loop in `search/src/match_engine.rs::match_spectra` into a `par_iter().enumerate().map(...).collect()` chain. All current inputs are `&` immutable after the upfront `enumerate_candidates` and `aa_set_for_gf` setup, so no shared `&mut` state is introduced. Each Rayon worker owns its `TopNQueue`, `scored_per_charge`, and per-bin GF state for the duration of its closure. CLI gains a `--threads N` flag (default `num_cpus::get()`) configuring a global Rayon pool.

**Tech Stack:** Rust 2021, Rayon 1.10, num_cpus 1.16. Existing search/scoring/output crates unchanged structurally.

**Spec:** [docs/superpowers/specs/2026-05-08-rayon-match-spectra-design.md](../specs/2026-05-08-rayon-match-spectra-design.md)

---

## File map

| File | Change |
|---|---|
| `rust/crates/search/Cargo.toml` | add `rayon = "1.10"` dependency |
| `rust/crates/search/src/match_engine.rs` | convert outer spectrum loop to `par_iter`; closure body returns `TopNQueue` per index |
| `rust/crates/search/tests/match_spectra_thread_invariance.rs` | new — assert `--threads 1` and `--threads 4` produce bit-identical PSM identity + spec_e_value |
| `rust/crates/msgf-rust/Cargo.toml` | add `num_cpus = "1.16"` and `rayon = "1.10"` |
| `rust/crates/msgf-rust/src/bin/msgf-rust.rs` | add `--threads N` CLI flag; configure global Rayon pool before `match_spectra` |
| `benchmark/parity/run_pxd001819_2arm.sh` | add `--threads 12` to RUST_ARGS to mirror Java's default 12 threads |
| `~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_rust_perf_gap_vs_java.md` | append post-Rayon measurement |

---

## Task 1: Add `rayon` dependency to the search crate

**Files:**
- Modify: `rust/crates/search/Cargo.toml`

- [ ] **Step 1.1: Add the dependency**

In the `[dependencies]` section of `rust/crates/search/Cargo.toml`, add:

```toml
rayon = "1.10"
```

The exact placement (alphabetical):

```toml
[dependencies]
model = { path = "../model" }
rayon = "1.10"
scoring_crate = { path = "../scoring", package = "scoring" }
suffix = { workspace = true }
thiserror = { workspace = true }
```

- [ ] **Step 1.2: Verify the dep resolves**

Run:
```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && cargo build -p search 2>&1 | tail -5
```

Expected: clean build (no behavior change yet).

- [ ] **Step 1.3: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add rust/crates/search/Cargo.toml rust/Cargo.lock
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "chore(search): add rayon dependency"
```

---

## Task 2: Convert the outer spectrum loop to `par_iter`

**Files:**
- Modify: `rust/crates/search/src/match_engine.rs:66-194` (the `for (spec_idx, spec) in spectra.iter().enumerate()` block in `match_spectra`)

This is the core change. The loop body (~130 lines) becomes a closure that returns a `TopNQueue` per spectrum, and `collect()` builds the result `Vec<TopNQueue>`.

- [ ] **Step 2.1: Add the rayon import**

At the top of `rust/crates/search/src/match_engine.rs`, add:

```rust
use rayon::prelude::*;
```

Place it among the existing crate imports (alphabetical or grouped — match the file's existing style).

- [ ] **Step 2.2: Replace the loop with `par_iter`**

Find the current loop (lines ~40-194):

```rust
    let mut queues: Vec<TopNQueue> = (0..spectra.len())
        .map(|_| TopNQueue::new(params.top_n_psms_per_spectrum))
        .collect();

    let candidates: Vec<Candidate> = enumerate_candidates(idx, params, decoy_prefix).collect();
    // ... bucket_index + aa_set_for_gf setup unchanged ...

    for (spec_idx, spec) in spectra.iter().enumerate() {
        // ~130 lines of per-spectrum work, mutating queues[spec_idx]
    }

    queues
}
```

Replace the upfront `let mut queues` allocation and the `for` loop with a single `par_iter` chain. The new structure:

```rust
    let candidates: Vec<Candidate> = enumerate_candidates(idx, params, decoy_prefix).collect();

    // Build mass-bucket index: nominal(peptide.mass() - H2O) → Vec<candidate_idx>.
    // (existing comment block preserved)
    let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (cand_idx, cand) in candidates.iter().enumerate() {
        let nominal = nominal_from(cand.peptide.mass() - H2O);
        bucket_index.entry(nominal).or_default().push(cand_idx);
    }

    // Build an aa_set clone with enzyme registered (for GF cleavage scoring).
    // (existing comment block preserved)
    let mut aa_set_for_gf: AminoAcidSet = params.aa_set.clone();
    if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
        aa_set_for_gf.register_enzyme(params.enzyme, 0.95, 0.95);
    }

    // Parallel per-spectrum work. Each Rayon worker owns its TopNQueue,
    // scored_per_charge cache, and per-bin GF state for the closure's duration.
    // All inputs (spectra, idx, params, scorer, candidates, bucket_index,
    // aa_set_for_gf, decoy_prefix) are `&` immutable, so no &mut shared state.
    let queues: Vec<TopNQueue> = spectra
        .par_iter()
        .enumerate()
        .map(|(spec_idx, spec)| {
            let mut queue = TopNQueue::new(params.top_n_psms_per_spectrum);

            // Skip spectra with too few peaks (mirrors Java's `-minNumPeaks` filter).
            if spec.peaks.len() < params.min_peaks as usize {
                return queue;
            }

            // Determine which charge states to try for this spectrum.
            let charges_to_try: Vec<u8> = match spec.precursor_charge {
                Some(z) if z > 0 => vec![z as u8],
                _ => params.charge_range.clone().collect(),
            };

            // Build (and cache) a ScoredSpectrum per charge to evaluate.
            let mut scored_per_charge: HashMap<u8, ScoredSpectrum<'_>> = HashMap::new();
            for &z in &charges_to_try {
                scored_per_charge.entry(z)
                    .or_insert_with(|| ScoredSpectrum::new(spec, scorer.param(), z));
            }

            // Compute per-charge candidate windows and union them.
            let mut window_cand_indices: HashSet<usize> = HashSet::new();
            for &z in &charges_to_try {
                let charge_f = z as f64;
                let neutral_mass = (spec.precursor_mz - PROTON) * charge_f - H2O;
                let nominal_center = nominal_from(neutral_mass);
                let iso_min = *params.isotope_error_range.start() as i32;
                let iso_max = *params.isotope_error_range.end() as i32;
                let tol_da_left  = params.precursor_tolerance.left.as_da(neutral_mass);
                let tol_da_right = params.precursor_tolerance.right.as_da(neutral_mass);
                let widen_left  = (tol_da_left  - 0.4999_f64).round() as i32;
                let widen_right = (tol_da_right - 0.4999_f64).round() as i32;
                let min_nominal = nominal_center - iso_max - widen_right;
                let max_nominal = nominal_center - iso_min + widen_left;
                for (_nm, idxs) in bucket_index.range(min_nominal..=max_nominal) {
                    for &ci in idxs {
                        window_cand_indices.insert(ci);
                    }
                }
            }

            // Per-candidate scoring loop — unchanged, but pushes to the local
            // `queue` (was `queues[spec_idx]` in the serial version).
            for &cand_idx in &window_cand_indices {
                let cand = &candidates[cand_idx];
                for &z in &charges_to_try {
                    let scored_spec = &scored_per_charge[&z];
                    let mut best_for_charge: Option<(MassError, f32)> = None;
                    for offset in params.isotope_error_range.clone() {
                        if let Some(err) = matches_precursor(spec, &cand.peptide, z, offset, &params.precursor_tolerance) {
                            let score = score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da);
                            if best_for_charge.as_ref().map_or(true, |(_, s)| score > *s) {
                                best_for_charge = Some((err, score));
                            }
                        }
                    }
                    if let Some((err, score)) = best_for_charge {
                        let features = compute_psm_features(scored_spec, &cand.peptide, fragment_tolerance_da);
                        queue.push(PsmMatch {
                            spectrum_idx: spec_idx,
                            candidate: cand.clone(),
                            charge_used: z,
                            mass_error_ppm: err.mass_error_ppm,
                            score,
                            spec_e_value: 1.0,
                            de_novo_score: i32::MIN,
                            activation_method: Some(scorer.param().data_type.activation),
                            e_value: 1.0,
                            features,
                            isotope_offset: err.isotope_offset,
                        });
                    }
                }
            }

            // Phase 6: compute SpecEValue for the PSMs in this queue.
            if !queue.is_empty() {
                let enzyme_opt = if params.enzyme != Enzyme::NoCleavage
                    && params.enzyme != Enzyme::NonSpecific
                {
                    Some(params.enzyme)
                } else {
                    None
                };
                let top_charge = queue
                    .iter_psms()
                    .max_by(|a, b| a.cmp(b))
                    .map(|p| p.charge_used)
                    .unwrap_or(charges_to_try[0]);
                let scored_spec_for_gf = &scored_per_charge[&top_charge];
                compute_spec_e_values_for_spectrum(
                    spec,
                    params,
                    &mut queue,
                    &aa_set_for_gf,
                    enzyme_opt,
                    scorer,
                    scored_spec_for_gf,
                    top_charge,
                    fragment_tolerance_da,
                    idx,
                );
            }

            queue
        })
        .collect();

    queues
}
```

The body of the closure is the EXACT same code as the existing for-loop body, with two mechanical substitutions:
- `queues[spec_idx].push(...)` → `queue.push(...)`
- `queues[spec_idx].is_empty()` and `&mut queues[spec_idx]` → `queue.is_empty()` and `&mut queue`
- `continue` → `return queue` (early exit on `min_peaks` filter)

The `let mut queues: Vec<TopNQueue> = (0..spectra.len()).map(...)` upfront allocation is removed — `collect()` builds the Vec in `spec_idx` order automatically thanks to `par_iter().enumerate()`.

- [ ] **Step 2.3: Build + verify**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && cargo build --workspace 2>&1 | tail -10
```

Expected: clean build. If the borrow checker complains about `params.charge_range.clone().collect()` inside the closure (because `params` is captured by reference but `clone` is called twice), the fix is to bind `let charge_range = params.charge_range.clone();` before the loops.

- [ ] **Step 2.4: Run all existing tests under default thread count**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test --workspace --lib 2>&1 | grep "test result:" | tail -10 && \
  cargo test --release -p msgf-rust --test cli_smoke 2>&1 | tail -5 && \
  cargo test --release -p search --test match_engine_java_parity rust_top1_matches_java_top1_for_majority_of_spectra -- --nocapture 2>&1 | tail -5
```

Expected:
- All workspace lib tests pass (~294)
- `cli_smoke` BSA test still produces 435/5760 PSMs
- Top-1 identity rate still 214/217 = 98.6%

If any test fails, the parallelization broke something. Most likely cause: a closure variable captured the wrong way (by-move when by-ref was intended). Inspect the compiler error.

- [ ] **Step 2.5: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add rust/crates/search/src/match_engine.rs
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "perf(search): parallelize match_spectra outer loop with rayon par_iter"
```

---

## Task 3: Add `--threads N` CLI flag

**Files:**
- Modify: `rust/crates/msgf-rust/Cargo.toml`
- Modify: `rust/crates/msgf-rust/src/bin/msgf-rust.rs`

- [ ] **Step 3.1: Add deps**

In `rust/crates/msgf-rust/Cargo.toml` `[dependencies]` (alphabetical):

```toml
num_cpus = "1.16"
rayon = "1.10"
```

- [ ] **Step 3.2: Add the CLI field**

In `rust/crates/msgf-rust/src/bin/msgf-rust.rs`, find the `Cli` struct (it has `#[derive(Parser)]`). Add:

```rust
    /// Number of worker threads for the search loop. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,
```

Place it after the existing search-param flags (e.g., after `top_n` or `ntt`).

- [ ] **Step 3.3: Configure the Rayon pool before `match_spectra`**

In the `run()` function (or wherever `match_spectra` is called), insert BEFORE the `match_spectra(...)` call:

```rust
    // Configure the global Rayon worker pool. Mirrors Java's -thread N flag.
    // build_global() panics if called twice, so guard with a OnceLock so
    // multiple invocations (e.g., in tests) don't blow up.
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .expect("build_global");
    });
    eprintln!("Using {} worker threads", cli.threads);
```

Note: `OnceLock::get_or_init` cannot capture a per-call `cli.threads` value if called twice with different values — it's keyed on first call. For the CLI, this is fine (one call per binary invocation). For tests, see Task 4.

- [ ] **Step 3.4: Build + manual smoke**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && cargo build --release -p msgf-rust 2>&1 | tail -5
```

Then quickly test the new flag parses:

```bash
./target/release/msgf-rust --help 2>&1 | grep -A1 threads
```

Expected: `--threads <THREADS>` line with the description.

- [ ] **Step 3.5: Re-run cli_smoke (default thread count)**

```bash
cargo test --release -p msgf-rust --test cli_smoke 2>&1 | tail -5
```

Expected: still passes; output still 435/5760 PSMs.

- [ ] **Step 3.6: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add \
  rust/crates/msgf-rust/Cargo.toml \
  rust/crates/msgf-rust/src/bin/msgf-rust.rs \
  rust/Cargo.lock
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "feat(msgf-rust): --threads N CLI flag configures rayon pool"
```

---

## Task 4: Thread-count invariance test

**Files:**
- Create: `rust/crates/search/tests/match_spectra_thread_invariance.rs`

The test runs `match_spectra` twice on a small fixture under different thread counts, asserting per-spectrum PSM identity and spec_e_value are bit-identical.

- [ ] **Step 4.1: Write the failing test (TDD)**

Create `rust/crates/search/tests/match_spectra_thread_invariance.rs`:

```rust
//! Thread-count invariance: match_spectra must produce bit-identical output
//! regardless of the Rayon thread count, because each spectrum's full pipeline
//! runs entirely on one worker (no FP-accumulation non-determinism across
//! threads).

mod common;
use common::*;

use std::fs::File;
use std::io::BufReader;
use input::{FastaReader, MgfReader};
use search::{
    match_spectra, Enzyme, PrecursorTolerance, SearchIndex, SearchParams,
};
use model::Tolerance;

fn run_search(thread_count: usize) -> Vec<search::TopNQueue> {
    // Build a Rayon pool scoped to this thread count using `install`,
    // which doesn't conflict with global pool initialization.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .expect("build pool");

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("src/test/resources/BSA.fasta")).unwrap()
    )).unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX_");
    let aa = aa_set();
    let scorer = rank_scorer();
    let mut params = SearchParams::default_tryptic(aa.clone());
    params.enzyme = Enzyme::Trypsin;
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    params.charge_range = 2..=3;
    params.isotope_error_range = -1..=2;

    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    pool.install(|| {
        match_spectra(&spectra, &idx, &params, &scorer, 0.5, "XXX_")
    })
}

#[test]
fn match_spectra_output_invariant_across_thread_counts() {
    let q1 = run_search(1);
    let q4 = run_search(4);

    assert_eq!(q1.len(), q4.len(), "queue count differs");

    let mut spectra_with_psms = 0;
    for (i, (qa, qb)) in q1.iter().zip(q4.iter()).enumerate() {
        let psms_a = qa.clone().into_sorted_vec();
        let psms_b = qb.clone().into_sorted_vec();
        assert_eq!(
            psms_a.len(), psms_b.len(),
            "spectrum {}: PSM count differs ({} vs {})",
            i, psms_a.len(), psms_b.len()
        );
        if !psms_a.is_empty() {
            spectra_with_psms += 1;
            for (j, (a, b)) in psms_a.iter().zip(psms_b.iter()).enumerate() {
                let pep_a = a.candidate.peptide.residues.iter()
                    .map(|aa| aa.residue as char).collect::<String>();
                let pep_b = b.candidate.peptide.residues.iter()
                    .map(|aa| aa.residue as char).collect::<String>();
                assert_eq!(pep_a, pep_b,
                    "spectrum {} PSM rank {}: peptide differs ({} vs {})",
                    i, j, pep_a, pep_b);
                assert_eq!(a.charge_used, b.charge_used,
                    "spectrum {} PSM rank {}: charge differs", i, j);
                assert_eq!(a.score.to_bits(), b.score.to_bits(),
                    "spectrum {} PSM rank {}: score differs ({} vs {})",
                    i, j, a.score, b.score);
                assert_eq!(a.spec_e_value.to_bits(), b.spec_e_value.to_bits(),
                    "spectrum {} PSM rank {}: spec_e_value differs ({} vs {})",
                    i, j, a.spec_e_value, b.spec_e_value);
            }
        }
    }
    assert!(spectra_with_psms > 0, "no spectra produced PSMs — fixture problem");
    eprintln!("Verified bit-identical output across thread counts on {} spectra with PSMs",
              spectra_with_psms);
}
```

Note the use of `.to_bits()` for f32/f64 equality. Bit-equality is the correct check here — we verified per-spectrum work runs entirely on one worker, so output is deterministic at the bit level.

- [ ] **Step 4.2: Run the test (release mode for speed)**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo test --release -p search --test match_spectra_thread_invariance 2>&1 | tail -10
```

Expected: PASS. Wall time ~30-60 seconds (one BSA search × 2 thread counts).

If it FAILS with `assertion 'left == right' failed: spectrum N PSM rank M: spec_e_value differs`: there is FP non-determinism somewhere. The most likely cause is a Rayon-internal accumulation that is order-dependent. Investigate `compute_spec_e_values_for_spectrum`'s `add_prob_dist` calls — these should all be on the worker's local data, but verify.

- [ ] **Step 4.3: Commit**

```bash
git -C /Users/yperez/work/msgfplus-workspace/astral-speed add \
  rust/crates/search/tests/match_spectra_thread_invariance.rs
git -C /Users/yperez/work/msgfplus-workspace/astral-speed commit -m "test(search): match_spectra output bit-identical across thread counts"
```

---

## Task 5: Update PXD001819 harness to use 12 threads

**Files:**
- Modify: `benchmark/parity/run_pxd001819_2arm.sh`

**NOTE:** `benchmark/` is gitignored. This is local-only tooling. No commit.

- [ ] **Step 5.1: Add `--threads` to RUST_ARGS**

Find the `RUST_ARGS=(` block in the harness (around line ~140). Add `--threads 12` to mirror Java's default thread count (the harness uses `MSGFPLUS_THREADS` for Java, defaulting to 4; PXD001819 fixtures used 12 historically per the project memory).

After the `--max-length 40` line, insert:

```bash
    --threads "${MSGFPLUS_THREADS:-12}"
```

(Reusing the same env var keeps the two arms symmetric.)

- [ ] **Step 5.2: Sanity check the script**

```bash
bash -n /Users/yperez/work/msgfplus-workspace/astral-speed/benchmark/parity/run_pxd001819_2arm.sh && echo "syntax ok"
```

Expected: `syntax ok`.

---

## Task 6: Run PXD001819 end-to-end and capture wall time

This is a measurement task, not a code task. Produces evidence for the post-Rayon perf claim.

- [ ] **Step 6.1: Build release binary**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust && \
  cargo build --release -p msgf-rust 2>&1 | tail -3
```

- [ ] **Step 6.2: Run the harness (Rust arm only — Java baseline already exists)**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/benchmark/parity && \
  bash run_pxd001819_2arm.sh --skip-java 2>&1 | tee /tmp/pxd001819_rust_rayon.log | tail -30
```

Expected wall time: 5-15 min (was 103+ min single-thread). If wall time is > 30 min, parallelization is not effective — stop and investigate.

- [ ] **Step 6.3: Capture metrics**

Compute and record:
- Wall time (from `time` output)
- PSM count (`wc -l rust.pin`)
- Top-1 identity vs Java reference: count rows in `diff_report.tsv` (or just compare PSM totals)
- Speedup ratio vs single-thread baseline (~80 min projected from the 200-spec slice profile)

- [ ] **Step 6.4: Update project memory**

Edit `~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_rust_perf_gap_vs_java.md`:

Append a new section after the existing content:

```markdown

## 2026-05-08 update — post-Rayon measurement

Commit: <commit SHA from Task 2 + Task 3>

PXD001819 wall-time after Rayon parallelization (`--threads 12`):

  - Java MS-GF+ (12 threads): 1m42s
  - Rust msgf-rust (12 threads): <NN>m<NN>s
  - Speedup factor: <Rust / Java>

The remaining gap is per-call GF cost (Fix B in spec
`docs/superpowers/specs/2026-05-08-rayon-match-spectra-design.md`).
Profile target: identify which sub-frames inside `PrimitiveAaGraph::new`
and `compute_inner` account for the 6-12× constant factor.
```

Replace `<NN>m<NN>s` and `<commit SHA>` with the actual values.

---

## Self-review checklist (run before handing off to executor)

- [x] Spec coverage: Tasks 1-3 implement Fix A's parallelization; Task 4 covers the determinism test from the spec; Task 5+6 cover the measurement promise (post-Rayon perf gate).
- [x] No placeholders: every task has exact file paths, code blocks, and commands.
- [x] Type consistency: `TopNQueue`, `PsmMatch`, `SearchParams`, etc. are referenced by their actual types throughout.
- [x] Iteration shipping model honored: 5 milestone commits (Tasks 1-5), no per-task PR. Single closing PR at end of the rewrite remains the user's iteration model.
- [x] Out-of-scope items in the spec (Fix B, streaming candidates, output threading) are NOT in any task here — deferred per the spec.

## Execution handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-08-rayon-match-spectra-plan.md`. Two execution options:

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

Which approach?
