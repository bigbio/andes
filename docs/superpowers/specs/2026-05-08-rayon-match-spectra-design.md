# Fix A — Rayon-parallel `match_spectra` design

**Status:** drafted 2026-05-08

**Goal:** Parallelize the per-spectrum scoring loop in `search/src/match_engine.rs::match_spectra` over a worker thread pool. Mirrors Java MS-GF+'s `-thread N` model. Targets the dominant 1/12 utilization gap measured on PXD001819 (Rust 103+ min single-thread vs Java 102 s on 12 threads).

## Motivation

Profiling on a 200-spectrum PXD001819 slice (commit `d4667d3`) showed:

- **GF construction = 50.7%** of single-threaded CPU time, 92.8 ms/spectrum
- **Per-PSM scoring = 18.6%** (`score_psm` + node/peak lookup)
- Full-dataset projection: ~81 min Rust single-thread vs ~5.85 min Java single-thread (≈14× per-thread)

Java's `DBScanner.computeSpecEValues` runs the same outer "for each spectrum" loop with the same per-spectrum work (mass-window GF construction, `setUpScoreThreshold`, `gf.accept`). Java's wall-time win comes from **threading**, not from caching or skipping. The Rust port already does the algorithmic equivalent of Java's per-spectrum work — it just runs single-threaded.

This design adds Rayon-based parallelism over the outer `match_spectra` loop. Java-aligned: same per-spectrum work, executed in parallel.

## Scope

**In scope:**
- Parallelize the per-spectrum loop in `match_spectra`.
- Add a `--threads N` CLI flag on `msgf-rust` (default: `num_cpus::get()`).
- Plumb the thread count from the CLI through to the search.
- Verify deterministic outputs across thread counts (PIN row order may differ but PSM identity must be invariant).

**Out of scope (Fix B, later):**
- Per-call GF cost optimization (constant-factor work inside `PrimitiveAaGraph::new` / `compute_inner`).
- Caching, empty-bin skipping, or any change that diverges from Java's per-spectrum, per-mass-index, fresh-build pattern.
- Threading inside the per-spectrum work (e.g., parallelizing the mass-window GF loop). Java doesn't do this either.

## Architecture

### Current (single-thread)

```
match_spectra(spectra, idx, params, scorer, ...) -> Vec<TopNQueue>
  ├─ build candidates: Vec<Candidate>          (one-time, before loop)
  ├─ build bucket_index: BTreeMap<i32, Vec<usize>>
  ├─ register enzyme on a cloned aa_set
  └─ for spec_idx in 0..spectra.len() {           ← parallel target
        ├─ build per-charge ScoredSpectrum
        ├─ compute candidate window via bucket_index
        ├─ per-candidate score loop → push to queues[spec_idx]
        └─ compute_spec_e_values_for_spectrum(...)
     }
```

### Proposed (parallel)

The outer loop becomes `par_iter().enumerate().map(...)` collecting back into `Vec<TopNQueue>`:

```rust
let queues: Vec<TopNQueue> = spectra
    .par_iter()
    .enumerate()
    .map(|(spec_idx, spec)| {
        let mut queue = TopNQueue::new(params.top_n_psms_per_spectrum);
        // ... existing per-spectrum body ...
        queue
    })
    .collect();
```

Rayon's `par_iter` distributes work across a thread pool (default = num CPUs). The pool is configured once via `rayon::ThreadPoolBuilder` from the `--threads` CLI flag.

### Sharing model

Each per-spectrum body reads:

| Resource | Sharing | Notes |
|---|---|---|
| `spectra: &[Spectrum]` | `&` shared, immutable | safe |
| `idx: &SearchIndex` | `&` shared, immutable | safe |
| `params: &SearchParams` | `&` shared, immutable | safe |
| `scorer: &RankScorer` | `&` shared, immutable | safe |
| `candidates: &[Candidate]` | `&` shared, immutable (built once before parallel loop) | safe |
| `bucket_index: &BTreeMap<i32, Vec<usize>>` | `&` shared, immutable (built once) | safe |
| `aa_set_for_gf` | `&` shared, immutable (cloned + enzyme-registered before parallel loop) | safe |
| `decoy_prefix: &str` | `&` shared, immutable | safe |

Per-thread (allocated inside the closure):

- `queue: TopNQueue` — owned, returned via `collect()`.
- `scored_per_charge: HashMap<u8, ScoredSpectrum<'_>>` — owned, lives for the closure's duration.
- `GeneratingFunctionGroup` and per-bin `PrimitiveAaGraph` / `GeneratingFunction` — built and dropped within the closure.

**No `&mut` shared state required.** This is a textbook Rayon parallelization: the inputs are read-only and outputs are owned per task.

### Output ordering

`par_iter().enumerate().map(...).collect::<Vec<_>>()` preserves index order — `queues[spec_idx]` corresponds to `spectra[spec_idx]`. Downstream PIN/TSV writers expect this; no change needed there.

The PIN row order may differ from Java's row order if Java's `-thread N` uses a different work-distribution scheme, but **PSM identity per spectrum is invariant** across thread counts. The `gf_bsa_parity` and `match_engine_java_parity` tests should pass with any thread count.

### CLI surface

Add to `msgf-rust/src/bin/msgf-rust.rs`:

```rust
#[arg(long, default_value_t = num_cpus::get())]
threads: usize,
```

In `run()`, before invoking `match_spectra`, configure the global Rayon pool:

```rust
rayon::ThreadPoolBuilder::new()
    .num_threads(cli.threads)
    .build_global()?;
```

(Or use `install` on a per-`match_spectra` pool if a global pool is undesirable for library use; for the CLI, global is fine.)

## Dependencies

Add to `search/Cargo.toml`:

```toml
rayon = "1.10"
```

Add to `msgf-rust/Cargo.toml`:

```toml
num_cpus = "1.16"
rayon = "1.10"
```

Both are stable, low-churn dependencies already in the broad Rust ecosystem.

## Testing strategy

### Determinism gates

The existing parity tests must pass under multiple thread counts:

- `cli_smoke` — re-run with `--threads 1` and `--threads 4`; output PSM count must be identical (435/5760 for BSA).
- `gf_bsa_parity` — histogram (1 OOM 68.2%, 2 OOM 92.6%, etc.) must be invariant across thread counts.
- `match_engine_java_parity` — top-1 identity rate must be invariant (214/217 = 98.6%).

A new test `match_spectra_thread_count_invariance` runs the BSA fixture once with 1 thread and once with 4, asserting per-spectrum PSM identity AND `spec_e_value` are bit-identical. Since each spectrum's full pipeline (scoring + GF + spec_e_value assignment) runs entirely on a single Rayon worker, there is no FP-accumulation non-determinism across thread counts — only wall time changes.

### Unit-test invariants

- A unit test asserting `match_spectra(...)` produces the same `Vec<TopNQueue>` regardless of `--threads`.
- A unit test that `rayon::ThreadPoolBuilder::new().num_threads(N).build_global()` is called exactly once (idempotent guard if `--threads` is parsed multiple times in tests).

### Performance gate (informational, not a test)

Run the PXD001819 harness post-Rayon. Expected wall time:

- Single-threaded baseline: ~80 min (projected from 200-spec slice profile)
- 12 threads: ~7-10 min target (if scaling is near-linear)
- Java 12-thread reference: 102 s

If wall time > 15 min on 12 threads, the parallelization isn't effective and we need to investigate (likely culprit: a serial bottleneck like the upfront `enumerate_candidates.collect()` taking ~10 s of the budget).

## Risks

1. **`enumerate_candidates` is called once before the parallel loop and takes ~10 s on PXD001819.** That's a serial bottleneck not addressed by Rayon. Mitigation: leave for now; it amortizes to <1% of total wall on the full dataset. If it dominates after Rayon, addressed by Fix B or a future fragment-index iteration.

2. **`build_global()` panics if called twice.** Mitigation: guard with `OnceLock` or `Once` so multiple test cases or repeated CLI invocations within a test process don't blow up. Use `try_*` or `is_initialized` checks.

3. **Output non-determinism from FP non-associativity in the GF DP.** Mitigation: tests that check spec_e_value compare with ε-tolerance; tests that check peptide identity are exact.

4. **Memory pressure scales with thread count.** Each worker holds a `GeneratingFunctionGroup` mid-search. With 12 workers × ~10-50 MB each, total memory could be 120-600 MB on top of the ~1.3 GB candidate index. Should still fit comfortably on standard hardware.

## Migration steps (for the writing-plans skill)

1. Add `rayon` + `num_cpus` deps to `search` and `msgf-rust` Cargo.toml files.
2. Convert the outer `for (spec_idx, spec) in spectra.iter().enumerate()` loop in `match_spectra` to `par_iter`.
3. Add `--threads N` CLI flag, configure global Rayon pool from it.
4. Add the determinism + invariance unit tests.
5. Run the full integration suite + cli_smoke + parity tests under both `--threads 1` and `--threads 4`.
6. Run PXD001819 harness, capture wall-time + PSM-count parity vs Java.
7. Update project memory `project_rust_perf_gap_vs_java.md` with the post-Rayon measurement.

## Out-of-scope follow-ups (Fix B and beyond)

- Per-call GF cost diagnosis (Fix B): targeted profile inside `PrimitiveAaGraph::new` and `compute_inner` to find the constant-factor 6-12× per-thread gap.
- Streaming candidate enumeration (eliminating the upfront `Vec<Candidate>` materialization).
- Output writer threading (PIN/TSV writes are sequential; not the bottleneck).
- `-thread` semantics matching Java's task subdivision (Java has `-tasks N` separate from `-thread N`; Rust uses Rayon's work-stealing which is conceptually different but functionally equivalent).
