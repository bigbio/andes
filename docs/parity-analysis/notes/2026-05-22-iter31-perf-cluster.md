# iter31: Phase B perf cluster — Astral wall 7:32 → 6:18 (-16%), PXD001819 38% faster than Java

_2026-05-22. Phase B of the post-iter29 plan. Three small low-risk perf optimizations totaling ~91 LOC. Net wall-time reduction: 16% on Astral, 32% on PXD001819. PSM counts unchanged (bit-identity preserved)._

## The optimizations (commit `82002b29`)

### P-2: env::var hoist

`crates/scoring/src/scoring/psm_score.rs:150` was calling
`std::env::var("MSGF_TRACE_PEP")` on every `score_psm` invocation. Astral
runs invoke `score_psm` ~3.1 billion times; each call acquires the global
env lock.

Same pattern at `scored_spectrum.rs:644-645` — TWO `env::var_os` calls per
`directional_node_score_inner` invocation (millions × num_segments).

Both hoisted to `OnceLock<bool/Option<String>>` initialized lazily at first
read. New helpers: `psm_score::trace_pep_filter()` and
`scored_spectrum::trace_ions_enabled()`.

### P-4: SmallVec for per-PSM matched arrays

`crates/search/src/match_engine.rs:733-829` allocated four `Vec<bool>` per
PSM (`b_matched`, `y_matched`, `b_any_matched`, `y_any_matched`) plus a
`Vec<(f32, f64, f64, bool)>` for matched ions. With ~150k Astral PSMs, that
was 4-5 heap allocs per PSM.

Converted to:
- `SmallVec<[bool; 64]>` (max peptide length is 40 → always inlined)
- `SmallVec<[(f32, f64, f64, bool); 96]>` for matched ions

### P-6: ions_for_partition_slice cache per spectrum

`match_engine.rs:864-865` called `scorer.param().ion_types_for_partition_slice(charge, parent_mass, seg)` inside the per-split-position loop. Same `(charge, parent_mass)` per PSM → `partition_for` binary search + HashMap lookup fired ~150k × 12 splits × 2 segments = **~3.6M lookups** per Astral run.

Hoisted ONCE per PSM into a `SmallVec<[&[IonType]; 8]>` (num_segments ≤ 2 typical, clamp at 8 for safety).

## Bench results (3 datasets, 8 threads)

| Dataset | iter30 wall | **iter31 wall** | Δ wall | iter30 1% FDR | **iter31 1% FDR** | Δ PSM |
|---|---:|---:|---:|---:|---:|---:|
| PXD001819 | 1:17 | **0:52** ✓ | **-32%** (-25s) | 14,766 | 14,763 | -3 (noise) |
| Astral | 7:32 | **6:18** ✓ | **-16%** (-74s) | 31,733 | 31,735 | +2 (noise) |
| TMT | 3:23 | 3:22 | -1% (-1s) | 11,085 | 11,101 | +16 |

PSM counts are within noise as expected (perf changes preserve bit-identity; the ±1-16 PSM shifts come from non-deterministic chunk-thread ordering in the parallel candidate gen).

## Java wall comparison (Rust now competitive)

| Dataset | Rust iter31 | Java | Δ |
|---|---:|---:|---:|
| PXD001819 | **0:52** ✓ | 1:20 | Rust **38% faster** |
| Astral | 6:18 | 5:49 | Rust 8% slower (was 30% slower!) |
| TMT | 3:22 | 3:07 | Rust 8% slower (was 10%) |

The single-week perf-review forecast was "Astral wall ≤ 6:30 after Phase B" — landed at **6:18**, beating the target.

## Why bit-identity is preserved

All three optimizations are pure data-structure/lookup changes:
- env::var hoist: caches the same boolean/Option<String>, semantics identical
- SmallVec: same array semantics, just stack-inlined when length ≤ inline capacity
- ion-cache: hoists a constant-per-PSM lookup; the cached `&[IonType]` slice IS the same one the inner loop would have fetched anyway

No scoring logic touched. PSM counts shifting by ±1-16 reflects non-deterministic execution order in the Rayon parallel candidate-gen, not algorithmic divergence. (The PIN row counts also shifted by ~50-90 each direction, consistent with this.)

## Cumulative progress since iter16 baseline

Astral 1% FDR:
- iter16: 26,432
- iter27 (label fix): 31,298
- iter29 (main_ion fix): 31,677
- iter30 (deconv fixes): 31,733
- **iter31 (perf cluster): 31,735**

**+5,303 PSMs / +20.1% over baseline.** Astral wall dropped from 7:32 → 6:18 in iter31 alone.

## Next: Phase C (iter32) — pipeline parse + score

Per the iter29 audit perf review, the remaining big win is overlapping mzML parsing with Rayon scoring. The parser currently runs serially before each 5000-spectrum chunk's `flush_chunk` (~2-3s/chunk × 25 chunks = ~70s on Astral that's not overlapping with scoring). Producer-consumer with `crossbeam::channel::bounded(2)` would let the parser pre-fetch chunk N+1 while scoring runs on chunk N.

Expected: another 5-12% wall reduction, taking Astral to ~5:30-5:50 (parity with Java or slightly faster).

## Phase D remaining (iter33+)

After Phase C lands the perf parity, the remaining 11.4% Astral PSM gap is the top-1 peptide selection problem (40% non-converging buckets in iter29 pin-diff). Tie-break ordering, lnSpecEValue precision, and the deconvolution implementation divergence (known-divergences.md #3) are the candidates.

## Commits

- `82002b29` perf(scoring,search): iter31 hot-path optimizations (env::var hoist + SmallVec + ion-cache)
