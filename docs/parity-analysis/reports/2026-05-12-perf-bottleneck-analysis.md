# PXD001819 Wall-Time Bottleneck Analysis (2026-05-12)

**Status:** Combined finding from instrumented PHASE timing + targeted code reading.

**Context:** Today's Rust↔Java wall comparison on full PXD001819:
- **Java MS-GF+:** 64.92s real (12 threads, Mac)
- **Rust (post-`507bcb1`):** 356.04s real (12 threads, Mac)
- Gap: **5.49× slower wall** — the iteration's 2× MUST target (≤ 32.5s; or 178s relative to Rust's 356s) is not met.

## Phase wall breakdown (post-`507bcb1`)

Instrumentation: `[PHASE …]` markers added to `msgf-rust.rs` (`fasta_load`,
`search_index_build`, `param_and_scorer`, `spectra_load`, `match_spectra`,
`pin_write`, `TOTAL`).

| Phase | Wall | % of total | Serial / parallel |
|---|---|---|---|
| `fasta_load` | 0.02s | <1% | serial |
| `search_index_build` (target+decoy + SA) | 0.85s | <1% | serial |
| `param_and_scorer` | 0.00s | <1% | serial |
| `spectra_load` (mzML) | 6.87s | 2% | serial |
| **`match_spectra`** | **189.40s** | **53%** | parallel (Rayon, 12 threads) |
| **`pin_write`** | **158.90s** | **45%** | serial |
| **TOTAL** | **356.04s** | | |

PSM-preservation gate: Percolator @ 1% FDR = **14,850** (PSM-identical to
the pre-iteration baseline). Pin rows: 37,113.

## Three load-bearing root causes (code reading)

### #1 — Rust is missing Java's `FastScorer` precomputation

**Java path:**
- `FastScorer.java:28` precomputes two arrays per spectrum/charge:
  - `prefixScore[nominal_mass]`
  - `suffixScore[nominal_mass]`
- `FastScorer.java:68` per-split scoring is then ARRAY LOOKUPS:
  `score += prefixScore[i] + suffixScore[peptideMass - i]`.

**Rust path:**
- `psm_score.rs:34` (`score_psm`) iterates every peptide split.
- Each split calls back into `ScoredSpectrum::node_score` →
  `directional_node_score` at `scored_spectrum.rs:399`.
- `directional_node_score` re-iterates segments, ion types, peak lookups
  per split per candidate.

**Cost magnitude:** for PXD001819 with 5,463,285 PSMs pushed, each PSM is
typically scored at ~10 splits → ~54M `directional_node_score` calls,
each doing per-segment partition lookups + per-ion peak rank lookups.
Java does the SAME total work ONCE per spectrum (precomputing the array)
and amortizes to O(1) per split. Estimated saving on `match_spectra`: 3–5×.

**Files:** `rust/crates/scoring/src/scoring/psm_score.rs`,
`rust/crates/scoring/src/scoring/scored_spectrum.rs`. Add a
`FastSpectrumScores` struct on `ScoredSpectrum` carrying
`prefix_score: Vec<i32>` / `suffix_score: Vec<i32>`, populated once at
construction (loop over `1..peptide_mass`); `score_psm` becomes the array
lookup.

### #2 — Rust runs `compute_psm_features` BEFORE top-N retention

**Rust path:**
- `match_engine.rs:231` calls `compute_psm_features(...)` for every
  candidate before pushing to the `TopNQueue`.
- Feature extraction (`match_engine.rs:501`): predicts b/y ions, walks
  peak vectors, allocates `Vec<bool>` / `Vec<(mz, intensity, ...)>`,
  partial-sorts matched ions, computes per-ion error stats.

**The waste:** 5,463,285 PSMs pushed → 37,112 spectra with non-empty
queue → ~37,112 survivors (top-N=1). **99.3% of feature work is on
candidates that get evicted.**

**Java equivalent:** features computed in the writer path, after top-N
selection (`MSGFPlusMatch.java`).

**Estimated saving:** ~30–80s on `match_spectra` (depending on what
fraction of the 189s is feature work vs scoring).

**Files:** `rust/crates/search/src/match_engine.rs` line ~231 (the
inner-loop call site) + a new "lazy features" pattern: defer to a second
pass that iterates only the retained top-N PSMs per spectrum.

### #3 — Rust rebuilds the target+decoy SA from FASTA every run

**Java path:**
- `CompactSuffixArray.java:89` reads cached `.cseq` / `.canno` /
  `.csarr` / `.cnlcp` files when they exist on disk; writes them on
  first run.

**Rust path:**
- `msgf-rust.rs:148` calls `SearchIndex::from_target_db(...)`.
- `search_index.rs:48` builds `CompactFastaSequence` + decoys + suffix array.
- `suffix_array.rs:41` runs SA-IS over the compact bytes.

**Cost magnitude:** measured **0.85s** in `search_index_build` on
PXD001819. At <1% of total wall, this is **not** a load-bearing
bottleneck despite the architectural difference. **Defer.**

## Other findings (already addressed or smaller)

| # | Finding | Status |
|---|---|---|
| A | PIN Label cache by peptide sequence (pin.rs L118, L257) | **Landed in `507bcb1`** — dropped full wall ~626s → 356s |
| B | One-pass distinct-peptide count folded into candidate collection (match_engine.rs L40, search_index.rs L132) | **Landed in `507bcb1`** |
| C | `ScoredSpectrum::new` clones `partition_ion_logs.to_vec()` per spectrum/charge (scored_spectrum.rs L161) | Open — ~Vec<f32> clone per (segment × ion). Smaller constant-factor win |
| D | `pin_write` per-row formatting at 158.90s (45% of wall) | Open — Label cache landed; remaining is TSV format / I/O. Worth investigating |

## Issue → phase mapping

```
TOTAL 356s
├── match_spectra 189s  ← Issue #1 (FastScorer), Issue #2 (feature hoist)
│                          and minor C (segment clone)
├── pin_write 159s     ← Major finding D (TSV format remains)
└── other 8s            ← Issue #3 (cached SA: <1s today; defer)
```

## Lever ranking (best → worst by phase data)

1. **Issue #1 — FastScorer precompute** — biggest single search-side
   win; estimated 3–5× on `match_spectra` → match_spectra 189s → ~40–60s
2. **Issue #2 — Hoist `compute_psm_features` to post-top-N** — kills
   ~99% of wasted feature work; estimated 30–80s on `match_spectra`
3. **D — `pin_write` further optimization** — 45% of remaining wall
   (159s). Investigate parallel write OR cheaper per-row format
4. **C — Drop `partition_ion_logs.to_vec()` cloning in ScoredSpectrum::new**
   — small constant-factor; pairs naturally with #1's refactor
5. **Issue #3 — Cached SA artifacts** — **<1s today**; defer

## 2× MUST gate (wall ≤ 178s) — feasibility

Combining levers #1 + #2 + D:
- match_spectra: 189s → ~30–60s if #1 + #2 land cleanly
- pin_write: 159s → ~20–60s if D lands (parallel/cheaper write)
- other: ~8s

Realistic total: ~60–130s = **2.5–6× speedup** vs current 356s →
**hits 2× MUST and likely beats 2.5× STRETCH**.

## Recommended next iteration

**Parallel work, 3 tracks, can land independently:**

1. **Track A — FastScorer port (Issue #1)**
   `prefix_score[]` / `suffix_score[]` precompute on `ScoredSpectrum`,
   `score_psm` becomes array lookup. Java parity via the existing
   `gf_java_parity` + `match_engine_java_parity` tests. Bit-identity
   gate (sums must match). Highest single-lever payoff.

2. **Track B — Hoist `compute_psm_features` to post-top-N (Issue #2)**
   Two-pass: enumerate + score + push to TopNQueue (cheap features
   only), then after Phase 6 SpecEValue, iterate retained PSMs to fill
   in the heavyweight features. Bit-identity gate trivial — same
   features, fewer calls.

3. **Track C — `pin_write` parallelization / per-row optimization (D)**
   Profile what's actually slow inside pin_write (Label cache is in
   place; remaining 159s could be TSV format, BufWriter flushes, or
   per-row Vec allocations). Investigate parallel write via Rayon
   per-chunk OR replace per-row String allocations with reusable
   buffer.

Defer Issue #3 (cached SA) since phase data shows it's <1% of wall.

Defer Finding C (ScoredSpectrum clone) until Tracks A+B done — it pairs
naturally with the FastScorer refactor.
