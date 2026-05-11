# GF Tails Reconciliation Iteration 1 — Report

**Date:** 2026-05-10  
**Track:** Track 1 (Diagnostic & Targeted Fixes)  
**Status:** Partial delivery. Gate tightened 0.5 OOM (below the 1.0 OOM DoD target).

## Results Summary

Cumulative fixes from 2026-05-05 through 2026-05-10 (theo_mz formula, cleavage credit, partition Ord, per-partition ions) improved 3 of 5 traced PSMs dramatically, but regressed scan 3353:

| Scan | Peptide | Charge | Pre-iter 1 | Post-iter 1 | Δ | Status |
|---|---|---|---|---|---|---|
| 3416 | KVPQVSTPTLVEVSR | 3 | 0.106 | 0.461 | +0.355 | mildly regressed |
| 3353 | KVPQVSTPTLVEVSR | 3 | 1.010 | **3.276** | **+2.266** | **NEW OUTLIER** |
| 5442 | MYVDPSSPK | 2 | 2.396 | 1.031 | -1.365 | improved |
| 1507 | IPPITPTPK | 2 | 2.862 | 0.192 | -2.670 | improved |
| 2693 | DAFSHGFK | 2 | 3.675 | 0.480 | -3.195 | improved |

**TOLERANCE_LOG10:** 4.0 → **3.5** (tightened 0.5 OOM)

## Commits

Track 1 Tasks 2–4 (diagnostic infrastructure):

- **c6fbb14** — msgf-trace: add MGF support to spectrum reader
- **1962cfb** — msgf-trace: implement --print-score-dist (per-node GF distribution)
- **9d70a9f** — msgf-trace: complete --print-score-dist feature
- **e918376** — Java: GeneratingFunction.getSpectralProbabilityWithTrace gated (matching Rust tracer output)
- **1256269** — Python: diff_gf_distribution.py harness (diff Rust ↔ Java GF distributions, localize first-divergence node)

Track 1 Task 7 (gate adjustment):

- **bc8ee39** — Test: TOLERANCE_LOG10 4.0 → 3.5; refreshed per-PSM divergence table

## What Was Delivered

### Diagnostic infrastructure (Track 1 Tasks 2–4)

1. **Rust tracer:** `msgf-trace --print-score-dist` outputs per-node spectral probability distributions for each candidate in a spectrum
   - Enables side-by-side inspection of Rust GF DP nodes
   - Supports MGF input (T1-2 output)

2. **Java tracer:** `GeneratingFunction.getSpectralProbabilityWithTrace` gate
   - Parallel Java DP trace output (matching Rust structure)
   - Gated for legacy-mode compatibility (T1-3 output)

3. **Python diff harness:** `benchmark/parity/diff_gf_distribution.py`
   - Ingests Rust + Java tracer output
   - Walks both GF trees in lockstep
   - Reports first node where distributions diverge (node type, score, probability, delta)
   - Supports per-candidate or full-spectrum filtering

### Root-cause identification (Track 1 Task 7)

**scan 3353 identified as the new SP-level outlier:**
- Same peptide (KVPQVSTPTLVEVSR ch3) appears on two spectra
- scan 3416: log10 divergence = 0.461 (near parity, ~0.29 OOM above gate)
- scan 3353: log10 divergence = 3.276 (outlier, ~3.3 OOM above gate)
- **Conclusion:** divergence is spectrum-dependent, not peptide-dependent

**Probable root:** spectrum-specific quality factors (peak rank distribution, noise profile) interact with the recent fixes in `compute_inner` probability accumulation. The per-partition ion improvements (commits 2026-05-10) likely changed the weight distribution of edges in the DP, making scan 3353's spectral topology expose a latent numerical or algorithmic difference.

## What's Open

### Tightening gate to DoD (1.0 OOM)

Current gate = 3.5 OOM; DoD requires ≤ 1.0 OOM additional progress. scan 3353's 3.3 OOM contribution is the blocker.

### scan 3353 detailed investigation

Next iteration (Track 1 Task 6a/b/c/d) will:

1. Run `diff_gf_distribution.py` on scan 3353 traces to pinpoint the first-divergence node
2. Inspect the node type and scoring formula
3. Trace back through `compute_inner` to identify the root cause:
   - Per-edge probability accumulation divergence?
   - Score threshold or pruning mismatch?
   - `f32` vs `f64` precision artifact?
   - Floating-point reduction order (Rayon parallelism)?
4. Apply targeted fix (formula adjustment, precision promotion, or threshold alignment)
5. Validate gate tightening on full PXD001819

## Next Steps

1. Invoke `benchmark/parity/diff_gf_distribution.py` with scan 3353 traces
2. Follow diagnostic output to the first-divergence node
3. Map node back to `compute_inner` code location
4. Propose and implement targeted fix per Track 1 Task 6a/b/c/d
5. Measure gate tightening and confirm scan 3353 moves within the 1.0 OOM DoD target
