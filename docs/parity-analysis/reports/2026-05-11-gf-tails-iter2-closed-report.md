# GF Tails Parity — Iteration 2 Closed Report

**Date:** 2026-05-11  
**Commit:** 3b0bef1 (test fix)  
**Result:** 1.0 OOM SP-vs-SP gate met on all 5 BSA fixture PSMs

## Results

| scan | peptide | charge | Java SP | Rust SP | log10 SP diff |
|---|---|---|---|---|
| 3416 | KVPQVSTPTLVEVSR | 3 | 3.005e-9 | 5.220e-9 | 0.240 |
| 3353 | KVPQVSTPTLVEVSR | 3 | 4.658e-10 | 3.473e-10 | 0.127 |
| 5442 | LGEYGFQNALIVR | 2 | 4.315e-7 | 2.752e-6 | 0.805 |
| 1507 | YLYEIAR | 2 | 5.246e-4 | 2.914e-4 | 0.255 |
| 2693 | SLGKVGTR | 2 | 1.392e-3 | 1.652e-3 | 0.074 |

**TOLERANCE_LOG10: 1.0 OOM** (all PSMs pass).

## Key Finding: Test Unit Mismatch

Iter 1 reported scan 3353 as an outlier at 3.276 OOM, suggesting a regression in GF DP accuracy. Iter 2 revealed this was a **test unit artifact**: the test was comparing Rust **spectral probability** (`psm.spec_e_value` field, historically misnamed) against Java **SEV** (SP × `num_distinct_peptides_at_length`).

After fixing the test to compare SP-vs-SP (using Java `spec_prob` values captured via `-Dmsgfplus.gftrace=true`), scan 3353's true divergence is **0.127 OOM**, well within the 1.0 OOM gate.

## Implication

The GF DP (dynamic programming) agreement between Rust and Java is much tighter than iter 1 measurements suggested. The original Design of Defects (DoD) target of **1.0 OOM tolerance on spectral probability** is now achieved.

The remaining **SEV-level gap** (observed in full SpecEValue comparisons) is correctly reattributed to the `num_distinct_peptides_at_length` factor, which has its own follow-up (item #2 in `known-divergences.md`: mod-aware distinct counting).

## Parity Status

- **SP-vs-SP parity:** 1.0 OOM gate ✓ closed (all 5 fixture PSMs)
- **SEV-level parity:** deferred to `num_distinct` follow-up (item #2)
- **Next work:** if needed, debug scan 5442 (new outlier at 0.805 OOM) via `diff_gf_distribution.py`

See the previous iter 1 report ([2026-05-10-gf-tails-iter1-report.md](2026-05-10-gf-tails-iter1-report.md)) for diagnostic infrastructure context.
