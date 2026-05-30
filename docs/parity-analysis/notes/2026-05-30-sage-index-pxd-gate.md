# Sage-style index — PXD gate (low-res CID): NO degeneration, but low recall

**Date:** 2026-05-30
**Branch:** `feat/chimeric-dda-plus` (B-T1..T3 + recall-insurance 246008694)

| arm | rows | wall | @1% PSMs |
|---|---:|---:|---:|
| brute (`--chimeric-frag-index off`) | 365,661 | 1:28 | 18,234 |
| SageIndex (`on`, TOP_K=128, tol 0.05) | 286,267 | 2:00 | 11,432 |

Entrapment FDP (SageIndex on): 0.83% raw / 1.66% combined — still controlled.

## Read
- **Speed: NO degeneration** (2:00 vs Approach A's 28:24). The local sub-ms
  microbenchmark held on real data. SageIndex is slightly SLOWER than brute on PXD
  only because PXD is tiny/low-res — brute is already 1:28, so index build+query
  overhead isn't repaid. Per design, PXD (low-res) gates to brute; the index targets
  Astral.
- **Recall 63%** (11,432/18,234) — fails the ≥99.5% gate. Low-res CID: sparse
  fragments → coarse matched-fragment-count voting → true peptide misses top-128.
- Decisive test = Astral high-res (fine bins → discriminative voting; brute 18:30
  with thousands of candidates/spectrum → index should win speed big).
