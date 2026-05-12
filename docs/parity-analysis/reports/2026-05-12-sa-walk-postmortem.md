# SA-walk Refactor Postmortem (2026-05-11/12)

**Status:** Tasks 1–5 + 6a–6d shipped; Task 6e+6f+6g REVERTED.

## Final state on `rust-implement`

| Commit | Item | Status |
|---|---|---|
| `9c56797` | Task 1 — PrimitiveAaGraph arena pool | ✅ shipped |
| `b95d348` | Task 2 — ScoreDist flat arena | ✅ shipped |
| `317d6bc` | Task 3 — per-segment partition cache | ✅ shipped |
| `e19293a` | Task 5 — chunked add_prob_dist | ✅ shipped |
| `e5361c6` | Task 6a+6b — DistinctPeptide + SaPeptideStream | ✅ infrastructure (unused in prod) |
| `4f927b7` | Task 6c+6d — Met-cleavage + diff_pin_psms harness | ✅ infrastructure |
| _(working tree)_ | Task 6e+6f+6g — production wiring | ❌ REVERTED — gates failed |

## Gate results

Iteration target: **2× wall (≤ 3m57s)** AND **PSM-identical** AND **Percolator @ 1% FDR ≥ 14,850**.

| Gate | Pre-iter (iter 3) | After items 2+3+5 (Mac 12t) | After SA-walk integration (VM 8t) | Final ship |
|---|---|---|---|---|
| Total wall | 7m53s | **7m18s (−7.4%)** | 7m07s VM (+32s vs legacy) | 7m18s |
| match_spectra wall | 199s | **160.64s (−19%)** | 123.82s | 160.64s |
| Pin row count | 37,113 | 37,113 | 37,113 | 37,113 |
| Top-1 PSM identity vs prior | identical | identical | **99.087%** (339 disagreements) | identical |
| Percolator @ 1% FDR | 14,850 | 14,850 | **14,698** (−152) | 14,850 |
| CPU% (parallelism) | ~200% | similar | 146% (regressed) | ~200% |

## What caused the gate failures

339 of 37,112 scans had different top-1 peptides between protein-walk and
SA-walk. The disagreement pattern is consistent: protein-walk picks
Met-cleaved variants (`M.X…` flanking notation), SA-walk picks
non-Met-cleaved alternatives at the same scan.

Five examples:

| scan | protein-walk top-1 | SA-walk top-1 |
|---|---|---|
| 24840 | `M.KYM+15.99491GSFLRK.A` | `R.NMGVHITFVK.S` |
| 17125 | `M.VAFTVDQMR.S` | `R.GTPVYNYPR.T` |
| 1954 | `M.TISSAHPETEPK.W` | `K.TLDFDYAVQPK.G` |
| 24658 | `M.KYM+15.99491GSFLRK.A` | `R.NMGVHITFVK.S` |
| 11401 | `M.RMTTELDDLR.R` | `K.NELFGPSFPNK.T` |

The Met-cleavage handling in `SaPeptideStream` (commit `4f927b7`) was
correct enough to PASS the in-test fixture (a single M-prefixed protein
with controlled residue overlap), but the integration into the production
candidate stream introduced enough subtle divergence to flip 0.91% of
top-1 selections.

The lab CPU profile contradicted the spec's "serial pre-loop bottleneck"
hypothesis too: the post-SA-walk run had MORE serial time on the VM
(304s) than pre-SA-walk (278s estimated), so even if PSM identity had
held, the wall would have regressed.

## Lessons (canonical list)

1. **Profile before assuming bottlenecks.** The spec assumed
   `enumerate_candidates` Vec materialization was a major serial cost. CPU
   accounting after the change suggests it wasn't.
2. **PSM-loss gates must be same-machine comparisons.** Cross-machine
   diffs introduce float-precision and threading-order artifacts that
   confound the analysis.
3. **Met-cleavage is the load-bearing edge case** for any SA-walk approach.
   Get it right on a controlled fixture FIRST, with bit-identity vs the
   protein-walk on the same fixture, before any production wiring.
4. **The subagent-driven workflow correctly caught the regression at the
   end-to-end gate.** The implementer's lib + parity tests all passed; only
   `diff_pin_psms.py` (Task 6d, which was scoped exactly for this) revealed
   the issue.

## What stays useful

The SA-walk infrastructure (`sa_walk.rs`, `distinct_peptide.rs`,
`diff_pin_psms.py`) is committed and works in isolation. Future iterations
that revisit candidate-generation refactors can build on it without
starting from scratch — but should first profile and confirm the
bottleneck hypothesis, and tackle Met-cleavage parity as the first
feature, not the last.

## Recommended next steps

- **Defer further SA-walk work** until profiling shows the bottleneck IS in
  candidate enumeration (it currently appears to be elsewhere).
- **Investigate the actual serial-time culprit** with `cargo flamegraph`
  or `samply` — pin writing? FASTA load? Param parse? Something else?
- **Treat the 6 shipped commits as a self-contained iteration:** ~7% wall
  improvement, PSM-identical, new infrastructure for future use.
