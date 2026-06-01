# GF/SpecE parity audit (external AI) — bench results (2026-06-01)

An external AI audit proposed 6 P0 "make SpecEValue more Java-faithful" fixes
(SinkUnreachable retry removal, round->truncate, max_score-1 guard, precursor-filter
ranking, features-on-deconv-peaks, sentinel PSMs) + safe P1/P2 perf/robustness items.

## Safe P1/P2 — SHIPPED (commit cb808ce3)
Perf micro-opts (drop per-spectrum clone, isotope range clone, O(k^2)->FxHashSet dedup)
+ robustness (get_probability OOB->0.0 defensive, DROPPED_NODES counter, test/doc
cleanups). Behavior-neutral (parity + schema tests pass). No FDR risk.

## P0.4 (precursor-filter parity) — BENCH-VALIDATED, REVERTED
Implemented Java semantics: precursor-filtered peaks rank LAST but remain matchable
(effective intensity 0), instead of Rust's u32::MAX "invisible". Pervasive change to
scored_spectrum.rs peak-pickers + deconv.
Bench (--chimeric, vs baseline Astral 71,839 / PXD 16,552 / TMT 9,671):
  Astral 71,907 (+68, noise; entrapment FDP 1.04%->1.08%)
  PXD    16,695 (+143; FDP 1.13%->1.03% cleaner)
  TMT     9,579 (-92, -0.95%)  <-- REGRESSES THE GATE BLOCKER (CID)
Reverted. n=9 confirmation: even the audit's strongest parity fix (with a supporting
per-scan trace) does not improve aggregate Percolator FDR and HURTS the TMT blocker.
The reviewer's single-scan trace (34306 RawScore 79->80) did not generalize.

## Remaining P0 (round/truncate, max_score guard, sink retry, features-deconv)
Same class, lower confidence than P0.4. Per n=9, expected neutral-to-regressive.
Recommendation: do NOT pursue the P0 parity tweaks; the TMT CID gap needs the deferred
per-ion CID scoring trace + Java instrumentation, not incremental SpecE parity edits.

## P0.6 (sentinel GF PSMs) — ALREADY IMPLEMENTED (reviewer finding stale)
The safe "count + per-run log" the audit asked for already exists:
`GF_SPECTRA_NO_GROUP` (release-safe static AtomicU64) is incremented at the
`!group.is_computed()` early return and printed in the per-run "GF diagnostics
(cumulative): ... {} spectra with no successful bin" line. No change needed.
The remaining options (fallback SpecE policy / reduce-failures) are scoring
changes and were not pursued.

## Decision (user, 2026-06-01)
STOP the P0 parity grind. Bank the safe P1/P2 batch (shipped, commit cb808ce3).
Do NOT pursue P0.2 (round->truncate) / P0.3 (max_score guard) / P0.5
(features-on-deconv) / P0.1 (sink-retry): same class as P0.4 which bench-regressed
the TMT blocker; P0.3 additionally requires undoing a deliberate anti-inversion
guard and risks spec_e=0 -> lnSpecE=-f64::MAX Percolator outliers. The TMT CID gap
needs the deferred per-ion CID scoring trace + Java instrumentation, not these
incremental SpecE-parity edits.
