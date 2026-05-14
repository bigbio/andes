# score_psm under-scoring — diagnosis

**Date:** 2026-05-13/14
**Bug:** Rust's pin column 7 RawScore is 1/3 to 1/5 of Java's for the same (peptide, scan, charge) on PXD001819.
**Canary case:** scan=28787, peptide IVNEEFDQLEEDTPVYK, charge=2 → Java RawScore=297, Rust HEAD=108.

---

## Phase 1: Bisect outcome — strategy invalidated

The original spec/plan assumed the bug was a regression introduced somewhere in the rust-implement commit window (May 4–12, 2026) and that a `git bisect run` over that range against a single-PSM oracle would surface the bad commit. **This assumption is false.**

| Commit | scan=28787 RawScore | DeNovoScore | Verdict |
|---|---:|---:|---|
| `ab28821` (fix/score-psm-undercount HEAD) | 108 | (not captured here) | bug present |
| `5d912fc` (2026-05-11 "GF tails iter 2 closed") | **61** | **97** | worse than HEAD |
| `<earlier commits>` | — | — | not tested (DeNovoScore ceiling at 5d912fc proves the floor) |

Critical observation at 5d912fc: **DeNovoScore=97**. DeNovoScore is the theoretical maximum RawScore that Rust's scoring model can produce for the (spectrum, peptide-mass-window) pair (by construction it's the highest-score path through the GF graph). Rust at 5d912fc literally cannot reach 297 even with a perfect path — the GF graph itself is producing a score range whose top is ~97. By contrast Java reports RawScore=297 for the same PSM, so Java's graph/scoring produces values in a different scale or with substantially more contributing edges.

This is consistent with the user's project memory entry "2026-05-10: Rust↔Java reached 1% FDR parity on PXD001819 (14,839 vs 14,798). RawScore agreement now exact" referring to the **GF-internal** RawScore (the integer score along the SpecEValue graph, used as the SEV index) — NOT the pin's RawScore column, which is computed differently (via `score_psm` → `ScoredSpectrum`'s `directional_node_score` summed across splits). The two scores share a name but are different quantities at different points in the pipeline.

## Phase 2: Pivot — static Java↔Rust comparison

Bisect is abandoned. Replace with side-by-side reading of:

- **Java pin RawScore source:** `src/main/java/edu/ucsd/msjava/.../FastScorer.java` plus its callers in the DBScanner pipeline that emit pin column 7. The PIN writer is `src/main/java/edu/ucsd/msjava/mzid/DirectPinWriter.java` (project CLAUDE.md confirms this). Trace from `DirectPinWriter.writeRow` back to where RawScore is computed for that row.
- **Rust pin RawScore source:** `rust/crates/scoring/src/scoring/psm_score.rs::score_psm` and `rust/crates/scoring/src/scoring/scored_spectrum.rs::directional_node_score`. The PIN writer is `rust/crates/output/src/pin.rs`. Trace from `write_psm_row` back.

The two sources are written by the same author but the Rust port is incomplete: some scoring contributions (ion types, prefix/suffix accounting, sequence-edge bonuses) appear to be missing or differently scaled.

### Likely root causes (Phase-3 hypotheses to test)

1. **Ion-type set mismatch.** Java iterates a wider partition of ion types (e.g. b, b2+, y, y2+, internal, immonium, neutral-loss variants) than Rust's `directional_node_score_inner`. Each missing ion class drops a chunk of score contribution. A 3× factor is consistent with Rust evaluating ~1/3 of Java's ion classes.
2. **Edge-score vs node-score double-counting.** Java may sum BOTH the directional node score AND a separate edge-score along the path; Rust may sum only one.
3. **Charge-state ladder.** Java may sum scores across multiple charge states of the precursor or across multiple isotope offsets per split; Rust may take a single value.
4. **Per-segment vs per-split aggregation.** Java's RawScore may aggregate over all segments of the PrimitiveAaGraph, Rust may aggregate only over the path used by GF.

## Phase 3: Surface of the fix (filled in after Phase 2)

[To be populated by the code-comparison investigation. Likely files:
`rust/crates/scoring/src/scoring/scored_spectrum.rs`,
`rust/crates/scoring/src/scoring/psm_score.rs`,
or both.]

## Phase 4: Validation (filled in after Phase 3)

Regression-test scan=28787 must hit 297 ± tolerance. Then PXD001819 + Astral + TMT Percolator @ 1% FDR must hit the gates from the original spec.
