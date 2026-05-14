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

## Phase 2.5: Static-comparison findings (code-explorer pass)

Side-by-side reading of Java and Rust score chains identified three concrete divergences. None of them on its own conclusively explains the 3× gap, but Divergence B is the strongest single-cause candidate.

### Rust function chain (pin column 7)
- `match_engine.rs:277-284` → `score_psm(scored_spec, peptide, scorer, z, tol) + compute_cleavage_credit(...)`
- `psm_score.rs:29-84` `score_psm` — split loop, accumulates `prefix_mass_acc` in f64, computes `prefix_nominal = nominal_from(prefix_mass_acc)` (round of sum)
- `scored_spectrum.rs:96-222` `ScoredSpectrum::new` — builds `prefix_score_cache`/`suffix_score_cache`
- `scored_spectrum.rs:510-558` `directional_node_score_inner` — per segment, per ion-type, lookup nearest-peak rank
- `pin.rs:347` `let raw_score = psm.score.round() as i32;`

### Java function chain (pin column 7)
- `DirectPinWriter.java:213` → `match.getScore()`
- `DBScanner.java:519-541` → `rawScore = scorer.getScore(...)` then `score = cleavageScore + rawScore`
- `FastScorer.java:59-76` `getScore` — split loop, reads `nominalPrefixMassArr[i]` (sum of rounds), accumulates `Math.round(prefixScore[m] + suffixScore[m])`
- `FastScorer.java:20-35` constructor — builds `prefixScore[m]`/`suffixScore[m]` via `scoredSpec.getNodeScore(NominalMass(m), isPrefix)`
- `NewScoredSpectrum.java:134-166` `getNodeScore` — per ion-type, computes theoMass in **f32**, finds `spec.getPeakByMass(theoMass, mme)`, gets rank
- `NewScoredSpectrum.java:43-56` constructor — `filterPrecursorPeaks()` sets intensity=0 but keeps the peak; `setRanksOfPeaks()` ranks **all** peaks including zero-intensity

### Divergence A — Nominal mass accumulation
- Java: `nominalPRM[i] = nominalPRM[i-1] + round(exactMass_k * SCALER)` — **sum of rounds**, one round per residue
- Rust `psm_score.rs:47-52`: `prefix_nominal = nominal_from(prefix_mass_acc)` — **round of sum**, accumulated in f64 before a single round
- Effect: ±1 Da index shift on some split positions → wrong cache bucket. Drift, not 3× gap.

### Divergence B — Precursor peak filtering (likely root cause)
- Java `NewScoredSpectrum.java:43-56`: precursor-window peaks get intensity=0 but **retain valid ranks**; fragment ions landing near precursor get an **observed-rank score**.
- Rust `scored_spectrum.rs:264-299, 656-677`: precursor-window peaks are tagged `ranks[i] = u32::MAX`; `nearest_peak_rank_in` skips them → fragment ions near precursor get **missing-ion log score** (typically much lower).
- Effect: every fragment ion m/z overlapping the precursor isolation window in the Rust path scores as missing where Java scores as observed. For a charge=2 precursor at m/z ~1027 + 17-residue peptide, several mid- and high-mass fragment ions (especially y2+/b2+ near precursor) will be affected. Each missed ion costs 5–20 log units. Across all ion types × splits, this could easily account for 100+ score units, consistent with a 297→108 gap.

### Divergence C — `IonType.mz()` float precision
- Java f32, Rust f64. Sub-Da differences in theoMass. Rarely changes which peak is selected. Negligible.

## Phase 2.6: Divergence-B hypothesis test — no-op

Modified Rust to rank ALL peaks (including precursor-filtered) so filtered peaks land at worst-tier rank instead of `u32::MAX` skip → score-as-missing. Re-ran the oracle on scan=28787.

Result: RawScore=**108 unchanged** (oracle exit 1).

Cause: `HCD_QExactive_Tryp.param` (the param used for PXD001819) contains **only `reduced_charge=0` entries** for charge=2: offsets {-1.0005, 0.0, +1.0005} with tol=Da(0.5). The filter formula `filter_mz = (neutral_mass + c·PROTON)/c + offset` collapses with `c = charge - reduced_charge = 2`. For a charge=2 precursor this puts filter zones at the precursor m/z (~1027) ± 1 Da only — **no charge-reduced zones in mid-mass fragment territory**. So almost no fragment ions of K.IVNEEFDQLEEDTPVYK.L land in those zones, and the rank-vs-missing distinction never fires for this PSM.

Divergence B remains a real Java↔Rust semantics gap that should be fixed eventually (it would matter for ETD/ECD param files with `reduced_charge=1` entries), but it is NOT the cause of the PXD001819 RawScore gap. Change was reverted; lib tests + `gf_java_parity` were green throughout.

## Phase 3: Confirm root cause with per-split / per-ion instrumentation (NEXT)

Static analysis identified candidates A and B; B is now ruled out for this PSM. Divergence A (sum-of-rounds vs round-of-sums for nominal accumulation) is plausible but cannot account for a systematic 3× gap — it only causes ±1 Da index jitter per split.

The remaining hypotheses are at the ion-iteration / scoring-table level:
- Ion-type set mismatch: Rust's `partition_ion_logs(part)` may return fewer ion types than Java's `scorer.getIonTypes(charge, parentMass, seg)`.
- Log-table value mismatch: Same (ion, rank) keys may have different log values in Rust vs Java because of param-parsing drift.
- Per-segment partition mismatch: `Partition` keys may differ between Java and Rust for the same (charge, parentMass, seg).
- Edge-score double-count or omission: pin RawScore in Java may include both node and edge contributions; Rust may include only one.

Before changing code, instrument both implementations and dump per-ion-match tuples for scan=28787 + K.IVNEEFDQLEEDTPVYK.L + charge=2. The first divergence in the dump identifies the root cause.

### Instrumentation strategy

Rust: add a feature-gated `eprintln!` in `directional_node_score_inner` printing `(seg_idx, ion, theo_mz, observed_rank_or_none, log_score_added)`.
Java: project memory mentions `-Dmsgfplus.trace=true -Dmsgfplus.trace.scan=N -Dmsgfplus.trace.pep=SEQUENCE` infrastructure is already wired. If per-ion contributions are not in the existing trace output, extend `NewScoredSpectrum.getNodeScore` to emit them when the trace flags are set.

Compare the two dumps for the same PSM. Outcomes:
- Same ions iterated, different ranks → rank-assignment divergence (Divergence A and/or peak-lookup bug).
- Same ions iterated, same ranks, different log scores → param-table parse divergence.
- Different ions iterated → ion-type set / partition divergence (most likely if scale is ~3×).
- Same ions, same ranks, same scores, but Java has additional contributions outside this loop → edge-score or out-of-graph contribution.

### Java instrumentation
Project memory and the parity-analysis history indicate Java already has `msgf-trace` wiring via `-Dmsgfplus.trace=true -Dmsgfplus.trace.scan=N -Dmsgfplus.trace.pep=SEQUENCE`. Confirm what gets emitted; if per-split `(split, prefixMass, suffixMass, prefixScore, suffixScore, contribution)` is not already in the trace, add it temporarily.

### Rust instrumentation
The msgf-rust workspace has a `msgf-trace` binary at `rust/crates/cli/`. Add a feature-gated `eprintln!` block in `psm_score.rs::score_psm`'s split loop emitting the same fields. Run on scan=28787 + IVNEEFDQLEEDTPVYK + charge=2.

### Expected output of the diff
First divergent split position identifies the root cause:
- If both sides have the same `prefix_nominal` and `suffix_nominal` but different `prefix_score`/`suffix_score` → Divergence B (the rank-vs-missing path).
- If the `nominal`s diverge first → Divergence A.
- If both agree on split contribution but a charge/segment-level multiplier differs → a fourth divergence not in the static review.

## Phase 4: Validation (after Phase 3 fix)

Regression test scan=28787 must hit 297 ± 10 tolerance. Then PXD001819 + Astral + TMT Percolator @ 1% FDR must hit the gates from the original spec (≥14,800 / ≥33,000 / ≥10,500).
