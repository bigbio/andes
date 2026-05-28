# ID-rate improvement for PXD001819 + TMT — design

**Date:** 2026-05-28
**Branch:** `feat/id-rate-pxd001819-tmt`
**Type:** Phased, bench-gated investigation (outcome of each phase decides the next)

## Goal

Increase peptide-spectrum-match (PSM) counts at 1% FDR (Percolator) for the two
datasets where msgf-rust currently trails or barely matches upstream Java
MS-GF+ v2024.03.26:

| Dataset | Current Rust @1% | Upstream Java @1% | Stretch target (+10% over Rust) |
|---|---:|---:|---:|
| PXD001819 (UPS1 yeast tryptic) | 14,755 | 14,974 | ~16,230 |
| TMT (a05058, PXD007683)        | 9,605  | 10,115 | ~10,565 |

The +10%-over-Rust figure means *beating* upstream Java by ~8% (PXD) / ~4% (TMT)
— the same kind of edge msgf-rust already holds on Astral DDA (+9.8%: 36,715 vs
33,425). It is a **stretch target**, not a hard gate (see "Success criteria").

Constraints:
- **No wall-time regression beyond noise (~3%)** on any of the three datasets.
- **No PSM regression on Astral** (currently 36,715 @1% FDR — already beating Java).

## Background

The I5 trace investigation (`docs/parity-analysis/notes/2026-05-26-score-psm-trace-findings.md`,
merged via PR #37) localized the Rust↔Java per-PSM scoring divergence on
PXD001819. The gap is dominated by **label flips**: Rust ranks a different
(usually wrong) top-1 peptide than Java because its scoring over-scores weaker
peptides by +5..+20 RawScore and under-scores Java's favored peptide by ~13.
Three roughly-equal contributors, all in the hot scoring path:

- **H1 — ion-type enumeration** (~10% of divergent ions): Rust scores some ions Java doesn't.
- **H2 — peak-rank assignment** (~40%): when several peaks fall inside the
  fragment-tolerance window, the selected peak / assigned rank differs from Java's.
- **H3 — per-rank log-probability lookups** (~40%): a different log-prob value is
  returned even when Rust and Java agree on the rank.

The project's empirical audit pattern (n≈9+, recorded in workspace memory):
- **Additive** changes (new PIN columns) are *safe* but historically *flat* —
  Percolator already extracts most of that signal via correlated columns.
- **Modifying-existing-distribution** changes usually *regress* Percolator @1% FDR
  even when individually "more correct" — UNLESS they move Rust's top-1 selection
  *toward Java's*. The three biggest historical gains were all top-1-restoring
  changes: iter20 high-res tolerance fix (+4,650 Astral), iter33 edge-in-ranking
  (+3,705), iter29 main-ion fix (+379).

This design exploits that pattern: pursue the cheap top-1-restoring config lever
first, then the hot-path label-flip fix (also top-1-restoring), keeping additive
features as a zero-risk safety net.

## Relevant code (from the 2026-05-28 code map)

- Peak rank assignment (H2): `crates/scoring/src/scoring/scored_spectrum.rs`
  - `ScoredSpectrum::new` (~L168–357): assigns ranks intensity-desc then m/z-asc;
    precursor-filtered peaks excluded before ranking.
  - `nearest_peak_rank_in` (~L897–918): binary search + linear scan in tolerance
    window, selects **highest intensity** (strict `>`, so first peak wins on ties).
  - `directional_node_score_inner` (~L715–750): tolerance is `mme.as_da(theo_mz)`.
- Log-prob tables (H3): `crates/scoring/src/scoring/rank_scorer.rs`
  - `node_score` / hot-path indexing: `idx = min(rank, max_rank).max(1) - 1`;
    missing ion uses the `max_rank` sentinel slot.
  - `error_score`, `ion_existence_score`.
- Ion-type enumeration (H1): `crates/scoring/src/param_model.rs`
  - `build_partition_ion_types_cache`, `partition_for` / `find_partition`.
- RawScore assembly: `crates/scoring/src/scoring/psm_score.rs::score_psm`;
  `pin_score` (= node + cleavage) vs `rank_score` (= node + cleavage + edge) split
  in `crates/search/src/match_engine.rs` (~L408–448).
- Instrument/tolerance resolution: `crates/msgf-rust/src/bin/msgf-rust.rs` (~L585–631);
  `detect_instrument_type` in `crates/input/src/mzml.rs`;
  `InstrumentType::is_high_resolution()` gates the 20-ppm vs 0.5-Da feature tolerance
  in `compute_psm_features` (`crates/search/src/match_engine.rs`).
- PIN columns: `crates/output/src/pin.rs::write_header` / `write_psm_row`.
- Feature struct: `crates/search/src/psm.rs::PsmFeatures`.

## Success criteria

- **Hard gate (per change):** ships only if it gains PXD or TMT @1% FDR, regresses
  *none* of PXD/TMT/Astral beyond Percolator noise (~±0.3%), and keeps wall within
  ~3% on all three. Otherwise revert.
- **Stretch target:** +10% over current Rust on PXD and TMT. Treated as a
  *direction we bench toward*, NOT a revert-everything gate — we ship every change
  that net-gains under the hard gate, and report the cumulative result even if it
  lands below +10%.
- **Speed:** Phase-2 hot-path changes must show no wall regression beyond noise;
  revert if they do, even if PSMs gain.

## Bench protocol (every phase)

Reuse the established VM harness (`/srv/data/msgf-bench/`):
- Build Rust release with `target-cpu=sandybridge` (the committed `.cargo/config.toml`).
- Run all 3 datasets, cal=auto, `--top-n 1 --threads 8`, per-dataset
  tolerance/instrument/protocol matching the README bench.
- Percolator 3.7.1 via `run_percolator_docker.sh` (`--seed 42 --only-psms`),
  parse 1%-FDR target count.
- Compare against the locked baseline (PXD 14,755 / TMT 9,605 / Astral 36,715).
- One commit per change; revert in place if the hard gate fails. Keep reverts on
  the branch as record (matches the project's iteration-shipping model).

## Phases

### Phase 0 — Diagnostic (measurement only; no code change)

Determine, for each of PXD001819 and TMT, what the *current* run actually resolves
to. Output a table:

| Dataset | Resolved instrument | Feature tolerance | Scoring tolerance | Calibration fired? | Dominant activation |
|---|---|---|---|---|---|

Sources:
- The Rust run's stderr ("Param resolver: auto-detected ...", "instrument = ...").
- `InstrumentType::is_high_resolution()` for the resolved instrument → 20 ppm vs 0.5 Da.
- The MassCalibrator log line (fired vs skipped, and the <confident-PSM guard count).
- The mzML `<analyzer>` cvParam (FTMS/Orbitrap/ITMS) to know the *true* instrument.

**Decision:** if PXD/TMT are high-res but resolving to the 0.5-Da low-res tolerance,
Phase 1a has the iter20 win available. If they're already on the high-res path,
Phase 1a yields ~0 and we proceed to Phase 1b / Phase 2.

### Phase 1 — Config levers (low speed risk)

**1a — Instrument / fragment tolerance.** If Phase 0 shows a low-res/high-res
mismatch on the *tolerance* path, fix instrument resolution so the 20-ppm
high-res branch engages in both `compute_psm_features` and the scoring tolerance.
Likely a fix to how `detect_instrument_type` feeds `InstrumentType::is_high_resolution()`,
or a default change. Bench-gate.

**1b — Calibration.** If Phase 0 shows `--precursor-cal auto` skipping on PXD/TMT
(the historical guard skips when <200 confident PSMs), investigate the guard
threshold and whether calibration would add PSMs (it did on Astral). Bench-gate.

### Phase 2 — Label-flip gap (hot-path; bench-gated, revert-ready)

**2a — H2 peak-rank assignment.** Using the I5 artifacts (scan 41522,
`R.DPANLPWASLNIDIAIDSTGVFK.E`), find the first RANK_DIFF ion (theo_mz, rust_rank,
java_rank). Trace both peak-selection paths:
- Rust `nearest_peak_rank_in`: intensity-max in window, strict `>`.
- Java `getPeakByMass`: identify its tie-break (first peak / closest m/z / index order).
Make Rust match Java's rule. Re-run the I5 trace harness; confirm RANK_DIFF count
drops and (by the H2→H3 coupling) LOGPROB_DIFF drops proportionally. Bench-gate on
PXD + TMT + Astral.

**2b — H3 log-prob indexing.** Only if 2a lands and the bench gate passes. Close
the residual same-rank/different-value cases (307 of 608 in the I5 data). Likely a
table-indexing or clamping difference in `rank_scorer.rs`. Bench-gate.

H1 (ion enumeration) is deferred unless 2a/2b leave a clear ion-enumeration residual
in the re-run trace.

### Phase 3 — Additive PIN features (safety net)

Only if Phases 1–2 fall short of the stretch target. Add new PIN columns (no change
to existing column values):
- `MeanMatchedRank` — average `nearest_peak_rank` over matched b/y ions.
- `ScoreFractionTop1Split` — max single-split contribution / RawScore (requires
  keeping the max split alongside the sum in `score_psm`).
- `ln(num_distinct)` — already computed for E-value; expose as a column.
Bench-gate; expected flat-to-small.

## Risks & mitigations

- **Phase 2 regresses Percolator** (n=9+ audit): bench-gate per dataset, revert in
  place. The I5 trace re-run is the leading indicator (RANK_DIFF should drop) before
  trusting the Percolator delta.
- **Phase 2 slows the hot path:** measure wall on all 3 datasets; revert if >3%.
- **+10% unreachable via scoring alone:** acknowledged up front; success criteria
  treat +10% as a direction, not a revert-all gate. Worst case we close the
  gap-to-Java and modestly beat it, and ship the net-positive subset.
- **Astral regression:** every change is checked against Astral, which already
  wins — a change that helps PXD/TMT but breaks Astral is reverted.

## Out of scope

- Astral (already beats Java; only a regression guard here).
- Algorithmic rewrites of the GF DP or candidate enumeration.
- New scoring models / `.param` retraining.
- The README bench-table PR (#39, separate).
