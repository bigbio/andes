# Chimeric two-pass cascade: optimized + multi-dataset gate summary (2026-05-31)

Branch `feat/chimeric-dda-plus`, clean HEAD `161d719e`. This note summarizes the
two-pass chimeric cascade, the speed optimizations landed this session, and the
same-machine entrapment-validated results vs Java MS-GF+ on all three benchmark
datasets. It is the reference for a future opt-in `--chimeric` PR.

## What the cascade is

Speed-correct chimeric search that recovers co-isolated second peptides WITHOUT the
wide-window cost (MaxQuant "second-peptide" model):

- **Pass 1** — narrow top-1 primary search per scan (the normal, fast path).
- **Pass 2** — MS1-gated targeted secondary search: detect co-isolated precursors in
  the MS1 isolation window via averagine-KL (`detect_coisolated`), then score a handful
  of candidates at each co-isolated mass on the **residual** spectrum (primary's matched
  charge-1 b/y peaks removed), with one single-bin GF SpecEValue DP per secondary
  (`search_secondary`). Secondaries are `force_push`ed as extra emissions, not
  competitors for the primary's top-1 slot.

Files: `crates/search/src/coisolation.rs` (core), driver `run_pass2_coisolation` in
`crates/search/src/match_engine.rs`, wiring in `crates/msgf-rust/src/bin/msgf-rust.rs`.
The entire path is gated on `--chimeric`; `--chimeric off` (default) is byte-identical
to the narrow engine.

## Speed optimizations landed this session (Astral 7:25 → 5:59)

| commit | optimization | effect |
|---|---|---|
| `3556bee8` | `search_secondary` uses a lightweight raw-peak primary match instead of a second `ScoredSpectrum::new` (drop 1 of 2 deconv builds per secondary) | Pass-2 ~14% cut, results bit-identical |
| `29d5e3d2` | build candidate index ONCE — `PreparedSearch::into_parts`/`from_parts` move the precursor-tolerance-independent 16.8M-candidate enumeration from the cal pre-pass into the main pass | main_prepare 14.2s → 0.00s |
| `3d941a02` | overlap `read_with_ms1` (~20s) on a background thread behind cal+enumerate | read_with_ms1(join) → 0.00s |

Earlier in the cascade work: single-bin GF per secondary (the original Pass-2 win,
9:42 → 7:22), `max_kl` 1.0 → 0.3 (FDP toward nominal). Profiling localized the Pass-2
floor: GF SpecEValue ~48% (NOT RawScore-gateable — passing PSMs reach RawScore −16),
score-loop 22%, res_ss build 15%, detect 14%.

## Same-machine A/B vs Java (entrapment-validated, FDRBench 1:1 shuffled-target)

| dataset | Rust @1% | Java @1% | PSMs | Rust wall | Java wall | speed | Rust entrapment FDP (combined) |
|---|---:|---:|:--:|---:|---:|:--:|:--:|
| **Astral** (LFQ DDA, HCD) | **55,581** | 35,818 | **+55.2%** | 5:59 | 5:52 | tied (−7s, < run noise) | 1.54% |
| **PXD001819** (UPS1 yeast) | **18,197** | 14,989 | **+21.4%** | 1:12.9 | 1:21.6 | **−8.7s faster** | 1.52% |
| **TMT** (PXD007683, CID) | 9,628 | 10,194 | −5.5% | 2:02.9 | 3:07.1 | **−64s faster (+34%)** | 0.80% |

All FDPs ~nominal (the chimeric gains are real co-isolated peptides, not coincidental
targets — both the reversed-decoy and entrapment rulers agree). Astral and PXD001819
beat Java on **both** axes; PXD001819 flipped from the prior −1.5% narrow-path deficit
to a +21.4% win.

## TMT: the lone remaining merge-gate blocker (deferred)

TMT has almost no co-isolation (CID, narrow isolation → Pass-2 = 2.66s), so the cascade
cannot help it. The −5.5% TMT PSM gap is REAL (entrapment: Java 9,224 vs Rust 8,436 real
PSMs) but is NOT a chimeric or GF-DP problem — it traces to a per-peptide CID
node-scoring divergence (Rust under-scores Java's winning peptides on CID spectra; 95%
of the 438 label-flips have Java's peptide outside Rust's RawScore top-10). Ruled out:
window/calibration (cal-off was worse), modifications (identical), top-1 selection
(`PsmMatch::Ord` already SpecE-first), aggregate T/D discrimination (Rust equal/better).
Full analysis: `2026-05-31-tmt-gap-diagnosis-not-gf-bug.md`. Deferred to a future
iteration — candidate strategies are additive Percolator features (e.g. DeltaRawScore)
or a per-ion CID node-scoring trace, not the chimeric cascade.

## Merge-gate status

[[merge-gate-beat-java]] requires Rust to beat Java on PSMs AND speed on all 3 datasets.
2/3 met (Astral, PXD001819). TMT PSMs (−5.5%) blocks a full merge. NOT merged — branch
parked as a clean, reviewable opt-in `--chimeric` feature pending the TMT strategy.
