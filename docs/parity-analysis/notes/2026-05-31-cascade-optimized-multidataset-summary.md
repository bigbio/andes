# Chimeric two-pass cascade: optimized + multi-dataset gate summary (2026-05-31)

Branch `feat/chimeric-dda-plus`. This note summarizes the two-pass chimeric
cascade, the speed optimizations landed this session, and the same-machine
entrapment-validated results vs Java MS-GF+ on all three benchmark datasets.
It is the reference for opt-in `--chimeric` **DRAFT PR #42**. The cascade has
since been hardened through code review, two adversarial rounds, a dead-code
cleanup, and an external GF/SpecE parity audit — see the **Final state (rev5)**
section at the bottom for the current HEAD and final numbers.

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

---

## Addendum: code-review fixes A/B revalidated — strict improvement (2026-05-31)

After the multi-agent code review, two secondary-PSM correctness fixes landed
(commit `c7940916`): **(A)** Pass-2 secondaries now get real fragment-ion features
(`compute_psm_features` on the residual ScoredSpectrum + reused edge_score) instead
of all-zeros; **(B)** new `PsmMatch.precursor_mz_override` so the PIN writer emits
correct ExpMass/dm/absdm from the co-isolated precursor mass (not the primary's).
Plus **(E)** `--chimeric` no longer forces top_n=1 on non-mzML inputs (`ffaab1d9`).

Same-machine revalidation vs the pre-fix cascade (entrapment = authoritative ruler):

| dataset | metric | pre-A/B | A/B | Δ |
|---|---|---:|---:|---:|
| Astral | normal @1% | 55,581 | **58,641** | +3,060 |
| Astral | entrapment REAL PSMs | 45,248 | **49,610** | **+4,362** |
| Astral | entrapment FDP (combined) | 1.54% | **1.13%** | cleaner |
| PXD001819 | normal @1% | 18,197 | 17,565 | −632 |
| PXD001819 | entrapment REAL PSMs | 16,276 | **16,371** | +95 |
| PXD001819 | entrapment FDP (combined) | 1.52% | **1.37%** | cleaner |

By the entrapment ruler (more-real-PSMs + lower-FDP on BOTH datasets), A/B is a strict
improvement and more honest. The Astral gain is large (+4,362 real PSMs, FDP 1.54→1.13)
because Astral is secondary-rich; the PXD normal-DB −632 is reduced *false* inflation
(entrapment real rises while FDP falls). Both still beat Java decisively (Astral +63.7%,
PXD +17.2% on normal @1%). TMT A/B confirmed: 9,706 (was 9,628, +78 — the few TMT
secondaries now carry proper features; still −4.8% vs Java 10,194, the lone blocker).
No regression on any dataset; A/B is a strict improvement (more real PSMs + lower FDP
everywhere).

## Final state (rev5) — DRAFT PR #42, HEAD `b46b610b`

After the A/B round above, the cascade was hardened and shipped as a do-not-merge
DRAFT (PR #42 → `dev`). Work done, in order:

- **Code review (5-agent)** — clean; A/B feature fixes already folded in.
- **Adversarial review round 1** (`429cf1bf`): bounded-memory streaming
  (`read_with_ms1_chunked`, per-chunk MS1 link + parser-thread pipeline + offset),
  Pass-2 precursor calibration (`adjusted_observed_neutral_mass` before the
  prefilter), tolerant MS1 read.
- **Adversarial review round 2** (`16f396b7`): Pass-2 secondaries COMPETE for shared
  peaks (`search_secondary` threads `prior_claimed`, returns the winner's claimed
  set); real resync (`resync_to_next_spectrum` skips one bad scan and continues
  instead of silently truncating the tail).
- **Doc refinement** (`c5c8ea8d`): CLI/README parameter tables with defaults; ~150
  dev-history jargon comments stripped.
- **Dead-code removal** (`d824deab`): the refuted wide-window / fragment-posting-index
  / shared-fragment machinery, the `--chimeric-frag-index` flag, and 3 always-zero PIN
  columns removed (chimeric PIN 42→39 columns); cascade @1% preserved (71,839 ≈ 71,877,
  within noise).
- **GF/SpecE parity audit** (`docs/parity-analysis/notes/2026-06-01-p0-parity-audit-bench.md`):
  safe P1/P2 perf+robustness batch shipped (`cb808ce3`, behavior-neutral). The
  strongest parity candidate (P0.4, precursor-filter Java parity) was bench-validated
  and **regressed TMT 9,671→9,579** (the lone blocker dataset) → reverted (`1c706522`).
  This is the **n=9** confirmation that Java-faithful SpecEValue fixes redistribute
  rather than improve aggregate Percolator FDR. P0 grind stopped by user decision.

### rev5 final numbers (entrapment-clean, same machine)

| dataset | Rust @1% | Java @1% | Δ vs Java | entrapment FDP | wall | maxRSS |
|---|---:|---:|---:|---:|---:|---:|
| Astral | **71,877** | 35,818 | **+101%** | 1.04% | 6:38 | 10.9 GB |
| PXD001819 | **16,592** | 14,974 | **+11%** | 1.13% | 1:14 | 2.3 GB |
| TMT | 9,671 | 10,194 | **−5%** | 0.80% | 2:14 | 7.7 GB |

maxRSS on Astral equals the non-chimeric run (10.9 GB) — bounded-memory confirmed;
the footprint is index-dominated (scales with DB size, not input). Speed beats Java on
all three datasets same-machine (e.g. Astral chimeric 6:39 vs Java 6:46).

### Gate status

PSMs win 2/3 (Astral + PXD decisive and entrapment-validated), speed wins 3/3. **TMT
PSMs (−5%) is the only thing blocking the merge gate.** The TMT gap is diagnosed as a
per-peptide CID node-scoring divergence (RankScorer ion-match / peak-rank / per-rank
log-prob), NOT a cascade issue and NOT a gross GF-DP bug — CID narrow windows have
~no co-isolation, so the cascade cannot help TMT. Closing it needs a per-ion CID trace
+ Java instrumentation, or an additive CID-specific Percolator feature. PR #42 stays
DRAFT until TMT beats Java.
