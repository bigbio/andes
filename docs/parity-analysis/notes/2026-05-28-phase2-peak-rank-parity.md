# Phase 2 finding: peak-rank assignment is bit-identical to Java (H2 is NULL)

**Date:** 2026-05-28
**Branch:** `feat/id-rate-pxd001819-tmt`
**Dataset:** PXD001819, scan 41522 (a documented I5 label-flip scan)
**Tooling:** `msgf-trace --dump-peaks` (commit `12c1839c`) — dumps the post-filter,
post-deconvolution active peak list (rank, m/z, intensity).

## Headline

**Rust and Java assign the identical rank to every peak.** Comparing Rust's
dumped active peak list against Java's per-ion trace ranks, matched by *observed
peak m/z* (not theoretical-ion m/z), the rank offset is **+0 for all 465 matched
peaks**. No exceptions.

This **debunks the I5 doc's central hypothesis** (`2026-05-26-score-psm-trace-findings.md`),
which attributed ~40% of the scoring divergence to H2 (peak-rank assignment) and
made it the primary fix target for this whole investigation.

## What went wrong in the original I5 analysis

The I5 `analyze.py` (and an even cruder ad-hoc script) aligned Rust↔Java ions by
`(ion_kind, round(theo_mz/1e-3))` — i.e. by *theoretical* ion m/z. The same
theoretical m/z can correspond to **different physical ions** in the two peptides
being compared (e.g. Rust `y/1` at theo_mz X vs Java `y/2+offset` at a
coincidentally-equal theo_mz X). Those spurious cross-matches produced the
RANK_DIFF=301 count and the apparent "+2 rank offset."

When you instead match **actual observed peaks by their m/z** (what
`--dump-peaks` enables), the ranks are identical. The peak ranker
(`ScoredSpectrum::new` intensity-desc + m/z-asc, precursor-filtered peaks excluded)
already matches Java's `Spectrum.setRanksOfPeaks()` (intensity-desc via
`IntensityComparator`, precursor peaks zeroed-but-ranked-at-bottom) for every peak
that matters.

Confirmed details:
- Java zeroes precursor peaks but keeps them in the ranked list (sorted to the
  bottom at intensity 0); Rust removes them. This makes Java's total peak count
  ~3 higher (max rank 489 vs Rust active 486) but does **not** shift any real
  peak's rank, because the zeroed peaks sort below every real peak either way.
- Rust's precursor-filter tolerance (`pof.tolerance`) already equals Java's `mme`
  (0.5 Da for CID_LowRes) for these offsets — a test swapping to `param.mme` left
  the active-peak count unchanged at 486 (reverted as a no-op).

## The label-flip still exists, but is not a rank bug

On current master, scan 41522's Rust top-1 is a decoy (`VVYGNIYEIEIDRLFLTDQR`,
score 13, SpecE 4.66e-4) — Java picks a real peptide. But:
- The scores are tiny (top non-decoy is 11), so this is a low-quality spectrum
  near the noise floor.
- Since peak ranks are identical, the divergence must be in **candidate
  enumeration (H1)** — whether Rust enumerates/locates the Java-favored peptide in
  this scan's mass window — and/or **log-prob table values (H3)**, not rank
  assignment.

## Implication for the +10% goal

- The single lever the whole Phase 2 plan was built around (fix peak-rank
  assignment) **does not exist** — Rust already matches Java bit-for-bit there.
- Combined with the external reviewer's "BSA 217/217 top-1 parity" and Phase 0's
  finding that instrument/tolerance/calibration are already correct, the picture
  is: **Rust's scoring is at or near parity with Java.** The residual gap
  (PXD001819 −1.1% after Phase 1b's +53; TMT −5%) is small, concentrated on
  low-quality / SpecE-tail spectra, and attributable to H1/H3 micro-differences
  that the n=9+ audit says regress Percolator when "fixed" individually.
- **+10% over current Rust is not reachable via scoring-path fixes** — the scoring
  path is essentially correct. The shippable result from this investigation is
  Phase 1b (+53 PXD, zero-risk calibration nudge).

## Recommended next directions (if ID-rate work continues)

1. **Stop chasing rank/score parity** — it's already there.
2. If a large ID gain is still wanted, it likely requires *algorithmic* changes
   (e.g. a different candidate-generation / scoring model), not parity fixes — a
   research project, not a bench-gated tweak.
3. H1 (does Rust enumerate the Java-favored peptide in-window?) is the only
   remaining parity question worth a *cheap* check, but the per-scan evidence
   suggests these are low-score spectra that don't cross 1% FDR anyway.
