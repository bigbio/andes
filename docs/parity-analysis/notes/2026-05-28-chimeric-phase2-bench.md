# Chimeric Phase 2 bench — additive MS1 isotope-KL feature does NOT control FDR

**Date:** 2026-05-28
**Branch:** `feat/chimeric-dda-plus` (rebased on dev incl. PR #40 Phase-1b cal)
**Binary:** chimeric-build with Phase 2 (averagine envelope + MS1 link + PrecursorIsotopeKL/SNR additive PIN columns), cal=auto, --chimeric.

## Result (PSMs @1% FDR, Percolator --only-psms)

| Dataset | off baseline | Phase-1 on | Phase-2 on (+isotope-KL) | Phase-2 vs off |
|---|---:|---:|---:|---:|
| PXD001819 | 14,808 | 17,015 | 17,749 | +20% |
| TMT       | 9,605  | 9,608  | 9,473  | -1.4% |
| Astral    | 36,715 | 71,347 | 72,457 | +97% |

(maxRSS: PXD 2.8 GB, TMT 8.4 GB, Astral 11.9 GB — the batch MS1-load path; fit in 27 GB.)

## Conclusion: additive KL feature is INSUFFICIENT — hypothesis refuted

Phase 2's premise was that an MS1 targeted-XIC isotope-KL feature, fed to
Percolator, would let the rescorer reject the spurious co-isolated PSMs that
Phase 1's multi-PSM emission over-counts. **It does not.** Astral (narrow
isolation windows → minimal real co-isolation) stays at ~+97% (72,457 vs 36,715)
— biologically implausible, the same inflation as Phase 1 (it even rose +1,110).
PXD rose +734, TMT fell -135. The KL column changed Percolator's ranking
slightly but did not shrink the FDR computed over the inflated
multi-PSM-per-scan set.

**Root cause:** PSM-level target-decoy FDR over ~5 PSMs/scan is inflated
*structurally* — the decoy model doesn't encode the "few real peptides per scan"
constraint. A soft Percolator feature cannot fix that; it can only re-rank within
the already-inflated set (and may even overfit, as the small increases show).

## Implication

- Chimeric counts remain **untrustworthy** through the additive-feature path.
  PXD's "17,749 > Java 14,974" is NOT a real win — the Astral control proves the
  method inflates.
- **Chimeric does NOT pass the merge gate** ([[merge-gate-beat-java]]): it does
  not deliver a trustworthy PXD/TMT-beats-Java gain.

## What would actually be required (not yet done)

1. **Hard pre-FDR filter**, not a soft feature: discard PSMs whose precursor
   isotope envelope is absent/poor in MS1 (high KL) BEFORE Percolator, so the
   FDR is computed over a credible set. (The Phase 2 plan listed this as
   "optional" — the bench shows it is mandatory.)
2. **Per-scan / peptide-level FDR** or rank-1-vs-rest competition that models the
   multi-PSM-per-scan structure (closer to how MSFragger-DDA+ + Philosopher
   actually control FDR).
3. **Greedy shared-fragment rescoring (Phase 3)** to down-weight peptides
   explained only by another's fragments.

Until at least (1)+(2) land and the Astral control returns to a plausible value,
chimeric is not mergeable. Decision: do NOT pursue more additive features;
either implement the hard-filter + per-scan FDR, or shelve chimeric.
