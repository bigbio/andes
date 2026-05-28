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

## Deeper investigation (2026-05-28, "is it an implementation problem?")

Investigated whether the inflation is an implementation/faithfulness gap. Findings:

1. **Real bug: top-N ran at 10, not 5.** The `--chimeric` top-N bump only fires
   when `cli.top_n == 1`, but the CLI default is 10, so it never triggered — the
   bench emitted up to 10 PSMs/scan. Over-emission, but NOT the main driver.

2. **Hard MS1 isotope filter is INSUFFICIENT (decisive test).** Post-filtered the
   Astral chimeric PIN to keep only rows with a real precursor envelope
   (PrecursorIsotopeKL < 1; dropped the 32% no-envelope KL=10 rows + 33% poor
   KL≥1 rows) and re-ran Percolator: **72,457 → 69,250** @1% FDR — still +89% over
   the 36,715 off baseline. The isotope filter (soft OR hard) does NOT deflate.

3. **Root cause — precursor-presence ≠ fragment-origin.** The isotope envelope
   only confirms a precursor EXISTS at that mass in MS1. On dense MS1 (Astral),
   almost every wide-window candidate has *some* real co-eluting precursor near
   its mass, so its envelope passes — even though the MS2 fragments did not come
   from it. The discriminator between true and spurious co-fragmentation is
   FRAGMENT-level (shared-fragment competition), which neither Phase 1 nor Phase 2
   provides.

4. **Wide search inflates scans-passing too.** 58,314 unique scans yield a passing
   target (vs 36,715 off) — searching a wider window with top-N gives each scan
   many more chances; PSM-level FDR doesn't catch this multiple-testing inflation.

## Verdict

Implementation gaps exist (top-N bug; soft-not-hard filter; no shared-fragment
rescoring), but the specific Phase-2 fix (MS1 isotope KL) is empirically the WRONG
lever — proven by the hard-filter test. A trustworthy chimeric search needs:
(a) **fragment-level shared-fragment competition** (Phase 3) as the primary
discriminator, AND (b) a **co-isolation-gated emission** (only report extra
peptides when MS1 shows a distinct co-isolated precursor) + proper per-scan/
peptide FDR — NOT top-N-always. That is a substantial build, and the gain is real
only on genuinely co-isolated (wide-window) data (PXD-like), not narrow-window
(Astral/TMT). Recommendation: shelve the isotope-filter path; revisit chimeric
only via fragment-competition if pursued at all.
