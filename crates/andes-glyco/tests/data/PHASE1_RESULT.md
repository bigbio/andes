# Phase-1 glyco backbone-solver gate — RESULT

**Verdict: FAIL** (searchable-overall 59.8% < 70% gate; sparse stratum 2.7% < 55% target).
Significant improvement over prior run (37.3% → 59.8%); ranking fix recovered +22.5 pp overall.
Sparse stratum remains the hard limiter — see failure analysis below.

## The three numbers (n = 542 confident glyco PSMs)

| Metric | Value | Prior (pre-ranking-fix) |
|---|---|---|
| Oxonium-fire rate | **100.0%** (542/542) | 100.0% |
| **Searchable-backbone OVERALL** | **59.8%** (324/542) | 37.3% (202/542) |
| **Searchable-backbone SPARSE stratum** (≤2 true core-Y rungs, n=148) | **2.7%** (4/148) | 0.0% (0/148) |

Gate: PASS requires overall ≥ 0.70 AND sparse ≥ ~0.55. **Both fail.**

## Run provenance

- Binary commit: see commit on branch `glyco-phase1` for this PHASE1_RESULT.md update.
- Branch: `glyco-phase1`.
- Gate window (non-negotiable): symmetric `in_window(solved, truth) = |solved - truth| <= max(truth*20e-6, 0.01)`.
- Solver call: `solve_backbone(&peaks, precursor_neutral, prec_z, 20.0, 5)` (tol 20 ppm, top_k = 5).
- Oxonium gate: `oxonium_gate(&peaks, 0.10, 20.0)`.
- All 8 unit tests green.

## Ranking improvements applied (this run)

All changes are in `crates/andes-glyco/src/backbone.rs`:

1. **Intensity-weighted votes** (primary lever): each (peak, charge, rung) vote contributes
   `sqrt(intensity / base_peak_intensity)` instead of flat +1. Per-rung best across charges
   is tracked via a two-pass `rung_best: HashMap<(bin_key, rung_idx), (best_w, mass_sum, vote_count)>`.
   The intensity_score of each candidate = sum of per-rung best-weights (charge-deduplicated,
   so a peak contributing to the same rung at z=2 and z=3 only counts once).

2. **Rung-specific weights**: Y0 (bare backbone) and Y1 (backbone+HexNAc) each multiplied by 2.0×,
   Y2 by 1.5×, Y3-Y5 by 1.0×. Y0+Y1 are more diagnostic (specific to the exact backbone mass)
   than Y3-Y5 (hexose additions common to many glycan structures).

3. **Sort order**: PRIMARY = core_y_hits DESC, SECONDARY = intensity_score DESC (rung-weighted,
   charge-deduplicated), TERTIARY = backbone_mass ASC (deterministic tiebreak).
   Old order was (core_y_hits DESC, flat_votes DESC, backbone_mass ASC).

4. **Minimum backbone mass** raised from 200 Da to 500 Da (MIN_BB): eliminates very common
   spurious cluster at ~203 Da (HexNAc oxonium m/z 204 at z=1 voting as Y0).

5. **Minimum glycan mass** added (MIN_GLYCAN = 406 Da): backbone candidates where
   precursor_neutral − backbone < 406 Da (< 2×HexNAc core) are rejected as implausible.

6. **New unit test**: `solve_backbone_intensity_weighted_prefers_bright_ladder` — verifies
   that a high-intensity true Y-ion ladder ranks above a dim spurious cluster with equal rung count.

7. **BackboneCandidate struct**: added `intensity_score: f64` field (public); probe and all
   existing tests updated accordingly.

## Top_k=500 ceiling (current algorithm)

Running `solve_backbone(.., 500)` (all candidates) on this dataset: **82.7% overall,
38.5% sparse**. The true backbone IS present in the candidate pool at rank >5 for 22.9%
of spectra — those spectra are where the ranking fix still fails to promote truth into top-5.

## Failure mode analysis

**Three failure classes** (among the 218 non-searchable spectra):

### Class A: Chimeric spectra / genuine competing signal (~40% of failures)
Spectra where a co-isolated precursor's Y-ladder completely dominates the MS2 signal.
The solver correctly ranks the dominant backbone (another peptide), but the truth file
assigns to the weaker signal. Not fixable by ranking changes without MS1 isolation info.
Example: scan 4014 (bb_true=1638.811, rungs=6 present, but top=1988.925 with h=6, v=55).
The 1988 Da backbone has 6 Y-ion hits from genuinely bright peaks — a second co-isolated
glycopeptide at 1988 Da dominates.

### Class B: Deamidation variants (+0.984 / +1.968 Da) (~20% of failures)
Spurious backbone at truth+0.984 Da or truth+1.968 Da. The solver finds the deamidated form
(Asn→Asp at the glycosylation site N-X-S/T) as the better-supported candidate because the
in-source deamidated peptide generates brighter Y-ions. Not fixable without modification-aware
hypothesis generation.
Examples: scan 7242 (nearest=753ppm ≈ +0.991 Da), scan 7735 (nearest=328ppm ≈ +1.972 Da).

### Class C: Low-rung (sparse) spectra (~40% of failures, all of sparse stratum)
Spectra where the true backbone has ≤2 core-Y rungs present. The solver generates a h=2 candidate
for the true backbone but a spurious h=6 cluster from coincidental multi-charge peak assignments
always ranks higher (core_y_hits PRIMARY). With ≤2 true rungs, the true backbone can never win
on core_y_hits. To fix this, the sort order would need to demote core_y_hits — but doing so
regresses the overall rate (intensity-first sort tested: 53.9% overall, worse than current 59.8%).

## Implication for the gate

The 70% / 55% gate is not achievable with the current voting+ranking architecture at top_k=5:
- Overall ceiling with all improvements at top_k=500: 82.7% → structural cap on ranking performance.
- Sparse ceiling at top_k=500: 38.5% → below even the 55% sparse gate; the sparse stratum is
  fundamentally limited regardless of ranking, because low-rung spectra (≤2 true rungs) cannot
  generate candidates competitive with spurious h≥3 clusters.

### Next steps to cross the 70% overall gate

1. **Complement-pair confirmation**: require that a peak at approx `precursor_neutral − B + PROTON`
   (the glycan complement Y-ion) is present before accepting backbone B as a candidate. This would
   eliminate many chimeric backbone candidates.
2. **Glycan mass grid filter**: check if `precursor_neutral − backbone_mass` snaps to a known
   glycan mass (within 5 ppm) from a pre-loaded composition table. This is a database-assisted
   filter, not purely signal-processing.
3. **Precursor isotope deconvolution**: correct the precursor_neutral before backbone voting to
   eliminate the +0.984/+1.968 Da deamidation-variant spurious clusters.
4. **Sparse stratum specifically**: allow top_k=20 with a second-pass filter that checks if the
   backbone mass is biochemically plausible given the observed glycan delta — this would help
   sparse spectra where the true h=2 candidate would be found at rank 6-20.

## Truth source + cutoff (unchanged from prior run)

- **Dataset PXD025455** (PRIDE 2021/05), file `HCC_pool_Late_Fc3_r1`.
- **Ground truth = the reference engine-3.2 labile-search-mode NGLYCAN pepXML**.
- **backbone_mass = `calc_neutral_pep_mass`** (Cam-C 57.0215, ox-M 15.9949, N-term acetyl 42.0106).
- **Confidence cutoff**: rank-1 target PSMs, glycan delta > 200 Da, delta snaps to known composition
  within 0.05 Da, the reference engine expect ≤ 1.0. → **542 scans**.
