# Phase-1 glyco backbone-solver gate — RESULT

**Verdict: FAIL** (searchable-overall 37.3% < 70% gate; sparse stratum collapses to 0%).

## The three numbers (n = 542 confident glyco PSMs)

| Metric | Value |
|---|---|
| Oxonium-fire rate | **100.0%** (542/542) |
| **Searchable-backbone OVERALL** | **37.3%** (202/542) |
| **Searchable-backbone SPARSE stratum** (≤2 true core-Y rungs, n=148) | **0.0%** (0/148) |

Gate: PASS requires overall ≥ 0.70 AND sparse not collapsing (≥ ~0.55). **Both fail.**

## Run provenance

- Binary commit: see the commit that adds `crates/andes-glyco/src/bin/glyco_probe.rs`
  (this file is committed in the same commit). Branch `glyco-phase1`.
- Centroid fix (non-negotiable #1) applied to `backbone.rs` in the same branch:
  `solve_backbone` now reports the **vote-weighted centroid** of each merged
  near-mass cluster (not the lowest-mass edge); test
  `solve_backbone_reports_cluster_centroid_not_low_edge` added; all 7 unit tests pass.
- Gate window (non-negotiable #2): symmetric
  `in_window(solved, truth) = |solved - truth| <= max(truth*20e-6, 0.01)`.
- Solver call: `solve_backbone(&peaks, precursor_neutral, prec_z, 20.0, 5)`
  (tol 20 ppm, top_k = 5), oxonium gate `oxonium_gate(&peaks, 0.10, 20.0)`.

## Truth source + cutoff

- **Dataset PXD025455** (PRIDE 2021/05). The exact file named in the plan
  (`Pool_HCC_early_Fc3_r1`) does **not exist** in this project; the closest
  intact-N-glycopeptide DDA run with per-PSM ground truth is
  **`HCC_pool_Late_Fc3_r1`** (same study, HCC pooled serum, Fc3 fraction) — used here.
- **Ground truth = the reference engine-3.2 labile-search-mode `NGLYCAN` pepXML** (NOT a commercial glyco engine;
  this dataset ships the reference engine pepXML, not a commercial glyco engine). The glycan is encoded as a
  precursor **mass-offset delta** (`mass_offsets` open-search list), so the bare
  peptide backbone mass is exactly `calc_neutral_pep_mass` — the reference engine already
  excludes the glycan from the peptide mass. No subtraction needed.
- **backbone_mass = `calc_neutral_pep_mass`** (sequence + non-glycan fixed/var mods:
  Cam-C 57.0215, ox-M 15.9949, N-term acetyl 42.0106). The glycan delta =
  `precursor_neutral_mass − calc_neutral_pep_mass`.
- **Confidence cutoff** (documented): rank-1, **target** (non-`rev_`) PSMs with
  (a) glycan delta > 200 Da, (b) delta snaps to a known `mass_offsets` glycan
  composition within 0.05 Da (1 of 182 catalogued offsets — strong composition
  confirmation), and (c) the reference engine `expect` (e-value) ≤ 1.0. → **542 scans**.
  (Tighter `expect ≤ 0.1` yields 431; the rate is insensitive to the cutoff —
  the failure is structural, not a low-confidence-truth artifact.)

## The single most important finding (failure mode)

**The signal is recoverable; the solver's candidate ranking/pruning throws it away.**

Diagnostic (`glyco_diag.rs`): for the **true** backbone mass, the core-Y ladder
(≥2 of {Y0,Y1..Y5}) is actually present in the spectrum in **85.6%** of scans.
When `top_k` is raised from 5 → 500, the true backbone is recovered within
20 ppm in **80.3%** of scans. So `solve_backbone` *can* form the correct cluster
~80% of the time — but at the spec'd `top_k = 5` it ranks it below ~5 spurious
higher-vote low-mass clusters (artifact of voting every peak at every charge
1..=z as a possible Y-ion, which piles votes into low-mass bins). Result:
only **37.3%** of true backbones survive in the top-5.

The **sparse stratum collapses to 0%**: with ≤2 true core-Y rungs the solver
never lands the backbone in the top-5 — it is entirely dependent on a rich
(≥3-rung) ladder.

### Implication for the gate

As specified (oxonium 0.10/20 ppm, `solve_backbone(.., 20.0, 5)`, symmetric
±20 ppm), this is a **FAIL** — glyco mode is NOT worth building on this solver
as-is. BUT the failure is an addressable *ranking/precision* problem
(candidate de-noising, complement-pair confirmation, intensity-weighted votes,
charge-aware Y-ion deisotoping), not a fundamental backbone-unrecoverability —
the 80% ceiling at top_k=500 is the upside if the ranking is fixed. A Phase-1.5
"fix candidate ranking, re-run the same gate" is the honest next step before
either greenlighting or killing glyco mode.
