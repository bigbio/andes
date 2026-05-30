# Chimeric two-pass cascade (MS1-localized second-peptide search) — design

**Date:** 2026-05-30
**Branch:** `feat/chimeric-dda-plus`
**Type:** Feature design (the speed solution for chimeric; phased, bench-gated)
**Supersedes the chimeric ON path:** the blind wide-window search (and both
fragment-index approaches A/B).

## Motivation (evidence-grounded)

Chimeric NO_RESCORE delivers real, entrapment-validated PSM gains (Astral +116%,
PXD +21.6% vs Java) but is ~3× slower than Java. A `perf` profile
(`2026-05-30-chimeric-cost-profile.md`) shows the cost is **per-candidate SCORING
(~65%)** — `score_psm` + `edge_score` + `psm_edge_score` + node-mass lookups +
roundf-heavy peak matching — paid for THOUSANDS of wide-window candidates per scan.
NOT emission (`compute_psm_features` absent), NOT MS1 features, NOT GF (~6%). Both
fragment-index approaches failed: A degenerates; B's top-K prefilter drops the
low-evidence co-isolated tail that IS the gain.

**Insight:** the wide window exists only to find co-isolated peptides, but their
precursors are visible in MS1. So instead of scoring the whole window on every
scan, do a narrow first pass, then a **targeted** second pass at the precursor
masses MS1 says are actually co-isolated. This cuts the candidate-scoring
bottleneck while keeping the real co-isolated peptides (MS1 envelopes ARE the real
co-isolated precursors — no recall loss from fragment-prefiltering).

## Goal

Recover most of the chimeric +116% Astral PSM gain at wall time **< Java 6:18**,
with FDR validated by the entrapment ruler.

## Design

### §1 Two-pass architecture
- **Pass 1 — narrow search** (the existing non-chimeric path, unchanged): one
  primary peptide per scan. Already faster than Java; emits the primary PSMs.
- **Pass 2 — MS1-gated second-peptide search:** ONLY for scans where MS1 shows a
  co-isolated precursor in the isolation window. Score a handful of candidates near
  the MS1-detected co-isolated mass(es), on the residual spectrum; emit secondary
  PSMs.
- **Combine:** primary + secondary PSMs → one PIN → Percolator → FDR (entrapment-validated).

The whole cascade is behind `--chimeric` (the flag is reused; the ON path now means
"two-pass cascade", not blind wide-window). `--chimeric off` = Pass-1-only =
bit-identical to the current narrow search.

### §2 Pass-2 detail (the new component)
New module `crates/search/src/coisolation.rs`:
1. **Co-isolated precursor detection.** For each MS2 scan, take its linked MS1
   (`Ms1Link`, already loaded under `--chimeric`), restrict to the isolation window
   `[selected_mz − lower, selected_mz + upper]`, peak-pick + find isotope envelopes
   (charge-deconvolved), EXCLUDE the envelope at the selected precursor m/z → a list
   of co-isolated `(neutral_mass, charge)` candidates (typically 0–2). Reuse the
   `precursor_isotope_match` / averagine machinery from Phase 2. If the list is
   empty → scan skipped (the MS1 gate).
2. **Residual spectrum.** Subtract the Pass-1 primary peptide's matched charge-1 b/y
   peaks (reuse `matched_peak_keys`) so shared peaks don't inflate the secondary's
   score (MaxQuant second-peptide convention).
3. **Targeted scoring.** For each co-isolated `(mass, charge)`: enumerate candidates
   within the NARROW precursor tolerance of that mass (small set — the candidate-count
   cut), `score_psm` + `psm_edge_score` + per-charge GF SpecEValue on the residual,
   keep the best distinct peptide; emit as a secondary PSM (SpecId `<scan>_2`).

### §3 Speed model (why it works)
Pass 1 cost ≈ narrow (beats Java). Pass 2 cost ≈ (#co-isolated scans) × (few
candidates each) × scoring — a small fraction, because (a) only co-isolated scans
run it, and (b) each scores a NARROW candidate set at the MS1 mass, not the whole
window. The proven ~65% candidate-scoring bottleneck collapses: instead of
thousands of wide-window candidates/scan, Pass 2 scores a handful only on
co-isolated scans. Target: total < Java 6:18.

### §4 Correctness & gates
- **Off (`--chimeric off`):** Pass-2 code unreached → narrow search bit-identical.
- **Recall:** Pass-1+2 PSMs @1% recover most of the blind-chimeric Astral gain
  (gate: a high fraction of the 77,287; the exact bar set after the first
  measurement, since MS1-gating may not catch envelope-less co-isolation). The
  primary set alone should already ≈ narrow/Java; the secondary set is the gain.
- **FDR honesty:** re-run the entrapment harness on the cascade output — secondary
  PSMs must not inflate entrapment FDP beyond nominal.
- **Speed:** total wall < Java 6:18 on Astral; ≤ current narrow + small margin on PXD.

### §5 Phase decomposition
- **P1** `coisolation.rs`: MS1 isolation-window co-isolated-precursor detector
  (envelope find + exclude selected + charge/mass output), unit-tested on synthetic MS1.
- **P2** Pass-2 targeted residual search: given a scan, primary peptide, and a
  co-isolated `(mass,charge)`, return the best secondary PSM (reuses score_psm/GF),
  unit-tested.
- **P3** two-pass driver in the binary/match_engine: run Pass 1 (narrow), then Pass 2
  on MS1-gated scans, combine into one PIN; `--chimeric off` bit-identical.
- **P4** PXD + Astral gates: recall vs blind-chimeric, wall < Java, entrapment FDP
  preserved; bench note + PR.

## Risks & mitigations
| Risk | Mitigation |
|---|---|
| MS1 gate misses envelope-less co-isolation → recall < +116% | measure the recovered fraction; if too low, fall back to all-scans-MS1-localized (the option-2 variant) — a `coisolation.rs`-local change |
| Co-isolated detection false positives → spurious secondary PSMs | residual scoring + entrapment FDP gate; Percolator down-weights |
| Residual subtraction removes shared real peaks of the secondary | subtract only the primary's MATCHED peaks (charge-1 b/y), not all; measured via recall |
| MS1 peak-picking cost | only on `--chimeric`; isolation-window-restricted (tiny m/z range); profile if needed |
| Two-pass FDR comparability | one combined Percolator run; entrapment ruler is the authority |

## Out of scope
- Higher-order (3rd+) co-isolation (N-pass) — start with one secondary peptide.
- TMT PSM gap (Lever-2a) — orthogonal.
- The roundf 10% micro-opt (a separate flat speedup for narrow+chimeric).
- Fragment-index approaches (A/B) — refuted; superseded by this cascade.

## References
- `2026-05-30-chimeric-cost-profile.md` (candidate scoring = ~65%, the cost this targets).
- `2026-05-29-entrapment-fdp-reversal.md` (the validated gains to preserve + the ruler).
- `2026-05-30-sage-index-astral-and-chimeric-speed-conclusion.md` (why prefiltering failed).
- `2026-05-28-chimeric-phase2-bench.md` (MS1/averagine machinery to reuse).
- MaxQuant second-peptide search; MSFragger-DDA+ MS1 refinement.
