# Benchmark methods

## Verdicts

Every compared measure carries a verdict for andes against each other tool, set by one rule. The
ratio is the worse value over the better one, in the direction of the metric (more PSMs or
precursors is better; less error, CV, false calls, time or memory is better).

| verdict | ratio between andes and the other tool |
|---|---|
| 🟢 well ahead / 🔴 well behind | 1.5x or more |
| 🟢 ahead / 🔴 behind | 1.1x up to 1.5x |
| 🟡 on par | under 1.1x, or the other value lies inside andes's seed or draw range, or both are counts below 10 |
| ⚪ not compared / not tested | no other tool measured, no common denominator, or no decision was made |

- **Ranges and noise.** A range over seeds, draws, runs or sets is judged at its end nearest the
  other tool, so a verdict holds for the whole range. A single run (marked ¹) has no measured
  spread; the 1.1x floor is wider than the andes seed spread of PSM counts (under 1.2%), so a gap
  under 10% is never called a lead or a loss.
- **Shares of a reference** (glyco spectra confirmed): the ratio is 100% over the confirmed share.
- **Rates with 95% intervals** (TMT match between runs): 🟡 on par when the intervals overlap,
  otherwise ahead or behind; no "well".
- **Error control** (entrapment true FDP at nominal 1%): 🟢 on target when every reading is at
  or below 1.1%, 🟡 near nominal when the range crosses 1.1%, 🔴 above nominal when every reading
  is above it. Natural entrapments with hits on conserved homologs are context (⚪), not a reading.
- Several sets or tools in one row get one verdict each. Caveats (depth at which an error is
  measured, opt-in modes, unmatched settings) go in the note and do not change the verdict.

## Measurement conditions

Identification rows were measured on 2026-10-23 at
`2051616b`, except UPS1 and timsTOF (2026-11-07 and 2026-11-09 at `2ee76dd6`). Label-free
quantification rows were measured on 2026-11-05 at `09e1c00d`, TMT rows at `0cfa03cc`. A row says
so when it differs. PSM counts, wall time and memory of Astral, TMT and phospho come from
the benchmark VM (8 threads) on the plain databases the Comet rows used. True FDP and all other
rows come from a cluster node (8 threads) with an entrapment component added to the database
(rules under [Method notes](#method-notes)). Percolator 3.7.1 `-Y --num-threads 8`, seeds 42, 1 and 2 unless a setup
line says otherwise. "At a true FDP of 1%" means the deepest score-sorted list whose
entrapment-estimated FDP stays at or below 1%.

## Method notes

- **True FDP by entrapment.** FDP = (entrapment hits / accepted) × (1 + T/E), with T/E the ratio of
  searchable target to entrapment peptides, measured for each database. A PSM is an entrapment hit
  only when all its proteins are entrapment and its peptide, with I folded to L, occurs in no target
  protein. `entrap_common.py` holds this rule for `run.sh`, `build_databases.sh` and the glyco
  evaluators. A natural entrapment (a foreign proteome) can be contaminated by conserved homologs;
  the shuffle is the primary reading.
- **Entrapment designs.** Peptide-level shuffles ([`build_shuffled_entrap.py`](build_shuffled_entrap.py);
  each tryptic peptide shuffled with its C-terminal residue fixed, after Wen et al., Nat. Methods
  2025) give every entrapment peptide a target of the same composition and mass. Three realisations
  (seeds 1, 2, 3) for Astral, phospho and glyco; one for UPS1, TMT, TKO and timsTOF. The earlier
  whole-protein shuffle left 89.6% of entrapment peptides on the HYE database without a
  same-composition target, and read 0.81% on Astral where peptide-level databases read 0.92–1.02% and
  *Arabidopsis* (T/E 2.228) 0.99% (#225); on phospho it agreed (1.05% against 0.95–1.04%). FDPs marked
  "whole-protein shuffle" may read low.
- **Conserved-protein filter.** An `ENTRAP_` protein is conserved when its best local alignment
  (BLOSUM62, gap open 11, extend 1, against target proteins preselected by shared 4-mers) covers at
  least 60% of it at 40% identity or more. For UPS1, 269 of 4,531 *E. coli* proteins qualify (719 at
  30% / 50%; 23 at 60% / 60%); the cut was fixed in advance. Hit counts were not re-run at the other
  cuts; T/E without those proteins is 2.872 at 30% / 50% and 2.456 at 60% / 60%.
- **Percolator.** Version 3.7.1, `--seed 42 -Y`, on plain FASTA with andes `XXX_` decoys, for every
  engine on the same host. The q-value has a floor of 1 / (number of top targets), so small or sparse
  runs return few identifications at 1%; pool glyco fractions for this reason.
- **ProteoBench terms.** "Min runs" is the number of the six runs in which a precursor
  (peptidoform and charge) must be quantified. Epsilon is the measured log2 A/B ratio minus the
  expected one; the headline is its median absolute value. CV is the median within-condition
  coefficient of variation, averaged over the two conditions. A transfer q-value scores the peak a
  precursor is matched to in a run where it was not identified; it is not an identification FDR.
- **Truth-design evaluation.** Every tool's table is reduced to raw precursor-by-run intensities, then:
  per-run median normalisation over precursors quantified in all runs, Tukey median polish per
  protein, proteins with ≥ 2 precursors tested with a moderated t, Benjamini–Hochberg 5% per
  condition pair. A true positive is a significant protein of a changing species with the expected
  sign; a false call is a significant protein of a constant species (yeast in UPS1 sets, human in the
  mixtures) or one of a changing species with the wrong sign. DE is scored on the common-tested family
  (rows that andes, MaxQuant and FragPipe all test; PXD030780: andes and FragPipe). Earlier versions
  tested each arm on its own proteins (4,699 against MaxQuant's 4,312 on QE) and counted wrong-sign
  calls as true positives; both are corrected. normalisation and roll-up are the same code for all
  tools. Accuracy is scored on cells every listed tool identifies, against the replicate-based
  expectation. UPS1 detection counts a precursor at a level when one value is within 2x of its own
  high-level reference.
- **Absence classes.** For a precursor, d is the median log2 intensity over runs that contain yeast
  minus that over human-only runs, normalised on human proteins. It is *enriched* when d exceeds a
  boundary learned from human control precursors by the number of human-only runs the precursor is
  seen in (6.67, 5.62, 3.94, 1.07 log2 for 1, 2–4, 5–9, 10 or more runs), *not specific* otherwise,
  *only absent* when seen in human-only runs only, and *no reference* when andes's table of every row
  (2,892,797 rows) lacks it. Enriched-only is a lower bound of species-specific errors.
- **Masked wrong-peak rate** (results in [ADDITIONAL_BENCHMARKS.md](ADDITIONAL_BENCHMARKS.md#masked-wrong-peak-rate-andes-only)). 5% of identified precursor-run pairs are removed from the library,
  quantified again as transfers, and counted when the apex lies further from the original apex than
  that fit's peak width, at `quant_q ≤ 0.01`.
- **TMT knockout layout.** MET6 (P05694) is knocked out in 126–128N, HIS4 (P00815) in 128C–130N and
  URA2 (P07259) in 130C–132N; 132C–134N is the parental strain. For each PSM of those proteins, divide
  the mean of its knockout channels by the mean of the other channels in `<stem>.tmt.parquet` (rows
  with `quant_gate_passed`), report the median over scans whose reporter window starts below 126 Th,
  and repeat on PSMs with `purity ≥ 0.9`. The MBR evaluation reads `run.tmt.parquet` (`evidence_type`,
  `quant_q`, `quant_gate_passed`, channel columns).
- **Comparisons.** Identification comparisons are runs on the same host with matched parameters and the
  same rescoring; ProteoBench comparisons are public leaderboard submissions. Wall times compare only
  within one host and session; counts are portable.
