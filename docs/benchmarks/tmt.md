# TMT quantification benchmarks

## Knockout interference and extraction (PXD017803)

Setup: yeast triple knockout, TMTpro 16-plex; two files each of hrMS2, SPS-MS3 and RTS-MS3, as
mzML, with the kit impurity sheet. andes defaults including the quantification gate (precursor purity
≥ 0.7, Percolator PEP ≤ 0.05), `--rescore`, 16 threads, `0cfa03cc`, 2026-11-05. MaxQuant 2.7.0.0
run 2026-10-18 on vendor files (reporter intensity corrected, 1% PSM and protein FDR). FragPipe 21.1;
quantms (OpenMS IsobaricWorkflow, `96f263c`). Interference is the median ratio of knockout to
non-knockout channel intensity over PSMs of the three knocked-out proteins (lower is better; hrMS2
only scans whose window starts below 126 Th). Reporter extraction and impurity correction match OpenMS
IsobaricAnalyzer exactly on the same scans.

| acquisition | andes PSMs at 1% (quantified) | andes interference | no gate / purity ≥ 0.9 | MaxQuant PSMs; interference | FragPipe | quantms | PSMs vs MaxQuant | interference vs MaxQuant / FragPipe / quantms |
|---|---:|---:|---|---|---:|---:|---|---|
| SPS-MS3 | 9,806 (8,635) | 0.062 | 0.117 / 0.023 | 6,883; 0.054 | 0.139 | 0.072 | 🟢 ahead | 🔴 behind / 🟢 well ahead / 🟢 ahead |
| RTS-MS3 | 14,270 (10,518) | 0.027 | 0.033 / 0.015 | 8,930; 0.038 | 0.055 | 0.037 | 🟢 well ahead | 🟢 ahead / 🟢 well ahead / 🟢 ahead |
| hrMS2 | 13,474 (3,481) | 0.192 | 0.202 / 0.135 | 9,677; 0.156 | 0.214 | 0.180 | 🟢 ahead | 🔴 behind / 🟢 ahead / 🟡 on par; with `--tmt-unmix on` 🟢 well ahead of all three |

- `--tmt-unmix on` on hrMS2 gives a median knockout ratio of 0.000 at a background CV of 0.123 (0.079
  without it).
- Few knockout PSMs stand behind each value (34–91 andes, 26–65 MaxQuant); read them as a direction.
  The tools do not quantify the same PSMs, and MaxQuant reports fewer PSMs at its own 1% FDP on all
  three sets.
- On hrMS2 only 4,502 of 13,474 andes PSMs have a reporter scan covering the reporter ions (3,481
  pass the gate; 388 fail purity, 633 fail PEP). On those reporter-bearing PSMs the evaluated counts
  are andes 3,278, FragPipe 3,038 and MaxQuant 2,241: 🟢 ahead of MaxQuant, 🟡 on par with FragPipe.

## Known ratios (PXD007683)

Setup: human + yeast TMT11 with known yeast ratios, 11 SPS-MS3 fractions, ion-trap CID MS2; andes
default gate, `--rescore`, 16 threads, `0cfa03cc`, 2026-11-05; MaxQuant 2.7.0.0 as above. Yeast log2
ratios against the human-normalised expectation (A/B 1.078; A/C 1.688, nominal 3:1 = 1.585; B/C 0.610).

| metric | andes | MaxQuant 2.7.0.0 | verdict | note |
|---|---:|---:|---|---|
| rank-1 target PSMs at 1% (quantified) | 149,121 (128,265) | 97,107 | 🟢 well ahead | +54% |
| yeast log2 ratio A/B / A/C / B/C | 1.191 / 1.589 / 0.397 | 1.263 / 1.700 / 0.424 | 🟢 well ahead / 🔴 well behind / 🔴 behind | distance to the expectation: andes 0.113 / 0.099 / 0.213, MaxQuant 0.185 / 0.012 / 0.186 |
| yeast ratios without the gate | 1.109 / 1.473 / 0.353 | — | ⚪ not compared | |
| protein-level MAE A/B / A/C / B/C | 0.230 / 0.213 / 0.213 | 0.222 / 0.159 / 0.186 | 🟡 on par / 🔴 behind / 🔴 behind | 1.3x on A/C |
| same, without the gate | 0.238 / 0.266 / 0.258 | — | ⚪ not compared | |
| proteins evaluated | 10,797 | 8,427 | 🟢 ahead | +28% |
| cost, one process, 16 threads | 9.2 min; 17.8 GB | — | ⚪ not compared | |

## Match between runs for reporter ions (PXD014750, hrMS2)

Setup: opt-in `--tmt --mbr` (hrMS2 only; refused for SPS-MS3 and ion-trap MS2). Six files (three
TMT11 plexes of human and yeast peptides at changing mixing, each acquired twice), human + yeast +
entrapment database, `q ≤ 0.01`, `--rescore`, 16 threads; andes `0cfa03cc`. Truth: human-only
peptides must be flat; yeast-only peptides must follow the plex's yeast profile. A value is wrong
when its factor-normalised, impurity-corrected profile deviates by more than 0.5 log2 (RMS); the
excess is the wrong rate of transferred values minus that of the tool's identified values at matched
species and signal. Only cells with a truth are scored, so denominators differ. MaxQuant and
SIMSI-Transfer rows are an earlier analysis of the same files with the same truth and threshold.

| tool | MBR | identified rows | transferred | wrong, transferred | wrong, identified (matched) | excess, points (95% CI) |
|---|---|---:|---:|---:|---:|---|
| andes `--tmt` | off | 199,730 | 0 | — | 1.96% (all scored rows) | — |
| andes `--tmt --mbr`, all candidates | on | 199,730 | 7,841 (6,260 scored) | 3.90% | 3.50% | +0.40 (−0.04 to +0.89) |
| andes `--tmt --mbr`, `quant_q ≤ 0.01` | on | 199,730 | 3,636 (3,267 scored) | 2.79% | 2.84% | −0.05 (−0.56 to +0.51) |
| MaxQuant 2.7.0.0 | on | 137,539 evidence rows | 7,499 cells | 3.77% | 2.24% | +1.53 (+1.11 to +1.98) |
| MaxQuant 2.7.0.0 | off | 137,539 evidence rows | 0 | — | — | — |
| SIMSI-Transfer 0.8.0 (MaRaCluster 1.04.1), p10 | on | MaxQuant scans | 16,511 cells | 2.25% | 1.85% | +0.40 (+0.18 to +0.64) |

Verdicts, by the 95% intervals of the excess (rule in [Methods](methods.md#verdicts)) and by cells
transferred:

| andes | vs MaxQuant 2.7.0.0 | vs SIMSI-Transfer | note |
|---|---|---|---|
| `--mbr`, all candidates: excess | 🟢 ahead | 🟡 on par | |
| `--mbr`, `quant_q ≤ 0.01`: excess | 🟢 ahead | 🟡 on par | |
| `--mbr`, `quant_q ≤ 0.01`: cells transferred | 🔴 well behind | 🔴 well behind | 3,636 against 7,499 and 16,511; opt-in mode |

- The entrapment share of accepted andes transfers is 1.6% (1.7% of all candidates).
- `--mbr` adds 3.9% rows (1.8% at `quant_q ≤ 0.01`). PSMs through the quantification gate: 131,161
  without it, 134,797 with it. Cost with the flag off 3.4k CPU-s, 8.0 min, 14.6 GB; with it 3.7k
  CPU-s, 8.4 min, 14.6 GB (one run each, one node type).
- The MaxQuant row counts all evidence rows of the types that carry an MS/MS scan, before reverse and
  contaminant hits are removed (isobaric MBR as in MaxQuant ≥ 1.6.12.0: `MULTI-MATCH-MSMS` rows); andes
  counts PSMs at `q ≤ 0.01`.
