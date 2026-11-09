# andes benchmarks: supplementary tables

Tables behind the rows of [README.md](README.md) that most readers do not need: per-run and per-level
breakdowns, the full ProteoBench leaderboard, and diagnostics that only andes exposes. Builds,
settings and verdicts are those of the matching README section.

## timsTOF per run

From README section 2.1, timsTOF DDA-PASEF.

Per run (mean over six draws):

| run | tool | PSMs | true FDP | PSMs at 1% | peptides | wall | memory |
|---|---|---:|---:|---:|---:|---:|---:|
| PXD053748 `..._4363` (plasma, 30 min) | andes | 15,642 | 0.79% | 15,880 | 3,910 | 55 s | 2.4 GB |
| | Comet | 14,143 | 1.14% | 13,964 | 3,493 | 159 s | 6.2 GB |
| PXD053748 `..._4367` (plasma, 50 min) | andes | 18,916 | 0.90% | 19,133 | 3,865 | 56 s | 2.5 GB |
| | Comet | 16,139 | 1.15% | 15,865 | 3,396 | 176 s | 6.4 GB |
| PXD064202 `..._2442` (astrocytes) | andes | 36,092 | 0.88% | 36,279 | 14,682 | 54 s | 2.3 GB |
| | Comet | 30,848 | 1.31% | 29,957 | 13,248 | 186 s | 6.0 GB |
| PXD064202 `..._2403` (astrocytes) | andes | 58,578 | 0.91% | 58,828 | 24,146 | 62 s | 2.7 GB |
| | Comet | 50,814 | 1.31% | 49,517 | 22,007 | 222 s | 6.1 GB |

andes alone, one shuffle draw and the natural yeast design (PSMs / true FDP; "no conserved" removes
hits on yeast proteins with a human homolog at 40% identity over 60% length, 907 of 6,733):

| run | shuffle PSMs | shuffle true FDP | shuffle PSMs at 1% | natural PSMs / true FDP / no conserved / PSMs at 1% |
|---|---:|---:|---:|---|
| `..._4363` | 15,411–15,474 | 0.52–0.60% | 15,938–15,952 | 16,062–16,108 / 0.86–1.01% / 0.80–0.95% / 16,076–16,219 |
| `..._4367` | 18,900–18,957 | 0.82–0.94% | 19,031–19,152 | 19,456–19,507 / 0.91–0.98% / 0.79–0.86% / 19,528–19,576 |
| `..._2442` | 36,052–36,125 | 0.94–1.01% | 36,047–36,178 | 36,873–36,982 / 1.44–1.54% / 0.93–1.02% / 36,040–36,071 |
| `..._2403` | 58,814–58,873 | 1.12% | 58,542–58,637 | 59,739–59,761 / 1.17–1.24% / 0.72–0.78% / 59,068–59,260 |

## ProteoBench

From README section 2.2.

### andes on the two ProteoBench modules


Setup: "Quant LFQ DDA ion" modules, six runs each: QExactive (PXD028735, Q Exactive HF-X) and
Astral (PXD070049, 15 min, 50 ng); human/yeast/E. coli at known A:B ratios. andes `09e1c00d`,
2026-11-05, defaults (`--predict rt`), Percolator through `--rescore` (seed 42), scored with
`proteobench` 0.18.6 ("Custom" format via [`proteobench/andes_to_proteobench.py`](proteobench/)).
Quantity: `quantity_normalized` (run loading and gradient drift removed) unless a row says raw.
Cells read precursors ≥ 3 runs / median |ε| / E. coli |ε| / CV.

| view | QExactive | Astral |
|---|---|---|
| default: transfer q ≤ 0.01, normalised | 63,796 / 0.161 / 0.263 / 0.117 | 36,675 / 0.186 / 0.252 / 0.112 |
| raw `quantity` (what quantms receives) | 63,796 / 0.222 / 0.309 / 0.245 | 36,675 / 0.204 / 0.254 / 0.113 |
| transfer q ≤ 0.05 | 71,299 / 0.175 / 0.441 / 0.132 | 43,286 / 0.205 / 0.349 / 0.127 |
| identified rows only (no match between runs) | 52,166 / 0.154 / 0.203 / 0.105 | 31,739 / 0.186 / 0.235 / 0.100 |
| precursors in ≥ 1 run / all 6 runs (default) | 79,492 / 41,370 | 49,012 / 24,393 |
| by abundance, quartile Q1 (lowest): \|ε\| / E. coli \|ε\| / CV | 0.279 / 0.461 / 0.195 | 0.319 / 0.406 / 0.194 |
| Q2 | 0.169 / 0.263 / 0.130 | 0.208 / 0.303 / 0.138 |
| Q3 | 0.136 / 0.205 / 0.105 | 0.161 / 0.212 / 0.103 |
| Q4 (highest) | 0.122 / 0.183 / 0.077 | 0.122 / 0.197 / 0.060 |
| opt-in `--predict rt,ms2` / `ms2,peptdeep-rt,deeplc-rt`: precursors (default on that build) | 67,906 / 68,077 (65,391) | 34,413 / 33,987 (33,876) |
| CPU-s, six runs: search + Percolator / library / extract / finalize per threshold | 5.8k / 4 / 0.2k / 0.2k | 5.8k / 4 / 30 / 0.1k |
| wall: search per run / extract per run / finalize | 3–4 min / 6 s / 34 s | 3–4 min / 2 s / 19 s |
| PSM FDP at `q ≤ 0.01` (entrapment, T/E 1.075) | 0.84–0.85% | — |
| entrapment-only share of precursors ≥ 3 runs, transfer q 0.01 / 0.05 | 0.20–0.21% / 0.26–0.27% | — |
| false-donor share of transfers, entrapment / identity-decoy estimate, q 0.01 | 0.22–0.26% / 0.28–0.29% | — |
| same, transfer q 0.05 | 0.37–0.40% / 0.42–0.47% | — |

- The error-control rows were measured on QExactive only, seeds 42, 1 and 2, at `09e1c00d`. The
  decoy estimate was above the measured share in all three seeds at both thresholds; three seeds do
  not make it an upper bound. The wrong-peak part has no entrapment check.
- The `--predict` rows are from 2026-10-21 at `a27a72c`, before the extraction change of #240.
- Raw `quantity` is less accurate than most leaderboard submissions; normalisation differs between
  tools, so the |ε| column does not compare extraction alone.

### Against the ProteoBench leaderboard

Setup: public results of both modules (last updated 2026-07-16, read 2026-10-21, scored by
`proteobench` 0.17.0; andes by 0.18.6). Per tool the submission with the most precursors and, when
different, the one with the lowest median |ε|, at ≥ 3 runs. "Transfer control" is what limits the
error of match-between-runs values, from the submission's parameters, comment and output; "not
checked" means the parameter file could not be retrieved.

QExactive:

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | transfer control |
|---|---|---:|---:|---:|---:|---|
| andes, default | yes | 63,796 | 0.161 | 0.263 | 0.117 | transfer q ≤ 0.01 |
| andes, transfer q ≤ 0.05 | yes | 71,299 | 0.175 | 0.441 | 0.132 | transfer q ≤ 0.05 |
| andes, raw quantity | yes | 63,796 | 0.222 | 0.309 | 0.245 | as default |
| andes, identified rows only | no | 52,166 | 0.154 | 0.203 | 0.105 | none |
| FragPipe 23.0, most precursors | yes | 63,701 | 0.200 | 0.306 | 0.156 | IonQuant MBR ion FDR 0.01 |
| FragPipe 23.0, lowest \|ε\| | yes | 54,374 | 0.190 | 0.257 | 0.154 | as above |
| MaxQuant 1.5.3.30, most precursors | yes | 44,299 | 0.210 | 0.386 | 0.219 | not checked |
| MaxQuant 2.4.13.0, lowest \|ε\| | no | 29,737 | 0.192 | 0.263 | 0.199 | none |
| i2MassChroQ 1.0.16, most precursors | yes | 74,697 | 0.222 | 0.460 | 0.132 | no transfer setting |
| i2MassChroQ 1.2.0, lowest \|ε\| | no | 49,433 | 0.198 | 0.263 | 0.120 | none |
| AlphaPept 0.5.0, most precursors | yes | 37,760 | 0.219 | 0.281 | 0.165 | not checked |
| AlphaPept 0.5.3, lowest \|ε\| | no | 36,176 | 0.214 | 0.255 | 0.165 | none |
| Sage 0.15.0-beta.2 | yes | 81,933 | 0.388 | 0.877 | 0.393 | peptide q on IDs only; peak q not applied |
| MSAngel 2.2.10 | yes | 60,206 | 0.216 | 0.514 | 0.227 | not checked |
| ProlineStudio 2.3.0-SNAPSHOT | yes | 55,639 | 0.210 | 0.508 | 0.226 | not checked |
| WOMBAT 0.9.11, most precursors | yes | 55,943 | 0.211 | 0.692 | 0.235 | not checked |
| WOMBAT 0.9.11, lowest \|ε\| | yes | 41,016 | 0.195 | 0.394 | 0.213 | not checked |
| PEAKS 13 | no | 72,741 | 0.289 | 0.741 | 0.216 | none |
| quantms v1.3.1dev (OpenMS ProteomicsLFQ) | no | 33,040 | 0.209 | 0.387 | 0.224 | none |
| quantms, our run 2026-10-17, default | no | 34,677 | 0.210 | 0.401 | 0.225 | none |
| quantms, our run, feature-linking MBR | yes | 39,257 | 0.221 | 0.563 | 0.234 | feature quality score |

Astral:

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | transfer control |
|---|---|---:|---:|---:|---:|---|
| andes, default | yes | 36,675 | 0.186 | 0.252 | 0.112 | transfer q ≤ 0.01 |
| andes, transfer q ≤ 0.05 | yes | 43,286 | 0.205 | 0.349 | 0.127 | transfer q ≤ 0.05 |
| andes, raw quantity | yes | 36,675 | 0.204 | 0.254 | 0.113 | as default |
| andes, identified rows only | no | 31,739 | 0.186 | 0.235 | 0.100 | none |
| FragPipe 23.0, most precursors | yes | 39,057 | 0.212 | 0.250 | 0.127 | not checked |
| FragPipe 22.1-build02, lowest \|ε\| | yes | 29,105 | 0.181 | 0.405 | 0.078 | not checked |
| MaxQuant 2.6.7.0 | yes | 20,283 | 0.160 | 0.158 | 0.092 | within a condition only; not checked |
| i2MassChroQ 1.2.6, most precursors | yes | 36,469 | 0.197 | 0.417 | 0.095 | not checked |
| i2MassChroQ 1.2.9, lowest \|ε\| | no | 24,803 | 0.178 | 0.245 | 0.077 | none |
| Sage 0.14.6 | yes | 30,319 | 0.167 | 0.338 | 0.070 | peptide q on IDs only; peak q not applied |
| DIA-NN 2.6.0, most precursors | yes | 93,237 | 0.277 | 0.691 | 0.159 | run and global precursor q ≤ 0.01 |
| DIA-NN 2.6.0, lowest \|ε\| | yes | 16,363 | 0.123 | 0.168 | 0.079 | as above plus MS1-level q ≤ 0.01 |

On the precursors andes shares with a submission (default, transfer q 0.01; median |ε| / E. coli |ε|
/ CV; QExactive only, Astral per-submission tables could not be retrieved):

| submission | shared | submission | andes raw | andes normalised | only in the submission |
|---|---:|---|---|---|---|
| FragPipe 23.0, most precursors | 51,722 | 0.189 / 0.274 / 0.155 | 0.207 / 0.276 / 0.241 | 0.147 / 0.225 / 0.108 | 11,979 at 0.272 / 0.708 |
| FragPipe 23.0, lowest \|ε\| | 48,938 | 0.187 / 0.243 / 0.155 | 0.203 / 0.268 / 0.241 | 0.144 / 0.207 / 0.108 | 5,436 at 0.220 / 0.510 |
| i2MassChroQ 1.0.16, most precursors | 50,910 | 0.207 / 0.348 / 0.125 | 0.209 / 0.291 / 0.243 | 0.152 / 0.241 / 0.112 | 23,787 at 0.269 / 0.948 |
| i2MassChroQ 1.2.4, 57,787 precursors | 49,853 | 0.202 / 0.348 / 0.123 | 0.209 / 0.291 / 0.242 | 0.151 / 0.242 / 0.112 | 7,934 at 0.279 / 0.835 |
| Sage 0.15.0-beta.2, our reproduction | 54,633 | 0.309 / 0.613 / 0.319 | 0.215 / 0.291 / 0.245 | 0.156 / 0.246 / 0.114 | 25,724 at 0.683 / 1.364 |

- QExactive: at least six submissions report more precursors than andes (Sage and the largest
  i2MassChroQ and PEAKS submissions among them), at a median |ε| of 0.214–0.388. FragPipe 23.0
  reaches 0.190 at 54,374 precursors; i2MassChroQ 1.2.0 has a CV of 0.120 without match between runs.
- Astral: MaxQuant 2.6.7.0 and Sage 0.14.6 report fewer precursors than andes at lower median
  |ε| (0.160, 0.167); i2MassChroQ 1.2.6 about as many; FragPipe 23.0 more (39,057, 0.212); DIA-NN
  2.6.0 reports 93,237 at 0.277, and its 16,363-precursor submission has the lowest median |ε|
  (0.123). andes's CV (0.112) is above most submissions. Normalisation changes little here.
- Sage's per-submission table was not retrievable; the reproduction (2026-10-21, same version, same
  six mzML files, submission parameters, per-charge output) gives 80,357 precursors, |ε| 0.388,
  E. coli 0.878, CV 0.394 (submission: 81,933, 0.388, 0.877, 0.393). Most Sage precursors that
  andes lacks are charge states Sage never identified or that fail its own peak q-value of 0.01,
  which it reports but does not apply.

## UPS1 per level

From README section 2.3, UPS1 detection at the lowest three levels.

Per level, MBR on (rows; within 2x; share within 2x; median |log2 error|):

| set | level (amol) | andes | MaxQuant | FragPipe |
|---|---:|---|---|---|
| PXD009740 | 170 | 295; 66; 33.6%; 1.29 | 323; 38; 14.2%; 1.97 | 189; 74; 52.9%; 0.96 |
| | 382 | 467; 153; 59.3%; 0.85 | 532; 127; 35.2%; 1.26 | 345; 133; 67.0%; 0.63 |
| | 859 | 838; 283; 70.3%; 0.72 | 787; 252; 51.6%; 0.99 | 633; 251; 78.8%; 0.58 |
| PXD009815 | 10 | 174; 0; 0.0%; 7.23 | 195; 0; 0.0%; 7.16 | 131; 1; 0.8%; 6.53 |
| | 50 | 194; 3; 2.1%; 4.71 | 189; 4; 2.1%; 4.39 | 140; 19; 17.1%; 3.97 |
| | 100 | 197; 21; 19.8%; 3.03 | 219; 19; 15.5%; 2.94 | 172; 30; 25.6%; 2.26 |

## Masked wrong-peak rate (andes only)

From README section 2.3.

At `quant_q ≤ 0.01`, pooled over runs; 5% of identified precursor-run pairs are masked from the
library, re-quantified as transfers, and counted when the apex lies further from the identified apex
than the fit's peak width. No other tool exposes this quantity.

| set | seeds 42 / 1 / 7 | MS1-only transfers, seed 42 |
|---|---|---:|
| PXD009740 | 0.50 / 0.50 / 0.49% | 0.82% |
| PXD009815 | 0.55 / 0.55 / 0.60% | 0.92% |
| PXD001819 | 0.46 / 0.41 / 0.43% | 1.01% |
| QE | 0.77 / 0.80 / 0.78% | 1.30% |
| Astral | 0.66 / 0.66 / 0.62% | 0.66% |
| PXD014415 | 0.49 / 0.50 / 0.48% | 1.09% |

- Transfers with MS2 evidence read 0.18–0.65% (seed 42). The blank-run set has 126–130 masked pairs
  and none was wrong. Among runs with at least 200 masked pairs, 0 to 5 per set (of 6 to 40) are
  above 1%.
