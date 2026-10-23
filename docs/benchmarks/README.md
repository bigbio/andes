# andes benchmarks

andes is benchmarked on public PRIDE datasets for identification (high-res, TMT-labelled and
low-res LFQ, phospho, glyco) and for quantification (label-free with match between runs, and
TMT). Identifications are counted after Percolator at `q ≤ 0.01`. Where the database carries an
entrapment component, the page also gives the **true** false discovery proportion (FDP), which
is what tells you whether the nominal 1% holds. Other tools are named with the version compared.
Their numbers come either from runs on the same host with the same rescoring, or from a public
leaderboard. Section 2 gives the date and the andes commit behind each row.

## 1. Summary

All andes rows: `2051616b` (origin/main, 2026-10-23), default settings plus the dataset's mods,
Percolator 3.7.1 `-Y`, seeds 42, 1 and 2. Comparison rows keep the date they were measured.

**Identification** (PSMs at Percolator `q ≤ 0.01`, 8 threads):

| benchmark | dataset | andes | comparison | true FDP at nominal 1% | andes wall |
|---|---|---:|---|---|---:|
| Standard, high-res | Astral HCD LFQ, PXD070049 (1 file) | 52,169–52,437 PSMs, 31,147–31,238 peptides | Comet 2025.01 31,435 PSMs, 217 s (2026-10-05) | 0.89–1.08% (peptide-level shuffle); 0.99–1.05% (*Arabidopsis*) | 114–117 s |
| TMT-labelled, low-res | a05058, PXD007683 (1 file, ion-trap CID) | 13,443–13,516 PSMs | Comet 2025.01 10,504, 77 s (2026-10-05) | 0.67–0.76% (peptide-level); 0.83–1.05% (*Arabidopsis*) | 40–47 s |
| TMT-labelled, high-res | TKO hrMS2 TMTpro, PXD017803 (1 file) | 6,628–6,686 PSMs (entrapment database) | — | 0.98–1.31% (peptide-level) | 4 s † |
| Low-res LFQ | UPS1, PXD001819 (1 file, ion-trap CID) | 16,049–16,242 PSMs | Comet 2025.01 14,734, 42 s (2026-10-05); Java MS-GF+ v20240326 15,904 (earlier session) | 2.07–2.19% (natural *E. coli*); 0.91–0.98% (peptide-level) | 30–48 s |
| Phospho-enriched | mouse liver, PXD007653 (1 file, Q Exactive) | 38,596–38,629 PSMs, 27,693–27,724 phospho-bearing | Comet 2025.01 33,984 at 1.77% FDP, 1,096 s (2026-10-06; whole-protein shuffle) | 0.80–0.99% (peptide-level); 1.22–1.36% (*Arabidopsis*) | 337 s |
| timsTOF DDA-PASEF | plasma, PXD053748 (2 files) | 14,336–14,393 and 17,491–17,550 PSMs | — | 0.62–0.90% (peptide-level); 0.78–1.14% (natural yeast) | 43–51 s † |
| timsTOF DDA-PASEF | astrocytes, PXD064202 (2 files) | 34,210–34,249 and 55,590–55,609 PSMs | — | 0.98–1.03% (peptide-level); 1.29–1.55% (natural yeast) | 43–64 s † |
| Glyco, quick tier | mouse liver, PXD005553 (1 fraction) | 6,992–7,072 glycoPSMs | 88.2–88.4% of pGlyco2's and 89.3–89.7% of MSFragger-Glyco's identifications confirmed | 0.61–1.35% (peptide-level) | 5.8–6.6 min † |
| Glyco, deep tier | mouse liver, PXD005553 (5 fractions) | 31,666 ± 9 glycoPSMs (2026-09-05) | 78.9% (pGlyco2) and 88.0% (MSFragger-Glyco) confirmed | 1.11% ± 0.03 (whole-protein shuffle) | 23–29 min per fraction, 16 threads |

† On a cluster node (8 threads); the other walls are on the benchmark VM (8 cores), the host the
Comet times come from. Wall times compare only within one host.

**Quantification** (andes defaults, 2026-10-23 at `2051616b`):

| benchmark | dataset | andes | comparison (leaderboard snapshot 2026-07-16) | error control | andes cost |
|---|---|---|---|---|---|
| Label-free, match between runs (ProteoBench "Quant LFQ DDA ion QExactive") | PXD028735, 6 runs, Q Exactive HF-X, human/yeast/E. coli mix | 65,012 precursors in ≥ 3 runs; median \|ε\| 0.163, CV 0.118 (transfer q 0.01) | median \|ε\| 0.190 (FragPipe 23.0); CV 0.120 (i2MassChroQ 1.2.0, no MBR); most precursors 81,933 (Sage 0.15.0-beta.2, median \|ε\| 0.388, transferred values not error-controlled) | PSM FDP 0.83–0.85%; entrapment-only quantified precursors 0.20–0.22% | 6.9k CPU-s for 6 runs |
| Label-free, match between runs (ProteoBench "Quant LFQ DDA ion Astral") | PXD070049, 6 runs, Astral, human/yeast/E. coli mix | 37,352 precursors in ≥ 3 runs; median \|ε\| 0.193, CV 0.113 (transfer q 0.01) | median \|ε\| 0.123 (DIA-NN 2.6.0, 16,363 precursors); CV 0.070 (Sage 0.14.6); most precursors 93,237 (DIA-NN 2.6.0, median \|ε\| 0.277) | not measured | 6.4k CPU-s for 6 runs |
| TMT reporter ions | yeast triple knockout, TMTpro 16-plex, PXD017803 | knockout-channel interference 0.064 (SPS-MS3), 0.023 (RTS-MS3), 0.182 (hrMS2) | 0.139 with OpenMS precursor purity ≥ 0.9 (SPS-MS3, 2026-10-10) | — | — |
| TMT reporter ions | human + yeast TMT11, PXD007683 (11 fractions, SPS-MS3) | yeast A/C log2 ratio 1.591 (expected 1.585), protein MAE 0.212 | — | — | 9.6 min, 17.5 GB, 16 threads |
| Large study, match between runs | PXD042233, 500 runs (Q Exactive HF/HF-X, Exploris 480) | 243,895 precursors in ≥ 3 runs, 9,208 protein groups (transfer q 0.01) | — | PSM, precursor and protein q ≤ 0.01 | 146 CPU-h end to end (§2) |

**Opt-in modes** (Astral, PSMs / distinct peptides, Percolator seeds 42, 1 and 2; measured
2026-10-05 to 2026-10-13 on earlier builds and not re-measured at `2051616b`):

| mode | andes | without the mode |
|---|---|---|
| `--chimeric` (second and later peptides per spectrum), at an entrapment true FDP of 1% | 66,377–66,424 / 29,384–29,400 | 43,076–43,145 / 26,723–26,746 |
| `--refine` (PTM discovery, high-res only) | 48,312 PSMs, 28,728 peptides at `q ≤ 0.01`; its second-pass PSMs at ~1.2–1.8% FDP | 47,080 / 28,932 |

## 2. The benchmarks

**Where and how.** Every andes number on this page was measured on 2026-10-23 at `2051616b`
unless a row says otherwise. PSM counts, wall time, CPU time and peak memory of the four standard
sets (Astral, TMT a05058, UPS1, phospho) come from the benchmark VM (8 cores, 8 threads,
`/usr/bin/time`, two runs except phospho), on the plain databases the Comet rows used. The true
FDP of every set and all other rows come from a cluster node (8 threads), with each target
database extended by an entrapment component: three peptide-level shuffled realisations
(`build_shuffled_entrap.py` seeds 1, 2 and 3) for Astral, phospho and glyco, one realisation for
UPS1, TMT, TKO and timsTOF, plus a foreign (reviewed *Arabidopsis*) or natural entrapment where
one was used before. Percolator 3.7.1 `-Y --num-threads 8`, seeds 42, 1 and 2. "At a true FDP
of 1%" is the deepest score-sorted list whose entrapment-estimated FDP stays at or below 1%.

**Standard, high-res (Astral).** One 15-minute DDA run of the ProteoBench human/yeast/E. coli
sample against the 31,889-entry HYE FASTA. 52,169–52,437 PSMs and 31,147–31,238 peptides at
`q ≤ 0.01`. With the peptide-level entrapment (T/E 1.075) the nominal 1% reads 0.89–1.08% true
FDP over nine database × seed combinations, and andes keeps 47,267–47,977 PSMs and
28,141–28,857 peptides at a true FDP of 1%; with the *Arabidopsis* entrapment (T/E 2.228),
0.99–1.05% and 49,965–50,185 PSMs. 114–117 s wall, 569–636 CPU-s, 4.8 GB. Comet 2025.01 was run
on the same VM with matched parameters and the same Percolator rescoring on 2026-10-05.

**TMT-labelled, low-res.** One ion-trap CID run of a TMT 10-plex sample against human + yeast
reviewed proteins (26,483 sequences): 13,443–13,516 PSMs and 11,852–11,922 peptides at
`q ≤ 0.01`; 40–47 s, 200–216 CPU-s, 1.7 GB. With a peptide-level entrapment (T/E 1.073) the
nominal 1% is 0.67–0.76% true FDP (13,325–13,371 PSMs at a true FDP of 1%); with *Arabidopsis*
(T/E 2.021), 0.83–1.05% (13,161–13,388).

**TMT-labelled, high-res (TKO hrMS2).** `Fig2_m23142_hrMS2` of PXD017803 (TMTpro, Orbitrap
MS2), `--protocol TMT`, against reviewed yeast + cRAP with a peptide-level entrapment (T/E 1.039):
6,628–6,686 PSMs at `q ≤ 0.01`, 0.98–1.31% true FDP, 6,607–6,662 PSMs at a true FDP of 1%.

**Low-res LFQ (UPS1).** One run of the UPS1 standard in yeast, searched against yeast plus an
E. coli entrapment proteome: 16,049–16,242 PSMs and 4,123–4,183 peptides at `q ≤ 0.01`;
30–48 s, 162–166 CPU-s, 1.3 GB. Against the natural *E. coli* entrapment (T/E 2.351) the
nominal 1% is **2.07–2.19% true FDP**, and andes keeps 14,457–14,605 PSMs at a true FDP of 1%.
With a peptide-level shuffle of the same targets instead (T/E 1.039) it reads 0.91–0.98%
(16,199–16,206 PSMs at 1%). The two entrapments disagree on this low-res set; the natural one is
the established reference here. Java MS-GF+ (v20240326, earlier session, not re-run) reported
15,904 PSMs on this file.

**Phospho-enriched.** One single-shot EasyPhos mouse liver file
(`20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw`, as mzML), not in any model's
training data, 80 Da index slices. Against mouse reviewed proteins plus the whole-protein shuffled
entrapment that the Comet row used: 38,596–38,629 PSMs (27,693–27,724 phospho-bearing) and
20,422–20,461 peptides, 337 s, 2,160 CPU-s, 7.9 GB. Comet 2025.01 reported 33,984 at 1.77% true
FDP on 2026-10-06 (on the whole-protein shuffle, which can read low). andes's true FDP on three
peptide-level shuffles (T/E 1.061) is 0.80–0.99%, with 38,604–38,843 PSMs at a true FDP of 1%;
on *Arabidopsis* (T/E 1.403), 1.22–1.36% and 38,298–38,412. Phosphosite localisation is
**not** scored. Against the depositors' MaxQuant identifications (PEP ≤ 0.01, 23,563 scans)
andes covered 82.4% with 95.3% same-scan backbone agreement on 2026-10-07 (`40774aca`); not
re-measured.

**timsTOF DDA-PASEF.** Native `.d` input, 20 ppm precursor tolerance, 2 missed cleavages,
`--top-n 1`, 80 Da index slices, Carbamidomethyl C fixed, Oxidation M and protein N-terminal
Acetyl variable, against reviewed human proteins with a reviewed yeast entrapment proteome
(natural, T/E 3.845) or a peptide-level shuffle of the human targets (T/E 1.064).

| run | PSMs at `q ≤ 0.01` (peptide-level) | true FDP | PSMs at a true FDP of 1% | natural yeast: PSMs / true FDP / at 1% |
|---|---:|---:|---:|---|
| PXD053748 `..._4363` (plasma, 30 min) | 14,336–14,393 | 0.62–0.63% | 14,772–14,877 | 14,923–14,959 / 0.94–1.14% / 14,788–15,033 |
| PXD053748 `..._4367` (plasma, 50 min) | 17,491–17,550 | 0.87–0.90% | 17,644–17,660 | 18,081–18,136 / 0.78–0.91% / 18,379–18,434 |
| PXD064202 `..._2442` (astrocytes) | 34,210–34,249 | 0.98–1.02% | 34,216–34,222 | 35,105–35,150 / 1.53–1.55% / 33,744–33,959 |
| PXD064202 `..._2403` (astrocytes) | 55,590–55,609 | 0.98–1.03% | 55,561–55,655 | 56,801–56,836 / 1.29–1.32% / 56,028–56,088 |

On the astrocyte runs the natural yeast entrapment reads above nominal (1.29–1.55%) where the
peptide-level shuffle reads 0.98–1.03%. 43–64 s per run, 2.8–3.9 GB.

**Glyco.** Mouse liver glycoproteome data (PXD005553), against mouse proteins plus an entrapment,
with `--glyco --glyco-species mouse --decoy-strategy sequon`. Two published reference
identification sets for the same spectra ship in [`glyco/truth/`](glyco/truth/): pGlyco2 (the
depositors') and MSFragger-Glyco (a later reanalysis, PXD031032). "Confirmed" is the share of a
reference's spectra that andes identifies with the same peptide.
- *Quick tier:* the `MouseLiver-Z-T-1` fraction as converted mzML, three peptide-level shuffled
  databases × Percolator seeds 42, 1 and 2: 6,992–7,072 glycoPSMs at 1%, true FDP 0.61–1.35%,
  pGlyco2 88.2–88.4% and MSFragger-Glyco 89.3–89.7% confirmed; 5.8–6.6 min and 5.4–5.6 GB on 8
  threads. Each run has only 20–45 entrapment hits, so the spread is mostly counting noise.
  Reading the `.raw` natively, as the commands below do, can give somewhat different counts.
- *Deep tier:* all five fractions, pooled before Percolator, five seeds. Measured 2026-09-05 at
  `14818d3e` on an older, smaller glycan list and converted mzML; not re-measured since.
- The open problem is selection, not generation: almost every reference spectrum produces a
  candidate row, but in the deep tier about a fifth of pGlyco2's spectra lose the per-scan
  choice.

**Label-free quantification with match between runs (ProteoBench).**
- *Data and measurement:* the two ProteoBench "Quant LFQ DDA ion" modules, six runs each, where
  the human, yeast and E. coli ratios between conditions A and B are known: QExactive
  (`LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}`, PXD028735, Q Exactive HF-X) and
  Astral (`LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}`, PXD070049). Measured
  2026-10-23 at `2051616b` on a cluster node (8 threads per run), andes defaults throughout
  (§3, `--predict rt`), Percolator 3.7.1 through `--rescore` (seed 42). Scored with the official
  `proteobench` package 0.18.6 ("Custom" format, converted by
  [`proteobench/andes_to_proteobench.py`](proteobench/)).
- *Which quantity is scored:* `quantity_normalized` of the MBR table, which removes each run's
  loading level and its drift along the gradient, unless the table says raw. The raw `quantity`
  is what andes hands to quantms.
- *Results* at a minimum of 3 runs; precursors / median |ε| / E. coli |ε| / CV:

  | module | transfer q ≤ 0.01, normalised | raw quantity | transfer q ≤ 0.05, normalised | identified only, normalised |
  |---|---|---|---|---|
  | QExactive | **65,012 / 0.163 / 0.275 / 0.118** | 0.228 / 0.323 / 0.246 | 71,243 / 0.176 / 0.442 / 0.132 | 52,037 / 0.154 / 0.189 / 0.105 |
  | Astral | **37,352 / 0.193 / 0.258 / 0.113** | 0.213 / 0.256 / 0.115 | 42,656 / 0.211 / 0.343 / 0.130 | 30,759 / 0.194 / 0.243 / 0.103 |

  At least 1 / 6 runs (transfer q 0.01): QExactive 79,126 / 40,258, Astral 47,774 / 22,033.
  The opt-in `--predict` settings `rt,ms2` and `ms2,peptdeep-rt,deeplc-rt` (AlphaPeptDeep and
  DeepLC models, `--predictor-models`, ONNX Runtime backend) were measured on 2026-10-21 at
  `a27a72c`, before the extraction change of #240: on QExactive 67,906 and 68,077 precursors
  against 65,391 for the default on that build, on Astral 34,413 and 33,987 against 33,876.
- *By abundance* (transfer q 0.01, normalised; quartiles of each module's own mean intensity,
  |ε| / E. coli |ε| / CV): QExactive Q1 0.292 / 0.535 / 0.198, Q2 0.175 / 0.273 / 0.130,
  Q3 0.140 / 0.200 / 0.104, Q4 0.117 / 0.164 / 0.078; Astral Q1 0.317 / 0.415 / 0.191, Q2 0.217 /
  0.320 / 0.136, Q3 0.164 / 0.219 / 0.103, Q4 0.135 / 0.185 / 0.070.
- *Cost* (six runs, CPU seconds): QExactive search with prediction and Percolator 5.8k, library
  under 10, extract 0.8k, finalize 0.4k per transfer threshold; Astral 6.1k, under 10, under 0.1k
  and 0.2k. Each search takes 3–4 minutes of wall time, each extract under 20 seconds (QExactive)
  or about 2 seconds (Astral), finalize under a minute.
- *Error control* is measured on the QExactive runs with the HYE database plus a peptide-level
  shuffled entrapment (T/E 1.075), andes defaults, Percolator 3.7.1 `-Y` seeds 42, 1 and 2
  outside andes. PSM FDP at q 0.01 is 0.83–0.85%; precursors quantified in at least 3 runs at
  transfer q 0.01 are 0.20–0.22% entrapment-only (0.26–0.30% at 0.05). Of the accepted
  transfers, the entrapment measures a false-donor share of 0.24–0.33% at transfer q 0.01 and
  0.37–0.49% at 0.05; the identity-decoy estimate is 0.33–0.34% and 0.45–0.47%. The estimate is
  at or above the measured share at 0.01 in all three seeds and below it at 0.05 in one (0.45%
  against 0.49%), so it tracks the false-donor share but is not an upper bound. The wrong-peak
  part has no entrapment check.

**Comparison with the ProteoBench leaderboard.** Public results of both modules
(`Results_quant_ion_DDA` and `Results_quant_ion_DDA_Astral`, both last updated 2026-07-16; read
2026-10-21; scored by `proteobench` 0.17.0). For each tool the submission with the most
precursors and, when different, the one with the lowest median |ε|, at a minimum of 3 runs.
The last column says what limits the error of transferred (match-between-runs) values, taken
from the submission's parameters, comment and the tool's own output; "not checked" means the
submission's parameter file could not be retrieved (the ProteoBench data server was unreachable
on 2026-10-21).

*QExactive:*

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | error control on transferred values |
|---|---|---:|---:|---:|---:|---|
| andes, default (2026-10-23) | yes | 65,012 | 0.163 | 0.275 | 0.118 | PSM q ≤ 0.01; transfer q ≤ 0.01 (decoy-estimated wrong peaks and wrong donors) |
| andes, default, transfer q ≤ 0.05 | yes | 71,243 | 0.176 | 0.442 | 0.132 | transfer q ≤ 0.05 |
| andes, default, raw quantity | yes | 65,012 | 0.228 | 0.323 | 0.246 | as the default |
| FragPipe 23.0, most precursors | yes | 63,701 | 0.200 | 0.306 | 0.156 | identifications 1%; transfers at IonQuant MBR ion FDR 0.01 (submitted workflow) |
| FragPipe 23.0, lowest \|ε\| | yes | 54,374 | 0.190 | 0.257 | 0.154 | as above |
| MaxQuant 1.5.3.30, most precursors | yes | 44,299 | 0.210 | 0.386 | 0.219 | PSM and protein 1%; transfers not checked |
| MaxQuant 2.4.13.0, lowest \|ε\| | no | 29,737 | 0.192 | 0.263 | 0.199 | no transfers |
| i2MassChroQ 1.0.16, most precursors | yes | 74,697 | 0.222 | 0.460 | 0.132 | identifications about 1%; no transfer error setting in the submitted parameters |
| i2MassChroQ 1.2.0, lowest \|ε\| | no | 49,433 | 0.198 | 0.263 | 0.120 | no transfers |
| AlphaPept 0.5.0, most precursors | yes | 37,760 | 0.219 | 0.281 | 0.165 | PSM and protein 1%; transfers not checked |
| AlphaPept 0.5.3, lowest \|ε\| | no | 36,176 | 0.214 | 0.255 | 0.165 | no transfers |
| Sage 0.15.0-beta.2 | yes | 81,933 | 0.388 | 0.877 | 0.393 | peptide q ≤ 0.01 on identifications only; every charge state 2–4 of those peptides is integrated in every run and its peak q-value is not applied |
| MSAngel 2.2.10 | yes | 60,206 | 0.216 | 0.514 | 0.227 | PSM 1%; transfers not checked |
| ProlineStudio 2.3.0-SNAPSHOT | yes | 55,639 | 0.210 | 0.508 | 0.226 | PSM 1%; transfers not checked |
| WOMBAT 0.9.11, most precursors | yes | 55,943 | 0.211 | 0.692 | 0.235 | identifications 1%; transfers not checked |
| WOMBAT 0.9.11, lowest \|ε\| | yes | 41,016 | 0.195 | 0.394 | 0.213 | as above |
| PEAKS 13 | no | 72,741 | 0.289 | 0.741 | 0.216 | no transfers |
| quantms v1.3.1dev (OpenMS ProteomicsLFQ) | no | 33,040 | 0.209 | 0.387 | 0.224 | no transfers |
| quantms (OpenMS ProteomicsLFQ), our run (2026-10-17, its default) | no | 34,677 | 0.210 | 0.401 | 0.225 | no transfers |
| quantms (OpenMS ProteomicsLFQ), our run, its feature-linking MBR | yes | 39,257 | 0.221 | 0.563 | 0.234 | feature-quality score; no transfer error rate |

andes's rows use its own normalisation (`quantity_normalized`: median |ε| 0.163, CV 0.118). Its
raw quantity reads |ε| 0.228 and CV 0.246, less accurate than most submissions; normalisation
differs between tools, so the |ε| column does not compare extraction alone. Among the
submissions, FragPipe 23.0 reaches median |ε| 0.190 at 54,374 precursors, and i2MassChroQ 1.2.0
a CV of 0.120 without match between runs. Six submissions report more precursors than andes
(Sage, i2MassChroQ's four largest and PEAKS 13), at a median |ε| of 0.214–0.388. The two quantms
runs marked "our run" are ours (2026-10-17, its development version, raw intensities).

*On the precursors andes shares* (default, transfer q 0.01; the submission's values | andes raw
| andes normalised, median |ε| / E. coli |ε| / CV):

| submission | shared | submission | andes raw | andes normalised | only in the submission |
|---|---:|---|---|---|---|
| FragPipe 23.0, most precursors | 51,940 | 0.190 / 0.273 / 0.155 | 0.211 / 0.280 / 0.241 | 0.147 / 0.216 / 0.109 | 11,761 at 0.266 / 0.725 |
| FragPipe 23.0, lowest \|ε\| | 48,916 | 0.187 / 0.243 / 0.155 | 0.207 / 0.269 / 0.241 | 0.144 / 0.204 / 0.108 | 5,458 at 0.214 / 0.490 |
| i2MassChroQ 1.0.16, most precursors | 51,156 | 0.208 / 0.353 / 0.125 | 0.214 / 0.300 / 0.243 | 0.153 / 0.245 / 0.113 | 23,541 at 0.266 / 0.935 |
| i2MassChroQ 1.2.4, 57,787 precursors | 50,058 | 0.202 / 0.352 / 0.123 | 0.214 / 0.302 / 0.243 | 0.152 / 0.246 / 0.112 | 7,729 at 0.270 / 0.835 |
| Sage 0.15.0-beta.2, our reproduction | 55,218 | 0.311 / 0.628 / 0.321 | 0.220 / 0.303 / 0.245 | 0.158 / 0.254 / 0.116 | 25,139 at 0.685 / 1.349 |

Sage's per-submission table was not retrievable, so its row uses a reproduction (2026-10-21): the
same Sage version on the same six mzML files with the submission's stated parameters and
per-charge output gives 80,357 precursors at |ε| 0.388, E. coli 0.878 and CV 0.394 (submission:
81,933, 0.388, 0.877, 0.393). Most of its precursors that andes does not report are charge
states that Sage never identified or are above its own peak q-value of 0.01, which it reports but
does not apply.

*Astral:*

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | error control on transferred values |
|---|---|---:|---:|---:|---:|---|
| andes, default (2026-10-23) | yes | 37,352 | 0.193 | 0.258 | 0.113 | PSM q ≤ 0.01; transfer q ≤ 0.01 (decoy-estimated) |
| andes, default, transfer q ≤ 0.05 | yes | 42,656 | 0.211 | 0.343 | 0.130 | transfer q ≤ 0.05 |
| andes, default, raw quantity | yes | 37,352 | 0.213 | 0.256 | 0.115 | as the default |
| FragPipe 23.0, most precursors | yes | 39,057 | 0.212 | 0.250 | 0.127 | identifications 1%; transfers not checked |
| FragPipe 22.1-build02, lowest \|ε\| | yes | 29,105 | 0.181 | 0.405 | 0.078 | as above |
| MaxQuant 2.6.7.0 | yes | 20,283 | 0.160 | 0.158 | 0.092 | PSM and protein 1%; transfers within a condition only (submission comment), not checked |
| i2MassChroQ 1.2.6, most precursors | yes | 36,469 | 0.197 | 0.417 | 0.095 | identifications about 1%; transfers not checked |
| i2MassChroQ 1.2.9, lowest \|ε\| | no | 24,803 | 0.178 | 0.245 | 0.077 | no transfers |
| Sage 0.14.6 | yes | 30,319 | 0.167 | 0.338 | 0.070 | peptide q ≤ 0.01 on identifications only; peak q-value not applied (as in QExactive) |
| DIA-NN 2.6.0, most precursors | yes | 93,237 | 0.277 | 0.691 | 0.159 | run and global precursor q ≤ 0.01 (submission comment) |
| DIA-NN 2.6.0, lowest \|ε\| | yes | 16,363 | 0.123 | 0.168 | 0.079 | as above, plus MS1-level q ≤ 0.01 |

On Astral, MaxQuant 2.6.7.0 (20,283 precursors, median |ε| 0.160), Sage 0.14.6 (30,319, 0.167,
CV 0.070) and i2MassChroQ 1.2.6 (36,469, 0.197, CV 0.095) report fewer precursors than andes
(37,352, 0.193) and FragPipe 23.0 more (39,057, 0.212); DIA-NN 2.6.0 reports 93,237 at 0.277, and
its 16,363-precursor submission has the lowest median |ε| (0.123). andes's CV (0.113) is above
most submissions. Normalisation changes little here (raw |ε| 0.213, CV 0.115). Shared-precursor
comparisons need the per-submission tables, which could not be retrieved; they are not given.

**TMT quantification.**
- *Data and measurement:* PXD017803 (yeast triple knockout, TMTpro 16-plex; two files each of
  hrMS2, SPS-MS3 and RTS-MS3) and PXD007683 (human + yeast TMT11 with known yeast ratios, 11
  SPS-MS3 fractions, ion-trap CID MS2), as mzML, with the kit impurity sheet. andes defaults,
  including the quantification gate (precursor purity ≥ 0.7 and Percolator PEP ≤ 0.05;
  §10 of DOCS.md), `--rescore`, 16 threads; measured 2026-10-23 at `2051616b`.
- *Extraction:* reporter extraction and impurity correction match OpenMS IsobaricAnalyzer
  exactly on the same scans.
- *Knockout interference* (median ratio of knockout to non-knockout channel intensity over the
  PSMs of the three knocked-out proteins, lower is better; hrMS2 only scans whose window starts
  below 126 Th):

  | acquisition | PSMs at 1% | PSMs quantified (gate passed) | interference | without the gate | purity ≥ 0.9 |
  |---|---:|---:|---:|---:|---:|
  | SPS-MS3 | 9,816 | 8,626 | 0.064 | 0.128 | 0.032 |
  | RTS-MS3 | 14,307 | 10,605 | 0.023 | 0.033 | 0.015 |
  | hrMS2 | 13,464 | 3,618 | 0.182 | 0.185 | 0.134 |

  Few knockout PSMs stand behind each value (36–91), so read them as a direction. OpenMS
  precursor purity at ≥ 0.9 gave 0.139 on SPS-MS3 (2026-10-10). On hrMS2 only 3,618 of the
  13,464 PSMs have a reporter scan whose range covers the reporter ions.
  `--tmt-unmix` on hrMS2 sets the median knockout ratio to 0.000 at a background CV of 0.124
  (0.080 without).
- *Known ratios (PXD007683):* 149,133 rank-1 target PSMs at 1% over the 11 fractions, 128,177 quantified;
  yeast log2 A/C 1.591 (expected 1.585), A/B 1.192, B/C 0.396, with protein-level MAE 0.212,
  0.228 and 0.215 (without the gate: A/C 1.464, MAE 0.271). 9.6 minutes and 17.5 GB in one
  process.

**Large study (PXD042233, 500 runs).** Label-free DDA from Q Exactive HF, HF-X and Exploris 480,
`.raw` read natively, human reviewed proteins + contaminants, 5 ppm, 2 missed cleavages,
`--top-n 1`, Carbamidomethyl C fixed, Oxidation M and protein N-terminal Acetyl variable. Each
run is searched as its own job (8 threads, `--quant-dir`); the per-run PINs share one header
(it depends only on the command line) and are concatenated into one Percolator run (16 threads,
`--subset-max-train 1000000`); then `quant library`, `quant extract` per run (8 threads) and
`quant finalize` (16 threads). Measured 2026-10-23 at `2051616b` on cluster nodes.

| stage | wall | CPU | peak memory |
|---|---|---:|---:|
| search, 500 jobs | 0.7 min per run (at most 1.0) | 19.5 h | 3.6 GB |
| PIN concatenation (53.1 million rows) | 17 s | — | — |
| Percolator | 12.9 min | 0.7 h | 26.0 GB |
| `quant library` | 6.6 min | 0.1 h | 5.6 GB |
| `quant extract`, 500 jobs | 2.1 min per run (two runs 45 and 72 min) | 121.5 h | 12.7 GB (median 4.7) |
| `quant finalize`, transfer q 0.01, with the QPX views and protein inference | 25.8 min | 4.1 h | 45.2 GB |
| `quant finalize`, transfer q 0.05, without QPX views | 18.0 min | 3.9 h | 11.7 GB |

Outputs: 30.6 million PSMs at `q ≤ 0.01`. At transfer q 0.01, 194,370 accepted transfers
(175,087 `transfer_ms1`, 19,283 `transfer_ms2`) and 294,362 / 243,895 / 127,315 precursors
quantified in at least 1 / 3 / 50 runs (identified only: 294,362 / 243,730 / 126,741); 9,208
protein groups at protein q ≤ 0.01. At transfer q 0.05, 4,622,956 accepted transfers and
294,362 / 250,847 / 141,734 precursors. Extract is 83% of the CPU, and most of its transfer fits
are rejected (121 million `below_threshold` rows at q 0.01).

## 3. How to reproduce

You need `andes` on `PATH`, Docker (for the pinned Percolator image), `curl` and `python3`. Build
andes with `--features thermo` and install the .NET 8 runtime to read `.raw` directly; this is
preferred over converting (see [`reproduce/README.md`](reproduce/README.md)).

### Standard sets: Astral, TMT, UPS1 (tested scripts)

```bash
DATA=~/andes-bench                            # ~10 GB for these three
cd docs/benchmarks/reproduce
./build_databases.sh "$DATA"                  # UniProt download, ~1 min; prints the UPS1 T/E factor
./fetch_spectra.sh   "$DATA" astral tmt ups1  # 4.8 GB from PRIDE
./run.sh             "$DATA" astral tmt ups1  # search + Percolator; prints wall, PSMs, entrapment hits
```

`run.sh` uses 8 threads (`THREADS=` to change), Percolator 3.7.1 `--seed 42 -Y`, and the
`configs/` mods files (Astral: `astral_mods.txt`, `--precursor-tol 10ppm`; TMT: `mods-tmt.txt`).
UniProt is versioned and the original HYE FASTA is no longer served, so counts will differ
slightly from the table; the scripts print a sha256 and sequence count for each input.

### Phospho (commands)

```bash
GLYCO=1 ./build_databases.sh "$DATA"          # also writes databases/mouse_entrap.fasta and prints its T/E
./fetch_spectra.sh "$DATA" phospho            # 3.7 GB
andes --spectrum "$DATA"/phospho/20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw \
      --database "$DATA"/databases/mouse_entrap.fasta \
      --mods ../configs/mods-phospho.txt --threads 8 \
      --fragment-index-slice-da 80 --output-pin phospho.pin
docker run --rm -v "$PWD":/r quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 percolator \
      --seed 42 -Y --results-psms /r/phospho.t.psms --decoy-results-psms /r/phospho.d.psms /r/phospho.pin
```

Count rows at `q-value ≤ 0.01` and their entrapment hits with `python3 ../entrap_common.py count
phospho.t.psms "$DATA"/databases/mouse_entrap.fasta 0.01` (FDP = hits / PSMs × (1 + T/E), with the
T/E that `build_databases.sh` printed), and peptides carrying `79.96`. Expect about 6 minutes on 8 threads
(337 s on the benchmark VM). The slice width is pinned because it sets peak memory (7.9 GB at 80 Da).

### Glyco (commands; scripts in [`glyco/`](glyco/))

```bash
GLYCO=1 ./build_databases.sh "$DATA"
./fetch_spectra.sh "$DATA" glyco-mouse        # 5 fractions, 12.6 GB; the quick tier needs only T-1
cd ..                                         # docs/benchmarks
DB="$DATA"/databases/mouse_entrap.fasta
# --decoy-strategy sequon is required: other generated decoys carry fewer sequons than their targets.
for f in MouseLiver-Z-T-1 MouseLiver-Z-T-2 MouseLiver-Z-T-3 MouseLiver-Z-T-4 MouseLiver-Z-T-5; do
  andes --spectrum "$DATA"/glyco-mouse/$f.raw --database "$DB" \
        --glyco --glyco-species mouse --decoy-strategy sequon \
        --threads 8 --output-pin $f.pin                      # writes $f.glyco.pin
done
python3 glyco/pool_pins.py MouseLiver-Z-T-*.glyco.pin > pooled.pin   # pool BEFORE Percolator
docker run --rm -v "$PWD":/r quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 percolator \
      --seed 42 -Y --results-psms /r/pooled.t.psms --decoy-results-psms /r/pooled.d.psms /r/pooled.pin
python3 glyco/eval_yield.py  pooled.pin pooled.t.psms
python3 glyco/eval_entrap.py pooled.pin pooled.t.psms 0.01 "$DB"
for ref in glyco/truth/*_mouse_liver.tsv.gz; do              # the two reference sets
  python3 glyco/score_vs_truth.py "$ref" pooled.pin pooled.t.psms
done
```

For the quick tier, search `MouseLiver-Z-T-1` alone, run Percolator on its `.glyco.pin` without
pooling, and pass `--run MouseLiver-Z-T-1` to `score_vs_truth.py` (a single-file run writes
SpecIds without a file name). About 6 minutes and 6 GB on 8 threads. Pool the deep tier's
fractions: one fraction alone has too few glyco decoys for a stable 1% threshold.

### Label-free quantification, ProteoBench (commands; converter tested)

Inputs: the six runs of a module and the HYE FASTA linked from the ProteoBench module
documentation (`ProteoBenchFASTA_MixedSpecies_HYE.fasta`, 31,889 entries).
- *QExactive:* `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.raw` of PXD028735 (about 21 GB).
- *Astral:* `LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw` of PXD070049 (about 15 GB).

The measured runs read mzML peak-picked from these files. `R` is each run's stem. Every flag not
shown is the andes default (20 ppm precursor tolerance, one missed cleavage, Carbamidomethyl C
fixed, Oxidation M and protein N-terminal Acetyl variable, `--predict rt`).

```bash
# 1. per run: search, Percolator (--rescore) and the run's MS1, MS2 and PSM tables in q/
andes --spectrum $R.mzML --database ProteoBenchFASTA_MixedSpecies_HYE.fasta \
      --rescore --output-pin $R.pin --quant-dir q/ --threads 8
# 2. once: the library, from every run's Percolator results (one --psms/--decoy-psms per run)
andes quant library --quant-dir q/ --psms $R=$R.percolator.target.psms.txt \
      --decoy-psms $R=$R.percolator.decoy.psms.txt ...
# 3. per run
andes quant extract --quant-dir q/ --run $R
# 4. once per transfer threshold (the default --mbr-fdr is 0.01)
andes quant finalize --quant-dir q/ -o pb
andes quant finalize --quant-dir q/ -o pb05 --mbr-fdr 0.05
# 5. ProteoBench "Custom" table (python3 + pyarrow); --module qe or astral
python3 proteobench/andes_to_proteobench.py --mbr pb.mbr.parquet --quant-dir q/ \
      --module qe --column quantity_normalized -o custom_input.tsv
```

6. Score with `proteobench` 0.18.6: `DDAQuantIonModuleQExactive` (or `DDAQuantIonAstralModule`)
   `.benchmarking(custom_input.tsv, "Custom", ...)` reports minimum 1 to 6 runs in one pass.

The converter keeps the `identified`, `transfer_ms2` and `transfer_ms1` rows of the MBR table,
takes each precursor's peptidoform and proteins from the quant directory, and writes Unimod
names in ProForma. `--column quantity` writes the raw quantity instead of `quantity_normalized`,
and `--evidence identified` writes the identification-only table. Its tests:
`python3 -I tests/test_proteobench_convert.py`.

`andes --spectrum ... --spectrum ... --rescore --lfq --mbr --quant-dir q/` runs the same stages
in one process, with one Percolator over all runs; it was not the configuration measured here.

### TMT quantification, knockout (commands, not a tested script)

The PXD017803 files `Fig2_m23142_hrMS2.raw` and `Fig2_m23145_SPSms3.raw`, reviewed yeast
proteins, and a mods file with Carbamidomethyl C and TMTpro (304.20715) on K and the peptide
N-terminus as fixed modifications, plus Oxidation M as variable. The measured runs read mzML
converted from these files (both files of each acquisition in one process), with 16 threads:

```bash
andes --spectrum Fig2_m23142_hrMS2.raw  --database yeast.fasta --mods mods_tmtpro.txt \
      --precursor-tol 20ppm --protocol TMT --rescore --fdr 0.01 \
      --tmt tmt16 --tmt-level 2 --tmt-correction <kit impurity sheet> --output-tsv hr.tsv
andes --spectrum Fig2_m23145_SPSms3.raw --database yeast.fasta --mods mods_tmtpro.txt \
      --precursor-tol 20ppm --protocol TMT --rescore --fdr 0.01 \
      --tmt tmt16 --tmt-level 3 --tmt-correction <kit impurity sheet> --output-tsv sps.tsv
```

Channel layout: MET6 (P05694) is knocked out in 126–128N, HIS4 (P00815) in 128C–130N and URA2
(P07259) in 130C–132N; 132C–134N is the parental strain. For each PSM of those proteins, divide
the mean of its knockout channels by the mean of the other channels in `<stem>.tmt.parquet`
(rows with `quant_gate_passed`). Report
the median, using only scans whose reporter scan window starts below 126 Th. Repeat on PSMs
with `purity ≥ 0.9`.

## 4. Method notes

- **True FDP by entrapment.** A target PSM that matches only entrapment sequences is false by
  construction: FDP = (entrapment hits / accepted) × (1 + T/E), where T/E is the ratio of
  searchable target to entrapment peptides, measured for each database and never assumed to be
  1. A PSM whose peptide also occurs in a target protein counts as a target: it is entrapment
  only when **all** its proteins are entrapment and its peptide, with I folded to L (the two
  are isobaric), occurs in no target protein sequence. `entrap_common.py` holds this rule
  for `run.sh`, `build_databases.sh` (the T/E factor) and the glyco evaluators.
- **Entrapment design.** Shuffled entrapment is built peptide by peptide
  ([`build_shuffled_entrap.py`](build_shuffled_entrap.py)): each tryptic peptide of each target
  protein is shuffled with its C-terminal residue fixed, as in the peptide-level entrapment
  construction of Wen et al., Nat. Methods 2025. Every entrapment peptide then has a target of
  the same composition and mass. The earlier builder shuffled whole proteins: 89.6% of its
  entrapment peptides on the HYE database had no same-composition target peptide, and on the
  same Astral search it read 0.81% where the peptide-level databases read 0.92–1.02% and a
  foreign proteome (reviewed *Arabidopsis*, T/E 2.228) 0.99% (#225). On phospho it agreed with
  the peptide-level databases (1.05% against 0.95–1.04%). FDPs marked "whole-protein shuffle"
  were measured that way and may read low.
- **Percolator.** Version 3.7.1, `--seed 42 -Y`, on plain FASTA with andes's `XXX_` decoys, for
  every engine compared on the same host. Percolator's q-value has a floor of 1 / (number of
  top targets), so small or sparse runs return few or no identifications at 1%; pool glyco
  fractions for this reason.
- **ProteoBench terms.** "Min runs" is the number of the six runs in which a precursor
  (peptidoform and charge) must be quantified to be counted. Epsilon is the difference between
  the measured log2 A/B ratio and the expected one; the headline is its median absolute value
  over all precursors. CV is the median coefficient of variation within a condition, averaged
  over the two conditions. A
  transfer q-value is an internal quant-quality score for the peak a precursor is matched to in
  a run where it was not identified; it is not an identification FDR.
- **Comparisons.** Identification comparisons are runs on the same host with matched parameters
  and the same rescoring. ProteoBench comparisons are submissions on its public leaderboard.
  Wall times compare only within one host and session; counts are portable.

## 5. Known gaps

- **Label-free quantification rests on two datasets and one Percolator seed** (three for the
  error control), and its best scores use andes's own normalisation (`quantity_normalized`; the
  raw `quantity` is the quantms contract). On Astral andes's CV (0.113) is above most
  submissions, and its match-between-runs entrapment error is measured on QExactive only. The
  leaderboard was scored by `proteobench` 0.17.0 and andes by 0.18.6; several engines' transfer
  error control could not be checked.
- **Entrapment designs can disagree.** On UPS1 the natural *E. coli* entrapment reads the
  nominal 1% as 2.07–2.19% and a peptide-level shuffle as 0.91–0.98%; on the timsTOF astrocyte
  runs natural yeast reads 1.29–1.55% and the shuffle 0.98–1.03%. The 500-run study has no
  entrapment in this measurement.
- **Not re-measured at `2051616b`:** the glyco deep tier, the compared phospho FDP (Comet), the
  `--chimeric` and `--refine` opt-in table, the non-default `--predict` settings and the
  phospho agreement with MaxQuant.
- **Single files.** Phospho is one file, peptide-level only; `--refine` is measured on Astral
  only.
- **Glyco composition on fucose-rich tissue.** On heart and lung the backbone agrees with
  pGlyco2 94–99% of the time but the composition only 68–78% (97% on liver).
- **Not benchmarked:** iTRAQ, phosphosite localisation, timsTOF quantification. Java MS-GF+ was
  not re-run under the current defaults; its counts are historical.

## Layout

    docs/benchmarks/
      README.md      this page
      reproduce/     database, download and run scripts (no hardcoded paths)
      glyco/         glyco harness: pooling, yield, entrapment, reference scoring
      configs/       mods files and per-engine parameter files
      proteobench/   andes MBR output -> ProteoBench "Custom" table (python3 + pyarrow)
      entrap_common.py  shared entrapment rules (ALL-ENTRAP_, I/L rescue, T/E)
      tests/         tests: python3 -I tests/test_entrap.py (stdlib), tests/test_proteobench_convert.py

Bulk spectra and third-party engine binaries stay outside git.
