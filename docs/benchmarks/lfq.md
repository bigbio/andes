# Label-free quantification benchmarks

## ProteoBench

Setup: ProteoBench "Quant LFQ DDA ion", QExactive (PXD028735) and Astral (PXD070049), six runs each;
andes `09e1c00d`, 2026-11-05, scored with `proteobench` 0.18.6; public leaderboard read 2026-10-21.
The summary compares andes with the MaxQuant and FragPipe submissions; every andes view and the whole
leaderboard are in [ADDITIONAL_BENCHMARKS.md](ADDITIONAL_BENCHMARKS.md#proteobench).

- andes' raw `quantity` (what quantms receives) is less accurate than most submissions (QExactive
  median |ε| 0.222, CV 0.245); the default `quantity_normalized` reads 0.161 and 0.117.
- QExactive: at least six submissions report more precursors than andes (up to Sage, 81,933), all
  at a higher median |ε| (0.214–0.388).
- Astral: some submissions report lower error or CV with fewer precursors (Sage 0.167 and CV 0.070
  at 30,319; DIA-NN 0.123 at 16,363); andes' CV (0.112) is above most submissions.

## Truth designs

Setup: eight sets with a known answer, quantified by andes, MaxQuant and FragPipe and read through
one evaluation ([`quant-eval`](quant-eval/README.md); rules in [Methods](methods.md#method-notes)). andes `09e1c00d`
defaults, run 2026-11-05 over the 2026-10-31 searches (`--predict rt,ms2 --chimeric --refine
--max-missed-cleavages 2`); `quant library` (`--mbr-donor-fdr 0.01`), `quant extract`, `quant
finalize` (`--quant-fdr 0.01`), quant-model seeds 42, 1 and 7 (PXD030780: 42); tables show seed 42.
Raw `quantity` column. MBR on = identified rows plus transfers with `transfer_q ≤ 0.01`; MBR off =
identified rows. MaxQuant 2.7.0.0 (Standard DDA, LFQ on, 1% PSM and protein FDR, 2 missed
cleavages, `Intensity` of the evidence table per precursor and run; MBR-off run re-run 2026-11-05,
MBR-on run from 2026-10). FragPipe 21.1 with IonQuant 1.10.12, MBR on (ion FDR 0.01), normalisation
off (`--normalization 0`), 2026-10 outputs, no MBR-off run.

| set | what it is | runs | truth |
|---|---|---:|---|
| PXD009740 | UPS1 at ten levels (170 amol–250 fmol) in constant yeast | 30 | UPS1 amount per run; yeast constant |
| PXD009815 | UPS1 at ten levels (10 amol–50 fmol) in yeast | 40 | as above |
| PXD001819 | UPS1 at nine levels in yeast; level and run order confounded | 27 | as above; counts and accuracy only (every tool reports 1,000+ yeast false positives in DE) |
| QE (ProteoBench) | human, yeast, *E. coli*; Q Exactive HF-X (PXD028735) | 6 | yeast 2:1, *E. coli* 1:4, human 1:1 (A vs B) |
| Astral (ProteoBench) | same mixture, Astral (PXD070049) | 6 | as above |
| PXD030780 | human, yeast, *E. coli*, five mixing levels A–E | 25 | species amounts per group |
| PXD014415 | human and human + yeast runs | 40 | yeast absent from 20 runs |
| PXD016921 | blank run among cell-input runs | 6 | no analyte in the blank |

- Settings are not matched: andes used `--chimeric`, `--refine` and PSM-level FDR; MaxQuant and
  FragPipe used protein-level 1% FDR. With matched settings (±10 ppm, isotope errors 0 to 2, no
  `--chimeric`, no `--refine`, seed 42) andes's precursors fall 31% on PXD009740 and PXD009815 (to
  27,226 and 29,600, level with the other tools), 12% on QE (81,987) and 6% on Astral (51,798),
  and DE true positives stay within seed spread (PXD009740 1,293, PXD009815 831, QE 1,049, Astral
  673). The depth lead on the UPS sets comes mostly from the search settings; the DE results come
  from the quantification. A protein-level 1% filter on andes removes 0.4–1.6% of precursors.
- MaxQuant 2.7 does not remove match rows with its MBR FDR option (it only fills the `Match
  q-value` column), so MaxQuant's MBR-on rows are unfiltered, while andes's transfers are filtered
  at `transfer_q ≤ 0.01`. Filtering MaxQuant at `Match q-value ≤ 0.01` cuts PXD001819 transfers from
  43,373 to 2,550 and blank-run transfers from 85 to 9. No MaxQuant run exists for PXD014415 or
  PXD030780.

### Depth

Cells are precursor-by-run values. Each entry reads precursors / proteins (≥ 2 precursors) /
identified cells + transferred cells.

| set | andes, MBR on | MaxQuant, MBR on | FragPipe, MBR on | result |
|---|---|---|---|---|
| PXD009740 | 39,595 / 3,031 / 573,349 + 332,626 | 29,293 / 2,339 / 464,069 + 175,450 | 25,888 / 2,259 / 414,258 + 141,367 | 🟢 +35% |
| PXD009815 | 42,686 / 3,362 / 798,525 + 463,945 | 27,459 / 2,453 / 586,274 + 230,204 | 28,286 / 2,691 / 587,049 + 193,057 | 🟢 +51% |
| PXD001819 | 10,696 / 1,004 / 146,103 + 65,796 | 9,314 / 881 / 116,299 + 43,373 | 10,406 / 937 / 118,961 + 63,595 | 🟡 +3% |
| QE | 92,828 / 6,797 / 350,579 + 71,126 | 51,609 / 4,810 / 186,178 + 25,465 | 71,377 / 5,938 / 249,172 + 12,401 | 🟢 +30% |
| Astral | 55,246 / 4,901 / 212,633 + 28,812 | 24,745 / 2,855 / 93,669 + 29,072 | 46,347 / 4,334 / 168,481 + 46,343 | 🟢 +19% |
| PXD030780 | 93,566 / 8,602 / 1,047,714 + 560,171 | — | 65,586 / 6,885 / 824,966 + 351,518 | 🟢 +43% |
| PXD014415 | 76,099 / 5,470 / 1,494,110 + 801,869 | — | 61,771 / 5,030 / 1,175,139 + 550,247 | 🟢 +23% |
| blank run | 17,625 / 1,728 / 29,755 + 1,641 | 11,010 / 1,229 / 17,749 + 1,305 | 14,153 / 1,403 / 24,638 + 4,015 | 🟢 +25% |

Result is the precursor gain over the best other tool, each tool at its own search settings (first
bullet above). With andes at matched settings (no `--chimeric`, no `--refine`), measured on four
sets, the gain is −7% on PXD009740, +5% on PXD009815, +15% on QE and +12% on Astral. With MBR off, precursor counts change by less than 0.2% for every tool and the transferred
cells are 0.

### Accuracy and precision

Cells identified by every tool listed for the set. Entries read andes / MaxQuant / FragPipe (MBR on;
andes MBR-on accuracy equals MBR off because the cells are identified cells). Bright decile = top
tenth of predicted intensity; gross error = share of cells off by more than 2x. Protein CV is the
median replicate CV of matched features.

| set | cells off by > 2x (%) | bright-decile gross error (%) | precursor SD, top decile | protein CV (%) | protein CV, MBR off (%): andes / MaxQuant | result |
|---|---|---|---|---|---|---|
| PXD009740 | 2.63 / 1.55 / 4.02 | 3.64 / 1.13 / 2.37 | 0.132 / 0.116 / 0.126 | 8.7 / 10.5 / 13.5 | 10.2 / 10.5 | 🔴 accuracy; 🟢 CV |
| PXD009815 | 1.42 / 0.35 / 0.65 | 5.97 / 0.55 / 0.68 | 0.056 / 0.066 / 0.051 | 4.4 / 7.0 / 6.1 | 5.6 / 7.5 | 🔴 accuracy; 🟢 CV |
| PXD001819 | 1.81 / 0.86 / 1.80 | 1.79 / 0.46 / 3.76 | 0.044 / 0.060 / 0.038 | 5.7 / 8.3 / 5.2 | 6.8 / 7.9 | 🔴 vs MaxQuant; 🔴 CV vs FragPipe (minor) |
| QE | 3.75 / 2.92 / 5.84 | 6.44 / 1.62 / 2.64 | 0.208 / 0.184 / 0.206 | 9.2 / 11.6 / 11.5 | 9.6 / 10.3 | 🔴 accuracy; 🟢 CV |
| Astral | 1.62 / 0.69 / 1.13 | 4.39 / 0.26 / 0.66 | 0.057 / 0.072 / 0.041 | 5.8 / 7.0 / 4.9 | 6.0 / 6.9 | 🔴 accuracy; 🔴 CV vs FragPipe (minor) |
| PXD030780 | 2.63 / — / 1.44 | 7.19 / — / 0.30 | 0.112 / — / 0.096 | 10.1 / — / 10.2 | 12.1 / — | 🔴 accuracy; 🟡 CV |

- Bright-decile gross error is 1.8–7.2% for andes against 0.3–1.6% for MaxQuant (MBR on) and
  0.3–3.8% for FragPipe; accuracy is behind on every set. Protein CV is 4.4–10.1% for andes,
  7.0–11.6% for MaxQuant and 4.9–13.5% for FragPipe. MaxQuant MBR-off values are close to its MBR-on
  ones (gross error up to 1.65%).

### Differential expression (BH 5% per pair)

Entries are true positives / false calls on the common-tested family, sign-aware (wrong-sign
calls, counted as false, in brackets). Seed 42. Across seeds 42, 1 and 7 andes MBR-on gives
1,288–1,302 / 305–350 on PXD009740, 841–847 / 153–188 on PXD009815, 1,047–1,050 / 110–120 on QE and
681–685 / 296–301 on Astral.

| set | andes, MBR on | andes, MBR off | MaxQuant, MBR on | MaxQuant, MBR off | FragPipe, MBR on | result |
|---|---|---|---|---|---|---|
| PXD009740 | 1,300 / 305 (1) | 1,068 / 180 (0) | 1,304 / 129 (0) | 1,047 / 188 (0) | 1,237 / 302 (0) | 🟡 true positives; 🔴 2.4x MaxQuant false calls (moderate) |
| PXD009815 | 844 / 153 (7) | 640 / 76 (0) | 844 / 92 (2) | 541 / 54 (2) | 849 / 301 (2) | 🟡 true positives; 🔴 1.7x MaxQuant false calls (moderate) |
| PXD001819 | 375 / 3,296 (0) | 266 / 1,833 (0) | 358 / 939 (0) | 240 / 861 (0) | 359 / 1,977 (0) | 🟡; not a decision set |
| QE | 1,047 / 110 (0) | 1,007 / 62 (0) | 857 / 37 (0) | 824 / 43 (0) | 946 / 43 (1) | 🟢 +22% over MaxQuant, +11% over FragPipe; 🔴 3.0x false calls (moderate) |
| Astral | 685 / 296 (1) | 640 / 251 (0) | 650 / 210 (2) | 465 / 141 (0) | 689 / 389 (1) | 🟡 true positives; 🔴 1.4x MaxQuant false calls (minor) |
| PXD030780 | 17,692 / 1,374 (23) | 14,416 / 915 (29) | — | — | 18,227 / 1,751 (22) | 🟡 |

MaxQuant matches filtered at `Match q-value ≤ 0.01` (true positives / false calls; PXD009740,
PXD009815, QE, Astral):

| tool | PXD009740 | PXD009815 | QE | Astral | result |
|---|---|---|---|---|---|
| andes, MBR on | 1,156 / 272 | 763 / 151 | 983 / 105 | 568 / 297 | 🟢 +5% to +15% true positives; 🔴 1.9–2.5x false calls (moderate) |
| MaxQuant, matches filtered | 1,077 / 131 | 725 / 62 | 851 / 42 | 518 / 154 | |

Spread of unchanged proteins within a condition (tool ÷ andes; PXD009740 / PXD009815 / QE): MaxQuant
1.15 / 1.55 / 1.26, FragPipe 1.40 / 1.30 / 1.25. True positives / false calls by minimal fold change
(seed 42 at `e841aae0` (2026-11-08, the #333 per-run q-values), measured 2026-11-09; the table above is
the mean over seeds at `09e1c00d`, e.g.
PXD009740 reads 1,287 / 320 here and 1,300 / 305 above):

| set | tool | all | \|log2 FC\| ≥ 0.3 | ≥ 0.5 | ≥ 1.0 |
|---|---|---|---|---|---|
| QE | andes | 1,049 / 113 | 1,049 / 83 | 1,046 / 11 | 476 / 0 |
| | MaxQuant | 857 / 37 | 857 / 35 | 852 / 18 | 303 / 2 |
| | FragPipe | 946 / 43 | 946 / 37 | 945 / 20 | 680 / 8 |
| PXD009740 | andes | 1,287 / 320 | 1,287 / 310 | 1,287 / 167 | 1,271 / 27 |
| | MaxQuant | 1,292 / 129 | 1,292 / 128 | 1,292 / 74 | 1,238 / 25 |
| | FragPipe | 1,224 / 302 | 1,224 / 302 | 1,224 / 259 | 1,220 / 86 |
| PXD009815 | andes | 853 / 183 | 852 / 60 | 846 / 39 | 800 / 21 |
| | MaxQuant | 842 / 92 | 842 / 66 | 838 / 36 | 778 / 6 |
| | FragPipe | 859 / 301 | 859 / 184 | 858 / 102 | 844 / 22 |

- At `|log2 FC| ≥ 0.5`, andes has 167 false calls against MaxQuant's 74 on PXD009740 (🔴 2.3x), 39
  against 36 on PXD009815 and 11 against 18 on QE; Astral was not recomputed. Still behind on
  PXD009815 at large fold changes (21 against 6 at ≥ 1.0) and in wrong-sign calls (18 against 2 at
  seed 42 on `e841aae0`; 7 against 2 in the table above).
- Why: andes's protein values are tighter, so small real shifts of "unchanged" proteins reach
  significance; the truth labels assume no change at all, which stops holding below about 1.4-fold
  (|log2 FC| 0.5). The calls are not made by transferred values: from identified cells alone the
  called proteins keep 80–95% of their fold change. On PXD009740, 296 of 320 false calls are against
  the highest UPS1 level; without that level andes has 955 true positives and 24 false calls against
  MaxQuant's 953 and 31.
- Match between runs raises true positives for every tool that offers it (andes 1,068 to 1,300 on
  PXD009740, 1,007 to 1,047 on QE, 640 to 685 on Astral, 14,416 to 17,692 on PXD030780; MaxQuant
  1,047 to 1,304, 824 to 857, 465 to 650). andes's false calls rise with it (180 to 305, 62 to 110);
  MaxQuant's fall on PXD009740 and QE (188 to 129, 43 to 37) and rise on Astral and PXD009815.

### UPS1 detection at the lowest three levels

The three lowest levels of each set (PXD009740 170, 382, 859 amol; PXD009815 10, 50, 100 amol;
PXD001819 50, 125, 250 amol), MBR on. A row is one UPS1 precursor in one run. The table has two
units:

- precursors: at each level, the distinct precursors with at least one row within 2x of their own
  high-level reference, summed over the three levels. "At FragPipe's row count" is the expected
  number of such precursors if the tool reported only as many rows as FragPipe, drawn at random.
- rows: the rows at the three levels, the share of them within 2x and their median error. The
  share is over rows, so it is not the first column divided by the second.

| set | tool | precursors within 2x | rows | rows within 2x (% of rows) | median \|log2 err\| of rows | precursors within 2x at FragPipe's row count | result |
|---|---|---:|---:|---:|---:|---:|---|
| PXD009740 | andes | 502 | 1,600 | 60.3 | 0.84 | 436.7 | 🟢 vs MaxQuant; 🔴 vs FragPipe at equal rows (minor) |
| | MaxQuant | 417 | 1,642 | 38.9 | 1.20 | 340.7 | |
| | FragPipe | 458 | 1,167 | 71.1 | 0.64 | 458.0 | |
| PXD009815 | andes | 24 | 565 | 7.6 | 4.75 | 22.1 | 🔴 2x fewer than FragPipe (minor); 🟡 vs MaxQuant |
| | MaxQuant | 23 | 603 | 6.3 | 4.63 | 19.1 | |
| | FragPipe | 50 | 443 | 15.6 | 3.73 | 50.0 | |
| PXD001819 | andes | 0 | 53 | 0.0 | 6.00 | 0.0 | 🟡 none detects |
| | MaxQuant | 1 | 68 | 1.5 | 4.25 | 0.5 | |
| | FragPipe | 3 | 29 | 10.3 | 3.13 | 3.0 | |

Per-level table: [ADDITIONAL_BENCHMARKS.md](ADDITIONAL_BENCHMARKS.md#ups1-per-level).

- With MBR off: PXD009740 andes 206 (546 rows, 62.6% within 2x, median error 0.84) against MaxQuant
  111 (360 rows, 48.3%, 1.01); PXD009815 andes 3 against MaxQuant 0; PXD001819 andes 0, MaxQuant 0.
- At equal rows FragPipe is ahead on both UPS1 sets (458 against andes 437; 50 against 22) and
  andes is ahead of MaxQuant on PXD009740 (437 against 341) and level on PXD009815 (22 against 19).
- At 10 and 50 amol every tool over-reads: median error +6.5 to +7.4 log2 at 10 amol and +4.0 to
  +4.7 at 50 amol (a factor of 2^4 to 2^7).

### Absence: false transfers and the blank run

PXD014415: yeast rows among the transferred rows of the 20 human-only runs. PXD016921: rows
reported in the blank run (andes on main of 2026-11-09, after #333 `e841aae0`, with per-run q-values). Classes (definition in
[Methods](methods.md#method-notes)) are enriched / not specific / only
absent / no reference.

| metric | andes | MaxQuant 2.7 | FragPipe | result |
|---|---|---|---|---|
| yeast rows, identified (MBR off or on) | 451 | — | 208 | ⚪ |
| yeast rows, transferred; share of transferred rows | 769; 0.20% of 389,423 (seeds 0.19–0.20%) | — | 243; 0.09% of 266,898 | 🔴 2x at default (minor) |
| same at `quant_q ≤ 0.01` | 134; 0.08% of 162,219 | — | as above | 🟡 level |
| transferred classes | 202 / 567 / 0 / 0 | — | 27 / 52 / 0 / 164 | ⚪ |
| rate, enriched rows only | 0.052% | — | 0.010% | 🔴 (moderate) |
| rate, enriched plus no-reference rows | 0.052% | — | 0.072% | 🟢 |
| identified classes | 38 / 412 / 1 / 0 | — | 8 / 56 / 5 / 139 | ⚪ |
| identified rate, enriched only / with no-reference | 0.0053% / 0.0055% | — | 0.0014% / 0.0270% | 🔴 / 🟢 |
| blank run, rows identified | 0 on main of 2026-11-09, after #333 `e841aae0` (85 at `09e1c00d`) | 20 | 8 | 🟢 fewest |
| blank run, transfers attempted | 0 (run not aligned) | not reported | 25 candidates | ⚪ |
| blank run, rows transferred | 0 (65 at `09e1c00d`) | 85; 9 at `Match q-value ≤ 0.01` | 0 (scorer not trained) | ⚪ not tested |

- The reading flips with the no-reference rows, which the intensity profile cannot judge, so both
  are given. A tool's precursors that andes does not report cannot be classed.
- Blank run, identification: at `09e1c00d` andes took q-values from all six runs pooled, and 85
  rows of the blank passed 1%. With per-run q-values (#333) none does. Most of those 85 were good
  spectra of carry-over peptides, so the 0 means the blank has too few identifications for a
  per-run 1% estimate, not that 85 identifications were wrong.
- Blank run, transfers: neither zero is a tested rejection.
  - andes aligns runs on shared identifications. The blank has none, so it is not aligned and no
    transfer into it is attempted.
  - FragPipe (IonQuant) found 25 candidates in the blank but needs 10 identified ions in a run to
    train its transfer scorer. The blank has 8, so it scores nothing and writes no transfer.
  - MaxQuant aligns on MS1 features, so the blank is aligned and every match is reported; its MBR
    FDR option only fills the `Match q-value` column.
- The blank is not empty. MaxQuant's matched rows read a median 1.5% of the same peptide's
  intensity in the other runs (carry-over), so a row in the blank is not by itself a wrong value.
- A diagnostic andes build (not released) gives the blank the average alignment of the other
  runs. It fits 15,568 transfer targets in the blank and accepts 6 at `transfer_q ≤ 0.01`
  (7 proteins), against 0 of 23,639 decoys; the six read a median 7% of the same precursor's
  intensity in the other runs. The build is not used: adding the blank's rows to the pooled
  transfer model lowers the accepted transfers of the four single-cell runs by 51–66% (815 to 296
  in the first).

### Depth at matched true error (PXD014415)

| metric | andes | FragPipe | result |
|---|---:|---:|---|
| transferred rows at FragPipe's 0.09% (`transfer_q ≤ 0.003`) | 487,849 (0.88–0.89 × FragPipe's, over seeds) | 550,247 | 🔴 0.88–0.89x of FragPipe (minor) |
| transferred rows through `quant_q ≤ 0.015` | 428,888 | — | ⚪ |
| transferred rows at its own `quant_q ≤ 0.01` (0.08%) | 317,272 | — | ⚪ |

MaxQuant has no PXD014415 run. No tool reports a masked wrong-peak rate to compare.

## Large study: PXD042233 (500 runs)

Setup: label-free DDA from Q Exactive HF, HF-X and Exploris 480, `.raw` read natively, human
reviewed proteins + contaminants, 5 ppm, 2 missed cleavages, `--top-n 1`; one search job per run
(8 threads, `--quant-dir`), per-run PINs concatenated into one Percolator run (16 threads,
`--subset-max-train 1000000`), then `quant library`, `quant extract` per run, `quant finalize` (16
threads); andes `2051616b`, 2026-10-23, cluster nodes. No entrapment in this measurement.

| stage | wall | CPU | peak memory | output |
|---|---|---:|---:|---|
| search, 500 jobs | 0.7 min per run (at most 1.0) | 19.5 h | 3.6 GB | |
| PIN concatenation | 17 s | — | — | 53.1 million rows |
| Percolator | 12.9 min | 0.7 h | 26.0 GB | 30.6 million PSMs at `q ≤ 0.01` |
| `quant library` | 6.6 min | 0.1 h | 5.6 GB | |
| `quant extract`, 500 jobs | 2.1 min per run (two runs 45 and 72 min) | 121.5 h | 12.7 GB (median 4.7) | 121 million `below_threshold` rows at q 0.01 |
| `quant finalize`, transfer q 0.01, QPX views, protein inference | 25.8 min | 4.1 h | 45.2 GB | 194,370 transfers (175,087 `transfer_ms1`, 19,283 `transfer_ms2`); precursors in ≥ 1 / 3 / 50 runs 294,362 / 243,895 / 127,315 (identified only 294,362 / 243,730 / 126,741); 9,208 protein groups at protein q ≤ 0.01 |
| `quant finalize`, transfer q 0.05, no QPX views | 18.0 min | 3.9 h | 11.7 GB | 4,622,956 transfers; precursors in ≥ 1 / 3 / 50 runs 294,362 / 250,847 / 141,734 |

Extract is 83% of the CPU and most of its transfer fits are rejected. No other tool was run on this
set (⚪).
