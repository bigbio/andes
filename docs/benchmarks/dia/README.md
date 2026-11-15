# DIA (`--dia`) on the ProteoBench Astral entrapment set

Experimental, Phase A: spectrum-centric direct search of each DIA scan across its isolation window, with
co-elution, MS1 and neighbouring-cycle features rescored by Percolator (DOCS.md, "DIA"). This is the
dataset and the evaluation OpenDIA uses, so the numbers below compare directly with OpenDIA's and DIA-NN's.

## Data

- ProteoBench `LFQ_Astral_DIA_15min_50ng_Human` runs 01 and 02 (Orbitrap Astral, 300 isolation windows of 2 Th
  over 380–980 m/z, 1.65 s cycle, 292,883 / 292,914 MS2 scans, 2,928 MS1 scans), mzML from OpenMS FileConverter.
- Database: the PAIRLIB FASTA (`pairlib.fasta`, 100,634 one-peptide entries: human targets and their
  `_p_target` entrapment pairs). Its library holds 302,480 target precursors, half of them entrapment (r = 1).
- andes: `--dia --enzyme nocleavage --model hcd_astral_tryp --charge 1..5 --min-length 7 --max-length 35
  --predict rt,ms2`, Carbamidomethyl C fixed, Oxidation M variable (max 1), andes `XXX_` decoys
  ([`run.sh`](run.sh), [`mods.txt`](mods.txt)). `nocleavage` searches each FASTA entry as one peptide, which is
  the library's search space. Percolator 3.7.1 `--seed 42 -Y`, 14 threads.
- DIA-NN 2.2 with its own decoys on the PAIRLIB targets (`diann-pairlib-20261113`, `diann-decoys-20261114/own`).

## Evaluation ([`eval_dia.py`](eval_dia.py))

The same quantities as OpenDIA's `diann-decoys-20261114/evaluate.py`. A precursor is (peptidoform, charge),
scored by its best PSM. The q-value is precursor-level target-decoy competition `(D + 1) / T`. Only precursors
of the PAIRLIB library are counted. Entrapment FDP_c = E (1 + 1/r) / N. N@FDP_c 1 % is the number of IDs at a
true (entrapment) FDP of 1 %, a ranking measure that does not depend on how q-values are calibrated.

## Results (2026-11-15)

| tool | q ≤ 0.01 N, raw01 / raw02 | FDP_c at q ≤ 0.01 | FDP_c at q ≤ 0.05 | N@FDP_c 1 % | N@FDP_c 2 % |
|---|---|---|---|---|---|
| andes `--dia` | 1404 / 1499 | 0.004 / 0.005 | 0.042 / 0.046 | 1559 / 1651 | 1707 / 1691 |
| andes `--dia`, no MS1 gate | 1425 / 1465 | 0.008 / 0.004 | 0.032 / 0.040 | 1489 / 1569 | 1662 / 1731 |
| DIA-NN 2.2, own decoys | 2361 / 2322 | 0.014 / 0.010 | 0.060 / 0.060 | 2067 / 2301 | 2522 / 2543 |
| OpenDIA arm E (from its RESULT.md) | 1667 / 1457 | 0.010 / 0.007 | | 1694 / 1628 | 1908 / 1903 |

- **Error control holds**: andes's q ≤ 0.01 measures 0.4–0.5 % true FDP, and q ≤ 0.05 measures 4.2–4.6 %.
- **Ranking**: at a matched true FDP of 1 %, andes reaches 75 % (raw01) and 72 % (raw02) of DIA-NN, and
  92 % / 101 % of OpenDIA arm E.
- **MS1 gate** (a candidate is scored only when its monoisotopic m/z has a survey peak within 5 ppm in the
  nearest survey scan or one either side): +70 / +82 at N@FDP_c 1 % and 25–42 % less search time. On the
  slice it keeps the MS1 peak of 93 % of DIA-NN's precursors and of 13 % of random library candidates.
- **Cost per run**: andes search 6–7 min (wall, 14 threads) plus 12–15 min of Percolator; DIA-NN 8.8–21 min
  (`diann-decoys-20261114/RESULT.md`).

## Where the gap is

Measured before the MS1 gate. On a 3-minute slice (RT 10–13 min, [`slice.py`](slice.py), [`gap.py`](gap.py)), DIA-NN reports 604 library
precursors at q ≤ 0.01. In andes's PIN:

- **166 (27 %) are never retrieved.** They have no row in any scan. Most are 2+ near the low end of the
  m/z range (median 491 m/z). Raising the fragment-index cap from 100 to 300 and 1,000 candidates per scan
  did not retrieve them (in-PIN 438 → 422 → 417; N@FDP_c 1 % 275 → 238 → 294), so the cap is not the
  limit. The likely cause is spectrum-level matching: in a scan of a 2 Th window holding ~2,700 candidate
  forms, a weak precursor's few fragments do not lift it into the scan's top 10 rows.
- **About 140 more are retrieved but never ranked first in any scan**; Percolator rescues some of them.
- At the scan nearest DIA-NN's apex, the missed precursors still match a median of 5 b/y ions (13 of 166
  match fewer than 3), so they reach the fragment vote and lose the per-scan ranking to other forms of the
  window. The MS1 gate removes most of those competitors: on the slice, DIA-NN precursors in the PIN go
  from 438 to 473 and ranked first somewhere from 295 to 341.

Hypotheses tested on the slice that did not help: `--charge 2..3` or `2..4` instead of `1..5` (within
slice noise). The ranking score matters: `--score auto` picks the strong (fused) score for high-resolution
data, and on the slice it beats `--score rank` by a wide margin (N@FDP_c 1 % 285 vs 206; DIA-NN precursors
in the PIN 473 vs 309), so `--dia` keeps the auto choice.

Tested and rejected on the full raws or the slice:

- **MS1 gate width** (N@FDP_c 1 %, raw01 / raw02): 5 ppm ±1 survey scan (kept) 1559 / 1651; nearest
  survey scan only, 5 ppm, 1490 on raw01; 3 ppm ±1 scan 1618 / 1520, with q ≤ 0.01 at 1.0 / 1.35 % true
  FDP. The raw01 gain of 3 ppm does not hold on raw02 and its q-values are no longer calibrated.
- **Precursor-level rescoring** (slice): one row per (peptidoform, charge) with its best PSM's features plus
  its number of rows, scans where it ranks first, second-best and top-3 RawScore, and scan span, rescored
  by Percolator: N@FDP_c 1 % 285 → 192. Those aggregates barely separate targets from decoys, because a
  decoy has its target's m/z and passes the MS1 gate in the same scans. Counting spectrum matches is not
  peak-group evidence.

Slice probes (N@FDP_c 1 %): baseline 264; predicted MS2 features (`--predict rt,ms2`) 275; `--top-n 20` 232;
the neighbouring-cycle support columns had no measurable effect.

Both losses come from searching one scan at a time. Closing them is Phase B: a peptide-centric pass that
extracts each library precursor's fragment chromatograms across its elution window and scores the peak
group (DIA-NN, AlphaDIA), instead of asking each spectrum which peptide it holds.

## Peptide-centric pass (Phase B, `--dia-pc`, experimental, negative result)

`--dia-pc` (hidden, with `--dia`) adds a DIA-NN / OpenSWATH-style pass
(`crates/andes/src/bin/andes/dia_pc.rs`) and writes a precursor-level PIN `<stem>.dia.pin`:

- library: every target and decoy peptidoform of the database, at every charge whose m/z is in an
  isolation window (647,382 precursors);
- RT from andes-rt, calibrated on the Phase A identifications; mass calibration of fragments and MS1
  from the apexes of Phase A-identified precursors (−3.3 ppm fragments, +1.4 ppm MS1 on raw01);
- an assay of the 12 most intense andes-ms2 b/y ions; their traces over the precursor's window within
  the RT window; up to 8 local maxima scored, the 3 best by co-elution × spectral agreement kept as rows
  sharing a `ScanNr` so target-decoy competition keeps one per precursor;
- 23 peak-group scores: fragment co-elution (leave-one-out correlation), apex and integrated spectral
  angle against the prediction, matched fragments and their mass error, MS1 co-elution, error and isotope
  fit, RT error, intensities, and sparse-peak scores (fragments peaking at the apex, the top-3 predicted
  present, matched fragments above the window's background, lag of best alignment).

On the slice it does **not** beat the spectrum-centric search (N@FDP_c 1 %, Phase A 285):

| variant | Percolator | best offline rescorer |
|---|---|---|
| first version | 258 | 270 (with Phase A evidence joined) |
| + mass calibration, peak groups chosen by quality | 256 | 284 (gradient boosting) |
| ±3-scan leave-one-out correlation | 241 | 264 (hand-built score) |
| + sparse-peak scores | 244 | 250 (semi-supervised LDA) |
| narrower RT window (2.5 widths) | 131 | 226 |

Combining Phase A and Phase B per precursor (best z-scored score) gives 241. Phase B's IDs are mostly a
subset of Phase A's (372 and 210, 186 shared; Phase A ∪ B holds 272 of DIA-NN's 517 slice precursors,
Phase A alone 256).

Why, from the measurements:

- **Extraction is not the problem.** All 517 DIA-NN precursors of the slice are extracted. At DIA-NN's
  RT their 6 most intense library fragments match 4.7–5.8 peaks against 0.8–1.5 at a random RT, and
  83 target precursors (0 decoys) have peak groups correlating above 0.8.
- **Peaks are too short for chromatographic scores.** DIA-NN's peaks have a 3.3 s FWHM on a 1.65 s
  cycle, so a weak precursor's fragments appear in one to three scans, often not the same ones. The
  median fragment correlation of DIA-NN's true precursors is 0.45, 0.31 in the lowest-intensity
  quartile. That is the scan-level evidence Phase A already scores.
- **The predicted assay is weaker than an empirical one.** 2.46 of andes-ms2's top three fragments are
  present at the best apex, against 2.84 of the DIA-NN library's top three at DIA-NN's apex.
- **The RT prior is broad, for both tools.** The calibrated andes-rt error against DIA-NN's apexes is
  0.42 min robust σ (about 15 cycles), and DIA-NN's own library RT against its apexes is 0.38 min, so
  the generic predictors are about equally broad. `pc_rt_err` has no power here (AUC 0.46) and a
  narrower window loses true peaks. A run-specific refinement narrows it (see "Run-specific RT").

What would close the gap, from the same measurements: an empirical or fine-tuned fragment and RT
library (AlphaDIA-style transfer learning on the run's own identifications), a DIA-NN-style nonlinear
classifier trained per run, and interference handling between co-eluting precursors. Each is a larger
project than this pass. The code is kept behind the hidden flag, with its tests, as the starting point.

## Rescoring and RT tested on the full raws (Phase A PIN)

N@FDP_c 1 %, raw01 / raw02, against Percolator 1559 / 1651:

- **A nonlinear per-run classifier** (semi-supervised gradient boosting on the same 95 features, 3-fold
  cross-validated by precursor, Percolator `-Y` competition): 1553 / 1583 after one round, 1578 / 1561
  after two, 1520 / 1507 after three. It does not beat the linear SVM, so the gap is in the evidence,
  not in the classifier.
- **Run-specific RT.** A residual model of observed minus andes-rt RT, fitted on the run's q ≤ 1 %
  identifications (residue composition, length, oxidation, predicted RT; cross-validated by peptide),
  narrows the RT error from 0.43 to 0.18 min robust σ (95th percentile 1.06 to 0.54 min). As an
  extra PIN column it gives 1510 / 1648: Percolator weights it, but it barely separates DIA-NN's
  precursors that andes misses from decoys of the same score (AUC 0.46, against 0.49 for the
  unrefined RT error). Those precursors are already at the right RT; their spectra are what is weak.
- **Shared-scan interference** is not where the errors are. Of the targets at precursor q ≤ 5 % on
  raw01, 1,884 rank first in their scan (2.0 % entrapment) and 36 rank second (8 %), so a filter on
  rank would remove at most three false IDs.

Where DIA-NN's raw01 IDs sit in andes: 2,195 of 2,361 (93 %) have a row in the PIN, and for 84 % the
best row is inside DIA-NN's chromatographic peak. andes reaches q ≤ 1 % for 25 % of DIA-NN's
lowest-intensity quartile and 80 % of its highest. Against decoys of the same score, the missed ones
differ most in `matchedIonRatio` (AUC 0.75) and the andes-ms2 agreement (0.71), features the PIN
already has. Closing the remaining gap needs more fragment evidence per precursor than one scan
gives: summing or deconvolving the apex scans (a pseudo-spectrum per precursor), or an empirical
fragment library from the run.

## Reproduce

```bash
# mzML (OpenMS FileConverter) of the two raws, then from this directory:
./run.sh full01 raw01.mzML --enzyme nocleavage --model hcd_astral_tryp --predict rt,ms2
./run.sh full02 raw02.mzML --enzyme nocleavage --model hcd_astral_tryp --predict rt,ms2
# slice for iteration
python3 slice.py raw01.mzML slice01.mzML 600 780   # then set spectrumList count to the printed number
```

`run.sh` and `eval_dia.py` read the OpenDIA benchmark root
`/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107` (PAIRLIB FASTA and `lib_targets.tsv`); edit `B=` to
move it.
