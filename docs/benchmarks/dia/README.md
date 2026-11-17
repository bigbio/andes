# DIA (`--dia`) on the ProteoBench Astral entrapment set

Spectrum-centric direct search of each DIA scan across its isolation window ("Phase A" below), with
co-elution, MS1 and neighbouring-cycle features rescored by Percolator (DOCS.md, "DIA"). The rows of
the results table are the steps that were kept, in order; the sections after it record every variant
that was tested and dropped, with its numbers. This is the
dataset and the evaluation OpenDIA uses, so the numbers below compare directly with OpenDIA's and DIA-NN's.

## Data

- ProteoBench `LFQ_Astral_DIA_15min_50ng_Human` runs 01 and 02 (Orbitrap Astral, 300 isolation windows of 2 Th
  over 380–980 m/z, 1.65 s cycle, 292,883 / 292,914 MS2 scans, 2,928 MS1 scans), mzML from OpenMS FileConverter.
- Database: the PAIRLIB FASTA (`pairlib.fasta`, 100,634 one-peptide entries: human targets and their
  `_p_target` entrapment pairs). Its library holds 302,480 target precursors, half of them entrapment (r = 1).
- andes: `--dia --enzyme nocleavage --model hcd_astral_tryp --charge 1..5 --min-length 7 --max-length 35
  --predict rt,ms2`, Carbamidomethyl C fixed, Oxidation M variable (max 1), andes `XXX_` decoys
  ([`run.sh`](run.sh)). `nocleavage` searches each FASTA entry as one peptide, which is
  the library's search space. Percolator 3.7.1 `--seed 42 -Y`, 14 threads.
- DIA-NN 2.2 with its own decoys on the PAIRLIB targets (`diann-pairlib-20261113`, `diann-decoys-20261114/own`).

## Evaluation ([`eval_dia.py`](eval_dia.py))

The same quantities as OpenDIA's `diann-decoys-20261114/evaluate.py`. A precursor is (peptidoform, charge),
scored by its best PSM. The q-value is precursor-level target-decoy competition `(D + 1) / T`. Only precursors
of the PAIRLIB library are counted. Entrapment FDP_c = E (1 + 1/r) / N. N@FDP_c 1 % is the number of IDs at a
true (entrapment) FDP of 1 %, a ranking measure that does not depend on how q-values are calibrated.

## Results (2026-11-15 to 2026-11-16)

| tool | q ≤ 0.01 N, raw01 / raw02 | FDP_c at q ≤ 0.01 | FDP_c at q ≤ 0.05 | N@FDP_c 1 % | N@FDP_c 2 % |
|---|---|---|---|---|---|
| andes `--dia` | 1404 / 1499 | 0.004 / 0.005 | 0.042 / 0.046 | 1559 / 1651 | 1707 / 1691 |
| andes `--dia`, mean of 5 Percolator seeds | 1394 / 1565 | 0.007 / 0.009 | | 1496 / 1594 (1421–1575 / 1524–1651) | 1685 / 1727 |
| andes `--dia`, `--trainFDR 0.05` (default with `--dia --rescore`), mean of 5 seeds | 1493 / 1602 | 0.005 / 0.011 | 0.049 / 0.044 | 1635 / 1611 (1568–1672 / 1557–1647) | 1793 / 1744 |
| andes `--dia` + summed-apex columns, `--trainFDR 0.05`, mean of 5 seeds | 1433 / 1589 | 0.004 / 0.007 | | 1634 / 1687 (1619–1646 / 1620–1725) | 1798 / 1819 |
| andes `--dia --rescore` + two-pass RT refinement, 20 ppm fragments, mean of 5 seeds | 1432 / 1612 | 0.004 / 0.005 | 0.044 / 0.038 | 1666 / 1747 (1627–1688 / 1648–1792) | 1872 / 1884 |
| andes `--dia --rescore` (+ 7 ppm fragments), mean of 5 seeds | 1589 / 1755 | 0.004 / 0.009 | 0.043 / 0.050 | 1801 / 1819 (1733–1828 / 1748–1882) | 1964 / 1962 |
| andes `--dia --rescore`, current (main's UNIMOD masses), mean of 5 seeds | 1606 / 1747 | 0.003 / 0.007 | 0.044 / 0.045 | **1800 / 1821** (1740–1838 / 1791–1853) | 1909 / 1953 |
| andes `--dia`, no MS1 gate | 1425 / 1465 | 0.008 / 0.004 | 0.032 / 0.040 | 1489 / 1569 | 1662 / 1731 |
| DIA-NN 2.2, own decoys | 2361 / 2322 | 0.014 / 0.010 | 0.060 / 0.060 | 2067 / 2301 | 2522 / 2543 |
| OpenDIA arm E (from its RESULT.md) | 1667 / 1457 | 0.010 / 0.007 | | 1694 / 1628 | 1908 / 1903 |

- **Error control**: at 20 ppm andes's q ≤ 0.01 measured 0.4–0.5 % true FDP (0.7–0.9 % averaged over
  seeds) and q ≤ 0.05 4.2–4.6 %. At 7 ppm with the old mod masses, q ≤ 0.01 averaged 0.4 % on raw01
  and 0.9 % on raw02 (0.6–1.3 %, above 1 % for seed 42). The current build measures 0.3 % / 0.7 %
  (0.1–0.5 % / 0.6–1.0 % per seed) and q ≤ 0.05 4.4 / 4.5 %.
- **Percolator seed noise.** The andes PIN is deterministic (byte-identical across runs), but the
  Percolator seed alone moves N@FDP_c 1 % over a range of about 150 (seeds 42, 1, 2, 3, 7). Seed 42,
  used for the first row and the single-seed comparisons below, is at the top of that range on both
  raws. Compare features on the mean over seeds, paired by seed and by thread count: the same seed
  with `--num-threads` 7 or 10 gave 1483 or 1559 on one raw01 PIN. Averaging the PSM scores of the
  five seeds before precursor competition gives 1482 / 1682, no steadier than one seed.
- **Ranking**: at a matched true FDP of 1 %, andes reaches 75 % (raw01) and 72 % (raw02) of DIA-NN
  with seed 42, 72 % / 69 % on the seed mean, and 87 % / 79 % with `--trainFDR 0.05`, the
  summed-apex columns, the two-pass RT refinement and 7 ppm fragments; it is above OpenDIA arm E
  (106 % / 112 %).
- **Modification masses.** main now reads modification masses from UNIMOD (Oxidation 15.994915,
  Carbamidomethyl 57.021464, against 15.99491 / 57.02146 in the `mods.txt` the earlier rows used).
  The PIN then writes oxidation as `+15.99492`, which `eval_dia.py` now maps to UniMod:35 as well
  (before the fix, those precursors fell out as "not in library" and the current row read 1777 /
  1797). On the slice, the old build with the UNIMOD masses writes a PIN byte-identical to the
  current build, so the merge and cleanup do not change the search.
- **Fragment tolerance.** andes measures each run's fragment mass error from 100 confident PSMs of a
  3,000-spectrum sample. A DIA scan's best candidate rarely reaches 1 % q (21 of 3,000 on raw01), so
  the measurement never fired and the search kept the 20 ppm default. Mean of 5 seeds, N@FDP_c 1 %
  with `--fragment-tol-ppm` 5 / 7 / 10 / 15 / 20 (retrieval and features): 1567 / 1801 / 1720 / 1695 /
  1666 on raw01, 1685 / 1819 / 1833 / 1731 / 1747 on raw02. 7 ppm is higher than 20 ppm on raw01 for
  5 of 5 seeds and on raw02 for 4 of 5. Measuring on a 30,000-spectrum sample instead (offset −3.3 /
  −2.6 ppm, spread 1.3 ppm, 5 ppm features, 20 ppm retrieval) gave 1747 / 1681, below a flat 7 ppm.
  `--dia` on Orbitrap Astral MS2 now uses 7 ppm unless `--fragment-tol-*` is set (other analyzers keep
  20 ppm: only Astral was measured). Matching the
  co-elution and summed-apex columns at 7 or 12 ppm too (they use 20 and 10 ppm) gave 1817 / 1755 and
  1821 / 1785 with q ≤ 0.01 at 1.0–1.1 % true FDP on raw02; not kept. At 7 ppm, `--trainFDR` 0.02 /
  0.05 / 0.10 give 1726 / 1801 / 1822 on raw01 and 1703 / 1819 / 1824 on raw02 (0.05 kept: 0.10 is
  inside seed noise), and andes-ms2 at NCE 16 gives 1767 / 1811 (not kept). `--charge 2..4` (DIA-NN's
  IDs are 75 % 2+ and 22 % 3+) gives 1751 / 1826, `--top-n 5` 1722 / 1689 and `--top-n 15` 1711 / 1683
  (q ≤ 0.01 at 1.4 % true FDP on raw02); none kept, the default 10 rows per scan is best.
- **Percolator training FDR.** A DIA PIN has about ten candidate rows per scan and only ~6.7k PSMs at
  1 %, so Percolator's default `--trainFDR 0.01` trains on few positives. Mean N@FDP_c 1 % over seeds
  1, 2, 3 (10 threads, paired): `--trainFDR` 0.01 / 0.02 / 0.05 / 0.10 give 1449 / 1552 / 1615 / 1633
  on raw01 and 1568 / 1579 / 1609 / 1596 on raw02. Over all 5 seeds 0.05 is higher than the baseline on
  raw01 for 5 of 5 seeds and on raw02 for 3 of 5. Error control is unchanged (q ≤ 0.01 at 0.3–1.3 %
  true FDP). `--dia --rescore` now passes `--trainFDR 0.05` unless `--percolator-args` sets it.
  `--train-best-positive` is worse (3 seeds: 1403 / 1512).
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

Both losses come from searching one scan at a time. The natural remedy is a peptide-centric pass that
extracts each library precursor's fragment chromatograms across its elution window and scores the peak
group (DIA-NN, AlphaDIA), instead of asking each spectrum which peptide it holds. That was tried next
and did not close them on this data (below).

## Peptide-centric pass (Phase B, negative result, not in andes)

A DIA-NN / OpenSWATH-style pass was built and measured, then removed; the code is on the
`archive/dia-pc` branch. It wrote a precursor-level PIN from:

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
project than this pass. The archived code, with its tests, is the starting point.

## Rescoring and RT tested on the full raws (Phase A PIN)

N@FDP_c 1 %, raw01 / raw02, against Percolator 1559 / 1651 (seed 42; the first three items were run on
that seed only, and differences below about 100 are inside seed noise, see Results):

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
  Re-run under `--trainFDR 0.05` on the summed-apex PINs (5 seeds, paired, Percolator on codon): N@FDP_c
  1 % 1700 / 1722 against 1634 / 1687 (raw01 higher on 5 of 5 seeds, raw02 on 4 of 5), N@FDP_c 2 %
  1899 / 1853 against 1798 / 1819, q ≤ 0.01 at 0.25 / 0.42 % true FDP. It helps once Percolator
  trains on enough positives; the seed-42, `--trainFDR 0.01` rejection above was noise. This offline
  model (gradient boosting, anchors from a previous Percolator run) is not yet in andes. Variants
  that andes could compute in one pass, same 5 seeds (N@FDP_c 1 %, raw01 / raw02, against 1634 / 1687):
  ridge on Percolator anchors 1634 / 1753, ridge on anchors from target-decoy competition on
  `RawScore` at q ≤ 1 % 1603 / 1745, boosting on those anchors 1635 / 1722, at q ≤ 5 % 1648 / 1717.
  raw02 gains 30–66 with every variant; raw01 gains only with boosting on Percolator anchors, though
  its RT error is no smaller than ridge's (median 0.24 min for both, against 0.30 unrefined).
  **Shipped as a second Percolator pass** (`dia::refine_rt`, andes's own gradient boosting, depth 3,
  200 rounds, out of fold by peptidoform, anchors = first-pass targets at q ≤ 1 %, column
  `dia_rt_refined_err`): run end to end with `andes --dia --rescore` on codon, then Percolator
  re-scored on its PIN for 5 seeds: 1666 / 1747 against 1634 / 1687 (raw01 higher on 5 of 5 seeds,
  raw02 on 4 of 5), N@FDP_c 2 % 1872 / 1884, q ≤ 0.01 at 0.44 / 0.47 % true FDP. It adds one
  Percolator run (3–6 min) and ~8 s of fitting per raw. A boosting refit with a ridge on
  hinge-basis RT (spline) gave 1649 / 1745, a little below.
- **Shared-scan interference** is not where the errors are. Of the targets at precursor q ≤ 5 % on
  raw01, 1,884 rank first in their scan (2.0 % entrapment) and 36 rank second (8 %), so a filter on
  rank would remove at most three false IDs.
- **Summed apex.** Three columns per row: of the peptide's b/y ions, how many occur in at least two
  of its scan and the previous and next cycle of the same window, that count as a fraction, and how
  many the neighbours add. On the 3-min slice this lifted N@FDP_c 1 % from 285 to 340 (Python, 10 ppm)
  and 309 / 362 / 268 (in `dia.rs` at 20 / 10 / 5 ppm). On the full raws at 20 ppm, seed 42: 1539 / 1564
  against 1559 / 1651. At 10 ppm, mean over 5 seeds paired with the baseline: N@FDP_c 1 % 1452 / 1681
  against 1496 / 1594 (raw01 lower on 5 of 5 seeds, raw02 higher on 5 of 5), N@FDP_c 2 % 1703 / 1788
  against 1685 / 1727. Mixed under `--trainFDR 0.01`. Under `--trainFDR 0.05` (5 seeds, paired):
  1634 / 1687 against 1635 / 1611 (raw01 neutral, raw02 +76, higher on 4 of 5 seeds), N@FDP_c 2 %
  1798 / 1819 against 1793 / 1744, with the raw01 seed spread narrowing from 104 to 27. Kept
  (`dia_sum_n2`, `dia_sum_frac2`, `dia_sum_gain`, 10 ppm). The slice is too small to rank features
  (±50 IDs between near-identical runs); judge on full raws and several Percolator seeds.
- **andes-ms2 collision energy** (run-specific fragment intensity by conditioning). The run's
  metadata gives NCE 25 (Exploris family). `andes --dia --rescore` output (seed 42) at NCE
  19 / 22 / 25 / 28 / 31: 1652 / 1803, 1633 / 1759, 1622 / 1792, 1543 / 1774, 1580 / 1711; the
  Astral instrument class gave 1628 / 1415. Over 5 seeds paired with NCE 25 (1666 / 1747):
  NCE 16 1694 / 1757 (raw01 higher on 5 of 5, raw02 on 2 of 5), NCE 19 1679 / 1760 (3 / 1 of 5),
  NCE 22 1665 / 1784 (3 / 3 of 5). Agreement with anchor spectra rises as NCE falls (median spectral
  angle 0.721 / 0.750 / 0.770 at NCE 28 / 25 / 22), but the ID gain is at most +28 and inside seed
  noise on raw02. Not shipped; a per-run NCE search on anchor agreement is not worth its cost.
- **Both raws in one run** (`--spectrum raw01 --spectrum raw02`, one pooled Percolator model,
  per-run q-values), 5 seeds, scored per raw: 1632 / 1745 against 1666 / 1747. Pooling does not help.

Where DIA-NN's raw01 IDs sit in andes: 2,195 of 2,361 (93 %) have a row in the PIN, and for 84 % the
best row is inside DIA-NN's chromatographic peak. andes reaches q ≤ 1 % for 25 % of DIA-NN's
lowest-intensity quartile and 80 % of its highest. Against decoys of the same score, the missed ones
differ most in `matchedIonRatio` (AUC 0.75) and the andes-ms2 agreement (0.71), features the PIN
already has. Counting repeats over adjacent scans adds a little on raw02 (above); what is left is
deconvolving the apex scans or an empirical fragment library from the run.

## Reproduce

```bash
# mzML (OpenMS FileConverter) of the two raws; B is the OpenDIA benchmark root (PAIRLIB FASTA and
# lib_targets.tsv), ANDES the binary, PY a python with pandas.
export B=/path/to/OpenMS-r1-benchmark ANDES=andes PY=python3
./run.sh full01 raw01.mzML --enzyme nocleavage --model hcd_astral_tryp --predict rt,ms2
./run.sh full02 raw02.mzML --enzyme nocleavage --model hcd_astral_tryp --predict rt,ms2
# slice for iteration
python3 slice.py raw01.mzML slice01.mzML 600 780   # then set spectrumList count to the printed number
```

`run.sh` runs Percolator once at `--trainFDR 0.05`, as the first pass of `--dia --rescore`. The
two-pass numbers in the table (RT refinement) are `andes --dia --rescore` itself: its
`<stem>.pin` and Percolator output feed `eval_dia.py` the same way.
