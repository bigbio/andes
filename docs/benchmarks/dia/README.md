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
