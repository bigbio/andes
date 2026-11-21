# PB-Astral-full: `--dia` on the full ProteoBench Astral entrapment FASTA

This benchmark runs the same runs and evaluation as [../README.md](../README.md), with one difference: the
database is ProteoBench's full paired entrapment FASTA, not the PAIRLIB subset. It is used to measure
per-scan scores for `--dia` (`dia_llr.rs`), which are now default features of `--dia --predict ms2`. The scripts take every path from variables and arguments,
so the same comparison runs on PAIRLIB too (see "PAIRLIB" below).

## Data

- ProteoBench `LFQ_Astral_DIA_15min_50ng_Human` runs 01 and 02 for A/B tests, and run 03 to confirm a
  kept variant. mzML from ThermoRawFileParser 1.4.5 (centroided).
- Database: `ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta`. It holds 2,841,736
  one-peptide entries, of which 1,420,868 are `_p_target` entrapment peptides paired with the human
  targets (r = 1).
- DIA-NN 2.3.2, run on the same FASTA with its own decoys: `--fasta-search --predictor`, peptide length
  7–30, precursor charge 1–4, precursor m/z 300–1800, fragment m/z 200–1800, `--cut ""` (each entry
  searched as one peptide), `--unimod4` (fixed carbamidomethyl C), no variable modifications,
  `--qvalue 0.01`, `--reanalyse --rt-profiling`, all three runs in one job. The per-run comparison uses
  `report-first-pass.parquet` (before match-between-runs). `report.parquet` is also listed.
- andes search space, matched to DIA-NN ([`run.sh`](run.sh)): `--enzyme nocleavage --fixed-mods
  "Carbamidomethyl (C)" --variable-mods "" --charge 1..4 --min-length 7 --max-length 30 --dia --model
  hcd_astral_tryp --predict rt,ms2 --rescore`, andes `XXX_` decoys, 14 threads. Fragment tolerance is the
  Astral `--dia` default of 7 ppm. andes has no fragment m/z floor, so it keeps fragments below 200 m/z.
  The isolation windows (380–980 m/z) cover less than DIA-NN's precursor m/z range.
- Percolator 3.7.1, `-Y --trainFDR 0.05 --num-threads 14`, re-run on andes's final PIN (after the RT
  refinement pass of `--rescore`) with seeds 42, 1, 2, 3 and 7.

## Evaluation ([`eval.py`](eval.py))

These are the quantities of `../eval_dia.py`. A precursor is (peptidoform, charge), scored by its best PSM.
q = (D + 1) / T is computed from precursor-level target-decoy competition. Entrapment means every
protein of the PSM carries `_p_target`. FDP_c = E (1 + 1/r) / N = 2E / N. N@FDP_c 1 % / 2 % counts targets
ranked by score until the combined FDP reaches 1 % / 2 %. Each number is the mean over the 5 seeds,
paired by seed, with the per-seed range. DIA-NN's report is filtered at q ≤ 0.01, so its N@FDP_c can
only be the count at q ≤ 0.01.

## PAIRLIB

`eval.py andes <pin> <t.psms> <d.psms> --library $B/diann-pairlib-20261113/lib_targets.tsv` reproduces
`../eval_dia.py`: it keeps library precursors only, takes entrapment from the library, and computes r
from the library. To run on PAIRLIB, use `../run.sh` with `--predict rt,ms2`; the pool scores need no
flag. PAIRLIB inputs were not available for the D1 measurement below.

## Results (2026-11-18)

The arms in the first table were measured with hidden experiment switches (`--dia-score-cols`,
`--dia-rerank-key`, `--dia-pool`, and others) that no longer exist. The default build is the "all five
pool scores as features" arm with a pool of 100 rows; the table after them gives its numbers. The
experiment build is kept on the branch `exp/pr417-combo`.

Columns: N@FDP_c 1 % as the mean over 5 seeds (per-seed range); the paired difference against the
baseline, with the number of seeds where it is higher; FDP_c at q ≤ 0.01 (range over seeds); FDP_c at
q ≤ 0.05 (mean). raw01 / raw02.

| andes | N@FDP_c 1 %, raw01 | raw02 | vs baseline (seeds higher) | FDP_c at q ≤ 0.01 | FDP_c at q ≤ 0.05 |
|---|---|---|---|---|---|
| baseline: main with #410 (pool re-ranked by RawScore + angle) | 54546 (54433–54736) | 55420 (55310–55463) | | 0.94–0.97 % / 0.97–1.01 % | 4.98 / 4.69 % |
| pool scores as features (`--dia-score-cols dia_llr,dia_apex_llr`) | 55982 | 57163 | +1436 (5/5) / +1744 (5/5) | 0.96–0.98 / 0.91–0.95 % | 4.96 / 4.77 % |
| `dia_llr` only as a feature | 55809 | 56930 | +1262 (5/5) / +1510 (5/5) | 0.93–0.96 / 0.91–0.98 % | 4.96 / 4.85 % |
| re-rank by `dia_llr` | 54305 | 55702 | −241 (0/5) / +282 (5/5) | 1.00–1.04 / 0.90–0.93 % | 5.05 / 4.81 % |
| re-rank by RawScore + angle + `dia_llr` | 54394 | 55211 | −152 (0/5) / −209 (0/5) | 0.94–0.98 / 0.97–1.02 % | 4.99 / 4.73 % |
| re-rank by `dia_llr` + `dia_apex_llr` | 54240 | 54985 | −306 (0/5) / −434 (0/5) | 0.92–1.00 / 0.96–1.03 % | 4.95 / 4.94 % |
| re-rank by RawScore + angle + `dia_llr` + `dia_apex_llr` | 54613 | 55469 | +67 (4/5) / +49 (3/5) | 0.94–0.96 / 1.02–1.05 % | 5.00 / 4.85 % |
| features `dia_llr`, `dia_apex_llr` + re-rank by `dia_llr` | 55576 | 56627 | +1029 (5/5) / +1207 (5/5) | 0.98–1.01 / 0.96–1.02 % | 4.96 / 4.80 % |
| re-rank by `dia_llr_greedy` | 55703 | 56627 | +1157 (5/5) / +1207 (5/5) | 0.94–1.01 / 0.91–0.99 % | 5.16 / 4.79 % |
| re-rank by `dia_llr_soft` | 54405 | 55889 | −141 (0/5) / +469 (5/5) | 0.98–1.00 / 0.88–0.91 % | 5.14 / 4.83 % |
| re-rank by RawScore + `dia_angle_soft` | 54366 | 55355 | −180 (0/5) / −64 (1/5) | 0.94–0.98 / 0.96–1.06 % | 4.97 / 4.66 % |
| **all five pool scores as features** (`--dia-score-cols dia_llr,dia_apex_llr,dia_llr_greedy,dia_llr_soft,dia_angle_soft`) | **58533** (58477–58602) | **59267** (59169–59395) | **+3987 (5/5) / +3848 (5/5)** | 0.94–0.97 / 0.94–0.99 % | 4.86 / 4.85 % |
| the same without `dia_apex_llr` | 58449 | 59277 | +3902 (5/5) / +3857 (5/5) | 0.91–0.92 / 0.94–0.99 % | 4.86 / 4.83 % |
| all five as features + re-rank by `dia_llr_greedy` | 58434 | 59281 | +3888 (5/5) / +3861 (5/5) | 0.97–1.00 / 0.89–0.95 % | 5.08 / 4.93 % |
| DIA-NN 2.3.2, first pass (per run) | 82068 at q ≤ 0.01 | 83424 | | 0.84 / 0.83 % | |
| DIA-NN 2.3.2, final (match-between-runs over 3 runs) | 87611 at q ≤ 0.01 | 88047 | | 0.54 / 0.53 % | |

Confirmation on raw03 (not used for any choice): baseline 54350 (54301–54426), q ≤ 0.01 at
1.05–1.09 % true FDP. All five pool scores as features: 58440 (58346–58557), +4090, higher on 5 of 5
seeds, q ≤ 0.01 at 0.99–1.02 %. Re-ranking by `dia_llr_greedy` on top of that: 58382, +4032 (5/5). DIA-NN
first pass on raw03: 81598 at 0.79 %.

- **Verdict (experiment arms).** The five pool scores as Percolator features are well ahead of the baseline on all three
  raws (+7.3 % / +6.9 % / +7.5 %). Error control is unchanged: q ≤ 0.01 stays at or below about 1 % true FDP,
  and q ≤ 0.05 at about 4.9 %. The gain is more than ten times the per-seed range of any arm. Dropping `dia_apex_llr`, or also
  re-ranking by `dia_llr_greedy`, gives the same result within seed noise. andes reaches 71 % of DIA-NN's
  per-run (first-pass) count at a matched true FDP, against 66 % before.
- **Re-ranking keys alone** move N@FDP_c by −434 to +469, except the greedy key (+1157 / +1207). A key
  helps only when it changes which rows reach the PIN; once the scores are Percolator features, the
  re-ranking key no longer matters.
- **Cost.** The pool scores take 40–60 s per run on 14 threads and add 0.4 GB to the peak (19.4 → 19.8
  GB). Search plus `--rescore` takes 17–23 min per raw in every arm, and the variation follows the node
  rather than the arm. Percolator per seed takes 6–9 min in both the baseline and the kept arm. `--dia-dump`
  adds 4 GB and writes about 330 MB of Parquet per table and run.
- **Dump and trace tooling.** The per-row and per-ion dumps (`--dia-dump`) and the retrieval trace
  (`--dia-trace`) were removed from the default code; they live on `exp/pr417-combo`.

### Default build: pool scores and a pool of 100 rows (2026-11-18)

Percolator seeds 42, 1, 2, 3 and 7, Percolator 3.7.1, 14 threads; main is commit 232dd79a. Search plus
`--rescore` wall time and peak memory are for seed-independent single runs; Percolator is the seed 42
run.

| build | N@FDP_c 1 %, raw01 | raw02 | raw03 | vs main (seeds higher), raw01 / raw02 / raw03 | FDP_c at q <= 0.01 (raw01 / raw02 / raw03) | search + rescore, raw01 | peak memory | Percolator, seed 42 |
|---|---|---|---|---|---|---|---|---|
| main | 57413 | 58167 | 57426 | | 0.92-0.97 / 1.04-1.07 / 1.04-1.07 % | 19.7 min | 20.1 GB | 7.3 min, 11.5 GB |
| pool scores + pool 100 | 65055 | 66366 | 65237 | +7642 / +8199 / +7811 (5/5 each) | 0.90-0.95 / 0.89-0.95 / 1.03-1.07 % | 23.8 min | 29.6 GB | 8.9 min, 12.1 GB |

DIA-NN 2.3.2 first pass is 82068 / 83424 / 81598 at q <= 0.01, so andes is at 79 % / 80 % / 80 % of it at a
matched true FDP of 1 %. On raw03 main's nominal q <= 0.01 is already at 1.04-1.07 % true FDP, and the
change does not move it. PAIRLIB was not run (its inputs are not on the cluster).

## Full stack: pass 2, learned score, study mode and transfer

Data as above: the full paired entrapment FASTA, the three runs `LFQ_Astral_DIA_15min_50ng_Human`
01, 02 and 03, treated as one study. andes: per run `--dia --predict rt,ms2 --rescore --quant-dir`,
then `andes quant extract --run <stem>` per run, then `andes quant finalize`, with finalize seeds
42, 1, 2, 3 and 4. The q reported by andes is class-wise (own and transferred precursors each have
their own competition). DIA-NN 2.3.2 on the same FASTA, three runs in one job. The columns are
raw01 / raw02 / raw03.

| measure | raw01 | raw02 | raw03 |
|---|---|---|---|
| N@FDP_c 1 %, mean ± sd over 5 finalize seeds | 84,271 ± 202 | 85,268 ± 102 | 84,609 ± 113 |
| N@FDP_c 1 %, pooled q instead of class-wise | 84,542 | 85,652 | 85,007 |
| N at reported q ≤ 0.01 | 84,804 | 85,591 | 84,935 |
| FDP_c at q ≤ 0.01, all rows | 1.06 % | 1.04 % | 1.03 % |
| FDP_c at q ≤ 0.01, own nominees | 1.06 % | 1.05 % | 1.04 % |
| FDP_c at q ≤ 0.01, transferred (N about 3,000) | 1.28 ± 0.14 % | 0.90 ± 0.13 % | 0.94 ± 0.06 % |
| FDP_c at q ≤ 0.005 / 0.02, all rows | 0.55 / 2.03 % | 0.52 / 2.08 % | 0.52 / 2.05 % |
| DIA-NN 2.3.2, own q ≤ 0.01 (FDP_c 0.53 %) | 85,127 | 85,238 | 85,158 |
| DIA-NN 2.3.2, N@FDP_c 1 % | at least 88.8k | | |

Protein groups at q 0.01 over the three runs: 82,121 ± 148. DIA-NN's N@FDP_c 1 % is a lower bound:
its report is capped at q 0.05.

Verdicts:

- **Against the build before this work** (about 57-58k at N@FDP_c 1 %, measured on the base of the
  stack; the stack's own base commit was not re-measured per raw): well ahead, about +27k per run.
  The layers add up as follows, each measured on its own and not additive in general: pool scores
  and recentring, learned score (+7 % / +8 %), joint columns (+1.6 % / +1.8 %), RT-aware cut
  (+3.1 % to +3.6 %), cross-run columns (about +5.0k / +4.2k), transfer (+5.5k to +6.3k).
- **Against DIA-NN at matched true FDP of 1 %:** behind, by at least 4 to 5 % (andes is at most
  94.9-96.0 % of DIA-NN's lower bound).
- **Against DIA-NN's own q count:** on par, about 99-100 %.
- **Found of the DIA-NN precursors that andes never nominated before this work** (8,934 on
  raw01): 4,865.
- **Error control, all rows:** slightly above nominal at q 0.01 (1.03-1.06 %). It comes from decoy
  versus entrapment non-exchangeability that is present before Percolator; N@FDP_c 1 % accounts
  for it.
- **Error control, transferred class:** above nominal on raw01 (1.28 %), at or below it on raw02
  and raw03. Under one pooled competition the transferred targets measured 2.2-2.7 % at q 0.01,
  which is why the q is class-wise. The cost of class-wise q is about 300-580 precursors per run.
- **Leakage:** entrapment-versus-decoy AUC of the learned score, the Percolator score among
  transferred rows, `xfer` and the cross-run columns is 0.499-0.501 (0.4993-0.5008 on every run).

Cost per Astral run (15 min gradient): search 19.7 min on one node (37.7-37.9 min on a slower node,
same code), peak 35.3-35.8 GB; extract with transfer 8.7 min, 13.8 GB; finalize over three runs
64.7 ± 1.2 min, 22.6 GB. The RT-aware cut adds about 20 % to the search time. The learned score
takes about 570 s per run.

Limits:

- The three runs are injections of one sample. Cross-run columns and transfer across different
  samples, with RT drift, are not measured. This is the best case for them.
- diaPASEF is untested with pass 2 and transfer.
- The raw01 transferred class is above nominal at every finalize seed. About 3,000 transferred
  precursors pass per run, so E is about 20 and the estimate is noisy. A minimum class size or a
  +1 decoy correction for the class has not been tested.
- The per-run search peaks at about 36 GB.
- The RT-aware cut falls back to votes when the RT prediction gives fewer than 50 anchors.
- The layers other than the full stack were each measured with one seed or one raw in places
  (joint columns, RT-aware cut); the full-stack numbers above use 5 finalize seeds.
- Behind DIA-NN at matched true FDP, as stated.

User description: [DOCS.md](../../../../DOCS.md), "DIA pass 2, learned score and study mode".

## Per-ion measurements (`--dia-dump`, raw01; raw02 agrees to ±0.01)

Anchors are Percolator q ≤ 1 % targets. Decoy rows are all decoy rows of the 50-row pool. Ions are
sampled from every 10th cycle (31 M ion records).

- **Chance presence** p0 (decoy rows), z = 1, ions of 4+ residues by m/z (200 Th bins from 200):
  0.32, 0.18, 0.14, 0.10, 0.07, 0.05, 0.03, 0.015, 0.005. Short ions are mostly present by chance (length 1:
  0.87, length 2: 0.62–0.75, length 3 below 400 m/z: 0.37–0.65). For z = 2, p0 is 0.000–0.055. By local
  peak density (peaks within ±5 Th), z = 1: 0.17 at ≤ 10 peaks, 0.37–0.39 above 10.
- **True presence** q (anchors) by andes-ms2 intensity rank, MS1 quartile Q1 → Q4: rank 1 0.89 → 0.99,
  rank 3 0.72 → 0.98, rank 7–8 0.44 → 0.90, rank 11–15 0.21 → 0.59, rank 21+ 0.03 → 0.16.
  DIA-NN precursors that andes misses, at their best row at DIA-NN's RT: rank 1 0.87, rank 4–6 0.68,
  rank 7–10 0.46, rank 11–20 0.15, against 0.95 / 0.77 / 0.60 / 0.28 for anchors.
- **Same apex** m/k (of the fragments present in s−1, s or s+1, the share peaking at s): anchors 0.35,
  decoy rows 0.53, decoys that win their scan 0.66, DIA-NN precursors found / missed 0.55 / 0.52. A
  random match is usually present in one scan only, so it peaks there by construction. The term
  therefore does not separate (AUC 0.44, DIA-NN-missed best rows against decoys that win their scan).
  Continuity (fragments present at s that are also present at s ± 1) does separate: anchors 0.81, decoys
  0.56, DIA-NN found 0.76, missed 0.66. That evidence is already in `dia_sum_*`.
- **Where DIA-NN's raw01 precursors are** (81,725 at 1 %, non-entrapment, first pass): andes finds 51,029.
  8,144 are in the pool at DIA-NN's RT (±0.1 min) but never rank first by the #410 key, 1,857 rank first
  there but are not accepted, 7,300 have pool rows only away from DIA-NN's RT, and 13,395 are not in the
  pool. For the 8,144 that never rank first, the winning row at their best scan is a decoy 6.6 % of the
  time (PAIRLIB: 45 %), 5.1 % under the greedy key and 8.2 % under the soft key. With the full FASTA
  the competitors are mostly other targets.
- **Separation** of DIA-NN-missed best rows from decoys that win their scan (AUC): `dia_llr` 0.905,
  `n_matched` 0.927, RawScore 0.894, `dia_llr_greedy` 0.887, `dia_llr_soft` 0.790, the #410 key 0.776,
  `dia_angle_soft` 0.739, andes-ms2 angle 0.596, `dia_apex_llr` 0.437.

## Reproduce

```bash
export FASTA=ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta \
       ANDES=andes PERCOLATOR=/path/to/percolator-3.7.1 PY=python3 THREADS=14
./run.sh search run01 LFQ_Astral_DIA_15min_50ng_Human_01.mzML
for s in 42 1 2 3 7; do ./run.sh perc run01 $s; done     # appends one line per seed to RESULTS.txt
python3 eval.py diann report-first-pass.parquet Human_01  # DIA-NN on the same evaluation
python3 scans.py LFQ_Astral_DIA_15min_50ng_Human_01.mzML scans01.parquet   # scan -> RT for gap.py
python3 gap.py run01.pin run01.s42.t.psms run01.s42.d.psms report-first-pass.parquet Human_01 scans01.parquet
```

`PERCOLATOR` can be a wrapper script that runs the Percolator 3.7.1 container
(`exec docker run ... percolator "$@"` or `exec singularity exec ... percolator "$@"`), and `--rescore`
uses it too.
