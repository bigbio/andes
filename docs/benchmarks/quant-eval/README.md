# Quantification evaluation (andes, MaxQuant, FragPipe)

Refs #299.

Python code behind the quantification scorecard of issue #299. It evaluates the **raw** precursor x run quantities of andes,
MaxQuant and FragPipe (IonQuant) through one downstream pipeline that is identical for every tool, so the numbers of the
issue can be reproduced from the three tools' output tables. It runs no search and no quantification: it reads finished
tables.

## Inputs per tool

| tool | file | quantity used | notes |
|---|---|---|---|
| andes | `quant finalize` parquet | `quantity` (raw) | quality columns `quant_q`, `transfer_q`, `library_q`, `measurement_q` (when present) and `evidence_type` of the same file; `quantity_normalized` is not used |
| MaxQuant | `evidence.txt` | `Intensity`, summed over the evidence rows of a precursor x run | reverse and contaminants removed; the LFQ columns are not used |
| FragPipe | IonQuant `combined_ion.tsv` | `<run> Intensity` | IonQuant must be run with normalisation off (`--normalization 0`); rows with `Match Type` MBR are transfers |
| FlashLFQ (PIP-ECHO) | `QuantifiedPeaks.tsv` (`--flashlfq NAME:PEAKS[:PIPQ]`) | `Peak intensity` | key = `Full Sequence` + `/` + `Precursor Charge`, so the identifications given to FlashLFQ must carry andes-form full sequences; `Peak Detection Type` MBR = transfer, its `PIP Q-Value` is kept (`PIPQ` keeps transfers up to that value); decoy-peptide and random-RT peaks are dropped. FlashLFQ sums the charge states of a peptide into one intensity unless run charge-specific |

Tool modified sequences are mapped to the andes key form (`PEPTIDE[+mass]/charge`); only Carbamidomethyl, Oxidation and Acetyl
are mapped, other modifications are dropped and counted. Run names of the tools are matched to the run names of `design.tsv`
(case, `-` and `_` ignored). Every table carries a run tag (hash of the producing file plus a lineage label, `runguard.py`):
joins refuse tables with different tags, duplicates of (precursor, run) are refused, and no column of one arm is joined to
another arm's table.

andes row rules (`--andes NAME:RULE[@LIBQ]:PATH`), all with `library_q <= 0.01` (`@LIBQ` sets another limit for that arm, for runs with a different library entry gate, e.g. `weights@0.05`) and no below-threshold rows:

| rule | rows |
|---|---|
| `ident` | identified rows only: the same run without match between runs |
| `default` | identified rows + transfers with `transfer_q <= 0.01` |
| `all` | identified rows + every transfer the finalize output carries |
| `hard-cut` | identified rows + transfers with `quant_q <= --hard-cut` (default 0.01) |
| `weights` | rows of `all`, with the per-transfer error probability recovered from `measurement_q` (running mean of p over a run's accepted transfers sorted by p, rows with equal p sharing one value; `quant_q` holds that mean in tables without the column); the roll-up then uses learned precision weights (the "V" view) |
| `weights-bias` | rows of `weights`; where the table carries `quantity_bias` and `quantity_bias_se` (log2) for a transfer, its quantity is divided by `2**quantity_bias` and the learned class variance of the transfer is scaled by its `quantity_bias_se**2` relative to its class mean. A table without these columns gives the `weights` rows |

MaxQuant and FragPipe rows are used as the tools report them (their own match-between-runs rows included).

## Pipeline (`quant_eval/`)

1. `readers.py`, `keys.py`: per-tool readers and the key mapping; `runguard.py`: lineage guard.
2. `rollup.py`: per-run median normalisation of log2 quantity over the precursors an arm quantifies in all runs (at least 90 %
   of runs if fewer than 200 such precursors); Tukey median polish per protein over precursors x runs; proteins with at least 2
   precursors are tested. The `weights` view uses weighted medians with weights learned truth-free from the polish residual
   variance by decile of the transfer error probability.
3. `de.py`, `stats.py`: moderated t (empirical Bayes variance prior), Benjamini-Hochberg 5 % per condition pair, true positive
   = significant protein of a species whose amount changes with the observed sign equal to the expected sign; a significant
   protein of a changing species with the opposite sign is a false call (wrong sign), as is a significant protein of the null
   species; counts with and without `|log2FC| >= 0.5`. Comparisons between tools use the common-tested family (the
   (protein, pair) rows that every compared arm tests), so that no arm is counted on more rows than another.
   Three BH families are scored, because the number of tested proteins differs between tools and changes the BH threshold:
   `own` (each arm's own proteins), `inter` (proteins tested by every shared arm; the whole test, variance prior and BH, is run on that set;
   the headline of the scorecard) and `union` (prior fitted on the arm's own proteins; an arm that does not test a protein contributes it with p = 1 to BH). For each family: true positives
   at matched pooled FDP (1, 5, 10 %: the ranked p values of all pairs are swept, hits in the null species count as false) and
   the per-pair AUC of changed versus null species by p value.
4. `metrics.py`: precursor accuracy (leave-one-condition-out log2 error, share with abs(error) > 1, on cells every shared arm
   quantifies and on each arm's own rows, by species and by UPS1 level); counts of precursors, peptides and proteins
   (identified and transfer-only); UPS1 detection and quantification per level; matched-intensity replicate SD at protein and
   precursor level; detection by species and intensity quartile; design-consistent UPS1 detection (a precursor counts at a
   level when one of its rows is within 2x / 3x of its own reference, the median over the levels above the three lowest of
   log2 intensity minus log2 level, after run offsets from the null-species precursors; raw detection is shown next to it);
   absence (rows of a species in runs where it is absent), reported in its own table; rows in blank runs are carry-over and
   are reported in a separate table, not as part of the absence score.
   `rates.py`: truth-measured false-transfer rates for a species that is absent from some runs (PXD014415: yeast in the 20
   human-only runs): rate = transferred rows of the absent species / all transferred rows in those runs; the ratio of those rows
   to the absent species' transferred rows in the runs that contain it (IonQuant-style); a rate scaled by the species' share of
   the identified precursors (an extension, not a measured rate); identified rows of the species in those runs. `depth_curve`
   sweeps `transfer_q` and `quant_q` of the andes transfers and `PIP Q-Value` of FlashLFQ, and the scorecard reads the number of
   transferred rows of an andes arm at the threshold with the most rows among those whose rate is not above a reference tool's rate (depth at a matched
   truth-measured false-transfer rate; rows beyond the finalize `--quant-fdr` limit do not exist, so the curve stops there).
   `rollup.normalise_rt_local` is an optional replacement of the run normalisation (`--rt-local` adds a second DE pass): the
   run shift is smoothed along retention time with a running median over 401 reference precursors, the same function and the
   same RT proxy (median andes apex RT per precursor key) for every arm; keys without an RT get the run's median shift.
   `--norm-null` adds a second view next to the default one: run factors from the precursors of the null species only (the
   species whose amount is constant by design, e.g. yeast in the UPS1 sets, human in the HYE mixes), the same step for every
   arm, then the same roll-up and DE. It uses the design, so it is an evaluation view; the default (all-species) normalisation
   stays the headline of the scorecard.
   `--rollup` picks the protein roll-up, the same for every arm (default `median-polish`). `celltype` removes, before the
   polish, a per-precursor offset of the transferred cells from the identified cells, measured inside the conditions where the
   precursor has both kinds and shrunk to a global offset (empirical Bayes; every scale learned from the arm's own rows);
   `balance` re-weights the transferred cells per (protein, condition) so that their share of the cell weight equals the
   protein's pooled share; `celltype+balance` does both; `censored` treats the transferred cells as right-censored at their
   value (EM inside the polish, E[y | y <= x] with the identified cells' learned residual SD). None is the default: on PXD009740, PXD009815 and ProteoBench QE /
   Astral they did not improve or hold the truth score of every tool (cell type: more null-species calls; balance: fewer true
   positives at matched FDP; censored: more transfers kept at the same error).
5. `downstream.py`: runs everything for one dataset and writes one JSON; `scorecard.py`: markdown tables from the JSONs.

The masked wrong-peak measurement (`masked.py`) needs masked runs of andes (5% of the identified pairs removed from the library); the suite below reads them from the manifest.

## Dataset design

Two small TSV files per dataset (examples in `examples/` for PXD009740, PXD009815 and ProteoBench QE / Astral; run names and
paths are placeholders to edit):

- `design.tsv`: `run`, `condition` (optional), `blank` (optional), `absent_species` (optional).
- `amounts.tsv`: `species`, `condition` (or `*`), `amount` (or `=condition`). The expected log2 fold change between two
  conditions is log2(amount_high) - log2(amount_low); species without an amount are not scored.

The database, the species classes (`HUMAN`, `YEAST`, `ECOLI`, `UPS`) and the null species are command-line arguments.

## Run

```
python3 -m quant_eval.downstream --name DS --design design.tsv --amounts amounts.tsv --fasta db.fasta \
    --classes UPS,YEAST --null YEAST \
    --andes andes:default:finalize.parquet --andes andes_V:weights:finalize.parquet \
    --maxquant MQ:evidence.txt --fragpipe FP:combined_ion.tsv --out DS.json
python3 -m quant_eval.scorecard --results DS1.json DS2.json > scorecard.md
```

Add `--rt-local` to the first command to also score DE after the RT-local normalisation. Absence-only datasets (single-cell blank
run, species absent from runs) need no `--amounts`; see `examples/sc` and `examples/p14`.

`examples/<dataset>/run.sh` holds the full command of each example. Needs Python 3 with numpy, pandas, scipy and pyarrow.
Smoke test on synthetic data (seconds): `python3 tests/test_smoke.py`.

## Metrics and verdict rule

Each dataset section compares the first andes arm with MaxQuant and FragPipe: margin = andes minus tool (absolute and % of
the tool), verdict per tool.

- Counts and errors: `equal` when the margin is within max(3 % of the tool value, 2 units for integer counts); 10 % for FP and
  absence counts; otherwise `better` or `worse` by the metric direction (lower is better for FP, FDP, error and SD).
- FDP rows (FP / (TP + FP)): absolute band of 2 points.
- DE rows and the summary use the intersection family; matched-FDP rows use `equal` within max(3 %, 2).
- DE summary table (TP and FDP of every andes arm against each tool): better = TP >= 0.98 x tool TP and FDP lower by more than
  2 points, or TP > 1.03 x tool TP and FDP not higher by more than 2 points; worse = the mirror; else equal.

## Scoreboard suite, matched true error, gate

One manifest drives the whole truth suite; nothing is searched or quantified by the package (it reads finished tables; building andes
runs stays with the cluster scripts of the user).

```
python3 -m quant_eval.suite --manifest manifest.json --out results/ [--only u9740,qe] [--seeds 42,1,7] [--jobs 3]
python3 -m quant_eval.board --results results/ > scoreboard.md
python3 -m quant_eval.compare --results results/ [--build cand] [--arm V] [--tools MQ,FP] [--seed 42] [--only ups27]
python3 -m quant_eval.gate  --results results/ --candidate cand [--baseline main]
python3 -m quant_eval.proteobench --finalize fin.parquet --design design.tsv --fasta hye.fasta --rule hard-cut --out-prefix sub/QE
```

- `suite.py`: manifest format in its docstring. Per dataset and quant-model seed it runs `downstream` with the arms of every andes build
  (`b`, `b_V`, `b_H` and the matched-cut arm `b_M`), the tool outputs (MaxQuant, FragPipe, FlashLFQ, Sage `lfq.tsv` via `--sage`), the
  identified-cell accuracy block (`accuracy.py`: |err| > 1 share, within-2x, bright-decile gross error, replicate SD), the norm-null view,
  the masked wrong-peak curve (`masked.py`) and the finalize CPU / wall / RSS from `/usr/bin/time -v` files.
- `masked.py`: masked wrong-peak truth of the transfers as a curve over the quant_q limit (pooled, level A / B, runs above 1 %).
- `matched.py`: matched true error. The matched cut of a build is the loosest quant_q limit whose masked wrong-peak rate and whose
  PXD014415 false-transfer rate are not above the reference build's at its own quant_q <= 0.01 plus the reference's own seed-to-seed spread of that rate
  (max - min over the quant-model seeds of the manifest; recorded as `matched_tolerance` in the result JSON, 0 with one seed); depth, accuracy and DE are read at that cut.
- Arms: `X` (transfer_q <= 0.01), `X_V`, `X_H`, `X_I` (identified rows only, the run without match between runs) and `X_M`; a tool's output with match between runs off is another entry of the manifest's `tools` (e.g. `"maxquant": {"MQ": ..., "MQ_noMBR": ...}`).
- `suite --jobs N`: every (dataset, seed) pair is one unit that writes its own result JSON; with `--jobs` above 1 the units run in N worker
  processes. The result files do not depend on `--jobs`; the progress lines of the units are interleaved on stderr.
- `compare.py`: one short table per dataset for one andes build against the tools (default: the last build of the results, and the MaxQuant
  and FragPipe arms with match between runs). One row per metric: depth, cell accuracy, replicate spread, protein CV, DE calls in the
  common family and on the arm's own rows, low spike levels. Per tool its value, the difference of andes to it (relative to the tool; a
  factor from 2x on) and a marker: `✅` andes better, `❌` worse, `≈` equal or the ranges over the quant-model seeds overlap. It reads the
  result JSONs as `board.py` does and computes nothing new. A manifest with one build (that build as `reference_build`) and the tool
  outputs is enough for it. The two rows on the brightest cells depend on the reference build (its quantities, with those of the tools,
  rank the cells), so they can differ between a run with the build as reference and a run with another reference.
- `board.py`: the one scoreboard table. `gate.py`: pass / fail per criterion of the adoption gate (identified rows, wrong-peak, UPS1
  within 2x, DE by the pairwise TP-and-FDP rule, depth at the matched cut, cost); exit status 1 on any failure.
- `proteobench.py`: Custom-input file (and metadata JSON) for the ProteoBench DDA precursor-ion modules.
- Held-out datasets: a manifest entry with `"held_out": true` is refused by `suite` unless `--held-out-eval` is given, and `gate` refuses
  any result flagged held-out; they never enter a gate or tuning run.
- `--hard-cut` accepts `hard-cut=<tau>` as an andes rule (`--andes NAME:hard-cut=0.004:finalize.parquet`).

## Board metrics: low-level UPS1, absence classes, blank split

`board.py` keeps its legacy columns (`--legacy` prints exactly those) and adds, to their right, columns from `lowlevel.py` and `absence.py`
(computed by `downstream.py` into the result JSON, so a rerun of the suite is needed to produce them; `--detail` adds the long tables):

- L1 with its levels named (the amol levels summed, flagged when one is also a reference level, `ref_skip` = 3 lowest levels excluded from the
  references), the UPS1 rows the arm reports at those levels, the share of them within 2x and their median |log2 error|, and the L1 count at a matched
  number of rows (per level, every arm limited to the smallest row count among andes, MaxQuant and FragPipe; the headline is the expected count for rows drawn at random, the best-first reading orders
  andes by its quality column and the other arms by intensity, so its order differs per arm); the rows of precursors without a design reference are printed (+ unscored) and counted as not within 2x in a second precision column. The L1 count alone grows with the number of rows an arm reports; the diagnostic andes rule `every-row`
  (all rows including the below-threshold ones) shows the count-only effect.
- A2: error of the transferred rows (median signed and absolute error, within-2x share) at each level, identified rows next to them.
- PXD014415: the yeast rows in the human-only runs by class of the precursor's own profile (yeast-enriched, not yeast-specific, only in absent runs),
  with the boundary learned from the control species' spread (per bin of absent-run rows, the control count of each bin printed; a bin with fewer than 100 control precursors uses the pooled cutoff) and the fixed 2x next to it. Profile and cutoffs are derived once from a reference table (the andes every-row table) and every arm is classified against it; a precursor missing from the reference is `not_in_reference`. The noise-removed rate is a lower bound of the species-specific errors; the board prints the enriched rate and, in brackets, enriched + the rows without a yeast-run reference.
- Blank run: rows split into identified and transferred, and the carry-over-like rows (below a fraction of the precursor's median in the other runs,
  the fraction learned from the ordinary runs).

- DE: the intersection family is a set of candidate proteins and each arm still tests only those with at least 2 finite values on both sides of a pair,
  so the tested (protein, pair) rows differ between arms. The `common` family (proteins testable in the pair in every match-between-runs arm; `common_all`
  over every shared arm) has identical tested rows; the board prints it next to the legacy column, with the tested-row counts. `TPs` (expected sign)
  is the sign-aware true-positive count, `WS` the wrong-sign calls, and the sign-aware false count is FP + WS (`tp_at_fdp_signed` for the pooled-FDP reading).
  The legacy DE column is unchanged (unsigned TP, intersection family, each arm's own testable rows) and is labelled so on the default board.
- Board cosmetics: `within 2x %` (equal to 100 - `err>1 %`) is dropped from the default board, an empty absence rate prints `-`, FlashLFQ is labelled
  peptide-level (charge states summed) and Sage as having no transfer flag (its MBR fills count as identified). `--legacy` prints the earlier columns and labels.

Tests: `python3 tests/test_board_metrics.py`, `python3 tests/test_compare.py`.

Refs #299.

## Roll-up sensitivity (mokume)

`rollup_sens.py` scores the protein level of every arm after several roll-ups of the same raw precursor intensities: the harness
median polish and the quantification methods of mokume (`rollup_mokume.py`: MaxLFQ, DirectLFQ, Top3, sum, median) under the harness
run normalisation and under mokume's own `globalMedian` normalisation. The tool's own protein intensity (MaxQuant `LFQ intensity`
of `proteinGroups.txt`) is added as a reference row. No roll-up is implemented in the package: mokume-py must be importable
(`pip install mokume-py directlfq`). The default of every other module is unchanged.

```
python3 -m quant_eval.rollup_sens --manifest manifest.json --dataset u9740 --seed 42 --out results/u9740_s42.json
python3 -m quant_eval.rollup_sens_table --results results/ > rollup_sensitivity.md
```
