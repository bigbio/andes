# Quantification evaluation (andes, MaxQuant, FragPipe)

Refs #299.

Python code behind the quantification scorecard of issue #299. It evaluates the **raw** precursor x run quantities of andes,
MaxQuant and FragPipe (IonQuant) through one downstream pipeline that is identical for every tool, so the numbers of the
issue can be reproduced from the three tools' output tables. It runs no search and no quantification: it reads finished
tables.

## Inputs per tool

| tool | file | quantity used | notes |
|---|---|---|---|
| andes | `quant finalize` parquet | `quantity` (raw) | quality columns `quant_q`, `transfer_q`, `library_q`, and `evidence_type` of the same file; `quantity_normalized` is not used |
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
| `default` | identified rows + transfers with `transfer_q <= 0.01` |
| `all` | identified rows + every transfer the finalize output carries |
| `hard-cut` | identified rows + transfers with `quant_q <= --hard-cut` (default 0.01) |
| `weights` | rows of `all`, with the per-transfer error probability recovered from `quant_q` (running mean of p over a run's transfers sorted by p); the roll-up then uses learned precision weights (the "V" view) |

MaxQuant and FragPipe rows are used as the tools report them (their own match-between-runs rows included).

## Pipeline (`quant_eval/`)

1. `readers.py`, `keys.py`: per-tool readers and the key mapping; `runguard.py`: lineage guard.
2. `rollup.py`: per-run median normalisation of log2 quantity over the precursors an arm quantifies in all runs (at least 90 %
   of runs if fewer than 200 such precursors); Tukey median polish per protein over precursors x runs; proteins with at least 2
   precursors are tested. The `weights` view uses weighted medians with weights learned truth-free from the polish residual
   variance by decile of the transfer error probability.
3. `de.py`, `stats.py`: moderated t (empirical Bayes variance prior), Benjamini-Hochberg 5 % per condition pair, true positive
   = significant protein of a species whose amount changes, with the observed sign equal to the expected sign counted
   separately; false positive = significant protein of the null species; counts with and without `|log2FC| >= 0.5`.
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
   transferred rows of an andes arm at the loosest threshold whose rate is not above a reference tool's rate (depth at a matched
   truth-measured false-transfer rate; rows beyond the finalize `--mbr-fdr` limit do not exist, so the curve stops there).
   `rollup.normalise_rt_local` is an optional replacement of the run normalisation (`--rt-local` adds a second DE pass): the
   run shift is smoothed along retention time with a running median over 401 reference precursors, the same function and the
   same RT proxy (median andes apex RT per precursor key) for every arm; keys without an RT get the run's median shift.
5. `downstream.py`: runs everything for one dataset and writes one JSON; `scorecard.py`: markdown tables from the JSONs.

Not included: the masked wrong-peak measurement, which needs masked search and finalize runs of andes.

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

Refs #299.
