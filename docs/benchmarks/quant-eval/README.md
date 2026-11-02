# Quantification evaluation (andes, MaxQuant, FragPipe)

Refs #299.

Python code behind the quantification scorecard of issue #299. It evaluates the **raw** precursor x run quantities of andes,
MaxQuant and FragPipe (IonQuant) through one downstream pipeline that is identical for every tool, so the numbers of the
issue can be reproduced from the three tools' output tables. FlashLFQ `QuantifiedPeaks.tsv` (e.g. with PIP-ECHO) can be added
as a reference arm. It runs no search and no quantification: it reads finished tables.

## Inputs per tool

| tool | file | quantity used | notes |
|---|---|---|---|
| andes | `quant finalize` parquet | `quantity` (raw) | quality columns `quant_q`, `transfer_q`, `library_q`, and `evidence_type` of the same file; `quantity_normalized` is not used |
| MaxQuant | `evidence.txt` | `Intensity`, summed over the evidence rows of a precursor x run | reverse and contaminants removed; the LFQ columns are not used |
| FragPipe | IonQuant `combined_ion.tsv` | `<run> Intensity` | IonQuant must be run with normalisation off (`--normalization 0`); rows with `Match Type` MBR are transfers |
| FlashLFQ (optional) | `QuantifiedPeaks.tsv` | `Peak intensity` of the most intense peak of a precursor x run | e.g. FlashLFQ with PIP-ECHO as a reference arm; `Peak Detection Type` MBR / MBR_IsoTrack rows are transfers; decoy peaks dropped; optional MBR q-value cut `--flashlfq-mbr-q` |

Tool modified sequences are mapped to the andes key form (`PEPTIDE[+mass]/charge`); only Carbamidomethyl, Oxidation and Acetyl
are mapped, other modifications are dropped and counted. Run names of the tools are matched to the run names of `design.tsv`
(case, `-` and `_` ignored). Every table carries a run tag (hash of the producing file plus a lineage label, `runguard.py`):
joins refuse tables with different tags, duplicates of (precursor, run) are refused, and no column of one arm is joined to
another arm's table.

andes row rules (`--andes NAME:RULE:PATH`), all with `library_q <= 0.01` and no below-threshold rows:

| rule | rows |
|---|---|
| `default` | identified rows + transfers with `transfer_q <= 0.01` |
| `all` | identified rows + every transfer the finalize output carries |
| `hard-cut` | identified rows + transfers with `quant_q <= --hard-cut` (default 0.01) |
| `weights` | rows of `all`, with the per-transfer error probability recovered from `quant_q` (running mean of p over a run's transfers sorted by p); the roll-up then uses learned precision weights (the "V" view) |

MaxQuant and FragPipe rows are used as the tools report them (their own match-between-runs rows included).

FlashLFQ reader (`--flashlfq NAME:PATH`, columns checked against `ChromatographicPeak.TabSeparatedHeader` and
`MbrChromatographicPeak.ToString` of mzLib `mzLib/FlashLFQ`, master branch, 2026-11):

- required columns: `File Name` (file name without extension = run), `Full Sequence`, `Precursor Charge`, `Peak intensity`,
  `Peak Detection Type`;
- `Full Sequence` in mzLib notation, modification in brackets after its residue (`M[Common Variable:Oxidation on M]`,
  `C[Common Fixed:Carbamidomethyl on C]`), an N-terminal one before residue 1 (`[Common Biological:Acetylation on X]S...`),
  a C-terminal one after `-` at the end. Names Carbamidomethyl / Carbamidomethylation, Oxidation, Acetyl / Acetylation and the
  UniProt `N-acetyl...` names are mapped; numeric tags (`[+15.9949]`), UNIMOD accessions and bare names are accepted too; any
  other modification drops the row (counted). Fixed modifications are printed by mzLib, so no Carbamidomethyl is added;
- `Peak Detection Type` starting with `MBR` (`MBR`, `MBR_IsoTrack`) = transfer; every other type (`MSMS`,
  `MSMSAmbiguousPeakfinding`, `IsoTrack_Ambiguous`) = identified;
- dropped and counted: `Decoy Peptide` True and `Random RT` True rows (FlashLFQ's decoy and random-RT MBR peaks, which
  `QuantifiedPeaks.tsv` carries), ambiguous rows (`|` in `Full Sequence`), rows with `Peak intensity` 0;
- optional MBR q-value column, matched by name ignoring case, spaces and dashes: `PIP Q-Value` (current mzLib), `MBR Q-Value`,
  `PIP-ECHO Q-Value`. With `--flashlfq-mbr-q Q` transfers with a q value above Q are dropped (transfers without a q value, e.g.
  `MBR_IsoTrack`, are kept); without it every MBR row is kept. `QuantifiedPeaks.tsv` lists MBR peaks of every q value (FlashLFQ
  applies its MBR q-value threshold, default 0.05, strict `<`, only in its peptide and protein tables), so pass the threshold
  FlashLFQ was run with to score what FlashLFQ reports;
- several peaks of one precursor x run (different apex): the most intense is kept with its transfer flag, as in FlashLFQ's own
  peptide table and in the FragPipe reader (MaxQuant evidence rows are summed, as before).

FlashLFQ arms count as tool arms: they join the shared cells and the intersection DE family, appear in the DE summary and family
tables and in the per-arm absence tables; the metric-by-metric verdict tables stay andes vs MaxQuant and FragPipe.

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
   absence (rows in blank runs, rows of a species in runs where it is absent, false-transfer rate, blank carry-over; see
   "Absence, false-transfer rate and carry-over" below), reported in its own tables.
   `rollup.normalise_rt_local` is an optional replacement of the run normalisation (`--rt-local` adds a second DE pass): the
   run shift is smoothed along retention time with a running median over 401 reference precursors, the same function and the
   same RT proxy (median andes apex RT per precursor key) for every arm; keys without an RT get the run's median shift.
5. `downstream.py`: runs everything for one dataset and writes one JSON; `scorecard.py`: markdown tables from the JSONs.

Not included: the masked wrong-peak measurement, which needs masked search and finalize runs of andes.

## Dataset design

Two small TSV files per dataset (examples in `examples/` for PXD009740, PXD009815 and ProteoBench QE / Astral; run names and
paths are placeholders to edit):

- `design.tsv`: `run`, `condition` (optional), `blank` (optional), `absent_species` (optional), `order` (optional: integer
  acquisition order, used only for the blank carry-over; empty cells allowed, repeated values refused).
- `amounts.tsv`: `species`, `condition` (or `*`), `amount` (or `=condition`). The expected log2 fold change between two
  conditions is log2(amount_high) - log2(amount_low); species without an amount are not scored.

The database, the species classes (`HUMAN`, `YEAST`, `ECOLI`, `UPS`) and the null species are command-line arguments.

## Run

```
python3 -m quant_eval.downstream --name DS --design design.tsv --amounts amounts.tsv --fasta db.fasta \
    --classes UPS,YEAST --null YEAST \
    --andes andes:default:finalize.parquet --andes andes_V:weights:finalize.parquet \
    --maxquant MQ:evidence.txt --fragpipe FP:combined_ion.tsv \
    --flashlfq FlashLFQ:QuantifiedPeaks.tsv --flashlfq-mbr-q 0.01 --out DS.json
python3 -m quant_eval.scorecard --results DS1.json DS2.json > scorecard.md
```

`--flashlfq` (and `--flashlfq-mbr-q`, set to the MBR / PIP-ECHO q-value FlashLFQ was run with) is optional. Add `--rt-local` to the
first command to also score DE after the RT-local normalisation. Absence-only datasets (single-cell blank run, species absent
from runs) need no `--amounts`; see `examples/sc` and `examples/p14`. Add an `order` column to their `design.tsv` to get the
blank carry-over split.

`examples/<dataset>/run.sh` holds the full command of each example. Needs Python 3 with numpy, pandas, scipy and pyarrow.
Smoke test on synthetic data (seconds): `python3 tests/test_smoke.py`.

## Absence, false-transfer rate and carry-over

Raw absence counts mix depth with error rate, and blank-run counts are dominated by carry-over that no decoy method can catch.
Two derived numbers are reported next to the raw counts (`metrics.absence`, JSON `absence.<arm>`):

- **False-transfer rate** (PIP-ECHO style; Solivais et al., "Improved Detection of Differentially Abundant Proteins through
  FDR-Control of Peptide-Identity-Propagation", J Proteome Res 2025, 24(9):4437-4449), per arm and per absent species S, over the
  runs where `design.tsv` says S is absent. Transfers are the rows each reader marks as transfers (andes `evidence_type` other
  than `identified*`, MaxQuant precursor x run groups made only of MULTI-MATCH / MULTI-MATCH-MSMS evidence rows, FragPipe
  `Match Type` MBR, FlashLFQ `Peak Detection Type` MBR*).
  - visible = transferred rows of S / all transferred rows (every species) in those runs; every transferred S row there is a
    false transfer;
  - f_S = S precursors / all precursors among the arm's transfer-eligible precursors (precursors with at least one identified,
    i.e. non-transfer, row in any run of the arm's own table);
  - scaled = visible / f_S. Assumption: false transfers land on species in proportion to the donor library composition, so
    the visible ones (onto S) are a share f_S of all false transfers and the wrong transfers onto present species are
    invisible. The paper itself counts foreign PIPs directly as observed errors in the pure runs and scales only its entrapment
    count (by database size); the f_S scaling here is the same idea applied to the library composition. The scaled value is
    not capped: above 100 % the assumption does not hold (e.g. transfers of S preferred).
  - counts reported: `transferred_s_rows`, `transferred_rows`, `identified_s_rows`, `eligible_precursors`,
    `eligible_s_precursors`, and `visible_rate`, `f_s`, `scaled_rate` (None when undefined).
  - andes arms: `false_transfer` is on the arm's rows (the reader's rule already selects them, e.g. `transfer_q <= 0.01` for
    `default`); `false_transfer_all` gives the same numbers on every transfer of the finalize output (rule `all`), read
    separately with its own run tag (not repeated for rules `all` and `weights`, whose rows are already every transfer).
- **Blank carry-over** (needs the `order` column): for each blank run, carry-over = every row of a precursor that the same arm
  quantifies (any row) in the run acquired immediately before the blank (the design run with the next lower `order`, which may
  itself be a blank). Reported: blank rows raw (`blank_rows`), `carry_over`, `blank_rows_excl_carry_over`,
  `blank_transferred_excl_carry_over` and the split per blank. Without `order` nothing changes: `carry_over` is null and
  `carry_over_note` (and a line of the scorecard) says it was not computed.

## Metrics and verdict rule

Each dataset section compares the first andes arm with MaxQuant and FragPipe: margin = andes minus tool (absolute and % of
the tool), verdict per tool.

- Counts and errors: `equal` when the margin is within max(3 % of the tool value, 2 units for integer counts); 10 % for FP and
  absence counts; otherwise `better` or `worse` by the metric direction (lower is better for FP, FDP, error and SD).
- FDP rows (FP / (TP + FP)): absolute band of 2 points.
- Absence table: blank rows without carry-over and the false-transfer rates (visible, scaled; margin in points) use the 10 %
  band of the absence counts, lower is better; the per-arm false-transfer and carry-over tables (every arm, FlashLFQ included)
  carry no verdict.
- DE rows and the summary use the intersection family; matched-FDP rows use `equal` within max(3 %, 2).
- DE summary table (TP and FDP of every andes arm against each tool): better = TP >= 0.98 x tool TP and FDP lower by more than
  2 points, or TP > 1.03 x tool TP and FDP not higher by more than 2 points; worse = the mirror; else equal.

Refs #299.
