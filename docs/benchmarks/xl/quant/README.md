# Cross-link quantification benchmark

The harness of phase P0 in [`docs/design/xl-quant.md`](../../../design/xl-quant.md). It reads
finished quantification tables of cross-links and scores them against designs with a known
answer; it runs no search and no quantification.

andes does not quantify cross-links yet (phase P2 writes that table). Until then the harness
scores the published results of the reference workflow, xiSEARCH + xiFDR + Skyline, which the
authors of two of the datasets deposited with their raw files.

## Datasets

`datasets.json` lists the runs, their design and the published reference outputs.

| Name | Data | Linker | Design | Scored |
|---|---|---|---|---|
| `hsa_dilution` | PXD011036, 15 DDA runs, Orbitrap Fusion Lumos | BS3 | cross-linked HSA at 0.1, 0.3, 0.5, 0.9 and 1×, three injections each | log2 ratio of each level against 1× against the expected one, slope, replicate CV, completeness, residue pairs called changed |
| `hsa_replicates` | PXD007250, 20 DDA runs, Orbitrap Fusion Lumos | BS3 | HSA: 10 injections of one reaction, and 10 separate reactions | CV and completeness per group; differential calls between two halves of a group, all of them false |
| `transfer_clasen` | `cla_b1` and `cla_b4` of [`../datasets.json`](../datasets.json) | DSSO | recombinant interaction groups defined per batch, so each run lacks the other batch's cross-links | false-transfer rate: transferred values of cross-links the run's batch does not contain |
| `transfer_matzinger` | `mat_dsso` and `mat_dsso_hek` of [`../datasets.json`](../datasets.json) | DSSO | peptide library alone, and the library spiked into HEK lysate, same day and instrument | wrong identities among transfers: transferred values whose peptide pair is not a library pair |

Notes on the data:

- **PXD011036 contains two dilution series.** The seven-protein series is DIA only, so the DDA
  series is the HSA one.
- **Diluent not stated.** The published Skyline document names the runs `..._HSA_BS3_10_90_DDA_...`
  for 0.1×, which suggests 10 parts cross-linked HSA in 90 parts of a diluent. The PRIDE record
  does not say what the diluent was, so `normalise` is off and the ratios are read from raw
  quantities.
- **Enzyme.** The community SDRF files of both projects say Lys-C, but the cross-linked peptides
  are tryptic, and the published Skyline document is set to trypsin. The manifest follows the
  document. (The PXD011036 SDRF also gives the organism as horse.)
- **Run order.** The runs of `hsa_replicates` are listed in injection order (the sequence
  numbers in their file names), which the null-call test relies on.
- **The two transfer designs measure different things.**
  - **Clasen:** the two runs hold different real cross-links, so a transferred value of a
    cross-link from the other batch is a false transfer in the usual sense.
  - **Matzinger:** the HEK background carries no cross-links. The truth of the existing
    benchmark is the library in both runs, and no real cross-link is in one run and missing
    from the other. A transferred value outside the library is therefore an identification
    error that match-between-runs carried into another run.
  - Both use cleavable linkers. An MS1 transfer does not depend on how the linker fragments.

## Running

```bash
cd docs/benchmarks/xl/quant
python3 -m xlq.fetch "$XLBENCH_DATA" --only hsa_dilution,hsa_replicates   # 35 raw files, about 92 GB
python3 -m xlq.fetch "$XLBENCH_DATA" --references-only                    # FASTA and published outputs, 14 MB
python3 -m xlq.evaluate --dataset hsa_dilution --data "$XLBENCH_DATA" \
    --arm skyline:skyline:ref/HSA_BS3_DDA_dilution.sky \
    --arm andes:table:/path/to/andes_xl_quant.tsv --out hsa_dilution.json
python3 -m xlq.published --dataset hsa_replicates --data "$XLBENCH_DATA" --out hsa_replicates_published.json
python3 tests/test_xlq.py                                                 # unit tests, seconds
```

The transfer designs read the raw files, FASTAs and truth tables that
[`../fetch.py`](../fetch.py) downloads into the same `$XLBENCH_DATA`. Needs Python 3 with numpy,
pandas, scipy and pyarrow, and openpyxl for the truth tables of the transfer designs. The
statistics come from [`../../quant-eval`](../../quant-eval).

## Input: one long table per tool

`xlq/table.py` defines the table every tool's output is brought to. It has one row per
precursor and run:

- **Required:** `run`, `precursor`, `residue_pair` (andes form, `ACC:pos -- ACC:pos`), `quantity`
  (raw) and `evidence` (`identified` or `transferred`).
- **Optional:** `quant_q`, and the bare peptides and 1-based link sites (`alpha_peptide`,
  `alpha_site`, `beta_peptide`, `beta_site`) that the peptide-library truth needs.

An arm is `NAME:KIND:PATH`:

- **Kinds.** `table` is the long table. `skyline` is a Skyline document whose cross-links are
  pseudo-peptides: peptide A, a spacer residue carrying the linker mass, then peptide B, the
  encoding used before Skyline's native cross-link support. A peak with an area of zero has no
  log2 value and is read as missing; the result counts these peaks.
- **Filters** follow the kind after `?`, joined by `&`:
  - `evidence=identified`: identified rows only;
  - `min_COL=X` and `max_COL=X`: rows whose column `COL` is at least or at most `X`. Rows
    without a value in `COL` are kept.
  - Example: `skyline?min_isotope_dotp=0.9`.

How rows are matched:

- **Runs.** Run names are compared without directory, extension, case, `-` and `_`.
  - Two run names of a table that stand for one run of the design are an error, and so is an
    arm that keeps no row.
  - Runs the design does not have, and runs of the design without a row, are named on stderr
    and in the result (`runs_dropped`, `runs_without_rows`).
- **Residue pairs.** Keys are compared with UniProt identifiers reduced to their accession.
  - **An arm's keys are never changed.** Its own scores and its truth read the positions and
    proteins it reported, so they do not depend on which other arms are scored with it.
  - **Matching is for the `common` view only.** Tools can list different protein positions
    for one peptide: andes lists the positions of the peptide's digestion occurrences, the
    Skyline reader every substring occurrence. Keys of different arms are treated as one
    residue pair when all of the following hold:
    - they share a pair of positions;
    - no arm has two keys in the set;
    - one pair of positions is common to every key of the set. Overlap that only chains
      (A with B, B with C, but A and C share nothing) is not a match.
  - The result says how many sets were matched and how many were ambiguous (`key_matching`).
    Ambiguous sets are not treated as shared.

## Metrics

Every arm goes through the same steps:

- **Roll-up.** Precursors become residue pairs by Tukey median polish of log2 quantities, the
  protein roll-up of quant-eval; `--rollup sum` sums them instead.
- **Normalisation:** only where the dataset allows it (`normalise`). The runs are then centred
  by the rule of quant-eval's `normalise`: run medians over the residue pairs with a value in
  every run, or in at least 90 % of the runs when fewer than 200 have all.
- **Statistics:** quant-eval's moderated t and Benjamini-Hochberg at 5 %.

Each arm is scored on its own residue pairs (`own`) and, with several arms, on the residue pairs
all of them quantify (`common`). The `own` scores of an arm are the same whatever other arms are
given.

- **Dilution.** For a residue pair with at least two values at a level and at 1×, the observed
  log2 ratio is the difference of the replicate means. A reference level without runs, or a
  level with fewer than two runs, is an error. The table reports:
  - its error against log2(level), as the median, the median absolute error and the share
    within 2×;
  - the least-squares slope of observed on expected. 1 means no compression; a shift common
    to every level goes into the intercept, not the slope;
  - the residue pairs called changed with the expected sign.
- **Precision.** Per group of replicate runs:
  - CV on the linear scale, over residue pairs with a value in every run and over those with at
    least two;
  - completeness: the share of residue pairs with a value in every run.

  Both are given at residue-pair and at precursor level.
- **Run offsets.** Per group of replicate runs, the offset of each injection: the median, over
  the precursors, of the run's log2 value minus the precursor's mean over the group. A large
  spread names an injection that differs from its replicates.
- **Null calls.** The runs of a replicate group are split into two halves, 1st, 4th, 5th,
  8th, … against 2nd, 3rd, 6th, 7th, …. A drift that is linear in the injection order cancels
  between them (exactly for a multiple of four runs; a fifth of one injection step remains for
  ten). Every residue pair called changed is a false call.
- **Transfers.** The rate is the transferred rows of residue pairs that the truth excludes
  from a run, divided by all transferred rows in the runs with a truth. The truth is that of
  [`../eval.py`](../eval.py):
  - Clasen run: anything outside that batch's interaction groups;
  - Matzinger run: anything that is not a library pair linked at its designated residues.

  Identified rows of such residue pairs are identification errors and are counted apart. With
  a `quant_q` column the rate is also given at quant_q ≤ 0.01 and ≤ 0.05.

## Reference: published xiSEARCH + xiFDR + Skyline results

Computed on 2026-11-28 from the published files with the commands in
[`../../reproduce.md`](../../reproduce.md#cross-link-quantification-commands-harness-tested-on-synthetic-data).
The result files, which hold the SHA-256 hashes of their inputs, are in
[`../../results/`](../../results): `xl_quant_hsa_dilution.json`,
`xl_quant_hsa_dilution_published.json` and `xl_quant_hsa_replicates_published.json`. The tables
below are generated from them by [`../../render.py`](../../render.py).

### HSA dilution

Median polish, no normalisation, each arm on its own residue pairs. The document is
`HSA_BS3_DDA_dilution.sky` (Skyline 4.2). Errors are in log2; the expected log2 ratios are
−3.32 (0.1×), −1.74 (0.3×), −1.00 (0.5×) and −0.15 (0.9×).

<!-- BEGIN generated: xl-quant-hsa-dilution (docs/benchmarks/render.py) -->
| Skyline values | residue pairs | values | without identification | at every level | slope | error 0.1× | error 0.3× | error 0.5× | error 0.9× | CV at 1× | changed at 0.1× |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| every peak (as published) | 90 | 8,912 | 90 % | 90 | 0.38 | +2.20 | +1.51 | +1.06 | +0.88 | 69 % | 29 of 90 |
| isotope dot product ≥ 0.7 | 90 | 6,435 | 87 % | 79 | 0.74 | +1.04 | +0.74 | +0.75 | +0.68 | 44 % | 47 of 81 |
| isotope dot product ≥ 0.9 | 89 | 3,405 | 82 % | 54 | 0.84 | +0.38 | +0.48 | +0.41 | +0.51 | 30 % | 43 of 57 |
| isotope dot product ≥ 0.95 | 88 | 2,169 | 80 % | 26 | 0.80 | +0.57 | +0.48 | +0.40 | +0.51 | 19 % | 28 of 36 |
| identified peaks only | 90 | 900 | 0 % | 0 | 1.02 | — | +1.39 | +0.77 | +0.48 | 31 % | 0 of 0 |
<!-- END generated: xl-quant-hsa-dilution -->

- **Every peak, as published.** Every residue pair has a value in every run, but the slope is
  far below 1.
  - Nine in ten values come from runs where that precursor was not identified (Skyline's
    `aligned` or `false` peaks). Skyline places these peaks without controlling their error,
    and filtering them by peak quality brings the slope closer to 1 (next point).
- **Isotope dot-product filter.** A filter on Skyline's isotope dot product trades depth for
  accuracy: the slope rises, and the residue pairs with values at every level fall by two
  fifths at 0.9 and by more than two thirds at 0.95.
- **Identified values alone** are too sparse for this design. No residue pair has two at 0.1×,
  and the slope rests on 7, 9 and 21 residue pairs at 0.3×, 0.5× and 0.9×.
- **Positive error at every level, including 0.9×.** An error above +0.15 at 0.9× means the
  0.9× runs read higher than the 1× runs, and they do in every view.
  - This is not one injection that went wrong: the injections of a level agree within a few
    tenths of a log2 (table below).
  - The 1× level as a whole reads low, which shifts every ratio by the same amount. Such a
    shift goes into the intercept of the fit, not into the slope.

<!-- BEGIN generated: xl-quant-hsa-injections (docs/benchmarks/render.py) -->
| level | offsets of its injections (log2) | spread |
|---|---|---:|
| 0.1× | -0.03, +0.16, -0.02 | 0.19 |
| 0.3× | -0.02, +0.20, +0.00 | 0.22 |
| 0.5× | +0.17, +0.17, +0.08 | 0.10 |
| 0.9× | +0.31, -0.02, -0.06 | 0.38 |
| 1× | +0.10, +0.06, -0.05 | 0.16 |
<!-- END generated: xl-quant-hsa-injections -->

What the document holds, and two checks of the reader:

<!-- BEGIN generated: xl-quant-hsa-document (docs/benchmarks/render.py) -->
| | count |
|---|---:|
| cross-linked pseudo-peptides | 447 |
| precursors | 638 |
| peaks | 9,516 |
| peaks with an area of zero (read as missing) | 604 |
| residue pairs in the document | 90 |
| residue pairs of the target-target matches of the xiFDR file | 90 |
| residue pairs in both | 90 |
| site-isomer sets | 32 |
| precursors in site-isomer sets | 64 |
| of these, with a peak that overlaps another site's peak in a run | 60 |
<!-- END generated: xl-quant-hsa-document -->

- **Mapping check.** The residue pairs mapped from the Skyline pseudo-peptides to protein
  positions are exactly those of the target-target matches of the published xiFDR PSM file
  (5 % residue-pair FDR).
- **Site isomers.** Precursors with the same two peptides, modifications and charge but another
  link site are isobaric. Where their peaks overlap in a run, MS1 measures one sum, yet each
  gets it as its own value. The reference numbers above still count such values under their
  residue pairs (design D2; the harness applies the exclusion to every arm in phase P4).

### HSA replicates

The published CV reports, per Skyline transition. The paper reports 14 % and 32 %
(Müller et al., J Am Soc Mass Spectrom 2018). The areas behind these CVs are not published.

<!-- BEGIN generated: xl-quant-hsa-replicates (docs/benchmarks/render.py) -->
| replicates | precursors | median CV, monoisotopic peak | M+1 | M+2 |
|---|---:|---:|---:|---:|
| injection | 518 | 13.8 % | 12.6 % | 11.9 % |
| reaction | 533 | 33.0 % | 32.0 % | 30.6 % |
<!-- END generated: xl-quant-hsa-replicates -->

What andes has to beat is therefore depth close to the every-peak view, accuracy at least that
of the filtered views, and transfers that are counted and error-controlled.

## Still to do

- **Reference workflow runs.** Run xiSEARCH + xiFDR + a current Skyline for `hsa_replicates`
  (its areas are not published) and for the transfer designs, from the raw files. Skyline
  runs on Windows, or in its Docker image under Wine.
- **andes arm.** Phase P2 writes the andes table in the contract above.
- **Site-isomer groups.** A `site_group` column and its rule for every arm (phase P4).
