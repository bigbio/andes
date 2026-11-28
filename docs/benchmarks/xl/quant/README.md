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
| `hsa_replicates` | PXD007250, 20 DDA runs, Orbitrap Fusion Lumos | BS3 | HSA: 10 injections of one reaction, and 10 separate reactions | CV and completeness per group; differential calls between interleaved halves of a group, all of them false |
| `transfer_matzinger` | `mat_dsso` and `mat_dsso_hek` of [`../datasets.json`](../datasets.json) | DSSO | peptide library alone, and the library spiked into HEK lysate, same day and instrument | false-transfer rate into the library-only run |
| `transfer_clasen` | `cla_b1` and `cla_b4` of [`../datasets.json`](../datasets.json) | DSSO | recombinant interaction groups defined per batch | false-transfer rate into each batch's run |

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
- **False-transfer designs use cleavable linkers.** An MS1 transfer does not depend on how the
  linker fragments, so they still test transfers.

## Running

```bash
cd docs/benchmarks/xl/quant
python3 -m xlq.fetch "$XLBENCH_DATA" --only hsa_dilution,hsa_replicates   # 35 raw files, about 92 GB
python3 -m xlq.fetch "$XLBENCH_DATA" --references-only                    # FASTA and published outputs, 14 MB
python3 -m xlq.evaluate --dataset hsa_dilution --data "$XLBENCH_DATA" \
    --arm skyline:skyline:ref/HSA_BS3_DDA_dilution.sky \
    --arm andes:table:/path/to/andes_xl_quant.tsv --out hsa_dilution.json
python3 -m xlq.published --dataset hsa_replicates --data "$XLBENCH_DATA" --out hsa_replicates.published.json
python3 tests/test_xlq.py                                                 # unit tests, seconds
```

The false-transfer designs read the raw files, FASTAs and truth tables that
[`../fetch.py`](../fetch.py) downloads into the same `$XLBENCH_DATA`. Needs Python 3 with numpy,
pandas, scipy and pyarrow. The statistics come from [`../../quant-eval`](../../quant-eval).

## Input: one long table per tool

`xlq/table.py` defines the table every tool's output is brought to. It has one row per
precursor and run:

- **Required:** `run`, `precursor`, `residue_pair` (andes form, `ACC:pos -- ACC:pos`), `quantity`
  (raw) and `evidence` (`identified` or `transferred`).
- **Optional:** `quant_q`, and the bare peptides and 1-based link sites (`alpha_peptide`,
  `alpha_site`, `beta_peptide`, `beta_site`) that the peptide-library truth needs.

Residue-pair keys are compared with UniProt identifiers reduced to their accession. Run names
are compared without directory, extension, case, `-` and `_`.

An arm is `NAME:KIND:PATH`:

- **Kinds.** `table` is the long table. `skyline` is a Skyline document whose cross-links are
  pseudo-peptides: peptide A, a spacer residue carrying the linker mass, then peptide B, the
  encoding used before Skyline's native cross-link support.
- **Filters** follow the kind after `?`, joined by `&`:
  - `evidence=identified`: identified rows only;
  - `min_COL=X` and `max_COL=X`: rows whose column `COL` is at least or at most `X`. Rows
    without a value in `COL` are kept.
  - Example: `skyline?min_isotope_dotp=0.9`.

## Metrics

Every arm goes through the same steps:

- **Roll-up.** Precursors become residue pairs by Tukey median polish of log2 quantities, the
  protein roll-up of quant-eval; `--rollup sum` sums them instead.
- **Normalisation:** only where the dataset allows it.
- **Statistics:** quant-eval's moderated t and Benjamini-Hochberg at 5 %.

Each arm is scored on its own residue pairs and, with several arms, on the residue pairs all of
them quantify (`common`).

- **Dilution.** For a residue pair with at least two values at a level and at 1×, the observed
  log2 ratio is the difference of the replicate means. The table reports:
  - its error against log2(level), as the median, the median absolute error and the share
    within 2×;
  - the least-squares slope of observed on expected (1 means no compression);
  - the residue pairs called changed with the expected sign.
- **Precision.** Per group of replicate runs:
  - CV on the linear scale, over residue pairs with a value in every run and over those with at
    least two;
  - completeness: the share of residue pairs with a value in every run.

  Both are given at residue-pair and at precursor level.
- **Null calls.** The runs of a replicate group are split alternately (1st, 3rd, … against
  2nd, 4th, …), so drift over the injection order is not a difference between the halves.
  Every residue pair called changed is a false call.
- **False transfers.** The rate is the transferred rows of residue pairs that cannot be in a
  run, divided by all transferred rows in the runs with a truth. The truth is that of
  [`../eval.py`](../eval.py):
  - library-only run: anything that is not a library pair linked at its designated residues;
  - Clasen run: anything outside that batch's interaction groups.

  Identified rows of such residue pairs are identification errors and are counted apart. With
  a `quant_q` column the rate is also given at quant_q ≤ 0.01 and ≤ 0.05.

## Reference: published xiSEARCH + xiFDR + Skyline results

Computed on 2026-11-28 from the published files, whose SHA-256 hashes are in the result files:

- [`results/hsa_dilution.skyline_published.json`](results/hsa_dilution.skyline_published.json):
  `xlq.evaluate` with the five arms below on `HSA_BS3_DDA_dilution.sky` (Skyline 4.2, 447
  cross-linked peptides, 638 precursors, 15 runs). It takes 3 s.
- [`results/hsa_replicates.published.json`](results/hsa_replicates.published.json):
  `xlq.published` on the two CV reports.

**HSA dilution.** Median polish, no normalisation, each arm on its own residue pairs:

| Skyline values | residue pairs | values | without identification | at every level | slope | error 0.1× | 0.3× | 0.5× | 0.9× | CV at 1× | changed at 0.1× |
|---|---|---|---|---|---|---|---|---|---|---|---|
| every peak (as published) | 90 | 8,912 | 90 % | 90 | 0.38 | +2.20 | +1.51 | +1.06 | +0.88 | 69 % | 29 of 90 |
| isotope dot product ≥ 0.7 | 90 | 6,435 | 87 % | 79 | 0.74 | +1.04 | +0.74 | +0.75 | +0.68 | 44 % | 47 of 81 |
| isotope dot product ≥ 0.9 | 89 | 3,405 | 82 % | 54 | 0.84 | +0.38 | +0.48 | +0.41 | +0.51 | 30 % | 43 of 57 |
| isotope dot product ≥ 0.95 | 88 | 2,169 | 80 % | 26 | 0.80 | +0.57 | +0.48 | +0.40 | +0.51 | 19 % | 28 of 36 |
| identified peaks only | 90 | 900 | 0 % | 0 | 1.02 | — | +1.39 | +0.77 | +0.48 | 31 % | 0 of 0 |

Errors are in log2. The expected log2 ratios are −3.32 (0.1×), −1.74 (0.3×), −1.00 (0.5×) and −0.15 (0.9×).

- **Every peak, as published.** Every residue pair has a value in every run, but ratios are
  compressed to 0.38 of the expected change.
  - 90 % of the values come from runs where that precursor was not identified (Skyline's
    `aligned` or `false` peaks). Skyline places these peaks without controlling their error,
    and filtering them by peak quality reduces the compression (next point).
  - The 90 residue pairs, mapped from the Skyline pseudo-peptides to protein positions, are
    exactly the 90 of the target-target PSMs in the published xiFDR PSM file (5 % residue-pair
    FDR). This checks the reader's mapping.
- **Isotope dot-product filter.** A filter on Skyline's isotope dot product trades depth for
  accuracy: at 0.9 the slope is 0.84, but residue pairs with values at every level fall from 90
  to 54.
- **Identified values alone** are too sparse for this design. No residue pair has two at 0.1×,
  and the slope rests on 7, 9 and 21 residue pairs at 0.3×, 0.5× and 0.9×.
- **Positive error at every level, including 0.9×.** On the five precursors with a dot product
  of at least 0.9 in all 15 runs, the median intensity of the three 1× injections spans
  2.6 log2. Run-to-run variation is large in this set, so these errors are not compression
  alone.

**HSA replicates.** The published CV reports give, for the monoisotopic peak, median CV
13.8 % over 518 precursors for the injection replicates and 33.0 % over 533 for the reaction
replicates. The paper reports 14 % and 32 % (Müller et al., J Am Soc Mass Spectrom 2018). The
areas behind these CVs are not published.

What andes has to beat is therefore depth close to the every-peak view, accuracy at least that
of the filtered views, and transfers that are counted and error-controlled.

## Still to do

- **Reference workflow runs.** Run xiSEARCH + xiFDR + a current Skyline for `hsa_replicates`
  (its areas are not published) and for the false-transfer designs, from the raw files. Skyline
  runs on Windows, or in its Docker image under Wine.
- **andes arm.** Phase P2 writes the andes table in the contract above.
- **Generated tables.** A `render.py` block for these tables once andes numbers sit next to
  the reference (phase P5).
