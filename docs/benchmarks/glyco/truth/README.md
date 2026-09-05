# Reference identification tables

Published glycopeptide identifications from other engines, converted to one canonical
format and committed here so a benchmark can be scored **without re-downloading the
multi-GB originals or re-parsing a proprietary format**.

| file | source | spectra | regime |
|---|---|---|---|
| `pglyco2_mouse_liver.tsv.gz` | pGlyco2, PRIDE PXD005553 | 17,855 | mouse liver, HCD |
| `msfragger_mouse_liver.tsv.gz` | MSFragger-Glyco (Philosopher-filtered `psm.tsv` deposited in PRIDE PXD031032, a re-analysis of the pGlyco2 raws) | 14,626 | the same five liver fractions, so andes can be scored against two independent engines on identical spectra |
| `pglyco2_mouse_lung.tsv.gz` | pGlyco2, PRIDE PXD005555 | 15,016 | mouse lung, HCD |
| `pglyco2_mouse_heart.tsv.gz` | pGlyco2, PRIDE PXD005413 | 5,383 | mouse heart, HCD |

372 KB for all four. Columns: `run, scan, charge, peptide, glycan, glycosite`.
`peptide` is the BARE backbone — uppercase, no modifications.

## The filter is in the file

Each table's header records the decoy and FDR filter that produced it, because getting
those wrong has produced wrong numbers here twice:

- a raw MSFragger `psm.tsv` is **pre-FDR**, and roughly a third of its rank-1 glyco rows
  are that engine's own decoys;
- a per-(PSM, protein) export counts a shared peptide several times — that inflated one
  earlier truth set by 70%;
- the deposited MSFragger table mixes N-glycans with chemical modifications and isotope
  labels in the same column (12,187 of the liver rows are `Label:15N`, `Deamidated`,
  `Phos`, `Kdn`...), and lists several candidate compositions for 7% of rows. Only rows
  parsing purely as HexNAc/Hex/dHex/NeuAc/NeuGc are kept, first candidate wins.

Both are handled: every table is **one row per (run, scan)**, decoy-filtered, FDR-filtered.

## Regenerating

```bash
python3 ../make_truth.py pglyco2 MouseLiver-Z-T-*-FDR.txt | gzip -9 > pglyco2_mouse_liver.tsv.gz
python3 ../make_truth.py strucgp *_result.xlsx            | gzip -9 > strucgp_*.tsv.gz
# PXD031032/Mouse_OpenSearch_6000Da_N-GlycanMode_psm.tsv (109 MB) covers five tissues; keep one
python3 ../make_truth.py msfragger Mouse_OpenSearch_6000Da_N-GlycanMode_psm.tsv MouseLiver | gzip -9 > msfragger_mouse_liver.tsv.gz
```

## Reading one

```python
import csv, gzip
with gzip.open("pglyco2_mouse_liver.tsv.gz", "rt") as fh:
    rows = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
```

## What these are and are not

They are another engine's **claims at its own stated FDR**, not ground truth. pGlyco2's 1% is
that tool's estimate; it has not been independently entrapment-validated here. Treat "% recovered" as agreement with a strong reference, not as
sensitivity. Disagreement can mean andes is wrong, the reference is wrong, or both are
defensible for an ambiguous spectrum.

Provenance: derived from data publicly deposited in PRIDE under the accessions above; cite
the original publications when quoting a comparison.
