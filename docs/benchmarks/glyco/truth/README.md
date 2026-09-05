# Reference identification tables

Published glycopeptide identifications from other engines, converted to one canonical
format and committed here so a benchmark can be scored **without re-downloading the
multi-GB originals or re-parsing a proprietary format**.

| file | source | spectra | regime |
|---|---|---|---|
| `pglyco2_mouse_liver.tsv.gz` | pGlyco2, PRIDE PXD005553 | 17,855 | mouse liver, HCD |
| `pglyco2_mouse_lung.tsv.gz` | pGlyco2, PRIDE PXD005555 | 15,016 | mouse lung, HCD |
| `pglyco2_mouse_heart.tsv.gz` | pGlyco2, PRIDE PXD005413 | 5,383 | mouse heart, HCD |
| `byonic_plasma.tsv.gz` | Byonic (depositors'), PRIDE PXD030622 | 629 | human plasma, sceHCD |

372 KB for all four. Columns: `run, scan, charge, peptide, glycan, glycosite`.
`peptide` is the BARE backbone — uppercase, no modifications.

## The filter is in the file

Each table's header records the decoy and FDR filter that produced it, because getting
those wrong has produced wrong numbers here twice:

- a raw MSFragger `psm.tsv` is **pre-FDR**, and roughly a third of its rank-1 glyco rows
  are that engine's own decoys;
- Byonic's `PQMs` table emits **one row per (PSM, protein)**, so a peptide in several
  proteins is counted several times — that inflated one truth set from 629 spectra to
  1,068 rows.

Both are handled: every table is **one row per (run, scan)**, decoy-filtered, FDR-filtered.

## Regenerating

```bash
python3 ../make_truth.py pglyco2 MouseLiver-Z-T-*-FDR.txt | gzip -9 > pglyco2_mouse_liver.tsv.gz
python3 ../make_truth.py byonic  *.byrslt                 | gzip -9 > byonic_plasma.tsv.gz
python3 ../make_truth.py strucgp *_result.xlsx            | gzip -9 > strucgp_*.tsv.gz
```

## Reading one

```python
import csv, gzip
with gzip.open("pglyco2_mouse_liver.tsv.gz", "rt") as fh:
    rows = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
```

## What these are and are not

They are another engine's **claims at its own stated FDR**, not ground truth. pGlyco2's 1%
and Byonic's 1% are each that tool's estimate; neither has been independently
entrapment-validated here. Treat "% recovered" as agreement with a strong reference, not as
sensitivity. Disagreement can mean andes is wrong, the reference is wrong, or both are
defensible for an ambiguous spectrum.

Provenance: derived from data publicly deposited in PRIDE under the accessions above; cite
the original publications when quoting a comparison.
