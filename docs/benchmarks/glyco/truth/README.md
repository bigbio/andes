# Reference identification tables

Published glycopeptide identifications from other engines in one canonical format, so a
benchmark can be scored without re-downloading multi-GB originals or parsing proprietary
formats.

| file | source | spectra | regime |
|---|---|---|---|
| `pglyco2_mouse_liver.tsv.gz` | pGlyco2, PRIDE PXD005553 | 17,855 | mouse liver, HCD |
| `msfragger_mouse_liver.tsv.gz` | MSFragger-Glyco (Philosopher-filtered `psm.tsv` deposited in PRIDE PXD031032, a re-analysis of the pGlyco2 raws) | 14,626 | the same five liver fractions, so andes can be scored against two independent engines on identical spectra |
| `pglyco2_mouse_lung.tsv.gz` | pGlyco2, PRIDE PXD005555 | 15,016 | mouse lung, HCD |
| `pglyco2_mouse_heart.tsv.gz` | pGlyco2, PRIDE PXD005413 | 5,383 | mouse heart, HCD |

Columns: `run, scan, charge, peptide, glycan, glycosite`; `peptide` is the bare uppercase
backbone. Each file is **one row per (run, scan)**, decoy- and FDR-filtered, and its header
records the filter. Two traps this avoids: a raw MSFragger `psm.tsv` is pre-FDR (about a third
of its rank-1 glyco rows are its own decoys), and a per-(PSM, protein) export counts shared
peptides several times (that inflated one earlier truth set by 70%). From the MSFragger table
only rows parsing purely as HexNAc/Hex/dHex/NeuAc/NeuGc are kept, first candidate wins.

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

Another engine's **claims at its own stated FDR**, not ground truth; they have not been
entrapment-validated here. Read "% recovered" as agreement with a strong reference, not as
sensitivity. Derived from data deposited in PRIDE under the accessions above; cite the original
publications when quoting a comparison.
