# PXD001819 benchmark inputs (public)

Used by `benchmark/ci/PXD001819/run_ci.sh` and documented in  
[`.claude/plans/pxd001819-ci-benchmark.md`](../../../../../.claude/plans/pxd001819-ci-benchmark.md).

## Public data (downloaded in CI, not committed)

| File | URL |
|------|-----|
| `UPS1_5000amol_R1.mzML.gz` | `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/lfq/LTQOrbitrapVelos/PXD001819/UPS1_5000amol_R1.mzML.gz` |
| `PXD001819_uniprot_yeast_ups.fasta` | `https://raw.githubusercontent.com/bigbio/quantms-test-datasets/quantms/databases/PXD001819_uniprot_yeast_ups.fasta` |

## Mods

`mods.txt` in this directory matches the LFQ yeast + UPS1 setup (fixed Carbamidomethyl C, variable Oxidation M, variable Acetyl protein N-term).

## Search flags (must match `benchmark/run_pxd001819_benchmark.sh`)

Trypsin/P (`-e 1`), CID (`-m 0`), LTQ Orbitrap (`-inst 0`), 5 ppm precursor (`-t 5ppm`), isotope `0,1`, TDA on, MS2 only (`-msLevel 2`), `-addFeatures 1`, `-n 1`, missed cleavages 2, charges 2–4.
