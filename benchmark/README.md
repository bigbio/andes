# Benchmarks

Public, reproducible benchmark **reports** (andes vs Java MS-GF+ and a comparison engine, uniform
Percolator, entrapment-FDP validation) live under [`docs/benchmarks/`](../docs/benchmarks/).

This directory holds only the local **sign-off harness**; heavy inputs (`data/`,
`results/`, prebuilt JARs) are gitignored and not distributed.

## Reference dataset

| Dataset | PXD | Instrument | Type | FASTA / SDRF |
|---------|-----|------------|------|--------------|
| LFQ (Yeast+UPS1) | PXD001819 | LTQ Orbitrap Velos, CID | LFQ | [quantms-test-datasets](https://github.com/bigbio/quantms-test-datasets) |

## VM sign-off harness (`vm/`)

The three-dataset andes-vs-field comparison harness lives under [`vm/`](vm/) — run it on
the self-hosted bench VM. See [`vm/README.md`](vm/README.md).
