---
hide:
  - navigation
  - toc
---

<div class="andes-hero" markdown>

![andes](assets/andes-logo.svg){ .andes-logo }

**The data-driven peptide search engine of the quantms ecosystem.**
Spectra and a FASTA go in; a Percolator-ready PIN comes out. The scoring model is chosen per
file from its metadata.

[Get started](getting-started/install.md){ .md-button .md-button--primary }
[Benchmarks](benchmarks/index.md){ .md-button }
[Parameters](parameters.md){ .md-button }

</div>

## What it does

<div class="grid cards" markdown>

-   **Reads your files as they are**

    mzML, MGF, native Thermo `.raw` and Bruker timsTOF `.d`. One of 17 bundled scoring models
    is picked from the file's activation, resolution and isobaric label.
    [How it works](how-it-works.md)

-   **Leaves the FDR to Percolator**

    andes computes no production FDR. It writes the PIN, and Percolator, run by you or by
    `--rescore`, decides what passes `q ≤ 0.01`.
    [Quick start](getting-started/quickstart.md)

-   **Optional modes**

    Co-isolated peptides (`--chimeric`), secondary modifications (`--refine`), and intact
    N-glycopeptides (`--glyco`).
    [Modes](modes.md)

-   **Quantification built in**

    TMT/iTRAQ reporter ions (`--tmt`) and label-free MS1 features with match between runs
    (`--lfq`).
    [Reference §10](reference.md#10-quantification)

</div>

## Identifications, measured the same way

PSMs at Percolator `q ≤ 0.01` on public reference datasets, rescored identically for every
engine on the same 8-core machine. The conditions, dates and commits behind each bar are on the
[benchmarks page](benchmarks/index.md).

![PSMs at q ≤ 0.01 for andes, Comet and Java MS-GF+](figures/psms.svg)

## Install

```bash
cargo install --git https://github.com/bigbio/andes --bin andes
```

Release archives for Linux, macOS and Windows are on the
[releases page](https://github.com/bigbio/andes/releases). Details: [Install](getting-started/install.md).
