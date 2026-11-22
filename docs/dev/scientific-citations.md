# Citing papers and ideas, not datasets

## Rule

Credit the **paper** that introduced a method and state the **idea** in one or two plain
sentences. Do not cite a dataset as the source of a method, a model or a design choice.

Good:

> The match-scoring network follows the transfer-learning recipe of Zeng et al.
> (Nat. Commun. 2022): a small MLP trained on run-pair features to tell true from false
> transfers.

Avoid:

> Trained on PXD012345 and PXD067890.

## Format

`Author et al., Journal Year` (or `Author & Author` for two authors), optionally with
volume/article number, as in NOTICE. Add a DOI link in user-facing docs when one exists.
Verify every reference (PubMed or the journal) before writing it; never cite from memory.

## Where things go

- Doc comments (`///`, `//!`): the paper and the idea the function implements, when the
  method is not obvious. No dataset names.
- README, DOCS.md, TRAIN.md, NOTICE: papers and ideas. Training data is described by its
  kind and openness ("annotated MS/MS spectra from public PRIDE reanalyses"), not by a list
  of accessions.
- Dataset accession numbers appear only where a benchmark must be reproduced (configs and
  results under docs/benchmarks/, benchmark tables). There they are identifiers, not
  citations. NOTICE lists no datasets.
- Commit messages and PRs: the design rationale, with the papers it rests on.

## Ideas without a paper

If a design comes from andes's own experiments, say so and describe the measurement that
justified it (in the commit or PR, per CLAUDE.md), rather than attributing it to a paper.

## Registry and markers

`docs/citations.yaml` lists the works the code and docs rest on, as a YAML list of
`{key, authors, title, year, journal, doi | url}`. Keys are lowercase letters, digits and
underscores, for example `kall2007_percolator`; each entry has a `doi` or a `url`.

A doc comment (`///` or `//!`) or a page under `docs/` credits a registry entry with the
marker `[cite:` followed by the key and `]`. A structure test fails when a marker names a key
that is not in the registry, or when an entry is incomplete or its key is duplicated. Add the
entry (checked against the journal or PubMed) before writing the marker.
