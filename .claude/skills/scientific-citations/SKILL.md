---
name: scientific-citations
description: How andes credits prior work in docs, NOTICE, doc comments, commit messages and PRs - cite papers and describe the idea, not datasets. Use when writing documentation, README/DOCS/TRAIN text, NOTICE, model descriptions, or a doc comment that explains where a method comes from.
---

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
