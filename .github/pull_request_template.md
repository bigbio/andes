## What and why

<!-- What the change does and the reason for it. -->

## Results changed?

- [ ] No: identifications, scores, quantities and output files are unchanged.
- [ ] Yes: describe the change below.

If yes, give the datasets, the seeds, the numbers before and after, and the evidence the rule in
[`docs/dev/benchmark-reproducibility.md`](../docs/dev/benchmark-reproducibility.md) asks for.

## Checklist

- [ ] Goldens are unchanged, or regenerated with the reason stated above.
- [ ] No ratchet baseline was raised.
- [ ] Every new parameter is declared with its doc, default and unit.
- [ ] New algorithms carry citations in `docs/citations.yaml`.
- [ ] No site identifiers (hostnames, cluster paths, home paths) in committed files.
- [ ] For a stacked PR: the body has a retarget note naming the base branch.
