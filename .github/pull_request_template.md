## What and why

<!-- What the change does and the reason for it. -->

## Results changed?

- [ ] No: identifications, scores, quantities and output files are unchanged.
- [ ] Yes: describe the change below.

If yes, give the datasets, the seeds, the numbers before and after, and the evidence the rules in
`docs/dev/benchmark-reproducibility.md` and `docs/dev/fdr-validation.md` (paths from the
repository root) ask for.

## Checklist

- [ ] Goldens are unchanged, or regenerated with the reason stated above.
- [ ] No ratchet baseline was raised.
- [ ] Every new parameter is declared with its doc, default and unit.
- [ ] New algorithms carry citations in `docs/citations.yaml`.
- [ ] No site identifiers (hostnames, cluster paths, home paths) in committed files.
- [ ] Logs: the clean-logs rules are followed and a default run's log fits one screen.
- [ ] Run report: `statistics.log` is updated if counts or parameters change.
- [ ] I/O and disk: no new temporary-file leaks; new outputs are versioned Parquet.
- [ ] Docs: `ARCHITECTURE.md` and `DOCS.md` are updated where the change touches them.
- [ ] Tests: end-to-end tests or goldens cover new behaviour.
- [ ] For a stacked PR: the body has a retarget note naming the base branch.
