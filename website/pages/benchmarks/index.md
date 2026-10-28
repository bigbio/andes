# Benchmarks

andes is measured on public PRIDE datasets for identification and quantification. Identifications
are counted after Percolator at `q ≤ 0.01`. Where the database carries an entrapment component,
we also report the **true** false discovery proportion (FDP), which says whether the nominal 1%
holds. Comparison engines are named with their version, and their numbers come from runs on the
same machine with the same rescoring, or from a public leaderboard.

This page is a visual summary. Dates, commits, database definitions and every row are in the
[full benchmark report](details.md); [how to reproduce](details.md#3-how-to-reproduce) lists
the commands.

## Identifications

![PSMs at q ≤ 0.01](../figures/psms.svg)

Java MS-GF+ was not re-run for this comparison: its counts are from an earlier session, so it
is shown only where one exists. The phospho set is searched against mouse proteins with a
whole-protein shuffled entrapment, the database both engines used.

## Time, CPU and memory

![Wall time, CPU time and peak memory](../figures/cost.svg)

Wall times vary 5–15% between sessions; CPU time is steadier. Wall bars for andes are the
midpoint of the reported range.

## Error control: does the nominal 1% hold?

![True FDP at nominal 1%](../figures/true_fdp.svg)

Each bar spans the entrapment-estimated true FDP over Percolator seeds and database
realisations. The peptide-level shuffle is the primary design; the foreign (*Arabidopsis*) and
natural (*E. coli*, yeast) entrapments are cross-checks. On UPS1 the two designs disagree, and
the natural *E. coli* entrapment is the established reference there, which reads about 2.1%.
Read the [method notes](details.md#4-method-notes) before comparing across rows.

## Optional modes

![Chimeric mode on Astral](../figures/chimeric.svg)

`--chimeric` reports second and later peptides per spectrum. Counts are at an entrapment true
FDP of 1% on Astral, measured on an earlier build.

## Also benchmarked

Quantification (label-free with match between runs against the ProteoBench leaderboard, TMT
reporter-ion accuracy), timsTOF DDA-PASEF, and glycopeptides are in the
[full report](details.md#1-summary).
