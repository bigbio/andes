# Benchmarks

andes is measured on public PRIDE datasets for identification and quantification. Identifications
are counted after Percolator at `q ≤ 0.01`. Where the database carries an entrapment component,
we also report the **true** false discovery proportion (FDP), which says whether the nominal 1%
holds. Comparison engines are named with their version, and their numbers come from runs on the
same machine with the same rescoring, or from a public leaderboard.

This page is a visual summary. Dates, commits, database definitions and every row are in the
[full benchmark report](details.md); [how to reproduce](reproduce.md) lists the commands.

## Identifications

![PSMs at q ≤ 0.01](../figures/psms.svg)

Java MS-GF+ was measured on UPS1 only, in an earlier session, and not re-run. UPS1 andes is the
`2ee76dd6` run of the report (2026-11-07, entrapment database, another host; bar at the midpoint
of its range). The phospho set is searched against mouse proteins with a whole-protein shuffled
entrapment, the database both engines used.

## Time, CPU and memory

![andes as a percentage of Comet](../figures/cost_relative.svg)

The same measurements in absolute terms:

![Wall time, CPU time and peak memory](../figures/cost.svg)

Wall times vary 5–15% between sessions; CPU time is steadier. Bars for andes are the midpoint
of the reported range. UPS1 is left out: its andes run was on another host.

## Error control: does the nominal 1% hold?

![True FDP at nominal 1%](../figures/true_fdp.svg)

Each bar spans the entrapment-estimated true FDP over Percolator seeds and database
realisations (timsTOF: per-run means over six shuffle draws). The peptide-level shuffle is the
primary design; the foreign (*Arabidopsis*) and natural (*E. coli*, yeast) entrapments are
cross-checks. On UPS1 the two designs disagree: the natural *E. coli* entrapment reads 1.95–2.07%,
about half of it on conserved proteins that may really be in the sample.
Read the [method notes](methods.md#method-notes) before comparing across rows.

## Optional modes

![Chimeric mode on Astral](../figures/chimeric.svg)

`--chimeric` reports second and later peptides per spectrum. Counts are at an entrapment true
FDP of 1% on Astral, measured on an earlier build.

## Also benchmarked

- **Glycoproteomics** against five published tools on PXD065680:
  [glyco benchmark](glyco.md).
- **Quantification** (label-free with match between runs against the ProteoBench leaderboard, TMT
  reporter-ion accuracy): [quantification](../quantification.md).
- timsTOF DDA-PASEF and the mouse-liver glyco set are in the [identification report](identification.md).
