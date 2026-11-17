# Known gaps in the benchmarks

| gap | size | status / where to fix |
|---|---|---|
| Label-free quantification rests on two ProteoBench datasets and one Percolator seed (three for error control) | best scores use `quantity_normalized`; raw is `quantity`; Astral CV 0.113 above most submissions; MBR entrapment error measured on QExactive only; leaderboard scored by `proteobench` 0.17.0, andes by 0.18.6; several engines' transfer control unchecked | open |
| Entrapment designs disagree | UPS1: natural 1.95–2.07% against shuffle 0.92–0.94%; timsTOF astrocytes: natural 1.17–1.54% against shuffle 0.94–1.12% on one draw (0.70–0.96% on five more) | the shuffle is the primary reading ([Methods](methods.md#method-notes)); the 500-run study has no entrapment |
| Not re-measured at `2ee76dd6` or `2051616b` | glyco deep tier, Comet phospho FDP, `--chimeric` / `--refine` table, non-default `--predict` settings, phospho agreement with MaxQuant | open |
| Single files | phospho is one file, peptide-level only; `--refine` measured on Astral only | open |
| Glyco composition on fucose-rich tissue | heart and lung: backbone agrees with pGlyco2 94–99%, composition 68–78% (97% on liver) | open |
| Missing tool arms | no MaxQuant run for PXD014415 or PXD030780; FragPipe has no MBR-off run; MaxQuant MBR-on rows unfiltered except in the filtered-matches table | open |
| Not benchmarked | iTRAQ, phosphosite localisation, timsTOF quantification; Java MS-GF+ not re-run under current defaults (counts historical) | open |
