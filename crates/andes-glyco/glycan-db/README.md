# Bundled pGlyco glycan databases

The pGlyco N-glycan structure databases
([Liu et al., *Nat. Commun.* 2017](https://www.nature.com/articles/ncomms15473)), bundled so
`--glyco-species` works without a pGlyco installation. Each file is a header line of
monosaccharide symbols, then one canonical S-expression per glycan (grammar:
`andes-glyco::glycan_db::load_glycan_gdb`).

| File | Species / scope | Structures | Provenance |
|---|---|---|---|
| `pGlyco-N-HighMannose.gdb` | high-mannose N-glycans (species-independent) | 30 | derived (curated subset) |
| `pGlyco-N-Human.gdb` | *Homo sapiens* N-glycans | 2922 | verbatim |
| `pGlyco-N-Human-multi.gdb` | human, multi-antennary | 5634 | verbatim |
| `pGlyco-N-Mouse.gdb` | *Mus musculus* N-glycans | 6662 | verbatim |
| `pGlyco-N-Mouse-large.gdb` | mouse, extended | 7878 | verbatim |

Plant databases are not bundled: xylose (`X`) is not yet represented in `GlycanComp`.

Attribution: redistributed from the pGlyco project (https://github.com/pFindStudio/pGlyco3),
Apache-2.0. If you redistribute andes with these files, keep this notice and the upstream
Apache-2.0 license.
