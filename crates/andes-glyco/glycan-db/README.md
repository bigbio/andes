# Bundled pGlyco glycan databases

These `.gdb` files are the pGlyco N-glycan structure databases
([Liu et al., *Nat. Commun.* 2017](https://www.nature.com/articles/ncomms15473)),
redistributed here so `--glyco-species` can select a species-specific glycan space
without a local pGlyco installation. Each file is a header line of monosaccharide
symbols followed by one canonical S-expression per glycan (see
`andes-glyco::glycan_db::load_glycan_gdb` for the grammar).

| File | Species / scope | Structures |
|---|---|---|
| `pGlyco-N-HighMannose.gdb` | high-mannose N-glycans (species-independent) | 30 |
| `pGlyco-N-Human.gdb` | *Homo sapiens* N-glycans | 2922 |
| `pGlyco-N-Human-multi.gdb` | human, multi-antennary | 5634 |
| `pGlyco-N-Mouse.gdb` | *Mus musculus* N-glycans | 6662 |
| `pGlyco-N-Mouse-large.gdb` | mouse, extended | 7878 |

Plant N-glycan databases are not bundled: plant glycans are xylosylated (`X` =
Xyl), which the search model (`GlycanComp`) does not yet represent.

Attribution: glycan structures sourced from the pGlyco project
(https://github.com/pFindStudio/pGlyco3). If you redistribute andes with these
files, keep this notice.
