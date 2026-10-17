# Bundled glycan databases

N-glycan structure databases (see Attribution), bundled so `--glyco-species` works without an
external glycan database. Each file is a header line of monosaccharide symbols, then one
canonical S-expression per glycan (grammar: `andes-glyco::glycan_db::load_glycan_gdb`).

| File | Species / scope | Structures | Provenance |
|---|---|---|---|
| `N-glycan-high-mannose.gdb` | high-mannose N-glycans (species-independent) | 30 | derived (curated subset) |
| `N-glycan-human.gdb` | *Homo sapiens* N-glycans | 2922 | verbatim (renamed) |
| `N-glycan-mouse.gdb` | *Mus musculus* N-glycans | 6662 | verbatim (renamed) |
| `N-glycan-mouse-large.gdb` | mouse, extended | 7878 | verbatim (renamed) |

The upstream `*-multi` databases are not bundled: each of their extra lines is a combination
of glycans for a peptide with several glycosites, and andes places one glycan per peptide.
Plant databases are not bundled: xylose (`X`) is not yet represented in `GlycanComp`.

Attribution: redistributed from the pGlyco project
([Liu et al., *Nat. Commun.* 2017](https://www.nature.com/articles/ncomms15473);
https://github.com/pFindStudio/pGlyco3), Apache-2.0. If you redistribute andes with these
files, keep this notice and the upstream Apache-2.0 license.
