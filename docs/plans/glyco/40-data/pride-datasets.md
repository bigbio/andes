# PRIDE / ProteomeXchange N-glycopeptide harvest — annotated-result datasets

*Compiled 2026-07-02 for the andes glyco SP-B / 2D-FDR campaign. Sources: PRIDE
Archive REST (v2/v3) via MCP, manuscripts (PMC/DOI), tool repos. andes fails at
**identification** (ranking + FDR), not backbone generation (see
`00-context/00-current-state.md`) — so we harvest datasets that ship
**engine-annotated glycopeptide result tables** as external truth/training labels.*

## Selection principle

We need (spectra + a downloadable per-PSM glycopeptide result table). Raw-only
datasets are training-useless without re-search. **PXD025455 is the eval holdout**
(the 523-scan truth already in `tests/data/`); it never enters TRAIN. To fight
species/matrix leakage, TRAIN mixes human + mouse + serpent and multiple matrices
(brain / saliva / serum / IgG / venom), EVAL is human serum only.

## Confirmed candidates (all PRIDE-verified)

| PXD | Content | Instrument | Activation | Species / matrix | Result files? | Engine | Size |
|---|---|---|---|---|---|---|---|
| **PXD005411** | Mouse brain glycoproteome (a glyco search engine2 landmark) | LTQ Orbitrap | Stepped-energy HCD | *M. musculus*, brain | **Yes** — 5× `MouseBrain-Z-T-*-FDR.txt` (a glyco search engine2 per-GPSM tables) | a glyco search engine 2.0 | ~5 raw + 5 txt |
| **PXD030670** | Saliva N-glyco, lung cancer | Q Exactive | HCD | *H. sapiens*, saliva | **Yes** — `HILIC-Intact_glycopeptides*.xlsx` (a commercial glyco engine) | a commercial glyco engine 3.10 / Byologic | ~30 raw + xlsx |
| **PXD011239** | Serum haptoglobin HCC vs cirrhosis | Q Exactive | **EThcD** (not HCD) | *H. sapiens*, serum (Hp only) | **Partial** — `*.prot.xml` (a commercial glyco engine/PD, protein-level) | a commercial glyco engine / PD 2.1 | ~40 raw + prot.xml |
| **PXD020254** | 2D/3D breast cancer + xenograft | Orbitrap Fusion Lumos | **Stepped HCD (±10%)** | *H. sapiens*, cell/tumor | **Yes** (archived) — `HepG2.rar`, `Breastcancercelllines.rar` | (multi-enzyme HILIC) | ~68 raw + 2 rar |
| **PXD016175** | Plasma IgG, PCa vs BPH | Orbitrap Fusion Lumos | HCD | *H. sapiens*, plasma IgG | **Yes** — `Results.zip` (a glyco search engine 2.0) | a glyco search engine 2.0 | 96 raw + zip |
| **PXD057219** | *Bothrops* venom N-glyco | Q Exactive HF | HCD | *Bothrops* (snake) | **No** — only `.raw` + `Rawfilesdescription.xlsx` + `checksum.txt`; **GlycReSoft results NOT deposited** | GlycReSoft / a quantitation tool | 190 raw, no glyco table |
| PXD025455 *(HOLDOUT)* | NASH-HCC serum panel (the eval truth) | Q Exactive HF | **Stepped HCD** | *H. sapiens*, serum | Yes — `.pepXML` (a commercial glyco engine) | a commercial glyco engine / Byologic | large, DDA+PRM |

Instrument/activation/engine fields are from each project's `dataProcessingProtocol`;
result-file presence is from the verified PRIDE file listing (`get_project_files`).

## Notes that change the ranking

- **PXD005411** is the highest-value TRAIN anchor: it is the a glyco search engine2 reference set
  (Liu 2017, *Nat Commun* 8:438, PMC5585273, doi:10.1038/s41467-017-00535-2), ships
  clean per-GPSM `-FDR.txt` tables with glycan+peptide+glycopeptide FDR already
  controlled, and adds **mouse** (anti-leakage vs human eval). Stepped-energy CID/HCD
  on ion-trap-Orbitrap — sparse-b/y regime that stresses SP-B exactly like eval.
- **PXD030670** and **PXD016175**/PXD020254 give **a commercial glyco engine-labelled** and
  **a glyco search engine2-labelled** human sets on *different instruments* (QE vs Lumos) — good
  cross-instrument diversity; xlsx/zip/rar tables need a light parse step.
- **PXD011239** is EThcD, not HCD — different fragmentation regime; and its
  `.prot.xml` is protein-level, so per-GPSM glyco labels are weak. **Demote to
  optional / regime-probe only.**
- **PXD057219** (venom) is attractive for species diversity but **has no downloadable
  glyco result table** (GlycReSoft output absent) — it becomes a re-search-only
  candidate, not a labelled-truth source. **Drop from labelled TRAIN.**

## Clean-room license map (respect: papers OK, code NOT)

andes reference algorithms must come from **published papers only**. License status
of the reference engines:

- **a glyco search engine 2.0 / 3.0** — binary, **registration-gated proprietary** license
  (`i.pfind.net/license`); GitHub `pFindStudio/a glyco search engine` hosts a compiled release, not
  reusable source. → Paper algorithms (glycan-first, Y-ion + peptide 2-stage FDR)
  are clean-room references; **do not vendor their code**.
- **a cross-spectrum glyco engine** — **Apache-2.0** (`DICP-1809/a cross-spectrum glyco engine`), but distributed
  as `.exe`. Apache = permissive; its **cross-spectrum peptide-matching + glycome
  network-smoothing** (Fang 2022, PMC8990002; doi:10.1038/s41467-022-29530-y) is the
  clean-room reference for andes G4 (RT-gated cross-spectrum transfer). ~11% direct
  b/y backbone-ID figure cited in our current-state doc comes from here.
- **O-Pair Search / an open-source glyco engine** — **GPL-3.0** (`smith-chem-wisc/an open-source glyco engine`).
  Copyleft: **do NOT copy code** into Apache andes. The O-Pair *graph localization +
  paired-scan* method (Lu 2020, *Nat Methods* 17:1133, PMC7606753) is paper-citable
  only. (O-glyco-focused; lower priority than N-glyco refs.)
- **a commercial glyco engine / Byologic / the reference glyco engine** — commercial / UM-proprietary. **Labels
  only** (their result TSVs are usable as external truth); **never** algorithm/code.

**2D-FDR reminder:** andes FDR is **Percolator-only**; the a glyco search engine separate glycan/
peptide-axis FDR (Liu 2017) is the *conceptual* reference for a thin Percolator
post-process (`G3′`), not a second FDR engine, and never Mokapot.

## Recommended TRAIN / EVAL split (anti-leakage)

**EVAL (frozen):** PXD025455 only — human serum, QE-HF stepped-HCD, a commercial glyco engine truth.
Never used for training or model selection.

**TRAIN (labelled truth, ranked):**
1. **PXD005411** (mouse brain, a glyco search engine2 `-FDR.txt`) — primary; different species,
   clean per-GPSM FDR labels, stepped-energy regime.
2. **PXD016175** (human plasma IgG, a glyco search engine2 `Results.zip`) — different instrument
   (Lumos), well-characterised IgG glycoforms; a glyco search engine2 labels homogeneous with #1.
3. **PXD030670** (human saliva, a commercial glyco engine xlsx) — a commercial glyco engine label diversity, QE.
4. **PXD020254** (human cell/tumor, Lumos stepped-HCD, `.rar`) — matrix diversity;
   parse archive to confirm per-GPSM table before trusting.

**Species/matrix mix achieved:** mouse-brain + human-plasma-IgG + human-saliva +
human-cell/tumor across LTQ-Orbitrap / Lumos / QE — no single organism or matrix
dominates, and **no serum stepped-HCD-QE-HF** leaks from EVAL into TRAIN (matrix +
instrument + activation of PXD025455 are held out as a combination).

**Optional / probe:** PXD011239 (EThcD regime check, weak labels); PXD057219
(species diversity only if re-searched — no deposited labels).

## Sources

- a glyco search engine2: Liu et al. *Nat Commun* 2017, PMC5585273, doi:10.1038/s41467-017-00535-2
- a cross-spectrum glyco engine: Fang et al. *Nat Commun* 2022, PMC8990002,
  doi:10.1038/s41467-022-29530-y; repo `github.com/DICP-1809/a cross-spectrum glyco engine` (Apache-2.0)
- a glyco search engine: Zeng et al. *Nat Methods* 2021, PMC8648562, doi:10.1038/s41592-021-01306-0;
  repo `github.com/pFindStudio/a glyco search engine` (registration-gated)
- O-Pair: Lu et al. *Nat Methods* 2020, PMC7606753, doi:10.1038/s41592-020-00985-5;
  `github.com/smith-chem-wisc/an open-source glyco engine` (GPL-3.0)
- Dataset detail: PRIDE Archive `get_project_details` / `get_project_files`
  (accessions above); OmicsDI mirror `omicsdi.org/dataset/pride/PXD005411`
