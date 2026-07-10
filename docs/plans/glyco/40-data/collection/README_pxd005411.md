# Glyco-PSM example collection — PXD005411 (a glyco search engine2 mouse brain)

**Status: REAL-HARVESTED** from a public PRIDE identification result file (not illustrative).

## Provenance

| Field | Value |
|---|---|
| Dataset | **PXD005411** — "Mouse glycoproteome: Brain — a glyco search engine 2" |
| Publication | Liu MQ, Zeng WF, Fang P, *et al.* **a glyco search engine 2.0 enables precision N-glycoproteomics with comprehensive quality control and one-step mass spectrometry for intact glycopeptide identification.** *Nat Commun* **8**, 438 (2017). PMID:[28874712](https://pubmed.ncbi.nlm.nih.gov/28874712/) · DOI:[10.1038/s41467-017-00535-2](https://doi.org/10.1038/s41467-017-00535-2) |
| Organism / tissue | *Mus musculus*, brain |
| Instrument | LTQ Orbitrap, **stepped-energy HCD** (stepped-collision-energy, the same regime andes targets) |
| Search engine | a glyco search engine 2.0 (Apache-2.0 lineage; clean-room reference OK per campaign constraints) |
| Source file | `MouseBrain-Z-T-2-FDR.txt` (PRIDE `SEARCH` output, 5.0 MB, 17,188 PSM rows) |
| Download (public FTP) | `https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005411/MouseBrain-Z-T-2-FDR.txt` |
| Harvest date | 2026-07-02 |
| Selection | Rank==1, `GlyDecoy==0` & `PepDecoy==0`, `TotalFDR ≤ 0.01`; then a **balanced 45-PSM sample** stratified across glycan classes (12 high-mannose / 12 fucosylated / 12 sialylated / 9 complex-hybrid), deduplicated by scan, ranked within class by `TotalScore`. |

The a glyco search engine2 mouse-brain series is the canonical benchmark for stepped-HCD N-glyco
(same fragmentation regime as andes' PXD025455 truth set), so its glycan
compositions + backbone masses are an authoritative notation/mass fixture.

## Mass computation (standardized, this repo's convention)

Monoisotopic residue masses used (from the campaign standard; identical to the
values in the harvest instruction / `30-standards/masses.md`):

```
Hex     162.0528     HexNAc  203.0794     Fuc     146.0579
NeuAc   291.0954     NeuGc   307.0903     (residue/dehydrated masses)
proton  1.0072764    water   18.0105646
```

- `glycan_mono_mass_calc = Hex·162.0528 + HexNAc·203.0794 + Fuc·146.0579 + NeuAc·291.0954 + NeuGc·307.0903`
- `peptide_mono_mass_file = PeptideMH(file) − proton` (bare backbone, glycan stripped; from a glyco search engine's `PeptideMH`)
- `glycopeptide_neutral_mono_calc = peptide_mono_mass_file + glycan_mono_mass_calc`
- `precursor_mz_calc = (glycopeptide_neutral_mono_calc + z·proton) / z`

**Validation.** Our `glycopeptide_neutral + proton` reproduces a glyco search engine's
`PrecursorMH` to **1–5 mDa** across the sample (the residual is standard-mass
rounding vs a glyco search engine's internal element table, well under the 20 ppm search
window). This confirms the mass pipeline and the glycan-vector decoding below.

## Columns (`psms_pxd005411.tsv`, tab-separated)

| Column | Meaning |
|---|---|
| `dataset` | `PXD005411` |
| `source_file` | raw/mzML base name (`MouseBrain-Z-T-2`) |
| `scan`, `charge` | precursor scan number and charge (parsed from a glyco search engine `GlySpec` `File.scan.scan.z.rank.dta`) |
| `peptide_a glyco search engine` | a glyco search engine backbone string; **`J` = the glycosylated Asn** (glycosite marker) |
| `peptide_plain` | same sequence with `J`→`N` (plain 20-AA alphabet) |
| `glycosite_pos` | 1-based position of the glycosite within the peptide (a glyco search engine `GlySite`) |
| `modifications` | a glyco search engine `Mod` (e.g. `Carbamidomethyl[C]`; `null` if none) |
| `glycan_a glyco search engine_vec` | the raw 5-integer a glyco search engine vector (**order Hex HexNAc NeuAc NeuGc Fuc** — see gotcha #1) |
| `Hex`,`HexNAc`,`Fuc`,`NeuAc`,`NeuGc` | **canonical** monosaccharide counts (re-ordered to this repo's tuple) |
| `canonical_glycan` | `HexNAc(n)Hex(n)Fuc(n)NeuAc(n)NeuGc(n)` string |
| `glycan_mono_mass_calc` | glycan mass computed from standardized masses above |
| `glycan_mass_file` | a glyco search engine `GlyMass` (cross-check; matches calc to <0.02 Da for all 688 vectors in the file) |
| `peptide_mono_mass_file` | bare-backbone neutral monoisotopic mass (`PeptideMH − proton`) |
| `glycopeptide_neutral_mono_calc` | backbone + glycan neutral mass |
| `precursorMH_file` | a glyco search engine `PrecursorMH` (observed glycopeptide [M+H]) |
| `precursor_mz_calc` | our computed precursor m/z at `charge` |
| `total_score`,`pep_score`,`gly_score` | a glyco search engine `TotalScore`/`PepScore`/`GlyScore` |
| `ppm_file` | a glyco search engine precursor mass error (ppm) |
| `total_fdr` | a glyco search engine `TotalFDR` (glycopeptide-level; all rows ≤ 0.01) |
| `protein`,`prosite` | source protein(s) and site(s) (`/`-joined for shared peptides) |

## Notes / caveats

- a glyco search engine2 `TotalFDR` is a **separate-axis glyco FDR** (glycan-axis × peptide-axis,
  combined) — exactly the 2D-FDR structure the andes roadmap wants as a Percolator
  post-process (G3′). These rows are the *target* stratum only (decoys excluded).
- 45 PSMs span **17 distinct backbones** — glycoform multiplicity (one backbone,
  many glycans) is the normal N-glyco pattern and is preserved here.
- Sequon note: several backbones end `...J[K/R]` (glycosite is the penultimate
  residue). The N-X-S/T sequon then spans the tryptic cleavage (the S/T sits on
  the next peptide); a glyco search engine validated these against the full protein sequence, so
  they are correct glycosites even though `N-X-[ST]` is not visible within the
  peptide string. All non-boundary sequons satisfy N-X-[ST] (X≠P).
- Clean-room: a glyco search engine is an acceptable published/permissive reference. This is
  **result harvesting** (public identification output), not code reuse.
