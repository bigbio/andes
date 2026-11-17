# Simulated substitutions: localization accuracy

This is the benchmark design of Yu et al., *Nat. Commun.* 2020, "Identification of
modified peptides using localization-aware open search": a **substituted** database is searched
with open search against **unmodified** spectra, so the mass offset is a known substitution's
mass difference and the substituted residue is the ground-truth site. It measures the two things
the paper is about — whether the modified peptide is found, and whether the offset is placed on
the right residue — without depending on any declaration of the modification.

## Data and method

| | |
|---|---|
| Spectra | PXD001468, `b1906_293T_proteinID_01A_QE3_122212.mzML` (HEK293T, Q Exactive), 41,820 MS2 |
| Plain database | `UP000005640_9606.fasta` (20,610 human entries) |
| Substituted database | `UP000005640_9606_modified.fasta`, same entries with **58,015 random substitutions over 8,832 proteins** (offset range −129..+129 Da) |
| andes | 0.4.0, branch `feat/open-search` (`7e6d71d4` … `c6953ead` and this commit) |
| Percolator | 3.7.1, `--seed 42 -Y --only-psms=false` |

Two searches of the same spectra:

- **reference** — narrow search against the plain database with the paper's variable mods
  (andes defaults: carbamidomethyl-C fixed; oxidation-M and protein-N-term-acetyl variable). Its
  PSMs at `q ≤ 0.01` are the truth peptides.
- **open** — `--open-search --open-search-delta=-150..150` against the substituted database,
  fixed mods only. `--precursor-mono auto` is active, and the localization-aware re-score runs
  on every offset candidate.

For every reference PSM whose span carries **exactly one** substitution, the truth is that
substituted peptide and the offset `reference_calc_mass − substituted_span_mass`
(original − substituted, so a mutation that makes the peptide lighter gives a positive offset).
Peptides with variable mods are excluded (the simulation is on unmodified spectra).

## Results

| metric | value |
|---|---:|
| reference PSMs at `q ≤ 0.01` | 23,076 |
| single-substitution truth cases | **8,563** |
| open-search PSMs at `q ≤ 0.01` | 18,645 |
| same scan + substituted peptide found | 3,726 (43.5% of truth) |
| … with the offset within 0.02 Da | **3,011 (35.2%)** |
| … whose site was called | 1,772 (58.8% of those) |
| … at the substituted residue | **1,595 (90.0% of called)** |
| median absolute offset error | **0.0008 Da** |

## Shifted-ion shortlist (A/B)

The shortlist votes on each observed peak and on its complement
(`b + y = peptide + 2 H+`, so a shifted fragment's complement is the unmodified peptide's own
ion). Measured during development with a regular-vote arm that a then-present development flag
selected (the flag has been removed; open search now always votes with the shifted ions):

| metric | off | on (default) |
|---|---:|---:|
| open-search PSMs at `q ≤ 0.01` | 18,403 | 18,645 |
| same scan + substituted peptide | 3,235 (37.8%) | **3,726 (43.5%)** |
| offset within 0.02 Da | 2,516 (29.4%) | **3,011 (35.2%)** |
| site correct / called | 1,395 / 1,531 (91.1%) | 1,595 / 1,772 (90.0%) |
| median absolute offset error | 0.00089 Da | 0.00078 Da |
| open-search wall | 120.5 s | 163.6 s |

The shifted vote recovers 20% more of the truth peptides at the correct offset; it costs 36% of
the open search's wall time and moves localization accuracy by one point (the extra cases are
slightly harder).

## Comparison with the published benchmark (Yu et al. 2020)

The paper's simulated data is the same HeLa/HEK293T set used here: 24 fractions of PXD001468,
spectra unmodified, half of the high-confidence peptides given a single random amino-acid
replacement, searched against the **modified** database. Published totals (1% spectrum FDR,
Philosopher for MSFragger, each tool's own FDR module for the others):

| tool (version) | substitution-containing PSMs | substitution-free PSMs | precision (Type 1+2) |
|---|---:|---:|---:|
| MSFragger 2.2, regular open search | 81,648 | 178,966 | 0.9294 |
| **MSFragger 2.2, localization-aware OS** | **134,413 (+65%)** | 178,958 | **0.9623** |
| pFind3 3.1.5 (with mixture spectra) | ~106,000 (~134,000) | — | — |
| MetaMorpheus 0.0.303 (with precursor deconvolution) | ~102,000 (~139,000) | — | — |
| TagGraph 1.8 | ~113,000 | — | — |

The non-MSFragger counts are read from Fig. 3a and are approximate; only the MSFragger rows are
stated numerically in the text. On the phosphorylation set (PXD004940) the paper reports 23,828
PSMs for LOS against 20,072 for regular OS (+19%).

**Same file, same spectra.** One of the paper's 24 files is `b1906_293T_proteinID_01A`, the file
this benchmark runs. Its Supplementary Data 4 lists MSFragger LOS's PSMs per scan, so the
comparison can be made on that file alone (though the two runs use different substituted
databases — the paper's replacements are confined to its high-confidence peptides, ours are
58,015 random substitutions across the proteome — and different FDR tooling):

| run, b1906 only | substitution-containing PSMs | Type 1+2 share of them |
|---|---:|---:|
| MSFragger LOS (paper, Supplementary Data 4) | **5,320** | 0.9623 (paper, over all PSMs) |
| andes, regular votes | 4,737 | 0.897 |
| andes, shifted-ion votes | **5,185** | 0.904 |

On this file andes' shifted-ion open search finds 97.5% as many substitution-containing PSMs as
MSFragger LOS (−2.5%) and 11% more than the same search with regular votes. The Type 1+2 share
is computed over containing PSMs only (the paper's is over all PSMs), so the two are not the same
statistic; the paper's Type 2 is judged against its Comet ∩ MS-GF+ list, ours against the andes
narrow search.

Parameter optimization: the paper reports an average +10% PSMs across four datasets from mass
calibration plus parameter optimization (its benchmark search used 20 ppm and 100 peaks); andes
`--optimize` on one file gives +5.1% (23,076 → 24,244 PSMs). The two procedures optimize
different knobs (the paper searches fragment tolerance and peak count after a 2-D m/z/RT
calibration; andes searches precursor tolerance, top-N and charge range, and measures the
fragment window from the run).

## Summary of configurations (all on this dataset)

| reference search | open search | reference PSMs | truth cases | open PSMs | same peptide | offset within 0.02 Da | site correct / called | median offset error |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| default 20 ppm | regular votes | 23,076 | 8,563 | 18,403 | 3,235 (37.8%) | 2,516 (29.4%) | 1,395 / 1,531 (91.1%) | 0.00089 Da |
| default 20 ppm | **shifted votes** | 23,076 | 8,563 | 18,645 | 3,726 (43.5%) | **3,011 (35.2%)** | 1,595 / 1,772 (90.0%) | 0.00078 Da |
| **`--optimize` (50 ppm)** | shifted votes | 24,244 | 8,872 | 18,645 | 3,773 (42.5%) | 3,025 (34.1%) | 1,597 / 1,778 (89.8%) | 0.00078 Da |

The rates use their own truth denominator (8,563 vs 8,872): the optimized reference search
finds 5.1% more identifications, and the 309 extra single-substitution cases are harder than
average, so the rate is slightly lower while the recovered count is higher (3,025 vs 3,011).
The open search is the same run in the last two rows; only the reference (truth) search differs.

Reference-search settings, same file:

| reference setting | reference PSMs at 1% | truth cases | note |
|---|---:|---:|---|
| default 20 ppm | 23,076 | 8,563 | |
| `--precursor-cal auto` | 21,082 | — | tightens 20 ppm to 4.24 ppm, **−8.6%** |
| `--optimize` | 24,244 | 8,872 | picks 50 ppm on a strided 5,000-spectrum sample, **+5.1%** |

## Parameter search on the reference search

`--optimize` (a strided 5,000-spectrum sample) chose **50 ppm** for the reference search, where
the tuned default is 20 ppm, and the reference's Percolator run then had **24,244 PSMs at 1%**
instead of 23,076 (+5.1%). With that larger truth set (8,872 single-substitution cases) the open
search recovers 3,025 (**34.1%**) with the offset within 0.02 Da, and 1,597 of 1,778 called sites
(**89.8%**) are the substituted residue.

The gain is not automatic: a sample of the *first* 1,000 spectra preferred 5 ppm, which lost 7%
of the full file's PSMs. The grid needs a sample spread across the file and the native
rescorer's features to rank settings; both are what `--optimize` does (see `DOCS.md`).

The open search's own PSMs are reported at the same 1% (target-decoy) threshold as the
reference; no entrapment component exists in the human database, so this is a nominal rather
than a measured FDR (contrast the mouse-liver benchmark, where entrapment puts the offset rows
at 0.72%).

**Reading.** Roughly one in three reference peptides is recovered as its substituted form with
the correct mass, and when the fragmentation pins the offset to one residue, that residue is the
substituted one **90% of the time** at a median offset error of 0.8 mDa. The recovery is lower
than the paper's MSFragger results because andes scores the substituted peptide directly instead
of indexing shifted fragment ions and jointly scoring them; the localization accuracy
conditional on a called site is the part this implementation reproduces.

On calibration: andes' **fragment**-error calibration is automatic and was active (it sets the
feature tolerance and node-score window from confident matches), and the open search also uses
the MS1 monoisotopic correction. Its **precursor** calibration (`--precursor-cal auto`) is a
file-wide ppm shift followed by a tightened window; it is off by default and was not used here —
enabling it on the reference search tightens 20 ppm to 4.24 ppm, costs 16 s, and **loses 8.6% of
the PSMs** (23,076 → 21,082), matching its documented behaviour on the UPS1 benchmark (−7.6%).
The paper's supervised m/z- and RT-dependent mass calibration and its parameter-grid
optimization have no andes equivalent.

The 2,516 recovered offsets are 29.4% of the 8,563 truth cases but only 14% of the 18,403 open
PSMs: most open-search PSMs are background matches, which is expected when the whole proteome is
eligible at every offset.

## Reproduce

```bash
cargo build --release -p andes --bin andes --features thermo
export DOTNET_ROOT=~/.dotnet

./run_substitution_bench.sh /path/to/pxd001468-parent
# one file by default; passes OUT=<dir>, THREADS=<n> and ANDES=<binary> through.
```

`run_substitution_bench.sh` runs both searches, both Percolator runs and
`analyze_substitution.py`, which writes `report.json` with the counts above. The scripts take no
hardcoded host paths.

## Limits

- One file of the four datasets in the paper; the other three are not on this host.
- The substituted database is the one already present in the PXD001468 download; its
  substitution rate (~1.3% of residues) is the design used here, not a re-run of the paper's
  generator.
- Variable-mod reference peptides are excluded, so 3,594 of the 23,076 reference PSMs are not
  scored; 10,919 more are excluded because their span has zero or several substitutions.
- The open search's FDR is nominal (target-decoy), not entrapment-measured.
