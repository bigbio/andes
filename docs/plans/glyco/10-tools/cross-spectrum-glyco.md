# a cross-spectrum glyco engine — clean-room reference for andes N-glyco

**Sources.** Paper: Yang *et al.*, "a cross-spectrum glyco engine enables glycan database-independent
peptide matching and in-depth characterization of site-specific N-glycosylation,"
*Nat Commun* **13**, 1900 (2022). DOI 10.1038/s41467-022-29530-y · PMC8990002
([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/)). Repo:
[github.com/DICP-1809/a cross-spectrum glyco engine](https://github.com/DICP-1809/a cross-spectrum glyco engine)
(Java + JavaFX GUI). **License: Apache-2.0** — clean-room-borrowable (algorithms from
the published paper; do not copy Java verbatim, but the recipes below are paper-derived).
Lab: Mingliang Ye, DICP/CAS. Headline claim: **+33.5%–178.5%** more glyco-PSMs than
a commercial glyco engine / the reference glyco engine / StrucGP / a glyco search engine 3.0.

## 1. Candidate generation — PEPTIDE-first, glycan-database-INDEPENDENT

The core inversion vs andes: a cross-spectrum glyco engine finds the **peptide backbone first, without
committing to any glycan database**.

1. **In-silico deglycosylation.** Strip glycan B/Y ions from the MS2 by removing peaks at
   the characteristic N-glycan core mass gaps (see §2 ladder). The residual spectrum is
   treated as a bare-peptide spectrum.
2. **Peptide DB search with MS-GF+** against sequon-containing tryptic peptides
   (N-X-S/T/C, X≠P). The precursor is set to a de-glycosylated backbone mass derived from
   the Y-ladder anchor. **They reuse MS-GF+ as the backbone scorer** — directly relevant
   to andes (own MS-GF+-lineage scorer).
3. **Glycan = precursor − backbone**, then annotated by composition (not looked up in a
   fixed DB), with **monosaccharide stepping** to explain unusual/modified glycans:
   iteratively add one monosaccharide (or a modified unit, e.g. +179 Da ammonium) and keep
   stepping while new B/Y ions match, "from Y-Hex(3)HexNAc(2) to the terminal."

There is **no fragment-index / open-mass-offset engine** (contrast the reference engine). Generation
is: Y-ladder anchor → backbone mass → MS-GF+ → glycan by mass difference + stepping.

## 2. Glycan handling

- **DB source:** GlyTouCan (≈1,766 unique compositions / 10,936 entries), **WURCS 2.0**
  encoding — but used only for *annotation*, not to constrain the search (the whole point).
- **Composition, not structure** for scoring; stepping can characterise structure-ish
  branch order via the sequence of matched Y ions.
- **Core-Y ladder** (the anchor set, N-glycan invariant core):
  `Y0, Y+HexNAc(1), Y+HexNAc(2), Y+Hex(1)HexNAc(2), Y+Hex(2)HexNAc(2), Y+Hex(3)HexNAc(2)`.
  Masses (backbone + adduct): HexNAc 203.0794, Hex 162.0528, Fuc 146.0579,
  NeuAc 291.0954, NeuGc 307.0903.
- **Oxonium ions** gate glyco-spectra; the Y1 (peptide+HexNAc) ion is *required* for a
  confident expansion match.

## 3. Scoring

Peptide and glycan scores are **separate additive terms**, both intensity-weighted by an
empirical fragment **frequency** raised to α:

```
Score_peptide = coeff · Σ_i  Intensity_i · frequency_i^α        (α = 0.3)
Score_core    = Σ_j  Intensity_j · frequency_j                  (core-Y ions, m terms)
Score_PSM     = Score_peptide + Score_core
Score_glycan  = Σ_i  Intensity_i · frequency_i^α
```

- `frequency_i` = empirical occurrence of fragment *i* across PSMs of that backbone
  (a **data-derived fragment prior**, not a deep predictor).
- `coeff` = **cosine similarity** between this spectrum's backbone-fragment pattern and the
  **average pattern for that peptide backbone** — the cross-spectrum term.
- **Spectrum expansion (the +48% idea):** backbone fragmentation is glycan- and
  charge-independent, so the fragment pattern learned from *confidently identified* PSMs of
  a backbone is **transferred** to score *unassigned* spectra of the same backbone that
  MS-GF+ alone could not identify. This is a **cross-spectrum prior**, not a per-spectrum
  score. Only ~11% of glyco-spectra yield a direct backbone ID; expansion recovers the rest.

## 4. FDR — 2D, decoy = random fragment m/z shift

- **Decoy construction (glycan axis):** for each GPSM, a decoy spectrum is built by
  **shifting each MS2 fragment m/z by a random 1–30 Da**. This is a *spectrum* decoy, not a
  sequence-reversal decoy — cheap and glycan-axis-specific.
- **Peptide axis:** standard MS-GF+ target/decoy on the deglycosylated search; keep q<0.01.
- **FDR estimator:** `FDR = 2·N_decoy / (N_decoy + N_target)`, PSMs ranked by expectation
  value via **linear tail-fit**.
- **2D gate = intersection of two independent axes:** peptide-axis q<0.01 **AND**
  glycan-axis FDR<0.01 **AND** ≥3 matched core-Y ions with **Y1 required**. The two axes are
  computed *separately* then combined — never one merged decoy pile.

## 5. Notation emitted

Composition strings: `Hex(9)HexNAc(2)`, `Hex(7)HexNAc(2) + 179 (NH3 adduct)`, with
Fuc/NeuAc/NeuGc and phospho/ammonium modifier annotations; structures in WURCS 2.0.

## 6. Single most valuable clean-room idea for andes — and what to avoid

**BORROW: the separate two-axis FDR with a random-m/z-shift glycan decoy.** andes already
learned (SPA2 / current-state §2) that a *unified* Percolator pile with glycan-decoy rows
crashed recovery 29.4%→4.4%. a cross-spectrum glyco engine's answer is exactly the fix andes flagged as
**G3′**: compute the peptide axis in Percolator (targets vs sequence decoys — the andes
hard-constraint: **Percolator only, never Mokapot**) and the **glycan axis as a thin
post-process** using an independent random-shift decoy + `FDR = 2·N_d/(N_d+N_t)`, then gate
on the **intersection** plus the ≥3-core-Y / Y1-required rule. This is a clean-room,
Percolator-compatible 2D-FDR — no code copied, recipe from the paper.

**Second borrow (aligned with andes G4):** cross-spectrum **backbone-pattern transfer** —
learn each backbone's fragment prior from confident PSMs, transfer via **cosine coeff** to
rank sparse-b/y spectra. This is the +48% lever and precisely what andes's sparse-stepped-HCD
stratum needs (single-spectrum scoring cannot rank it). andes differentiates by making the
prior a **learned regime-matched model + RT gating**, not a flat empirical frequency table.

**AVOID / differentiate:**
- **Do not adopt peptide-first as andes's primary generation.** andes is glycan-Y-first
  (verified G1: +7–10 pts backbone-findability) and near-ceiling (~90%). Peptide-first
  reintroduces the deglycosylation-anchor fragility a cross-spectrum glyco engine depends on; keep it only
  as a *fallback* branch, not the spine.
- **Do not copy the flat `frequency^0.3` empirical prior** as the scorer — andes's edge is a
  *learned* peptide-axis model (SP-B), which should subsume it.
- **Do not merge glycan-decoy rows into the Percolator target/decoy pile** (already refuted).
- MS-GF+ reuse is validated by them, but andes's own scorer + Percolator is the productionised
  path; treat a cross-spectrum glyco engine as algorithmic confirmation, not an implementation to mirror.
