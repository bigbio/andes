# the reference glyco engine / FragPipe (N-glyco) — clean-room study for andes

*Sources:* the reference glyco engine, Polasky et al., **Nat Methods 17, 1125 (2020)** —
[PMC7606558](https://pmc.ncbi.nlm.nih.gov/articles/PMC7606558/) /
[doi:10.1038/s41592-020-0967-9](https://www.nature.com/articles/s41592-020-0967-9).
Glycan FDR = an open-search PTM tool, Polasky et al., **MCP 21, 100205 (2022)** —
[PMC8933705](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705/) /
[doi:10.1016/j.mcpro.2022.100205](https://www.mcponline.org/article/S1535-9476(22)00013-5/fulltext).
Repo [github.com/Nesvilab/the reference engine](https://github.com/Nesvilab/the reference engine);
[FragPipe glyco tutorial](https://fragpipe.nesvilab.org/docs/tutorial_glyco.html).

> **License — AVOID copying code.** the reference engine, IonQuant, diaTracer are **UM-proprietary**
> (academic-only, commercial via Fragmatics/UM OTT). an open-search PTM tool itself is
> open-source, but the glyco pipeline runs *inside* proprietary the reference engine.
> **Clean-room rule holds: borrow published formulas only, never source.**

## 1. Candidate generation — **peptide-first, mass-offset (NOT glycan-Y-first)**
the reference glyco engine extends the peptide-first **fragment-index** engine. The glycan is a
**restricted mass-offset list** (labile / "mass offset" search), not an open search and
not a variable mod. Workflow:
1. Enumerate peptides; keep those with an N-X-S/T sequon (X≠P) for N-glyco.
2. For each such peptide, add to the fragment index: unmodified b/y, the **Y-ion ladder**
   (peptide backbone + retained glycan stubs: peptide+0, +HexNAc, +2HexNAc, +HexNAc·Hex,
   … up to intact), and optionally **b/y+HexNAc**. No b/y carrying the *intact* glycan.
3. Precursor delta = `precursor_neutral − peptide_mass`; must match one supplied glycan
   offset (they use the **~182-mass N-glycan offset list** from Riley et al.).
Composition is *inferred from the delta mass* + Y/oxonium evidence — the search itself is
composition-agnostic (one offset = one mass, many isomeric compositions collapse).

## 2. Glycan handling
- **Oxonium gate:** spectrum accepted as glyco only if summed oxonium intensity (204.087
  HexNAc, 366.14 HexHexNAc, 138.055, 163.06 Hex, 292/274 NeuAc…) ≥ **10 % of base peak**.
- **Y-ion ladder** = the discriminating backbone evidence: peptide+{0, HexNAc, 2HexNAc,
  2HexNAc+Hex (core), …}. These give the *peptide* mass anchor even when b/y are dead —
  exactly andes's Y0/Y1 anchor idea.
- **DB**: a fixed internal composition list (mammalian N-glycans; Hex/HexNAc/Fuc/NeuAc/
  NeuGc/etc.), **composition not structure**. a glyco search engine uses a structure DB — clean-room OK.

## 3. Scoring
- Single **hyperscore** (∝ Σ matched-fragment count × intensity), Y and b/y+HexNAc ions
  scored *equally* with normal b/y — no separate glycan sub-score in the primary match.
- **No learned/spectral-predictor** in the core search (contrast: andes's own strong model).
- Glycan composition scoring is deferred to **an open-search PTM tool** (§4), not the search.

## 4. FDR — **sequential two-stage, NOT joint 2D** (this is the key andes lesson)
Stage A: the reference engine search → **peptide** FDR via PeptideProphet's *extended mass model*
(each mass-offset modeled independently) + Philosopher, 1 % PSM/protein.
Stage B (an open-search PTM tool): on peptide-FDR-passing PSMs, assign the **glycan** with its own
target/decoy FDR. Reported as "1 % peptide **and** glycan FDR" — two thin filters in
series, **not a single joint likelihood**. Directly validates andes's §2-refuted-unified-pile
finding: keep the glycan axis a *separate post-process*.

**Glycan decoy recipe** (per target, seeded/reproducible):
shift intact target mass by a random Δ within the mass tolerance, assign a random isotope
error ∈{−1,0,+1,+2,+3}, and **randomize each theoretical Y/oxonium mass by ±1–20 Da**
(keeps composition & ion *count*, destroys mass coincidence). ≈1 decoy per target.

**Candidate pick (pairwise, Eq. 1):** `S = S_Y + S_oxo + β·log|Δm₂/Δm₁| + α_iso`, β=1.
- `S_Y = Σ_t 2·[α_hit,t·(U₁−U₂) + α_miss,t·(V₁−V₂)]`, t∈{HexNAc, Hex-only, Fuc}; Y counts
  **√-normalized** (don't reward big glycans); U=unique-found, V=unique-missing.
- `S_oxo = Σ_t 5·[α_hit,t·(I_obs/I_exp)(U₁−U₂) + α_miss,t·(V₁−V₂)]`, t∈{NeuAc/NeuGc,Fuc,
  Phos,Sulf}, intensity-weighted (capped), **not** √-normalized. Log prob-ratios (Table S1:
  NeuAc 1.5, Fuc 1.3, Phos 2.0).
Best candidate is then **absolute-scored** vs σ of unmodified-peptide mass error;
**glycan q = #decoy≥t / #target≥t**, threshold at 1 %.

## 5. Glycan notation
a commercial glyco engine-style condensed composition, e.g. `HexNAc2Hex7`, `Hex5HexNAc4Fuc1NeuAc2` — count
per monosaccharide, no structure/linkage. andes should emit the same for interop.

## 6. Single most valuable clean-room idea — **and what to AVOID**
**BORROW: the an open-search PTM tool separate-axis glycan FDR** — a per-PSM **absolute glycan score**
from Y+oxonium evidence with mass/isotope penalties, its own target/decoy list
(random-mass-shift + ±1–20 Da fragment randomization decoys), q-valued *after* peptide FDR.
This is exactly the fix for andes's crash (unified Percolator pile → 4.4 %): make it a
**thin Percolator post-process on the glycan axis only**. Reuse the four score terms
(S_Y √-normalized, intensity-weighted S_oxo, `β·log(Δm)`, isotope prior) as andes glycan
features, and the decoy recipe verbatim (it's published, not code).

**AVOID:** (a) the reference engine's peptide-first offset architecture — andes is deliberately
**glycan-Y-first + own learned models**; don't regress to composition-agnostic offsets.
(b) Folding glycan evidence into the *peptide* hyperscore/one FDR — the spectrum-level
Y/oxonium features can't separate two peptides at one backbone mass (andes SPA2 finding);
keep peptide-axis discrimination in the learned strong model + Y0/Y1 anchor, glycan-axis in
the separate FDR. (c) Any code reuse (proprietary).
