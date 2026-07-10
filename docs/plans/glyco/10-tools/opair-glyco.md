# an open-source glyco engine / O-Pair — N-glyco reference for andes

*Sources: O-Pair, Nat Methods 2020 ([PMC7606753](https://pmc.ncbi.nlm.nih.gov/articles/PMC7606753/), [doi:10.1038/s41592-020-00985-5](https://doi.org/10.1038/s41592-020-00985-5)); the multi-attribute glycan-score + glycan-FDR model, Mol Cell Proteomics 2022 ([PMC8933705](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705/), [doi:10.1016/j.mcpro.2022.100205](https://doi.org/10.1016/j.mcpro.2022.100205)); repo [github.com/smith-chem-wisc/an open-source glyco engine](https://github.com/smith-chem-wisc/an open-source glyco engine). The 2022 scoring/FDR paper is the same Smith/Nesvizhskii-adjacent lineage shared by an open-search PTM tool; its equations are the published, clean-room-usable core.*

## 1. Candidate generation — PEPTIDE-FIRST, ion-indexed open search
O-Pair builds a **fragment-ion index**: all theoretical target+decoy peptide fragment masses from the whole protein DB are indexed once. Per MS2, it assembles every peptide with ≥1 fragment matching a peak, ranks by match count, and takes top candidates — an **open (mass-offset) search** where the glycan is a delta mass on the precursor, not enumerated up front. Rationale: HCD O-glycans are labile and "rarely retain glycan mass," so the backbone is found glycan-agnostic, then the delta is assigned. Reported ~2,000× faster than prior O-glyco tools.
> **Contrast with andes:** this is the the reference engine/a comparison search engine lineage (peptide-first + mass-offset). andes is deliberately **glycan-Y-first** (G1: precursor − known glycan → backbone; VERIFIED 59→70% findability). Do NOT clone peptide-first; borrow their *scoring/FDR*, keep our generation.

## 2. Glycan handling
Composition-only (no topology): an internal curated list of mammalian N-glycan compositions (≈disk DB; user-replaceable). Residue categories: HexNAc, Hex, Fuc, NeuAc, NeuGc, Phosphate, Sulfate. Uses **B/oxonium ions** (204.087 HexNAc marker, sialic markers 274/292, etc.) and a **Y-ion ladder** (peptide+peptide+HexNAc+…+intact) as the two orthogonal evidence axes.

## 3. Scoring — the multi-attribute glycan score (the valuable part)
Pairwise glycan score between two competing compositions (Eq. 1–3):
```
S_pairwise = S_Y + S_oxo + β_mass·log|Δm₂/Δm₁| + α_isotope
S_Y   = Σ_t ( α_Yt_hit·(U₁−U₂) + α_Yt_miss·(V₁−V₂) ) · √(norm)      # t ∈ {HexNAc, Hex-only, Fuc}
S_oxo = Σ_t ( α_ot_hit·(I_obs/I_exp)·(U₁−U₂) + α_ot_miss·(V₁−V₂) )   # t ∈ {NeuAc/NeuGc, Fuc, Phos, Sulf}
```
- **U** = # unique matched fragments of that type, **V** = # theoretical-but-missed. **α_hit = log(P_hit ratio) > 0, α_miss = log(miss ratio) < 0** — empirical per-type log-probability-ratios (values in their Table S1). This is a **learned/empirical likelihood-ratio**, not a raw count.
- **√-normalization on Y-ions only** — stops large glycans (more theoretical Y) from being over-rewarded. Oxonium is intensity-weighted (`I_obs/I_exp`) with a low-intensity floor (no negative for weak oxonium).
- **Mass term** `β_mass·log|Δm₂/Δm₁|` (β_mass=1.0 default), isotope removed first; **isotope term** α_isotope is a log-ratio over allowed errors {−1,0,+1,+2,+3}.
- **Absolute score** (Eq. 4–6) drops the pairwise Δ and uses `log|Δm/σ_unmod|` (σ = empirical unmodified-peptide mass error) — this is the score thresholded for FDR.
- Peptide b/y (Morpheus: #matched fragments + matched-intensity fraction) is scored/FDR'd **separately, first**; glycan score is a second, independent stage.

## 4. FDR — separate glycan-axis target-decoy (thin, two-stage; NOT joint 2D)
**Decoy-glycan recipe (borrowable, clean):** per target glycan build ONE decoy that keeps the *same nominal composition/# fragments* but (a) shifts intact glycan mass by a random value within tolerance, (b) shifts **each Y and oxonium fragment by a unique random 1–20 Da**, (c) assigns a random isotope error, (d) **fixed RNG seed** → deterministic. So a decoy can only win by chance, never by real mass agreement. FDR = cumulative-decoys/cumulative-targets on the **absolute glycan score**; q per PSM. Applied as **two stages**: peptide-FDR (1%) then glycan-FDR (1%) on the peptide-passing set — cumulative stringency, matching a glyco search engine's separate-axis philosophy, **not** one joint pile.
> This is exactly andes G3′: our unified-Percolator-pile crashed recovery 29→4%; the fix is this **separate-axis post-process**. For andes: keep Percolator for the peptide axis, add a thin glycan-axis target/decoy q on the random-shifted-fragment decoy — no Mokapot.

## 5. Notation
`H`=Hex, `N`=HexNAc, `F`=Fuc, `A`=NeuAc, `G`=NeuGc, `P`=Phos, `S`=Sulf; e.g. `HexNAc(4)Hex(5)NeuAc(1)` = N4H5A1.

## 6. License
**MIT** ([LICENSE.txt](https://github.com/smith-chem-wisc/an open-source glyco engine/blob/master/LICENSE.txt)) — permissive; algorithms and even code are clean-room OK to reference (contrast a commercial glyco engine commercial / the reference engine UM-proprietary — algorithms only from *papers*, which these are).

## 7. Single most valuable borrow — and what to AVOID
**BORROW:** the **empirical log-probability-ratio glycan score** (per-fragment-type α_hit/α_miss on Y + intensity-weighted oxonium + mass/isotope terms) *as a per-composition PIN feature*, plus their **random-fragment-shift decoy-glycan recipe** to drive a **separate-axis glycan q-value** (SP-B + G3′ in one move). Crucially this score is a **likelihood-ratio between competing (peptide,glycan) candidates at one backbone mass** — precisely the target/decoy-*discriminating* signal SPA2_RESULT found missing (spectrum-level oxonium/YLadder were identical across competitors → 0 IDs@1%). The √-normalized Y term and intensity floor are directly transplantable.
**AVOID:** (a) their **peptide-first open search** — conflicts with andes's glycan-Y-first thesis and gives up our precursor-anchored generation edge; (b) folding glycan evidence into the *peptide* FDR pile (their two-stage separation is the point — our unified attempt already failed); (c) treating the α constants as fixed — andes should **re-learn** them regime-matched (stepped-HCD, own models) rather than copy Table S1 values, staying differentiated and license-clean.
