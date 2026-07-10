# a glyco search engine — clean-room study for andes N-glycopeptide search

Sources: a glyco search engine *Nat Methods* 2021 [PMC8648562](https://pmc.ncbi.nlm.nih.gov/articles/PMC8648562/) / [DOI 10.1038/s41592-021-01306-0](https://www.nature.com/articles/s41592-021-01306-0); a glyco search engine2 [PMC5585273](https://pmc.ncbi.nlm.nih.gov/articles/PMC5585273/); a glyco search engine decoy/FDR [PMC4853738](https://pmc.ncbi.nlm.nih.gov/articles/PMC4853738/); repo [github.com/pFindStudio/a glyco search engine](https://github.com/pFindStudio/a glyco search engine).

## 1. Candidate generation — glycan-FIRST via ion-indexing (this is the crown jewel)
a glyco search engine is **glycan-first**: it identifies the glycan portion *before* the peptide, then filters, then scores peptides only for surviving backbones. The trick is a **glycan ion-indexing table** that runs in **linear O(#peaks) time per spectrum**. Instead of matching Y-ions directly (each Y = peptide + partial glycan, so peptide-mass-dependent), a glyco search engine indexes the **Y-complementary mass**:

```
precursor_mass − peak_mass  ≈  glycan_mass − glycan_Y_mass  (+ mass error)
```

The RHS is **peptide-independent** — it depends only on the glycan composition and which sub-glycan the Y ion lost. So an index keyed on `glycan_mass − subglycan_mass` lets every observed peak vote, in one pass, for **all** glycans whose Y-ladder it is consistent with. Each index entry carries an **extra bit** marking whether that Y corresponds to a **trimannosyl-core** ion. Output per spectrum: matched-ion count and matched-**core**-ion count for *every* glycan, in O(#peaks). Glycans are retained if **≥ n core-Y matched (n=2 for N-glycan, n=1 for O-glycan)**. This is the exact structure andes calls its "glycan-Y index / two-axis retention" (G1) — a glyco search engine validates it and adds the peptide-independence + core-bit refinements andes should adopt.

## 2. Glycan handling
Built-in N-/O-glycan **composition** DBs (canonicalization-based; GlycoWorkbench `.gwp` import for structures). Search is at **composition** level; `a glyco search engineSite` (separate DP algorithm) does site-specific structure localization post-hoc. Oxonium ions gate glyco spectra; **core Y-ladder** (peptide+HexNAc, +HexNAc₂, +Hex, i.e. the trimannosyl core Y0/Y1/…) is the discriminative backbone signal. Notation: `Hex(9)HexNAc(2)Fuc(1)NeuAc(2)` (monosaccharide-count composition) — dense, engine-friendly; andes already emits this.

## 3. Scoring (a glyco search engine = a glyco search engine2 scheme, re-tuned)
Three scores. Peaks weighted by **log intensity × a quartic mass-error term** `(1 − (merr/tol)⁴)`:

- **Glycan:** `ScoreG = Σᵢ log(intᵢ)·(1−(merrᵢ/tol)⁴) · ratio_ion^α · ratio_core^β`  (α≈0.22–0.56, β≈0.42–0.45; re-tuned in v3)
- **Peptide:** `ScoreP = Σᵢ log(intᵢ)·(1−(merrᵢ/tol)⁴) · ratio_ion^γ`  (γ≈0.94)
- **Combined GPSM:** `ScoreGP = w·ScoreG + (1−w)·ScoreP`, **w≈0.35** (peptide weighted 0.65).

`ratio_ion` = matched/theoretical ions; `ratio_core` = matched/theoretical **trimannosyl-core** ions. Note: **no learned/spectral-predictor model** — a glyco search engine is entirely a hand-tuned intensity·error·coverage function. This is exactly the seam andes differentiates on (own learned peptide-axis model + cross-spectrum), so **borrow the decomposition, not the hand weights**.

## 4. FDR — 2D / multi-dimensional (glycan ⟂ peptide), Percolator-compatible as post-process
**Glycan decoy recipe (from PMC4853738):** for each theoretical Y-ion of the deduced backbone, **add a random mass 1–30 Da** to the Y-ion mass → a decoy glycan Y-spectrum; compete target vs decoy Y-ladders. (Backbone/peptide decoys are ordinary reversed/shuffled sequences.) FDR is estimated with a **finite mixture model** to de-bias the spectrum-based decoy.

**2D-FDR** treats a glycopeptide as false if **glycan OR peptide is false**:
```
FDR_GP = P(G=false ∪ P=false | X ≥ x) = FDR_G + FDR_P − FDR_{G∩P}
```
→ glycan-, peptide-, **and** glycopeptide-level control, computed on **separate axes** then combined. **This directly explains andes's G3 crash:** dumping glycan-decoy rows into one Percolator `Label` pile made YLadder dominate the −1 pile. The fix is a glyco search engine's shape — **two independent decoy axes**, each a *thin post-process* of Percolator q-values (peptide-axis TDC on backbone decoys; glycan-axis on the +1–30 Da Y-decoy), intersected per the union formula. Percolator stays the FDR engine; 2D-FDR is arithmetic on top.

## 5. License
Repo code is **Apache-2.0** (clean-room–usable, cite it). BUT the *distributed binary* requires a per-user license from `i.pfind.net/license/a glyco search engine` — operational/runtime restriction, **not** an algorithm-IP restriction. **Algorithms are from published papers ⇒ clean-room re-implementation is fine.** Do not vendor a glyco search engine binaries or GUI; re-derive from the formulas above.

## 6. Single most valuable borrow + what to avoid
**BORROW:** the **glycan-first ion-indexing on Y-complementary mass** (`precursor − peak = glycan − subglycan`, peptide-independent, O(#peaks), core-bit, ≥2-core-Y retention) — it is the exact algorithmic core of andes's G1 and it makes candidate generation both fast *and* discriminative before any peptide scoring. Second: the **2D separate-axis FDR union formula** as andes's G3′ Percolator post-process.
**AVOID:** (a) a glyco search engine's **hand-tuned static weights** (w=0.35, α/β/γ) — andes's edge is *learned* regime-matched scoring + cross-spectrum, so keep the score *decomposition* but learn the combiner; (b) the **finite-mixture-model FDR** — it competes with Percolator; use TDC + the union formula instead; (c) treating scoring as the fix for sparse stepped-HCD b/y — a glyco search engine has no cross-spectrum transfer, which is precisely andes's differentiator (G4).
