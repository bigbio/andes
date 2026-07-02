# a commercial glyco engine (Protein Metrics) — N-glycopeptide search, distilled for andes

**License: COMMERCIAL, closed-source.** No code, no clean-room borrowing of implementation.
Papers/patents/docs only. Everything below is from *published* sources; treat it as
*prior-art awareness*, not a template to copy. Where a clean-room analogue exists in an
Apache/permissive tool (a glyco search engine, a cross-spectrum glyco engine, O-Pair, MCP multiattribute-FDR paper),
that is the safe reference to actually implement against.

## 1. Candidate generation — PEPTIDE-FIRST, brute enumeration
a commercial glyco engine is **peptide-first / glycan-as-modification**, *not* glycan-Y-first. It in-silico
digests the protein DB, then combines each candidate peptide with **every glycan composition
in a user glycan DB**, predicts a theoretical intact-glycopeptide spectrum, and compares it to
any observed spectrum whose precursor mass lands within tolerance ([Bern 2012, PMC3545648](https://pmc.ncbi.nlm.nih.gov/articles/PMC3545648/);
[community eval, Nat Methods 2021](https://www.nature.com/articles/s41592-021-01309-x)).
The glycan is treated as a variable modification on the N of the N-X-S/T sequon
(N-linked). There is **no fragment-index and no Y-ladder-first backbone deduction** in the
core DDA engine — the peptide-mass anchor comes from `precursor − glycanDB[i]`, exactly andes's
SP-A DB-branch. (a commercial glyco engine's MS3/Y1 method US9484193 is a *separate* instrument workflow, not the
DDA engine.) **Contrast for andes:** andes's glycan-Y-first index (G1, +7–10 pts backbone-
findability) is a genuine differentiator a commercial glyco engine lacks; keep it.

## 2. Glycan handling
- **DB source:** curated composition lists, not structures. a commercial glyco engine ships e.g. **132 human
  N-glycans, 57 human plasma N-glycans, 9 O-glycans** ([Nat Methods 2021](https://www.nature.com/articles/s41592-021-01309-x));
  users can supply their own.
- **Composition, not structure** (like andes's 2510-comp DB). Notation §5.
- **Y-ladder + oxonium** both scored as fragment evidence: Y0 (bare peptide), Y1 (peptide+HexNAc),
  and glycan-loss Y-ions; oxonium (204.087 HexNAc, 366.14 HexHexNAc, 274/292 NeuAc, etc.) used
  as glyco-class evidence ([Bern 2018 peak-filtering, MCP](https://www.sciencedirect.com/science/article/pii/S1535947620351252)).

## 3. Scoring — probabilistic fragment-match, hand-tuned weights
- **a commercial glyco engine score** = a −10·log10(p-value)-style fragment-match score, 0–~1000; **300 good, 400
  very good, >500 near-certain** ([Protein Metrics docs](https://support.proteinmetrics.com/hc/en-us/articles/11608104962836-Delta-Mod-Score)).
  It counts/weights matched b/y (and c/z under EThcD) against a chance model of a search of that
  size (the "Log Probability" = log p-value).
- **Weighting rule (published guidance, not learned):** *peptide backbone fragments weighted
  MORE than Y-ions, and much more than oxonium, regardless of oxonium abundance* ([Bern 2018](https://www.sciencedirect.com/science/article/pii/S1535947620351252)).
  So a commercial glyco engine's combined score ≈ w_pep·S_bly + w_Y·S_Y + w_oxo·S_oxo with hand-set w_pep ≫ w_Y ≫ w_oxo.
- **Delta Mod Score** = drop in a commercial glyco engine score from top peptiform to 2nd-best; a *localization/
  ranking* confidence, 0…Score. High Δ ⇒ confident site/assignment. This is a commercial glyco engine's answer to
  the exact problem andes has (competing candidates at one backbone mass): it is a **top-1-minus-
  top-2 separation feature**, not a new absolute score.
- **No learned spectral predictor.** Weights are static. **This is andes's opening:** replace
  hand-tuned w_pep/w_Y/w_oxo + a commercial glyco engine's chance model with andes's *own regime-matched learned
  peptide-axis model* (SP-B) → Percolator features. Differentiate, don't clone.

## 4. FDR — 1D peptide FDR only (a commercial glyco engine); glycan-axis FDR came later, elsewhere
- **a commercial glyco engine itself does PEPTIDE FDR only** (PEP/1D-FDR). Practical glyco cutoffs used in the
  field: **Score > 300, |Log Prob| > 4, 1D-PEP ≤ 0.001, mass error < 10 ppm** ([community usage](https://www.nature.com/articles/s41592-021-01309-x)).
  a commercial glyco engine has **no glycan-specific decoy** — a known gap.
- **Glycan-axis / 2D-FDR** is published by the *same lineage* (Bern-adjacent, Protein Metrics)
  in the **Multiattribute Glycan Identification** paper ([MCP 2022, PMC8933705](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933705/)),
  and **this is the clean, citable, formula-level source andes should mirror**:
  - **Decoy-glycan recipe:** shift a target glycan's intact mass by a *random value within the
    glycan mass-error tolerance* + assign a *random isotope error*; each glycan **fragment ion
    randomly shifted by a unique 1–20 Da**. Decoys stay near target mass but have distinguishable
    mass/isotope-error distributions → separable by Y+oxonium evidence. *(This is the recipe
    andes's G3′ should adopt — cleaner than andes's current same-backbone-only decoy.)*
  - **Pairwise glycan score (Eq. 1):**
    `S_pairwise = S_Y + S_oxo + β_mass·log|Δm₂/Δm₁| + α_isotope`
    `S_Y  = Σ_t 2( α^hit_{Y,t}(U₁−U₂) + α^miss_{Y,t}(V₁−V₂) )`, t ∈ {HexNAc, Hex-only, Fuc-containing}
    `S_oxo = Σ_t 5( α^hit_{oxo,t}(I_obs/I_exp)(U₁−U₂) + α^miss_{oxo,t}(V₁−V₂) )`, t ∈ {NeuAc/NeuGc, Fuc, phosphate, sulfate}
  - **Absolute score** (for FDR): same terms but Δm/isotope compared to *typical* errors, not a
    rival candidate. **FDR** = collect absolute scores of all target & decoy best candidates, pick
    the threshold giving the target decoy ratio.
- **HARD andes constraint:** FDR = **Percolator only**. Do NOT implement a commercial glyco engine-style internal
  target-decoy counting. Use the *decoy-glycan recipe* + score *terms* as **PIN features** and a
  **thin separate-axis 2D post-process** (a glyco search engine style) on Percolator output. Feeding glycan-decoy
  rows under one `Label` already crashed andes 29%→4% — keep the axes separate.

## 5. Glycan notation emitted
`HexNAc(4)Hex(5)Fuc(1)NeuAc(2)` — monosaccharide-count string. Keywords: HexNAc, Hex, Fuc/dHex,
NeuAc, NeuGc, Pent, GlcA, IdoA, DiNAcBac, Acetyl, Sulfo, Phospho, Na ([Protein Metrics N-glyco docs](https://support.proteinmetrics.com/hc/en-us/articles/17137163463316-a commercial glyco engine-N-Linked-Glycopeptide-Analysis)).
andes already uses HexNAc/Hex/Fuc/NeuAc/NeuGc — **adopt this exact parenthesized string for
interoperability** (a glyco search engine/GlyGen-compatible).

## 6. License
**a commercial glyco engine = commercial, closed (Protein Metrics / Dotmatics).** Clean-room boundary: only the
*published algorithm ideas & formulas* above are usable, and prefer the **Apache/permissive
analogues** for actual implementation — a glyco search engine (GPL-adjacent academic; ideas only, check
license), a cross-spectrum glyco engine (cross-spectrum, Apache), O-Pair/an open-source glyco engine (permissive), and the
**MCP multiattribute-FDR paper** (formulas are publishable prior art). Do **not** read/port a commercial glyco engine
or the reference glyco engine (UM-proprietary) source.

## 7. Single most valuable clean-room idea + what to AVOID
- **BORROW:** the **Delta-Mod-Score concept** — a *top1-minus-top2 backbone-score separation*
  feature — as a **Percolator PIN feature**. It directly attacks andes's #1 failure (ranking loses
  71/154; competing peptides share all spectrum-level glyco features). Compute Δ = (best backbone
  b/y score) − (2nd-best *different-peptide* backbone b/y score) per scan; it discriminates
  target-vs-decoy at one backbone mass where OxoniumScore/YLadderScore cannot. Pair it with SP-B's
  learned peptide-axis model and the MCP decoy-glycan recipe for G3′.
- **AVOID:** (1) a commercial glyco engine's **hand-tuned static score weights + internal 1D target-decoy** — andes
  must stay learned (SP-B) + Percolator-only. (2) **Peptide-first brute enumeration as the primary
  path** — it is the combinatorial-false-match generator that gave andes ~20 candidates/scan and a
  1.2:1 target:decoy pile; keep glycan-Y-first (G1) as the pruning front-end. (3) a commercial glyco engine's **lack of
  glycan-axis FDR** — don't inherit the gap; do G3′.
