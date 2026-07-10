# Fragmentation physics: why stepped-HCD N-glycopeptide backbone ID is hard

*Theory note grounding the "andes ranking fails" diagnosis in dissociation physics.*
*Companion to `00-context/00-current-state.md`, `PHASE1_RESULT.md`, `SPA2_RESULT.md`.*

## The energy-partition problem

In collisional activation the internal energy deposited into a glycopeptide precursor
does not distribute uniformly across bonds — it flows to the **lowest-barrier
dissociation channels first**. Glycosidic bonds (HexNAc–peptide, Hex–HexNAc,
sialic linkages) are far more labile than the amide (b/y) backbone bonds. So at the
collision energies that fragment *anything*, the glycan cleaves preferentially and
the ion current is spent producing:

- **Oxonium ions** (low-mass glycan reporters: HexNAc 204.087, Hex 163.060,
  NeuAc 292.103, 274.092, 138.055 …) — diagnostic of glycosylation but carrying
  **zero peptide-sequence information**;
- **Y-ions** (peptide + progressively stripped glycan; Y0 = bare peptide+H,
  Y1 = peptide+HexNAc) — a glycan ladder anchored on the intact backbone.

a cross-spectrum glyco engine states it plainly: "the fragmentation of peptide backbones is
typically **suppressed** in HCD-MS/MS analysis due to the **labile nature of
glycosidic bonds**" ([Fang et al., Nat Commun 2022, PMC8990002](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/)).
The dissociation-methods survey quantifies it: "the **majority of signal** in HCD and
sceHCD spectra is in **glycan-related channels**, i.e., Y-type fragments and oxonium
ions" ([Riley, Malaker et al., JPR 2020, PMC7425838](https://pmc.ncbi.nlm.nih.gov/articles/PMC7425838/)).
The peptide b/y series — the only fragments that identify *which peptide* — receive
the residual, and are often too sparse to sequence.

## What stepped collision energy actually buys

Fragment yield is strongly energy-dependent and the two channels peak at different
energies. Glycan cleavage / Y-ions dominate at **low NCE (~20–30%)**; Y-series
intensity is maximal near NCE 30% and the Y-ions **vanish by NCE 50%**, where
backbone b/y finally appear ([JPR 2020](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00218)).
No single energy gives both. Stepped-HCD (sceHCD, e.g. **20/30/40%**) sums 2–3
energies into one spectrum specifically to "produce cleavages of **both** the peptide
backbone **and** the glycosidic bonds" — a deliberate compromise. It *broadens* the
b/y coverage but does not remove suppression: the composite spectrum still allocates
most of its intensity to oxonium + Y, and sialylated glycans make it worse —
"negatively charged sialic residues … hamper the fragmentation of peptide backbones
and lead to **insufficient peptide fragments**" ([PMC8990002](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/)).
Only electron-based methods (EThcD/AI-ETD) distribute signal evenly enough to
reliably sequence the backbone — and those aren't in the target stepped-HCD data.

## The hard ceiling: ~11% direct backbone ID

The decisive number: in a cross-spectrum glyco engine's in-silico deglycosylation stage, peptide
backbones were directly identified in only **140,250 of 1,282,263** oxonium-containing
glycopeptide spectra — **≈11%** ([PMC8990002](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/)).
The other ~89% *fired oxonium* (are glycopeptide spectra) but **lack backbone b/y
sufficient for direct sequencing**. The authors are explicit: "in-silico
deglycosylation … only permits identification of peptide backbones for the
glycopeptide spectra **with sufficient peptide fragment ions**." This is a property
of the *physics of the spectrum*, not of the search engine.

## Implication for single-spectrum scoring — and why andes ranking fails

andes already solves **generation**: `precursor_neutral − known_glycan_mass` places
the backbone within ±20 ppm for ~90% of spectra (`PHASE1_RESULT.md`: DB-branch =
87.3% of searchable hits). The failure is **identification** — choosing the *correct*
peptide among ~20 sequon candidates per scan and controlling FDR
(`SPA2_RESULT.md`, `00-current-state.md`: 154 recovered → 83 top-1 correct → 66 true).

The physics makes this a **hard ceiling for any single-spectrum peptide-axis scorer**,
andes's learned models included. Discrimination between two candidate peptides at the
same backbone mass can come *only* from peptide-specific fragments — the **b/y
series**. But b/y is exactly the channel the fragmentation suppressed: on ~89% of
spectra there simply **aren't enough backbone ions** to rank the true peptide above a
decoy. Oxonium, Y-ladder, glycan-mass, core-Y features are **backbone/spectrum-level**
— identical for the true peptide and a same-mass decoy — so they add **no
target/decoy discrimination** (`SPA2_RESULT.md`: top-1 T:D ≈ 1.24:1 → Percolator
returns 0 @1% FDR). A better-calibrated model reshapes a likelihood over the *same
sparse peaks*; **it cannot manufacture b/y ions the dissociation never deposited.**

Two physics-grounded consequences for the roadmap:

1. **A peptide-mass-anchor feature (Y0/Y1)** is the exception worth training: Y0/Y1
   are high-intensity *even when the interior b/y ladder is dead*, and they *are*
   peptide-mass-specific, so they discriminate competing peptides where oxonium/Y
   cannot. Pair this with the regime-matched strong model (SP-B).
2. **Cross-spectrum transfer is the escape from the ceiling, not scoring.**
   a cross-spectrum glyco engine's +48.4% (91,535 → 135,840) and +53.3% more backbone spectra came
   from **spectrum expansion** — transferring the fragmentation fingerprint of a
   confidently-IDed backbone to co-eluting spectra sharing that backbone within an RT
   window — *not* from a better per-spectrum score. This directly motivates andes G4
   (RT-gated cross-spectrum, in-process). It recovers the sparse-b/y stratum that
   single-spectrum ranking provably cannot.

**Bottom line:** andes ranking fails because ~89% of stepped-HCD glyco spectra are
physically b/y-starved; single-spectrum peptide scoring is capped near the ~11%
direct-ID ceiling. The differentiator is not to out-score the field per spectrum but
to (a) exploit the always-present Y0/Y1 anchor and (b) borrow backbone evidence
across spectra — andes's glycan-Y-first + own learned models + in-process
cross-spectrum design.

## Clean-room provenance

All algorithmic references are from **published papers with permissive licenses**;
no a commercial glyco engine (commercial) or the reference engine (UM-proprietary) code is used.

- **a cross-spectrum glyco engine** — spectrum-expansion / cross-spectrum matching — Apache-2.0
  ([github.com/DICP-1809/a cross-spectrum glyco engine](https://github.com/DICP-1809/a cross-spectrum glyco engine),
  [Nat Commun 2022, DOI 10.1038/s41467-022-29537-5](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990002/)).
- **a glyco search engine / a glyco search engine** — glycan-first + separate-axis (2D) FDR — published methods
  ([Sci Rep 2016 srep25102](https://www.nature.com/articles/srep25102);
  [a glyco search engine 2.0, Nat Commun 2017, PMC5585273](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5585273/)).
- **O-Pair / an open-source glyco engine** — paired-dissociation localization — **MIT** license
  ([github.com/smith-chem-wisc/an open-source glyco engine](https://github.com/smith-chem-wisc/an open-source glyco engine),
  [Nat Methods 2020, PMC7606753](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7606753/)).
- Dissociation physics — [Riley/Malaker et al., JPR 2020 (PMC7425838)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7425838/).

Reuse is **conceptual** (physics + published algorithms). andes differentiates:
glycan-Y-first candidate selection, own learned peptide-axis models, in-process
cross-spectrum transfer — not a re-implementation of the reference engine/a comparison search engine. **FDR remains
Percolator-only** (never Mokapot); 2D-FDR is a thin separate-axis Percolator
post-process (`00-current-state.md` G3′), consistent with the a glyco search engine separate-axis
scheme.
