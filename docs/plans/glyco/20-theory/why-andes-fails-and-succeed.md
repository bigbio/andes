# Why andes glyco identification fails — and the mechanism of success

*Theoretical study. Branch `glyco-phase1`. Truth set: PXD025455 `HCC_pool_Late_Fc3_r1`,
523 re-searched N-glyco scans. Source-cited; clean-room references only.*

## 1. The paradox: generation is solved, identification is not

andes recovers the true peptide **backbone** in ~90% of glyco spectra
(`PHASE1_RESULT.md`: 90.4% searchable via `precursor − known_glycan`) yet only
**29.4% pass @1% FDR (154/523)**, of which **83/154 are top-1 correct** and only
**66 are true @1% FDR (111 false pass)**. The loss is entirely in **ranking +
FDR**, not candidate generation. This is the diagnostic decomposition to keep:
`IDs ≈ Coverage × Separation`. Coverage ≈ 0.90; the deficit is **Separation**.

## 2. Root cause — the SPA2 insight: glyco features are backbone-level

The mechanistic finding (`SPA2_RESULT.md`) is decisive. For a spectrum with
precursor-neutral mass *M*, andes enumerates ~20 candidate backbones as
*M − glycan_i* over 2,510 clean-room compositions, each matched to sequon peptides.
The glyco-specific PIN features — **OxoniumScore, YLadderScore, CoreYHits,
GlycanMass** — are computed from oxonium ions (204.087, 366.14, …) and the Y-ladder
(peptide+glycan-remnant series). **These peaks depend only on the glycan and the
intact precursor mass, not on the peptide sequence.** Two different peptides
competing at the *same* backbone mass window produce *identical* oxonium/Y-ladder
features. So these features contribute **zero target/decoy discrimination on the
peptide axis**. Percolator is left with only the backbone **b/y RankScore**, which
on stepped-HCD glyco spectra is physically sparse — collision energy is spent
cleaving the labile glycosidic bonds, depositing little backbone-fragment intensity.
Result: top-1 collapse and, when the raw PIN was fed as one unified pile, **0 IDs
@1% FDR** at a ~1.24:1 target:decoy ratio. The features that fire are the ones that
cannot separate; the feature that could separate barely fires.

This mirrors field consensus: direct b/y interpretation succeeds in only **~11%** of
glycopeptide spectra — the basis of a cross-spectrum glyco engine's design
([Yang et al., *Nat Commun* 13:1900, 2022](https://www.nature.com/articles/s41467-022-29530-y)).
A better-calibrated model reshapes a likelihood over the *same* sparse peaks; it
cannot manufacture b/y ions the fragmentation never deposited.

## 3. Recoverable signal — the Y0/Y1 peptide-mass anchor

There is one peptide-axis-discriminating signal hiding in the glyco peak forest.
The **Y0** (bare peptide backbone, +H) and **Y1** (peptide+GlcNAc) ions encode the
**peptide monoisotopic mass directly**, and in stepped-HCD they are typically
**high-intensity** even when interior b/y is dead. Because Y0/Y1 are a function of
the *peptide* mass (not the glycan), an observed-Y0-mass-conditioned feature *does*
distinguish competing peptides at one backbone window — unlike oxonium/YLadder. This
is exactly the anchor a glyco search engine's coarse-scoring exploits: for each glycan it computes
backbone = precursor − glycan and scores the Y-ion ladder anchored on that backbone
mass ([Liu et al., *Nat Commun* 8:438, 2017](https://www.nature.com/articles/s41467-017-00535-2)).
The fix is a **learned, regime-matched peptide-axis glyco scorer (SP-B)** that adds
a **Y0/Y1 anchor feature** as an *additive* PIN column (never modifying existing
features — additive-only is the parity-safe rule here). This is andes's
differentiator: **glycan-Y-first candidate selection + own learned fragment models**,
not a re-implementation of any existing engine.

## 4. Ceiling-breaker — cross-spectrum transfer

SP-B raises Separation but is bounded by the ~11% of spectra with usable backbone
evidence. The lever that breaks that ceiling is **cross-spectrum transfer**:
glycoforms of the *same* peptide backbone share a correlated backbone-fragmentation
pattern, so a confident ID on one glycoform transfers backbone evidence to sibling
spectra whose own b/y is dead. a cross-spectrum glyco engine's peptide-code / shared-backbone
matching delivered its **33.5%–178.5%** PSM increase from this transfer, *not* from
scoring ([Yang et al. 2022](https://www.nature.com/articles/s41467-022-29530-y);
fragmentation+elution+glycome-connectivity modeling independently adds ~9.5%,
[Klein et al., *Nat Commun* 15, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11263600/)).
andes implements this **in-process** (`ANDES_GLYCO_CROSSSPECTRUM`), which must be
**RT-gated** to avoid transferring across co-eluting confounders.

## 5. FDR — 2D, and Percolator-only

Glyco FDR is **two-dimensional**: peptide axis and glycan axis fail independently.
a glyco search engine 2.0 estimates the glycan axis with a **decoy Y-ladder** (a random 1–30 Da mass
added to each Y-ion) scored against a *separate* database
([Liu et al. 2017](https://www.nature.com/articles/s41467-017-00535-2)). andes's
`glycan_y_intensity_decoy` scorer is sound (TDD: target > decoy), but pushing ~352K
glycan-decoy rows into Percolator under **one `Label`** crashed recovery 29.4%→4.4%:
they differ from targets only in YLadder, dominate the −1 pile, and Percolator
over-weights YLadder. **Constraint: 2D-FDR must be a thin Percolator *post-process*
on separate axes — never a unified pile, never Mokapot.**

## 6. Cause → Fix

| # | Cause (mechanism) | Fix | Clean-room source (license) |
|---|---|---|---|
| C1 | Oxonium/YLadder/GlycanMass are **backbone-level** → 0 peptide-axis discrimination; only sparse b/y separates | **SP-B** learned regime-matched peptide-axis glyco scorer; glycan-stripped backbone training rows | a glyco search engine coarse/fine scoring, [Liu 2017](https://www.nature.com/articles/s41467-017-00535-2) (a glyco search engine Apache-2.0) |
| C2 | Interior b/y physically sparse (~11% usable) under stepped-HCD | **Y0/Y1 peptide-mass anchor** as *additive* PIN feature (high-intensity, peptide-specific) | a glyco search engine Y-anchor; O-Pair paired-dissociation anchoring, [Lu et al., *Nat Methods* 17:1133, 2020](https://www.nature.com/articles/s41592-020-00985-5) (an open-source glyco engine, permissive) |
| C3 | Single-spectrum scoring cannot rank the sparse-b/y stratum at all | **RT-gated in-process cross-spectrum transfer** (G4) — the ceiling-breaker | a cross-spectrum glyco engine shared-backbone, [Yang 2022](https://www.nature.com/articles/s41467-022-29530-y) (Apache-2.0) |
| C4 | Glycan-decoy rows in one Percolator pile → over-weighted YLadder, recovery crash | **Separate-axis 2D-FDR** as thin Percolator post-process (peptide axis + glycan-Y decoy axis) | a glyco search engine 2.0 two-DB decoy, [Liu 2017](https://www.nature.com/articles/s41467-017-00535-2) |
| C5 | ~20 coincidental (peptide,glycan) pairs/scan inflate false matches | **Glycan-Y-first** candidate selection (`ANDES_GLYCO_YINDEX`; +7.2 pts findability) prunes the pair space before scoring | andes-native (differentiator) |

**Do not copy code** from a commercial glyco engine (commercial) or the reference glyco engine (UM-proprietary,
[academic-only license](https://available-inventions.umich.edu/product/the reference engine-ultrafast-and-comprehensive-identification-of-peptides-from-tandem-mass-spectra)).
Reference algorithms only from a glyco search engine (Apache-2.0), a cross-spectrum glyco engine (Apache-2.0), and
O-Pair/an open-source glyco engine (permissive).

## 7. Prediction

SP-B (C1+C2) is necessary and lifts top-1 from 83/154 toward the ~90% generation
ceiling, but is bounded by the ~11% b/y stratum; C4 converts the raised Separation
into *valid* q-values (kills the 111 false passers). The identification-rate
ceiling-break past ~30% requires **C3 cross-spectrum transfer** — the field's own
+33.5–178.5% came from transfer, not scoring. Gate on **decoy-separated ranking**
(true peptide's b/y-score vs a same-backbone decoy), not find-rate.
