# Crosslinking (XL-MS) support in andes — Research Plan

**Date:** 2026-06-13
**Status:** Research / feasibility plan — not yet scoped for implementation
**Author basis:** feasibility study of the andes codebase + XL-MS literature/engine survey

> This is a **research roadmap**, not an implementation spec. It records the
> landscape, a feasibility verdict against andes's actual architecture, a phased
> approach, and what to prototype first. Implementation of any phase requires its
> own brainstorm → spec → plan cycle.

## 0. Verdict up front

Crosslinking is a **large** feature that does **not** fit cleanly into andes's
single-peptide architecture. The scoring core assumes one peptide per spectrum
(`prefix_mass + suffix_mass == peptide_nominal`); a crosslinked spectrum is
co-explained by **two** peptides joined by a linker, which breaks that invariant
at the root, along with the candidate index, the PIN output schema, and the
decoy/FDR model.

- A **cleavable-crosslinker (DSSO/DSBU) MS2–MS3 MVP** is a realistic,
  independence-clean, mostly-reuse extension: detect reporter doublets → deduce
  each peptide's mass → run **two ordinary linear searches** (stub as a residue
  mod) → pair + FDR. This rides on andes's existing engine.
- A **non-cleavable (DSS/BS3) pair search** is effectively a **different engine**
  (O(n²) pairing + a two-peptide co-scorer + a new index) and should live in a
  separate `andes-xl` crate/mode if pursued at all.

**Recommendation:** pursue **Phase 0 + Phase 1** if XL is strategic; gate
non-cleavable pair search on a Phase-1 success and a scope/IP review.

## 1. Motivation

XL-MS maps protein structure and interactions by covalently bridging spatially
close residues, digesting, and identifying the crosslinked peptide pairs.
Supporting it would extend andes from linear PSM identification into
structural/interactomics proteomics. The question this plan answers: *what would
it take, and is it the right shape for andes?*

## 2. Background — XL-MS for an engine designer

- **Crosslinker** = a bivalent reagent bridging two residues (usually K–K / K–Nterm
  for NHS-ester reagents). The species reaching the MS is a **covalent dipeptide**
  (peptide α + peptide β + linker residual), linked at specific internal residues.
- **Product types:** *inter-peptide* (two distinct peptides — the target of XL-MS),
  *intra-peptide loop-link* (two residues in one peptide), *dead-end/mono-link*
  (one end hydrolyzed — behaves like a single-residue mass mod).
- **Non-cleavable (DSS, BS3):** linker stable in MS2; the dipeptide fragments as
  one entity — both peptides' b/y ions appear in one spectrum, and any fragment
  spanning the linked residue carries *the entire partner peptide + linker* as a
  mass adduct. No diagnostic ions. This is the hard combinatorial problem.
- **MS-cleavable (DSSO ~158, DSBU/BuUrBu ~196, CDI, DSBSO):** a labile bond in the
  spacer cleaves in MS2 *before* backbone fragmentation, releasing each peptide
  with a small characteristic **stub** and producing **reporter doublets** (peak
  pairs at a fixed Δ-mass; DSSO's classic Δ≈31.97 Da stub family). The doublet (a)
  flags a crosslink spectrum, (b) lets the engine *measure* each peptide's mass,
  (c) reduces the problem to sequencing two near-linear peptides carrying a stub.
- **Acquisition:** **MS2-only** (stepped-HCD; doublets + both backbone ladders in
  one scan) vs **MS2–MS3** (MS2 cleaves the linker + measures doublets; each
  released peptide gets an MS3 sequenced like an ordinary linear peptide). MS2–MS3
  is the most andes-friendly (each MS3 ≈ a linear spectrum + stub mod).

## 3. State of the art (engines, algorithms, licensing)

**Candidate strategy — taming the O(n²) pair explosion:**
- **Top-K α-then-β:** score single peptides, keep top-K (pLink ~500, ProteinProspector
  ~1000, Kojak ~250), pair only among those. ([pLink2, Nat Commun 2019](https://www.nature.com/articles/s41467-019-11337-z))
- **pLink2 fragment-index + two-stage open search**; **Kojak 2.0** two-stage with
  modified Comet XCorr; both prune aggressively. ([Kojak 2.0, PMC10234491](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10234491/))
- **Cleavable-reporter narrowing (XlinkX, MeroX, MS Annika, MaXLinker):** detect
  doublets → measure each mass → collapse to **two linear searches**. ([XlinkX, Nat Commun 2016](https://www.nature.com/articles/ncomms15473); [MaXLinker, MCP 2020](https://www.mcponline.org/content/19/3/554))
- **Exhaustive references:** StavroX, **OpenPepXL (OpenMS, BSD-3)**.

**Scoring two-peptide spectra:** predict both peptides' b/y ions (linked-residue
fragments carry partner/stub mass); combine per-peptide evidence. **MaxLynx**
uses a dipeptide Andromeda score with **partial + total** scores (handles the
"one peptide strong, one weak" asymmetry). **Kojak** sums modified XCorr; **pLink**
uses a joint matched-ion p-value; **Xolik** shows linear-time paired scoring. ([MaxLynx, Anal Chem 2021](https://pubs.acs.org/doi/10.1021/acs.analchem.1c03688); [Xolik, Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/2/251/5047755))

**FDR:** each ID has two peptides → three classes **TT / TD / DD**; standard
`FDR = (TD − DD) / TT`, with **separate inter- vs intra-protein** FDR. **MS Annika**
and **MaxLynx** score best on correct FDR + ID count in benchmarks. XL decoys are
typically **reverse-between-cleavage-sites** to preserve linked-residue mass. ([MS Annika, PMC8155564](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8155564/); [ribosomal benchmark, Nat Commun 2022](https://www.nature.com/articles/s41467-022-31701-w))

**Licensing (relevant to andes's Apache-2.0 goal):** **pLink2 is license-gated**
(per-user time-limited licenses) — *not* a code/algorithm source to borrow from. ([pLink2 GitHub](https://github.com/pFindStudio/pLink2)) **Kojak** (Apache-friendly) and
**OpenPepXL** (BSD-3) are the permissive references. The general methods (top-K
α-then-β, reverse-between-cleavage decoys, TT/TD/DD FDR, cleavable-doublet →
two-linear-searches) are widely published across independent engines and
implementable clean-room — **but a patent check is required before committing to
any specific method.**

## 4. Feasibility against andes's architecture

| Concern | andes today | XL impact |
|---|---|---|
| Crosslinker representation | `Modification` = single-residue scalar mass delta ([modification.rs:28-43](../../crates/model/src/modification.rs)) | Dead-end & loop-link fit as mods; **inter-peptide link is bivalent/relational — needs a new first-class `Crosslinker`/`CrosslinkSite` concept**, cannot be a `Modification`. |
| Candidate generation | enumerates **single** peptides; no mass-sorted index, no pairing ([candidate_gen.rs:40-203](../../crates/search/src/candidate_gen.rs)) | Non-cleavable pairing is O(n²) and **needs infrastructure andes lacks** (mass index, top-K, "find β of mass `precursor−α−linker`"). |
| Scoring core | `score_psm` walks one peptide with `suffix = peptide_nominal − prefix` ([psm_score.rs:234-311](../../crates/scoring/src/scoring/psm_score.rs)); per-(partition,IonType) rank LLR; single-peptide node-cache ([scored_spectrum.rs:790-877,416-426](../../crates/scoring/src/scoring/scored_spectrum.rs)) | **Structural blocker.** Two peptides co-explain one spectrum; `suffix = total − prefix` is false; node-cache doesn't transfer; peak-ownership is new. Per-ion LLR tables are *reusable for unshifted ions* but uncalibrated for XL spectra. |
| Output / FDR | one `Peptide`/`Proteins`/mass per PIN row ([pin.rs:11-20,342](../../crates/output/src/pin.rs)); whole-protein-reversal decoys; Percolator FDR | PIN can carry XL *features* additively, but needs **two peptides per row** and a **TT/TD/DD inter/intra FDR layer** beyond Percolator + **reverse-between-cleavage decoys**. |
| Cautionary prior | the **fragment-index candidate generator was built and abandoned** (git `7095e061`, `206d1ae0`; CLAUDE.md "Abandoned") with an "irreducible recall/speed tension — a top-K fragment prefilter drops exactly the secondary co-isolated peptides" | This is the **same failure mode** as XL pairing (β is the "secondary" peptide a top-K α prefilter discards), and andes deliberately *deleted* the only mass-index-like structure. Re-introducing it must re-validate that finding. |
| Interaction w/ neutral-loss work | `loss_class: u8` IonType discriminant + per-mod fragment emission (in progress) | **Helps modestly:** the cleavable stub is a residue mass-delta mod; the per-mod fragment-emission plumbing is the same *shape* needed for stub-shifted linked fragments; the IonType discriminant proves the model store tolerates extra per-ion classes. Not the same feature; no collision. |

## 5. Phased roadmap

- **Phase 0 — mono-links & loop-links (days, no architecture change).** Confirm/document
  that **dead-end mono-links** search today as `mods.txt` deltas. **Loop-links** additionally
  need the b/y ladder to skip fragments *between* the two linked residues — a small
  fragment-prediction tweak. Ship as "limited crosslink-adduct support."
- **Phase 1 — DSSO/DSBU MS-cleavable MS2–MS3 MVP (the real first project; medium-large,
  isolated in `andes-xl` crate / `--xl` mode).** (1) reporter-doublet detection; (2) per-peptide
  mass deduction from doublets; (3) **two stub-modified linear searches reusing `score_psm` +
  `candidate_gen`**; (4) pairing + mass-consistency check (`α + β + linker ≈ precursor`);
  (5) crosslink output schema + **TT/TD/DD inter/intra FDR** + reverse-between-cleavage decoys.
  Reuses the linear engine; avoids the n² explosion and the two-peptide-scoring rewrite.
- **Phase 2 — MS2-only cleavable (stepped-HCD).** One scan holds doublets *and* both backbone
  ladders → requires two-peptide co-scoring. Only after Phase 1 proves the FDR/pairing/output
  layer.
- **Phase 3 — non-cleavable DSS/BS3 pair search (large; reconsider whether to do at all).**
  Needs a mass-sorted peptide index, top-K α-then-β candidate generation, and the **two-peptide
  co-scorer with peak ownership** — i.e. re-introducing the indexing andes deleted and a scorer
  the per-(partition,IonType) model only partially supports. Honest framing: **a separate XL
  engine.**

## 6. Prototype & benchmark first (make-or-break, before committing to Phase 1)

1. **Doublet detector** on a public DSSO dataset — the synthetic peptide library
   ([Beveridge et al., Nat Commun 2020](https://www.nature.com/articles/s41467-020-14608-2.pdf))
   is the standard benchmark. Measure detection sensitivity/specificity. This is the make-or-break
   of Phase 1.
2. **Two-stub-linear-search + naive pairing** on MS2–MS3 data, scored against **MS Annika / XlinkX**
   on the same file: does andes's linear scorer + simple pairing recover known crosslinks at a
   controlled FDP? Tells you whether the linear core suffices before investing in two-peptide scoring.

## 7. Hardest parts (ranked)

1. **Two-peptide co-scoring with peak ownership** (Phase 2/3) — the structural rewrite of the scoring core.
2. **Re-building the mass index / α-then-β prefilter** andes abandoned, avoiding the known
   "top-K drops the secondary peptide" failure (Phase 3).
3. **Correct XL FDR** (TT/TD/DD, inter/intra, XL decoys) — a new layer beyond Percolator (all phases).
4. **Robust doublet detection** across instruments/linkers (Phase 1 gate).

## 8. Independence & licensing

XL can fit the Apache-2.0 independence goal **if** andes uses its own algorithms and
permissively-licensed references only. **Avoid pLink/pLink2** (license-gated; do not copy
algorithm details). Implement clean-room from published methods used across multiple
independent engines (Kojak/Apache-friendly, OpenPepXL/BSD, MS Annika). **Required: a patent
review** of the specific cleavable-doublet / pairing / FDR methods before implementation — no
blocking patent on the *general* cleavable method was found, but no exhaustive IP search was done.

## 9. Open questions / decisions before any implementation

- Is XL strategically in scope for andes, or out of mission (linear ID + independence + glyco)?
- Cleavable-only (Phase 0–1) acceptable as the product, or is non-cleavable required?
- Separate `andes-xl` crate vs `--xl` mode in the main binary?
- Which crosslinker(s) and instrument method(s) to target first (DSSO MS2–MS3 recommended)?
- Datasets available for the Phase-1 prototype/benchmark?

## 10. References

- pLink2 (license-gated): https://www.nature.com/articles/s41467-019-11337-z · https://github.com/pFindStudio/pLink2
- Kojak / Kojak 2.0: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4428575/ · https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10234491/
- XlinkX (MS-cleavable, MS2–MS3): https://www.nature.com/articles/ncomms15473
- MaXLinker: https://www.mcponline.org/content/19/3/554
- MaxLynx (dipeptide Andromeda, partial+total scores): https://pubs.acs.org/doi/10.1021/acs.analchem.1c03688
- MS Annika (FDR) / 2.0 (MS2-MS3): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8155564/ · https://pubs.acs.org/doi/10.1021/acs.jproteome.3c00325
- Xolik (linear-time paired scoring): https://academic.oup.com/bioinformatics/article/35/2/251/5047755
- XL benchmark (ribosome): https://www.nature.com/articles/s41467-022-31701-w
- DSSO synthetic-library benchmark: https://www.nature.com/articles/s41467-020-14608-2.pdf
- Thermo MS-cleavable crosslinkers PI: https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0016303_MSCleavableCrosslinkers_PI.pdf

### andes code references
- single-peptide candidate enumeration (no pairing/mass index): `crates/search/src/candidate_gen.rs:40-203`
- `Modification` scalar single-residue delta: `crates/model/src/modification.rs:28-43`
- single-peptide scoring invariant `suffix = peptide_nominal − prefix`: `crates/scoring/src/scoring/psm_score.rs:234-311`
- per-(partition,IonType) rank model + node DP cache: `crates/scoring/src/scoring/rank_scorer.rs`, `scored_spectrum.rs:790-877,416-426`
- canonical b/y only: `crates/scoring/src/scoring/fragment_ions.rs:1-3,98-144`
- single-peptide PIN row: `crates/output/src/pin.rs:11-20,342`
- abandoned fragment index (the prefilter XL needs) + failure mode: git `7095e061`, `206d1ae0`; `.claude/CLAUDE.md`
