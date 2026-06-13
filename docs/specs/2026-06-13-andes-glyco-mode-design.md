# andes `--glyco` Mode — Design Spec (N-glyco, stepped-HCD MVP)

**Status:** design. Built on the feasibility analysis in
`docs/specs/2026-06-13-andes-glycopeptide-search-research-plan.md` and the
neutral-loss scoring primitive shipped on `feat/enzyme-support`
(`docs/specs/2026-06-13-andes-neutral-loss-glyco-design.md`).

> **⚠ READ ALONGSIDE — v2 addendum is authoritative where they differ:**
> `docs/specs/2026-06-13-andes-glyco-mode-design-v2-improvements.md` (from a
> 6-agent SOTA research + synthesis pass) supersedes specific decisions below.
> Key replacements: oxonium trigger → calibrated two-part AND gate (not
> count-only); backbone solver → Y-complementary voting + core-Y quorum + bounded
> top-K fallback + backbone-b/y arbiter (not a fragile direct Y0/Y1 read); 2D FDR
> → **sequential TDC** with `max(q_pep,q_gly)≤α` (not inclusion-exclusion);
> scoring → the neutral-loss primitive scores only **glycosite-spanning,
> partial-glycan-retaining** backbone b/y with **per-rung trim masses** (NOT the
> full glycan — a full-glycan loss goes `theo_mz≤0` and is silently dropped on
> singly-charged ions), with the Y0/Y1/Y2 whole-backbone ladder + oxonium as a
> separate cheap matcher feeding **additive PIN features** (numeric, `--glyco`-
> gated; strings TSV-only); Phase-0 adds a precision arm + a *searchable-backbone*
> gate + a DB-coverage check, and a precondition (the Byonic per-scan truth must
> be obtained/verified, with a pGlyco3/MSFragger fallback). Plan from the addendum.

**Goal:** add an intact-N-glycopeptide search mode to andes for stepped-HCD/HCD
data — reusing the existing peptide engine — that identifies the peptide
backbone AND its N-glycan composition with glyco-aware FDR. Gated on a
measurement-only prototype that must pass first.

**Scope (MVP):** N-linked glycosylation only; sequon-restricted (N-X-S/T, X≠P);
stepped-HCD or HCD; curated N-glycan composition list (not de-novo, not a full
glycan tree). **Out of scope (MVP):** O-glyco, ETD/EThcD, de-novo glycan,
glycan structure/linkage. Deferred to a follow-on.

**Why a mode, not a new engine:** an intact glycopeptide is still ONE peptide —
the glycan is a single labile precursor adduct on a normal backbone, so andes's
core invariant `suffix = peptide_nominal − prefix` holds *at the backbone mass*.
No two-peptide co-scoring, no n² pairing (contrast: crosslinking). The candidate
engine, `score_psm`, reverse-decoy FDR, Percolator PIN, and the neutral-loss
primitive are all reusable.

---

## Phase 0 — Prototype (go/no-go; measurement only; NO engine change)

The make-or-break question: **can we detect glyco scans and recover the peptide
backbone mass from the spectrum alone, accurately enough to drive a search?**
Answer it on real data before building anything in the engine.

- **Data (already staged):** `/tmp/glyco_andes/Pool_HCC_early_Fc3_1.mzML`
  (PXD025455, stepped-HCD, Q Exactive HF). Ground truth: the dataset's own Byonic
  results — download `Pool_HCC_early_Fc3_1.pepXML` from the same PRIDE folder
  (`.../2021/05/PXD025455/`). It gives, per scan, the Byonic glycopeptide
  sequence + glycan, i.e. the true backbone mass and the true glyco/non-glyco label.
- **Tool:** a standalone **Python** script (`pyteomics` for mzML + pepXML). No
  andes/Rust changes — keeps Phase 0 purely diagnostic.
- **Detector 1 — oxonium classifier:** flag an MS2 as glyco when ≥2 of the
  canonical oxonium ions are present above an intensity floor (fraction of base
  peak). Ion set (m/z): 126.0550, 138.0550, 144.0655, 168.0655, 186.0761,
  204.0867 (HexNAc series); 274.0921, 292.1027 (NeuAc); 366.1395 (HexNAc+Hex);
  657.2349 (HexNAc-Hex-NeuAc). Tune the floor + count on the data.
- **Detector 2 — Y0/Y1 backbone-mass reader:** in the Y-ion region, find the
  Y0 (bare peptide) / Y1 (peptide+HexNAc, +203.0794) pair and the Hex ladder
  (+162.0528). Derive backbone neutral mass from Y0 (or Y1−HexNAc), then
  `glycan_mass = precursor_neutral − backbone_neutral`; snap to the composition
  list (below).
- **Validation metrics (the go/no-go):**
  1. **Trigger recall/precision:** of Byonic-confirmed glyco scans, what fraction
     does the oxonium classifier flag? (target recall ≥ ~0.85) and how many
     non-glyco scans does it wrongly flag?
  2. **Backbone-mass accuracy:** for flagged scans with a detectable Y0/Y1, does
     the derived backbone mass match Byonic's peptide mass within tolerance
     (≤ ~0.02 Da)? (target ≥ ~0.7 of flagged glyco scans)
  3. **Glycan-mass accuracy:** does `precursor − backbone` snap to the correct
     Byonic glycan composition?
- **Decision:** if recall + backbone accuracy clear the bar → proceed to Phase 1.
  If Y0/Y1 is too sparse to pin the backbone, reconsider (e.g. require the
  oxonium trigger + fall back to Route (ii) bounded composition search). Record
  the numbers in this doc.

---

## Phase 1+ — the `--glyco` search mode (outline; gated on Phase 0)

Per glyco-flagged MS2 (oxonium-triggered):

1. **Backbone mass** from Y0/Y1 (Phase-0 reader, ported to Rust).
2. **Peptide search at the backbone mass** — reuse `candidate_gen` +
   `score_psm` UNCHANGED. The sequon (N-X-S/T) restricts candidates to those
   with a glycosite. This is an ordinary andes search window; the GF mass-indexed
   DP works as-is because the backbone is a normal peptide.
3. **Glycan assignment** — `glycan_mass = precursor − backbone`; match against the
   curated N-glycan composition list (~182 human N-glycans; see "Glycan database
   & decoy strategy"), disambiguated by Y-ladder consistency.
4. **Glycan-fragment scoring** — reuse the **neutral-loss primitive**: the Y-ion
   ladder rungs and glycan-retaining backbone ions ARE the loss-shifted b/y ions
   that `predict_by_ions_with_losses` emits and `score_psm` scores against the
   pooled per-class glyco loss table. (ETD activation-gating already correct.)
5. **2D FDR** — peptide-level (existing reverse decoy) × glycan-level (decoy
   glycan list / FragPipe-style), reported jointly.
6. **Output** — additive glyco PIN/TSV columns (glycan composition, glycan score,
   Y-ion count); the existing `Modifications` TSV column carries the glyco mod.

**Reuse map** (detail in the research plan §"architecture map"):
- REUSE as-is: candidate generation, `score_psm`, reverse-decoy, Percolator PIN, neutral-loss primitive.
- NET-NEW: oxonium trigger, Y0/Y1 backbone solver, glycan composition DB + assignment, glycan decoy + 2D FDR, glyco output columns.
- Likely housing: an `andes-glyco` crate (DB + FDR) + a `--glyco` mode flag wiring it into the search loop.

---

## Independence
Clean-room implementable (the method is multiply, independently published). Read
O-Pair/MetaMorpheus (MIT) + Glyco-Decipher (Apache-2.0) freely; papers-only for
pGlyco3 / MSFragger-Glyco / Byonic / StrucGP. Steer clear of Protein Metrics'
wildcard-modification patent — the curated-composition MVP does.

## Glycan database & decoy strategy (resolved from SOTA)

Every leading tool follows the same pattern — **a curated default composition
list + a user-overridable custom file, searching glycan *compositions* (not
structures), with mass-based glycan decoys.** Sizes: Byonic ships 57 / 182 / 309
human N-glycan lists; MSFragger-Glyco's default is a **182-mass** N-glycan list
(loads Byonic/MetaMorpheus/pGlyco files for custom); pGlyco3 ships 1,234
compositions. So the MVP adopts that consensus:

- **Default DB:** ship a curated **~182 human N-glycan composition list** under
  `resources/glycans/human_nglycan_182.txt` (the de-facto-standard size — broad
  enough to cover most serum/tissue N-glycans, small enough to stay fast). Each
  entry = composition (`HexNAc(n)Hex(m)NeuAc(k)Fuc(j)…`) + its monoisotopic mass.
- **Override:** `--glycan-db <file>`; accept a simple composition/mass format and
  (stretch) Byonic-style composition strings for interoperability.
- **Search unit:** glycan **composition** (mass + monosaccharide counts), matched
  by `precursor − backbone` mass + Y-ladder consistency. Not glycan structure/linkage.
- **Glycan decoys (glycan-level FDR):** **mass-based**, à la MSFragger-Glyco /
  PTM-Shepherd (a decoy glycan whose mass sits within the match tolerance of a
  target glycan but is not a real composition). Orthogonal to andes's existing
  reverse-peptide decoy → the two compose into the 2D FDR. Simplest effective
  strategy and matches SOTA. (pGlyco's "Y-complementary mass" = precursor − Y is
  the same backbone-mass idea as our Y0/Y1 reader.)

## Open questions (resolve during planning)
- Exact oxonium floor + count and Y0/Y1 tolerances — set empirically from Phase 0.
- Stepped-HCD only, or also single-energy HCD (Y-ladder may be sparser)?
- Whether to support O-glyco composition files later (the DB format should not preclude it).
