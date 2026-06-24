# Glycosylation / neutral-loss scoring in andes — and the MaxSBM (sequence-based-modifier) extension

Status: reference + design note. The neutral-loss/glyco scoring described in §1 is
implemented and on-branch. The MaxSBM-informed extension in §3–§4 is a design direction,
not yet built. Any ID-count claim from it must clear **1% true entrapment-FDP**, not
target-decoy (see §5).

---

## 1. The current andes glyco / neutral-loss algorithm

andes models a modification that **fragments off the backbone** (a glycan, a labile PTM)
as a residue that carries one or more **neutral losses** plus a **`loss_class`** tag, and
scores the resulting loss ions as an **additive, peptide-aware** term on top of the normal
b/y rank score. The design goal is that the term is exactly zero for every standard
peptide, so the integer `RawScore` stays byte-identical to the pre-feature engine.

Key pieces (all in `crates/scoring/src`):

- **`IonType::Prefix { charge, offset_bits, loss_class }`** — the ion vocabulary carries a
  per-class `loss_class: u8`. `loss_class == 0` is the normal (no-loss) series; non-zero
  classes are distinct loss families (e.g. a glycan oxonium/Y-ion family). This lets the
  model keep **separate trained frequency tables per loss class** instead of folding all
  losses into the backbone series.
- **`Modification { mass_delta, neutral_losses: Vec<f64>, loss_class, .. }`** — a mod
  declares the neutral losses it can shed (e.g. hexose `162.0528`) and which `loss_class`
  they belong to.
- **Activation gate** — `data_type.activation.predicts_neutral_losses()`. ETD/EThcD strip
  no labile groups, so loss ions are not predicted there; HCD/CID do. Loss scoring only
  runs when the activation *preserves* the loss.
- **Additive scoring** (`psm_score.rs`, the split loop ~L286–L334). For each fragment split
  position `s`, after the normal `node_score`, a **second** contribution is added **only
  when** `score_losses` holds — i.e. (a) the activation predicts losses, (b) the model has
  trained loss tables (`has_loss_tables()`), and (c) the peptide actually carries a
  loss-bearing mod (`loss_class != 0 && !neutral_losses.is_empty()`). For standard
  peptides this gate is false → no allocation, no score change.
  - `span_losses(peptide, 0..s)` / `span_losses(peptide, s..n)` collect which losses the
    prefix and suffix fragments can carry (peptide-specific: a loss only applies to a
    fragment that spans the modified residue).
  - `ScoredSpectrum::loss_node_score(...)` adds the per-class loss-ion log-likelihood using
    the model's pooled per-class loss tables.
- **Why peptide-aware and not cached** — the normal split score is mass-indexed and cached
  (`cached_split_score`), but the loss contribution is `−L/z`-shifted and **peptide-specific**
  (depends on which residue carries the mod), so it cannot flow through the mass-keyed cache;
  it is added per split outside the cache. The hot path for standard peptides is untouched.
- **Training** — `ion_match_facts` derives loss facts from the intact-ion vocabulary × the
  declared losses; the rank-distribution estimator was patched to pick up loss keys that are
  absent from the no-loss `frag_off_table`. Loss keys live in the rank-dist table, not the
  frag-off table.

Net: glyco/neutral-loss support is a **per-class additive ion family**, gated three ways,
trained from the same own-data corpus, inert for non-modified peptides.

---

## 2. MaxSBM (MCP, June 2026) — the relevant prior art

Lennartsson, Kyriakidou, Nielsen, Olsen, Cox, Hendriks, *"Improved Peptide Search for
Identification of SUMO and Sequence-Based Modifiers, in MaxSBM"*, Mol. Cell. Proteomics
25(6):101589, 2026. DOI 10.1016/j.mcpro.2026.101589 (bioRxiv 2025.08.27.672604).

**The problem they name.** Standard engines model a PTM as a single static mass shift on a
residue. That is wrong for **sequence-based modifiers (SBMs)** — protein-derived modifiers
(SUMO2/3, ubiquitin) that leave a **multi-residue peptide remnant** on the lysine after
digestion. In HCD that remnant **itself fragments**, scattering the modifier's mass across
many product ions instead of sitting cleanly on the backbone, so a fixed-delta model
mis-scores the spectrum.

**Their method (inside MaxQuant/Andromeda).** Add two modifier-aware ion series to the
candidate's predicted spectrum:
- **d-ions ("diagnostic"):** product ions from the **cleaved-off modifier** sequence
  (for SUMO2/3, ~d2–d8/d9, plus a few "double-fragmented" internal d-ions for subsequences
  like FQ/QQ/FQQ). These are characteristic of the modifier *class*.
- **p-ions ("partial-remnant"):** backbone ions where the peptide **retains a partial
  modifier** — a partial neutral-loss ladder of the remnant. Optimum was a *small* set
  (p2, p3, p7), not the full ladder.

These extra theoretical ions are matched by Andromeda's existing probability score against a
reversed-decoy DB; an optional Percolator stage rescatters on top.

**Their results / caveats.** SUMO sites +17.7–25.5%, median score +9%; **+70% with
Percolator**. But: HCD DDA legacy instruments only (no Astral, no DIA), and **no entrapment /
no ground truth** — the authors explicitly flag that Percolator "could amplify false
positives" and that endogenous SUMO data has no objective validation. Ubiquitin (small LysC
remnant) gained only ~2% — the method only pays when the remnant is large enough to fragment
informatively, and **each modifier needs a hand-tuned d/p ion set**.

---

## 3. Why MaxSBM maps almost directly onto §1

andes's `loss_class` abstraction already expresses "this modification has its own family of
extra ions that appear under the right activation." MaxSBM's two series are exactly that:

| MaxSBM concept | andes mechanism it maps to |
|---|---|
| **p-ions** (backbone retaining a partial remnant) | a **per-class neutral-loss series** — the *exact* shape of `span_losses` + `loss_node_score`. A p-ion ladder = a set of `neutral_losses` on a `loss_class`, scored additively for fragments spanning the modified K. |
| **d-ions** (ions from the cleaved-off modifier itself) | a **new ion family** keyed by a `loss_class` whose masses come from the *modifier* sequence, not the peptide backbone — a small extension of the `IonType`/vocabulary, scored with its own per-class table. Unlike p-ions, d-ions are **peptide-position-independent** (they are the modifier's own b/y series), so they can be cached more aggressively than p-ions. |
| **per-modifier d/p tuning** (SUMO ≠ ubiquitin) | a **PTM-class profile keyed by (modifier, protease)** — consistent with how andes already keys models by (activation, instrument, enzyme, protocol). The remnant sequence (and thus the d/p masses) is protease-dependent, so the key must include the protease. |

So extending andes from "glycan/labile neutral losses" to "fragmenting protein-remnant
modifiers (SUMO, ubiquitin, and large glycans)" is mostly: (a) let a `loss_class` also carry
a **modifier-derived ion series** (d-ions), and (b) curate the d/p mass sets per
(modifier, protease). The additive, three-way-gated, inert-for-standard-peptides structure
of §1 is reused unchanged.

This is also the right frame for **large/branched glycans**: the glycan oxonium/Y-ion
"diagnostic" ions are andes's analogue of d-ions, and stepwise glycan losses are the p-ion
analogue. The same machinery serves glyco and SBM.

---

## 4. Proposed extension path (design, not yet built)

1. **d-ion family.** Allow a `loss_class` to reference a short **modifier ion series**
   (the remnant's own b/y masses for SBMs; the oxonium/Y set for glycans). Score it with a
   dedicated per-class table, like the existing loss tables. Keep the series **small** —
   MaxSBM found p2/p3/p7 and d2–d8 sufficient; the full ladder over-predicts peaks.
2. **Per-(modifier, protease) profile.** A keyed config object carrying the remnant
   sequence, the d/p mass sets, and the activation gate. Curated per modifier (SUMO2/3,
   ubiquitin-GG, large glycans), exactly as MaxSBM hand-tunes each.
3. **`--refine` integration.** Today `--refine` discovers PTMs as delta masses. Add an
   SBM/glyco-remnant tier so a refinement pass can apply a remnant profile (with its d/p
   ions) on top of the base PSM — the same Pass-2 cascade already used for the entrapment-
   validated +7.4% refine gain, but with loss-aware scoring instead of plain delta masses.

---

## 5. The integration that survives entrapment-FDP (the important part)

MaxSBM's headline **+70% is Percolator-amplified, target-decoy, on endogenous data with no
ground truth** — the same class of "add more theoretical peaks → match more" change that
inflates IDs unless validated against truth. andes's own refine "+6%" had exactly this trap
(an upper bound until grouped/entrapment validation). So:

- **Gate any d/p-ion change at 1% TRUE entrapment-FDP**, never target-decoy. If extra ions
  only raise the count by diluting decoys, entrapment-FDP exposes it.
- **Prefer the additive-feature integration over expanding the match-set.** The safest
  pattern in this codebase's history is: *additive PIN/GBDT features help; modifying the
  existing match-set / score regresses.* So rather than (or in addition to) injecting d/p
  ions into the matched-peak set, expose **"d-ion intensity fraction"** and **"p-ion ladder
  completeness"** as **PIN features** (and native-GBDT features). The FDR model then sees the
  diagnostic signal directly, which both discriminates SBM/glyco PSMs and sidesteps the
  match-set-inflation risk.
- **FDR stays Percolator** (production) / the single-pass native GBDT (labeled fallback).
  d/p-ions are an andes *scoring/feature* change; they do not touch the FDR boundary.

**Bottom line:** MaxSBM validates that andes's `loss_class` design is the right substrate for
fragmenting modifiers, and points at a concrete extension (modifier d-ions + small p-ion
ladders, keyed per modifier/protease). The discipline is unchanged: add it as additive,
gated, own-data-trained scoring + diagnostic PIN features, and prove every ID gain at true
entrapment-FDP before banking it.
