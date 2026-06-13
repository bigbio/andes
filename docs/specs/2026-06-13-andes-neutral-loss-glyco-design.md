# Neutral-loss-aware scoring for labile mods (glycopeptides, e.g. Unimod 393)

**Date:** 2026-06-13
**Status:** Design — approved, pending implementation plan
**Scope:** `crates/model`, `crates/scoring`, `crates/model-train`, `crates/andes`, docs
**Spec type:** unified design for a phased feature (SP1→SP4), one branch + milestone commits

## Goal

Let andes score **labile-modification glycopeptides** well by modeling the
**neutral losses** of a modification. Today fragment prediction is "canonical
b/y only, no neutral losses" ([fragment_ions.rs:3](../../crates/scoring/src/scoring/fragment_ions.rs)),
and the node-score `IonType` is only `Prefix`/`Suffix`
([param_model.rs:495](../../crates/scoring/src/param_model.rs)) — so for a
glycopeptide whose glycan falls off in MS/MS, the dominant deglycosylated
backbone ions go unmatched and scoring is degraded.

Concretely, andes must search **Unimod 393** (Glucosylgalactosyl,
glucosylgalactosyl-hydroxylysine, +340.100562 on K) and score its
characteristic stepwise sugar losses.

## Literature basis (informs the design)

- Glycopeptide HCD/CID fragmentation is a **stepwise Y-ion ladder** (peptide +
  progressively smaller retained glycan) plus oxonium markers — not a single
  clean loss. Glyco engines model a **set** of fragment masses: MSFragger-Glyco
  / -Labile use a user-specified Y-ion list; O-Pair/MetaMorpheus and pGlyco
  compute Y-ion ladders. General engines (Comet) offer a **single** per-mod
  fragment neutral loss.
- **Unimod-393-specific:** the observed backbone losses are **−162.05 Da (one
  Hex)** and **−324.11 Da (two Hex)** — the sugar ladder that **retains the
  hydroxylysine**. Note `324 ≠ 340`: the Unimod "neutral loss = 340.100562"
  field includes the hydroxyl O, but real backbone fragments lose only the
  sugars (keep +O). So loss masses are **user-specified sugar units**, NOT
  auto-derived from the Unimod composition. Collagen glycopeptides also notably
  *retain* glycan on many backbone ions, so predicting both intact and loss
  ions matters. (JASMS 10.1007/s13361-013-0624-y; PMC4235766.)

## Decisions (settled in brainstorming)

| # | Decision |
|---|----------|
| A | **Full route**: a dedicated loss `IonType` in the GF node-score model + a retrained `glyco` model. (Not feature-only.) |
| B | Glyco training data is available; retraining a `glyco` model is acceptable now. |
| C | **Multiple, user-specified losses per mod** (a list), e.g. `162.0528;324.1056`. Predict the stepwise ladder + the intact ion. **No oxonium-marker modeling in v1.** |
| D | **Pooled loss key**: `IonType` gains one 2-state discriminant (intact vs lost). ALL loss ions share one learned distribution family per `(partition, base series, charge)`; the loss mass only shifts m/z, never the table key. (Loss masses are per-search, so they cannot be stable model keys.) |
| E | One unified design; **phased implementation** SP1→SP4 as milestone commits. |
| F | **Configuration must be easy** (mods.txt is hand-edited): clear grammar, validation errors, common-loss reference, copy-paste glyco example, shipped template. Optional mod attributes use an **extensible `key=value` tail** (not positional fields). |
| G | **Unimod accession interop**: parse `accession=UNIMOD:393` into `Modification.accession` and emit it in a NEW **additive** PIN/TSV column (standard CURIE). Existing mass-delta `Peptide` column unchanged; ProForma notation deferred. |

## Phase SP1 — Loss-aware representation

### Data model — `crates/model/src/modification.rs`
Add `neutral_losses: Vec<f64>` to `Modification` (default empty ⇒ no loss ⇒
current behavior). Update the 4 in-file constructors/tests.

### Configuration & grammar (decision F — first-class)
`from_mods_txt_line` currently does `splitn(5, ',')` and requires exactly 5
fields; the 5th (name) absorbs trailing commas. Change to **5 positional core
fields followed by an optional, extensible `key=value` attribute tail** (each
attribute is its own comma-separated field):

```
<mass>,<residue>,<fix|opt>,<location>,<name>[,<key>=<value>]...
```

v1 attribute keys:
- `loss=<m1;m2;…>` — semicolon-separated neutral-loss masses in Da.
- `accession=<CURIE>` — modification CV accession, conventionally `UNIMOD:393`
  (generic, so `RESID:AA0153` / `MOD:00000` also parse). Stored in
  `Modification.accession`. See *Interoperability* below.

Rationale for labeled over positional: losses and accession are the first of
several optional mod attributes; a `key=value` tail is order-independent,
self-documenting, and extensible without further grammar changes.

- **Backward-compatible:** a 5-field line parses exactly as today (no attrs).
  Parsing: take the first 5 comma fields positionally; treat every subsequent
  comma field as `key=value`.
- **Caveat (documented):** when any attribute is present, the **name must not
  contain a comma** (it would be misread as an attribute field). Mod names with
  commas are pathological; documented as a v1 limitation. Pure 5-field lines
  retain comma-in-name behavior.
- **Validation / errors:** unknown key → `ModParseError::UnknownModAttr { key }`;
  malformed `loss` value → `BadNeutralLoss { value }` (each parses as `f64`,
  `> 0`, `< 2000`); malformed attribute (no `=`) → `BadModAttr { field }`. All
  messages name the offending text.
- **Ergonomics deliverables (DOCS §2):**
  - A **common-loss reference table**: Hex `162.0528`, HexNAc `203.0794`,
    NeuAc `291.0954`, phospho (H₃PO₄) `97.9769`, sulfo `79.9568`.
  - A **copy-paste Unimod-393 line**:
    `340.100562,K,opt,any,Glucosylgalactosyl,loss=162.0528;324.1056,accession=UNIMOD:393`
  - A complete **glyco `mods.txt` template** shipped alongside the existing
    examples, plus a one-paragraph "how losses are scored" note.
  - Explicit statement that attributes are **orthogonal** to `fix|opt`,
    `location`, and `NumMods` — they describe *fragment* behavior / identity,
    not residue placement; they work with any of those settings unchanged.

### Interoperability — Unimod accession (decision: capture + additive output)
Today mods reach the PIN/TSV output only as **inline mass deltas** in the
`Peptide` column (`Peptide`'s `Display`, [peptide.rs:105](../../crates/model/src/peptide.rs);
emitted at [pin.rs:514](../../crates/output/src/pin.rs)), so downstream tools map
mods by mass — ambiguous for isobaric mods. To improve quantms / SDRF / PSI
interoperability:
- Parse `accession=` into the existing `Modification.accession` field (currently
  always `None`).
- **Emit it in a NEW additive column** in PIN and TSV — e.g. `Modifications`
  listing `pos:CURIE` entries (`6:UNIMOD:393`) in standard CV form. **Additive
  only** — the existing mass-delta `Peptide` column is unchanged, so current
  quantms scripts keep working and accession-aware tools gain unambiguous
  mapping.
- ProForma notation (`K[UNIMOD:393]`) in the `Peptide` column is **deferred**
  (out of scope): it would modify an existing column, risking parser breakage
  and violating the additive-only rule.

### IonType — `crates/scoring/src/param_model.rs`
Add a 2-state loss discriminant to `Prefix`/`Suffix` (the pooled key, decision
D). Preferred shape: a `loss: bool` field on each variant (default `false`).
Touches: every `match` on `IonType`, `Hash`/`Eq`/ordering, the
`partition_ion_types_cache`, and serialization (SP2). Intact ions keep
`loss: false` ⇒ unchanged keys ⇒ existing tables unaffected.

### Fragment prediction — `crates/scoring/src/scoring/fragment_ions.rs`
For each predicted b/y ion at charge `z` whose span includes a residue carrying
`neutral_losses`, additionally emit one loss-shifted ion per declared loss `L`
at `mz_intact − L/z` with the same series/charge but `loss: true`, alongside the
intact ion. **v1 simplification:** if a fragment spans multiple loss-bearing
residues, emit each residue's losses independently — no cross-products (the
common case is one glyco site per peptide). **Inert guard:** a peptide with no
loss-bearing residue predicts zero loss ions ⇒ byte-identical to today.

## Phase SP2 — Model format + scoring

- **Store/Param:** the rank/error/existence/ion-existence tables already key on
  `IonType`; with the `loss` flag they gain pooled `loss: true` entries.
  Extend the **parquet store** IonType encoding to carry the flag.
  **Backward-compatible:** models with no `loss: true` entries load unchanged;
  a loss ion with no table scores as **absent** (the missing-ion path), never
  panics. The legacy `.param` reader never produces `loss: true`.
- **Scorer — `rank_scorer.rs`:** node-score lookup for a `loss: true` ion uses
  its pooled table; absent ⇒ missing-ion score. Matched-ion counting / DP
  include loss ions only when the peptide declares losses AND the model has
  loss tables.
- **Byte-identical guard:** standard searches (no loss mods) predict no loss
  ions; non-glyco models have no loss tables ⇒ existing benchmarks unchanged.

## Phase SP3 — Training

- **Estimator — `crates/model-train/src/estimate.rs`:** when accumulating from
  confident glyco PSMs (whose `mods.txt` declares losses), observe the predicted
  `loss: true` ions and accumulate their rank/error/existence stats into the
  **pooled** loss distributions — same machinery as intact ions. No new
  estimator math; just the extra IonType key flows through.
- Produces a **`glyco` model** (existing catalog slug
  [catalog.rs:95](../../crates/model-train/src/catalog.rs)) with loss tables,
  trained on the available glyco corpus. Wires into the model store like any
  trained model.

## Phase SP4 — Validation (benchmark-gated)

- **Glyco dataset** (the available corpus): search with the `glyco` model +
  glyco `mods.txt` (losses declared), FDP/entrapment-controlled → demonstrate
  **PSM gain on glycopeptides** vs the no-loss baseline.
- **Standard 3** (Astral / UPS1 / a05058 TMT): **byte-identical** (no loss mods)
  — the regression guard, runnable as the existing golden/parity tests + the VM
  benchmark.
- Per project rules, the scoring change is **gated on the VM benchmark**; the
  in-repo deliverable is implemented + unit-tested, handed off for the VM run.
  Do not claim a benchmark pass without it.

## Cross-cutting guarantees

- **Additive only:** loss tables and the `loss` flag are purely additive;
  existing ion types/tables/PIN columns are untouched (respects "additive PIN
  features only / never regress").
- **Inert when unused:** no loss mod ⇒ no loss ions ⇒ byte-identical; non-glyco
  model ⇒ loss ions score absent.
- **Independence-aligned:** the loss `IonType` is andes's own design, trained on
  andes's own glyco data — consistent with the MS-GF+ independence program.

## Testing

- **Unit (grammar):** `loss=` parses a `;`-list and `accession=` a CURIE;
  5-field lines unchanged; unknown key → `UnknownModAttr`; malformed loss →
  `BadNeutralLoss`; out-of-range loss rejected; attribute without `=` →
  `BadModAttr`; comma-in-name caveat covered.
- **Unit/Integration (accession):** `accession=UNIMOD:393` reaches
  `Modification.accession`; the additive PIN/TSV `Modifications` column emits
  `pos:UNIMOD:393`; the existing `Peptide` column is byte-identical (additive
  guard); a mod with no `accession=` emits an empty entry.
- **Unit (prediction):** a loss-bearing residue emits intact + one ion per loss
  at `mz − L/z` with `loss:true`; no loss-bearing residue ⇒ zero loss ions
  (inert); multi-charge correct.
- **Unit (model):** `IonType{loss:true}` round-trips through the parquet store;
  a model without loss tables scores a loss ion as absent (no panic).
- **Unit (scorer):** pooled loss table looked up for loss ions; intact-ion
  scoring byte-identical to pre-change for a no-loss peptide.
- **Unit (training):** estimator accumulates loss-ion stats into the pooled
  distribution from a synthetic glyco PSM.
- **Integration:** a glyco peptide with declared losses yields loss-ion matches
  in the PIN/score path; a standard search produces byte-identical goldens.
- **Validation (VM):** glyco benchmark gain + standard-3 parity.

## Out of scope (v1)

- Oxonium / glycan marker ions as spectrum-level evidence (decision C).
- Per-loss-mass or per-tier distributions (decision D: pooled only).
- Cross-product losses across multiple glyco residues in one fragment.
- Unimod composition-string parsing / auto-deriving losses from the mod
  composition (losses are explicit, user-specified).
- Glycan-structure database search / localization (this is labile-mod scoring,
  not a glyco search engine).

## File-touch summary

- SP1: `crates/model/src/modification.rs` (struct + grammar + `accession=`), `crates/scoring/src/param_model.rs` (IonType), `crates/scoring/src/scoring/fragment_ions.rs`, `crates/output/src/{pin.rs,tsv.rs}` (additive `Modifications` accession column), `DOCS.md` §2 + glyco template.
- SP2: `param_model.rs` (store schema), the model-store read/write, `crates/scoring/src/scoring/rank_scorer.rs`.
- SP3: `crates/model-train/src/estimate.rs`, `catalog.rs`.
- SP4: benchmark scripts/docs; `andes` search (no code change beyond SP1–SP3 — it just loads the glyco model + glyco mods.txt).
