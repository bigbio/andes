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
| D | **Per-class pooled loss key** (revised after the phospho/generality review): `IonType` gains a `loss_class: u8` discriminant (0 = intact; 1.. = a stable per-mod-class pool, e.g. Glyco=1, Phospho=2, Sulfo=3, Generic=255). Loss ions pool **within a class** (the −162/−324 glyco ladder shares one distribution — the data-efficiency win) but **not across classes** (glyco and phospho do not contaminate each other; phospho −98 is often dominant, glyco losses are modest ladder rungs). The specific loss *mass* still only shifts m/z, never the key. Rationale: no production engine pools losses globally (Comet/Andromeda/Mascot/MSFragger scope per-mod/residue/activation); widening `bool`→`u8` is ~free now but a breaking model-format change if deferred. |
| H | **Activation-gated loss prediction**: electron-based dissociation (ETD/EThcD) *preserves* labile mods, so loss-shifted ions are suppressed for ETD spectra; predicted only for collisional activation (CID/HCD/PQD). andes already detects activation (`SpecDataType.activation`). |
| E | One unified design; **phased implementation** SP1→SP4 as milestone commits. |
| F | **Configuration must be easy** (mods.txt is hand-edited): clear grammar, validation errors, common-loss reference, copy-paste glyco example, shipped template. Optional mod attributes use an **extensible `key=value` tail** (not positional fields). |
| G | **Unimod accession interop**: parse `accession=UNIMOD:393` into `Modification.accession` and emit it in a NEW **additive TSV-only** column (standard CURIE). PIN omits the column: PIN's `Proteins` is rest-of-line, so a trailing `Modifications` field is parsed by Percolator as a phantom extra protein — Percolator already receives modification positions via the `Peptide` inline mass-delta column; the CURIE accession column is emitted in andes's **TSV only** for quantms/SDRF interop. Existing mass-delta `Peptide` column unchanged; ProForma notation deferred. |

## Phase SP1 — Loss-aware representation

### Data model — `crates/model/src/modification.rs`
Add `neutral_losses: Vec<f64>` to `Modification` (default empty ⇒ no loss ⇒
current behavior) and `loss_class: u8` (default `0`; set from the `class=`
attribute via the fixed registry — Glyco=1/Phospho=2/Sulfo=3/Generic=255). A
`LossClass` name→id registry lives in `model` (a small `const`/`fn`).
**(Note: Task 1 already shipped `neutral_losses` + `loss=`/`accession=`; the
`class=` attribute + `loss_class` field are an additive follow-up — Task 1b.)**

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
- `class=<name>` — loss-class pool selecting which trained distribution the loss
  ions score against. Names map to stable `u8` ids via a fixed registry:
  `glyco`=1, `phospho`=2, `sulfo`=3 (extensible); when `loss=` is present but
  `class=` is omitted, default `generic`=255. Mods in the same class share one
  pooled loss distribution; different classes are scored independently
  (decision D). Stored on `Modification.loss_class`.
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
    `340.100562,K,opt,any,Glucosylgalactosyl,loss=162.0528;324.1056,class=glyco,accession=UNIMOD:393`
  - A complete **glyco `mods.txt` template** shipped alongside the existing
    examples, plus a one-paragraph "how losses are scored" note.
  - Explicit statement that attributes are **orthogonal** to `fix|opt`,
    `location`, and `NumMods` — they describe *fragment* behavior / identity,
    not residue placement; they work with any of those settings unchanged.

### Interoperability — Unimod accession (decision: capture + TSV-only additive column)
Today mods reach the PIN/TSV output only as **inline mass deltas** in the
`Peptide` column (`Peptide`'s `Display`, [peptide.rs:105](../../crates/model/src/peptide.rs);
emitted at [pin.rs:510](../../crates/output/src/pin.rs)), so downstream tools map
mods by mass — ambiguous for isobaric mods. To improve quantms / SDRF / PSI
interoperability:
- Parse `accession=` into the existing `Modification.accession` field (currently
  always `None`).
- **Emit it in a NEW additive column in TSV only** — `Modifications`,
  listing `pos:CURIE` entries (`6:UNIMOD:393`) in standard CV form. **Additive
  only** — the existing mass-delta `Peptide` column is unchanged, so current
  quantms scripts keep working and accession-aware tools gain unambiguous
  mapping.
- **PIN does NOT get this column.** PIN's `Proteins` column is
  **rest-of-line**: the row writer emits one tab-separated accession per protein
  in a loop, and the line ends only with `writeln!` after the last protein.
  Appending any field after the proteins loop would be parsed by Percolator as
  an additional protein accession (a phantom protein), corrupting FDR
  computation. Percolator already receives modification positions via the
  `Peptide` inline mass-delta column, so no information is lost; the CURIE
  accession column is emitted in andes's **TSV only** for quantms/SDRF interop.
- ProForma notation (`K[UNIMOD:393]`) in the `Peptide` column is **deferred**
  (out of scope): it would modify an existing column, risking parser breakage
  and violating the additive-only rule.

### IonType — `crates/scoring/src/param_model.rs`
Add a `loss_class: u8` discriminant to `Prefix`/`Suffix` (the per-class pooled
key, decision D). `0` = intact; `1..` = a loss-class pool (Glyco=1, Phospho=2,
Sulfo=3, Generic=255). Touches: every `match` on `IonType`, `Hash`/`Eq`/ordering,
the `partition_ion_types_cache`, and serialization (SP2). Intact ions keep
`loss_class: 0` ⇒ unchanged keys ⇒ existing tables unaffected. Add accessors
`is_loss()` (`loss_class != 0`) and `loss_class()`.

### Fragment prediction — `crates/scoring/src/scoring/fragment_ions.rs`
For each predicted b/y ion at charge `z` whose span includes a residue carrying
`neutral_losses`, additionally emit one loss-shifted ion per declared loss `L`
at `mz_intact − L/z` with the same series/charge and `loss_class` = that mod's
class id, alongside the intact ion. **Activation gate (decision H):** loss ions
are predicted only for collisional activation (CID/HCD/PQD); for electron-based
methods (ETD) — which preserve labile mods — loss prediction is suppressed. The
caller supplies the activation (from `scorer.param().data_type.activation`); a
helper `ActivationMethod::predicts_neutral_losses()` returns `false` for ETD.
**v1 simplification:** if a fragment spans multiple loss-bearing residues, emit
each residue's losses independently — no cross-products (the common case is one
glyco site per peptide). **Inert guard:** a peptide with no loss-bearing residue
(or an ETD spectrum) predicts zero loss ions ⇒ byte-identical to today.

## Phase SP2 — Model format + scoring

- **Store/Param:** the rank/error/existence/ion-existence tables already key on
  `IonType`; with `loss_class` they gain per-class entries. Extend the **parquet
  store** IonType encoding to carry the `loss_class` byte (added to the flat
  ion-type encoding; reader defaults `0` when absent → back-compat).
  **Backward-compatible:** models with no `loss_class != 0` entries load
  unchanged; a loss ion whose `(partition, ion_type-with-class)` table is absent
  scores as **absent** (the missing-ion path), never panics. The legacy
  `.param` reader never produces `loss_class != 0`.
- **Scorer — `rank_scorer.rs`:** node-score lookup for a `loss_class != 0` ion
  uses its per-class pooled table; absent ⇒ missing-ion score. Matched-ion
  counting / DP include loss ions only when the peptide declares losses AND the
  model has the matching loss-class tables.
- **Byte-identical guard:** standard searches (no loss mods) predict no loss
  ions; non-glyco models have no loss tables ⇒ existing benchmarks unchanged.

## Phase SP3 — Training

- **Estimator — `crates/model-train/src/estimate.rs`:** when accumulating from
  confident glyco PSMs (whose `mods.txt` declares `loss=`/`class=glyco`),
  observe the predicted `loss_class=Glyco` ions and accumulate their
  rank/error/existence stats into the **Glyco-class** pooled distributions —
  same machinery as intact ions. No new estimator math; the extra IonType key
  flows through.
- Produces a **`glyco` model** (existing catalog slug
  [catalog.rs:95](../../crates/model-train/src/catalog.rs)) populating the
  **Glyco loss class only**, trained on the available glyco corpus. **v1 trains
  and validates the Glyco class only**; other classes (Phospho, Sulfo, …) are
  *configurable* but require their own training data + model and are explicit
  fast-follows (do not ship a glyco-trained table as a phospho solution).

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
- **Inert when unused:** no loss mod (or ETD spectrum) ⇒ no loss ions ⇒
  byte-identical; a model lacking the declared loss class ⇒ loss ions score
  absent.
- **Independence-aligned:** the loss `IonType` is andes's own design, trained on
  andes's own glyco data — consistent with the MS-GF+ independence program.

## Testing

- **Unit (grammar):** `loss=` parses a `;`-list and `accession=` a CURIE;
  5-field lines unchanged; unknown key → `UnknownModAttr`; malformed loss →
  `BadNeutralLoss`; out-of-range loss rejected; attribute without `=` →
  `BadModAttr`; comma-in-name caveat covered.
- **Unit/Integration (accession):** `accession=UNIMOD:393` reaches
  `Modification.accession`; the additive **TSV** `Modifications` column emits
  `pos:UNIMOD:393`; the existing `Peptide` column is byte-identical (additive
  guard); a mod with no `accession=` emits an empty entry. **PIN is unchanged**
  (the Modifications column is TSV-only; PIN's Proteins column is rest-of-line
  and cannot be followed by a trailing column without corrupting Percolator's
  protein parsing).
- **Unit (prediction):** a loss-bearing residue emits intact + one ion per loss
  at `mz − L/z` tagged with the mod's `loss_class`; no loss-bearing residue ⇒
  zero loss ions (inert); **ETD activation ⇒ zero loss ions (gate, decision H)**;
  multi-charge correct.
- **Unit (model):** `IonType{loss_class:1}` round-trips through the parquet
  store; a model without that class's tables scores a loss ion as absent (no
  panic).
- **Unit (scorer):** per-class pooled table looked up for loss ions; intact-ion
  scoring byte-identical to pre-change for a no-loss peptide.
- **Unit (training):** estimator accumulates loss-ion stats into the pooled
  distribution from a synthetic glyco PSM.
- **Integration:** a glyco peptide with declared losses yields loss-ion matches
  in the PIN/score path; a standard search produces byte-identical goldens.
- **Validation (VM):** glyco benchmark gain + standard-3 parity.

## Out of scope (v1)

- Oxonium / glycan marker ions as spectrum-level evidence (decision C).
- Per-loss-mass or per-tier distributions (decision D: per-**class** pooled, not
  per-mass).
- **Trained non-glyco classes:** Phospho/Sulfo/etc. are configurable
  (`class=phospho`) and predicted, but v1 does not train/validate them; a
  glyco-trained model has no Phospho loss tables ⇒ phospho loss ions score
  absent until a phospho model exists (fast-follow). Docs must warn: declare
  pY's loss as `−79.9663` (HPO₃), **not** `−97.9769`, since pY rarely loses
  H₃PO₄ (false-localization trap).
- Cross-product losses across multiple loss-bearing residues in one fragment.
- Unimod composition-string parsing / auto-deriving losses from the mod
  composition (losses are explicit, user-specified).
- Glycan-structure database search / localization (this is labile-mod scoring,
  not a glyco search engine).

## File-touch summary

- SP1: `crates/model/src/modification.rs` (struct + grammar + `accession=`), `crates/scoring/src/param_model.rs` (IonType), `crates/scoring/src/scoring/fragment_ions.rs`, `crates/output/src/tsv.rs` (additive `Modifications` accession column — **TSV only**; PIN unchanged because `Proteins` is rest-of-line), `DOCS.md` §2 + glyco template.
- SP2: `param_model.rs` (store schema), the model-store read/write, `crates/scoring/src/scoring/rank_scorer.rs`.
- SP3: `crates/model-train/src/estimate.rs`, `catalog.rs`.
- SP4: benchmark scripts/docs; `andes` search (no code change beyond SP1–SP3 — it just loads the glyco model + glyco mods.txt).
