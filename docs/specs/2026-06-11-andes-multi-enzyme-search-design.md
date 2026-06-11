# Andes multi-enzyme search support — design (Sub-project A)

Date: 2026-06-11
Status: approved (brainstorming) → ready for implementation plan
Scope: **Sub-project A only** — the andes search engine. Sub-project B (harvest per-file
enzyme routing + re-harvest) is described at the end for context but specced separately.

## Problem

Andes ships scoring models for many proteases (Trypsin, Lys-C, Arg-C, Asp-N, Glu-C, Lys-N,
ALP/aLP, NoCleavage; plus Chymotrypsin and NonSpecific in the cleavage engine), but the
**search path can only ever digest with Trypsin**. Two hardcodes block everything else:

- `params.enzyme` is fixed to `Enzyme::Trypsin` in `default_tryptic` (`crates/search/src/search_params.rs:97`) and never reassigned in the search/run path.
- The model-selection key hardcodes the enzyme to `"Trypsin"` in `build_selection_key` (`crates/andes/src/bin/andes.rs:3311`).

There is **no `--enzyme` search flag** (`SearchArgs`, `crates/andes/src/bin/andes.rs:109-284`).
Forcing a non-tryptic model with `--model <id>` only swaps the *scoring* model; candidate
peptides are still generated tryptically — a mismatched, unusable combination. The N-terminal
cutters (Asp-N, Lys-N) are the most affected. The code self-documents the gap:
`resolve_bundled_param` comment *"Only Tryp is supported as the enzyme component for now"*
(`andes.rs:2948`).

Additionally, real-world digests are frequently **combinations** (e.g. Trypsin+Lys-C is a
standard co-digestion; the PXD000900 multi-protease study contains GluC-LysC, ArgC-AspN,
LysC-elastase, etc.). Andes has no way to express a combined digest.

## Goals

1. A search-time `--enzyme` flag accepting **one or several** enzymes.
2. Combined digestion = **union of cut sites** in a single search (one candidate set, one score).
3. Model selection by **canonical cut-set signature**, so equivalent digests collapse
   (Trypsin ≡ Trypsin+Lys-C → reuse the trypsin model) and only genuinely-orthogonal digests
   (Trypsin+Glu-C → cut after K,R,E) resolve to a distinct model.
4. Full back-compat: existing invocations (no `--enzyme`) behave exactly as today (trypsin).

Non-goals (Sub-project B / later): per-file enzyme detection in the harvest pipeline,
MSFragger combined-cleavage config, minting orthogonal-combo *models*, re-harvesting PXD000900.

## Cleavage semantics (decided)

- **Union of cut sites, one search.** A peptide bond is cleavable if ANY enzyme in the set
  cleaves it. `Trypsin+Lys-C` → cut after K,R; `Trypsin+Glu-C` → cut after K,R,E. This matches
  combined/co-digestion as implemented by MaxQuant / MSFragger "multiple enzymes".

## Model identity (decided)

- **Canonical cut-set, collapse equivalents.** A digest is identified by its effective cut-set
  signature, not by enzyme names. Equivalent sets share a model; distinct cut-sets get distinct
  models. "Lys-C is usually combined with Trypsin" falls out automatically — their signatures
  are identical.

## Design

### 1. `CleavageSpec` (new type, `crates/model/src/enzyme.rs`)

Represents the union of cleavage rules from a **set** of `Enzyme`s.

- Constructed from `&[Enzyme]` (or an `EnzymeSet`). Internally unions each enzyme's
  `rules()` output: cut-after residue set, cut-before residue set, any restrict/not-before
  residues, and the `universal` flag (NonSpecific/aLP cleave everywhere; NoCleavage never).
- Exposes the interface the search already consumes:
  - `is_cleavable_after(prev_residue, next_residue) -> bool`
  - `is_cleavable_before(prev_residue, next_residue) -> bool`
- Exposes `signature() -> CutSetSignature`: a canonical, order-independent encoding of the
  effective cut rules (normalized so that, e.g., union ordering and duplicate enzymes don't
  change it). Equality of signatures defines digest equivalence.

The existing `Enzyme` enum stays as the per-enzyme rule source. `CleavageSpec` is the only new
abstraction. Because `compute_cleavage_positions` (`crates/search/src/candidate_gen.rs:120,430-446`)
and `sa_walk.rs:211-223` already call `is_cleavable_after/before` generically, the digestion
loop changes only the *type* it reads (`CleavageSpec` instead of a single `Enzyme`) — no
algorithmic change.

### 2. Search params (`crates/search/src/search_params.rs`)

- `SearchParams` carries a `CleavageSpec` (or an `EnzymeSet` from which the spec is derived)
  in place of the single `enzyme: Enzyme` field.
- `default_tryptic` builds the spec from `[Enzyme::Trypsin]` — preserving today's default.
- NTT / `num_tolerable_termini` logic (`andes.rs:1173-1177`) is unchanged; it already operates
  generically against the cleavage interface.

### 3. CLI (`crates/andes/src/bin/andes.rs`, `SearchArgs`)

- Add `--enzyme <LIST>`: one or more enzyme names, `+` or `,`-separated
  (`--enzyme trypsin`, `--enzyme trypsin+lysc`, `--enzyme gluc,lysc`).
- Parse each token via the existing `Enzyme::from_name` (case-insensitive, aliases incl. `aLP`,
  `NoCleavage`, `Chymotrypsin`; `elastase` added as an alias to `NonSpecific` — see below).
  Build a `CleavageSpec` from the resulting set.
- Default when omitted: `[Trypsin]` (back-compat).
- Unknown name → hard error (consistent with the train path's
  `from_name(...).ok_or_else(... "unknown --enzyme '{enz}'")`, `andes.rs:2530-2531`).
- `--model <id>` continues to force an exact model, bypassing selection.

### 4. Model selection (`crates/andes/src/bin/andes.rs:3311`, `build_selection_key`)

- Replace the hardcoded `enzyme: "Trypsin".to_string()` with the **signature** of the requested
  `CleavageSpec`.
- `load_param_from_store` / the `select()` ladder (`crates/model-train/src/select.rs:83`)
  matches a candidate model by computing the signature of that model's stored `enzyme` value and
  comparing to the requested signature. Single-enzyme store values (`Trypsin`, `LysC`, `AspN`,
  …) keep working unchanged; the collapse (Trypsin+LysC → trypsin model) is automatic because
  signatures match.
- No store rewrite required for A. (Persisting orthogonal-combo enzyme strings is Sub-project B.)

### 5. Enzyme vocabulary

In scope for first pass: Trypsin, Lys-C, Arg-C, Asp-N, Glu-C, Lys-N, ALP/aLP, NoCleavage,
NonSpecific, **Chymotrypsin** (rules already exist: cut after F,Y,W,L — `enzyme.rs:33`).
**Elastase**: no dedicated rules and very low specificity → aliased to `NonSpecific` for the
first pass (it will not get a hand-tuned rule set or a dedicated high-quality model). Revisit
if a real elastase model is wanted.

## Testing (TDD)

Unit (`crates/model` / `crates/search`):
- `signature()` equivalence: `Trypsin` == `Trypsin+LysC`; `Trypsin+GluC` distinct from both;
  N-terminal cutters Asp-N (before D) and Lys-N (before K) produce the expected signatures;
  order-independence (`gluc,lysc` == `lysc,gluc`); duplicate enzyme is idempotent.
- `CleavageSpec` union: `is_cleavable_after` true for K,R **and** E under `Trypsin+GluC`;
  NoCleavage never cleaves; NonSpecific/aLP cleave everywhere.

Unit (`crates/andes`):
- CLI parse: `--enzyme trypsin+lysc` → 2-enzyme set; unknown token → error; omitted → trypsin.
- Selection: a requested combo whose signature collapses selects the single-enzyme model;
  a distinct signature with no matching model surfaces the expected behavior (error or
  documented fallback — pick one explicitly in the plan).

End-to-end:
- A small fixture search with `--enzyme gluc` yields a measurably different candidate set
  (and different PSMs) than the trypsin default, proving digestion actually changed.
- A no-`--enzyme` run is byte-identical to the pre-change trypsin baseline (back-compat).

## Risks / open implementation choices (resolve in the plan)

- **Combo with no matching model**: when a requested signature matches no store model (e.g. a
  novel orthogonal combo), decide between a hard error vs. a documented nearest-signature
  fallback. Default recommendation: hard error in A (keeps behavior honest); add fallback only
  if needed.
- **Signature canonical form**: must be stable across enzyme ordering, duplicates, and the
  universal/NoCleavage edge cases. Pin it with the `name_round_trips`-style test discipline
  already used in `enzyme.rs:299-310`.
- **Branch hygiene**: the working tree currently sits on `feat/train-data-type-override` with
  unrelated WIP. Implement Sub-project A on its own branch off a clean base.

## Sub-project B (context only — separate spec)

Harvest pipeline (`harvest_model.sh` + python on Codon) routes enzyme **per dataset** (one
route per MS2 analyzer), which mislabels multi-protease datasets — confirmed on PXD000900 where
GluC/AspN/ArgC/elastase files were all routed to `lysc` (wrong-cleavage PSM counts of 2–264).
Those 20 flats are quarantined. Sub-project B will: detect enzyme(s) **per file** (filenames
carry them: `GluC-LysC`, `ArgC-AspN`…), configure the MSFragger label-search to the combined
cleavage, route each flat to the slug matching its cut-set signature (reusing A's signature
primitive), and re-harvest PXD000900 / re-examine PXD010154 to rebuild the now-empty enzyme
corpora.
