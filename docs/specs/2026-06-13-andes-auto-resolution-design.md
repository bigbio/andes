# Auto-detect instrument; remove `--instrument`; make fragmentation/tolerance MGF-only

**Date:** 2026-06-13
**Status:** Design — approved (Option B), pending implementation plan
**Scope:** `andes` search CLI + `crates/search` + `crates/scoring` + docs

## Goal

Remove `--instrument` from the search command and make `--fragmentation` plus a
new `--fragment-tol-*` override **MGF-only extended parameters**. The driving
insight: **MGF is the only systematically metadata-less input format.** mzML,
native Thermo `.raw`, and Bruker `.d` all carry the activation method and
analyzer in their metadata, so for those formats both the instrument resolution
class *and* the activation method are already auto-detected with no flags.

Therefore:

- **Metadata formats (mzML / `.raw` / `.d`): zero-config.** Instrument and
  activation are read from metadata. `--fragmentation` / `--fragment-tol-*` are
  ignored for model selection (auto-detection wins).
- **MGF (and rare metadata-poor mzML): the user describes it.** Activation comes
  from `--fragmentation`; resolution/tolerance from `--fragment-tol-ppm` /
  `--fragment-tol-da`. If neither is given, a **safe warned default** applies.
- **No data-driven inference pre-pass.** (Option C, the fragment-mass-error
  classifier, was considered and dropped — too complex/risky for the payoff once
  MGF is the only metadata-less case and the extended params cover it.)

Andes is a new tool with no released users, so there are **no back-compat
shims, aliases, or deprecation notices** — removed flags are simply gone.

## Background (why the flags existed)

`--instrument` had two jobs: (1) it was part of the `(activation, instrument,
protocol)` model-store lookup key; (2) the selected model's instrument drives
the fragment-matching tolerance — `is_high_resolution()` picks 20 ppm (high-res)
vs 0.5 Da (low-res) at several sites (`match_engine::matched_peak_keys`,
`compute_psm_features`, `coisolation`, `scored_spectrum`).

Crucially the tolerance is read from the **selected model's** instrument, not the
CLI flag directly — so determining the instrument once (from metadata, or from
an explicit override) makes the tolerance follow automatically. The model store
already accepts `instrument: None` (`build_selection_key` defaults it).

Auto-detection already runs today for mzML/`.raw`/`.d` when `--fragmentation
auto` (the default). The flags only ever bite on MGF (no metadata), on detection
returning `None`, or when `--fragmentation` is set non-auto (which currently
*disables* detection entirely — a wart this design removes).

## Decisions

| # | Decision |
|---|----------|
| A | `--instrument` is **hard-removed** from the search command. No alias. |
| B | `--fragmentation` becomes a **hidden, MGF-only** extended param; it sets the activation method for metadata-less input only. For mzML/`.raw`/`.d` it no longer disables detection — metadata always wins. |
| C | New hidden **`--fragment-tol-ppm` / `--fragment-tol-da`** (mutually exclusive) MGF-only override sets the resolution class + matching tolerance. |
| D | **No inference pre-pass.** Metadata-less input is handled by the extended params or a warned default. |
| E | Metadata-less input with no override → **CID / LowRes / 0.5 Da** default + a loud warning telling the user to set `--fragmentation` / `--fragment-tol-*`. LowRes is the non-catastrophic choice (a tight 20 ppm window on genuinely low-res data matches almost nothing; a 0.5 Da window on high-res data is merely noisier). |
| F | README + DOCS must state explicitly that `--fragmentation` and `--fragment-tol-*` are **MGF-only** (no effect on metadata-bearing formats). |

## Design

### 1. CLI surface (`crates/andes/src/bin/andes.rs`)

- **Remove** the `instrument: Instrument` field from the search `Cli`. Remove
  now-dead helpers (`parse_instrument`, `cli_instrument_to_instrument_type`,
  and the `Instrument` enum if unused on the search path). Old invocations
  passing `--instrument` fail with clap's standard "unexpected argument" error.
- **`--fragmentation`**: add `hide = true`; help text states it applies to MGF
  input only. It sets the activation method for metadata-less input; for
  metadata-bearing input the detected activation always wins (no longer
  disables detection).
- **New hidden** `--fragment-tol-ppm <f64>` / `--fragment-tol-da <f64>`,
  mutually exclusive (clap `conflicts_with`); help states MGF-only. Reuses the
  names already present on the `train-from-msnet` subcommand.
- The `train` subcommand's separate `--instrument` (`String`, for tagging
  trained models) is **untouched**.
- Update the `--param-file` doc comment, which references the removed flag.

### 2. Resolution / activation resolution order (in `main`)

The instrument + activation fed to `load_param_from_store` resolve in this
precedence:

1. **Metadata detection** (mzML peek / `.raw` vendor / `.d` = CID+TimsTOF) →
   if it yields an instrument, **use it** (and the detected activation). This is
   the unchanged hot path; the extended params are ignored here.
2. **Metadata-less input** (MGF, or detection returned `None`):
   - resolution/tolerance: `--fragment-tol-ppm` → high-res (`QExactive`);
     `--fragment-tol-da` → low-res (`LowRes`); neither → `LowRes` (decision E).
   - activation: `--fragmentation` if set; else the class-consistent default
     (`LowRes → CID`, `HighRes → HCD`).
   - emit the decision-E warning when both overrides are absent.

**Key interaction:** `build_selection_key` normalizes `(HCD, LowRes) → (HCD,
QExactive)`. Pairing the low-res default with `CID` activation (not `HCD`) keeps
the resolved class low-res. An explicit `--fragmentation` is still respected.

### 3. Effective-tolerance centralization (refactor)

Today 20 ppm / 0.5 Da is recomputed independently from `is_high_resolution()` at
~4 sites (`match_engine.rs:996`, `:1069`, `coisolation.rs:123`,
`scored_spectrum.rs:1328`). Resolve the **effective fragment tolerance once**
into a single `Tolerance` value (= `--fragment-tol-*` override if present, else
the selected model's `is_high_resolution()`-derived default) and thread it to
those consumers.

- The override replaces the *matching* tolerance only; the model's trained
  `mme` (rank-table binning) is left untouched (overriding `mme` miscalibrates
  the model — the scoring `mme` must equal the training `mme`).
- **Regression guard:** with no override, the resolved value must equal *exactly*
  today's 20 ppm / 0.5 Da, so model selection and PIN features stay
  **byte-identical** on the Astral / TMT / UPS1 benchmarks. Gated by existing
  golden tests.

### 4. Documentation (decision F)

- **README:** the options table must label `--fragmentation` and
  `--fragment-tol-ppm` / `--fragment-tol-da` as **MGF-only** extended
  parameters, with a one-line note that mzML/`.raw`/`.d` are fully auto-detected
  and ignore them.
- **DOCS.md:** update §1 (remove `--instrument`; mark the two params MGF-only),
  §4 Auto-detection (state metadata formats are zero-config; MGF is described by
  the extended params or gets the warned default), and any `--instrument`
  references.

## Behavior matrix

| Input | Today | After |
|---|---|---|
| mzML / `.raw` / `.d` with metadata | metadata detect | **unchanged**, zero-config |
| mzML with detect = `None` (metadata-poor) | falls back to `--instrument` default | extended params, else warned `CID/LowRes` default |
| MGF + `--fragment-tol-*` and/or `--fragmentation` | flag-based | uses the extended params |
| MGF + nothing | `--instrument` default → `hcd_qexactive` | warned **`cid_lowres`** default (decision E) |

## Error handling & edge cases

- **Both `--fragment-tol-ppm` and `--fragment-tol-da`** → clap rejects at parse.
- **Metadata-bearing input + extended params** → params ignored for model
  selection (detection wins); warn **once** that they had no effect (not
  per-file, to avoid noise on mixed batches).
- **Metadata-less + no override** → decision-E default + warning (not an error;
  unattended MGF runs still proceed).

## Testing

- **Unit:** effective-tolerance resolver — override present vs absent yields the
  expected `Tolerance`; absent reproduces 20 ppm / 0.5 Da exactly.
- **Integration:** MGF + `--fragment-tol-da 0.6` selects the low-res family,
  applies 0.6 Da matching, leaves model `mme` untouched. MGF + `--fragment-tol-ppm`
  selects the high-res family. MGF + nothing → warns + selects `cid_lowres`.
- **Integration:** MGF + `--fragmentation ETD` selects an ETD model.
- **Regression:** mzML/`.raw`/`.d` select byte-identical models vs today and
  ignore the extended params (protects Astral/TMT/UPS1 benchmarks).
- Update/remove existing CLI tests that pass `--instrument`; update the CI
  benchmark script + README/DOCS references.

## Out of scope

- A data-driven fragment-mass-error inference pre-pass (Option C) — dropped.
- Inferring the activation method from data — `--fragmentation` covers MGF;
  metadata formats auto-detect.
- Overriding the model's trained `mme` binning.
- Dropping MGF support entirely (Option A) — rejected; the benchmark infra and
  README depend on MGF.
