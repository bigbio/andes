# Porting an externally trained model into andes

How a model trained elsewhere reaches a search, what exists today, and what has to be
built. The answer depends on the kind of model: a tree ensemble fits a slot that already
exists, a neural network goes through the ONNX predictor path. This page describes the code
as it is; "not built" marks the missing pieces.

## What andes loads

Scoring models live in one Parquet store, `resources/models/protocol=<P>/models.parquet`,
read by `crates/model-train/src/store/read.rs`. A model is keyed by `model_id`. Its manifest
row carries the selection tuple `(activation, instrument, enzyme, protocol)`, the model's own
fragment tolerance (`mme_val`, `mme_is_ppm`) and three nullable binary columns, each an
`AGBD` blob (`crates/model-train/src/store/schema.rs`):

| Column | `Param` field | Role | Output used |
|---|---|---|---|
| `gbdt_model_bytes` | `gbdt_peak_model` | per-peak signal/noise classifier | `predict_logit` (LLR term) |
| `frag_intensity_model_bytes` | `frag_intensity_model` | fragment-intensity regressor | `predict_value` (raw tree sum) |
| `rich_ion_model_bytes` | `rich_ion_model` | per-ion classifier | `predict_value` (raw tree sum) |

A null blob is legal, and a store written before a column existed reads as `None`
(`read.rs`, optional column reads), so adding a model kind does not break older stores.

`GbdtPeakModel::from_bytes` (`crates/scoring/src/gbdt_eval.rs`) decodes the blob. Layout
(little endian): magic `AGBD`, `u32` version (2), `u32 n_features`, `u32 flags`
(bit 0: apply sigmoid to the tree sum), `u32 n_trees`; per tree `u32 n_nodes` followed by the
arrays `feature` (i32, -1 marks a leaf), `threshold` (f32), `left` (i32), `right` (i32),
`value` (f32, leaf output) and `default_left` (u8, 1 sends NaN left); then the calibration
footer `u32 n_iso`, `iso_x` and `iso_y` (f32, empty means identity); then the feature names
`u32 n_names` (must equal `n_features`) and per name a `u32` byte length and UTF-8 bytes. Nothing
may follow the names. `to_bytes` writes the same layout, so the format round-trips. Version 1
blobs (no names) are rejected. Decoding validates each tree: all arrays share one length,
and every internal node's children are in range and have a larger index than the node
(pre-order layout).

## Path A: a tree ensemble into an existing slot

LightGBM, XGBoost and scikit-learn ensembles map onto `Tree` directly. No new runtime or
dependency is involved.

What exists:

- The blob format and its decoder/encoder (above).
- Store writers that take a blob per model: `write_models_with_gbdt` (`store/write.rs`),
  `write_all_models_with_sources_and_gbdt_pub` and `attach_to_model` (`store/update.rs`).
  `attach_to_model` replaces components of one model in a single-file store and keeps the rest
  (this is how `train_intensity.rs` attaches the intensity and rich-ion models).
- The opposite direction: `crates/onnx-infer/scripts/export_agbd_to_onnx.py` converts an
  `AGBD` blob to an ONNX `TreeEnsembleRegressor`.

Not built: a converter from a LightGBM/XGBoost dump to `AGBD`, and a command that attaches an
imported blob to a store row (attaching is a Rust API call today). A converter has to get
these details right, none of which the decoder can detect:

- **Split direction.** `Tree::eval` sends a value left when `v <= threshold`. LightGBM and
  scikit-learn use the same rule. XGBoost sends left when `v < threshold`, so its thresholds
  need the next lower f32 value.
- **Missing values.** `default_left` is read only for NaN. Zero-as-missing modes
  (LightGBM `missing_type=Zero`) and categorical splits have no representation.
- **Node order.** Children must come after their parent; LightGBM's dump (negative leaf ids)
  has to be renumbered into pre-order.
- **Output.** The score is the sum of leaf values. Fold the base score into the leaves,
  set the sigmoid flag for a binary-logistic objective, and provide `iso_x`/`iso_y` only if the
  source model had a calibration step. A classifier slot uses `predict_logit`; the
  regression slots use the raw sum.
- **Feature contract (below).**

### The feature contract

The blob carries the ordered names of its features (`GbdtPeakModel::feature_names`; index `i`
is the feature the tree splits call `i`). `validate_features` (`gbdt_eval.rs`) is called when a
store row is loaded (`read.rs`) and compares them, position by position, with the extractor
for the slot:

| Slot | Names | Count |
|---|---|---|
| `gbdt_model` | `peak_features::FEATURE_NAMES` | 18 |
| `frag_intensity_model` | `frag_features::FRAG_FEATURE_NAMES` | 19 |
| `rich_ion_model` | `ion_features::ION_FEATURE_NAMES` | 23 (the 19 fragment names, then `ion_rank`, `mass_err_ppm`, `complement_present`, `complement_ppm_delta`) |

A different count is reported as such; otherwise the error names the first differing position,
the name the extractor produces and the name the blob has. An ensemble whose features are in
a different order therefore fails at load. Names are required: the trainer writes them from the
dataset it was given (`Dataset::feature_names`), and every bundled blob carries them. A
converter must write the extractor's names in the extractor's order; the contract checks the
labels, not what a feature means, so the fixture score comparison below is still the check that
a ported model computes the same function.

A golden test on a fixture, pinning the ported model's output, belongs with any importer.

## Path B: a new tree-ensemble kind

Use this when the model answers a different question than the three existing slots. It is
Path A plus a nullable column and the plumbing that follows it:

- the column in `store/schema.rs`, and a `None` default in the reader for stores without it
  (follow `rich_ion_model_bytes`: `read.rs` optional column read, `ManifestRow` field, decode
  and `validate_features` against the new extractor's name array);
- a field on `scoring::Param` (`param_model.rs`);
- the writers, which carry one blob per model per column (`store/write.rs`, `store/update.rs`);
- the scoring code that consumes the model, and its feature extractor with a fixed feature
  count.

Bundled stores need rewriting only if they should carry the new model.

## Path C: a neural model (ONNX)

Part of this exists. The ONNX substrate is the `onnx-infer` crate (`crates/onnx-infer`), and
`andes search --predict` uses it (`crates/andes/src/bin/andes/predicted.rs`). It is not a
scoring-model slot: predictors run after the search, on the rows of the written PIN, and add
Percolator features. andes computes no FDR from them.

What exists:

- **Engines.** `tract` (pure Rust, Apache-2.0, default feature of the `andes` crate) and
  ONNX Runtime through `--features ort` (loaded from a shared library at run time; see
  `crates/onnx-infer/README.md`). `--predictor-backend tract|ort` picks one.
- **Manifests.** An `.onnx` file travels with a sidecar JSON, `ModelManifest`
  (`crates/onnx-infer/src/manifest.rs`): model id, input/output tensor specs (fixed or dynamic
  axes), task, `Provenance` (origin, SPDX licence of the weights, notice, notes) and free-form
  `extras`.
- **Wired predictors.** `--predict rt` (the bundled native andes-rt model,
  `resources/predictors/andes_rt.artn`, no engine needed), `ms2` and `peptdeep-rt`
  (AlphaPeptDeep graphs), `deeplc-rt` (DeepLC). The AlphaPeptDeep and DeepLC graphs are not
  bundled: `--predictor-models <dir>` points at exported graphs and manifests
  (`predictors::PEPTDEEP_FILES`, `DEEPLC_FILE`). Exporters for those, MS2PIP and Prosit are in
  `crates/onnx-infer/scripts/`; `crates/onnx-infer/MODELS.md` records versions, licences and
  parity checks.
- **Where inference runs.** On retained rows after the search, batched per predictor
  (`--predictor-batch`), never per candidate (`rescoring::TopNPolicy`).
- **A native alternative.** `onnx_infer::ms2_native` loads a fragment-intensity network
  without an ONNX engine; it is exercised by the examples under `crates/onnx-infer/examples/`
  and is not used by `andes search`.

Not built:

- A predictor that is neither AlphaPeptDeep, DeepLC nor andes-rt needs a `Predictor`
  implementation (`predictors/mod.rs`): the input packing for its graph, the modifications it
  can encode, and its output layout. The export scripts and `ms2pip`/`prosit` input modules
  in `onnx-infer` are the starting point, but only the three predictors above are reachable
  from `--predict`.
- An ONNX model inside the per-candidate scoring loop. That loop evaluates on the order of
  1e5 to 1e6 candidates per run, so a network there removes the speed that the tree-ensemble
  slots keep. Anything that must influence the search itself should be a tree ensemble
  (Path A or B), or be distilled into one.
- Resource discovery for weights shipped beside the binary. The exe-relative probe used for
  `unimod.obo` (`model::unimod::bundled_obo_path`) falls back to a compile-time source path,
  and `unimod::bundled()` then returns an empty database if the file is missing. A new weight
  file needs an error on a missing file instead.

## Traps that apply to every path

- **Declare the regime.** Auto-selection picks a model by the selection tuple. A ported
  model that does not state the activation, instrument, enzyme and protocol it was trained for
  is chosen for data it has never seen. Two cases have occurred: a "low-res" model trained on
  high-resolution Orbitrap data, and an input that skipped instrument detection and was
  searched with a model of the wrong resolution. For predictors, `--predictor-nce` and
  `--predictor-instrument` are part of the regime too.
- **Check the training corpus against every benchmark accession before quoting a number.**
  The training accessions of each bundled model are in the README "Supported models" table
  and in the store's source ledger (`ModelStore::load_sources`). A ported model needs the same
  statement, and overlap with a benchmark dataset makes the benchmark meaningless rather
  than only optimistic. `crates/onnx-infer/MODELS.md` lists what could and could not be
  checked for DeepLC.
- **The serve window is part of the model.** Each model carries its fragment tolerance.
  Serving high-resolution models at their training window cost 21% of Astral identifications
  once. A ported model must bring the tolerance it expects, and that value must be measured
  on a search, not assumed.
- **Licence and provenance.** Record origin, licence of the weights and whether the file is
  verbatim or derived: a NOTICE entry for bundled files (see how `NOTICE` treats the glycan
  databases), and `Provenance` in the manifest for ONNX graphs. Weights under a licence that
  is not Apache-2.0-compatible stay out of the Apache core, like the third-party predictor
  graphs.
- **Merge gate.** A ported model ships as a default only if it identifies more at the same
  true false-discovery proportion, or is faster at equal identifications, than the current
  model. Use an entrapment database where one exists; raw Percolator counts at 1% do not
  compare across models whose score distributions differ. For stores updated with
  `andes train --update`, the held-out acceptance gate of [`TRAIN.md`](../TRAIN.md) section 5
  applies on top.

## The other direction

andes writes only its own format: `GbdtPeakModel::to_bytes`. `export_agbd_to_onnx.py` turns
a blob into an ONNX `TreeEnsembleRegressor` for other tools, including its sigmoid flag and
calibration footer (as manifest `extras`, not applied in the graph). A LightGBM/XGBoost
text export of an andes ensemble is not built.

## Suggested first step

Path A end to end on one model whose answer is known: export an andes ensemble to LightGBM
text, write the converter that turns it back into `AGBD`, and assert that scoring on a fixture
reproduces the original output. That exercises the converter, the split and missing-value
rules and the store attach, with no new dependency, and leaves the converter that Paths A and
B both need.
