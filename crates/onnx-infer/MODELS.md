# Models

## Committed fixtures (always in git)

Small, license-clean ONNX / AGBD files under [`fixtures/`](fixtures/) used by CI:

| Artifact | Role |
|---|---|
| `double.onnx` | tract smoke (`y = 2x`) |
| `intensity_mlp.onnx` | tract + `RescoringBatch` smoke |
| `tree_leq.agbd` / `.onnx` | AGBD ↔ TreeEnsemble parity |
| `peak_gbdt.agbd` / `.onnx` | multi-tree peak-GBDT shape |

Rebuild: `python scripts/make_local_fixtures.py`  
Parity: `python scripts/test_fixtures.py` (needs `onnxruntime`)

## PeptDeep / AlphaPeptDeep (not in git)

Pretrained `.pth` and exported `.onnx` are **gitignored** under `models/` (tens of MB; separate license/NOTICE). Export when you have ~2 GB free disk:

```bash
curl -L -o /tmp/pretrained_models_v3.zip \
  https://github.com/MannLabs/alphapeptdeep/releases/download/pre-trained-models/pretrained_models_v3.zip
mkdir -p models/peptdeep/pretrained
unzip /tmp/pretrained_models_v3.zip -d models/peptdeep/pretrained

# Python 3.11 recommended; pin numpy<2
python scripts/export_peptdeep_to_onnx.py \
  --pretrained-dir models/peptdeep/pretrained/generic \
  --out-dir models/peptdeep \
  --tasks rt ms2 ccs
```

Then optionally: `cargo test -p onnx-infer -- --ignored`

For `andes --predict ms2 / peptdeep-rt`, the model directory also needs
`peptdeep_mods.json`, written by
[`peptdeep_reference.py`](scripts/peptdeep_reference.py) next to the graphs. By
default it carries every modification andes parses (carbamidomethyl, oxidation,
acetyl, phospho S/T/Y, TMT6plex, deamidation N/Q, pyro-Glu); a peptide with a
modification missing from the file is a missing row for AlphaPeptDeep.

## Andes production AGBD

In-search peak GBDT stays native `AGBD` (`scoring::gbdt_eval`). Convert for quantms / ort:

```bash
python scripts/export_agbd_to_onnx.py --agbd path/to/model.agbd --out model.onnx
```

The graph outputs the raw tree sum (`GbdtPeakModel::predict_value`). The blob's
sigmoid flag and isotonic calibration footer are copied into the manifest
`extras` (`apply_sigmoid`, `iso_x`, `iso_y`) and are not applied in the graph;
`predict_proba` is `isotonic(sigmoid(raw))` (sigmoid only when the flag is set).

## Third-party predictors exported to ONNX (not in git)

Weights stay outside git (`models/`, gitignored); the export scripts below
rebuild them. Every graph runs under tract (the `tract` feature) and ONNX
Runtime (the `ort` feature); parity with the reference Python package is
checked by
`cargo run --release -p onnx-infer --features tract --example model_parity -- --model <m> --reference <tsv> --models DIR [--backend ort]`
(build with `--features "tract ort"` for `--backend ort`).

| Model | Version / source | Licence (code / weights) | Export | Inputs (Rust module) |
|---|---|---|---|---|
| MS²PIP | ms2pip 4.2.0; XGBoost files from the ms2pip model server (Zenodo record 13270668 fallback): `HCD2021` (v20210416), `CID` / `TMT` (v20190107), `CID-TMT` (v20220104) | Apache-2.0 / Apache-2.0 (distributed by the ms2pip project) | [`export_ms2pip_to_onnx.py`](scripts/export_ms2pip_to_onnx.py) | [`ms2pip`](src/ms2pip/) — 139 features per cleavage site, as `ms2rescore_rs.ms2pip_compute_features` |
| DeepLC | deeplc 4.5.0, bundled `multitask_flexcnn_model.pt`, head `PXD005573_mcp` (the uncalibrated `deeplc.predict` output) | Apache-2.0 / Apache-2.0 (weights ship inside the package) | [`export_deeplc_to_onnx.py`](scripts/export_deeplc_to_onnx.py) | [`deeplc`](src/deeplc/) — atom composition, positional and terminal blocks, one-hot |
| Prosit | Prosit 2020 intensity HCD + Prosit 2019 iRT, Koina SavedModels (Zenodo 7681452 / 7671314) | Apache-2.0 (Koina) | [`export_prosit_to_onnx.py`](scripts/export_prosit_to_onnx.py) | [`prosit`](src/prosit/) — residue codes, charge one-hot, NCE/100 |

Reference scripts: [`ms2pip_reference.py`](scripts/ms2pip_reference.py),
[`deeplc_reference.py`](scripts/deeplc_reference.py),
[`peptdeep_reference.py`](scripts/peptdeep_reference.py).

### MS²PIP

- tract implements `ai.onnx.ml.TreeEnsembleClassifier` but not
  `TreeEnsembleRegressor`, and its leaf lookup walks the leaf table from the
  start (O(#leaves) per tree visit), which does not finish on the 9 M-leaf
  `HCD2021` y model. The exporter therefore writes the trees as flat node
  tables and an unrolled `Gather` / `GatherElements` / `Where` traversal (one
  step per tree level). `HCD2021` y: 400 trees, 18.9 M nodes, depth 18, 379 MB
  ONNX, 0.5 s tract load. Exporting it needs the 888 MB `.xgboost` file and
  ~7.5 GB RAM.
- `--format ensemble` writes the same trees as one `ai.onnx.ml`
  `TreeEnsembleRegressor` node (same input, output and manifest names), which
  only ONNX Runtime runs (`--features ort`, `--backend ort`). On ONNX Runtime it
  is about 4x faster than the traversal graph (MS²PIP CID on UPS1: 127 vs
  502 µs CPU per (peptidoform, charge)); parity with ms2pip is unchanged.
- The features ignore modifications (a modified residue encodes as its base
  residue, I = L), as ms2pip does; modifications only enter fragment m/z.
  ms2pip predicts 4..=100 residues; other lengths are rejected.
- Output: `log2(intensity + 0.001)` per cleavage site for b and y (z = 1).
- Training data ([ms2pip docs](https://github.com/compomics/ms2pip/blob/main/docs/source/prediction-models.rst)):
  CID = NIST CID Human library (eval NIST CID Yeast); HCD2021 = combined set
  incl. MassIVE-KB; TMT = Peng Lab TMT library; CID-TMT = PXD041002.

### DeepLC

- Trained across 6,543 LC setups; 5,553 setup names carry a PRIDE/MassIVE
  accession. None names PXD001819, PXD070049 or PXD007683 (the other ~990
  names are run labels and cannot be checked).
- The graph is the model's own forward with three ops rewritten into
  equivalent tract-compatible ones (one-hot matmul for the embedding, one-hot
  sum for `scatter_add` counts, a −1e30 offset instead of a `-inf` fill for
  the masked max pool).
- Modification compositions come from unimod via psm_utils. If psm_utils
  cannot open its unimod database (for example a Python whose `sqlite3` fails
  to import), DeepLC drops every modification with only a warning; the
  reference script must run in an environment where
  `Peptidoform("M[Oxidation]").parsed_sequence[0][1][0].composition` works.

### Prosit — blocked

Koina's Prosit SavedModels (MD5 as listed in Koina's `.zenodo` pointers) are
built on the GPU-only `CudnnRNN` GRU kernel: TensorFlow on a CPU cannot
execute them, so no local TensorFlow reference exists. tf2onnx 1.17 converts
them to ONNX `GRU` without error, but the result disagrees with Koina's served
predictions (`prosit_reference.py`): iRT max |Δ| ≈ 48 iRT units (Pearson 0.91
over the 8-peptide panel), normalised intensities max |Δ| = 1.0. Setting
`linear_before_reset = 1` does not fix it, and the GRU weights are sliced from
the cuDNN parameter blob inside the graph, so gate-block reordering of the
initialisers has no effect. tract also cannot run the converted graphs: the
bidirectional encoder's backward GRU is fed through `ReverseSequence`, which
tract 0.21 does not implement (`Unimplemented(ReverseSequence)` at optimise).
The Rust side (`prosit` module, `Prosit` predictor in the examples,
`model_parity --model prosit` against Koina) is in place; what is missing is a faithful graph — likely a re-implementation of the
Prosit architecture (e.g. in PyTorch) from the SavedModel's canonical
weights, or a CPU-compatible Keras export.
