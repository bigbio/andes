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

## Andes production AGBD

In-search peak GBDT stays native `AGBD` (`scoring::gbdt_eval`). Convert for quantms / ort:

```bash
python scripts/export_agbd_to_onnx.py --agbd path/to/model.agbd --out model.onnx
```
