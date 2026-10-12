# onnx-infer

Tool-agnostic **ONNX inference substrate** for proteomics: database search
rescoring, de novo sequencing, MS2/RT predictors.

**Future quantms role:** Rust home for
[quantms-rescoring](https://github.com/bigbio/quantms-rescoring) inference
(AlphaPeptDeep, DeepLC, spectrum features → PIN / idparquet), sharing feature
names and PeptDeep I/O with OpenMS. Developed in the andes workspace today
([#92](https://github.com/bigbio/andes/issues/92)); extractable later.

## Features

| Feature | Backend | Notes |
|---|---|---|
| *(none)* | — | Manifest / tensor / rescoring / peptdeep API (always on) |
| `tract` | [tract](https://github.com/sonos/tract) | Pure Rust, Apache-2.0 — preferred for Apache-clean builds |
| `ort` | reserved | Broader ops (`ai.onnx.ml` TreeEnsemble, many Transformers) |

Default features are empty so Apache-core builds do not pull an engine until opted in.

## Modules

| Module | Purpose |
|---|---|
| `OnnxSession` / `ModelManifest` | Load + run ONNX with sidecar provenance |
| [`rescoring`](src/rescoring/) | quantms-rescoring feature names, `CandidateKey`, `TopNPolicy`, `RescoringBatch` |
| [`peptdeep`](src/peptdeep/) | PeptDeep ONNX tensor packing (OpenMS export contract) |

## Speed + data-model contract

1. **Top-N only** — use [`rescoring::TopNPolicy`](src/rescoring/psm.rs); never call ONNX inside the full candidate loop.
2. **Feature names** — short forms matching `quantmsrescore.constants` / AlphaPeptDeep generator (`spec_pearson_norm`, `rt_diff`, …) so andes PIN, idparquet, and quantms-rescoring stay aligned.
3. **Search-engine agnostic keys** — [`CandidateKey`](src/rescoring/psm.rs) (`spectrum_ref`, `rank`, `charge`, `peptide`); map from andes `PsmMatch`, Comet, or Sage at the boundary.
4. **Native AGBD stays hot** — andes in-search peak GBDT keeps the `AGBD` blob evaluator; ONNX TreeEnsemble is the interchange / quantms path (`scripts/export_agbd_to_onnx.py`).

## Model export / porting

| Script | From → To |
|---|---|
| [`scripts/export_agbd_to_onnx.py`](scripts/export_agbd_to_onnx.py) | andes `AGBD` v1 → ONNX `TreeEnsembleRegressor` |
| [`scripts/export_peptdeep_to_onnx.py`](scripts/export_peptdeep_to_onnx.py) | PeptDeep `.pth` → ONNX (+ manifests), OpenMS-compatible names |

Fixtures (committed):

| File | Kind | Runner |
|---|---|---|
| `double.onnx` | `y=2x` | tract CI |
| `intensity_mlp.onnx` | MatMul+ReLU rescoring smoke | tract CI |
| `tree_leq.agbd` / `.onnx` | 1-feature TreeEnsemble | ort (`ai.onnx.ml`) |
| `peak_gbdt.agbd` / `.onnx` | 3-feature multi-tree (andes peak shape) | ort |

```bash
# Rebuild local fixtures (no torch)
python scripts/make_local_fixtures.py

# Round-trip AGBD → ONNX
python scripts/export_agbd_to_onnx.py \
  --agbd fixtures/peak_gbdt.agbd --out /tmp/peak_gbdt_rt.onnx

# PeptDeep (needs ~2GB free disk + peptdeep/torch; gitignored under models/)
curl -L -o /tmp/pretrained_models_v3.zip \
  https://github.com/MannLabs/alphapeptdeep/releases/download/pre-trained-models/pretrained_models_v3.zip
unzip /tmp/pretrained_models_v3.zip -d models/peptdeep/pretrained
python scripts/export_peptdeep_to_onnx.py \
  --pretrained-dir models/peptdeep/pretrained/generic \
  --out-dir models/peptdeep --tasks rt ms2 ccs
```

## Quick start

```toml
onnx-infer = { path = "../onnx-infer", features = ["tract"] }
```

```rust
use onnx_infer::rescoring::{TopNPolicy, ALPHAPEPTDEEP_CORE_FEATURES};
use onnx_infer::peptdeep::{pack_rt_inputs, NCE_DEFAULT};
use onnx_infer::{BackendKind, OnnxSession};

let policy = TopNPolicy::new(5)?;
let _ = ALPHAPEPTDEEP_CORE_FEATURES;
let batch = pack_rt_inputs(&["PEPTIDE", "SAMPLER"], 30)?;
// let mut session = OnnxSession::from_manifest("models/peptdeep_rt_dynamic.json", BackendKind::Ort)?;
// session.run(&batch.as_inputs())?;
let _ = (policy, batch, NCE_DEFAULT, BackendKind::Tract);
```

## Tests

```bash
# Always-on contracts (no tract): feature-name freeze, Top-N, manifests, PeptDeep packing
cargo test -p onnx-infer --test contracts

# Tract + committed fixture models (also run in CI)
cargo test -p onnx-infer --features tract

# AGBD ↔ ONNX Runtime parity + MLP reference (local; needs onnxruntime)
python scripts/test_fixtures.py
```

See [`MODELS.md`](MODELS.md) for committed fixtures vs PeptDeep weights.

## Measured: AlphaPeptDeep rescoring

`andes` → AlphaPeptDeep RT and MS2 predictions (`examples/rescore_pin.rs`, tract backend) →
extra PIN features (`scripts/peptdeep_pin_arms.py`) → Percolator 3.7.1 `-Y`, seeds 42/1/2.
AlphaPeptDeep `pretrained_models_v3` generic. Counts are at entrapment true FDP ≤ 1%: an
entrapment version of each database, FDP = entrapment hits / accepted × (1 + T/E).

| Astral (high-res HCD) | PSMs | Unique peptides |
|---|---|---|
| andes | 43,076–43,145 | 26,723–26,746 |
| andes + AlphaPeptDeep (MS2 + RT) | **48,565–48,629 (+12.7%)** | **29,638–29,665 (+10.9%)** |
| andes `--chimeric` | 60,403–60,435 | 26,696–26,717 |
| andes `--chimeric` + AlphaPeptDeep (MS2 + RT) | **63,323–63,372 (+4.9%)** | 27,849–27,876 (+4.3%) |

- On chimeric searches, AlphaPeptDeep adds 4.9% PSMs, against 12.7% on the plain search.
- Chimeric + AlphaPeptDeep gives the most PSMs, but **plain andes + AlphaPeptDeep gives the
  most unique peptides** (29,648 vs 27,861). If peptides or proteins are the target metric,
  plain andes + AlphaPeptDeep is the best setup.
- With the added features, Percolator's nominal 1% sits above 1% true FDP: 1.18–1.22% on the
  plain search and 1.28–1.35% on the chimeric one. The counts above are compared at equal true
  FDP.

On UPS1 (ion-trap CID) the gain is within seed noise. The PeptDeep MS2 model is mostly trained
on HCD, and its RT fits poorly there: the median residual is 8.2 min for confident targets,
against 0.28 min on Astral.

The prediction step is expensive: tract runs MS2 at about 16 ms per peptide per thread, against
onnxruntime's 3.1 ms. The plain Astral PIN (1.02M rows) took 31.6 min and 14,180 CPU-s on 8
threads. The chimeric PIN (288k rows) took 11.7 min.

## Roadmap

1. Offline #92 A/B with PeptDeep ONNX → five additive PIN columns.
2. Wire `ort` for TreeEnsemble + PeptDeep graphs tract cannot run.
3. Port quantms-rescoring feature generators onto this crate incrementally.
4. Optional de novo Transformer export (decode may stay outside ONNX).
