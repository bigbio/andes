# andes-ms2: a native fragment-intensity predictor

andes-ms2 v1 is a small CNN (about 0.36M parameters by default) that predicts b/y
fragment intensities. It runs in plain Rust (`onnx-infer::ms2_native`) and is
trained on our own reanalysis spectra. It uses andes-rt's peptide encoding:
per-residue embeddings, with I encoded as L, plus the elemental composition of
each modification.

## Model

**Conditioning.** Each prediction is conditioned on:
- precursor charge;
- NCE (or "unknown");
- activation (`HCD-FTMS`, `HCD-ITMS`, `CID-ITMS`, `CID-FTMS`, `HCD-ASTMS`);
- instrument family.

Conditioning enters through FiLM in every block.

**Outputs.** There are 12 channels per cleavage site:
- b/y at 1+ and 2+;
- the same four with −H3PO4, valid only on fragments that hold pS/pT;
- the same four with −CH4SO, valid only on fragments that hold oxidised M.

The first four channels follow the PeptDeep order, so `rescore_pin --model
andes-ms2` writes the same features as the PeptDeep and MS²PIP arms.

**Loss.** 1 − spectral angle over valid channels. An optional KL term to
AlphaPeptDeep predictions can be added.

## Files

| File | Role |
|---|---|
| `ms2_common.py` | Ion channels, fragment m/z, annotation, filter-string parsing, shard schema |
| `mzml_reader.py` | Streaming mzML reader for selected scans (stdlib + numpy) |
| `ms2_manifest.py` | Datasets → fragmentation group (`HCD-LFQ`, `HCD-TMT`, `CID-TMT`, `CID-ITMS-LFQ`); drops benchmark accessions |
| `ms2_select.py` | Confident PSMs (PEP ≤ 0.01, supported mods, ≤ 4 per precursor) → `selected/*.parquet` + `tasks.tsv` |
| `ms2_extract_raw.py` | One array task: raw → mzML (ThermoRawFileParser) → annotated shard |
| `ms2_extract_flats.py` | The same from andes training flats (peak lists attached; e.g. PXD009875) |
| `ms2_teacher.py` | Adds a `teacher` column with AlphaPeptDeep predictions (optional distillation) |
| `ms2_model.py` | Network, loss, `AMS20001` export |
| `ms2_train.py` | Trainer: sequence-hash split, held-out datasets, stratified test report |
| `ms2_export.py` | Weights plus the Rust parity reference (`--random` writes the test fixture) |
| `ms2_compare.py` | Median spectral angle per predictor on confident PSMs |
| `ms2_run_info.py` | Instrument / activation / NCE of an mzML, for `rescore_pin` |
| `ms2_selftest.py` | End-to-end run on synthetic mzML (about a minute on CPU) |
| `slurm/extract_array.sh`, `slurm/train_gpu.sh` | SLURM wrappers; all paths are arguments |

The Rust side has three parts:
- `src/ms2_native/` holds the forward pass and its unit tests, including parity on `fixtures/ms2_native/tiny.ams2`.
- `examples/ms2_native_check.rs` checks parity and speed for a trained model.
- `examples/common/models.rs` provides `AndesMs2` (`rescore_pin --model andes-ms2`).

## Shard schema

One row per annotated PSM:

| Column | Content |
|---|---|
| `dataset`, `run`, `scan` | Where the spectrum came from |
| `peptidoform` | ProForma with Unimod accessions |
| `sequence`, `charge` | Stripped sequence, precursor charge |
| `nce` | NaN when unknown |
| `activation`, `instrument` | Instrument is the model name from the mzML header |
| `pep` | Posterior error probability of the PSM |
| `n_matched` | Number of plain b/y ions matched |
| `explained` | Matched intensity ÷ total intensity |
| `intensities` | (n − 1) × 12, site-major. −1 = invalid channel; otherwise 0–1, normalised to the most intense valid ion |

Annotation details:
- Tolerance is 20 ppm for FTMS/ASTMS and 0.5 Da for ITMS.
- The precursor m/z must match the peptidoform (isotope −1..3, 20 ppm), which catches run/scan join errors.
- Fewer than 4 matched b/y ions drops the PSM.

## Checks before a cluster run

```bash
python ms2_selftest.py --work WORKDIR                  # synthetic: manifest → select → extract → train → export
cargo test -p onnx-infer --features tract              # includes ms2_native parity on the tiny fixture
cargo run --release -p onnx-infer --example ms2_native_check -- \
    --weights WORKDIR/andes_ms2.ams2 --reference WORKDIR/ms2_reference.tsv
```

## Cluster run plan

Placeholders:
- `$PY`: a Python with numpy/pyarrow/pandas, plus torch for training.
- `$BASE`: the project directory.
- `$QPX`: the qpx `psm.parquet` tables.
- `$RAW`: raw files per dataset.
- `$SIF`: the ThermoRawFileParser image.
- `$FLATS`: the andes training flats.
- `$SCRATCH`: node-local scratch.

Training must never see the benchmark accessions (PXD070049, PXD001819,
PXD007683, PXD007653, PXD005553, PXD002047). `ms2_manifest.py` drops them.

### 1. Manifest and PSM selection (one CPU job, about 1 h, 64 GB)

```bash
S=crates/onnx-infer/scripts/student/ms2
$PY $S/ms2_manifest.py --out $BASE/ms2/datasets.tsv --raw-template "$RAW/{dataset}" $QPX/*/psm.parquet
$PY $S/ms2_select.py --manifest $BASE/ms2/datasets.tsv --out-dir $BASE/ms2/sel \
    --max-per-precursor 4 --psms-per-task 150000          # prints "N tasks (array 1-N)"
```

Check the manifest first:
- every `raw_dir` exists (no `MISSING` in the output);
- the groups match the table in `NOTES.md`.

### 2. Spectrum extraction (SLURM array)

Pilot two tasks and read `part-*.json`. Expect:
- `found` ≈ `selected`;
- `precursor_mismatch` ≈ 0;
- `activation_skipped` ≈ 0 except on EThcD runs.

A high `precursor_mismatch` means the run/scan join is wrong. Fix it before
launching the full array.

```bash
TRFP="singularity exec --bind $BASE $SIF ThermoRawFileParser.sh"
sbatch --array=1-2 --cpus-per-task=2 --mem=8G --time=12:00:00 \
    $S/slurm/extract_array.sh $PY $BASE/ms2/sel/tasks.tsv $BASE/ms2/shards $SCRATCH "$TRFP"
# then the rest; finished tasks are skipped, so resubmitting the same array resumes it
sbatch --array=1-N%40 --cpus-per-task=2 --mem=8G --time=12:00:00 \
    $S/slurm/extract_array.sh $PY $BASE/ms2/sel/tasks.tsv $BASE/ms2/shards $SCRATCH "$TRFP"
```

Add `--trfp-args "-i {raw} -o {outdir} -f 1 -L 2"` to write MS2 only when the
image's ThermoRawFileParser supports `-L` (1.4+).

**Size.** About 4 PSMs per (peptidoform, charge) across 33 datasets is about
10–15M spectra, around 1 KB each before zstd.

### 3. PXD009875 flats (optional quick start, one CPU job)

```bash
$PY $S/ms2_extract_flats.py --dataset PXD009875 --out $BASE/ms2/shards_flats/PXD009875 \
    --activation HCD-FTMS --instrument "Q Exactive" --tol 20 --unit ppm $FLATS/hcd_qexactive_tryp/*.parquet
```

- PXD009875 is also a qpx dataset. Train on the flats **or** on its raw extraction, not both.
- The flats have no NCE, so `--nce` stays unknown.

### 4. Teacher labels (optional, one GPU job)

```bash
$PY $S/ms2_teacher.py --device cuda --out-dir $BASE/ms2/shards_teacher $BASE/ms2/shards
```

Use only HCD rows for distillation. PeptDeep is wrong on ion-trap CID (UPS1 SA
0.376 vs GBDT 0.575, `NOTES.md`). Train CID on real spectra only.

### 5. Training (one GPU)

```bash
sbatch --gres=gpu:1 --cpus-per-task=8 --mem=128G --time=24:00:00 \
    $S/slurm/train_gpu.sh $PY $BASE/ms2/shards $BASE/ms2/runs/v1 --steps 200000 --batch 512 --balance-activation
```

- **Default held-out datasets:** PXD019263 (label-free), PXD005235 (TMT) and PXD012255 (phospho). These are the andes-rt ones.
- **Ablation arms:**
  - `--ch 64 --dilations 1 2 4 8 --n-heads 0` (≈ 0.12M parameters, faster);
  - `--teacher-weight 0.5 --shards $BASE/ms2/shards_teacher`;
  - no `--balance-activation`.
- **Throughput:** the batch is built in Python, so check steps/s in the first log lines and set `--steps` to the time budget.

### 6. Export and parity

`train_gpu.sh` already writes `andes_ms2.ams2` and `ms2_reference.tsv`.

```bash
cargo run --release -p onnx-infer --example ms2_native_check -- \
    --weights RUN/andes_ms2.ams2 --reference RUN/ms2_reference.tsv --threads 8
```

### 7. Benchmark arms (Astral, UPS1, TMT, phospho)

Per dataset, take the PIN and mzML of the existing PeptDeep runs. Read the
conditions with `ms2_run_info.py run.mzML`.

```bash
cargo run --release -p onnx-infer --features tract --example rescore_pin -- --model andes-ms2 \
    --models RUN/andes_ms2.ams2 --pin search.pin --mzml run.mzML --out ams2.tsv \
    --frag-tol 20 --frag-unit ppm --nce 28 --activation HCD-FTMS --instrument "Orbitrap Exploris 480"
# PeptDeep, MS²PIP: the same example with --model peptdeep / ms2pip
# andes GBDT: model-train/examples/frag_model_baseline.rs (writes gbdt_sa)
$PY $S/ms2_compare.py --psms conf.tsv ams2=ams2.tsv pdeep=pdeep.tsv ms2pip=ms2pip.tsv gbdt=gbdt.tsv
$PY crates/onnx-infer/scripts/rescore_pin_arms.py --pin search.pin --baseline-psms base_seed42.psms \
    --features ams2=ams2.tsv pdeep=pdeep.tsv ms2pip=ms2pip.tsv rt=andes_rt.tsv --out-prefix arms \
    --arm ams2=ams2:ms2 pdeep=pdeep:ms2 ms2pip=ms2pip:ms2 ams2_rt=ams2:ms2,rt:rt
```

Percolator 3.7.1, seeds 42 / 1 / 2. Report PSMs **and** peptides at true
FDP ≤ 1% (entrapment) for every arm.

**Pass criteria:**
- spectral angle ≥ PeptDeep on Astral;
- ≥ GBDT on UPS1;
- the `ams2` arm ≥ `pdeep` at equal true FDP.

## Known limits of v1

- Fragment charge stops at 2+, so z ≥ 4 precursors lose their 3+ fragments.
- No a/c/z ions and no glycans: the composition table has none.
- Neutral-loss channels are trained but not yet used by `rescore_pin`, whose `Ion` layout has no loss slot.
