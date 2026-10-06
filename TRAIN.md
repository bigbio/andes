# Training scoring models with andes

`andes train` builds scoring models from your data into a Parquet model store
(`--out-store`). Search flags are in [`DOCS.md`](DOCS.md).

---

## 1. What a scoring model is, and when to train one

A scoring model holds, per `(charge, parent-mass, fragment-segment)` partition, the
intensity-rank and mass-error statistics of fragment ions. andes ships 17 models in
`resources/models/`, which is also the default seed source. Train your own when your
**instrument** (e.g. Astral, timsTOF) or **experiment class** (TMT/iTRAQ, phospho,
immunopeptidomics, glyco …) is not well covered.

Training is **bootstrap-supervised**: andes searches with a seed model, keeps PSMs at
target-decoy q ≤ a threshold, and learns the statistics from them.

## 2. Quick start

```bash
andes train \
  --spectra mydata.mzML \
  --database mydb.fasta \
  --out-store models.parquet \
  --seed-model hcd_qexactive_tryp \
  --train-fdr 0.01 \
  --instrument OrbitrapAstral \
  --protocol Automatic \
  --model-id astral_tryp
```

This keeps PSMs at ≤ 1% FDR, estimates a model (Laplace smoothing + partition backoff for thin
partitions) and writes it as `astral_tryp` into `models.parquet` (created or appended), with
its **source statistics** for later updates (§5). Search with it:

```bash
andes --spectrum mydata.mzML --database mydb.fasta --output-pin out.pin \
  --model-store models.parquet --model astral_tryp
```

**Key flags** (full list: `andes train --help`):

| Flag | Meaning | Default |
|---|---|---|
| `--spectra` | training spectra (mzML/MGF; `.raw`/`.d` with the native features) | *(required)* |
| `--database` | target FASTA (decoys auto-generated) | *(required)* |
| `--out-store` | Parquet store to create/append | *(required)* |
| `--seed-model` | seed slug or `.param` path for the first-pass search | `hcd_qexactive_tryp` |
| `--train-fdr` | q-value threshold for confident labels | `0.01` |
| `--instrument` | instrument tag for the model | `QExactive` |
| `--protocol` | experiment-class tag(s) (§4) | `Automatic` |
| `--model-id` | id written to the store | `trained_<instrument>_<protocol>` |
| `--mods` | mods.txt (same format as search) | Cam-C + Ox-M |
| `--date` | ISO-8601 date in the source ledger | today |

Use a lenient `--train-fdr` (e.g. `0.1`) on small datasets and `0.01` on full runs.

## 3. Training-data sources

Use data matching your instrument and experiment class, from **PRIDE**
(<https://www.ebi.ac.uk/pride/>; native `.raw`/`.d` work), **ProteoBench**
(<https://proteobench.readthedocs.io/>) or **MassIVE** (<https://massive.ucsd.edu/>). A few
thousand confident PSMs make a usable model; more helps high-charge / high-mass partitions.
Match the enzyme and `--mods` of your searches, and keep a held-out or entrapment set (§7).

## 4. Experiment-class catalog (`--protocol`)

"Protocol" is the sample-prep regime that reshapes fragment statistics. Built-in classes:
`standard`, `phospho`, `tmt`, `itraq`, `acetyl`, `ubiquitin`, `glyco`, `immuno` (aliases fold,
e.g. `Phosphorylation` → `phospho`). Combine them with a comma list (`--protocol tmt,phospho`,
stored as `phospho+tmt`); new classes need no code change. At search time the class comes
from a flag or is inferred from the mods (a TMT mass + a phospho mass → `tmt+phospho`).

## 5. Incremental training (add / remove / reweight / decay)

Updates are exact count arithmetic; no spectra are re-read. Each update makes a **candidate**
that must pass an **acceptance gate** (held-out yield ≥ the current model).

```bash
# Add a new dataset to an existing model
andes train --update astral_tryp --out-store models.parquet \
  --add --spectra more.mzML --database mydb.fasta --source-id batch2 \
  --validate heldout.mzML

# Remove a source
andes train --update astral_tryp --out-store models.parquet \
  --remove-source batch2 --validate heldout.mzML

# Down-weight a source, or decay stale sources by age
andes train --update astral_tryp --out-store models.parquet --reweight batch1=0.5 --validate heldout.mzML
andes train --update astral_tryp --out-store models.parquet --decay 180 --validate heldout.mzML
```

The candidate is committed only if it identifies at least as many target PSMs at 1% FDR on
`--validate` (`--force` commits anyway; without `--validate` the gate is skipped with a
warning). `--decay <days>` weights sources exponentially by age. `add` then `remove` of the
same source restores the model exactly.

## 6. Model selection at search time

Without overrides, search selects by **detected instrument** × **experiment class**, backing
off from an exact match to the largest matching class subset, the instrument family, then a
generic model. `--model-store <path>` uses another store; `--model <id>` forces one model.

## 7. Evaluation & validation

- **Acceptance gate** (§5) for updates.
- **Yield non-regression:** `cargo test -p model-train --test yield_nonregression` with
  `MSGF_TRAIN_BENCH=<dir>` trains a model and asserts its 1% FDR yield ≥ the bundled fallback
  on held-out spectra.
- Judge FDR with an entrapment or held-out set, not raw counts.
