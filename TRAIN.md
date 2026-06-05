# Training scoring models with andes

andes can **generate its own scoring models** from your data and store them in a single
Parquet model store (`resources/ionstat/models.parquet` by default). This guide covers training a
model from scratch, where to get training data, the experiment-class catalog, incremental
updates, and how a model is selected at search time.

For the full CLI/parameter reference see [`DOCS.md`](DOCS.md).

---

## 1. What a scoring model is, and when to train one

A scoring model captures, per `(charge, parent-mass, fragment-segment)` partition, the
intensity-rank and mass-error statistics of fragment ions — the numbers the generating-function
scorer turns into per-peak scores. andes ships 39 models (consolidated into one Parquet
store) covering common fragmentation × instrument × enzyme × protocol combinations.

**Train your own when:**
- your **instrument** isn't well covered (e.g. Orbitrap Astral, Bruker timsTOF), or
- your **experiment class** differs systematically from the bundled models (TMT/iTRAQ labeling,
  phospho-enrichment, immunopeptidomics, glyco …), or
- you simply want a model tuned to your own acquisition.

Training is **bootstrap-supervised**: andes searches your data with a seed model, keeps the
confident PSMs (target-decoy q-value ≤ a threshold), and learns the statistics from them.

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

This searches `mydata.mzML` with the seed model, keeps PSMs at ≤ 1% FDR as labels, accumulates
the per-partition histograms, estimates a model (Laplace smoothing + partition backoff so thin
partitions still get sensible scores), and writes it as `astral_tryp` into `models.parquet`
(created if absent, appended otherwise). The model carries its training **source statistics** so
it can be updated later (§5).

Search with the new model:

```bash
andes --spectrum mydata.mzML --database mydb.fasta --output-pin out.pin \
  --model-store models.parquet --model astral_tryp
```

(Without `--model-store`/`--model`, search auto-selects from the bundled store by detected
instrument — see §6.)

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

> Use a **lenient `--train-fdr`** (e.g. `0.1`) on small datasets so enough labels are collected;
> use `0.01` on full runs.

## 3. Training-data sources

Train on a dataset that matches the instrument/experiment class you're targeting:

- **PRIDE** (<https://www.ebi.ac.uk/pride/>) — public proteomics repository; pick projects on your
  instrument/protocol. Native `.raw`/`.d` work directly.
- **ProteoBench** (<https://proteobench.readthedocs.io/>) — curated LFQ/DDA reference datasets
  (e.g. the Orbitrap Astral DDA set used in this repo's benchmarks).
- **MassIVE** (<https://massive.ucsd.edu/>) — another large public source.

Guidance:
- A few thousand confident PSMs are enough for a usable model thanks to partition backoff; more
  is better, especially for high-charge / high-mass partitions.
- Match the **enzyme** and **modifications** of your downstream searches (pass `--mods`).
- For trustworthy FDR, keep an **entrapment**/held-out portion for validation (§7).

## 4. Experiment-class catalog (`--protocol`)

"Protocol" is the sample-prep regime that reshapes fragment statistics — not arbitrary PTMs.
Classes are **combinable** (a model is tagged with a set, canonicalized like `phospho+tmt`). The
built-in catalog: `standard`, `phospho`, `tmt`, `itraq`, `acetyl`, `ubiquitin`, `glyco`,
`immuno` (with alias folding, e.g. `Phosphorylation` → `phospho`). New classes are added by
training + tagging — no code change. Combine with a comma list, e.g. `--protocol tmt,phospho`.

At search time the experiment class is taken from a flag or **inferred from the configured mods**
(e.g. a TMT mass + a phospho mass → `tmt+phospho`).

## 5. Incremental training (add / remove / reweight / decay)

Because training accumulates **counts**, a model is updated with exact count arithmetic — no
spectra are re-read. Every update produces a **candidate** that must pass an **acceptance gate**
(yield on held-out data ≥ the current model) before it's committed.

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

- The candidate is searched against `--validate` and compared to the current model at 1% FDR; it
  is committed **only if it identifies at least as many** target PSMs (`--force` commits anyway).
- If `--validate` is omitted, the gate is skipped (a warning is printed) and the candidate is
  committed.
- `--decay <days>` applies exponential age weighting using each source's recorded date.

`add` then `remove` of the same source restores the model **exactly** (the per-source statistics
are retained losslessly).

## 6. Model selection at search time

With no `--model`/`--model-store` override, search auto-selects from the bundled store by the
**detected instrument** (mzML cvParams, Thermo `.raw` vendor metadata, Bruker `.d` → timsTOF) ×
the **experiment class** (flag or mod-inferred), with a backoff ladder: exact match → largest
matching experiment-class subset → instrument-family fallback → generic. Train an
instrument-specific model into the store and it is selected automatically for that instrument.

Overrides: `--model-store <path>` (use a different store), `--model <id>` (force a specific
model), `--param-file <path>` (load a binary `.param` directly).

## 7. Evaluation & validation

- **Acceptance gate** (§5) — the built-in held-out yield check for updates.
- **Yield non-regression** — the repo ships an env-gated harness
  (`cargo test -p model-train --test yield_nonregression`, set `MSGF_TRAIN_BENCH=<dir>`) that
  trains a model and asserts its 1% FDR yield ≥ the bundled fallback on held-out spectra.
- Judge FDR with an **entrapment**/held-out set, not just raw counts, consistent with how this
  project validates search quality.
