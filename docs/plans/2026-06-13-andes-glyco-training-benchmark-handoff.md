# Glyco Neutral-Loss Model: Real-Corpus Training + Benchmark Handoff

> **Status:** the *code* for neutral-loss-aware glyco scoring is complete, reviewed,
> and byte-identical for standard searches (Tasks 1–7 + SP3 + SP4-mechanism, on
> branch `feat/enzyme-support`). This document is the remaining, **dataset-gated**
> work: training a real glyco model and benchmarking the PSM gain. It needs a
> glyco MS/MS dataset (mzML + FASTA) that the local checkout does not have.

## What is already done (no dataset required)

The full neutral-loss machinery is implemented and unit-proven end-to-end on a
synthetic corpus (`crates/model-train/tests/neutral_loss_glyco.rs`):

- **Config** (Tasks 1–1b): mods.txt `loss=<m1;m2>` / `class=<glyco|phospho|sulfo|generic>` /
  `accession=<CURIE>` attributes; `Modification.{neutral_losses, loss_class}`.
- **Prediction** (Task 3): activation-gated, per-class loss-shifted b/y ions.
- **Scoring** (Task 7): `score_psm` adds a peptide-aware loss contribution —
  for each fragment spanning a loss-bearing residue, the model's pooled per-class
  loss rank table is probed at `intact_mz − loss/z`. Byte-identical when the
  peptide has no loss mod OR the model has no loss table (every bundled model).
- **Training** (SP3): `ion_match_facts` derives loss-ion facts from the intact ion
  vocabulary × active losses; `build_rank_dist_table` materialises the per-class
  loss rank tables; serialization round-trips them (`ion_loss_class` column).
- **Output**: TSV `Modifications` column emits `pos:UNIMOD:393` (PIN unchanged —
  its `Proteins` column is rest-of-line).

**Byte-identical guarantees verified:** full workspace suite (only the documented
pre-existing `fragment_tolerance_override_changes_model` failure), the
`precursor_cal_off` golden BSA search (RawScore unchanged; only the additive TSV
column differs), and `score_psm` Java parity.

## SP3 — Train a real glyco model

**Inputs you must provide:**
1. A glyco MS/MS dataset: `glyco.mzML` (HCD or CID — the model's activation must
   predict losses; ETD does not). Stepped-HCD glyco runs are ideal.
2. A target FASTA (`glyco.fasta`) for the same sample (decoys auto-generated).
3. A `glyco_mods.txt` declaring the glyco mod with its losses + accession. Start
   from `resources/mods/glyco_example.txt`:
   ```
   NumMods=2
   57.021464,C,fix,any,Carbamidomethyl
   340.100562,K,opt,any,Glucosylgalactosyl,loss=162.0528;324.1056,class=glyco,accession=UNIMOD:393
   ```
   Adjust the mod mass / residue / `loss=` list to the real glycan(s) in the
   sample. The `loss=` masses are the **labile-bond losses** (e.g. −Hex 162.0528,
   −Hex2 324.1056), NOT the composition delta.

**Command** (the HCD seed makes the trained model HCD ⇒ `predicts_neutral_losses`):
```bash
cargo build --release -p andes
./target/release/andes train \
  --spectra   glyco.mzML \
  --database  glyco.fasta \
  --mods      glyco_mods.txt \
  --seed-model hcd_qexactive_tryp \
  --protocol  glyco \
  --model-id  glyco_hcd \
  --train-fdr 0.01 \
  --out-store models_glyco.parquet
```

**Acceptance checks on the produced store:**
- The trained `glyco_hcd` model has at least one `loss_class != 0` entry in its
  `rank_dist_table` (i.e. `RankScorer::has_loss_tables()` is true). If not, no
  confident PSM produced a matched loss ion — check activation, the `loss=`
  masses, and that confident glyco PSMs exist at `--train-fdr`.
- The 39 bundled models are untouched (this writes a *new* store / model id).

There is no `--activation` flag on `train`; activation is inherited from the
seed model, so **always seed from an HCD/CID model for glyco** (the default
`hcd_qexactive_tryp` is HCD).

## SP4 — Benchmark / validation

**A. Glyco PSM gain (the win):** search the glyco dataset twice, FDP/entrapment-
controlled, and compare PSM/peptide counts at 1% FDR:
1. Baseline: search with a no-loss model (or the glyco model with `loss=` removed
   from the mods file) — loss ions are never predicted/scored.
2. Loss-aware: search with `glyco_hcd` + the loss-declaring `glyco_mods.txt`.
   Expect a glyco-PSM gain from the recovered loss-ion signal.

Use the same entrapment-FDP protocol as the other andes benchmarks (target-only
FASTA + entrapment; see `reference_andes_infra_layout` — benchmark on the VM, not
Codon). Report the gain and the FDP, not just raw counts.

**B. Standard-3 byte-identical regression (the guarantee):** run the standard
Astral / UPS1 / a05058-TMT searches with a **bundled** (non-glyco) model and no
loss mods; the PIN + RawScore must be byte-identical to the pre-feature baseline
(loss ions are never predicted or scored). The synthetic test already proves the
code path is inert; this confirms it on real data.

Do **not** claim a benchmark pass without running A on the real dataset — the
synthetic test proves the *mechanism*, not the *magnitude* of the gain.

## Gotchas

- **Activation gating is real:** an ETD glyco run will (correctly) predict/score
  no loss ions. Use HCD/CID data.
- **`loss=` ≠ mass delta:** declare the labile-bond loss masses, not the Unimod
  composition. pY phospho is the trap — declare HPO₃ (−79.9663), not H₃PO₄.
- **Per-class pooling:** the loss table is pooled per `loss_class`, split per
  `(charge, offset)`. Multiple loss masses in one class (Hex, Hex2) pool into the
  same per-class table — by design (matches how training counts them).
- **Loss keys are not in `frag_off_table`:** they live only in `rank_dist_table`
  (intact `frag_off` layout is unchanged). `RankScorer` builds the loss pass from
  `rank_dist_table`, and the writer serializes all `rank_dist_table` entries, so
  this is handled — just don't expect loss ions in `frag_off_table`.

## Carry-forward minor cleanup
- `IonType::loss_class()` returns 0 for `Noise` too; tighten the doc + add a
  `Noise.is_loss() == false` assertion when next editing `param_model.rs`.
