# VM benchmark harness — Phase V strong-score gate

Run on the self-hosted bench VM (`/srv/data/msgf-bench`). These scripts implement
the Phase V validation gate from
`internal-docs/2026-06-06-andes-strong-score-replacement-plan.md`.

## Prerequisites

- Rust toolchain `1.95.0` on VM (`cargo +1.95.0`)
- `thermo` feature build for Astral RAW input
- Docker (Percolator biocontainer)
- Python 3 + duckdb + pyarrow (intensity aggregation)
- 53+ GB free under `/srv/data` for MSNet RAW staging

## Quick start

```bash
# 1. Sync repo from your workstation (or git pull on VM)
bash benchmark/vm/sync_repo_to_vm.sh

# 2. Train intensity model (Phase T) — ~30–90 min depending on downloads
bash /srv/data/msgf-bench/phase_v_train_intensity.sh

# 3. Run the rank vs strong A/B gate (Astral + TMT + optional PTM)
bash /srv/data/msgf-bench/phase_v_strong_score_gate.sh

# 4. Summarize pass/fail
python3 /srv/data/msgf-bench/summarize_phase_v_gate.py /srv/data/msgf-bench/phase-v-*
```

## Gate criteria (all must pass)

| Check | Rule |
|-------|------|
| PSMs@1% | `strong` ≥ `rank` on NORMAL db for each dataset |
| Entrapment FDP | `strong` ≤ `rank` combined FDP on ENTRAP db for each dataset |
| Speed | `strong` search wall ≤ 110% of `rank` search wall (Pass-1) |

Default `--score rank` must remain byte-identical when no `--intensity-model` is set.
The gate runs both arms **with** `--intensity-model` so the numerator is live.

## Scripts

| Script | Purpose |
|--------|---------|
| `sync_repo_to_vm.sh` | Rsync local `msgf-rust` → VM repo path |
| `phase_v_train_intensity.sh` | Download MSNet RAW, aggregate, `andes train-intensity` |
| `phase_v_strong_score_gate.sh` | A/B `--score rank` vs `--score strong` on 3 datasets |
| `summarize_phase_v_gate.py` | Emit PASS/FAIL table from gate output dir |

## CID decision benchmark (before strong-score design changes)

Run MSFragger + Comet + andes (`rank` / `strong`) on **3 CID datasets × 4 raw files**,
after training CID-aware intensity + scoring models from MSNet:

```bash
# Full pipeline (train + benchmark, ~hours depending on downloads)
nohup bash /srv/data/msgf-bench/phase_v_cid_decision_pipeline.sh \
  > /srv/data/msgf-bench/phase_v_cid_decision.log 2>&1 &

# Or stepwise:
bash /srv/data/msgf-bench/phase_v_cid_train_models.sh
bash /srv/data/msgf-bench/phase_v_cid_benchmark.sh
python3 /srv/data/msgf-bench/summarize_phase_v_cid_benchmark.py /srv/data/msgf-bench/cid-bench-*
```

| Dataset | Class | 4 raw files |
|---------|-------|-------------|
| PXD016999 | TMT tissue CID | Instrument1 Fr01/Fr02 + SecondInstrument Fr01/Fr03 |
| PXD007683 | TMT Lumos CID | a05058–a05061 |
| PXD001819 | UPS1 standard CID | UPS1_5000amol R1–R4 |

MSNet training sources (intensity + `train-from-msnet`): `PXD000865-CID`, `PXD016986`,
`PXD009449-CID`. Outputs: `intensity_model_cid.parquet`, `models_cid_msnet.parquet`.

**Strong score is parked** for CID/TMT (see `internal-docs/2026-06-07-cid-rank-lumos-plan.md`).
Focus: expand rank training + close Lumos gap (PXD007683). List candidates:
`python3 benchmark/vm/list_cid_msnet_catalog.py`.

## PTM dataset

No phospho search fixture is pre-staged on the VM. Set these env vars before
running the gate (or export in the script header):

```bash
export PTM_MZML=/srv/data/msgf-bench/ptm-data/sample.mzML
export PTM_FASTA=/srv/data/msgf-bench/ptm-data/database.fasta
export PTM_ENTRAP=/srv/data/msgf-bench/ptm-data/entrapment.fasta
export PTM_MODS=/srv/data/msgf-bench/ptm-data/mods.txt
export PTM_EXTRA="--fragmentation HCD --instrument QExactive --protocol phospho"
```

Suggested MSNet catalog row for intensity training: `PXD000865`, `IPX0002031000`,
`PXD024364-Trypsin` (HCD tryptic, high-res).
