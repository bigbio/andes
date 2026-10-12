# Student MS2/RT predictor — working notes

Goal: a fast CPU fragment-intensity + RT predictor native to andes, distilled
from AlphaPeptDeep and refined on our own spectra. Scope so far is **step 1 only**
(data inventory and in-house baseline). No teacher set up, nothing trained.

## Step 1a — labelled PSM corpus (2026-10-12)

Two sources exist on the compute cluster.

**1. andes training flats.** One slug remains: `hcd_qexactive_tryp`. It holds 40
parquet files from PXD009875 (HeLa, Q Exactive, FTMS HCD) with 526,234 PSMs
(MSFragger expect ≤ 0.01).

Schema: `seq, charge, prec_mz, res_mod_pos, res_mod_delta, nterm/cterm_delta, mz[], intensity[]`.
- Each row carries the full centroided MS2 peak list, so observed b/y intensities
  can be extracted directly.
- There is no RT, NCE or run/scan.

The other slug flats and stores listed in older notes are gone. PXD009875 also
trained the shipped GBDT, so it is not clean held-out data for that model.

**2. quantms reanalysis PSM tables (qpx `psm.parquet`).** There are 43 tables
covering 33 datasets.
- Columns: ProForma peptidoform, charge, PEP, q-value, `is_decoy`, `rt` (s),
  `run_file_name`, `scan`.
- One row is one (run, scan); decoys are about 0.5% of rows, so the tables are
  pre-filtered.
- `mz_array` / `intensity_array` are empty in every file.
- RT is present on 100% of rows.
- The source `.raw` files are still on disk for every dataset in the table below.
  A ThermoRawFileParser pass plus the scan numbers recovers observed spectra and
  NCE.

Target PSMs with PEP ≤ 0.01, latest run per dataset (M = millions):

| Group | Datasets | PSMs |
|---|---|---|
| HCD, label-free, Orbitrap | PXD042233 (26.0M; Exploris/QE mix), PXD005354 (9.3M; Velos HCD), PXD000612 (5.6M; phospho), PXD002395, PXD030406, PXD009875 (GBDT training data), PXD000396, PXD036643, PXD004940, PXD014058, PXD018450 (phospho), PMID23933261-FTMS, PXD000269, PXD013658, PXD001305, PXD005481, PXD010150, PXD019263, PXD012255 (phospho), PXD003209 | ~62.8M |
| HCD, TMT10 | PXD023662 (6.1M), PXD011896, PXD005235 | ~10.3M |
| CID, TMT10 | MSV000085836 (5.1M), PXD026581 | ~6.0M |
| CID, ion trap, label-free | PMID23933261-ITMS | ~5.4M |

- Unique peptidoforms per dataset run from 30k to 290k. All datasets are human,
  so the union is far smaller than the sum, probably ~1–1.5M (not measured).
- None of the benchmark accessions (PXD070049, PXD001819, PXD007683, PXD007653,
  PXD005553, PXD002047) appear.
- Raw-only datasets with no qpx table are not counted (e.g. PXD014525, PXD017199,
  PXD048709 timsTOF `.d`).
- Instrument, label and dissociation come from each run's SDRF config.

## Step 1b — in-house predictors vs PeptDeep on benchmark confident PSMs

Setup:
- PSM set: Percolator q ≤ 0.01 targets from the plain Astral and UPS1 PeptDeep
  rescoring runs (seed 42, base arm).
- PeptDeep features: copied unchanged from that run (`rescore_pin`, Astral
  20 ppm / NCE 30; UPS1 0.5 Da / NCE 35).
- andes features: `model-train/examples/frag_model_baseline.rs` scores the shipped
  frag-intensity GBDT (`hcd_astral_tryp` for Astral, `cid_lowres_tryp` for UPS1,
  300 trees) with the same b/y 1+/2+ layout and the same similarity code.
- Precursor charge comes from the precursor m/z (30 ppm, isotopes 0–3).
- "Flat" = every ion predicted equal.
- Analysis script: `compare_frag_baseline.py`.

| Spectral angle, median (IQR) | Astral (HCD), n = 41,234 | UPS1 (ion-trap CID), n = 15,598 |
|---|---|---|
| PeptDeep generic | **0.660** (0.512–0.788) | 0.376 (0.247–0.515) |
| andes frag GBDT | 0.555 (0.445–0.658) | **0.575** (0.477–0.672) |
| flat | 0.271 | 0.292 |
| Pearson, PeptDeep / GBDT | 0.832 / 0.725 | 0.449 / 0.733 |
| Fraction of PSMs with SA > 0.7, PeptDeep / GBDT | 0.42 / 0.15 | 0.03 / 0.19 |

**Astral.** The GBDT falls further behind as charge and length grow:

| Stratum | PeptDeep SA | GBDT SA |
|---|---|---|
| z = 2 | 0.686 | 0.609 |
| z = 3 | 0.597 | 0.435 |
| z = 4 | 0.607 | 0.333 |
| length 6–10 | 0.717 | 0.632 |
| length 21+ | 0.576 | 0.382 |

The GBDT sees only flanking residues, position fraction, length and charge
(19 features). It has no sequence context.

**UPS1.** The GBDT wins at every charge and length. PeptDeep generic is an HCD
model, and this is ion-trap CID data.

**RT.** Calibration is fitted on half the peptides and scored on the other half
(minutes):

| | Astral MAE / median / p90 | UPS1 (linear) MAE / median / p90 |
|---|---|---|
| PeptDeep iRT, binned-median calibration | **0.42 / 0.28 / 0.86** | n/a (see defect below) |
| PeptDeep iRT, linear calibration | 0.69 / 0.51 / 1.27 | **7.56 / 6.39 / 11.9** |
| andes RT seed index (shipped) | 2.80 / 2.30 / 5.86 | 20.8 / 18.1 / 40.8 |
| andes RT index refit in-run (best case) | 0.88 / 0.66 / 1.89 | 10.9 / 8.7 / 22.2 |

**Defect found.** The binned-median calibration in `peptdeep_pin_arms.py` is not
robust to a single outlier bin.
- On UPS1, one iRT bin has a median of 81 min while its neighbours sit at 37–40 min.
- `np.maximum.accumulate` then carries 81 min forward across about 10 bins. The
  held-out MAE is 21.7 min against 7.6 min for a plain linear fit.
- The UPS1 "+RT" arm of the earlier PeptDeep experiment used this curve.
- Fix to try: isotonic regression (PAV) or LOWESS in place of the cumulative max.

## Not done (scope cut to step 1)

No teacher set-up, teacher predictions, student training, Rust inference or
rescoring arms.

## Step A — RT calibration fix (2026-10-12)

`peptdeep_pin_arms.py` now calibrates RT with isotonic regression (PAV).
- The fit uses confident targets and is applied unchanged to every row.
- Below 20 anchors it falls back to a straight line.
- The code is in `rt_calibration.py`.

Setup:
- Same andes build (fa715f98) and the same PeptDeep features as the VM run.
- The PINs reproduce the VM ones: UPS1 base 7,998 / 2,328 and Astral base 43,076 / 26,723 at seed 42.
- Percolator 3.7.1, seeds 42 / 1 / 2.
- Counts are PSMs / peptides at true FDP ≤ 1%. Entrapment T/E is 2.311 for UPS1 and 0.900 for Astral.

| Arm | old (binned median + cumulative max) | new (isotonic) |
|---|---|---|
| UPS1 rt, s42 / s1 / s2 | 7,781 / 8,283 / 7,790 | 8,337 / 8,082 / 7,830 |
| UPS1 both | 8,541 / 7,805 / 7,865 | **9,017 / 9,238 / 8,485** |
| UPS1 confident-target residual, median / p90 | 8.24 / 34.3 min | **4.69 / 13.8 min** |
| Astral rt | 47,157 / 46,843 / 47,035 | 47,048 / 46,892 / 46,942 |
| Astral both | 48,612 / 48,565 / 48,629 | 48,798 / 48,677 / 48,808 |
| Astral confident-target residual, median / p90 | 0.29 / 0.85 min | 0.29 / 0.86 min |

- **Astral:** unchanged within seed noise. Its binned curve had no outlier bin.
- **UPS1:** the old curve was the problem. "Both" now beats base (7,998–8,281) by about 8%; before it was within noise.
