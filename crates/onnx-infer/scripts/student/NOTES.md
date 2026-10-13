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

## Step B — andes-rt v1, a native RT model trained on our own RTs (2026-10-13)

**Data** (`rt_extract.py`):
- Source: quantms reanalysis PSM tables, 33 datasets, latest run each. Targets with PEP ≤ 0.01.
- Reduced to the exact median RT per (run, peptidoform): 58.8M rows, 1.28M peptidoforms, 6,754 runs.
- No benchmark accessions are present.
- Splits: train / val / test by an md5 hash of the stripped sequence (90 / 5 / 5). Held-out datasets never trained on: PXD019263 (label-free), PXD005235 (TMT) and PXD012255 (phospho).
- Unmapped mods dropped: 13.8k rows, all `UNIMOD:99913` (OpenMS phospho-decoy placeholder on Ala).
- **The quantms sequences carry no isoleucine (I is written as L).** The encoder now maps I→L. Before that fix the I embedding was untrained, and benchmark MAE was 4× worse (Astral 1.72 vs 0.43 min) while the qpx test looked excellent.

**Model** (`rt_model.py`, 77,794 parameters):
- Per position: AA embedding(32) concatenated with the mod composition [C,H,N,O,S,P,other,mass/100] from a built-in Unimod table (`mod_composition.py`, 35 entries, masses checked).
- N-/C-term tokens carry the terminal mods.
- 1x1 conv to 64 channels, then 4 residual dilated k3 convs (1/2/4/8), then attention + mean pooling, then an MLP to a scalar index.

**Training** (`rt_train.py`):
- Loss: pairwise logistic ranking within a run plus a per-run affine MSE head; batches of 16 runs × 64 peptidoforms.
- Hardware: one A100, 30k steps in about 6 min (CUDA 12.8 driver, torch 2.8 cu128). Validation per-run Pearson 0.984; 150k steps gives no gain (0.9839).

**Calibration:** per-run isotonic on confident PSMs (`rt_calibration.py`).

**Native inference:**
- `onnx-infer::rt_native` (plain Rust), `ARTN0001` weight file; weights are not in git.
- Parity with PyTorch: max |Δ| 9.5e-6 over 4,850 peptidoforms in both ProForma and andes PIN notation, including phospho / TMT / pyro-Glu / deamidation / acetyl.
- Speed: 114 µs per peptide on one thread (8 threads: 65k peptides/s). In `rt_pin` on 16 threads, Astral's 801,552 forms take 7.7 s (153 µs per form per thread).
- For comparison on tract: DeepLC 399 s on 8 threads (4.0 ms per form per thread); PeptDeep RT about 1.9 ms per form per thread. andes-rt is about 12× faster than PeptDeep RT and 26× faster than DeepLC.

### Accuracy: qpx held-out sequences
291 runs (≤ 12 per dataset), 187,533 rows. Calibration uses that run's non-test rows; error is on test-split sequences.

| MAE, min (R²) | PeptDeep | andes-rt full | no-phospho | no-TMT |
|---|---|---|---|---|
| all | 5.91 (0.921) | 3.29 (0.980) | 3.21 | 3.69 |
| unmodified/CAM | 3.98 | 3.14 | 3.05 | 3.15 |
| held-out datasets | 11.0 | 3.38 | 3.35 | 4.66 |
| TMT rows | 19.9 (0.41) | 2.89 (0.983) | 2.82 | 6.49 zero-shot |
| phospho rows (n = 607, mostly PXD012255, where every model has R² ≈ 0.3) | 7.12 | 7.04 | 8.33 zero-shot | 6.73 |

Few-shot fine-tuning (3k steps, lr 5e-4, half PTM batches):
- TMT: 100 peptides give 6.47 (no gain), 500 give 4.67. The full model reaches 2.89.
- Phospho: 100 give 7.99 and 500 give 7.67. The phospho slice is too noisy to separate the models.

### Accuracy: benchmark confident PSMs
Calibration on half the peptides, error on the other half. MAE in min; phospho peptides only in brackets.

| | andes-rt | no-phospho | no-TMT | PeptDeep | DeepLC | andes RT index (shipped / in-run refit) |
|---|---|---|---|---|---|---|
| Astral (span 18.7) | **0.43** | 0.42 | 0.41 | 0.43 | 0.44 | 2.75 / 0.89 |
| UPS1 (span 123) | **6.40** | 6.25 | 6.29 | 7.16 | 6.90 | 20.6 / 10.1 |
| TMT a05058 (span 159) | **4.90** | 4.75 | 10.1 zero-shot | 27.6 | 4.75 | 32.4 / 8.73 |
| Phospho PXD007653 (span 249) | **8.67** (8.73) | 9.98 (10.2 zero-shot) | 8.22 (8.32) | 11.2 (11.5) | 9.32 (9.47) | 52.5 / 16.3 |

Zero-shot through composition:
- **Phospho:** the no-phospho model scores 10.2 min on phospho peptides. That beats PeptDeep generic (11.5), but DeepLC is better (9.5) and the full model better still (8.7).
- **TMT:** the no-TMT model scores 10.1 min. That beats PeptDeep (27.6, which has no TMT), but trails DeepLC (4.75) and the full model (4.9).
- **Composition alone gets roughly half-way;** a few hundred labelled peptides are needed to close the gap.

### Rescoring, RT feature only
Features: `rt_diff`, calibrated RT and a missing flag, with isotonic calibration on base confident targets for every predictor. Percolator 3.7.1 `-Y`, seeds 42 / 1 / 2. Counts are PSMs at true FDP ≤ 1%, with peptides in brackets.

| | base | +andes-rt | +PeptDeep RT | +DeepLC |
|---|---|---|---|---|
| Astral (T/E 0.9) | 43,076 / 43,145 / 43,132 | **47,389 / 47,429 / 47,369** (≈ 28,980) | 47,024 / 46,906 / 46,934 (≈ 28,820) | 47,175 / 46,997 / 46,856 (≈ 28,780) |
| Phospho (T/E 1.0) | 37,110 / 37,046 / 37,038 (≈ 14,535) | 37,839 / 37,871 / 37,836 (≈ 14,720) | **37,908 / 37,874 / 37,876** (≈ 14,755) | 37,731 / 37,763 / 37,658 |
| UPS1 (T/E 2.311) | 7,998 / 8,281 / 8,181 | 7,855 / 8,008 / 7,863 | 8,331 / 8,198 / 8,045 | **8,635 / 9,019 / 8,500** |

- **Astral:** andes-rt is about 400 PSMs ahead of PeptDeep, and its nominal 1% sits at 0.95–0.97% true FDP. The no-phospho and no-TMT models give the same counts.
- **Phospho:** andes-rt and PeptDeep are tied. The no-TMT model gives 37,963–38,009.
- **UPS1:** despite the best UPS1 MAE, andes-rt raises the nominal count while true FDP rises to 3.6–3.8%, so there is no gain at true FDP. The cause is not understood (ion-trap CID; T/E-noisy dataset). DeepLC helps there.

### Defects found on the way
- **I→L:** the training sequences carry no isoleucine; fixed by encoding I as L.
- **Peptide parsing:** `K+42.01` at residue 1 has to be read as N-terminal acetyl.
- **Few-shot pool:** the first version filtered out the very PTM peptides it picked.
- **Extraction memory:** Arrow `approximate_median` builds one t-digest per group (125 GB); replaced by an exact sort-based median.
- **Environment:** the CUDA 13 torch wheel cannot use the CUDA 12.8 GPU driver.

Glycopeptides are out of scope for v1. The composition table has no glycans, and none are in the training data.
