# Pre-PR completion plan — everything except glyco (2026-06-14)

User directive: before the closing PR, **validate multiple/non-tryptic enzymes on real
datasets** (Codon-harvested, MSFragger-searched) **and fix every deferred hardening issue**.
**Only glyco stays parked.** This plan is grounded in an actual Codon recon (datasets, models,
raw availability) — not hand-waving.

## 0. Grounding facts (from Codon recon, 2026-06-14)
- **Non-tryptic / non-HCD datasets WITH MSFragger gold PSMs already harvested:**
  | Enzyme/type | Dataset(s) | harvested PSMs |
  |---|---|---|
  | **LysC** | PXD000900 (HCD-highres + CID-lowres) | 75,683 + 51,966 |
  | ETD | PXD004732, PXD010595 | 233k, 214k |
  | UVPD | PXD003109 (human), PXD018176 (E. coli), PXD065289 | 40k/242k/11k |
  | iTRAQ | PXD002214, PXD003702 | 4k/14k |
  | phospho | PXD007740, PXD014525, PXD019646 | up to 453k |
  | LFQ hybrid (HYE) | PXD028735 (ProteoBench) | 354k |
- **GAP: non-tryptic MODELS are NOT trained.** Stores hold only tryptic variants (+tmt/itraq/
  phospho/astral). `cid_lowres_lysc` / `hcd_highres_lysc` flats exist but no store → `--enzyme
  lysc` currently can't load a LysC model (falls back to tryptic). **Training step required.**
- **Raw spectra NOT kept on Codon** (deleted post-harvest) → benchmarking re-downloads the raw.
- andes enzyme is **single** (`params.enzyme: Enzyme`) → **multi-enzyme (2+ proteases in one
  search) is unimplemented.**

---

## Workstream 1 — Multiple + non-tryptic enzymes (the headline)

### 1a. Train the non-tryptic own models (Codon; flats already exist)
Train from the harvested flats, same `train-from-msnet` path as the astral model:
- `cid_lowres_lysc`, `hcd_highres_lysc` (PXD000900 flats).
- (already trained, reuse) `etd_highres_tryp`, `uvpd_qexactive_tryp`.
Ship the new model stores Codon→local→VM. **Verify provenance = own (train-from-msnet).**

### 1b. Non-tryptic head-to-heads (matched, 1% true entrapment-FDP)
For each: re-download raw → convert (docker TRFP) → build a `<enzyme>` entrapment FASTA →
andes (`--enzyme <e>` + matching own model) vs MSFragger (`search_enzyme = <e>`, dt=0) →
percolate → 1% true entrapment-FDP. Datasets: **LysC PXD000900** (primary), then **UVPD
PXD018176** and **ETD PXD004732** if time. This is the real proof `--enzyme` works end-to-end.

### 1c. Multi-enzyme search support (the branch's namesake)
andes currently digests with ONE enzyme. Implement **multi-protease search**:
- Add `params.enzymes: Vec<Enzyme>` (or `--enzyme` accepting a comma list, e.g.
  `--enzyme trypsin,lysc`); `candidate_gen` unions the cleavage sites of all listed enzymes
  (a residue is a cleavage site if ANY enzyme cuts there). Default = single trypsin
  (bit-identical). Model selection uses the primary (first) enzyme.
- Validate: a trypsin+LysC combined search vs the single-enzyme runs (more peptides, FDP-
  controlled). If no native multi-protease dataset is staged, validate on PXD000900 by
  searching `trypsin,lysc` and confirming it ≥ each single-enzyme run at matched FDP.

---

## Workstream 2 — Fix every deferred issue (implement + validate)

| ID | Fix | Validation dataset/approach |
|---|---|---|
| **H5** low-res-HCD routing | `build_selection_key`: `HCD/LowRes → CID/LowRes` (b/y, 0.5 Da) not high-res QExactive; **update equivalence test** to assert the corrected mapping (intentional divergence, documented). | Unit: model_id for HCD/LowRes input. Full PSM validation needs an **HCD-ion-trap dataset** — harvest one (or simulate by forcing the config on a05058) and A/B the model. |
| **H6** EThcD | Detect ETD+supplemental HCD/CID in one activation block (mzML) / Thermo codes 5/6/9/10 → route to **HCD** (b/y) not pure ETD, log it. | Needs an **EThcD dataset** (harvest a PTM/phospho EThcD PXD); unit-test the detection on a synthetic mixed-activation mzML. |
| **L6** match-tier WARN | `select()` returns its match tier (exact/family/empty/last-resort); binary **WARNs** on family/last-resort fallback. | Integration: feed a CID/QExactive input → family fallback → assert WARN. No new data. |
| **M1** Sec / non-standard residues | Add **selenocysteine `U` = 150.95364** to the amino-acid mass table so `candidate_gen` no longer drops U-spans; count+log truly-non-standard drops. | **Data exists**: human FASTA has selenoproteins (GPX1-4, SELENOP). Run PXD009875 (human) → assert U-containing peptides now identified; parity for non-U unchanged. |
| **H1** profile-mode mzML | Detect `MS:1000128` (profile); **centroid** in `build_peaks` (local-max + intensity-weighted) or fail loud. | **Creatable**: convert a `.raw` with ThermoRawFileParser `--noPeakPicking` → profile mzML; assert PSMs recover to ~centroided level. (qe-holdout HeLa.raw on Codon is a source.) |
| **decoy** none/shuffle | Refactor `SearchIndex` to store an explicit `target_count` (drop the `len()/2` invariant); add `--decoy-strategy {reverse,shuffle,none}`; seeded shuffle for reproducibility; pre-existing-decoy detection already warns. | Unit: target-only count, shuffle determinism (fixed seed), reverse bit-identical (default). |
| **cal** threshold tune | The cal-auto pre-pass self-gated on all 3 campaign sets (~220 confident < threshold). Lower `MIN_CONFIDENT_PSMS` / widen the sample so it fires; expose as a flag. | A/B on the **campaign datasets** (cal-off vs cal-auto-tuned) at 1% true entrapment-FDP — keep only if it doesn't raise FDP. |

---

## Workstream 3 — Re-bench + PR
- **Full re-bench** with the final binary: the 3 campaign datasets (TMT/UPS1/Astral, top-1 +
  chimeric) **+** the new non-tryptic (LysC, UVPD, ETD) **+** iTRAQ (C1) — all matched vs
  MSFragger at 1% true entrapment-FDP. Confirm no regression + capture the non-tryptic wins.
- **Parity** after every change; the equivalence test updated where H5 intentionally diverges.
- Then the closing PR (scope per your earlier call).

## Sequencing (each step gated on the prior; parity-tested)
1. **W2 safe/local** first (no datasets): L6 WARN, Sec mass-table, decoy refactor+strategy, profile detect+centroid — implement + unit-test + parity.
2. **W1a** train non-tryptic models (Codon) → ship.
3. **W1b** non-tryptic head-to-heads (LysC first) → real `--enzyme` validation.
4. **W1c** multi-enzyme search support → validate.
5. **W2 dataset-dependent**: H5 (HCD-IT), H6 (EThcD), cal-tune — implement + validate on harvested/created data.
6. **W3** full re-bench → PR.

## Risk / honesty
- H5/H6 still need datasets that aren't staged (HCD-ion-trap, EThcD) — Workstream 1's harvest
  pipeline produces them, but it's real Codon compute + download time. If a dataset can't be
  obtained, the fix ships with unit-test coverage + a clear "validated on synthetic, full
  PSM A/B pending data" note (not silently unvalidated).
- Multi-enzyme `candidate_gen` change is parity-sensitive → default single-trypsin must stay
  bit-identical (gated by the parity suite).
- **Park: glyco mode only** (specs + dormant primitive untouched).
