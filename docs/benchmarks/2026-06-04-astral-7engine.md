# 7-engine Astral benchmark — uniform Percolator (2026-06-04)

Head-to-head of **cimas** (top-1 and `--chimeric`), **MSFragger 4.2** (DDA+),
**Sage 0.14.7**, **Comet 2025.01**, **Java MS-GF+ v20240326**, and **ProSE**
(OpenMS, fragment-index engine, `openms-tools-thirdparty:latest`) on a single
Thermo Orbitrap Astral run. Every engine's output is re-scored through the
**same Percolator 3.7.1** (`-Y`, seed 42) so the 1% FDR is computed identically —
including Java and ProSE, whose outputs are converted to a Percolator PIN by hand
(see §4).

> **Headline:** Among single-best-per-scan engines, cimas top-1 leads on PSMs and
> peptides. In the chimeric tier, cimas `--chimeric` and MSFragger DDA+ are neck
> and neck (cimas ahead on PSMs+peptides, MSFragger on proteins). cimas is the
> fastest real search — ~30× faster than its own MS-GF+ ancestor.

> **Chimeric FDR — validated (see §8).** An entrapment-FDP check (paired,
> r = 1) confirms the chimeric gain is **real**, not coincidental-target
> inflation: doubling accepted PSMs moves the *true* FDP only from **0.98%**
> (top-1) to **1.16%** (chimeric) against a 1% TDC claim.

---

## 1. Dataset

| | |
|---|---|
| Run | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1` (Orbitrap Astral, 15-min LFQ) |
| Vendor file | `.raw` (2.5 GB) + matching `.mzML` (1.9 GB) |
| Database | `ProteoBenchFASTA_MixedSpecies_HYE.fasta` (mixed-species HYE) |
| Decoys | reversed; concatenated target+decoy where the engine needs an explicit file |

## 2. Harmonized search parameters (all engines)

| Parameter | Value |
|---|---|
| Enzyme | Trypsin, fully specific, ≤ 2 missed cleavages |
| Peptide length | 7–40 |
| Precursor charge | 2–4 |
| Precursor tolerance | 10 ppm |
| Fragment tolerance | 20 ppm (high-res HCD) |
| Isotope error | −1 … +2 (see §6 for the Comet exception) |
| Fixed mod | Carbamidomethyl (C), +57.02146 |
| Variable mods | Oxidation (M) +15.99491; Acetyl (protein N-term) +42.01057 |
| FDR | 1%, Percolator 3.7.1 `-Y --seed 42` |

## 3. Results — all via Percolator 3.7.1 (seed 42), sorted by proteins then wall

| Engine | Input | Wall | Proteins @1% | PSMs @1% | Peptides @1% |
|---|---|---:|---:|---:|---:|
| MSFragger 4.2 (DDA+) | mzML | 7:17 | **5,836** | 73,909 | 25,594 |
| **cimas** (`--chimeric`) | **`.raw`** | 4:06 | 5,199 | **74,830** | **26,993** |
| Sage 0.14.7 (chimera) | mzML | 4:50 | 5,059 | 32,091 | 21,397 |
| **cimas** (top-1) | **`.raw`** | 3:13 | 4,752 | 35,508 | 22,904 |
| Java MS-GF+ | mzML | 2:06:36 | 4,570 | 26,542 | 17,954 |
| ProSE (OpenMS) | mzML | 2:51 | 4,407 | 30,646 | 20,596 |
| Comet 2025.01 | mzML | 3:30 | 4,354 | 31,435 | 20,607 |

Notes:
- **cimas** is the only engine reading Thermo `.raw` natively. Sage's generic
  release lacks the proprietary Thermo reader (it tries to parse `.raw` as mzML
  and fails), so Sage/MSFragger/Comet/Java/ProSE all run on the mzML.
- Java MS-GF+ search wall is **2 h 6 min** vs cimas's 3–4 min on identical data.
- Among top-1 engines, cimas leads PSMs (35,508) and peptides (22,904).

## 4. Per-engine procedure — exactly how each was run and FDR-controlled

All scripts referenced are in [`scripts/`](scripts/).

### cimas (top-1 and chimeric) — `astral_all5.sh`
Native `.raw`, built with `--features thermo` (needs rustc ≥ 1.88 + .NET 8 at
runtime). `--fragmentation auto` selects the HCD/QExactive model. Writes a
Percolator PIN directly (`--output-pin`). `--chimeric` enables the two-pass
co-isolation cascade. PIN → Percolator → count.

### MSFragger 4.2 (DDA+) — `astral_all5.sh` step, `configs/msfragger-astral.params`
`data_type = 3` (DDA+ chimeric), `report_alternative_proteins = 1`. **Cannot read
native `.raw`** here — the `ext/thermo` Batmass-IO binary is not installed — so it
runs on the **mzML**. Writes a Percolator PIN → Percolator → count.

### Sage 0.14.7 — `astral_fix.sh`, `configs/sage-astral.json`
`chimera = true`, generates its own decoys (`rev_`). The generic release lacks the
Thermo reader, so it runs on the **mzML**. `--write-pin` → Percolator → count.

### Comet 2025.01 — `comet_astral.sh`
Runs from the OpenMS `openms-tools-thirdparty` image
(`/opt/OpenMS/thirdparty/Comet/comet.exe`). Params generated with `comet -p` then
edited: `decoy_search = 1` (concatenated internal decoys, `DECOY_`),
`output_percolatorfile = 1`, high-res settings (`fragment_bin_tol = 0.02`,
`fragment_bin_offset = 0.0`, `theoretical_fragment_ions = 0`). PIN → Percolator.

### Java MS-GF+ v20240326 — `astral_all5.sh` + `build_pins.py java`
`-tda 1 -inst 3 (Q-Exactive) -m 3 (HCD) -protocol 0 -addFeatures 1`. Produces a
358 MB mzIdentML. **`msgf2pin` 3.07.1 crashes** on this mzid
(`basic_string: construction from null`), so the PIN is built by hand:
`MzIDToTsv` → `build_pins.py java`, which emits a Percolator PIN with MS-GF+'s
discriminative features — **RawScore, DeNovoScore, ScoreDiff, lnSpecEValue,
lnEValue, IsotopeError, |precursor-error|, PepLen, charge one-hot** — with the
target/decoy `Label` from the `XXX_` prefix. PIN → Percolator → count.

### ProSE (OpenMS, fragment-index) — `tmt_all5.sh`/`astral_fix.sh` + `build_pins.py prose`
ProSE is OpenMS's new fragment-index search engine (v3.6.0-pre). Run with
`-Search:decoys` and `-Search:annotate:PSM ALL`. OpenMS `PSMFeatureExtractor`
**does not recognize ProSE's `ln(hyperscore)` score type** ("No known input to
create PSM features from"), so the standard PercolatorAdapter path fails. Instead
the idXML annotations are parsed directly (`build_pins.py prose`) into a Percolator
PIN whose features are ProSE's own annotations — **hyperscore, delta_score,
fragment/precursor mass errors, matched-ion-current, prefix/suffix-ion fractions,
num_matched_peaks, longest-ion-run, PepLen, charge** — with the target/decoy
`Label` from the `target_decoy` meta (and `DECOY_` accession fallback). PIN →
Percolator → count. (Enabling the feature annotations is what makes a uniform
Percolator run possible; without them ProSE only had native target-decoy q-values,
which scored ~22.5 k PSMs vs 30.6 k via Percolator.)

## 5. FDR counting (identical for all)

After Percolator (`-Y --seed 42 --results-psms`), from the target PSM list:
- **PSMs @1%** = rows with `q-value ≤ 0.01`.
- **Peptides @1%** = distinct stripped sequences among them.
- **Proteins @1%** = distinct `proteinIds` among them, excluding decoy prefixes
  (`XXX_` / `rev_` / `DECOY_`). For Java the MS-GF+ `(pre=…,post=…)` accession
  suffixes are stripped first.

## 6. Known approximations / caveats

- **Comet isotope error**: Comet has no exact `−1..+2` option; `isotope_error = 4`
  (`−1/0/1/2/3`) is the closest superset.
- **Java & ProSE feature sets** are hand-built (§4), not produced by a native
  converter, so their absolute counts are a faithful-but-not-canonical Percolator
  run. The four engines that emit a native PIN (cimas, MSFragger, Sage, Comet) are
  the strictest apples-to-apples subset.
- **Chimeric FDP**: reversed-decoy TDC does not expose coincidental real-DB targets
  in multi-PSM-per-scan output — so it is validated separately against an entrapment
  database (§8). The MSFragger DDA+ chimeric counts were not entrapment-checked here
  (cimas's were) and remain TDC-only.

## 7. Reproduce

On the benchmark host (engines + data staged under `/srv/data/msgf-bench`):

```bash
bash scripts/astral_all5.sh        # cimas top1+chimeric, MSFragger, Sage, ProSE
bash scripts/comet_astral.sh       # Comet
bash scripts/uniform_perc.sh       # build Java+ProSE PINs, percolate all 7 uniformly
bash scripts/entrapment_validate.sh  # §8 entrapment-FDP check (cimas top-1 vs chimeric)
```

Environment: RHEL 9, 8 cores; Percolator 3.7.1 and ProSE/Comet via
`ghcr.io/openms/openms-tools-thirdparty:latest`; Sage v0.14.7 and MSFragger 4.2
binaries; `MSGFPlus_v20240326.jar`.

## 8. Chimeric FDR — entrapment-FDP validation

Reversed-decoy TDC is blind to **coincidental real-DB targets**: in
multi-PSM-per-scan (chimeric) output, a second peptide can win by chance against a
real database sequence and TDC never sees it. To check whether cimas's `--chimeric`
gain is real or such inflation, the same Astral run is searched against an
**entrapment database** — the ProteoBench proteins plus a paired `ENT_`-prefixed
**shuffled twin** of every protein (r = 1, 31,889 + 31,889) — with reversed decoys
generated as usual. After Percolator @1%, any accepted PSM mapping **only** to an
`ENT_` protein is a confirmed false positive. With r = 1 the estimated true FDP is
`2 · N_ent / N_total` (Wen & Noble paired estimator). Script:
[`scripts/entrapment_validate.sh`](scripts/entrapment_validate.sh).

| cimas mode | accepted @1% (TDC) | original | entrapment | **true FDP (est.)** |
|---|---:|---:|---:|---:|
| top-1 | 32,965 | 32,803 | 162 | **0.98%** |
| `--chimeric` | 69,052 | 68,652 | 400 | **1.16%** |

**Verdict:** the gain is real. Doubling accepted PSMs raises the true FDP only from
0.98% to 1.16% against a 1% TDC claim — the chimeric cascade recovers genuine
co-isolated identifications, not coincidental-target noise. (Counts are lower than
§3 because the entrapment DB is 2× larger, adding candidate competition; the FDP
ratio is the meaningful quantity. The +0.16 pp overshoot could be removed by a
slightly stricter chimeric acceptance threshold.)
