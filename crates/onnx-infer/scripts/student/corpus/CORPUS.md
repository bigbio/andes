# Training corpus for andes-rt and andes-ms2

This directory holds the design, the catalog and the quality control for the corpus
that trains the native predictors: andes-rt (retention time) and andes-ms2 (fragment
intensities). It is a design and survey step. Nothing here trains a model.

Status on 2026-10-13:
- The scripts work.
- One worked QC example has been run (PXD009875 flats).
- The PRIDE metadata survey is complete.
- The cluster steps (MsNet inventory, qpx profiling, raw headers) wait for the
  cluster to come back.

| File | Role |
|---|---|
| `corpus_catalog.py` | Catalog schema (`FIELDS`, `spec` → `SCHEMA.tsv`). Fillers per source: `msnet`, `qpx`, `proteometools`, `pride-ids`, `pride` (metadata only), `other`. Also `merge`. |
| `corpus_qc.py` | QC metrics. Input kinds: `flat` (MsNet flats with peaks), `psm_table` (qpx), `msms` (MaxQuant msms.txt), `shard` (andes-ms2 shards). Subcommands `peptidoforms` and `overlap` give cross-source duplicates. |
| `msnet_inventory.py` | Walks an MsNet collection. Writes per-file and per-dataset inventories plus `qc_commands.sh` (one QC command per dataset). |
| `pride_survey.py` | PRIDE survey. Subcommands: `fetch` (paged, with a count check), `classify`, `summary`, `rank`, `files` (paged, with a `total_records` check), `peek-zip` (reads a remote zip's member list by HTTP Range and can extract one member). |
| `corpus_sampling.py` | Stratum weights, per-dataset caps, dataset and PTM-class hold-outs. |
| `slurm/run_line.sh` | Array wrapper: task *i* runs line *i* of a commands file. |
| `SCHEMA.tsv` | The catalog fields: name, type, level and description. |
| `results/qc_PXD009875/` | Worked QC example: `dataset_qc.json` and `run_qc.tsv`. |
| `results/catalog/` | Catalog rows (dataset + 40 runs) for PXD009875. |
| `results/survey/` | `strata_summary.tsv`, `candidates_files_checked.tsv` (1,004 projects), `usability_by_stratum.tsv`, `shortlist.tsv` (56), `massive_hla_search.tsv` (147). |

All paths are command-line arguments. Benchmark accessions are hard-excluded in
`corpus_catalog.EXCLUDED`, `ms2_manifest.BENCHMARK` and `pride_survey`:
PXD070049, PXD001819, PXD007683, PXD007653, PXD005553, PXD002047.

---

## 1. Sources

| Source | What it gives | Spectra | Status |
|---|---|---|---|
| **MsNet collection** (cluster long-term storage) | MsNet-format annotated PSMs | Peak lists attached | Not readable until the cluster's data-mover/compute nodes return. Inventory and QC are runbook step 1. |
| **MsNet-format flats, PXD009875** (benchmark VM) | 526,234 PSMs, 40 runs, MSFragger expect ≤ 0.01 | Attached | QC done (§3) |
| **qpx reanalysis** (`psm.parquet`, 33 human datasets) | ~85M target PSMs at PEP ≤ 0.01: RT, run, scan, PEP, q | None (raw on disk, TRFP + scan join) | Profiling is runbook step 2 |
| **PRIDE submissions with their own IDs** | mzIdentML / mzTab / PRIDE XML, or engine output (MaxQuant msms.txt, PD .msf, pepXML, Mascot .dat) | mzML/MGF/apl in about half; else raw | Survey done (§4). Readers: `corpus_qc.py msms` now; mzIdentML/mzTab later |
| **ProteomeTools** (synthetic) | MaxQuant txt per pool (zip) | raw; msms.txt carries matched-peak m/z and intensities | Accessions and formats verified (§4.3) |
| **MassIVE HLA** (mono-allelic sets) | Varies | Varies | Accessions verified through PROXI; formats not checked |

---

## 2. Catalog schema

There are 58 fields (`SCHEMA.tsv`; `python corpus_catalog.py spec`). There are two
tables, each written as parquet with a TSV twin.

**`datasets`**, one row per (source, accession):
- Provenance: `source`, `accession`, `title`, `year`, `licence`, `doi`, `pubmed`.
- Biology: `organism`, `taxid`, `sample_type`, `enrichment`.
- Chemistry: `labelling`, `enzyme`, `enzyme_specificity`.
- Acquisition: `instrument_model`, `analyser`, `fragmentation`, `nce`,
  `ms2_resolution`, `gradient_min`, `acquisition`.
- Identification provenance: `id_engine`, `id_fdr_method`, `id_format`,
  `spectra_available`.
- Totals: `n_runs`, `n_psms`, `n_unique_peptidoforms`, `n_unique_peptidoforms_il`.
- Content: `ptm_inventory` (JSON of Unimod → PSMs; unmapped deltas kept as
  `mass:+x`), `n_ptm_types`, `length_hist`, `charge_hist`, `tryptic_cterm_frac`,
  `isoleucine_written_as_L`.
- QC: `qc_ann_frac_median`, `qc_ann_frac_p10`, `qc_low_ann_frac`,
  `qc_n_peaks_median`, `qc_prec_err_ppm_median`, `qc_prec_fit_frac`,
  `qc_frag_err_da_median`, `qc_charge_alt_fit_frac`, `qc_decoy_auc`, `qc_pep_p99`,
  `qc_outlier_runs`.
- Governance: `excluded_benchmark`, `in_model_store_ledger`, `predictor_heldout`,
  `contamination_note`, `priority`, `notes`.

**`runs`**, one row per (source, accession, run):
- The acquisition fields as read from the run itself: `activation`
  (e.g. `HCD-FTMS`, `CID-ITMS`), `analyser` and `nce`.
- Counts and the run-level QC, plus `qc_outlier_flags`.

**Provenance of each value.** `analyser_src` and `fragmentation_src` say where a
value came from: `cv`, `text`, `qc` or `manual`. Missing values stay null, never
guessed.

**Precedence.** A run-level value read from the scan header or inferred by QC beats
the PRIDE CV term, which beats free text.

---

## 3. QC metrics and the PXD009875 worked example

### Metric definitions (`corpus_qc.py`)

**Annotated-intensity fraction.** The fraction of total MS2 intensity explained by
plain b/y 1+/2+ ions, with each peak counted once.
- It is computed at 20 ppm and at 0.5 Da.
- The headline `ann_frac` uses the tolerance of the analyser *inferred* for the run.

**Analyser inference.** Uses the median |error| of b/y 1+ ions matched in a 0.5 Da
window.
- < 0.01 Da with the 20 ppm annotation intact → **FTMS**.
- > 0.05 Da, or 20 ppm loses half the annotation → **ITMS**.
- Anything else → `mixed?`.

This catches mislabelled data even when the metadata is wrong.

**Peaks per spectrum and TIC.**

**Precursor error.** Best isotope offset in −1..3 at the stated charge.
- `prec_fit` is the share within 20 ppm. A low value means a run/scan join error.
- `charge_alt_fit` is the share that fits better at z ± 1.

**Fragment error.** Median signed error in Da and ppm per PSM.

**Deamidation/isotope artefact.** The share of deamidated PSMs whose precursor needs
isotope −1, compared with the same rate on other PSMs. Deamidation (+0.984 Da) plus
an isotope −1 pick (−1.003 Da) is just the unmodified peptide.

**Composition.**
- Length and charge histograms.
- PTM inventory: Unimod via `mod_composition`; MaxQuant names via `MQ_MODS`.
- Tryptic C-terminus share, missed cleavages.

**Isoleucine handling.** `isoleucine_written_as_L` is true when no sequence has I but
L is present. Cross-source keys always merge I and L.

**Scores (`psm_table`, `msms`).**
- PEP and q quantiles of retained targets.
- Decoy share.
- Target/decoy AUC of the best score column (Mann–Whitney), where decoys exist.

**Duplicates across sources.** `peptidoforms` writes the I/L-merged keys per
dataset; `overlap` reports shared and unique counts and the union.

**Run-level outliers.**
- Robust z (median, 1.4826·MAD) of PSM count, median `ann_frac`, peaks, precursor
  error and fragment error.
- Computed *within runs of the same inferred analyser*, with |z| > 3.5 flagged.
- Also flagged: `analyser_mixed`, and `precursor_join` when `prec_fit` < 0.95.

### PXD009875 flats

Command:

```bash
python corpus_qc.py flat --dataset PXD009875 --out-dir QC/PXD009875 --workers 6 FLATS/*.parquet
```

It took about 1 minute on 6 cores. Results are in `results/qc_PXD009875/`.

**Basic shape**

| | |
|---|---|
| PSMs / runs | 526,234 / 40 |
| Unique peptidoforms | 66,455 (45,518 sequences); I is kept, so the I/L-merged count is the same |
| Precursor fits the peptidoform (20 ppm, isotope −1..3) | **100.0%**. No join errors, no charge inconsistency (0.0%), CAM present in the flats (0 rows needed it added). |
| Precursor error quantiles (ppm) | 1%: −11.2, 50%: +0.2, 99%: +15.2 |
| Precursor isotope chosen | 0: 77.4%; +1: 13.9%; −1: 6.6%; +2: 2.1% |
| Tryptic C-terminus | 98.8% |
| Missed cleavages | 0: 95.9%; 1: 3.9%; 2: 0.2% |
| Charge | z2 72.0%, z3 26.0%, z4 1.9%, z5–6 0.07% |
| Length | median 14 (1%–99%: 7–34) |
| Modifications (PSMs) | Carbamidomethyl 110,517; **Deamidated 61,204 (11.6%)**; Oxidation 21,128. Nothing else. |

**By inferred analyser**

| | FTMS (16 runs) | ITMS (24 runs) |
|---|---|---|
| PSMs | 330,650 (62.8%) | **195,584 (37.2%)** |
| `ann_frac` median (p10–p90) | 0.284 (0.147–0.462) at 20 ppm | 0.183 (0.098–0.338) at 0.5 Da |
| PSMs with `ann_frac` < 0.2 | 23.9% | 56.9% |
| Peaks/spectrum median (p10–p90) | 151 (60–331) | 796 (409–1,246) |
| Fragment error, median (5–95%) | −1.3 ppm (−3.4 to +0.6) | −10 ppm (−121 to +81); median abs 0.066 Da |
| `ann_frac` by charge (z2/z3/z4) | 0.297 / 0.258 / 0.234 | 0.186 / 0.179 / 0.151 |
| `ann_frac` by length (≤10 / 11–15 / 16–20 / 21–30 / >30) | 0.23 / 0.29 / 0.33 / 0.33 / 0.34 | 0.25 / 0.20 / 0.18 / 0.16 / 0.15 |

### Findings that change how these flats can be used

**1. 37% of the "`hcd_qexactive_tryp`" flats are ion-trap MS2.** The affected runs are
the 21 MudPIT runs and the 3 `_FL_` runs.

- `ms2/README.md` step 3 extracts *all* flats with `--activation HCD-FTMS --tol 20 --unit ppm`.
- That annotates ion-trap spectra at 20 ppm, where they keep only ~2% of their 0.5 Da
  annotation, so they are dropped as too-few-ions or, worse, mislabelled.
- **Fix:** split the flats by inferred analyser using `run_qc.tsv`
  `analyser_inferred`, and extract the ITMS runs with `--tol 0.5 --unit da`.
- The ITMS activation is not in the flats. The PRIDE record lists Velos/Elite/Lumos.
  Treat it as `CID-ITMS?` (unverified) until the raw scan headers are read
  (runbook step 3). On a Lumos, ion-trap MS2 can also be HCD.

This also re-confirms the 2026-06 note that PXD009875 is a mixed high-res/low-res
dataset.

**2. Most "Deamidated" labels are isotope mis-picks.**
- Of the 61,204 deamidated PSMs, **45.5% need isotope −1**, against 1.5% of other
  PSMs. Overall, 63.0% use a non-zero isotope, against 17.3% of the rest.
- By analyser, the isotope −1 rate is 31% (FTMS) and 55% (ITMS).
- These PSMs are the unmodified peptide carrying a wrong +0.984 label: a wrong
  peptidoform for both RT and MS2 training. The paper searched only C+57 and M-ox;
  Deamidated came from the SDRF-driven reharvest.
- **Fix:** drop deamidated PSMs whose isotope ≠ 0, or re-search without Deamidated.
  Apply the same gate to every MSFragger-harvested source. This is the same defect
  class as the PR #101 deamidation fix.

**3. The run outliers are size outliers, not quality outliers.**
- Within each analyser group, only the three `_FL_` runs are flagged, on PSM count
  (35–40k) and peak count.
- The MudPIT ITMS runs show a fragment-error drift across salt steps: median
  −0.087 Da (step 1) to +0.026 Da (step 10). This is harmless at 0.5 Da but is a
  calibration signal worth keeping in the run table.

**4. Not computable on flats.** The flats carry targets only, with no scores or
q-values, so decoy separation and the q/PEP distribution are N/A. Those metrics are
exercised on `msms` input instead. ProteomeTools citrullination (one pool, 4,854
rows):
- target/decoy Score AUC 0.81; decoy share 0.02%;
- MaxQuant intensity coverage median 0.56 for HCD-FTMS and 0.43 for CID-ITMS.

---

## 4. PRIDE survey (metadata only, no data downloaded)

### 4.1 Method

**`fetch`.** Pages `/v3/projects` at 100 per page. Larger page sizes are silently
clamped to 100, so the script always pages through.
- Two full independent passes returned the same **41,465** projects, while
  `/projects/count` says **42,027**.
- 40 randomly sampled unlisted PXD numbers all return 404 from `/projects/{acc}`, so
  the 562-project gap is most likely non-public records. This was not verified for
  all 562 (`--allow-short` documents the choice).
- Submitter and PI contact fields are dropped at fetch time.

**`classify`.** Assigns strata from:
- PRIDE CV terms: instruments → family/analyser; organisms; identified-PTM CV;
  quantification method.
- Text patterns over the title, description and both protocols: fragmentation, MS2
  analyser, labelling, enrichment, HLA, enzymes, DIA, synthetic, plasma.

Precision was spot-checked by reading 10–20 sampled titles per stratum:
- HLA ≈ 85% after tightening. "Neoantigen" or "antigen presentation" alone no longer
  counts.
- Ubiquitin ≈ 90%. The earlier `k-?e-?gg` pattern matched every KEGG analysis.
- Acetyl ≈ 60%.
- ETD ≈ 60% after removing "Orbitrap Fusion **ETD**" instrument names.

Every text-derived stratum is a *heuristic*. Verify each shortlisted dataset against
its paper or run headers.

**`files`.** For the top 900 by `rank` plus named datasets (1,004 projects), pages
`/v3/projects/<acc>/files`. All 1,004 listings matched the `total_records` header.
Formats are then classified:
- `standard_ids`: mzIdentML, mzTab or PRIDE XML.
- `engine_ids`: MaxQuant txt, PD .msf/.pdResult, pepXML, Mascot .dat/xml, FragPipe,
  Sage, or a generic `search*.zip` / per-pool zip that is MaxQuant (verify with
  `peek-zip`).
- `research_needed` otherwise.
- `+peaks` when mzML/MGF/apl are present, `+raw` when only raw files are.

### 4.2 Counts per stratum (projects)

All values are from `results/survey/strata_summary.tsv`. COMPLETE = PRIDE COMPLETE
submission (standard IDs inside).

| Stratum | All | COMPLETE | COMPLETE & not DIA |
|---|---|---|---|
| Total | 41,465 | 6,949 | 6,713 |
| Human | 19,874 | 2,978 | 2,865 |
| Q Exactive family / Exploris | 15,799 / 2,947 | 2,464 / 402 | 2,421 / 378 |
| Tribrid (Fusion/Lumos/Eclipse/Ascend) | 9,862 | 1,658 | 1,630 |
| LTQ-Orbitrap hybrid (Velos/Elite/XL) | 7,037 | 1,624 | 1,616 |
| **Astral** | 435 | **12** | **6** |
| **timsTOF** | 2,467 | 174 | 165 |
| Other Q-TOF (SCIEX/Waters/Agilent) | 3,720 | 642 | 520 |
| Ion trap only / FT-ICR | 569 / 210 | 147 / 41 | 147 / 41 |
| MS2 in ion trap (text) | 2,719 | 570 | 565 |
| HCD / CID (text) | 11,423 / 4,699 | 1,782 / 941 | 1,750 / 925 |
| **ETD / EThcD / ECD+UVPD** (text) | 529 / 338 / 43 | 87 / 42 / 4 | 86 / 42 / 4 |
| TMT (any) / **TMTpro** | 5,888 / 1,418 | 863 / 256 | 846 / 253 |
| **iTRAQ** / SILAC | 1,327 / 2,151 | 180 / 159 | 176 / 154 |
| Phospho / glyco | 3,435 / 1,056 | 555 / 169 | 544 / 165 |
| **Ubiquitin-GG / acetyl / methyl** | 374 / 478 / 163 | 37 / 69 / 19 | 37 / 65 / 18 |
| **Acyl (succ/cro/mal/lac…) / citrullination / ADP-ribosyl** | 486 / 88 / 132 | 46 / 14 / 21 | 43 / 12 / 21 |
| Histone / crosslink | 1,347 / 2,237 | 194 / 339 | 187 / 338 |
| **HLA / immunopeptidome** (HLA-I / HLA-II) | **509** (402 / 139) | 137 (112 / 36) | 131 (110 / 33) |
| **GluC / AspN / chymotrypsin / ArgC / LysN** | 414 / 200 / 562 / 208 / 28 | 73 / 39 / 129 / 43 / 0 | 69 / 38 / 122 / 43 / 0 |
| Non-specific (text) | 1,847 | 327 | 318 |
| Synthetic peptides | 497 | 92 | 91 |
| Plasma / serum | 5,897 | 1,058 | 1,001 |
| Licence CC0 (rest: "EBI terms of use") | 34,347 | 5,370 | 5,172 |

COMPLETE submissions run at 700–930 per year since 2021. Astral data exists almost
only as PARTIAL.

**Usable without re-searching**, among the 1,004 files-checked candidates
(`usability_by_stratum.tsv`):

| Stratum | Checked | Standard IDs | Engine outputs | Re-search needed |
|---|---|---|---|---|
| HLA (any) | 410 | 136 | 143 | 131 |
| HLA-II | 123 | 38 | 56 | 29 |
| ETD/EThcD | 115 | 54 | 26 | 35 |
| timsTOF | 157 | 39 | 60 | 58 |
| Astral | 10 | 2 | 2 | 6 |
| Other Q-TOF | 108 | 42 | 28 | 38 |
| Non-tryptic enzyme (non-HLA) | 193 | 89 | 65 | 39 |
| Rare-PTM enrichment | 176 | 57 | 75 | 44 |
| TMTpro / iTRAQ | 43 / 25 | 18 / 10 | 13 / 7 | 12 / 8 |
| Synthetic | 79 | 22 | 32 | 25 |
| **All checked** | 1,004 | 341 | 367 | 296 |

About 70% of the gap-relevant candidates can be used without a search.

**Formats we need to parse, in order of yield:**
1. MaxQuant msms.txt (done: `corpus_qc.py msms`).
2. mzIdentML 1.1/1.2 (PEAKS, Comet/MHCquant, Mascot, Byonic exports).
3. PD .msf (SQLite).
4. mzTab.
5. pepXML.

### 4.3 Verified details

Zip contents were read without downloading, from the central directory by HTTP Range
(`peek-zip`):

**ProteomeTools.** Every pool is a zip of MaxQuant txt (msms.txt, evidence.txt,
msmsScans.txt). Peeked members:
- PXD004732 `TUM_first_pool_61_..._3xHCD-1h-R2-tryptic.zip`
- PXD010595 `02208a_GC3-..._ETD-1h-R1.zip`
- PXD021013 `TUM_HLA_41_..._DDA-1h-R1-unspecific.zip`
- PXD023119 `TMT_TUM_isoform_75_...zip`
- PXD056559 (citrullination)
- PXD009449 `SEARCH_Kmod_Methyl.zip`

msms.txt carries per-scan `Fragmentation` and `Mass analyzer` (HCD-FTMS and CID-ITMS
in the same pool), `PEP`, `Score`, `Reverse` and the matched peaks (`Matches`,
`Intensities`, `Masses`). Full spectra and NCE still need the raw files.

The ProteomeTools family, by title from the PRIDE records (all CC0 except PXD004732
and PXD009449, which are "EBI terms of use"):

| Accession | Content |
|---|---|
| PXD004732 | Part I |
| PXD010595 | Part II, Prosit |
| PXD021013 | Part III: HLA I/II and non-tryptic |
| PXD023119 | Part IV: TMT tryptic |
| PXD023120 | Part V: TMT non-tryptic |
| PXD056559 | Part VI: citrullination |
| PXD009449 | 21 PTMs (Zolg 2018) |
| PXD010871 | Prosit figure data |

PXD004732 and PXD010595 are in the model-store ledger.

**HLA MaxQuant zips.** The following were checked:
- PXD000394 (Bassani-Sternberg 2015). This is a 12 GB zip; extract only msms.txt with
  `peek-zip --extract`.
- PXD004894 and PXD006939, both **in the model-store ledger**: train only.
- PXD012308.
- PXD010154 (29 tissues, multi-protease, 2 TB raw). msms.txt was checked inside the
  per-tissue `*_txt.zip`.

**MassIVE** (PROXI API; formats not checked):
- MSV000084172: Sarkizova 2020 mono-allelic HLA-I, Lumos + QE Plus.
- MSV000080527: Abelin 2017 mono-allelic.
- PXD014293: HLA-II ligand processing rules, a MassIVE-hosted PXD.

The 147 HLA/MHC MassIVE hits are in `massive_hla_search.tsv`.

### 4.4 Shortlist

`results/survey/shortlist.tsv` has 56 datasets. Each row carries the stratum, a
value note, the ID format and how it was verified, licence, raw GB, DOI/PubMed, URL
and contamination flag. Ranked by value to the gaps in §5:

| # | Accession | Stratum | Why | IDs |
|---|---|---|---|---|
| 1 | PXD021013 | Synthetic HLA-I/II, non-tryptic | ProteomeTools III; HCD NCE ladder, IT-CID and ETD per pool | MaxQuant zip ✔ |
| 2 | PXD019643 | HLA-I + II | HLA Ligand Atlas: 29 benign tissues; Comet + Percolator | mzIdentML + mzML ✔ (COMPLETE) |
| 3 | MSV000084172 | HLA-I | 95 mono-allelic lines | MassIVE (unchecked) |
| 4 | PXD023120 | Synthetic TMT non-tryptic/HLA | ProteomeTools V | zip (pattern) |
| 5 | PXD009449 | Synthetic PTMs | 21 PTM classes (acetyl, methyl, succinyl, crotonyl, GG, nitro, citrulline …) | MaxQuant zip ✔ |
| 6 | PXD056559 | Synthetic citrullination | ProteomeTools VI | MaxQuant zip ✔ (QC'd) |
| 7 | PXD000394 / PXD012308 | HLA-I / HLA-II | Bassani-Sternberg 2015 / Racle 2019 | MaxQuant zip ✔ |
| 8 | PXD022930, PXD012437, PXD012498 | HLA + **EThcD** | Lumos/Fusion EThcD vs HCD of the same ligandomes | mzid+mzML / PD .msf |
| 9 | PXD055268, PXD054498, PXD038165, PXD029648 | HLA on **timsTOF** | Mouse/human HLA-I; HLA-II on timsTOF | mzid + peaks (COMPLETE) |
| 10 | PXD014754, PXD025877, PXD000450, PXD002951 | HLA on **TripleTOF** | Q-TOF fragmentation; PXD000450 has phosphorylated ligands | mzid + peaks / PD |
| 11 | PXD057177 | HLA + **drug-hapten PTM**, ZenoTOF 7600 | Benzylpenicillin-modified ligands: an unexpected PTM | mzid + peaks (COMPLETE) |
| 12 | PXD004732, PXD010595, PXD023119 | Synthetic tryptic / TMT | NCE ladder; the Prosit sets | MaxQuant zip ✔ (ledger: train only) |
| 13 | PXD058160, PXD061848 | **Astral** TMTpro | The only COMPLETE Astral DDA sets with IDs + peaks; mouse | mzid / mzTab + apl |
| 14 | PXD010154 | **Multi-protease** | 29 tissues | MaxQuant zip ✔ |
| 15 | PXD022403, PXD007058, PXD028432, PXD057316 | **ETD/EThcD** + PTMs | Arg-methyl/phospho, multi-protease (PXD007058 is in the ledger) | PD / MaxQuant / mzid |
| 16 | PXD057588, PXD047750 | HLA + **TMTpro** | Labelled immunopeptidomes | mzid + peaks |
| 17 | PXD052972, PXD043291, PXD057347, PXD018758, PXD021214, PXD022077, PXD037092, PXD034992 | **Rare PTMs** | Histone GG, yeast ubiquitylome (ion trap), acyl-K validation, succinylome (timsTOF), N-terminal acetylome (ArgC) | mzid + peaks mostly |
| 18 | PXD001574 | **iTRAQ** on TripleTOF/Synapt/QE, ETD | Phospho-iTRAQ | mzid + peaks |
| 19 | PXD013112, PXD066497, PXD000900 | Non-tryptic | LysC/GluC/AspN; endogenous neuropeptides (timsTOF); Confetti (ledger, re-search needed) | |
| 20 | PXD022949, PXD022884, PXD009531, PXD029029, PXD015240 | More HLA | Melanoma, ex-vivo tissues, HLA-C/G, citrullinated MHC-II, bovine BoLA | |

The full list with licences and URLs is the TSV.

**Licences.** 42 are CC0 and 11 are "EBI terms of use". The MassIVE ones are
unchecked. Treat "EBI terms of use" as usable for research and training, but record
it in the catalog for the Apache-release audit
(`project_independence_license_status`).

**Contamination.** None of the 56 is a benchmark. Seven are in the bundled
model-store ledger: PXD004732, PXD010595, PXD004894, PXD006939, PXD000900, PXD007058
and PXD014502. Use them to train; never use them to test a bundled model.

**PXD000138.** The synthetic phospho library is also a candidate phospho
ground-truth benchmark (`reference_benchmark_training_contamination`). It can play
one role only.

---

## 5. Coverage-gap analysis

**What we hold today:**
- qpx: 33 human datasets, ~85M PSMs:
  - HCD-LFQ 62.8M (4 phospho);
  - HCD-TMT10 10.3M;
  - CID-TMT 6.0M;
  - ion-trap CID LFQ 5.4M (one dataset);
  - all tryptic, all Orbitrap families, no spectra yet.
- MsNet flats: PXD009875 (HeLa), 331k FTMS-HCD + 196k ITMS.
- The cluster MsNet collection: unknown until inventoried.
- The `mod_composition` table: 35 Unimod entries.
- andes-ms2 v1 channels: b/y 1+/2+ and the H3PO4 and CH4SO losses. No c/z, no glycans.

Value = (model need) × (current scarcity) × (availability without search). ★ marks
the user's priorities.

| Rank | Gap (dimension) | Current corpus | What the models need | Available (projects, usable/checked) | Value |
|---|---|---|---|---|---|
| 1 | ★ **HLA-I peptides** (8–12-mers, C-term not K/R, z1–2, no tryptic y-ion bias) | **0** | RT and MS2 on non-tryptic short peptides; the tryptic C-terminal K/R prior is wrong for them | 402 projects (110 COMPLETE DDA); 279/410 HLA candidates usable; ProteomeTools III synthetic HLA; MassIVE mono-allelic | **Highest** |
| 2 | ★ **Unexpected / rare PTMs** (GG, acetyl-K, methyl-K/R, succinyl/crotonyl/malonyl, citrulline, nitro, drug adducts) | Phospho (4 datasets) + fixed/common only | Composition-encoded mods generalise only if many compositions are seen; zero-shot needs held-out classes | Synthetic PXD009449 (21 PTMs) and PXD056559; 176 enrichment candidates (132 usable) | **Highest** |
| 3 | ★ **HLA-II** (12–25-mers, nested sets, higher z) | 0 | Longer non-tryptic, z3–4 | 139 projects; 94/123 usable | Very high |
| 4 | **Non-tryptic enzymes** (GluC, AspN, chymotrypsin, ArgC, LysN, semi) | 0 | C-terminal residue and basic-residue position drive y/b balance and RT | ProteomeTools III AspN/LysN pools; PXD010154; 154/193 usable | High |
| 5 | **timsTOF** (TOF CID with ion mobility) | 0 | Different fragmentation regime (beam CID, collision energy ramps with 1/K0) | 2,467 projects; 99/157 usable; HLA on timsTOF exists | High |
| 6 | **Astral** (ASTMS HCD) | 0 for training (the benchmark is excluded) | The headline benchmark instrument; HCD-ASTMS is a separate activation token | Only 6 COMPLETE DDA; PXD058160 and PXD061848 (TMTpro, mouse); PARTIAL Astral needs our own search | High (scarce) |
| 7 | **ETD / EThcD** | 0 | Needs c/z channels (andes-ms2 v1 has none); catalog now, train later | ETD 529 / EThcD 338 projects (text); 80/115 usable; ProteomeTools ETD pools | Medium now, high once c/z exists |
| 8 | **Ion-trap CID** | 5.4M (1 qpx dataset) + 196k flats | More labs and instruments for the CID-ITMS head (PeptDeep fails on CID; GBDT wins) | ProteomeTools IT-CID (every pool); 565 COMPLETE DDA with IT MS2 | Medium–high |
| 9 | **Q-TOF** (SCIEX, Waters) | 0 | Separate instrument family | 3,720 projects; 70/108 usable | Medium |
| 10 | **TMTpro / iTRAQ** | TMT10/11 only | Different reporter and label masses; RT shift | ProteomeTools TMT (IV, V); 253 + 176 COMPLETE DDA | Medium |
| 11 | **High charge (z ≥ 4) and long peptides** | 1.9% z4+ in flats; the GBDT is weakest here | Oversample; HLA-II and missed cleavages | Inside HLA-II and long-gradient sets | Medium (sampling, not sourcing) |
| 12 | **NCE diversity** | NCE unknown on flats; qpx has 25–35 HCD | NCE conditioning in FiLM | ProteomeTools NCE ladders (20–35 per pool) | Medium |
| 13 | **Non-human organisms** | 0 | Composition generalisation; mostly redundant for RT and MS2 | Most of PRIDE | Low |
| 14 | Glyco (intact) | 0 | Not representable in the current encoder | 1,056 projects | Defer |
| 15 | DIA | – | DDA PSMs only | 5,487 projects | Excluded |

The gaps by dimension, as a matrix:

```
                         current      gap fill (shortlist)
instrument   Orbitrap    ████████     -
             Astral      ·            PXD058160, PXD061848 (+ PARTIAL via search)
             timsTOF     ·            PXD055268/054498/038165/029648/066497/042613/022077
             Q-TOF       ·            PXD014754/025877/000450/002951/001574/057177(Zeno)
activation   HCD-FTMS    ████████     ProteomeTools NCE ladders
             CID-ITMS    ██           ProteomeTools IT pools, PXD043291, PXD057316
             ETD/EThcD   ·            PXD022930/012437/012498/022403/007058/028432
labelling    LFQ, TMT10  ████         TMTpro: PXD023119/023120/058160/057588; iTRAQ: PXD001574
enzyme       trypsin     ████████     PXD010154, PXD021013 (AspN/LysN), PXD013112, PXD034992
peptides     HLA-I/II    ·            PXD021013, PXD019643, MSV000084172, PXD000394, PXD012308 …
PTM          phospho     ██           PXD009449, PXD056559, PXD052972, PXD057347, PXD057177 …
```

---

## 6. Sampling plan (`corpus_sampling.py`)

**1. Unit.** One unit per dataset, split per inferred analyser when its runs differ
(`--runs`). PXD009875 becomes `PXD009875:FTMS` and `PXD009875:ITMS`.

**2. Stratum.** fragmentation × analyser × labelling × enrichment class, where the
enrichment class is HLA-I, HLA-II, PTM, nonspecific or none.

**3. Temperature sampling over strata.** The share is p_s ∝ N_s^α × boost_s, with
α = 0.5. Boosts:

| Class | Boost |
|---|---|
| HLA-I | 4 |
| HLA-II | 3 |
| PTM | 2.5 |
| Non-specific enzyme | 2 |
| Other | 1 |

With α = 0.5, a 26M-PSM HCD-LFQ stratum gets ~√ its natural weight. On the toy
catalog, HLA strata move from 4% to 15% and from 1.4% to 6.7% of draws.

**4. Per-dataset cap.** Within a stratum, no unit contributes more than 25% of the
stratum's draws. Today PXD042233 alone is 26M of 62.8M HCD-LFQ PSMs.

**5. Per-precursor cap** (selection time, `ms2_select.py --max-per-precursor`).
- MS2: at most k = 4 spectra per (peptidoform, charge, activation, NCE bin of 5),
  lowest PEP first.
- RT: one median per (run, peptidoform).
- Cross-source duplicates (`corpus_qc.py overlap`, I/L-merged) count against the same
  cap, so a HeLa peptide seen in 30 datasets cannot dominate.

**6. Activation balance.** Keep `ms2_train.py --balance-activation`. Each activation
token gets a floor share of 5% while the data lasts.

**7. Quality gates before sampling.**
- Drop `prec_fit` = false.
- Drop deamidated PSMs with isotope ≠ 0 (finding 2 in §3).
- MS2: drop `ann_frac` below the dataset's own p5 within each analyser.
- Drop runs with `precursor_join` or `analyser_mixed` flags until resolved.
- Mined-PTM labels (engine outputs with many variable mods) need localisation
  probability ≥ 0.75 where available, and the PTM must be in `mod_composition` or be
  added with an elemental composition. No mass-only labels in training.

**8. Hold-outs.** Whole datasets and whole PTM classes, never PSM-level.

**Datasets:**
- The andes-rt v1 three: PXD019263 (LFQ), PXD005235 (TMT), PXD012255 (phospho).
- One per new stratum (`--holdout-per-stratum 1` picks the second-smallest dataset
  by PSMs). Proposed explicit picks:

| Stratum | Held-out dataset |
|---|---|
| HLA-I | PXD022949 (melanoma, QE HF) |
| HLA-II | PXD029648 (timsTOF) or PXD025877 (TripleTOF) |
| EThcD | PXD012498 |
| Astral | PXD061848 (keep PXD058160 for training) |
| Non-tryptic | PXD013112 |
| PTM | PXD052972 (histone GG) |

**PTM classes** (`--holdout-ptms`). Every PSM carrying the class is test-only in every
dataset, for the zero-shot evaluation the composition encoder is meant to win. Pick
classes with synthetic ground truth in PXD009449/PXD056559 and enough natural data:
- UNIMOD:121 GlyGly;
- UNIMOD:64 succinyl;
- citrullination: UNIMOD:7 on R only, which needs a site-aware filter (the same
  Unimod id as N/Q deamidation).

Report every metric per held-out PTM class: seen-composition vs unseen.

**Benchmarks** (the six above) are dropped, not held out. Datasets in the
model-store ledger may train the predictors but are never used to evaluate a bundled
andes model.

---

## 7. Cluster runbook (when the cluster returns, ~2026-10-20)

Placeholders: `$PY` (Python with pyarrow/pandas), `$S` (this directory), `$MSNET`
(the MsNet collection root, on long-term storage that only compute and data-mover
nodes mount), `$QPX`, `$RAW`, `$OUT`. Use `--comment` tags and log to project
storage, not `$HOME`.

### Step 1: MsNet inventory, then QC (first job)

```bash
sbatch --cpus-per-task=4 --mem=16G --time=4:00:00 -o $OUT/logs/inv_%j.log --wrap \
  "$PY $S/msnet_inventory.py --root $MSNET --out-dir $OUT/msnet_inv \
       --qc-script $S/corpus_qc.py --qc-out $OUT/qc_msnet"
# review $OUT/msnet_inv/datasets.tsv: files, GB, rows, flat-compatible, peaks, columns per PXD
N=$(wc -l < $OUT/msnet_inv/qc_commands.sh)
sbatch --array=1-$N%20 --cpus-per-task=8 --mem=32G --time=6:00:00 -o $OUT/logs/qc_%A_%a.log \
  $S/slurm/run_line.sh $OUT/msnet_inv/qc_commands.sh
for d in $OUT/qc_msnet/*/; do a=$(basename $d)
  $PY $S/corpus_catalog.py msnet --accession $a --qc-dir $d --out-dir $OUT/catalog; done
```

What to read from it:
- Datasets and PXDs, per-file formats and sizes.
- Whether peak lists are attached.
- Instruments and activations: from columns if present, else the QC `analyser_inferred`.
- PSM counts.

The PXD009875 example finished in about 1 min per 0.5M PSMs on 6 cores, so a 50M-PSM
collection is about 1.5 h of array time.

Then:
- Re-run `corpus_sampling.py` and §5 with the new rows.
- Flag any MsNet dataset that duplicates a qpx dataset. Train on one of them.
- Check every MsNet dataset against the model-store ledger.

### Step 2: qpx profiling

```bash
for t in $QPX/*/psm.parquet; do a=$(basename $(dirname $t))
  echo "$PY $S/corpus_qc.py psm_table --dataset $a --out-dir $OUT/qc_qpx/$a $t"; done > $OUT/qpx_cmds.sh
sbatch --array=1-$(wc -l < $OUT/qpx_cmds.sh)%10 --cpus-per-task=4 --mem=64G -o $OUT/logs/qpx_%A_%a.log \
  $S/slurm/run_line.sh $OUT/qpx_cmds.sh
# then: corpus_catalog.py qpx (with --meta from the SDRF/paper for instrument/NCE/labelling)
$PY $S/corpus_qc.py peptidoforms --out $OUT/keys.parquet $OUT/qc_*/*/psm_qc.parquet
$PY $S/corpus_qc.py overlap $OUT/keys.parquet     # cross-source duplicate table
```

qpx writes I as L. `corpus_qc.py` reports `isoleucine_written_as_L`, and every key
merges I and L.

### Step 3: raw-header pass

This resolves activation, NCE and analyser per run:
- `ms2/ms2_run_info.py run.mzML` on one mzML per run;
- or on the TRFP metadata of `$RAW`.

Order:
1. **PXD009875 MudPIT / `_FL_` runs.** CID-ITMS or HCD-ITMS? (§3 finding 1.)
2. The qpx datasets: NCE is missing from the catalog for all of them.

Fold the results into `runs.activation` and `runs.nce`.

### Step 4: PRIDE ID harvest (light; IDs first, spectra later)

```bash
$PY $S/pride_survey.py peek-zip --accession PXD021013 --file <pool>.zip --extract msms.txt --out $OUT/pt/<pool>.msms.txt
$PY $S/corpus_qc.py msms --dataset PXD021013 --out-dir $OUT/qc_pride/PXD021013 $OUT/pt/*.msms.txt
$PY $S/corpus_catalog.py proteometools --accession PXD021013 --qc-dir $OUT/qc_pride/PXD021013 --out-dir $OUT/catalog
```

- Work down `shortlist.tsv`.
- Pull only the identification members: MaxQuant msms.txt by range extraction;
  mzIdentML/mzTab files.
- Profile them first. Download raw files only for datasets that pass QC, and only the
  runs selected by `ms2_select.py`.

Still to build:
- an mzIdentML reader (PSM, score, PEP/q, decoy flag) for `corpus_qc.py`;
- a PD .msf (SQLite) reader.

### Step 5: merge and weights

```bash
$PY $S/corpus_catalog.py merge --out-dir $OUT/catalog_merged $OUT/catalog
$PY $S/corpus_sampling.py --catalog $OUT/catalog_merged/datasets.parquet --runs $OUT/catalog_merged/runs.parquet \
    --out $OUT/weights.tsv --holdout-ptms UNIMOD:121 UNIMOD:64
```

---

## 8. Caveats and next steps

- **Survey strata are text heuristics** (§4.1 precision). Verify the shortlist
  against papers and run headers before training.
- **The usability classes come from file names.** `search_zip?` must be confirmed
  with `peek-zip`. The 11 shortlisted MaxQuant zips listed in §4.3 were confirmed.
- **`--activation` in `ms2_extract_flats.py` is one value per invocation.** Split the
  PXD009875 flats by `run_qc.tsv` before extraction (§3).
- **The `mod_composition` table needs entries** for the target PTM classes:
  - drug haptens (benzylpenicilloyl);
  - ADP-ribose;
  - lactyl, β-hydroxybutyryl, glutaryl;
  - citrulline: currently the same id as deamidation; needs a site-aware feature.
  Without them, PTM-labelled PSMs are dropped by `supported()`.
- **c/z channels** (ETD/EThcD) and **fragment charge ≥ 3** are model work. They
  decide when gaps 7 and 11 become trainable.
