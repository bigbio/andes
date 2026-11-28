# Design: label-free quantification of non-cleavable cross-links

Status: design proposal, not implemented.

Decided: label-free MS1 quantification with match-between-runs comes first. Isotope-labelled
linkers, TMT and DIA come later (§7).

This page is the plan the implementation PRs refer to. Each PR records its measurements in
its own description.

## 1. Goal

Quantify cross-links made with non-cleavable linkers (DSS/BS3, BS2G, EGS, sulfo-SDA,
EDC/DMTMM, PhoX, …) across runs and conditions. Report quantities and changes at four
levels: cross-link precursor, cross-linked peptide pair, residue pair and protein pair (PPI).

Error is controlled at every step:

- the identification FDR of each level;
- a measured false-transfer rate for match-between-runs;
- the false-discovery proportion (FDP) of differential calls.

MS1 quantification does not depend on how the linker fragments, so cleavable linkers get the
same label-free path.

## 2. Where andes stands

### Cross-link search

Code: `crates/xl`, `crates/search/src/crosslinks.rs`, `crates/andes/src/bin/andes/xl.rs` and
`xl/`.

- **Search.** Non-cleavable linkers are searched by open-mass retrieval on a second fragment
  index. The partner peptide is found by the remaining precursor mass, at isotope offsets
  0..2.
- **FDR.** Target-target, target-decoy and decoy-decoy counts, split into intra, inter and
  homotypic links and computed per level (`crates/xl/src/fdr.rs`, `rollup.rs`). This is the
  xiFDR scheme. Residue pairs and protein pairs are rolled up by maximum score. Protein-pair
  q-values exist for inter links only.
- **Missing for quantification.** A cross-link spectrum match (`Csm`, `XlSpectrumResult` in
  `crates/search/src/crosslinks.rs`) carries:
  - no retention time;
  - no precursor m/z, only the neutral mass, the charge and the isotope offset;
  - no elemental composition.

  `Crosslinker` keeps only the intact mass, although `resources/crosslinkers.yaml` has the
  compositions.
- **Linear search skipped.** `--xl` runs without the linear search
  (`crates/andes/src/bin/andes/search/files.rs`, `search/finish.rs`).
- **Quant flags rejected.** `--xl` rejects `--lfq`, `--mbr`, `--quant-dir`, `--tmt` and
  `--label` (`crates/andes/src/bin/andes/search/resolve/checks.rs`).
- **Cohort.** Cross-link runs get only the `xl-fdr` stage, and their quantification setting
  is cleared (`crates/andes/src/bin/andes/cohort.rs`). Outside the isobaric `RunDesign`
  there is no notion of sample or condition.

### Quantification stack

`crates/quant` works on any precursor; nothing in it assumes a linear peptide:

- `LfqTarget` (`lfq.rs`) holds a monoisotopic m/z, a charge, an anchor retention time and an
  isotope envelope.
- `evidence.rs` measures with a Gaussian template, with M−1 and interference checks.
- `align.rs` fits a monotone piecewise-linear alignment.
- `decoys.rs` makes the transfer decoy probes: mass shift, retention-time shift and random
  retention time. The binary scores transfers with a bagged gradient-boosted model and
  k-decoy q-values (`crates/andes/src/bin/andes/transfers.rs`).
- `claims.rs` detects a peak claimed by another precursor; `normalize.rs` normalises
  RT-locally.

The linear-peptide assumptions live in the binary: peptidoform keys, the Percolator join and
the b/y ladders used for MS2 support. Glyco quantification
(`crates/andes/src/bin/andes/glyco/quantities.rs` and `quantities/`) already reuses
`quant::lfq` for a precursor that is not a linear peptide, with its own targets, transfers and
tables. It is the template to follow.

Gaps:

- `quant::formula::Formula` has C, H, N, O and S only. PhoX/DSPP (phosphorus) cannot be
  modelled.
- Protein roll-up is total intensity or Top3. There is no MaxLFQ and no differential testing
  in Rust.

DOCS.md §10 leaves normalisation, roll-up and statistics to quantms/MSstats. Those tools have
no notion of residue pairs or protein pairs, so cross-link roll-up has to live in andes. Glyco's
site and protein tables are the same kind of exception.

## 3. The tools to compete with

| Workflow | Licence | Quantification | Transfer error control | Roll-up | Statistics |
|---|---|---|---|---|---|
| xiSEARCH + xiFDR, then Skyline (MS1 filtering) | Apache-2.0 | label-free MS1, d0/d4 | none; peak boundaries corrected by hand | by hand or spreadsheet | none |
| xiSEARCH or MS Annika, then Spectronaut (DIA) | commercial quantifier | DIA MS1/MS2 | Spectronaut q-values; cross-link decoys and entrapment in Birklbauer et al., bioRxiv 2026 | residue pair by top-N or mean | ANOVA + Benjamini-Hochberg |
| xTract (xQuest/xProphet) | open Perl, unmaintained | label-free, DSS-d0/d12 | mProphet decoy peak groups | sum per cross-link | t-test + Benjamini-Hochberg |
| MaxQuant 1.5.4 quantification-only mode | free, closed | BS3-d0/d4 only, no match-between-runs | none | median per residue pair | none |
| MS Annika 3 + apQuant (Proteome Discoverer) | freeware, closed | label-free MS1 | not documented | spectrum match, residue pair | Proteome Discoverer |
| MaxLynx | free, closed | none: identification only | — | — | — |
| OpenPepXL (OpenMS) | BSD-3 | none (OpenPepXLLF was removed in OpenMS 3.5) | — | — | — |
| iqPIR | BSD-3 | isobaric; cleavable PIR linkers only | identification FDR only | cross-link and protein pair, monolink-normalised | t-test, Bonferroni |

Sources:

- Skyline: its documentation of cross-linked peptide support (skyline.ms).
- DIA cross-link quantification: Müller et al., Mol Cell Proteomics 2019; Birklbauer et al.,
  bioRxiv 2026.
- xTract: Walzthoeni et al., Nat Methods 2015.
- MaxQuant mode: Chen et al., Mol Cell Proteomics 2016.
- MaxLynx: Yılmaz et al., Anal Chem 2022.
- iqPIR: Chavez et al., Anal Chem 2021.

Full references are in §9.

What follows:

1. **No open tool does the whole job.** None identifies, quantifies, transfers between runs
   with a controlled error, rolls up to protein pairs and tests for differences in one run.
2. **Transfer error is unmeasured.** No cross-link match-between-runs false-transfer rate, and
   no FDP for differential cross-links, has been published. andes already controls transfer
   error for linear peptides in the manner of IonQuant (Yu et al., Mol Cell Proteomics 2021).
   Measuring and controlling both for cross-links is the clearest way to come out ahead.
3. **Abundance and interaction change are not separated.** Only iqPIR separates
   protein-abundance change from cross-link change, by monolink normalisation and with its
   own cleavable linkers. andes can do it label-free from the linear peptides of the same
   runs.

## 4. Design (label-free)

### D1. Linear co-search in quantified cross-link runs

When quantification is on, an `--xl` run also searches the same spectra for linear and
monolinked peptides. The path exists already for `--xl-linear-only`. The co-search provides:

- **Alignment anchors.** Thousands of them, where cross-links alone give hundreds.
- **Normalisation** from the same injection.
- **Protein abundances** for D7.
- **Linear precursors in the claims library**, so a cross-link envelope that sits on a
  co-eluting linear peptide's isotopes is detected (D4).

Its cost relative to the cross-link search is measured in P0. The cross-link search's own
time and memory are not documented today.

### D2. Precursor identity

A cross-link precursor is keyed by
`(linker, α peptidoform, α site, β peptidoform, β site, charge)`, in a canonical order.

- **Site isomers.** Two link sites on the same peptide pair are isobaric and can co-elute.
  They are measured at `(linker, α peptidoform, β peptidoform, charge)` and attributed to the
  site of the best spectrum match, with a `site_ambiguous` flag, as `IsomerClass` does for
  linear peptides.
- **Conflicts between runs.** When runs assign one precursor to different peptide pairs, the
  best-scoring assignment wins, and the conflicts are counted in the log.

### D3. Theoretical envelope

- **Formula.** α formula + β formula + the linker composition from `crosslinkers.yaml`. The
  composition is kept on `Crosslinker`, and `Formula` gains phosphorus.
- **Consistency test.** A chemistry test checks that every reagent's composition mass equals
  the intact mass the search uses.
- **Number of isotopes** is chosen by envelope coverage, as glyco does (4 to 14 isotopes
  covering 95 %). A fixed four is not enough for 3–6 kDa precursors at charge 3 to 7.
- **Monoisotope checks.** The fixed `min_mono_share` of `LfqParams` (0.3) cannot apply: at
  these masses the monoisotopic peak is a small share of the envelope. Monoisotope errors
  are caught by the envelope cosine and the M−1 trace of `evidence.rs` instead.
- **m/z used.** The search accepts isotope offsets 0..2, so quantification uses the
  theoretical monoisotopic m/z, never the reported precursor m/z.

### D4. Per-run measurement

- **Targets.** Every accepted cross-link precursor, measured with `lfq::quantify_at` or the
  `evidence.rs` measurement and the D3 envelope.
- **Per-run q-values.** Each run has its own decoys (m/z-shifted traces) and picked
  target-decoy q-values.
- **Interference.** Linear and cross-link precursors form one library for `claims.rs`.

### D5. Match between runs

- **Library.** The cross-link precursors whose residue pair passes the identification
  threshold, on the pooled cohort q-values, each with a consensus retention time.
- **Alignment.** `align.rs` on the joint anchors (linear and cross-link). Cross-links follow
  the maps the linear peptides define.
- **Transfers.**
  - Measurement: `decoys.rs` probes and the `evidence.rs` measurement with the cross-link
    envelope.
  - Model: the transfer model is trained together with the linear transfers.
  - q-values: stratified by precursor class (cross-link or linear), because cross-link
    precursors are a different population, of higher charge and lower intensity.
  - Fallback: once a run has enough decoys, a cross-link-only competition takes over.
- **Reporting.** Transfers are reported apart from identifications and gated by
  `--quant-transfer-fdr`.

### D6. Roll-up

- **MaxLFQ.** A shared, deterministic `quant::rollup` module, after Cox et al., Mol Cell
  Proteomics 2014: pairwise median log-ratios between runs, solved by least squares. It
  rolls precursors up to residue pairs, and residue pairs up to protein pairs.
- **Top-N baseline.** A top-N value is computed too, as the residue-pair quantity of
  Müller et al. (Mol Cell Proteomics 2019), who compared N from 1 to 10. The benchmark picks
  the default between MaxLFQ and top-N, which is then an algorithm constant.
- **Protein pairs.** A protein-pair value comes from its inter-protein residue pairs. It
  carries `n_residue_pairs` and a `heterogeneity` value, the spread of the residue-pair
  log-ratios. A high spread flags interfaces that change differently, which points to a
  change of conformation rather than of abundance.
- **What is quantified.** Only entities that pass their level's identification FDR. Every
  value records how many of its precursors were identified and how many transferred.

### D7. Abundance-aware change

Protein abundances come from the linear co-search. Each residue pair reports:

- `log2fc_xl`, `log2fc_protein_a` and `log2fc_protein_b`;
- an adjusted change:
  - intra link: `xl − A`, a change of conformation or accessibility;
  - inter link: `xl − mean(A, B)`.

Monolinks on the same lysine are an optional normaliser, as in iqPIR; the water monolink is
already searched.

This adjustment is andes's own proposal, not taken from a paper. Above all the inter-link
model is a modelling choice, so the adjusted change is labelled experimental until the
benchmark designs show it helps.

### D8. Statistics and conditions

andes writes an MSstats-format table: residue pair (or protein pair) as `ProteinName`, the
cross-link precursor as `PeptideSequence`. MSstats then handles missing values and testing
unchanged. This matches DOCS.md §10.

Native testing is not part of the first phases. The benchmark harness measures the FDP of
differential calls with the shared evaluation in `docs/benchmarks/quant-eval` (median polish,
moderated t-test, Benjamini-Hochberg). SDRF `factor value[...]` columns are passed through to
a `condition` column, which needs a small addition to the cohort's planned run.

### D9. Outputs

New tables are Parquet with a `TableVersion`, written atomically.

- **`<base>.xl.csms.parquet`**, the input contract of quantification, as `RunPsmRow` is for
  linear peptides. Columns: run, scan, `rt_s`, the precursor m/z after the MS1 snap, charge,
  isotope offset, theoretical m/z, the keys, and the spectrum-match, residue-pair and
  protein-pair q-values.
- **`<base>.xl.features.parquet`**, one row per precursor and run. Columns: evidence
  (`identified` or `transfer_ms1`), apex, area, `quantity_normalized`, envelope cosine,
  feature q-value and `quant_q`.
- **`<base>.xl.residue_pairs.quant.parquet`** and **`<base>.xl.ppi.quant.parquet`**, long
  format. Columns: entity, run, condition, quantity (MaxLFQ and top-N), `n_precursors`,
  `n_identified`, `n_transferred`, the protein log2 fold changes and the adjusted change.
- **Wide TSV matrices**, as glyco's `.matrix.tsv` already does, and `<base>.xl.msstats.csv`.
- **A `statistics.log` section:**
  - quantified counts per level, identified against transferred, share missing, CV per
    condition;
  - every inferred value (retention-time shift, tolerance), with the amount of evidence
    behind it.

### D10. Command line and cohort

- **Flags allowed with `--xl`.** `checks.rs` stops rejecting `--lfq`, `--mbr`, `--quant-dir`,
  `--quant-fdr` and `--quant-transfer-fdr`.
  - `--xl … --lfq` quantifies per run, with the linear co-search.
  - `--mbr` adds match between runs.
- **New parameter: `--xl-fdr`.** The identification threshold of the quantification library
  encodes the user's intent, so it is a parameter: residue-pair level, default 0.01, as in
  the cross-link benchmark. It is the only new flag.
- **Staged path.** The staged `andes quant` path accepts cross-link rows.
- **Cohort.** Cross-link runs get `xl-fdr`, then `xl-quant`, when their SDRF label is
  `label free sample` or the configuration turns on label-free quantification.

## 5. Phases

There is one reviewable PR per phase, and each is measured on the benchmark before merge.

| # | Change | Results changed? | Measured on |
|---|---|---|---|
| P0 | Benchmark harness: fetch, truth and evaluation for §6; a scripted xiSEARCH + xiFDR + Skyline reference; results JSON and a `render.py` block; cross-link search time and memory. | no | — |
| P1 | Retention time, MS1-snapped precursor m/z, run and theoretical m/z on spectrum matches; `.xl.csms.parquet`; linker composition on `Crosslinker`; phosphorus in `Formula`; the composition test. | no; columns added | tests |
| P2 | Linear co-search in quantified runs; per-run cross-link quantification (D3, D4); `.xl.features.parquet`; `--xl --lfq` allowed. | new outputs | replicate CV and completeness |
| P3 | Match between runs (D5): joint library and alignment, class-stratified transfer q-values, `--mbr`, the staged `andes quant` path, the cohort `xl-quant` stage. | new outputs | false-transfer rate (§6) |
| P4 | Roll-up and reporting (D6 to D9): `quant::rollup`, residue-pair and protein-pair tables, protein fold changes, MSstats export, SDRF condition, `statistics.log`, DOCS.md. | new outputs | dilution accuracy, FDP of differential calls |
| P5 | `docs/benchmarks/xl/quant.md`, README row, `reproduce.md` steps. | docs | — |

## 6. Benchmarks and acceptance targets

The rules are those of `docs/dev/fdr-validation.md` and
`docs/dev/benchmark-reproducibility.md`:

- several datasets, with tuning and held-out sets kept apart;
- seed ranges reported;
- full commands, versions, wall time and peak memory recorded;
- negative results kept.

| Dataset | Linker | Design | Measures | Role |
|---|---|---|---|---|
| PXD011036 | BS3 | 7 proteins, dilution 0.1–1× in triplicate; *E. coli* matrix series | ratio accuracy and compression, CV, share quantified, true positives of differential calls | tuning |
| PXD007250 | BS3 | HSA, 10 reaction and 10 injection replicates (every ratio is 1) | CV, completeness; every differential call is false | held out |
| Clasen batches (in `docs/benchmarks/xl`) | DSSO | recombinant interaction groups defined per batch | false transfers: a cross-link transferred into a run of a batch whose interaction groups do not contain its protein pair | false-transfer truth |
| Matzinger library and library in HEK (in `docs/benchmarks/xl`) | DSSO, DSBU | library alone and library spiked into HEK | false transfers: a HEK cross-link transferred into a library-only run | false-transfer truth |
| PXD059096 | PhoX | Cas9 dilution, Orbitrap Astral | accuracy and CV, once the data are public; needs phosphorus in `Formula` | held out |

Published figures to compare against:

- PXD011036: 292 of 414 residue pairs quantified; replicate CV 42–76 % by DDA and 10–15 %
  by DIA.
- PXD007250: 106 of 196 residue pairs complete across replicates; median CV 14 % for
  injection and 32 % for reaction replicates.

The false-transfer designs use cleavable linkers, but an MS1 transfer does not depend on the
linker's chemistry. They would give the first measured false-transfer rate for cross-link
match-between-runs.

No public two-proteome spike-in exists for a non-cleavable linker (for example
DSS-crosslinked *E. coli* spiked 1:1, 1:2 and 1:4 into HEK). It would be the decisive
benchmark, and it needs a wet lab.

Targets, to be measured and not promised:

- **Coverage.** At least as many quantified residue pairs as xiSEARCH + xiFDR + Skyline, at
  the same residue-pair FDR, on every dataset.
- **Precision.** Median CV no higher than the Skyline workflow on PXD007250.
- **Accuracy.** Median absolute log2 error of 0.2 or less on PXD011036 at 0.3–1×, with the
  compression at 0.1× reported.
- **Transfer error.** Measured false-transfer rate at or below `--quant-transfer-fdr` on the
  false-transfer designs. Transfers are always counted separately.
- **Differential calls.** FDP at or below the nominal 5 % on PXD007250 and on the known
  designs.

## 7. Later

- **Isotope-labelled linkers (BS3-d4, DSS-d12, BS2G-d4).**
  - Reagents: YAML entries with a ²H composition and `label_of`. A light reagent listed
    together with its label (`--xl "BS3; BS3-d4"`) is searched as two channels of one
    reagent, not as competitors.
  - Quantification: the `label/fit.rs` joint fit (`GroupProbe`, `ChannelEnvelope` with 4 or
    12 ²H). `learn_rt_shifts` learns the deuterium shift from the data.
  - Why a joint fit: with d4 the light M+4 isotope sits on the heavy monoisotope, which is
    where separate integration lost accuracy (Fischer et al., J Proteomics 2013; Chen et
    al., Mol Cell Proteomics 2016).
  - Benchmarks: PXD004107 (HSA, d0:d4 at 1:1, 1:2, 1:4) and PXD001675 (C3 against C3b,
    curated residue pairs).
- **MS2-supported transfers.** An MS2 scan at a transfer's apex scored with the cross-link
  pair scorer.
- **TMT reporters on cross-link spectra**, reusing `isobaric.rs`.
- **DIA**, with a library from cross-link identifications or predicted spectra (Prosit-XL,
  Kalhor et al., Nat Commun 2025), and the cross-link decoy design of Birklbauer et al.,
  bioRxiv 2026.

## 8. Risks

- **Alignment and normalisation anchors too sparse:** handled by the linear co-search (D1).
- **Global normalisation can be wrong.** In a pure dilution series, normalising on all
  peptides removes the signal. Raw and normalised quantities are both written, and each
  evaluation uses the one its design calls for.
- **Cross-link transfers miscalibrated** by a model trained mostly on linear transfers:
  handled by class-stratified q-values, a cross-link-only competition when there are enough
  decoys, and the measured false-transfer rate.
- **High charge and chimeric spectra** cause monoisotope errors and interference: handled by
  the coverage-based envelope, the M−1 check and joint claims.
- **Constants fitted to small benchmarks** (seven proteins, HSA): handled by keeping tuning
  and held-out sets apart.
- **Pooled cohort q-values do not control each project's share** (DOCS.md, "Crosslink FDR
  over a cohort"). The library uses pooled q-values, and that is documented.

## 9. References

- Birklbauer et al., bioRxiv 2026. A DIA-based quantitative cross-linking MS framework, with
  a four-state target/decoy library and entrapment.
  [doi:10.64898/2026.06.16.732730](https://doi.org/10.64898/2026.06.16.732730)
- Chavez et al., Anal Chem 2021. Six-plex isobaric quantitative protein interaction reporter
  (iqPIR). [doi:10.1021/acs.analchem.1c02209](https://doi.org/10.1021/acs.analchem.1c02209)
- Chen et al., Mol Cell Proteomics 2016. Quantitative cross-linking/MS with isotope-labelled
  cross-linkers and MaxQuant.
  [doi:10.1074/mcp.M115.056481](https://doi.org/10.1074/mcp.M115.056481)
- Cox et al., Mol Cell Proteomics 2014. MaxLFQ: delayed normalisation and maximal peptide
  ratio extraction. [doi:10.1074/mcp.M113.031591](https://doi.org/10.1074/mcp.M113.031591)
- Fischer et al., J Proteomics 2013. Quantitative cross-linking/MS with isotope-labelled
  cross-linkers (XiQ). [doi:10.1016/j.jprot.2013.03.005](https://doi.org/10.1016/j.jprot.2013.03.005)
- Kalhor et al., Nat Commun 2025. Prosit-XL: fragment-intensity prediction for cross-linked
  peptides.
  [doi:10.1038/s41467-025-61203-4](https://doi.org/10.1038/s41467-025-61203-4)
- Lenz et al., Nat Commun 2021. Protein-protein interaction FDR in cross-linking MS (xiFDR).
  [doi:10.1038/s41467-021-23666-z](https://doi.org/10.1038/s41467-021-23666-z)
- Müller et al., J Am Soc Mass Spectrom 2018. Reproducibility of label-free quantitative
  cross-linking/MS. [doi:10.1007/s13361-017-1837-2](https://doi.org/10.1007/s13361-017-1837-2)
- Müller et al., Mol Cell Proteomics 2019. DIA improves quantitative cross-linking/MS.
  [doi:10.1074/mcp.TIR118.001276](https://doi.org/10.1074/mcp.TIR118.001276)
- Walzthoeni et al., Nat Methods 2015. xTract: quantitative cross-linking MS of conformational
  changes. [doi:10.1038/nmeth.3631](https://doi.org/10.1038/nmeth.3631)
- Yılmaz et al., Anal Chem 2022. MaxLynx: cross-link identification in MaxQuant.
  [doi:10.1021/acs.analchem.1c03688](https://doi.org/10.1021/acs.analchem.1c03688)
- Yu et al., Mol Cell Proteomics 2021. IonQuant match-between-runs with FDR control.
  [doi:10.1016/j.mcpro.2021.100077](https://doi.org/10.1016/j.mcpro.2021.100077)

The method references go into `docs/citations.yaml`, with citation markers, in the PRs that
implement those methods.
