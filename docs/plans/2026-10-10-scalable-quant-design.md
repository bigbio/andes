# Scalable LFQ and TMT quantification with match-between-runs: design

Status: proposed, 2026-10-10. Builds on PR #119 (`--tmt`, `--lfq`).

## 1. Goal

Make andes quantification more accurate and faster than the current DDA quantification
tools, for label-free (LFQ) and isobaric (TMT) data, at any number of runs. The design
keeps one tool (andes), but splits the work into stages: per-run stages run in parallel
(one quantms task per raw file), and global stages run once on compact tables.

Success means, on public ground-truth data:

- LFQ: fewer missing values than the best DDA tools at a measured, controlled
  transfer FDR, with fold-change accuracy and CVs at least as good.
- TMT: less ratio compression than purity filtering alone, measured on the knockout
  channels of the TKO standard.
- Speed: one raw-file parse per run, and global-stage memory that grows with the
  library, not with the number of runs.

### Out of scope

Protein inference, protein roll-up, normalisation and statistics stay in quantms. Reported
identification FDR stays PSM-level and comes from Percolator.

## 2. Starting point (PR #119, measured 2026-10-10)

- TMT extraction and impurity correction reproduce OpenMS IsobaricAnalyzer exactly.
  On TKO (PXD017803, TMTpro16), the median knockout/non-knockout ratio after correction
  is 0.210 (hrMS2) and 0.127 (SPS-MS3). With andes purity at least 0.9 it drops to 0.040.
- LFQ on UPS1 in yeast (PXD001819, 6 runs): yeast median CV 3.9% and fold-change
  median absolute error 0.16. But 32% of yeast and 51% of UPS1 run cells are missing,
  because there is no match-between-runs; a comparison tool with MBR has 2-3% missing.

The per-run measurement is precise. What is missing is a value in every run.

What the PR #119 per-run output lacks for this design:

- the MS1 cache is a temporary file
- only identified precursors are extracted, and no decoy features are written
- the feature score and feature q-value are run-local
- no isotope traces, theoretical envelope, FAIMS voltage, isolation-window contents
  or reporter noise are written

## 3. Principle

Separate quantification tools see only an identification table plus MS1. andes runs the
search, so it can treat every (precursor, run) pair as one question: is the precursor
present in this run, how much is there, and how confident are we? It answers that with
all the evidence at once: MS1 features, MS2 fragments, chimeric co-identifications and
agreement across runs.

## 4. Pipeline

```
 PER RUN (parallel)          GLOBAL (once, m/z-shardable)      PER RUN (parallel)        GLOBAL / PER RUN
 andes search            ->  Percolator (merged, quantms)  ->  andes quant extract   ->  andes quant finalize
   PIN / idparquet           andes quant library                 reads run.cache only      LFQ: transfer FDR,
   run.cache                   global confidence                 MS1 traces + MS2 check         signal selection
   run.quant.parquet           joint RT alignment                consensus boundaries      TMT: reporter unmixing
                               library + decoys                                            -> QPX feature parquet
```

`andes run <files...>` runs all stages in one process on one machine, using the same code
and the same intermediate files. A laptop run and a 1,000-task quantms run give the same
result.

## 5. Stage 1: `andes search` (per run)

The search is unchanged. While spectra are in memory, it writes two more files.

### 5.1 `run.cache`

- Centroided MS1, indexed by retention time, one series per FAIMS compensation voltage.
- The reduced MS2 peaks that the search already computes, with each scan's isolation
  window, survey-scan reference, acquired m/z range and FAIMS voltage.
- Peaks stored as f32 m/z and f32 intensity in RT-ordered blocks, zstd-compressed, with
  a block index for random access by RT range.

No later stage opens the raw file again.

### 5.2 `run.quant.parquet` (versioned schema)

- One row per rank-1 PSM, targets and decoys, with no FDR filter.
- For each row, the MS1 feature: apex RT, boundaries, per-isotope areas, envelope cosine,
  mass error, FAIMS voltage, and the short isotope traces themselves (4 isotopes, about
  30 points each).
- The theoretical isotope envelope.
- RT anchors: PSM id, RT and score.
- TMT: per MS2 or MS3 scan, the reporter intensities, the reporter noise where the
  `.raw` provides it, and the isolation-window contents. The window contents are the MS1
  isotope clusters inside the window with their intensity fractions, interpolated between
  the flanking MS1 scans, plus any chimeric co-identifications in that scan.

The schema carries a version field, and readers reject unknown major versions.

## 6. Stage 2: `andes quant library` (global)

### 6.1 Global confidence

A precursor enters the library based on Percolator's peptide-level q-value from one
Percolator run over all files. This is an internal entry criterion, not a reported
result. OR-ing per-run PSM q-values across N runs would let false library entries grow
with N, and each one is transferred into every run.

### 6.2 Joint RT alignment

It fits one consensus RT per precursor and one monotone spline per run, using robust
(Huber) least squares and alternating between the two. There is no reference run, so no
reference bias, and the cost grows linearly with the number of runs. It also fits each
run's residual width as a function of RT; this sets the extraction windows.

### 6.3 Library entries

- precursor and charge
- consensus RT and consensus integration boundaries, on the aligned scale
- the observed isotope envelope
- the best fragment spectrum, from the best PSM

### 6.4 MBR decoys

Two candidate families: RT-displaced (outside the precursor's window) and m/z-shifted
(offsets that avoid isotope spacings). The entrapment benchmark (section 10) decides
which family is calibrated; it is not assumed.

### 6.5 Sharding

The library and the per-run target lists are partitioned by m/z, so this stage and the
next can be split across tasks.

## 7. Stage 3: `andes quant extract` (per run)

It covers every library precursor in that run, identified there or not, so each
(precursor, run) pair gets one consistent measurement.

### 7.1 MS1

A single sweep along RT over the cached MS1 extracts all active targets at once
(targets sorted by m/z, binary search per scan). It reports:

- isotope traces within the consensus boundaries mapped into the run
- envelope cosine and mass error
- RT deviation divided by the local residual width
- how well the isotope traces correlate in shape

### 7.2 MS2 check

Every MS2 scan whose isolation window covers the precursor within its elution peak is
scored against the library spectrum, with chimeric scoring allowed. Only scans inside
transfer windows are touched.

### 7.3 Output

`run.extract.parquet`: one row per (precursor, run) with the MS1 measures, the MS2
evidence and the traces.

## 8. Stage 4: `andes quant finalize`

### 8.1 LFQ transfer FDR

It reuses the native rescorer: a single-pass GBDT with 3-fold cross-validation. Folds are
grouped by precursor, so a precursor's runs never span folds. It trains target versus
decoy transfers and outputs a q-value and posterior error probability per transfer.
Every value carries an evidence type: identified, MS2-backed transfer, or MS1-only
transfer.

### 8.2 Signal selection

Per precursor, it chooses which isotopes and RT sub-ranges to trust from how
consistently the signals agree with each other across runs. It never uses how
consistent the totals are across runs. A real fold change scales every signal of a
precursor together, while interference makes one signal's ratio to its siblings drift.
Ranking by between-run variance of totals would quietly suppress real fold changes.

Output: the quantity plus a variance estimate, so quantms can weight precursors.

### 8.3 TMT reporter unmixing (per run)

Each scan's reporter vector is modelled as `y_s = sum_k a_sk * x_k + b_s`:

- `a_sk` is the fraction of isolated ion current from co-isolated peptide k, taken from
  the isolation-window record
- `x_k` is peptide k's unknown reporter profile
- `b_s` is a background term for unidentified interference

Every peptide appears both in its own scans and as an interferer in others. Solving all
scans together by non-negative least squares, one group of scans linked by shared
peptides at a time, recovers clean profiles. A per-run scale factor absorbs differences
in fragmentation efficiency. Reporter S/N weighting is applied when the `.raw` provides
noise; mzML falls back to intensity.

### 8.4 Output

The QPX feature parquet (OpenMS QPXFeatureSchema), extended with evidence type, transfer
q-value, posterior error probability and variance.

## 9. Speed

- One raw-file parse per run.
- Extraction is one RT sweep per run with m/z-sorted targets, not one chromatogram
  lookup per target.
- MS2 checks only touch scans inside transfer windows.
- Every global stage streams column-pruned Parquet and can be sharded by m/z.

Throughput is benchmarked against a fast integrated search-plus-quant engine, and
match-between-runs runtime against the established DDA MBR tools, at 10 and at 100+ runs.

## 10. Validation

- LFQ accuracy and missing values: ProteoBench DDA LFQ ion-level (HYE), plus UPS1 in yeast.
- Transfer FDR: an entrapment design in which one species is absent from some runs.
  Every transfer of that species into those runs is false, giving the true transfer FDP.
- TMT: TKO knockout-channel signal (interference), plus a known-ratio TMT mixture when
  one is available.
- Every change is a one-variable A/B on the same inputs.

## 11. Phases (one PR each; a phase that does not win is dropped)

| Phase | Delivers | Gate |
|---|---|---|
| P0 | `run.cache`, versioned `run.quant.parquet`, decoy features written, global library confidence | PR #119 numbers reproduce with no raw re-read |
| P1 | MS1-only match-between-runs: library, joint alignment, extract, transfer FDR | missing values near the best DDA MBR tools; true transfer FDP at or below nominal on the entrapment design |
| P2 | MS2-checked transfers | more transfers at the same true FDP |
| P3 | signal selection with variance | better CV and fold-change accuracy; no fold-change compression |
| P4 | TMT reporter S/N and unmixing | lower TKO knockout-channel signal than purity filtering at matched PSM counts |
| P5 | quantms wiring (scatter, gather, scatter, gather) | same results as `andes run` on the same files |

## 12. Risks and open points

- **Cache size.** About 290 MB per Velos run for MS1 alone in PR #119's cache; the
  compressed block format must be measured on Astral before P0 is accepted.
- **Library FDR at very large N.** One merged Percolator run over 1,000 files is heavy;
  Percolator's subset training is the expected answer, to be measured.
- **Decoy realism for MBR.** If neither decoy family is calibrated on the entrapment
  design, P1 does not ship.
- **Interferer contribution.** In TMT unmixing, MS1 intensity is only a proxy for an
  interferer's contribution to the reporters; the per-run scale factor may not be enough
  on its own.
- **S/N on mzML.** Converted mzML usually has no noise arrays, so S/N weighting is
  available only on `.raw` input.
- **Thermo paths.** The `.raw` paths added in PR #119 are compile-checked only; runtime
  verification needs a host with the .NET runtime.
