# andes benchmarks

andes is benchmarked on public PRIDE datasets for identification (high-res, TMT-labelled and
low-res LFQ, phospho, glyco) and for quantification (label-free with match between runs, and
TMT). Identifications are counted after Percolator at `q ≤ 0.01`. Where the database carries an
entrapment component, the page also gives the **true** false discovery proportion (FDP), which
is what tells you whether the nominal 1% holds. Other engines are anonymised (engine A, B, ...).
Their numbers come either from runs on the same host with the same rescoring, or from a public
leaderboard. Section 2 gives the date and the andes commit behind each row.

## 1. Summary

**Identification** (PSMs at Percolator `q ≤ 0.01` unless stated; 8 threads, one host):

| benchmark | dataset | andes | best comparison | true FDP at nominal 1% | andes wall |
|---|---|---:|---|---|---:|
| Standard, high-res | Astral HCD LFQ, PXD070049 (1 file) | **47,080** PSMs | engine K 31,435 (217 s) | 0.95–1.03% | 160 s |
| TMT-labelled, low-res | a05058, PXD007683 (1 file, ion-trap CID) | **12,347** PSMs | engine K 10,504 (77 s) | not measurable | 28 s |
| Low-res LFQ | UPS1, PXD001819 (1 file, ion-trap CID) | **15,848** PSMs | engine K 14,734 (42 s); the Java predecessor 15,904 | **2.3–2.5%** (above nominal; 2026-10-20, see below) | 18 s |
| Phospho-enriched | mouse liver, PXD007653 (1 file, Q Exactive) | **37,190** PSMs, ~26,800 phospho-bearing | engine K 33,984 at 1.77% FDP (1,096 s) | 0.97–1.08% | 637 s |
| Glyco, quick tier | mouse liver, PXD005553 (1 fraction) | **6,922–6,959** glycoPSMs | 89.0% of engine L's and 89.6% of engine A's identifications confirmed | 0.85–1.36% | 6 min |
| Glyco, deep tier | mouse liver, PXD005553 (5 fractions) | **31,666 ± 9** glycoPSMs | 78.9% (engine L) and 88.0% (engine A) confirmed | 1.11% ± 0.03 | 23–29 min per fraction, 16 threads |

**Quantification:**

| benchmark | dataset | andes | best comparison | error control | andes cost |
|---|---|---|---|---|---|
| Label-free, match between runs (ProteoBench "Quant LFQ DDA ion") | PXD028735, 6 runs, Q Exactive HF-X, human/yeast/E. coli mix | **54,932** precursors in ≥ 3 runs; median \|ε\| **0.155**, CV **0.110** (transfer q 0.01) | median \|ε\| 0.190 (engine A, submission 2); CV 0.127 (engine C, submission 7); most precursors 81,933 (engine E, median \|ε\| 0.388) | PSM FDP 0.97%; entrapment-only quantified precursors 0.26% | 8.6k CPU-s for 6 runs |
| TMT reporter ions | yeast triple knockout, TMTpro 16-plex, PXD017803 | knockout-channel interference **0.127** (SPS-MS3), **0.040** at purity ≥ 0.9 | 0.139 with a reference purity filter at ≥ 0.9 | — | — |

**Opt-in modes** (Astral, PSMs / distinct peptides, Percolator seeds 42, 1 and 2):

| mode | andes | without the mode |
|---|---|---|
| `--chimeric` (second and later peptides per spectrum), at an entrapment true FDP of 1% | 66,377–66,424 / 29,384–29,400 | 43,076–43,145 / 26,723–26,746 |
| `--refine` (PTM discovery, high-res only) | 48,312 PSMs, 28,728 peptides at `q ≤ 0.01`; its second-pass PSMs at ~1.2–1.8% FDP | 47,080 / 28,932 |

## 2. The benchmarks

**Standard, high-res (Astral).** One 15-minute DDA run of the ProteoBench human/yeast/E. coli
sample against the 31,889-entry HYE FASTA, measured 2026-10-07 at `40774aca`. The true FDP comes
from a 1:1 shuffled entrapment version of the same database. Engine K was run on the same host
with matched parameters and the same Percolator rescoring.

**TMT-labelled, low-res.** One ion-trap CID run of a TMT 10-plex sample against human + yeast
reviewed proteins (26,483 sequences), measured 2026-10-10 at `6a61f61`. The database has no
entrapment component, so only Percolator's own q-value is available.

**Low-res LFQ (UPS1).** One run of the UPS1 standard in yeast, searched against yeast plus an
E. coli entrapment proteome (PSM count 2026-10-10, `6a61f61`). Error re-measured 2026-10-20 at
`8c587fd` from the mzML file, Percolator seeds 42, 1 and 2: 16,446–16,491 PSMs at `q ≤ 0.01`,
and the nominal 1% is **2.30–2.50% true FDP** (T/E 2.351); at a true FDP of 1% andes keeps
14,142–14,420 PSMs. Ten or eleven accepted PSMs map only to entrapment proteins but match a
yeast sequence once I is folded to L; they count as target. Without that rule the same runs
read 2.49–2.67% and 13,754–13,854. This is the regime where the Java predecessor is
strongest; on the 2026-10-10 count andes trails it by 56 PSMs.

**Phospho-enriched.** One single-shot EasyPhos mouse liver file
(`20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw`), not in any model's training data,
against the 1:1 shuffled mouse entrapment database, 80 Da index slices, 2026-10-07 at
`40774aca`. andes finds 9.6% more PSMs than engine K at a lower measured error. Against the
depositors' engine B identifications (PEP ≤ 0.01, 23,563 scans) andes covers 82.4%, with 95.3%
same-scan backbone agreement. Phosphosite localisation is **not** scored.

**Glyco.** Mouse liver glycoproteome data (PXD005553), against the 1:1 shuffled mouse entrapment database with
`--glyco --glyco-species mouse --decoy-strategy sequon-reverse`. Two published reference
identification sets for the same spectra ship in [`glyco/truth/`](glyco/truth/): engine L (the
depositors') and engine A (a later reanalysis, PXD031032). "Confirmed" is the share of a
reference's spectra that andes identifies with the same peptide.
- *Quick tier:* the `MouseLiver-Z-T-1` fraction as converted mzML, measured 2026-10-07; the
  range is over three Percolator seeds; 19.9 GB peak memory. Reading the `.raw` natively, as
  the commands below do, can give somewhat different counts.
- *Deep tier:* all five fractions, pooled before Percolator, five seeds. Measured 2026-09-05 at
  `14818d3e` on an older, smaller glycan list and converted mzML; not re-measured since.
- The open problem is selection, not generation: almost every reference spectrum produces a
  candidate row, but in the deep tier about a fifth of engine L's spectra lose the per-scan
  choice.

**Label-free quantification with match between runs (ProteoBench).**
- *Data and measurement:* the six `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}`
  runs, where the human, yeast and E. coli ratios between conditions A and B are known. Measured
  2026-10-17 at `47d2844` on a cluster node, scored with the official `proteobench` package
  0.18.6.
- *Which quantity is scored:* the `quantity_normalized` column of the MBR table, which removes each
  run's loading level and its drift along the gradient. The raw `quantity` column, which is what
  andes hands to quantms, scores median |ε| 0.225 and CV 0.240 on the same precursors.
- *Transfer threshold:* at `--mbr-fdr 0.05` (the default is 0.01), andes quantifies 67,035 precursors in
  at least 3 runs, at median |ε| 0.164, E. coli |ε| 0.317 and CV 0.124.
- *Where andes stands:* its median |ε| and CV are lower than those of every match-between-runs
  submission on the leaderboard (snapshot 2026-07-16), and that also holds on the precursors
  andes shares with engine A's and engine C's largest submissions. Its precursor count is not
  the highest: 12 of the 30 submissions report more at transfer q 0.01. Engine labels follow
  #155.
- *Error control* is measured on the same runs with a shuffled entrapment database. Wrong-peak
  transfers are 1.02% of accepted transfers at transfer q 0.01 and 1.85% at 0.05.

**TMT quantification (yeast triple knockout).**
- *Data and measurement:* the hrMS2 and SPS-MS3 runs of PXD017803 (TMTpro 16-plex). Three
  knocked-out proteins must read near zero in their own channels, so the median ratio of
  knockout to non-knockout channel intensity measures interference (lower is better). Measured
  2026-10-10 on the reporter-ion code merged as `0a458a2`; not re-measured since.
- *Extraction:* reporter extraction and impurity correction match a reference isobaric
  implementation exactly on the same scans.
- *Results:* after correction the ratio is 0.210 on hrMS2 and 0.127 on SPS-MS3. Keeping PSMs
  whose precursor isolation purity is ≥ 0.9 (3,487 PSMs) brings SPS-MS3 to 0.040; a reference
  purity filter at the same threshold gives 0.139. Few knockout PSMs stand behind this
  comparison, so read it as a direction.

## 3. How to reproduce

You need `andes` on `PATH`, Docker (for the pinned Percolator image), `curl` and `python3`. Build
andes with `--features thermo` and install the .NET 8 runtime to read `.raw` directly; this is
preferred over converting (see [`reproduce/README.md`](reproduce/README.md)).

### Standard sets: Astral, TMT, UPS1 (tested scripts)

```bash
DATA=~/andes-bench                            # ~10 GB for these three
cd docs/benchmarks/reproduce
./build_databases.sh "$DATA"                  # UniProt download, ~1 min; prints the UPS1 T/E factor
./fetch_spectra.sh   "$DATA" astral tmt ups1  # 4.8 GB from PRIDE
./run.sh             "$DATA" astral tmt ups1  # search + Percolator; prints wall, PSMs, entrapment hits
```

`run.sh` uses 8 threads (`THREADS=` to change), Percolator 3.7.1 `--seed 42 -Y`, and the
`configs/` mods files (Astral: `astral_mods.txt`, `--precursor-tol 10ppm`; TMT: `mods-tmt.txt`).
UniProt is versioned and the original HYE FASTA is no longer served, so counts will differ
slightly from the table; the scripts print a sha256 and sequence count for each input.

### Phospho (commands)

```bash
GLYCO=1 ./build_databases.sh "$DATA"          # also writes databases/mouse_entrap.fasta and prints its T/E
./fetch_spectra.sh "$DATA" phospho            # 3.7 GB
andes --spectrum "$DATA"/phospho/20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw \
      --database "$DATA"/databases/mouse_entrap.fasta \
      --mods ../configs/mods-phospho.txt --threads 8 \
      --fragment-index-slice-da 80 --output-pin phospho.pin
docker run --rm -v "$PWD":/r quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 percolator \
      --seed 42 -Y --results-psms /r/phospho.t.psms --decoy-results-psms /r/phospho.d.psms /r/phospho.pin
```

Count rows at `q-value ≤ 0.01` and their entrapment hits with `python3 ../entrap_common.py count
phospho.t.psms "$DATA"/databases/mouse_entrap.fasta 0.01` (FDP = hits / PSMs × (1 + T/E), with the
T/E that `build_databases.sh` printed), and peptides carrying `79.96`. Expect about 11 minutes on 8 threads.
The slice width is pinned because it sets peak memory (13.7 GB at 80 Da).

### Glyco (commands; scripts in [`glyco/`](glyco/))

```bash
GLYCO=1 ./build_databases.sh "$DATA"
./fetch_spectra.sh "$DATA" glyco-mouse        # 5 fractions, 12.6 GB; the quick tier needs only T-1
cd ..                                         # docs/benchmarks
DB="$DATA"/databases/mouse_entrap.fasta
# --decoy-strategy sequon-reverse is required: plain reversal lets decoys pass the sequon gate.
for f in MouseLiver-Z-T-1 MouseLiver-Z-T-2 MouseLiver-Z-T-3 MouseLiver-Z-T-4 MouseLiver-Z-T-5; do
  andes --spectrum "$DATA"/glyco-mouse/$f.raw --database "$DB" \
        --glyco --glyco-species mouse --decoy-strategy sequon-reverse \
        --threads 8 --output-pin $f.pin                      # writes $f.glyco.pin
done
python3 glyco/pool_pins.py MouseLiver-Z-T-*.glyco.pin > pooled.pin   # pool BEFORE Percolator
docker run --rm -v "$PWD":/r quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 percolator \
      --seed 42 -Y --results-psms /r/pooled.t.psms --decoy-results-psms /r/pooled.d.psms /r/pooled.pin
python3 glyco/eval_yield.py  pooled.pin pooled.t.psms
python3 glyco/eval_entrap.py pooled.pin pooled.t.psms 0.01 "$DB"
for ref in glyco/truth/*_mouse_liver.tsv.gz; do              # the two reference sets
  python3 glyco/score_vs_truth.py "$ref" pooled.pin pooled.t.psms
done
```

For the quick tier, search `MouseLiver-Z-T-1` alone, run Percolator on its `.glyco.pin` without
pooling, and pass `--run MouseLiver-Z-T-1` to `score_vs_truth.py` (a single-file run writes
SpecIds without a file name). About 6 minutes and 20 GB on 8 threads. Pool the deep tier's
fractions: one fraction alone has too few glyco decoys for a stable 1% threshold.

### Label-free quantification, ProteoBench (commands, not a tested script)

The six `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.raw` files of PXD028735 (about
21 GB) and the HYE FASTA linked from the ProteoBench module documentation. The measured run used
peak-picked mzML; `R` is each run's stem.

```bash
# 1. per run: search (writes the run's MS1, MS2 and PSM tables into q/), then Percolator
andes --spectrum $R.mzML --database ProteoBenchFASTA_MixedSpecies_HYE.fasta \
      --mods configs/astral_mods.txt --precursor-tol 10ppm --max-missed-cleavages 2 \
      --threads 6 --output-pin $R.pin --quant-dir q/
percolator --seed 42 --only-psms --results-psms $R.target.psms \
      --decoy-results-psms $R.decoy.psms $R.pin
# 2. once: the library, from every run's Percolator results (one --psms/--decoy-psms per run)
andes quant library --quant-dir q/ --psms $R=$R.target.psms --decoy-psms $R=$R.decoy.psms ...
# 3. per run
andes quant extract --quant-dir q/ --run $R --lfq-tol 10ppm
# 4. once per transfer threshold
andes quant finalize --quant-dir q/ -o pb --mbr-fdr 0.01
```

5. Convert `pb.mbr.parquet` to ProteoBench's "Custom" format. Keep the `identified`,
   `transfer_ms1` and `transfer_ms2` rows, and write one row per precursor with `Sequence`,
   ProForma `Modified sequence`, `Proteins`, `Charge` and the six
   `LFQ_Orbitrap_DDA_Condition_*` columns filled from `quantity_normalized`:

   ```python
   import pandas as pd
   m = pd.read_parquet("pb.mbr.parquet")
   m = m[m.evidence_type.isin(["identified", "transfer_ms1", "transfer_ms2"])]
   wide = m.pivot(index="precursor_key", columns="run", values="quantity_normalized")
   ```
6. Score with `proteobench` 0.18.6:
   `DDAQuantIonModuleQExactive(...).benchmarking(file, "Custom", ...)`. It reports minimum
   1, 3 and 6 runs in one pass.

`andes --spectrum ... --spectrum ... --rescore --lfq --mbr` runs the same four stages in one
process, with one Percolator over all runs. It was not the configuration measured here. On a
cluster node (6 threads) the six searches take about 1 minute each, Percolator about 7 minutes
per run on one CPU (under 3 minutes with its default three threads on three CPUs), extract
about 20 seconds per run and finalize under a minute. CPU time for the six runs: search 1.1k s,
Percolator 2.4k s, library and extract 0.7k s, finalize 0.1k s. Percolator's
`--subset-max-train 300000` cuts its share to 0.8k s; see the match-between-runs section of
`DOCS.md` for what it changes.

### TMT quantification, knockout (commands, not a tested script)

The PXD017803 files `Fig2_m23142_hrMS2.raw` and `Fig2_m23145_SPSms3.raw`, reviewed yeast
proteins, and a mods file with Carbamidomethyl C and TMTpro (304.20715) on K and the peptide
N-terminus as fixed modifications, plus Oxidation M as variable. The measured run read mzML
converted from these files, with 6 threads:

```bash
andes --spectrum Fig2_m23142_hrMS2.raw  --database yeast.fasta --mods mods_tmtpro.txt \
      --precursor-tol 20ppm --protocol TMT --rescore --fdr 0.01 \
      --tmt tmt16 --tmt-level 2 --tmt-correction <kit impurity sheet> --output-tsv hr.tsv
andes --spectrum Fig2_m23145_SPSms3.raw --database yeast.fasta --mods mods_tmtpro.txt \
      --precursor-tol 20ppm --protocol TMT --rescore --fdr 0.01 \
      --tmt tmt16 --tmt-level 3 --tmt-correction <kit impurity sheet> --output-tsv sps.tsv
```

Channel layout: MET6 (P05694) is knocked out in 126–128N, HIS4 (P00815) in 128C–130N and URA2
(P07259) in 130C–132N; 132C–134N is the parental strain. For each PSM of those proteins, divide
the mean of its knockout channels by the mean of the other channels in `<stem>.tmt.parquet`
(rows with `quant_gate_passed`). Report
the median, using only scans whose reporter scan window starts below 126 Th. Repeat on PSMs
with `purity ≥ 0.9`.

## 4. Method notes

- **True FDP by entrapment.** A target PSM that matches only entrapment sequences is false by
  construction: FDP = (entrapment hits / accepted) × (1 + T/E), where T/E is the ratio of
  searchable target to entrapment peptides, measured for each database and never assumed to be
  1. A PSM whose peptide also occurs in a target protein counts as a target: it is entrapment
  only when **all** its proteins are entrapment and its peptide, with I folded to L (the two
  are isobaric), occurs in no target protein sequence. `entrap_common.py` holds this rule
  for `run.sh`, `build_databases.sh` (the T/E factor) and the glyco evaluators.
- **Percolator.** Version 3.7.1, `--seed 42 -Y`, on plain FASTA with andes's `XXX_` decoys, for
  every engine compared on the same host. Percolator's q-value has a floor of 1 / (number of
  top targets), so small or sparse runs return few or no identifications at 1%; pool glyco
  fractions for this reason.
- **ProteoBench terms.** "Min runs" is the number of the six runs in which a precursor
  (peptidoform and charge) must be quantified to be counted. Epsilon is the difference between
  the measured log2 A/B ratio and the expected one; the headline is its median absolute value
  over all precursors. CV is the median coefficient of variation within a condition, averaged
  over the two conditions. A
  transfer q-value is an internal quant-quality score for the peak a precursor is matched to in
  a run where it was not identified; it is not an identification FDR.
- **Comparisons.** Identification comparisons are runs on the same host with matched parameters
  and the same rescoring. ProteoBench comparisons are submissions on its public leaderboard.
  Wall times compare only within one host and session; counts are portable.

## 5. Known gaps

- **The match-between-runs numbers predate #166 and #167**, which since changed how transfers
  are scored and which donors are used; they have not been re-measured on current `main`.
  Quantification rests on one label-free dataset and one Percolator seed, and its best scores
  use andes's own normalisation (`quantity_normalized`; the raw `quantity` is the quantms
  contract). The TMT figures are from 2026-10-10.
- **Only some identification sets measure error.** TMT has no entrapment component. UPS1's
  entrapment is not 1:1, and its nominal 1% is 2.3–2.5% true FDP.
- **Single files.** Phospho is one file, peptide-level only; `--refine` is measured on Astral
  only.
- **Glyco composition on fucose-rich tissue.** On heart and lung the backbone agrees with
  engine L 94–99% of the time but the composition only 68–78% (97% on liver).
- **Not benchmarked:** iTRAQ, timsTOF `.d`, phosphosite localisation. The Java predecessor was
  not re-run under the current defaults; its counts are historical.

## Layout

    docs/benchmarks/
      README.md      this page
      reproduce/     database, download and run scripts (no hardcoded paths)
      glyco/         glyco harness: pooling, yield, entrapment, reference scoring
      configs/       mods files and per-engine parameter files
      entrap_common.py  shared entrapment rules (ALL-ENTRAP_, I/L rescue, T/E)
      tests/         stdlib tests for those rules: python3 -I tests/test_entrap.py

Bulk spectra and third-party engine binaries stay outside git.
