# How to reproduce the benchmarks

You need `andes` on `PATH`, Docker (for the pinned Percolator image), `curl` and `python3`. Build andes
with `--features thermo` and install the .NET 8 runtime to read `.raw` directly (preferred over
converting; see [`reproduce/README.md`](reproduce/README.md)).

## Standard sets: Astral, TMT, UPS1 (tested scripts)

```bash
DATA=~/andes-bench                            # ~10 GB for these three
cd docs/benchmarks/reproduce
./build_databases.sh "$DATA"                  # UniProt download, ~1 min; prints the UPS1 T/E factor
./fetch_spectra.sh   "$DATA" astral tmt ups1  # 4.8 GB from PRIDE
./run.sh             "$DATA" astral tmt ups1  # search + Percolator; prints wall, PSMs, entrapment hits
```

`run.sh` uses 8 threads and Percolator 3.7.1 `--seed 42 -Y` with the modification lists it names; UniProt is
versioned and the original HYE FASTA is gone, so counts differ slightly (the scripts print a sha256
for each input).

## Phospho (commands)

```bash
GLYCO=1 ./build_databases.sh "$DATA"          # also writes databases/mouse_entrap.fasta and prints its T/E
./fetch_spectra.sh "$DATA" phospho            # 3.7 GB
andes --spectrum "$DATA"/phospho/20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw \
      --database "$DATA"/databases/mouse_entrap.fasta \
      --variable-mods "Oxidation (M), Acetyl (Protein N-term), Phospho (STY)" --max-mods 4 \
      --threads 8 \
      --fragment-index-slice-da 80 --output-pin phospho.pin
docker run --rm -v "$PWD":/r quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 percolator \
      --seed 42 -Y --results-psms /r/phospho.t.psms --decoy-results-psms /r/phospho.d.psms /r/phospho.pin
```

Count rows at `q-value ≤ 0.01` with `entrap_common.py` (FDP = hits / PSMs × (1 + T/E), T/E from
`build_databases.sh`); it takes about 6 minutes on 8 threads, and the slice width is pinned because it
sets peak memory.

## Glyco (commands; scripts in [`glyco/`](glyco/))

```bash
GLYCO=1 ./build_databases.sh "$DATA"
./fetch_spectra.sh "$DATA" glyco-mouse        # 5 fractions, 12.6 GB; the quick tier needs only T-1
cd ..                                         # docs/benchmarks
DB="$DATA"/databases/mouse_entrap.fasta
# --decoy-strategy sequon is required: other generated decoys carry fewer sequons than their targets.
for f in MouseLiver-Z-T-1 MouseLiver-Z-T-2 MouseLiver-Z-T-3 MouseLiver-Z-T-4 MouseLiver-Z-T-5; do
  andes --spectrum "$DATA"/glyco-mouse/$f.raw --database "$DB" \
        --glyco --glyco-species mouse --decoy-strategy sequon \
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

For the quick tier search `MouseLiver-Z-T-1` alone, skip pooling and pass `--run MouseLiver-Z-T-1` to
`score_vs_truth.py`; pool all fractions for the deep tier, because one fraction has too few decoys.

## Label-free quantification, ProteoBench (tested script; converter tested)

Inputs: the six runs of a module (QExactive `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.raw`
of PXD028735, about 21 GB; Astral `LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw` of
PXD070049, about 15 GB) and `ProteoBenchFASTA_MixedSpecies_HYE.fasta` (31,889 entries). The measured
runs read peak-picked mzML.

Run [`reproduce/run_proteobench.sh`](reproduce/run_proteobench.sh) from a scratch directory; it
writes `q/`, `pb/`, `pb05/` and `custom_input.tsv` to the current directory, and `--module` is `qe`
or `astral`:

```bash
docs/benchmarks/reproduce/run_proteobench.sh qe "$DATA"/ProteoBenchFASTA_MixedSpecies_HYE.fasta \
      "$DATA"/LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.mzML
```

Per run it searches and rescores (`--rescore --percolator-docker`, which runs the pinned Percolator 3.7.1 image) and writes the run's MS1, MS2 and PSM tables to `q/`; then it builds the library from every run's
Percolator results (one `--psms`/`--decoy-psms` per run), extracts each run and finalizes twice
(the default `--quant-fdr` is 0.01, plus the 0.05 run in `pb05/`), and writes the ProteoBench
"Custom" table with `proteobench/andes_to_proteobench.py --column quantity_normalized`. It needs
python3 with pyarrow and prints the same provenance line as `run.sh`.

Score with `proteobench` 0.18.6 (`DDAQuantIonModuleQExactive` or `DDAQuantIonAstralModule`
`.benchmarking(custom_input.tsv, "Custom", ...)`).

The converter keeps `identified`, `transfer_ms2` and `transfer_ms1` rows; `--column quantity` writes
the raw quantity and `--evidence identified` the identification-only table (tests:
`python3 -I tests/test_proteobench_convert.py`). A one-process variant, `andes --spectrum ... --rescore
--lfq --mbr --quant-dir q/`, was not the configuration measured.

## Truth-design quantification suite (commands; harness tested on synthetic data)

[`quant-eval/`](quant-eval/README.md) needs Python 3 with numpy, pandas, scipy and pyarrow, reads
finished tables and runs no search; write one manifest per set (format in `quant_eval/suite.py`).

```bash
# andes (one search per run, Percolator over the set, then the three quant stages)
andes --spectrum $R.mzML --database db.fasta --output-pin $R.pin --quant-dir q/ --threads 8 \
      --predict rt,ms2 --chimeric --refine --max-missed-cleavages 2
andes quant library  --quant-dir q/ --psms target.psms.tsv --decoy-psms decoy.psms.tsv
andes quant extract  --quant-dir q/ --run $R
andes quant finalize --quant-dir q/ -o set --quant-fdr 0.01 --seed 42      # also --seed 1 and 7
# masked copy of the library for the wrong-peak measurement: the same three stages on a
# quant directory whose library lacks 5% of the identified pairs (mask_pairs.parquet)
# MaxQuant 2.7.0.0: one parameter file per set; the MBR-off run changes only
#   <matchBetweenRuns>False</matchBetweenRuns> (and <matchBetweenRunsFdr>False</matchBetweenRunsFdr>)
# FragPipe: IonQuant with --normalization 0
python3 -m quant_eval.suite --manifest manifest.json --out results/ --seeds 42,1,7
python3 -m quant_eval.board --results results/ > scoreboard.md
```

`X_I` on the board is andes with MBR off and `MQ_noMBR` a second MaxQuant manifest entry; ProteoBench
files come from `python3 -m quant_eval.proteobench` or `proteobench/andes_to_proteobench.py` and
nothing is uploaded.

## TMT match between runs (commands, not a tested script)

PXD014750 hrMS2 (six files, TMT11 impurity sheet, human + yeast + entrapment database); the two runs
differ by one flag, and the truth evaluation reads `run.tmt.parquet` as described in [Methods](methods.md#method-notes).

```bash
andes --spectrum f1.mzML ... --spectrum f6.mzML --database hy_entrap.fasta \
      --fixed-mods "Carbamidomethyl (C), TMT6plex (K), TMT6plex (N-term)" \
      --protocol TMT --rescore --psm-fdr 0.01 --tmt tmt11 --tmt-level 2 --tmt-correction lot.txt \
      --threads 16 --output-tsv run.tsv [--mbr]
```

## TMT knockout (commands, not a tested script)

PXD017803 `Fig2_m23142_hrMS2.raw` and `Fig2_m23145_SPSms3.raw`, reviewed yeast proteins,
Carbamidomethyl C and TMTpro (304.20715) on K and the peptide N-terminus fixed, plus Oxidation M
variable; the measured runs read mzML (both files of each acquisition in one process, 16
threads).

```bash
andes --spectrum Fig2_m23142_hrMS2.raw  --database yeast.fasta \
      --fixed-mods "Carbamidomethyl (C), TMTpro (K), TMTpro (N-term)" --variable-mods "Oxidation (M)" \
      --precursor-tol 20ppm --protocol TMT --rescore --psm-fdr 0.01 \
      --tmt tmt16 --tmt-level 2 --tmt-correction <kit impurity sheet> --output-tsv hr.tsv
andes --spectrum Fig2_m23145_SPSms3.raw --database yeast.fasta \
      --fixed-mods "Carbamidomethyl (C), TMTpro (K), TMTpro (N-term)" --variable-mods "Oxidation (M)" \
      --precursor-tol 20ppm --protocol TMT --rescore --psm-fdr 0.01 \
      --tmt tmt16 --tmt-level 3 --tmt-correction <kit impurity sheet> --output-tsv sps.tsv
```

Channel layout and scoring are in [Methods](methods.md#method-notes).
