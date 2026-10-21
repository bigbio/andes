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
| Label-free, match between runs (ProteoBench "Quant LFQ DDA ion QExactive") | PXD028735, 6 runs, Q Exactive HF-X, human/yeast/E. coli mix | **65,391** precursors in ≥ 3 runs; median \|ε\| **0.164**, CV **0.120** (defaults, transfer q 0.01, 2026-10-21) | median \|ε\| 0.190 (engine A); CV 0.120 (engine C, no MBR); most precursors 81,933 (engine E, median \|ε\| 0.388, transferred values not error-controlled) | PSM FDP 0.97%; entrapment-only quantified precursors 0.26% | 6.3k CPU-s for 6 runs |
| Label-free, match between runs (ProteoBench "Quant LFQ DDA ion Astral") | PXD070049, 6 runs, Astral, human/yeast/E. coli mix | **33,876** precursors in ≥ 3 runs; median \|ε\| **0.186**, CV **0.096** (defaults, transfer q 0.01, 2026-10-21) | median \|ε\| 0.123 (engine P, 16,363 precursors); CV 0.070 (engine E); most precursors 93,237 (engine P, median \|ε\| 0.277) | not measured | 5.4k CPU-s for 6 runs |
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
- *Data and measurement:* the two ProteoBench "Quant LFQ DDA ion" modules, six runs each, where
  the human, yeast and E. coli ratios between conditions A and B are known: QExactive
  (`LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}`, PXD028735, Q Exactive HF-X) and
  Astral (`LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}`, PXD070049). Measured
  2026-10-21 at `a27a72c` on a cluster node (8 threads per run), andes defaults throughout
  (§3), Percolator 3.7.1, one seed. Scored with the official `proteobench` package 0.18.6
  ("Custom" format, converted by [`proteobench/andes_to_proteobench.py`](proteobench/)).
- *Which quantity is scored:* `quantity_normalized` of the MBR table, which removes each run's
  loading level and its drift along the gradient, unless the table says raw. The raw `quantity`
  is what andes hands to quantms.
- *Results* at a minimum of 3 runs; precursors / median |ε| / E. coli |ε| / CV:

  | module, `--predict` | transfer q ≤ 0.01, normalised | raw quantity | transfer q ≤ 0.05, normalised | identified only, normalised |
  |---|---|---|---|---|
  | QExactive, `none` | 62,651 / 0.163 / 0.274 / 0.119 | 0.226 / 0.319 / 0.246 | 68,315 / 0.175 / 0.436 / 0.133 | 50,500 / 0.154 / 0.197 / 0.105 |
  | QExactive, `rt` (default) | **65,391 / 0.164 / 0.278 / 0.120** | 0.228 / 0.330 / 0.246 | 71,649 / 0.177 / 0.439 / 0.134 | 52,110 / 0.154 / 0.198 / 0.107 |
  | QExactive, `rt,ms2` | 67,906 / 0.165 / 0.285 / 0.121 | 0.229 / 0.331 / 0.247 | 74,464 / 0.179 / 0.450 / 0.136 | 53,182 / 0.155 / 0.199 / 0.107 |
  | QExactive, `ms2,peptdeep-rt,deeplc-rt` | **68,077 / 0.165 / 0.276 / 0.121** | 0.229 / 0.322 / 0.247 | 74,776 / 0.179 / 0.445 / 0.136 | 53,289 / 0.155 / 0.203 / 0.107 |
  | Astral, `none` | 31,622 / 0.183 / 0.294 / 0.093 | 0.212 / 0.272 / 0.094 | 35,717 / 0.206 / 0.426 / 0.110 | 25,707 / 0.175 / 0.291 / 0.079 |
  | Astral, `rt` (default) | **33,876 / 0.186 / 0.296 / 0.096** | 0.215 / 0.278 / 0.097 | 37,991 / 0.210 / 0.424 / 0.112 | 26,604 / 0.176 / 0.308 / 0.079 |
  | Astral, `rt,ms2` | **34,413 / 0.185 / 0.314 / 0.095** | 0.214 / 0.289 / 0.097 | 38,985 / 0.210 / 0.427 / 0.113 | 26,822 / 0.177 / 0.311 / 0.079 |
  | Astral, `ms2,peptdeep-rt,deeplc-rt` | 33,987 / 0.184 / 0.316 / 0.094 | 0.213 / 0.282 / 0.096 | 38,892 / 0.210 / 0.442 / 0.113 | 26,853 / 0.176 / 0.310 / 0.079 |

  With the default, at least 1 / 6 runs: QExactive 79,707 / 40,628, Astral 44,667 / 21,185.
  `ms2`, `peptdeep-rt` and `deeplc-rt` load the AlphaPeptDeep and DeepLC models
  (`--predictor-models`, ONNX Runtime backend); `rt` is andes's built-in model.
- *By abundance* (default, transfer q 0.01, normalised; quartiles of each module's own mean
  intensity, |ε| / E. coli |ε| / CV): QExactive Q1 0.292 / 0.511 / 0.201, Q2 0.176 / 0.278 / 0.131,
  Q3 0.140 / 0.238 / 0.106, Q4 0.117 / 0.155 / 0.079; Astral Q1 0.320 / 0.411 / 0.174, Q2 0.210 /
  0.282 / 0.119, Q3 0.165 / 0.273 / 0.086, Q4 0.116 / 0.257 / 0.054.
- *Cost* (default, six runs, CPU seconds): QExactive search with prediction and Percolator 5.2k,
  library under 10, extract 0.7k, finalize 0.4k per transfer threshold; Astral 5.1k, under 10,
  0.1k and 0.2k. Each search takes about 3 minutes of wall time, each extract under 20 seconds,
  finalize about a minute. Prediction alone costs 0.4k CPU seconds for six runs with `rt`, 4.2–4.3k
  with `rt,ms2` and 8.3–8.9k with `ms2,peptdeep-rt,deeplc-rt`.
- *Error control* is measured on the QExactive runs with a shuffled entrapment database
  (2026-10-17, before this measurement). Wrong-peak transfers are 1.02% of accepted transfers at
  transfer q 0.01 and 1.85% at 0.05.

**Comparison with the ProteoBench leaderboard.** Public results of both modules
(`Results_quant_ion_DDA` and `Results_quant_ion_DDA_Astral`, both last updated 2026-07-16; read
2026-10-21; scored by `proteobench` 0.17.0). For each engine the submission with the most
precursors and, when different, the one with the lowest median |ε|, at a minimum of 3 runs.
The last column says what limits the error of transferred (match-between-runs) values, taken
from the submission's parameters, comment and the engine's own output; "not checked" means the
submission's parameter file could not be retrieved (the ProteoBench data server was unreachable
on 2026-10-21). Engine labels follow #155; engine P is new.

*QExactive:*

| engine (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | error control on transferred values |
|---|---|---:|---:|---:|---:|---|
| **andes, default** | yes | **65,391** | **0.164** | **0.278** | **0.120** | PSM q ≤ 0.01; transfer q ≤ 0.01 (decoy-estimated wrong peaks and wrong donors) |
| andes, `--predict ms2,peptdeep-rt,deeplc-rt` | yes | 68,077 | 0.165 | 0.276 | 0.121 | as above |
| andes, default, transfer q ≤ 0.05 | yes | 71,649 | 0.177 | 0.439 | 0.134 | transfer q ≤ 0.05 |
| andes, default, raw quantity | yes | 65,391 | 0.228 | 0.330 | 0.246 | as the default |
| engine A, most precursors | yes | 63,701 | 0.200 | 0.306 | 0.156 | identifications 1%; transfers at MBR ion FDR 0.01 (submitted workflow) |
| engine A, lowest \|ε\| | yes | 54,374 | 0.190 | 0.257 | 0.154 | as above |
| engine B, most precursors | yes | 44,299 | 0.210 | 0.386 | 0.219 | PSM and protein 1%; transfers not checked |
| engine B, lowest \|ε\| | no | 29,737 | 0.192 | 0.263 | 0.199 | no transfers |
| engine C, most precursors | yes | 74,697 | 0.222 | 0.460 | 0.132 | identifications about 1%; no transfer error setting in the submitted parameters |
| engine C, lowest \|ε\| | no | 49,433 | 0.198 | 0.263 | 0.120 | no transfers |
| engine D, most precursors | yes | 37,760 | 0.219 | 0.281 | 0.165 | PSM and protein 1%; transfers not checked |
| engine D, lowest \|ε\| | no | 36,176 | 0.214 | 0.255 | 0.165 | no transfers |
| engine E | yes | 81,933 | 0.388 | 0.877 | 0.393 | peptide q ≤ 0.01 on identifications only; every charge state 2–4 of those peptides is integrated in every run and its peak q-value is not applied |
| engine F | yes | 60,206 | 0.216 | 0.514 | 0.227 | PSM 1%; transfers not checked |
| engine G | yes | 55,639 | 0.210 | 0.508 | 0.226 | PSM 1%; transfers not checked |
| engine H, most precursors | yes | 55,943 | 0.211 | 0.692 | 0.235 | identifications 1%; transfers not checked |
| engine H, lowest \|ε\| | yes | 41,016 | 0.195 | 0.394 | 0.213 | as above |
| engine I | no | 72,741 | 0.289 | 0.741 | 0.216 | no transfers |
| engine J | no | 33,040 | 0.209 | 0.387 | 0.224 | no transfers |
| engine J, our run (2026-10-17, its default) | no | 34,677 | 0.210 | 0.401 | 0.225 | no transfers |
| engine J, our run, its feature-linking MBR | yes | 39,257 | 0.221 | 0.563 | 0.234 | feature-quality score; no transfer error rate |

On `quantity_normalized` andes's median |ε| is lower than that of every submission, with or
without match between runs, and its CV equals the lowest (engine C without match between runs,
0.120). Its raw quantity (|ε| 0.228, CV 0.246) is less accurate than most submissions. Six
submissions report more precursors (engine E, engine C's four largest and engine I), at a median
|ε| of 0.214–0.388. The two engine J runs marked "our run" are ours (2026-10-17, its
development version, raw intensities).

*On the precursors andes shares* (default, transfer q 0.01; the submission's values | andes raw
| andes normalised, median |ε| / E. coli |ε| / CV):

| submission | shared | submission | andes raw | andes normalised | only in the submission |
|---|---:|---|---|---|---|
| engine A, most precursors | 52,066 | 0.190 / 0.273 / 0.156 | 0.211 / 0.283 / 0.241 | 0.147 / 0.225 / 0.110 | 11,635 at 0.264 / 0.801 |
| engine A, lowest \|ε\| | 48,958 | 0.188 / 0.239 / 0.155 | 0.206 / 0.272 / 0.241 | 0.143 / 0.207 / 0.109 | 5,416 at 0.213 / 0.517 |
| engine C, most precursors | 51,309 | 0.208 / 0.348 / 0.125 | 0.213 / 0.304 / 0.244 | 0.153 / 0.255 / 0.114 | 23,388 at 0.265 / 0.950 |
| engine C, 57,787 precursors | 50,217 | 0.203 / 0.349 / 0.123 | 0.214 / 0.305 / 0.243 | 0.152 / 0.256 / 0.113 | 7,570 at 0.269 / 0.850 |
| engine E, our reproduction | 55,486 | 0.311 / 0.627 / 0.322 | 0.220 / 0.312 / 0.246 | 0.158 / 0.257 / 0.117 | 24,871 at 0.688 / 1.370 |

Engine E's per-submission table was not retrievable, so its row uses a reproduction: the same
engine version on the same six mzML files with the submission's stated parameters and per-charge
output gives 80,357 precursors at |ε| 0.388, E. coli 0.878 and CV 0.394 (submission: 81,933,
0.388, 0.877, 0.393). Of its 24,871 precursors that andes does not report, 13,852 are charge
states that engine never identified and 19,183 are above its own peak q-value of 0.01, which it
reports but does not apply.

*Astral:*

| engine (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | error control on transferred values |
|---|---|---:|---:|---:|---:|---|
| **andes, default** | yes | **33,876** | **0.186** | **0.296** | **0.096** | PSM q ≤ 0.01; transfer q ≤ 0.01 (decoy-estimated) |
| andes, `--predict rt,ms2` | yes | 34,413 | 0.185 | 0.314 | 0.095 | as above |
| andes, default, transfer q ≤ 0.05 | yes | 37,991 | 0.210 | 0.424 | 0.112 | transfer q ≤ 0.05 |
| andes, default, raw quantity | yes | 33,876 | 0.215 | 0.278 | 0.097 | as the default |
| engine A, most precursors | yes | 39,057 | 0.212 | 0.250 | 0.127 | identifications 1%; transfers not checked |
| engine A, lowest \|ε\| | yes | 29,105 | 0.181 | 0.405 | 0.078 | as above |
| engine B | yes | 20,283 | 0.160 | 0.158 | 0.092 | PSM and protein 1%; transfers within a condition only (submission comment), not checked |
| engine C, most precursors | yes | 36,469 | 0.197 | 0.417 | 0.095 | identifications about 1%; transfers not checked |
| engine C, lowest \|ε\| | no | 24,803 | 0.178 | 0.245 | 0.077 | no transfers |
| engine E | yes | 30,319 | 0.167 | 0.338 | 0.070 | peptide q ≤ 0.01 on identifications only; peak q-value not applied (as in QExactive) |
| engine P, most precursors | yes | 93,237 | 0.277 | 0.691 | 0.159 | run and global precursor q ≤ 0.01 (submission comment) |
| engine P, lowest \|ε\| | yes | 16,363 | 0.123 | 0.168 | 0.079 | as above, plus MS1-level q ≤ 0.01 |

On Astral andes is in the middle: more precursors than engines B and E at a higher median |ε|,
fewer than engine A's and engine C's largest submissions at a lower one, and a CV above most.
Normalisation helps little here (raw |ε| 0.215, CV 0.097). Shared-precursor comparisons need
the per-submission tables, which could not be retrieved; they are not given.

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
GLYCO=1 ./build_databases.sh "$DATA"          # also writes databases/mouse_entrap.fasta (1:1 shuffled)
./fetch_spectra.sh "$DATA" phospho            # 3.7 GB
andes --spectrum "$DATA"/phospho/20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw \
      --database "$DATA"/databases/mouse_entrap.fasta \
      --mods ../configs/mods-phospho.txt --threads 8 \
      --fragment-index-slice-da 80 --output-pin phospho.pin
docker run --rm -v "$PWD":/r quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 percolator \
      --seed 42 -Y --results-psms /r/phospho.t.psms --decoy-results-psms /r/phospho.d.psms /r/phospho.pin
```

Count rows at `q-value ≤ 0.01`, rows whose proteins are all `ENTRAP_` (FDP = 2 × hits / PSMs on
this 1:1 database), and peptides carrying `79.96`. Expect about 11 minutes on 8 threads.
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

### Label-free quantification, ProteoBench (commands; converter tested)

Inputs: the six runs of a module and the HYE FASTA linked from the ProteoBench module
documentation (`ProteoBenchFASTA_MixedSpecies_HYE.fasta`, 31,889 entries).
- *QExactive:* `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.raw` of PXD028735 (about 21 GB).
- *Astral:* `LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw` of PXD070049 (about 15 GB).

The measured runs read mzML peak-picked from these files. `R` is each run's stem. Every flag not
shown is the andes default (20 ppm precursor tolerance, one missed cleavage, Carbamidomethyl C
fixed, Oxidation M and protein N-terminal Acetyl variable, `--predict rt`).

```bash
# 1. per run: search, Percolator (--rescore) and the run's MS1, MS2 and PSM tables in q/
andes --spectrum $R.mzML --database ProteoBenchFASTA_MixedSpecies_HYE.fasta \
      --rescore --output-pin $R.pin --quant-dir q/ --threads 8
# 2. once: the library, from every run's Percolator results (one --psms/--decoy-psms per run)
andes quant library --quant-dir q/ --psms $R=$R.percolator.target.psms.txt \
      --decoy-psms $R=$R.percolator.decoy.psms.txt ...
# 3. per run
andes quant extract --quant-dir q/ --run $R
# 4. once per transfer threshold (the default --mbr-fdr is 0.01)
andes quant finalize --quant-dir q/ -o pb
andes quant finalize --quant-dir q/ -o pb05 --mbr-fdr 0.05
# 5. ProteoBench "Custom" table (python3 + pyarrow); --module qe or astral
python3 proteobench/andes_to_proteobench.py --mbr pb.mbr.parquet --quant-dir q/ \
      --module qe --column quantity_normalized -o custom_input.tsv
```

6. Score with `proteobench` 0.18.6: `DDAQuantIonModuleQExactive` (or `DDAQuantIonAstralModule`)
   `.benchmarking(custom_input.tsv, "Custom", ...)` reports minimum 1 to 6 runs in one pass.

The converter keeps the `identified`, `transfer_ms2` and `transfer_ms1` rows of the MBR table,
takes each precursor's peptidoform and proteins from the quant directory, and writes Unimod
names in ProForma. `--column quantity` writes the raw quantity instead of `quantity_normalized`,
and `--evidence identified` writes the identification-only table. Its tests:
`python3 -I tests/test_proteobench_convert.py`.

`andes --spectrum ... --spectrum ... --rescore --lfq --mbr --quant-dir q/` runs the same stages
in one process, with one Percolator over all runs; it was not the configuration measured here.

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

- **Label-free quantification rests on two datasets and one Percolator seed**, and its best
  scores use andes's own normalisation (`quantity_normalized`; the raw `quantity` is the quantms
  contract). On Astral andes is mid-table, and its match-between-runs entrapment error is
  measured on QExactive only (2026-10-17, an earlier build). The leaderboard was scored by
  `proteobench` 0.17.0 and andes by 0.18.6; several engines' transfer error control could not
  be checked. The TMT figures are from 2026-10-10.
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
      proteobench/   andes MBR output -> ProteoBench "Custom" table (python3 + pyarrow)
      entrap_common.py  shared entrapment rules (ALL-ENTRAP_, I/L rescue, T/E)
      tests/         tests: python3 -I tests/test_entrap.py (stdlib), tests/test_proteobench_convert.py

Bulk spectra and third-party engine binaries stay outside git.
