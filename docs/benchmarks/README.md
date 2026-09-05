# andes benchmarks

## Start here (contributors)

```bash
DATA=~/andes-bench                                  # ~200 GB free for the standard sets

./reproduce/build_databases.sh "$DATA"              # databases from UniProt, ~1 min
GLYCO=1 ./reproduce/build_databases.sh "$DATA"      # ...plus the mouse glyco database
./reproduce/fetch_spectra.sh   "$DATA"              # astral tmt ups1 from PRIDE (4.8 GB)
./reproduce/run.sh             "$DATA"              # search + Percolator + results table
```

**Four rules that are not optional.** Each exists because breaking it produced a wrong
published number in this project, and each is explained where it applies below.

1. **Read `.raw` natively** (`--features thermo`), do not convert. A converter version
   silently cost 30% of identifications on one file; native reading is the reference and
   costs nothing. If you must convert, use ThermoRawFileParser **1.4.3**.
2. **Measure the entrapment ratio; never assume 1:1.** `FDP = hits/total x (1 + T/E)`,
   and `T/E` is a property of the database you built.
3. **Pool at least three files for glyco, all from ONE acquisition regime.** A single file
   returns zero at 1% FDR because Percolator's floor is `1/T_top`.
4. **Report the metric you computed.** `q <= 0.01` is a claim; an entrapment FDP is a
   measurement. On these datasets they differ by 2-3x.

Reference identifications for scoring ship in [`glyco/truth/`](glyco/truth/) — no
re-download, no proprietary parser.


Below: current results, how to reproduce them, the methodology, and the known gaps.
This folder was consolidated from eight documents in September 2026; the superseded ones
are in git history.

**One variable per comparison — the rule every inconsistency here violated.** Over one
week of benchmarking, five results had to be withdrawn or corrected. Every single one came
from comparing runs that differed in more than one respect, where each difference looked
harmless on its own:

| what looked like | what it actually was |
|---|---|
| a 30% engine regression | two different raw-converter versions |
| `--chimeric` running *faster* than baseline | it forces `top_n=1`, baseline uses 10 |
| andes 2.07x slower than Comet | a stale figure predating a default change |
| a glyco dataset finding nothing | Percolator's q floor on a single file |
| 21% more PSMs in an outside reproduction | a differently-built entrapment database |

So: run every arm you intend to compare in **one session, on one host, with one binary**,
and record binary commit, converter (or native), database build, thread count and date
beside the number. A result measured otherwise is not comparable to these and should not
be published next to them.

**Provenance rule.** Every figure below names the commit, host, thread count and date that
produced it. A number without that is not a result. Where something has *not* been
re-measured, this document says so rather than carrying the old value forward silently.

---

## 1. Current results

### Everything measured, in one table

Every row names what produced it. Standard and opt-in rows: commit `1b8520f8`, benchmark
VM (8-thread Xeon Gold 6238), Percolator 3.7.1 `--seed 42 -Y`, `q ≤ 0.01`, **2026-09-04,
one session**. Glyco rows: commit `6b37bb2c` (see the caveat in §2), 5 Percolator seeds, native `.raw`
(or TRFP 1.4.3, which is byte-identical on these files).

| benchmark | dataset | what is measured | andes | reference | measured error | wall (andes) |
|---|---|---|---:|---:|---|---:|
| **Standard, high-res** | Astral, PXD070049, HCD LFQ | PSMs @ q≤0.01 | **38,394** | Comet 31,435 (+22.1%) · Java MS-GF+ 26,542 † | not measurable (no entrapment component) | 244 s |
| **Standard + TMT labels** | a05058, PXD007683, ion-trap CID | PSMs @ q≤0.01 | **12,281** | Comet 10,504 (+16.9%) · Java 10,651 | not measurable | 97 s |
| **Standard, low-res LFQ** | UPS1, PXD001819, ion-trap CID | PSMs @ q≤0.01 | 15,838 | Comet 14,734 (+7.5%) · **Java 15,904** | 166 entrapment hits ⇒ **3.6% true FDP** at nominal 1% | 50 s |
| **Chimeric** (`--chimeric`) | Astral | PSMs @ q≤0.01 | **65,028** (+69%) | baseline 38,394 | not measurable | 322 s |
| | TMT a05058 | PSMs @ q≤0.01 | 12,540 (+2.1%) | baseline 12,281 | not measurable | 72 s |
| | UPS1 | PSMs @ q≤0.01 | **17,112** (+8.0%) | baseline 15,838 | 167 entrapment hits — flat against 166 | 48 s |
| **PTM discovery** (`--refine`) | Astral | PSMs @ q≤0.01 | **43,929** (+14.4%) | baseline 38,394 | not measurable by entrapment (pass 2 is protein-anchored) | 345 s |
| | TMT, UPS1 | — | skipped | high-res only, by design | | |
| **Glyco, deep tier** | pGlyco2 mouse liver PXD005553, 5 fractions | glycoPSMs @1% | **31,658 ± 34** | pGlyco2 17,855: **78.8% confirmed**, 92.2% of its peptides found | **1.02% ± 0.03 true FDP** (1:1 shuffled database) | ~25 min / fraction, 16 cores |
| | | same-scan agreement | 83.8% peptidoform · 99.1% backbone | (`glyco/agreement.py`) | | |
| **Glyco, quick tier** | one pGlyco2 liver fraction (`MouseLiver-Z-T-1`) | glycoPSMs @1% | 6,554 | — | 1.04% true FDP (1:1 database) | 4,592 s, 8 threads |

† Java MS-GF+ v20240326 was not re-run in the 2026-09 session; its counts are historical
(same protocol, earlier session) and it remains ~10-40x slower than andes.

**Not benchmarked yet, and therefore not claimed:** a phospho-enriched (or any
PTM-enriched) dataset, iTRAQ, timsTOF `.d`, MSFragger on the standard sets, and Comet's
fragment-index mode. The four bundled phosphorylation models have never been scored against
a reference. See §5.

### Standard search, in detail


Percolator image `quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2`; the three datasets
were measured in one session so they are mutually comparable.

| dataset | regime | wall | PSMs @ q≤0.01 |
|---|---|---:|---:|
| Astral (PXD070049) | high-res HCD LFQ | **244 s** | **38,394** |
| TMT a05058 (PXD007683) | low-res ion-trap CID, TMT | **97 s** | **12,281** |
| UPS1 (PXD001819) | low-res CID LFQ | **50 s** | **15,838** |

### Effect of the tree-count default (`--gbdt-max-trees`, now 100)

Same binary, same host, same session; the default is the only variable:

| dataset | all trees (previous) | 100 trees (current) | speedup | ΔPSMs |
|---|---:|---:|---:|---:|
| Astral | 400 s / 38,402 | 244 s / 38,394 | **1.64x** | −8 (−0.02%) |
| TMT | 111 s / 12,278 | 97 s / 12,281 | **1.14x** | +3 (+0.02%) |

Identification-neutral on both. The gap between 1.64x and 1.14x is expected: Astral runs in
`strong` score mode where the ensembles both rank candidates and build features, while
low-res TMT runs in `rank` mode where they only build features. A repeat Astral run the
same day gave 262 s, so treat wall times as ±8%, not single-second figures.

### Against Comet — measured head-to-head, 2026-09-04

Comet 2025.01 rev 1 (`4181df6`) was reinstalled and re-run on the same host, same 8 threads,
same Percolator protocol, the same day as the andes numbers above. So this is a real
head-to-head rather than two figures from different sessions:

| dataset | andes | Comet 2025.01 | PSM gain | andes speed |
|---|---|---|---|---|
| Astral | 38,394 / 244 s | 31,435 / 209 s | **+22.1%** | 1.17x slower |
| TMT a05058 | 12,281 / 97 s | 10,504 / 77 s | **+16.9%** | 1.26x slower |
| UPS1 | 15,838 / 50 s | 14,734 / 48 s | **+7.5%** | 1.04x slower |

**andes finds 7.5–22.1% more PSMs for 1.04–1.26x the wall time.** Note how much smaller
that speed gap is than the previously published "2.07x slower on Astral" — that figure
predates the tree-count default, and the gap on UPS1 is now within run-to-run noise.

Validation that the re-run reproduces the original: Comet's Astral count of 31,435 matches
the stored artifact from the earlier benchmark exactly, and its 209 s is close to the 217 s
published on 2026-08-28.

### Matched output depth — the fair speed comparison

The comparison above is **not** like-for-like on I/O: andes `--top-n` defaults to 10 while
Comet's params set `num_output_lines = 5`, so andes was writing roughly twice the PIN rows.
Re-run with `--top-n 5` (measured 2026-09-04, same session):

| dataset | andes `--top-n 5` | Comet 2025.01 | PSM gain | andes speed |
|---|---|---|---|---|
| Astral | 38,607 / 259 s (4.9 rows/spec) | 31,435 / 209 s (4.9) | **+22.8%** | 1.24x slower |
| TMT a05058 | 12,402 / 80 s (6.0) | 10,504 / 77 s (4.9) | **+18.1%** | 1.04x slower |
| UPS1 | 15,723 / 40 s (5.3) | 14,734 / 48 s (4.9) | **+6.7%** | **0.83x — faster** |

At matched depth andes ranges from 1.24x slower to **1.2x faster** than Comet depending on
regime, while finding 6.7-22.8% more PSMs. Output depth barely moves Astral (PIN write is
~4 s of ~250 s) but is worth ~20% on the two smaller, low-res sets.

Two honest caveats. Astral and TMT reuse the original Comet parameter files verbatim, with
only `database_name` repointed; **UPS1 had no stored Comet parameters**, so they were
derived here from the TMT low-res CID parameters with the TMT-specific fixed modifications
removed — a defensible derivation, but ours rather than the original benchmark's. And these
are single runs on a host with ~8% measured run-to-run variance, so read the times as
approximate.

---

### Opt-in modes

All nine arms in one session, 2026-09-04:

| dataset | baseline | `--chimeric` | `--refine` |
|---|---|---|---|
| Astral | 38,394 / 269 s | **65,028** / 322 s | 43,929 / 345 s |
| TMT | 12,281 / 90 s | 12,540 / 72 s | — *(skipped: high-res only)* |
| UPS1 | 15,838 / 52 s, 166 entrap | **17,112** / 48 s, 167 entrap | — *(skipped)* |

**`--chimeric` gains are real and conservatively stated.** On UPS1, the one dataset here
with a measurable error rate, +1,274 PSMs came with entrapment hits flat (166 -> 167). And
the mode **forces `top_n = 1`** (against a baseline default of 10) to avoid blind
multi-emission inflating FDR — so it wins while retaining a tenth of the pass-1 candidates.
That also makes its *wall time* non-comparable to baseline: it is not "baseline plus a
second pass", it is a shallower pass 1 plus a second pass, which is why it can finish
faster on the smaller sets.

**`--refine` is gated to high-resolution data** and skips on both low-res sets by design,
logging `refine is high-res-only and the data is low-res; skipping refinement`. At low
resolution a deamidation (+0.984) is not separable from a C13 isotope error, so identical
counts there are correct behaviour, not a silent no-op. Its Astral gain is **not**
entrapment-validatable: pass 2 is anchored to already-accepted proteins, so it cannot land
on an entrapment sequence by the mechanism the metric relies on.

## 2. Glyco

Two tiers of one dataset, pGlyco2 mouse liver (PXD005553). Two independent references for
the same spectra ship in `glyco/truth/`: the depositors' pGlyco2 identifications (17,855)
and MSFragger-Glyco's, from the Philosopher-filtered table deposited in PXD031032
(14,626). The MSFragger comparison has not been measured yet; the recipe below scores
against both. Quick is one fraction on a VM; Deep is all five
fractions on a cluster. The earlier human-plasma set was retired: its reference was a
proprietary Byonic `.byrslt` export, which cannot be rebuilt from public artifacts.

### Quick tier — one pGlyco2 liver fraction, VM-local, 76 min

`MouseLiver-Z-T-1.raw` from PXD005553 (2.7 GB) against the 1:1 shuffled mouse database, on
a VM at 8 threads: **4,592 s search**, **6,554 glycoPSMs @1%**, **1.04% true FDP**. It
needs no cluster, which is the point, but 76 minutes is a pre-merge check, not an inner
loop. One fraction clears Percolator's q floor here because a liver fraction carries
thousands of confident targets; see the floor rule below before assuming that of any
other single file.

### Deep tier — pGlyco2 mouse liver, 5 fractions, cluster-scale

PXD005553, ~25 min per fraction on 16 cores, **ThermoRawFileParser 1.4.3**, 5 Percolator seeds,
against a 1:1 shuffled-mouse entrapment database built by `glyco/build_shuffled_entrap.py`
(17,537 targets, UniProt reviewed mouse as served on 2026-09-03). Measured at commit
`6b37bb2c`, which is **not an ancestor of `main`** (a development branch whose glyco output
was later shown byte-identical to `main` on the mouse-brain set, but that equivalence has
not been re-checked on this dataset). Treat the row as provisional until re-measured on a
`main` commit; see §5.

| | measured |
|---|---|
| glycoPSMs @1% | **31,658 ± 34** |
| **true FDP** | **1.02% ± 0.03** (1:1 database, factor exactly 2) |
| confirmed vs pGlyco2 | **14,072 / 17,855 = 78.8%** |
| wrong target won | 1,692 (9.5%) |
| decoy won | 1,972 (11.0%) |
| rejected at FDR | 98 (0.5%) |
| never emitted | **21 (0.1%)** |

**The FDR is right where it claims to be** — 1.02% measured against a nominal 1%, on a
database where the correction factor is exactly 2 and there is nothing to get wrong.
**Generation is not the bottleneck here**: 0.1% of reference spectra produce no row.
**Selection is**, at 20.5%, and decoys win 1,972 of those.

andes accepts 31,602 spectra of which 17,530 are not in the pGlyco2 reference. At 1.02%
measured true FDP those are mostly not false positives, but the two searches used different
databases (UniProt mouse reviewed + shuffled twin here; pGlyco2's own there), so read this
as "accepts substantially more at comparable measured error", not as a clean gain.

### ⚠ Why a single glyco file cannot be benchmarked

Percolator's smallest attainable q is `1/T_top`, so a run with too few confident targets
cannot reach `q = 0.01` **at all** — the answer is 0 regardless of data quality. Measured
(2026-09, on a human-plasma set since retired from this benchmark; the mechanism is general):

| configuration | PIN rows | glycoPSMs @1% |
|---|---:|---:|
| sceHCD, 1 file | 7,257 | **0** |
| EThcD, 1 file | 14,349 | **0** |
| sceHCD-EThcD, 1 file | 12,701 | 112 |
| three different regimes pooled | 34,307 | 143 |
| **three sceHCD replicates pooled** | 22,590 | **385** |

Two rules follow. **Pool at least three files** — one is not a benchmark, it is a zero. And
**pool replicates of ONE acquisition regime, not different regimes**: mixing sceHCD with
EThcD gave 143 against 385 for the same number of files, because Percolator fits a single
model and the feature distributions differ between fragmentation chemistries.



### Refuted — do not re-try without new evidence

Each was measured, not argued: the matched-ion selector term `--glyco-gp-m` (every weight
worse, and the selection buckets do not move); the two-stage split election
(fewer correct identifications at higher error); generation-side expansion in general —
wider glycan box, two-axis Y retention, isobar resolution all moved yield **down**; and the
oxonium gate as an explanation for unemitted spectra (it fires for 33 of the 34).

---

## 3. How to reproduce

```bash
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
perc () { docker run --rm --platform linux/amd64 -v "$PWD":/r $PIMG percolator \
            --seed 42 -Y --only-psms=false \
            --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms /r/$1.pin; }
count () { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="q-value") q=i; next} $q<=0.01{c++} END{print c+0}' "$1"; }
```

| dataset | file | database |
|---|---|---|
| Astral | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` | ProteoBench HYE, 31,889 seqs |
| TMT a05058 | `a05058.raw` | `tmt_db.fasta`, human + yeast reviewed, 26,483 seqs |
| UPS1 | `UPS1_5000amol_R1.raw` | `yeast_entrap.fasta` (yeast + E. coli entrapment) |

```bash
andes --spectrum LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw --database hye.fasta \
      --mods configs/astral_mods.txt --precursor-tol 10ppm --enzyme trypsin \
      --threads 8 --output-pin astral.pin
perc astral && count astral.t.psms

andes --spectrum a05058.raw --database tmt_db.fasta \
      --mods configs/mods-tmt.txt --threads 8 --output-pin tmt.pin

andes --spectrum UPS1_5000amol_R1.raw --database yeast_entrap.fasta \
      --threads 8 --output-pin ups1.pin
```

andes auto-detects activation, analyser resolution and labelling from the file, so
tolerances and the model usually need no flags. It prints the parameters it actually
resolved and writes them to `statistics.log` — **quote those, not the ones you intended**,
since precursor calibration can tighten a window mid-run.

### Glyco

```bash
# 0. Build the database FIRST. It is 1:1 SHUFFLED-SELF entrapment, not a foreign proteome.
GLYCO=1 ./reproduce/build_databases.sh "$DATA"     # writes databases/mouse_entrap.fasta

# One search per fraction (the quick tier is the same recipe on MouseLiver-Z-T-1 alone).
# 1. --decoy-strategy sequon-reverse is REQUIRED, not optional:
#    plain reversal maps an N-X-S/T sequon to S/T-X-N, so reversed decoys sail through the
#    glyco sequon gate and q-values come out anti-conservative.
for f in MouseLiver-Z-T-1 MouseLiver-Z-T-2 MouseLiver-Z-T-3 MouseLiver-Z-T-4 MouseLiver-Z-T-5; do
  andes --spectrum $f.raw --database "$DATA/databases/mouse_entrap.fasta" --glyco \
        --decoy-strategy sequon-reverse \
        --threads 8 --output-pin $f.pin            # writes $f.glyco.pin
done

# 2. Pool BEFORE Percolator (one fraction has 0-2 glyco decoys; see the rules below).
python3 glyco/pool_pins.py *.glyco.pin > pooled.pin
perc pooled

# 3. Evaluate. eval_yield.py takes TWO arguments: the pooled PIN and the psms;
#    score_vs_truth.py attributes every miss to a stage against a committed reference.
python3 glyco/eval_yield.py  pooled.pin pooled.t.psms
python3 glyco/eval_entrap.py pooled.pin pooled.t.psms 0.01 "$DATA/databases/mouse_entrap.fasta"
python3 glyco/score_vs_truth.py glyco/truth/pglyco2_mouse_liver.tsv.gz   pooled.pin pooled.t.psms
python3 glyco/score_vs_truth.py glyco/truth/msfragger_mouse_liver.tsv.gz pooled.pin pooled.t.psms  # 2nd engine, same spectra
python3 glyco/agreement.py      glyco/truth/msfragger_mouse_liver.tsv.gz pooled.t.psms             # peptidoform agreement
#    Quick tier (ONE fraction, not pooled): a single-file run writes SpecIds with no file
#    name, so tell the scorers which reference run it is:
#      python3 glyco/score_vs_truth.py --run MouseLiver-Z-T-1 glyco/truth/pglyco2_mouse_liver.tsv.gz liver1.glyco.pin liver1.psms
```

**The entrapment database is 1:1 shuffled-self, and swapping it changes the answer.**
`mouse_entrap.fasta` is 17,537 mouse targets plus a shuffled twin of each — same length,
same amino-acid composition, scrambled order — tagged `ENTRAP_` *inside* the accession:
`>sp|ENTRAP_Q99JY4|ENTRAP_TRABD_MOUSE`. Build it with `glyco/build_shuffled_entrap.py`, not
`glyco/build_entrap.py`; the latter appends a *foreign* proteome, which is a different
experiment. An independent reproduction that used mouse + E. coli measured **~21% more
glycoPSMs**, because a foreign proteome changes three things at once: the ratio (~3.9:1, so
the FDP factor is ~4.9 rather than 2), the sequon density of the entrapment space, and the
total search-space size. Detect entrapment by SUBSTRING, not prefix — the tag is inside the
accession, not at the start of the header.

`score_vs_truth.py` works on **any** dataset with a committed reference in `truth/`; it
replaced a plasma-only script that hardcoded one dataset's paths and could not run
elsewhere. `make_truth.py` builds those references from pGlyco2 TSV or StrucGP xlsx.

**To reproduce any of this from scratch**, see [`reproduce/`](reproduce/) — three scripts that
pull the data from PRIDE, rebuild the databases from UniProt, and run the whole thing with no
hardcoded paths. It also documents what will *not* reproduce exactly, and why.

## Layout

    docs/benchmarks/
      README.md      this file - current results, method, gaps
      reproduce/     the maintained, path-independent way to run everything
      glyco/         the glyco harness (pooling, yield, entrapment, gap decomposition)
      configs/       per-engine parameter files

Benchmark material used to be spread across five directories in and around the repository;
it was consolidated here on 2026-09-04. Bulk spectra and third-party engine binaries
deliberately stay outside git; the workspace-root README says what remains and how to
re-fetch it.


---

## 4. Methodology, and the traps

**Read `.raw` natively; do not convert.** andes does this with `--features thermo`, at no
speed cost and with output identical to a correct conversion. Conversion is where a silent
30% error entered: on one pGlyco2 file, ThermoRawFileParser 2.0.0 wrote 33,893 MS2 where
1.4.3 wrote 45,905, and native reading — Thermo's own RawFileReader, hence the reference —
confirms 45,905. So 1.4.3 is right and 2.0.0 dropped 26% of the scans. It is
FILE-DEPENDENT: on the plasma file both versions agreed exactly, so no converter can be
assumed safe for your data. All figures here used native reading or 1.4.3.

**One rescorer for every engine.** Percolator 3.7.1, pinned, `--seed 42 -Y`. Comparing one
engine's own score against another's rescored q-value is not a comparison. Percolator also
auto-detects concatenated vs separate target-decoy *from the PIN's shape*, so two engines
can silently end up rescored under different modes — check the mode line in each
`.perc.log`.

**Pool glyco fractions before Percolator.** A single fraction yields on the order of 0–2
glyco decoys, so a per-fraction 1% q-value is estimated from almost nothing and swings
between runs. Differences measured that way are noise.

**Percolator's q has a floor of `1/T_top`.** Identifications at 1% are therefore a *step
function* of any threshold you sweep, and plateaus in such a sweep are artifacts. Below
~101 targets the answer is always zero.

**Replicate over seeds.** Single-seed glyco differences are routinely inside seed noise;
the 5-seed design here has a floor of about 117 PSMs, below which an effect is *not
demonstrable* rather than refuted.

### Entrapment FDP

An entrapment database adds a proteome the sample cannot contain; a target PSM matching
only those sequences is false by construction, which makes true error measurable rather
than assumed.

```
FDP = (entrapment hits / total accepted) x (1 + T/E)
```

`T/E` is the ratio of **searchable space** and must be measured for the database in front
of you — never assumed to be 1:

| database | T : E | factor |
|---|---|---:|
| `yeast_entrap.fasta` (UPS1) | 734,280 : 303,537 tryptic peptides | **3.42** |
| a genuine 1:1 database | 1 : 1 | 2.00 |

Assuming 1:1 has produced wrong published numbers here **twice** — understating UPS1's true
error, roughly fivefold on a foreign-proteome database. Peptide space is the better basis: it is
what the search samples, and entrapment proteomes usually have a different length
distribution from the target.

**Measured on UPS1** (2026-09-04): 166 entrapment hits on 15,838 PSMs ⇒ **3.58% true FDP**
at a nominal 1% (3.42 factor), or 2.61% on the cruder protein-count basis. Either way the
nominal 1% is not 1%.

**Know the estimator's resolution.** At ~380 accepted PSMs with 1–2 entrapment hits, one
hit moves the estimate by ~2.6 points — such a design cannot distinguish 1% from 5%, and a
reported "FDP 0.00" means *too few hits to measure*, not *clean*.

**Not every database here has an entrapment component.** The Astral HYE database has none,
so no entrapment FDP is computable from it and its counts are rescored `q ≤ 0.01` only.

---

## 5. Known gaps

- **The Astral database is ProteoBench's own file and is no longer served.** The spectra
  fetch (an earlier version of this document said they did not; that was a pagination bug
  in our script, fixed 2026-09-05), but `ProteoBenchFASTA_MixedSpecies_HYE.fasta` returns
  404 from every URL it was ever at, so `build_databases.sh` reconstructs a Human/Yeast/
  E. coli equivalent from UniProt (30.9k vs 31.9k sequences). Expect the Astral count to
  move by a few hundred PSMs on the reconstruction.
- **Java MS-GF+ has not been re-run** under the current default; every Java figure is
  historical. Comet 2025.01 *was* re-run head-to-head on 2026-09-04 (above). Neither
  Comet's newer fragment-index mode nor MSFragger has been benchmarked here at all.
- **Two of three databases cannot support an entrapment claim.** Astral has no entrapment
  component; UPS1's is not 1:1. Rebuilding both near 1:1 is the fix.
- **The deep glyco tier was measured off-`main`** (commit `6b37bb2c`) and must be
  re-measured on a `main` commit before it is quoted as current. The quick tier is measured
  on this branch with the scripted database (above), so it no longer has that problem.
- **andes has not yet been scored against the MSFragger liver reference**; the table is
  committed and the recipe scores it, but no run has been done since it was added.
- **No PTM-enriched benchmark exists.** `--refine` is measured only on Astral, and the four
  bundled phosphorylation models have never been scored against a reference dataset. A
  public phospho-enrichment set with a deposited identification list is the next dataset
  to add, following the same rules as the glyco tiers (pool, measure T/E, native `.raw`).
- **Glyco selection is the open problem**, and whether those 37% are recoverable by scoring
  at all is unknown — it needs a candidate-pool dump taken at retention time under
  production settings, which does not exist yet.
