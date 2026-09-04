# How to run the benchmarks

Every command here has been executed as written. Where a figure is quoted it names the
commit, host and date that produced it. If you change a dataset, a database or a
Percolator flag, you have changed the benchmark — say so when reporting the number.

## 0. The invariants

- **One rescorer for every engine.** Percolator 3.7.1, pinned by digest, `--seed 42 -Y`.
  Comparing an engine's own score against another engine's rescored q-value is not a
  comparison.
- **Report the metric you actually computed.** `q ≤ 0.01` is a *claim* the rescorer makes;
  an entrapment FDP is a *measurement*. They are not interchangeable, and on these datasets
  they differ by 2-3x (§4).
- **One variable per comparison.** Same binary, same host, same day, or it is not an A/B.
- **Speed is not portable.** Wall times are only comparable within one host and one
  session; identification counts are portable.

```bash
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
perc () {  # perc <stem>  -> <stem>.t.psms / <stem>.d.psms
  docker run --rm --platform linux/amd64 -v "$PWD":/r $PIMG percolator \
    --seed 42 -Y --only-psms=false \
    --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms /r/$1.pin
}
count () { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="q-value") q=i; next} $q<=0.01{c++} END{print c+0}' "$1"; }
```

## 1. The three standard datasets

| dataset | regime | accession | file | database |
|---|---|---|---|---|
| Astral | high-res HCD LFQ | PXD070049 | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` | ProteoBench mixed-species HYE, 31,889 seqs |
| TMT a05058 | low-res ion-trap CID, TMT | PXD007683 | `a05058.raw` | UniProt human + yeast reviewed, 26,483 seqs |
| UPS1 | low-res CID LFQ | PXD001819 | `UPS1_5000amol_R1.mzML` | `yeast_entrap.fasta` (yeast + E. coli entrapment) |

andes reads Thermo `.raw` natively; the Astral figures below were produced from converted
mzML, which is noted because it is a difference.

```bash
# Astral — 262 s / 244 s and 38,394 PSMs on commit 1b8520f8, benchmark VM, 8 threads, 2026-09-04
andes --spectrum astral.mzML --database hye.fasta \
      --mods astral_mods.txt --precursor-tol 10ppm --enzyme trypsin \
      --threads 8 --output-pin astral.pin
perc astral && count astral.t.psms

# UPS1 — 50 s and 15,838 PSMs, same binary and session
andes --spectrum UPS1_5000amol_R1.mzML --database yeast_entrap.fasta \
      --threads 8 --output-pin ups1.pin
perc ups1 && count ups1.t.psms

# TMT — same shape; mods in configs/mods-tmt.txt, and set --protocol TMT
andes --spectrum a05058.raw --database human_yeast.fasta \
      --mods ../configs/mods-tmt.txt --protocol TMT \
      --threads 8 --output-pin tmt.pin
perc tmt && count tmt.t.psms
```

andes auto-detects activation, analyser resolution and labelling from the file, so
tolerances and the scoring model usually need no flags. It prints the parameters it
actually resolved, and writes them to `statistics.log` — **quote those, not the ones you
intended**, since precursor calibration can tighten a window mid-run.

## 2. The glyco datasets

Glyco has its own harness in [`../../benchmarks/glyco/`](../../benchmarks/glyco/) — read
its README first; it encodes two rules that are not optional.

| dataset | regime | accession | truth available |
|---|---|---|---|
| mouse brain, 6 fractions | HCD + EThcD | PXD011533 | Byonic `.byrslt` (depositors') |
| human plasma, R1-R3 | sceHCD | PXD030622 | Byonic `.byrslt` (depositors') |

```bash
# one search per fraction/replicate, THEN pool, THEN one Percolator run
for f in Frac1 Frac2 Frac3 Frac4 Frac5 Frac6; do
  andes --spectrum $f.mzML --database mouse_entrap.fasta --glyco \
        --threads 8 --output-pin $f.pin          # writes $f.glyco.pin
done
python3 ../../benchmarks/glyco/pool_pins.py *.glyco.pin > pooled.pin
perc pooled
python3 ../../benchmarks/glyco/eval_yield.py  pooled.t.psms      # no truth needed
python3 ../../benchmarks/glyco/eval_entrap.py pooled.t.psms ...  # entrapment FDP
python3 ../../benchmarks/glyco/gap_decompose.py                  # vs Byonic, per scan
```

**Pooling is mandatory.** One fraction produces on the order of 0-2 glyco decoys, so a
per-fraction 1% q-value is estimated from almost nothing and swings wildly between runs.
Differences measured per-fraction are noise.

**Percolator's q has a floor of `1/T_top`.** Identifications at 1% are therefore a *step
function* of any threshold you sweep — sweeping a parameter and reading IDs@1% will show
plateaus that are artifacts. Below ~101 targets the answer is always zero.

**Truth must be filtered.** Byonic's `PQMs` table emits one row per (PSM, protein), so
peptides in several proteins are counted several times unless deduped by `(fraction, scan)`.
A raw MSFragger `psm.tsv` is pre-FDR and roughly a third of its rank-1 glyco rows are that
engine's own decoys. Both mistakes have been made here and both inflate the comparator.

## 3. Reproducing a competitor

Configs for the other engines are in [`configs/`](configs/) and the driver scripts in
[`scripts/`](scripts/). Two cautions from experience:

- Engines differ in what they call a PSM. Filter to the same thing before comparing —
  e.g. an open/offset search will report unmodified peptides alongside modified ones.
- Percolator auto-detects concatenated vs separate target-decoy **from the PIN's shape**,
  so two engines can end up rescored under different modes without saying so. Check the
  mode line in each `.perc.log`.

## 4. Entrapment FDP — how to compute it correctly

An entrapment database adds a proteome the sample cannot contain. A target PSM matching
only those sequences is false *by construction*, which makes the true error measurable
instead of assumed.

```
FDP = (entrapment hits / total accepted) x (1 + T/E)
```

`T/E` is the ratio of **searchable space**, and it must be measured for the database in
front of you:

| database | T : E | correct factor |
|---|---|---|
| `yeast_entrap.fasta` (UPS1) | 6,733 : 4,531 proteins, or 734,280 : 303,537 tryptic peptides | 2.49 (proteins) / **3.42 (peptides)** |
| plasma `human_entrap.fasta` | 20,411 : 4,531 proteins | **9.81** |
| a true 1:1 database | 1 : 1 | 2.00 |

**Assuming 1:1 has produced wrong published numbers twice in this project** — understating
UPS1's true FDP and understating plasma's roughly fivefold. Peptide space is the better
basis: it is what the search samples, and entrapment proteomes often have a different
length distribution from the target.

**Know the estimator's resolution.** At ~380 accepted PSMs with 1-2 entrapment hits, one
hit moves the estimate by ~2.6 points — such a design cannot distinguish 1% from 5%, and a
reported "FDP 0.00" means *too few hits to measure*, not *clean*. Check that the hit count
can support the precision you are claiming.

**Not every database here has an entrapment component.** The Astral HYE database has none,
so no entrapment FDP is computable from it — see
[`2026-09-04-tree-count-default-and-entrapment-audit.md`](2026-09-04-tree-count-default-and-entrapment-audit.md).

## 5. Before you publish a number

- Name the commit, host, thread count and date.
- State the metric: rescored `q ≤ 0.01`, or a measured entrapment FDP with its factor.
- For yield claims, run several Percolator seeds and report the spread — single-seed glyco
  differences are routinely inside seed noise.
- A flat or byte-identical result is a red flag: check the arm actually differed before
  reporting it as neutral.
