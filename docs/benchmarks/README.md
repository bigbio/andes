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
3. **Check Percolator's q floor before trusting a single glyco file.** The smallest
   attainable q is `1/T_top`; a file with too few confident targets returns zero at 1%
   regardless of quality, and must be pooled with others from ONE acquisition regime. A
   rich fraction (the liver quick tier) clears the floor on its own; the plasma files did not.
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
one session**. Glyco rows: 5 Percolator seeds; deep tier at `main` commit `14818d3e`, TRFP
1.4.3 mzML (byte-identical to native reading on these files), NeuGc ≤ 1 glycan list; quick tier
re-measured 2026-09-06 at the commit introducing the gated NeuGc bound, native `.raw`.

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
| **Glyco, deep tier** | pGlyco2 mouse liver PXD005553, 5 fractions, TRFP 1.4.3, `main` `14818d3e` | glycoPSMs @1% | **31,666 ± 9** | pGlyco2 **78.9% confirmed** · MSFragger **88.0% confirmed**, 95.8% peptidoform agreement | **1.11% ± 0.03 true FDP** (1:1 database) | 23–29 min / fraction, 16 cores |
| **Glyco, quick tier** | one pGlyco2 liver fraction (`MouseLiver-Z-T-1`), native `.raw`, gated NeuGc bound (2026-09-06) | glycoPSMs @1% | **7,122** (6,532 with the previous NeuGc ≤ 1 list, same binary) | pGlyco2 **86.7% confirmed** (was 77.9%) · MSFragger 87.8% confirmed, 96.3% / 95.6% peptidoform agreement | **1.13% true FDP** (CI 0.80–1.55; 1.10% before) | 8,145 s, 8 threads (WSL2 host) |

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
(14,626). The recipe below scores against both; the quick tier is measured against both. Quick is one fraction on a VM; Deep is all five
fractions on a cluster. The earlier human-plasma set was retired: its reference was a
proprietary Byonic `.byrslt` export, which cannot be rebuilt from public artifacts.

### Quick tier — one pGlyco2 liver fraction, VM-local, ~2 h

`MouseLiver-Z-T-1.raw` (PXD005553, 2.70 GB, sha256 `2f0142b7…`) read natively, against
`mouse_entrap.fasta` (34,554 sequences = 17,277 UniProt reviewed mouse + shuffled twins,
sha256 `5ee15d8d…`), `--glyco --decoy-strategy sequon-reverse`, 8 threads, Percolator
`--seed 42 -Y`. **Measured 2026-09-06 at the commit that introduced the gated NeuGc bound**
(`--glyco-max-neugc`, parent `fdf1f689`) on an 8-thread WSL2 host with 47 GB. The gate
raised the default glycan list from 612 to 852 compositions on the run's own evidence (log
line `glycan list: NeuGc <= 4 per composition`). The same binary with `--glyco-max-neugc 1`
reproduces the previous default byte-for-byte — 6,532 glycoPSMs, 34 entrapment hits, the
numbers measured 2026-09-05 at `d085c0fb` on the benchmark VM — so the two columns below are
a controlled A/B: one binary, one database, one fraction, one Percolator container.

| | NeuGc ≤ 1 (previous default) | **gated bound (current default)** |
|---|---:|---:|
| glycan compositions searched | 612 | **852** |
| search wall (WSL2, 8 threads) | 8,717 s ‡ | **8,145 s** |
| glyco rows | 41,929 | 42,108 |
| glycoPSMs @1%, seed 42 | 6,532 (2,843 glycopeptides, 902 compositions) | **7,122** (3,047 glycopeptides, 932 compositions) |
| 5 seeds | 6,475 – 6,565 | **7,078 – 7,122** |
| **true FDP** | 1.10% (95% CI 0.76–1.54%, 34 entrapment hits) | **1.13%** (95% CI 0.80–1.55%, 38 hits, 1:1 database) |
| pGlyco2 confirmed (3,877 spectra) | 77.9% | **86.7%** |
| pGlyco2 spectra carrying NeuGc ≥ 2 confirmed (of 501) | 58 | **393** |
| accepted PSMs whose glycan carries NeuGc ≥ 2 | 0 | 770 |
| MSFragger confirmed (3,040 spectra) | 87.8% | 87.8% |

‡ Run at the memory ceiling of a 24 GB guest (peak RSS 23.3 GB); the gated run peaked at
27.0 GB with headroom. Not a runtime comparison. The benchmark VM measured 6,329 s for the
NeuGc ≤ 1 list.

**Why the bound matters here.** The default glycan list is human-tuned and allowed at most
one NeuGc per composition. Mouse is CMAH-competent, and 501 of the 3,877 pGlyco2 reference
spectra on this fraction (12.9%) carry two or more NeuGc; under the old list not one
accepted identification did. A spectrum whose true composition cannot be enumerated is not
left unidentified — it is scored anyway and the best available wrong answer wins — which is
why the loss surfaced as `wrong target` and `decoy won`, not as `never emitted`. Raising the
bound moved 335 of those 501 spectra to `confirmed` and the entrapment FDP did not move
(1.10% → 1.13%, inside one CI); this is the same check that caught the 4,034-composition
list inflating error 5.4×. Against the MSFragger reference, which carries no NeuGc ≥ 2
calls, the two lists are identical (2,669 confirmed either way). `--glyco-taxon human` or
`--glyco-no-neugc` keep the previous behaviour, and on human samples the gate never fires:
the PIN is byte-identical.

Scored against both deposited references for this fraction (gated default, seed 42):

| reference (fraction 1) | spectra | confirmed | wrong target | decoy won | never emitted | peptide coverage | same-scan backbone | same-scan **peptidoform** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| pGlyco2 (depositors) | 3,877 | **86.7%** | 5.6% | 7.3% | 0.1% | 95.1% | 99.4% | **96.3%** |
| MSFragger-Glyco (PXD031032) | 3,040 | **87.8%** | 5.3% | 6.8% | 0.1% | 94.1% | 99.7% | **95.6%** |

Where both engines identify a scan, andes agrees on the full peptidoform 96.3% of the time
with pGlyco2 and 95.6% with MSFragger; the residual disagreements are the isobaric
compositions (Hex + NeuAc ≡ Fuc + NeuGc exactly; Hex + Fuc vs NeuGc within 1.02 Da). An
earlier version of this table reported 83.3% against pGlyco2. That was an artefact:
`make_truth.py` read pGlyco2's five-column glycan vector in the wrong monosaccharide order,
rotating the Fuc/NeuAc/NeuGc labels, and the committed `pglyco2_mouse_*.tsv.gz` tables
carried the rotation. The reader now checks the order against pGlyco2's own `GlyMass`
column on every row of every file and refuses a file that does not fit; the tables were
regenerated 2026-09-06 (on the previous list the corrected figure is 95.2%). Backbone-based
columns — confirmed, coverage, same-scan backbone — were never affected. andes accepts 7,122
spectra of which roughly 3,800–4,500 are in neither reference at a measured 1.13% FDP.

### Precursor mono-correction A/B (`--precursor-mono`, issue #64)

The firmware records a class of wide, high-mass glycopeptide envelopes on an M+3..M+6
isotopologue instead of the monoisotope. On this fraction 153 of 3,824 mass-resolvable
pGlyco2 reference scans sit above the default `0..2` window (88 of them at exactly +4),
and widening the window is not a fix: every extra offset is mass-degenerate with a glycan
composition change, so arms C and D of #64 piled hundreds of arbitrary compositions onto
the largest allowed offset. `--precursor-mono auto` instead reads the preceding MS1,
fits the observed envelope at the reported charge against a glycopeptide isotope model
under "the recorded precursor is M+k" (k = 0..6), and moves the precursor down by
k−1 isotopes when a k > 0 clearly wins (the sweep's +1 takes the last step, so an
overshoot cannot lose the true mass); the window stays `0..2`.

**Measured 2026-09-07, one session, one binary, one database, one Percolator seed** (this
commit; 4-thread, 16 GB sandbox VM; native `.raw`; Percolator 3.7.1 `--seed 42 -Y`;
the gated NeuGc default, 852 compositions). The 16 GB host cannot hold the full
candidate index (~27 GB, killed at the cgroup limit twice), so **both arms ran with
`--glyco-index-sequon-only`** (3.2 M sequon-bearing candidates, 3.4 GB; peak RSS
11.5 GB). The glyco scorer never reads a non-sequon candidate, and arm B reproduces the
full-index quick-tier row above to within Percolator seed noise (42,108 glyco rows in
both; 7,109 vs 7,122 glycoPSMs; 3,361 vs 3,362 pGlyco2 confirmed; 2,669 vs 2,669
MSFragger confirmed; 96.3% vs 96.3% peptidoform agreement; 37 vs 38 entrapment hits),
so the two arms are comparable with each other and, within seed noise, with the table
above. Measured directly on a 1,500-protein subset (8,000 spectra, same binary): the
sequon-only index writes the same 7,113 rows with the same peptides and labels as the
full index; 16 rows (0.2%) differ in `RawScore`/`CandidateRankEntropy` only, and the
flag-off binary is byte-identical to `main`.

| | B: `--precursor-mono off` | **E: `--precursor-mono auto`** | #64 acceptance |
|---|---:|---:|---|
| search wall (4 threads) | 9,561 s | 9,730 s (+1.8%; envelope fit itself is ~20 s) | — |
| MS2 shifted before search | 0 | **447 of 45,905** (shift 1:106 2:37 3:202 4:74 5:28) | — |
| glycoPSMs @1% | 7,109 | **7,225** (+116) | — |
| **true FDP** (1:1 database) | 1.10% (37 hits; CI 0.78–1.52%) | **0.91%** (31 hits; CI 0.62–1.29%) | inside B's CI ✔ |
| pGlyco2 confirmed (3,877) | 3,361 (86.7%) | **3,444 (88.8%)** | ≥ 3,386 ✔ |
| reference spectra gained / lost vs B | — | **+89 / −6** (all 6 lost are FDR-rejected, right answer emitted) | fewer than arm D's 64 flips ✔ |
| MSFragger confirmed (3,040) | 2,669 (87.8%) | **2,746 (90.3%)** | — |
| same-scan peptidoform agreement vs pGlyco2 / MSFragger | 96.3% / 95.7% | **96.8% / 95.9%** | ≥ 96.3% ✔ |
| the 84 target spectra | 0 confirmed (46 wrong target, 38 decoy) | **69 confirmed**, 5 FDR-rejected, 6 wrong target, 4 decoy | — |
| the 62 targets at +4 | 0 | **55 confirmed** + 4 FDR-rejected (59 emitted right) | ≥ 58 confirmed: **55**, short by 3 |
| accepted PSMs by effective offset (shift + `isotope_error`) | +0 5,332 · +1 1,347 · +2 430 | +0 5,325 · +1 1,337 · +2 402 · **+3 16 · +4 80 · +5 50 · +6 15** | no pile-up at the largest offset ✔ |
| HexNAc3 share at +0 / at the corrected tiers | 5.4% / — | 5.5% / **1.2% at +4, 0% at +5 and +6** (arm D: 79% at +4) | matches the offset-0 population ✔ |
| reference coverage of the corrected tiers | — | +4: 71 of 80 in pGlyco2 (89%); +5: 15 of 50; +6: 4 of 15; **0 entrapment hits** in all of +3..+6 | — |

**What the +4 tier is.** 80 accepted PSMs, 71 of them pGlyco2 reference spectra and 47
MSFragger's, no entrapment hit, 1.2% HexNAc3 — the real firmware-failure population, not
the offset-0 population shifted. The +5 and +6 tiers are small (50 and 15) and clean
(0 entrapment, 0% HexNAc3). There is no pile-up: the largest allowed shift carries 15
PSMs, against 613 and 735 for arms C and D.

**Where it falls short.** The strict "≥ 58 of 62 at +4 confirmed" criterion lands at 55:
four more have the right peptidoform emitted but q slightly above 1% (they are searched at
`isotope_error = 1` because of the back-off, and a corrected scan's mass-error columns are
computed one isotope from the searched mono), and three lose the collapse. Arms C and D
reached 58 by widening the window, at the cost of the degenerate tiers above. The
offline envelope fit alone reaches 85 of the 88 reference scans at +4 (below).

**Arm F: the `0..1` window, measured 2026-09-07** (same commit and database, a second
4-thread 16 GB VM, native `.raw`, Percolator 3.7.1 `--seed 42 -Y`, `--glyco-index-sequon-only`).
Once the corrector is on, the `+2` step of the sweep is the one that carries the
(k, X) ≡ (k−1, X + Hex + Fuc − NeuGc) degeneracy, so this arm runs
`--precursor-mono auto --isotope-error 0..1`. The MS1 pass is bit-identical to arm E
(447 of 45,905 MS2 shifted, 1:106 2:37 3:202 4:74 5:28). Search wall 11,740 s and peak
RSS 13.4 GB on this host; not a runtime comparison with the rows above (different VM,
arm B not re-run here, so the per-scan gained/lost column is not available).

| | B: off, `0..2` | E: auto, `0..2` | **F: auto, `0..1`** | #64 acceptance |
|---|---:|---:|---:|---|
| glycoPSMs @1% (seed 42) | 7,109 | 7,225 | **7,149** (seeds 1–5: 7,117–7,158) | — |
| **true FDP** (sequon-corrected, 1:1 database) | 1.10% (37; CI 0.78–1.52%) | 0.91% (31; CI 0.62–1.29%) | **0.98%** (33; CI 0.67–1.38%) | inside B's CI ✔ |
| pGlyco2 confirmed (3,877) | 3,361 (86.7%) | 3,444 (88.8%) | **3,466 (89.4%)** | ≥ 3,386 ✔ |
| pGlyco2 wrong target / decoy / FDR-rejected / not emitted | — | — | 162 / 226 / **21** / 2 | — |
| MSFragger confirmed (3,040) | 2,669 (87.8%) | 2,746 (90.3%) | **2,749 (90.4%)** | — |
| same-scan peptidoform agreement vs pGlyco2 / MSFragger | 96.3% / 95.7% | 96.8% / 95.9% | **96.9% / 96.0%** | ≥ 96.3% ✔ |
| the 84 target spectra | 0 confirmed | 69 confirmed, 5 FDR-rejected, 6 wrong target, 4 decoy | **71 confirmed**, 0 FDR-rejected, 9 wrong target, 4 decoy | — |
| the 62 targets at +4 | 0 | 55 confirmed + 4 FDR-rejected | **59 confirmed** | ≥ 58 ✔ |
| accepted PSMs by effective offset | +0 5,332 · +1 1,347 · +2 430 | +0 5,325 · +1 1,337 · +2 402 · +3 16 · +4 80 · +5 50 · +6 15 | +0 5,419 · **+1 1,552 · +2 21** · +3 14 · +4 87 · +5 48 · +6 8 | decreasing beyond +1 ✔ |
| HexNAc3 share at +0 / +1 / +2..+6 | 5.4% / — / — | 5.5% / — / 1.2% at +4 | 5.5% / 4.1% / **0%** | matches offset-0 ✔ |
| entrapment hits at +0 / +1 / +2..+6 | — | 0 in +3..+6 | 26 / 7 / **0** | — |
| reference coverage of +4 / +5 / +6 | — | 71 of 80 / 15 of 50 / 4 of 15 | **78 of 87** / 13 of 48 / 3 of 8 | — |

**What the `+2` tier was.** Arm B accepted 430 PSMs at `isotope_error = 2` and arm E 402;
with the window closed to `0..1` only 21 survive (all through a corrector shift), yet the
total goes up against B (+40 PSMs, +105 pGlyco2 confirmed, +80 MSFragger confirmed) and
the `+1` tier grows by about 200. Those are the same scans re-accepted one isotope lower
with the other composition of the degenerate pair, and the peptidoform agreement rising
from 96.3% to 96.9% against pGlyco2 says the lower one is the right one more often. Against
E the arm gives up 76 PSMs at 1% (the collapsed `+2` tier) and gains 22 pGlyco2 and 3
MSFragger confirmations; the four `+4` targets E emitted correctly but lost at q are
confirmed here, because their competitors at `isotope_error = 2` no longer exist.

**The three `+4` targets still missed** (of 62): scans 13348 (+3.920 Da) and 30044
(+3.963 Da) sit 40–80 mDa off an isotope multiple, so they are not clean firmware picks
and no envelope correction can place them; scan 15490 (HexNAc4Hex7Fuc1 on YHHYSSNFSIPK,
+4.011 Da) is corrected but a different sequon peptide wins the collapse. The other ten
misses among the 84 are the +5/+6 and non-integer cases (19835 at +5.852, 26334 at +5.105,
27562 at +2.876, 35708 at +6.147), unchanged from arms B–E.

**Arm F on all five fractions, measured 2026-09-08** (same binary, database and VM as the
arm F row above; each fraction searched alone with `--precursor-mono auto --isotope-error 0..1
--glyco-index-sequon-only`, 4 threads; Percolator `--seed 42 -Y` per fraction for the table,
and once over the pooled five PINs below). The question was whether the T-1 result was a
property of that file; it is not. Every fraction carries the same firmware population
(360–432 MS2 shifted, 84–118 reference scans truly at +3..+6), the corrector confirms
89–94% of it on each, and no fraction shows a pile-up, an entrapment excess or a HexNAc3
excess in the corrected tiers.

| | T-1 | T-2 | T-3 | T-4 | T-5 |
|---|---:|---:|---:|---:|---:|
| search wall (s) | 11740 | 9217 | 11429 | 10019 | 11373 |
| peak RSS (GB) | 13.44 | 13.60 | 13.62 | 13.85 | 13.77 |
| glycoPSMs @1% | 7149 | 6616 | 6960 | 6511 | 6919 |
| entrapment FDP | 0.98% (33; CI 0.67% - 1.38%) | 1.15% (36; CI 0.81% - 1.60%) | 0.98% (32; CI 0.67% - 1.38%) | 0.75% (23; CI 0.48% - 1.12%) | 0.86% (28; CI 0.57% - 1.24%) |
| pGlyco2 confirmed | 3466 of 3877 (89.4%) | 3063 of 3362 (91.1%) | 3270 of 3638 (89.9%) | 3068 of 3364 (91.2%) | 3251 of 3614 (90.0%) |
| pGlyco2 FDR-rejected | 21 | 17 | 11 | 7 | 7 |
| MSFragger confirmed | 2749 of 3040 (90.4%) | 2696 of 2975 (90.6%) | 2588 of 2873 (90.1%) | 2579 of 2818 (91.5%) | 2672 of 2920 (91.5%) |
| peptidoform agreement pGlyco2 / MSFragger | 96.9% / 96.0% | 97.3% / 95.8% | 96.1% / 95.9% | 97.7% / 97.2% | 96.6% / 95.7% |
| MS2 shifted | 432 | 360 | 407 | 376 | 372 |
| accepted by effective offset | +0 5,419 · +1 1,552 · +2 21 · +3 14 · +4 87 · +5 48 · +6 8 | +0 4,989 · +1 1,476 · +2 23 · +3 12 · +4 78 · +5 34 · +6 4 | +0 5,105 · +1 1,670 · +2 17 · +3 17 · +4 87 · +5 56 · +6 8 | +0 4,942 · +1 1,411 · +2 9 · +3 12 · +4 87 · +5 42 · +6 8 | +0 5,104 · +1 1,640 · +2 21 · +3 9 · +4 83 · +5 53 · +6 9 |
| firmware scans (true +3..+6) confirmed | 102 of 118 | 79 of 84 | 98 of 106 | 92 of 103 | 93 of 104 |
| true +2 scans confirmed | 10 of 25 (9 corrected) | 6 of 9 (5 corrected) | 13 of 21 (4 corrected) | 4 of 8 (4 corrected) | 7 of 13 (3 corrected) |
| entrapment hits in +2..+6 | 0 | 2 | 0 | 0 | 2 |
| HexNAc3 share in +2..+6 | 0.0% | 1.3% | 1.6% | 2.5% | 0.6% |
| pGlyco2 coverage of +2..+6 | 108 of 178 | 84 of 151 | 104 of 185 | 96 of 158 | 96 of 175 |

*Search wall is this VM's; peak RSS 13.4–13.9 GB with the sequon-only index. "Firmware
scans" are pGlyco2 reference scans whose recorded precursor sits an integer 3–6 isotopes
above peptide + Cam-C + Ox-M + glycan (the true-offset cross-tab used for the T-1 table).
"pGlyco2 coverage of +2..+6" counts accepted PSMs in those tiers that are reference scans.*

**Five fractions summed, by true offset of the pGlyco2 reference scan** (17,855 scans, 17,562
mass-resolvable; "confirmed" is the backbone, "peptidoform-correct" adds the composition):

| true offset | reference scans | corrector shifted | confirmed | peptidoform-correct |
|---|---:|---:|---:|---:|
| 0 | 16,057 | 0 | 14,840 (92.4%) | 14,574 |
| +1 | 787 | 3 | 689 (87.5%) | 634 |
| **+2** | **76** | **25** | **40 (52.6%)** | **18** |
| +3 | 29 | 23 | 23 | 19 |
| **+4** | **395** | **386** | **368 (93.2%)** | **345** |
| +5 | 79 | 69 | 66 | 58 |
| +6 | 12 | 9 | 7 | 7 |
| **+3..+6** | **515** | **487** | **464 (90.1%)** | **429** |

The +4 tier is confirmed at the same rate as the untouched offset-0 population (93.2% vs
92.4%) and is peptidoform-correct in 94% of the confirmed cases; the 515 firmware scans are
2.9% of the reference and were unreachable in every arm before the corrector.

**The cost of the `0..1` window is the true +2 population**, and it is now measured: 76
reference scans over five fractions, of which the corrector shifts only 25 (the k = 2
hypothesis rarely clears the 0.15 fit-gain threshold), 40 are confirmed on the backbone
and only 18 with the right composition. The other 22 "confirmed" are the degeneracy in
reverse: accepted at `isotope_error = 1` with Hex + Fuc in place of NeuGc. Under `0..2` these
76 scans are reachable directly. So the choice is 76 scans mostly-wrong-composition against
the ~400 PSMs per fraction that the `+2` step accepts with an arbitrary composition; the
agreement figures above say the window wins. A cheaper fix than re-widening is to let the
corrector take the k = 2 hypothesis at a lower gain threshold, which would move those scans
into the corrected tiers; not measured.

**Pooled, deep-tier style** (five `.glyco.pin` pooled with `pool_pins.py`, one Percolator run):

| | deep tier, 2026-09-05 (NeuGc ≤ 1 list, full index, mzML, 16 threads) | **arm F, 2026-09-08** (gated NeuGc list, sequon-only index, native `.raw`, 4 threads) |
|---|---:|---:|
| glycoPSMs @1% | 31,666 ± 9 (5 seeds) | **34,410** (seeds 1–5 and 42: 34,373–34,451) |
| **true FDP** (1:1 database) | 1.11% ± 0.03 | **0.98%** (159 hits; CI 0.83–1.15%) |
| pGlyco2 confirmed (17,855) | 78.9% | **90.1%** (16,086; wrong target 4.0%, decoy 5.2%, FDR-rejected 0.5%, never emitted 0.1%) |
| MSFragger confirmed (14,626) | 88.0% | **90.5%** (13,231) |
| same-scan peptidoform vs pGlyco2 / MSFragger | n/a / 95.8% | **96.9% / 96.1%** |

Not a controlled comparison — the glycan list changed in between (the quick tier measured
the gate alone at +9% PSMs on T-1) — but the deep-tier row is now re-measured at the current
defaults plus the corrector, and every column moved the right way at a lower measured error.

**Offline validation of the corrector** (`--precursor-mono-dump`, 23 s for the file, no
search). For every pGlyco2 reference scan the true offset is the integer k that makes
recorded − (peptide + Cam-C + Ox-M + glycan) an isotope multiple. Best-fitting
hypothesis vs true offset, 3,824 resolvable scans:

| true offset | n | fit picks 0 | picks 1 | picks 2 | picks 3 | picks 4 | picks 5 | picks 6 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| −1 | 30 | 25 | 3 | 2 | 0 | 0 | 0 | 0 |
| 0 | 3,488 | 3,459 | 26 | 1 | 0 | 1 | 0 | 1 |
| 1 | 158 | 9 | 137 | 11 | 0 | 1 | 0 | 0 |
| 2 | 25 | 0 | 1 | 22 | 1 | 1 | 0 | 0 |
| 3 | 10 | 2 | 0 | 0 | 7 | 1 | 0 | 0 |
| 4 | 88 | 1 | 0 | 0 | 1 | 86 | 0 | 0 |
| 5 | 17 | 2 | 0 | 0 | 0 | 1 | 14 | 0 |
| 6 | 8 | 0 | 1 | 0 | 0 | 0 | 1 | 6 |

With the shipped rule (fit ≥ 0.90, gain ≥ 0.15 over the recorded hypothesis, mono SNR
≥ 3, back-off 1) the scans the `0..2` sweep can reach go from 3,671 to 3,783 (+112) with
**one** wrong shift (a true −1 scan). Sweeping the thresholds moved this by ±4; back-off 0
(apply the best fit verbatim) reaches 3,777 with 8 wrong shifts, which is why the default
backs off. The 30 true −1 scans (recorded one isotope *below* the monoisotope) are
unreachable in every arm and are not touched.

**Step 1 of #64, settled.** The Thermo reader's precursor m/z equals both the isolation
window target and the trailer `Monoisotopic M/Z` on all 45,905 MS2 of this file (probe
over `thermorawfilereader` 0.7.0), so the trailer holds nothing the reader is not already
using; the +4 scans are firmware picks.

Reproduce: the quick-tier recipe with `--precursor-mono auto --isotope-error 0..1` (add
`--glyco-index-sequon-only` on a < 32 GB host); score with `score_vs_truth.py --run
MouseLiver-Z-T-1 --buckets ...` and read the `MonoShift` / `isotope_error` PIN columns for
the per-offset table.

This number needs no cluster, which is the point, but two hours is a pre-merge check,
not an inner loop. One fraction clears Percolator's q floor here because a liver fraction
carries thousands of confident targets; see the floor rule below before assuming that of
any other single file.

### Deep tier — pGlyco2 mouse liver, 5 fractions, cluster-scale

PXD005553 `MouseLiver-Z-T-{1..5}.raw` (12.6 GB) as **ThermoRawFileParser 1.4.3** mzML
(byte-identical to native reading on these files), against `mouse_entrap.fasta` (34,554
sequences, sha256 `5ee15d8d…`, 1:1 shuffled, factor exactly 2.0), `--glyco
--decoy-strategy sequon-reverse`, 16 threads per fraction on the SLURM cluster,
pooled before Percolator 3.7.1, 5 seeds. **Measured 2026-09-05 at `main` commit `14818d3e`**
(binary sha256 `6de3c8db…`, rustc 1.85); the earlier off-`main` figure of 31,658 ± 34
reproduces within seed noise. **These numbers pre-date the gated NeuGc bound** (NeuGc ≤ 1
list, 612 compositions). The single-fraction A/B in the quick tier shows +9% glycoPSMs at flat
FDP from the bound alone, so the deep tier is expected to move; it has not been re-measured.

| | measured |
|---|---|
| search wall | 1,350–1,731 s per fraction (16 cores); fraction 1 = 45,905 MS2, 41,929 glyco rows |
| glycoPSMs @1% | **31,666 ± 9** (5 seeds: 31,653 / 31,674 / 31,659 / 31,673 / 31,669) |
| **true FDP** | **1.11% ± 0.03** (1.08–1.15%; 1:1 database) |

Scored against both deposited references (seed 1):

| reference | spectra | confirmed | wrong target | decoy won | FDR-rejected | never emitted | peptide coverage | same-scan backbone | same-scan **peptidoform** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pGlyco2 (depositors) | 17,855 | **78.9%** | 9.5% | 11.0% | 0.5% | 0.1% | 92.3% | 99.1% | n/a § |
| MSFragger-Glyco (PXD031032) | 14,626 | **88.0%** | 5.3% | 6.2% | 0.5% | 0.1% | 96.4% | 99.9% | **95.8%** |

§ The 83.8% previously in this cell was computed against the mislabelled pGlyco2 tables (see
the quick tier) and the deep-tier Percolator outputs were not retained, so it cannot be
re-scored without re-running the tier. Every other column of this row is backbone-based
and unaffected.

**The FDR is right where it claims to be**, and **generation is not the bottleneck**: 0.1%
of reference spectra produce no row. **Selection is**: 20.5% of pGlyco2's spectra (11.5%
of MSFragger's) are generated and scored but lose the per-scan collapse, and decoys win
about half of those. andes accepts 31,653 spectra of which ~17.5k are in neither
reference at 1.11% measured FDP; the two searches used different databases (UniProt
reviewed + shuffled twin here; each engine's own there), so read that as "accepts
substantially more at comparable measured error", not as a clean gain.

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
error by 1.7x, and a foreign-proteome glyco database's (factor ~4.9) by ~2.5x. Peptide space is the better basis: it is
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
- **No PTM-enriched benchmark exists.** `--refine` is measured only on Astral, and the four
  bundled phosphorylation models have never been scored against a reference dataset. A
  public phospho-enrichment set with a deposited identification list is the next dataset
  to add, following the same rules as the glyco tiers (pool, measure T/E, native `.raw`).
- **Glyco selection is the open problem**, and whether those 37% are recoverable by scoring
  at all is unknown — it needs a candidate-pool dump taken at retention time under
  production settings, which does not exist yet.
