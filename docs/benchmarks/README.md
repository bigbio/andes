# andes benchmarks

Current results, how to reproduce them, the methodology, the known gaps and a short history.
Every figure names the commit, host, thread count and date that produced it; anything not
re-measured is marked as such.

```bash
DATA=~/andes-bench                                  # ~200 GB free for the standard sets

./reproduce/build_databases.sh "$DATA"              # databases from UniProt, ~1 min
GLYCO=1 ./reproduce/build_databases.sh "$DATA"      # ...plus the mouse glyco database
./reproduce/fetch_spectra.sh   "$DATA"              # astral tmt ups1 from PRIDE (4.8 GB)
./reproduce/run.sh             "$DATA"              # search + Percolator + results table
```

Reference identifications for scoring ship in [`glyco/truth/`](glyco/truth/). Read
[§3](#3-methodology-and-the-traps) before comparing your numbers with these.

---

## 1. Current results

### Everything measured, in one table

**Refreshed 2026-10-05** on the benchmark VM (8-thread Xeon Gold 6238), Percolator 3.7.1
`--seed 42 -Y`, `q ≤ 0.01`.

- **Counts** (standard, opt-in, quick-glyco, phospho) are from `main` `7d1e4565`, one session.
  The standard, `--chimeric` and Comet counts reproduce the 2026-09-04 values (`1b8520f8`)
  exactly. #105 (`5e7e6bf1`, current `main`) leaves the three standard counts unchanged; its
  PINs are byte-identical on Astral and UPS1, and TMT gains 9 of 409k rows.
- **Wall times** for the standard rows are #105's: two runs each, alternated with `main` in a
  follow-up session the same day. Comet was re-run in the first session; `main` measured
  within 5% in both sessions.
- **#112** (`40774aca`, 2026-10-07) makes high-res searches retrieve candidates with the
  fragment-ion index. Astral moves to 46,774 PSMs; TMT, UPS1 and phospho counts are unchanged.
  The standard and phospho rows are timed at #112 (standard sets: two runs each).
- **Glyco deep tier** was not re-run: it is still `14818d3e`, with TRFP 1.4.3 mzML and the
  NeuGc ≤ 1 glycan list.

| benchmark | dataset | what is measured | andes | reference | measured error | wall (andes) |
|---|---|---|---:|---:|---|---:|
| **Standard, high-res** | Astral, PXD070049, HCD LFQ | PSMs @ q≤0.01 | **46,774** | Comet 31,435 (+48.8%), 215 s · Java MS-GF+ 26,542 † | entrapment version of the database: at an equal true FDP of 1.00%, +22.2% over the previous in-RAM retrieval; the nominal 1% is **1.07% true** | **160–161 s** |
| **Standard + TMT labels** | a05058, PXD007683, ion-trap CID | PSMs @ q≤0.01 | **12,281** | Comet 10,504 (+16.9%), 76 s · Java 10,651 | not measurable | **53–55 s** |
| **Standard, low-res LFQ** | UPS1, PXD001819, ion-trap CID | PSMs @ q≤0.01 | 15,838 | Comet 14,734 (+7.5%), 46 s · **Java 15,904** | 166 entrapment hits ⇒ **3.6% true FDP** at nominal 1%; Comet: 154 hits on 14,734, the same rate | **34–35 s** |
| **Chimeric** (`--chimeric`, in-RAM path) | Astral | PSMs @ q≤0.01 | **65,028** (+69%; 55,616 distinct scans) | in-RAM baseline 38,394 | not measurable | 337 s ¶ |
| | TMT a05058 | PSMs @ q≤0.01 | 12,540 (+2.1%) | baseline 12,281 | not measurable | 73 s ¶ |
| | UPS1 | PSMs @ q≤0.01 | **17,112** (+8.0%) | baseline 15,838 | 167 entrapment hits — flat against 166 | 50 s ¶ |
| **PTM discovery** (`--refine`, in-RAM path) | Astral | PSMs @ q≤0.01 | **40,028** (+4.3%) at `main` with #106; 40,283 at `7d1e4565` | in-RAM baseline 38,394 | Pass 2 (2,388 modified PSMs): **~2.8% FDP** under the merged 1% threshold, by matched entrapment anchors (`--refine-entrapment`: 27 hits on 1,935; 3.3% before #106) | 363 s ¶ |
| | TMT, UPS1 | — | skipped | high-res only, by design | | |
| **Glyco, deep tier** | pGlyco2 mouse liver PXD005553, 5 fractions, TRFP 1.4.3, `main` `14818d3e` | glycoPSMs @1% | **31,666 ± 9** | pGlyco2 **78.9% confirmed** · MSFragger **88.0% confirmed**, 95.8% peptidoform agreement | **1.11% ± 0.03 true FDP** (1:1 database) | 23–29 min / fraction, 16 cores |
| **Glyco, quick tier** | one pGlyco2 liver fraction (`MouseLiver-Z-T-1`), native `.raw`, bundled pGlyco mouse glycan database (`--glyco-species mouse`, 1,833 compositions) | glycoPSMs @1% | **7,162** (seed 42; seeds 1–5: 7,094–7,137) | pGlyco2 **90.1% confirmed** · MSFragger **90.7% confirmed**, 97.8% / 95.7% peptidoform agreement | **0.89–1.24% true FDP** over 6 seeds (seed 42: 1.24%, CI 0.90–1.68; 1:1 database) | 4,888 s, 8 threads |
| **Phospho-enriched** (first PTM benchmark) | PXD007653 mouse liver EasyPhos, one file (`control2`), Q Exactive HCD, mzML via TRFP 1.4.3 (sha256 `03cdb585…`) | PSMs @ q≤0.01 (Percolator seeds 42, 1–3) | **37,179–37,280** with the fragment-ion index · **26,806–26,883 phospho-bearing** | **Comet 2025.01, same VM, file and settings: 33,928–34,025 (+9.6% for andes), 23,874–23,942 phospho-bearing, 1.74–1.78% FDP, 1,166 s** · MaxQuant (PEP≤0.01, 23,563 scans): **82.4% covered** (19,420 scans), 95.3% same-scan backbone agreement | **1.11–1.16%** true FDP (1:1 database) | **637 s** at #112, 8 threads (949 s at `7d1e4565`; 2026-09 cluster run: 164 s at 32 threads) |

† Java MS-GF+ v20240326 was not re-run in 2026-09 or 2026-10; its counts are historical
(same protocol, earlier session) and it remains ~10-40x slower than andes.

¶ The opt-in and quick-glyco rows are timed at `7d1e4565`, before the speed changes of #105,
#110 and #111, which also apply to those paths; they were not re-timed after them.

**How to read it.**

- **Only UPS1 measures error among the standard sets.** Its nominal 1% is **~3.6% true FDP**
  (166 entrapment hits on 15,838 PSMs, factor 3.42), and Comet sits at the same rate (154 hits
  on 14,734). Astral was checked separately against a 1:1 shuffled entrapment version of its
  database (see [retrieval strategy](#choosing-the-retrieval-strategy)); TMT has no entrapment
  component, so its counts are rescored `q ≤ 0.01` only.
- **Speed against Comet.** #105 (exact-mass candidate lookup, memoised fragment predictions,
  parallel setup, no edge scoring on low-res models), #110 (precursor-calibration pre-pass off
  by default), #111 (strong-mode features for survivors only) and #112 (index retrieval with a
  dense vote buffer) took andes from 274–297 s, 93–102 s and 51–57 s to 160–161 s, 53–55 s and
  34–35 s on Astral, TMT and UPS1: 0.70–0.76x Comet's wall time.
- **Not benchmarked, therefore not claimed:** iTRAQ, timsTOF `.d`, MSFragger on the standard
  sets, Comet's fragment-index mode, and phospho *site localisation*.

### Time, CPU and memory against Comet — 2026-10-07

Same VM and inputs, 8 threads, `/usr/bin/time -v`. andes is #112 (`40774aca`), two runs on the
standard sets and one on phospho, 2026-10-07; Comet 2025.01 is one run per dataset, 2026-10-06.
Phospho PSMs are Percolator seed 42.

| dataset | engine | PSMs @ q≤0.01 | wall | CPU time | peak memory |
|---|---|---:|---:|---:|---:|
| Astral | **andes** | **46,774** | **160 / 161 s** | **787 / 796 s** | **4.1 GB** |
| | Comet | 31,435 | 217 s | 1,594 s | 8.1 GB |
| TMT a05058 | **andes** | **12,281** | **53 / 55 s** | **295 / 303 s** | 5.9 GB |
| | Comet | 10,504 | 77 s | 568 s | 2.9 GB |
| UPS1 | **andes** | **15,838** | **34 / 35 s** | **199 / 205 s** | **2.5 GB** |
| | Comet | 14,734 | 42 s | 309 s | 2.9 GB |
| Phospho (PXD007653) | **andes** | **37,179** (1.11% FDP) | **637 s** | **3,608 s** | 16.4 GB |
| | Comet | 33,984 (1.77% FDP) | 1,096 s | 8,399 s | 6.0 GB |

- **CPU time:** andes uses 36–57% less than Comet on every dataset.
- **Wall time:** 17–42% less than Comet on every dataset.
- **Memory:** lower than Comet on Astral (4.1 vs 8.1 GB) and UPS1; 2x on TMT (5.9 vs 2.9 GB) and
  2.7x on phospho (16.4 vs 6.0 GB). Phospho rose from 13.7 GB at `7d1e4565`: #112 indexes each
  chunk's whole mass interval, which on a phospho search space holds many more forms.

### Opt-in modes

- **`--chimeric` forces `top_n = 1`** in pass 1 (baseline: 10), so its wall time is not
  comparable to baseline. On UPS1, +1,274 PSMs came with entrapment hits flat (166 → 167).
- **`--refine` is high-res only** and skips both low-res sets by design: at low resolution a
  deamidation (+0.984) cannot be told from a C13 isotope error.
- **`--refine` on Astral.** Three fixes landed before the refresh: #100 pairs each Pass-2 decoy
  with its anchor peptide, #102 drops deamidations better explained by a precursor one isotope
  high, and #103 adds `--refine-entrapment` (one shuffled entrapment anchor per real anchor, so
  Pass-2 false discoveries can be counted; the run-level entrapment database never reaches Pass
  2). At `7d1e4565`: 40,283 PSMs (+4.9%), Pass-2 FDP 3.3% (35 entrapment hits on 2,127). #106
  then stopped offering protein-N-terminal Acetyl on internal peptides: acetylations fell from
  512 to 125, PSMs to **40,028 (+4.3%)**, and Pass-2 FDP to **2.8%** (27 hits on 1,935). That is
  under one threshold over the merged PIN; thresholded on their own, Pass-2 rows measured
  0.8–1.1% entrapment FDP in an earlier run.

### Phospho-enriched (PXD007653)

Krahmer et al., mouse liver EasyPhos, Q Exactive HCD, one file
(`20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw`, 102,820 MS2, TRFP 1.4.3), not in any
model's training ledger. Reference: the depositors' MaxQuant PEP ≤ 0.01 scans (23,563;
`glyco/truth/maxquant_mouse_liver_phospho.tsv.gz`). Search: the 1:1 mouse entrapment database,
`configs/mods-phospho.txt` (sha `bdc523fe…`), defaults, 32 threads. The fragment-ion index
(#76) is the default here because the search goes out-of-core.

| arm | model | PSMs @ q≤0.01 (seeds 1–5) | entrapment @1% | true FDP | phospho-bearing PSMs (seed 1) | wall |
|---|---|---:|---:|---:|---:|---:|
| default, **fragment index** (today's default) | `hcd_qexactive_tryp` | 37,280 · 37,271 · 37,258 | 214–216 | **1.15–1.16%** | **26,883** | **164 s** |
| default, enumeration path | `hcd_qexactive_tryp` | 36,817 · 36,872 · 36,922 · 36,862 · 36,865 | 225–244 | 1.22–1.32% | 26,629 | 5,954 s |
| `--protocol phospho`, enumeration path | `hcd_qexactive_tryp_phosphorylation` | 36,835 · 36,819 · 36,826 · 36,780 · 36,884 | 202–210 | 1.10–1.14% | 26,767 | 6,031 s |

**Against Comet** on the same node (Comet 2025.01, OpenMS container, trypsin ≤1 missed
cleavage, 20 ppm precursor, 0.02 Da fragment bins, the same mods, `decoy_search = 1`): 226 s,
33,888–34,025 PSMs, 1.73–1.78% FDP, 23,888 phospho-bearing, 77.9% of MaxQuant scans covered
(andes with the index: 82.5%). andes finds 9.6% more PSMs and 12.5% more phospho-bearing PSMs at a third less measured error,
27% faster (peak RSS 25 GB). The phospho model is identification-neutral against the general
model. Sites are **not** scored.

**Reproduce.** `reproduce/fetch_spectra.sh phospho`, build the mouse entrapment database as for
the glyco tiers, then per arm:

```text
andes --spectrum <file>.mzML --database mouse_entrap.fasta \
      --mods docs/benchmarks/configs/mods-phospho.txt --threads 32 \
      [--protocol phospho] --output-pin <arm>.pin
percolator --seed <1..5> -Y --only-psms=false --results-psms <arm>_s<seed>.psms <arm>.pin
```

Count `q-value ≤ 0.01` rows, `ENTRAP_` hits (FDP = 2 × hits / PSMs) and `79.96`-bearing
peptides. Expect ~3 min on 32 cores, ~17 min on 8.

### Glyco

One dataset, pGlyco2 mouse liver (PXD005553), in two tiers. Two references for the same
spectra ship in `glyco/truth/`: the depositors' pGlyco2 identifications (17,855) and
MSFragger-Glyco's Philosopher-filtered table from PXD031032 (14,626).

**Quick tier** — `MouseLiver-Z-T-1.raw` (2.70 GB, sha256 `2f0142b7…`) read natively, against
`mouse_entrap.fasta` (34,554 sequences = 17,277 UniProt reviewed mouse + shuffled twins, sha256
`5ee15d8d…`), `--glyco --decoy-strategy sequon-reverse`, 8 threads, Percolator `--seed 42 -Y`.
Re-measured 2026-10-05 at `7d1e4565`. The glycan source changed in between: `--glyco` now
requires a glycan database and this run uses `--glyco-species mouse` (1,833 compositions),
where 2026-09-06 used the former built-in list with the gated NeuGc bound (852 compositions).
**Read the columns as two configurations, not a code A/B.**

| | 2026-09-06 (built-in list, gated NeuGc) | **2026-10-05 (bundled mouse list)** |
|---|---:|---:|
| glycoPSMs @1%, seed 42 | 7,122 | **7,162** (3,008 glycopeptides, 896 compositions) |
| seeds 1–5 | 7,078 – 7,122 | 7,094 – 7,137 |
| true FDP (1:1 database) | 1.13% | 0.89 – 1.24% over 6 seeds |
| pGlyco2 confirmed | 86.7% | **90.1%** |
| MSFragger confirmed | 87.8% | **90.7%** |
| same-scan peptidoform agreement, pGlyco2 / MSFragger | 96.3% / 95.6% | **97.8% / 95.7%** |
| search wall, 8 threads | 8,145 s (WSL2 host) | 4,888 s (benchmark VM) |

Selection losses against pGlyco2 fell from 12.9% (wrong target 5.6% + decoy won 7.3%) to 9.4%
(4.0% + 5.3%). #95 (the collapse prefers the monoisotopic hypothesis) landed in between; with
the glycan list also different, the gain cannot be apportioned.

**Deep tier** — all five fractions (`MouseLiver-Z-T-{1..5}.raw`, 12.6 GB) as TRFP 1.4.3 mzML,
same database, 16 threads per fraction, pooled before Percolator, 5 seeds; measured 2026-09-05
at `main` `14818d3e` on the NeuGc ≤ 1 list (612 compositions) and **not re-measured since**.
Result: 31,666 ± 9 glycoPSMs at 1.11% ± 0.03 true FDP (main table). Generation is not the
bottleneck (0.1% of reference spectra produce no row); selection is: 20.5% of pGlyco2's
spectra (11.5% of MSFragger's) are generated but lose the per-scan collapse, and decoys win
about half of those. Its pGlyco2 peptidoform agreement was computed against mislabelled
reference tables and cannot be re-scored without re-running the tier.

**Precursor mono-correction (`--precursor-mono`, issue #64, now the default).** On T-1, 153 of
3,824 reference scans were recorded on an M+3..M+6 isotopologue (88 at exactly +4). Widening the
window cannot fix this (each offset is mass-degenerate with a composition change), so the
corrector refits the MS1 envelope. T-1, 2026-09-07, one binary, 4 threads:

- arm B (off, `0..2`): 7,109 glycoPSMs, 1.10% FDP, pGlyco2 confirmed 3,361 (86.7%);
- arm E (auto, `0..2`): 7,225, 0.91%, 3,444 (88.8%);
- arm F (auto, `0..1`, shipped): 7,149, 0.98%, 3,466 (89.4%); 59 of the 62 firmware-mispicked
  +4 targets confirmed (55 with `0..2`), peptidoform agreement 96.9% vs 96.8%, and the `+2`
  tier falls from 402 to 21 PSMs.

The shipped thresholds reach 112 more reference scans with one wrong shift (back-off 0: 106,
with 8 wrong). Over five fractions, 464 of the 515 scans truly at +3..+6 (90.1%) are confirmed;
pooled arm F gives 34,410 glycoPSMs at 0.98% FDP (2026-09-08; the glycan list also changed, so
not a controlled comparison with the deep tier). On heart and lung the corrector transfers
(84–91% of the firmware population confirmed, against 17% for the heart baseline).

**Isotope-aware collapse (issue #79).** NeuGc and Hex+Fuc differ by 1.020401 Da, close to a
neutron (1.003355), so both fit the window. Penalising the isotope offset (`--glyco-gp-iso`,
default 1.0) cut the swap on heart from 124 to 5 scans (composition agreement 76.2% → 80.8%);
liver, the control, is unchanged.

**Defaults (2026-09-10).** `--precursor-mono auto` became the default (flat yield; pGlyco2
confirmations 3,375 → 3,471). `--glyco-min-core-y 2` (about 500 fewer glycoPSMs) and
`--glyco-pin-curated` (neutral) did not.

### Choosing the retrieval strategy

andes either enumerates every candidate in each precursor window and scores them all, or takes
a shortlist of the 100 candidates whose singly-charged b/y ions best match the spectrum (at
least 3 matched) from a fragment-ion index. Since #112 the index is the default on high-res
data; low-res data, `--chimeric`, `--refine` and `--glyco` keep enumeration.

**Low-res data must not use the index.** Forcing it (same binary, 8 threads, seed 42, 2026-09):
TMT a05058 12,281 → 3,613 PSMs at 1%, UPS1 15,838 → 10,312. At 0.5 Da the bins do not
discriminate.

**On high-res data the index finds more, at the same true error rate.** A September run flagged
the Astral gain (38,394 → 46,774) as possibly optimistic: decoy wins fell 9.4% against 1.8% for
targets, and Astral had no entrapment component. Measured on 2026-10-06 against a 1:1 shuffled
entrapment version of the Astral database (measured T/E = 0.900; FDP = entrapment hits /
accepted × (1 + T/E)), Percolator `-Y`, seeds 42/1/2:

| retrieval | PSMs @ q≤0.01 | true FDP at q≤0.01 | PSMs at a true FDP ≤ 1.00% |
|---|---:|---:|---:|
| enumeration, strong re-ranking of the top 25 (previous default) | 34,429–34,464 | 1.00–1.01% | 34,411–34,457 |
| enumeration, top 50 (diagnostic) | 36,287–36,313 | 0.99–1.02% | 36,254–36,367 |
| enumeration, top 100 (diagnostic) | 38,018–38,147 | 0.87–0.96% | 38,314–38,448 |
| **fragment-ion index (default)** | **42,374–42,610** | 1.06–1.08% | **41,997–42,139 (+22.2%)** |

The index's q-values are slightly optimistic (1.07% true at a nominal 1%), which costs about
one point; at an equal true FDP the gain is +22.2%. Per scan (seed 42): of the 10,838 scans only
the index identifies, **96%** had a peptide the enumeration path never emitted among its 10 rows.
Enumeration ranks the whole window by the rank score, which on high-res data separates targets
from decoys poorly (in a 484-PSM sample, 206 were decoys), so the true peptide often misses the
25 candidates strong mode re-ranks; widening that pool recovers part of the gain, slowly. The
other 4% are the same peptide, accepted because the index's smaller pool sharpens the
competition features (median DeltaRankScore 0 → 2, rank entropy 1.21 → 0.19). Charge 3 gains
most. A wider index pool (200 or 300 candidates, or 2 matched ions) does not recover the 2,751
scans only enumeration identifies.

### Refuted — do not re-try without new evidence

- **Charge-neighbour search** (`--charge-expand`, default off; 2026-10-03). On UPS1, accepted
  PSMs at reported z ≥ 4 fell from 197–210 to 84–91, because nothing tells Percolator which
  charge the MS1 envelope supports. Restore a precursor-envelope charge feature
  (`PrecursorIsotopeKL` is always 0.0) before retesting.
- **Glyco selector and generation changes:** the matched-ion term `--glyco-gp-m`, the two-stage
  split election, and generation-side expansion (wider glycan box, two-axis Y retention, isobar
  resolution) all lost; the oxonium gate does not explain unemitted spectra (it fires for 33 of
  34).
- **Filtering low-information rows out of the glyco PIN before Percolator.** Keeping scans with
  ≥ 40 peaks gave 270.6 ± 139.2 PSMs against 384.6 ± 19.9 ungated (one seed returned zero);
  Percolator needs the low-scoring rows to place a threshold. Use `--trainFDR 0.05` instead
  ([`DOCS.md` §9](../../DOCS.md#9-glycopeptide-search-experimental--advanced-knobs)).

Pending: `--fragment-index-intensity-tiebreak` (default off) breaks the index's top-100 vote
ties by matched intensity. On Astral forced out-of-core it changes 46% of assigned peptides
and raises distinct peptides ~0.30% at flat PSMs and an unchanged decoy share; it waits for
an FDP measurement on phospho.

---

## 2. How to reproduce

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

andes prints the parameters it actually resolved and writes them to `statistics.log`. **Quote
those, not the ones you intended**, since calibration can tighten a window mid-run.

### Glyco

```bash
# 0. Build the database FIRST. It is 1:1 SHUFFLED-SELF entrapment, not a foreign proteome.
GLYCO=1 ./reproduce/build_databases.sh "$DATA"     # writes databases/mouse_entrap.fasta

# One search per fraction (the quick tier is the same recipe on MouseLiver-Z-T-1 alone).
# 1. --decoy-strategy sequon-reverse is REQUIRED, not optional:
#    plain reversal maps an N-X-S/T sequon to S/T-X-N, so reversed decoys sail through the
#    glyco sequon gate and q-values come out anti-conservative.
for f in MouseLiver-Z-T-1 MouseLiver-Z-T-2 MouseLiver-Z-T-3 MouseLiver-Z-T-4 MouseLiver-Z-T-5; do
  andes --spectrum $f.raw --database "$DATA/databases/mouse_entrap.fasta" \
        --glyco --glyco-species mouse \
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
`mouse_entrap.fasta` is the mouse targets plus a shuffled twin of each, tagged `ENTRAP_` inside
the accession (`>sp|ENTRAP_Q99JY4|ENTRAP_TRABD_MOUSE`), so detect it by substring. Build it
with `glyco/build_shuffled_entrap.py`, not `glyco/build_entrap.py`, which appends a foreign
proteome: an independent reproduction with mouse + E. coli measured ~21% more glycoPSMs,
because the ratio (~3.9:1, factor ~4.9), the sequon density and the search-space size all
change at once.

`score_vs_truth.py` works on any dataset with a reference in `truth/`; `make_truth.py` builds
those from pGlyco2 TSV or StrucGP xlsx. To run everything from scratch, see
[`reproduce/`](reproduce/), which also lists what will not reproduce exactly.

## Layout

    docs/benchmarks/
      README.md      this file - current results, method, gaps
      reproduce/     the maintained, path-independent way to run everything
      glyco/         the glyco harness (pooling, yield, entrapment, gap decomposition)
      configs/       per-engine parameter files

Bulk spectra and third-party engine binaries stay outside git.

---

## 3. Methodology, and the traps

**Setup.** 8-thread Intel Xeon Gold 6238 VM, Linux x86_64, for every engine: andes, Java MS-GF+
[v20240326](https://github.com/MSGFPlus/msgfplus/releases/tag/v2024.03.26) and Comet 2025.01
(via OpenMS), with parameters matched per dataset. All are rescored by
`quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2` (`--seed 42 -Y`) on plain FASTA with
andes `XXX_` decoys. Java MS-GF+ PINs come from `MzIDToTsv` + `build_pins.py`, and its Astral
count reuses a prior run. Protein counts are omitted (they need uniform parsimony grouping).

**One variable per comparison.** Five results here had to be withdrawn, each from comparing
runs that differed in more than one respect: a "30% engine regression" was two converter
versions; `--chimeric` "running faster" was `top_n=1` against 10; "2.07x slower than Comet"
predated a default change; a glyco dataset "finding nothing" was Percolator's q floor; "21%
more PSMs" elsewhere was a different entrapment database. Run every arm in **one session, on
one host, with one binary**, and record commit, converter (or native), database build, thread
count and date beside the number.

**Read `.raw` natively** (`--features thermo`); output equals a correct conversion at no speed
cost. ThermoRawFileParser 2.0.0 dropped 26% of the MS2 scans on one pGlyco2 file, and the
effect is file-dependent. If you must convert, use 1.4.3. Details:
[`reproduce/`](reproduce/README.md).

**One rescorer for every engine.** Percolator auto-detects concatenated vs separate
target-decoy from the PIN's shape, so check the mode line in each `.perc.log`.

**Percolator's q has a floor of `1/T_top`.** Identifications at 1% are a step function of any
threshold you sweep, and a run with too few confident targets returns zero at 1% regardless of
quality. On a since-retired human-plasma set, single files gave 0, 0 and 112 glycoPSMs; three
different regimes pooled gave 143; three sceHCD replicates pooled gave 385. So **pool glyco
fractions before Percolator, at least three files, from one acquisition regime.** A rich
fraction (the liver quick tier) clears the floor alone.

**Replicate over seeds.** The 5-seed glyco design has a floor of about 117 PSMs; smaller
effects are *not demonstrable*, not refuted.

### Entrapment FDP

A target PSM matching only entrapment sequences is false by construction, which makes true
error measurable:

```
FDP = (entrapment hits / total accepted) x (1 + T/E)
```

`T/E` is the ratio of **searchable space** and must be measured for your database, never
assumed to be 1:

| database | T : E | factor |
|---|---|---:|
| `yeast_entrap.fasta` (UPS1) | 734,280 : 303,537 tryptic peptides | **3.42** |
| a genuine 1:1 database | 1 : 1 | 2.00 |

Assuming 1:1 has understated error here twice (UPS1 by 1.7x; a foreign-proteome glyco database
by ~2.5x). **UPS1** (2026-09-04): 166 entrapment hits on 15,838 PSMs ⇒ **3.58% true FDP** at a
nominal 1% (2.61% on the cruder protein-count basis). At ~380 PSMs with 1–2 hits, one hit moves
the estimate by ~2.6 points, so "FDP 0.00" means *too few hits to measure*. The Astral HYE
database has no entrapment component.

---

## 4. Known gaps

- **The fragment-ion index retrieves at a much tighter window than the scorer matches at.**
  `RankScorer::feature_match_tolerance()` (a constant 20 ppm on high-res) drives retrieval;
  scoring uses the model's `mme` (0.5 Da in every bundled model). Unmeasured; A/B any change
  to the retrieval window against enumeration on identifications.
- **The Astral database is no longer served.** `ProteoBenchFASTA_MixedSpecies_HYE.fasta`
  returns 404 everywhere, so `build_databases.sh` reconstructs a Human/Yeast/E. coli
  equivalent (30.9k vs 31.9k sequences); expect the Astral count to move by a few hundred PSMs.
- **Java MS-GF+ has not been re-run** under the current defaults. Comet 2025.01 was re-run on
  2026-10-05; Comet's fragment-index mode and MSFragger have not been benchmarked.
- **Two of three standard databases cannot support an entrapment claim.** Astral has none;
  UPS1's is not 1:1.
- **The fragment-ion index does not reproduce the enumeration path's multiplicity copies** in
  the PIN `Proteins` column.
- **Phospho is one file, peptide-level only;** `--refine` is measured only on Astral.
- **Glycan composition on fucose-rich tissue** (#79): on heart and lung the backbone agrees with
  pGlyco2 94–99% of the time but the composition only 68–78% (97% on liver); the core-fucose
  ion Y1+Fuc is not a ladder rung, which is the open lead.
- **Glyco selection is the open problem**; testing it needs a candidate-pool dump under
  production settings.

---

## 5. History

- **2026-09-04:** first same-session head-to-head with Comet: +7.5–22.1% PSMs at 1.04–1.26x
  Comet's wall time (andes 244 / 97 / 50 s). With `--top-n 5` to match Comet's output depth,
  andes ranged from 1.24x slower (Astral) to 0.83x (UPS1).
- **2026-09:** `--gbdt-max-trees` default set to 100: Astral 400 s → 244 s (1.64x, −8 PSMs),
  TMT 111 s → 97 s (1.14x, +3 PSMs).
- **2026-09-04:** `--refine` on Astral gave 43,929 PSMs, before the #100/#102/#103/#106 fixes.
- **2026-09-06:** gating the NeuGc bound on mouse raised the quick tier from 6,532 to 7,122
  glycoPSMs at flat FDP (1.10% → 1.13%); superseded by the bundled glycan databases.
- **2026-09:** the human-plasma glyco set was retired: its reference was a proprietary Byonic
  export that cannot be rebuilt from public artifacts.
- **2026-09:** TRFP 2.0.0 was found to drop 26% of MS2 on one pGlyco2 file; all figures here use
  native reading or 1.4.3.
- **2026-09-04:** benchmark material from five directories was consolidated here.
