# Public DIA benchmark (`--dia`) on slices of the PXD070049 Astral runs

The benchmark in [`../README.md`](../README.md) needs a private paired-entrapment library, DIA-NN
reports and about 16 GB of memory. This one uses public inputs only and runs on 4 cores and 15 GB,
so a change to `--dia` can be measured by anyone. It was used to measure the merged DIA commits
against each other and to test six further changes (2026-11-19 to 2026-11-20). The sections after
the results record every variant that was tested, kept or dropped, with its numbers.

The numbers here do not compare with [`../README.md`](../README.md): the library is 9 times
larger, the input is a 70 s slice at the densest part of the gradient, and there is no DIA-NN
reference.

## Data

- PRIDE PXD070049, Orbitrap Astral, 300 isolation windows of 2 Th, 1.65 s cycle.
  - Tuning runs: `LFQ_Astral_DIA_15min_50ng_Human_01` and `_02` (292,883 / 292,914 MS2 scans,
    2,928 MS1 scans, 26.85 min).
  - Held-out runs, used only to confirm a variant chosen on the tuning runs:
    `LFQ_Astral_DIA_15min_50ng_Human_03` and `LFQ_Astral_DIA_5min_50ng_Human_01` (15.4 min).
- mzML: ThermoRawFileParser 1.4.5 (bioconda, mono 6.12), `-f 1 -l 2` (centroided).
- Slices ([`slice_mzml.py`](slice_mzml.py)): RT 950–1020 s of runs 01, 02 and 03 (12,825 /
  12,823 / 12,826 MS2 scans, 128 MS1 scans each, about 42 cycles) and RT 545–615 s of the 5-min
  run (12,844 MS2 scans; the same relative start in the gradient, chosen before any measurement).
- Database ([`build_library.py`](build_library.py), seed 20261119): every fully tryptic peptide of
  the reviewed human UniProt proteome (UP000005640, release 2026_04, 20,414 entries) of 7–35
  residues without a missed cleavage, and for each one shuffled partner that keeps the C-terminal
  residue and is neither a human tryptic peptide nor a substring of a human protein (I = L).
  473,370 targets and 473,370 `_p_target` entrapments; 1,405,100 target and 1,405,100 entrapment
  precursors at charges 1–5 within 370–990 m/z, so r = 1. Sequences that start with M are left
  out, because andes clips an initiator Met from each FASTA entry.
- sha256 of the slices, the FASTA and the library table: [`SHA256SUMS`](SHA256SUMS).

## Protocol ([`bench.sh`](bench.sh))

Per slice: `andes --dia --enzyme nocleavage --model hcd_astral_tryp --charge 1..5 --min-length 7
--max-length 35 --predict rt,ms2 --threads 4 --rescore`, Carbamidomethyl (C) fixed, Oxidation (M)
variable (at most 1), andes's own decoys. `--dia --rescore` runs Percolator 3.7.1 (bioconda build
`h3b5f4bd_2`; `--seed 42 -Y --trainFDR 0.05 --num-threads 4`), the per-run RT refinement, and
Percolator again. Percolator is then re-run on the RT-refined PIN for seeds 1, 2, 3 and 7 with
the same arguments; seed 42 is andes's own second pass.

Each seed is evaluated with [`../eval_dia.py`](../eval_dia.py): the best row per precursor,
precursor-level target-decoy q-values, entrapment FDP_c = E (1 + 1/r) / N with r = 1, and
N@FDP_c 1 % / 2 %, the number of precursors at an entrapment FDP of 1 % / 2 %. Every number below
is the mean of the five seeds with (min–max); differences are paired by seed
([`compare.py`](compare.py)).

Host: 4 cores of an Intel Xeon at 2.10 GHz, 15 GB, Linux 6.18; shared with builds running under
`nice`, so wall times are inflated and CPU time is the comparable figure. Timings are from
`/usr/bin/time -v`.

"Native" below means the `percolator-native` crate (`--rescore-native`) through the small CLI in
[`percnative/`](percnative), used for triage; every kept or final number is from Percolator 3.7.1.

## Results

### The merged DIA commits

N@FDP_c 1 %, tuning slices, Percolator 3.7.1. Per-seed values, FDP_c, time and memory of these
runs are in [`baselines.tsv`](baselines.tsv), those of every variant below in
[`experiments.tsv`](experiments.tsv).

| commit | slice 01 | slice 02 | FDP_c at q ≤ 0.01 | search wall / CPU / peak RSS |
|---|---|---|---|---|
| 2e3b8f8 (#402, window search, MS1 gate, co-elution columns) | 4983 (4957–5010) | 5084 (5061–5109) | 0.0102 / 0.0090 | 70–76 s / 224 s / 1.4 GB |
| d2088cc (#410, 50-row pool, predicted-MS2 re-rank) | 5501 (5456–5562) | 5418 (5375–5462) | 0.0082 / 0.0108 | 142 s / 483 s / 2.4 GB |
| 8c69c7a (#417, summed-spectrum and elution-profile columns) | 5657 (5583–5700) | 5821 (5801–5852) | 0.0108 / 0.0071 | 146 s / 498 s / 2.4 GB |
| 3f2b66c (main after #436) | 5669 (5657–5678) | 5798 (5751–5851) | 0.0110 / 0.0076 | 148–151 s / 510 s / 2.4 GB |

- #402 → #410 is +10.4 % / +6.6 % and #410 → #417 is +2.9 % / +7.4 %, each higher on 5 of 5
  seeds on both slices. On the full runs [`../README.md`](../README.md) reports +13.6 % / +9.2 %
  and +5.6 % / +7.5 %.
- #417 → main is +11 / −23, inside seed noise; nothing between them changes the Astral search.
- Held-out slices, main: run 03 5770 (5714–5809), FDP_c at q ≤ 0.01 0.0079; 5-min run
  13,625 (13,540–13,688), FDP_c 0.0112.

### Percolator 3.7.1 against the native engine, same PINs (main)

| rescoring | slice 01 | slice 02 | FDP_c at q ≤ 0.01 | per run |
|---|---|---|---|---|
| Percolator 3.7.1 | 5669 (5657–5678) | 5798 (5751–5851) | 0.0110 / 0.0076 | 111–115 s wall, 0.45 GB |
| native, same RT-refined PINs | 5670 (5649–5700) | 5825 (5812–5861) | 0.0111 / 0.0074 | 2.6–3.5 s wall, 0.21 GB |
| native for both in-andes passes too | 5685 (5668–5707) | 5811 (5794–5821) | 0.0113 / 0.0078 | about 160 s per slice in total |

The two engines agree within seed noise. The native output is byte-identical between two runs
and between 1 and 4 threads. With Percolator 3.7.1 the two in-andes passes take 224–241 s per
slice, more than the search.

## What was kept

Three changes are on the branch; none of them changes the number of identifications beyond
noise. Percolator 3.7.1, N@FDP_c 1 %, final build (bdb64ee, flag off) against main (3f2b66c):

| | slice 01 | slice 02 | held-out 03 | held-out 5-min |
|---|---|---|---|---|
| main | 5669 (5657–5678) | 5798 (5751–5851) | 5770 (5714–5809) | 13,625 (13,540–13,688) |
| final build | 5676 (5633–5748) | 5825 (5760–5862) | 5814 (5786–5843) | 13,600 (13,467–13,731) |
| paired difference, seeds higher | +8, 3/5 | +27, 4/5 | +44, 4/5 | −24, 2/5 |
| N@FDP_c 2 %, main → final | 6040 → 6048 | 6104 → 6144 | 6096 → 6095 | 14,548 → 14,568 |
| FDP_c at q ≤ 0.01, main → final | 0.0110 → 0.0110 | 0.0076 → 0.0078 | 0.0079 → 0.0078 | 0.0112 → 0.0113 |
| search CPU, main → final | 504 → 459 s | 516 → 463 s | 509 → 462 s | 578 → 537 s |
| peak RSS, main → final | 2.44 → 2.04 GB | 2.46 → 2.06 GB | 2.42 → 2.04 GB | 2.69 → 2.28 GB |

1. **`RawScoreCal` over every charge of a charge-less spectrum, and no leave-one-out z-score from
   fewer than four other scores** (5be0a4e, a82536c). A DIA scan has no recorded charge, so the
   search runs one unit per charge; each unit calibrated the strong score against its own
   charge's candidates only, often two or three, and a leave-one-out z-score over two or three
   other scores has no finite variance. On slice 01, 23 of 123,474 PIN rows had
   |`RawScoreCal`| > 1000 (largest 2.5e5); after the fix the largest is 67 (slice 02: 1.1e4 →
   69). One entrapment precursor that had `RawScoreCal` 151,473 was accepted on every seed and
   is no longer accepted. The value changes on nearly every DIA row because the pool now spans
   the charges. The minimum pool size applies to every mode: `RawScoreCal` becomes 0 on 8 of
   1,718 rows of the chimeric golden PIN, 18 of 1,548 of the refine golden and 25 of 366 of a
   DDA golden against one protein (values of −14.75 to 28.33 from two or three other
   candidates), and nothing else changes in them; the four E. coli goldens are unchanged. The
   effect on DDA identifications was not measured.
2. **`rerank_pin` streams the PIN** (2227602) and **the predicted-MS2 row features are computed
   inside the prediction workers** (b2aad7b). The PIN is byte-identical before and after. Peak
   RSS falls by 16 % on a slice and by 20 % on a four times longer slice (RT 950–1230 s, 51,185
   MS2: 5.11 → 4.07 GB, CPU 1332 → 1229 s). What scales with the number of scans is the spectra
   and the 50-row queues (about 40 KB per MS2 scan on a fixed 1.2 GB), the PIN before the
   re-rank and the predictions.

Where the time goes on a slice (main, CPU): search 56 % (PSM features are 71 % of it, the GBDT
alone 31 %), andes-ms2 inference and its features 30 %, everything else under 15 %.

## Experimental: recentring the fragment window (`--dia-recentre-fragments`, hidden, off)

On Orbitrap Astral `--dia` matches fragments within 7 ppm of the theoretical m/z. The runs have
a systematic offset: −3.25 / −2.27 ppm on slices 01 / 02 and −1.93 / −1.43 ppm on the held-out
slices, with a spread of 1.1–1.3 ppm, independent of fragment m/z and constant over a slice.
The per-run measurement that DDA uses never fired under `--dia` (a scan's best candidate rarely
reaches 1 % q, and the 7 ppm default counted as a user override). The flag measures the offset
on 3,000 spectra spread over the run, with the best 10 candidates of each scan competing, and
moves the window by it for retrieval, scoring, the PIN features, the predicted-spectrum columns
and the co-elution columns.

Percolator 3.7.1, flag on against flag off, with the 95 % interval of a bootstrap over
precursors ([`bootstrap.py`](bootstrap.py)):

| | slice 01 | slice 02 | held-out 03 | held-out 5-min | pooled |
|---|---|---|---|---|---|
| N@FDP_c 1 %, off | 5676 (5633–5748) | 5825 (5760–5862) | 5814 (5786–5843) | 13,600 (13,467–13,731) | |
| N@FDP_c 1 %, on | 5906 (5856–6006) | 5840 (5803–5878) | 5737 (5698–5780) | 13,723 (13,654–13,830) | |
| paired difference, seeds higher | +230, 5/5 | +16, 4/5 | −77, 0/5 | +122, 4/5 | +291 |
| 95 % interval | [+73, +402] | [−127, +150] | [−182, +77] | [−190, +344] | [−81, +641] |
| N@FDP_c 2 %, difference | +201 [+88, +330] | +33 [−67, +157] | +31 [−63, +121] | +11 [−136, +175] | +276 [+63, +542] |
| FDP_c at q ≤ 0.01, off → on | 0.0110 → 0.0095 | 0.0078 → 0.0101 | 0.0078 → 0.0092 | 0.0113 → 0.0117 | |
| search CPU, off → on | 459 → 540 s | 463 → 557 s | 462 → 542 s | 537 → 585 s | |

- The rule set before the measurement was: higher on both tuning slices beyond seed noise, FDP_c
  at q ≤ 0.01 not higher, confirmed on the held-out slices. It fails: slice 01 gains, the other
  three are inside the noise of a slice, and run 03 is lower on every seed. So it is off by
  default.
- Pooled over the four slices it is +0.9 % at 1 % FDP with an interval that includes zero, and
  +0.8 % at 2 % with one that does not. The gain follows the size of the offset. This needs the
  full runs, where the counts are five to ten times larger, to decide.
- The measurement costs about 25 s wall per slice (60–75 s CPU) on a fixed 3,000 spectra; peak
  RSS is unchanged.
- Recentring also changes what the mass-error columns mean: the signed mean fragment error loses
  its target/decoy separation and its Percolator weight, the absolute error and the Gaussian
  mass-error score gain. Percolator refits per run, so this is not a mismatch with a fixed model.
- Without the co-elution columns (native): slice 01 +222 (5/5), slice 02 −1 (3/5), run 03 −130
  (0/5), 5-min +99 (5/5). Removed from the code.

## What was tested and dropped

### Rescoring variants (native engine on main's PINs; reference 5670 / 5825, FDP_c 0.0111 / 0.0074)

| variant | slice 01 (paired, seeds higher) | slice 02 (paired, seeds higher) | FDP_c at q ≤ 0.01 |
|---|---|---|---|
| `--trainFDR 0.01` | 5701 (+30, 4/5) | 5745 (−80, 0/5) | 0.0105 / 0.0080 |
| `--trainFDR 0.02` | 5677 (+7, 4/5) | 5776 (−48, 1/5) | 0.0114 / 0.0077 |
| `--trainFDR 0.10` | 5693 (+23, 4/5) | 5814 (−11, 3/5) | 0.0112 / 0.0077 |
| best row per precursor by RawScore before Percolator (89.5k of 123k rows) | 5599 (−72, 1/5) | 5690 (−134, 0/5) | 0.0093 / 0.0082 |
| best two rows per precursor (105k rows) | 5721 (+50, 4/5) | 5851 (+27, 4/5) | 0.0105 / 0.0084 |
| `RawScoreCal` set to 0 | 5696 (+26, 5/5) | 5825 (+0, 3/5) | 0.0108 / 0.0075 |
| mean score of three seeds | 5671 (+1) | 5834 (+9) | 0.0112 / 0.0080 |

- `--trainFDR` stays 0.05. Its original reason (few PSMs at 1 %) no longer holds, the first pass
  has about 21,600 PSMs at q ≤ 0.01, but lower values lose on slice 02.
- Reducing the PIN to one row per precursor loses on both slices. Two rows gain less than seed
  noise and raise FDP_c at q ≤ 0.01 on slice 02 (also on the fixed PINs, 0.0078 → 0.0086).
- Averaging seeds neither raises the mean nor narrows the spread.
- FDP_c above q on slice 01 (0.0110 at q ≤ 0.01) is within what about 31 entrapments allow (a
  binomial standard deviation of about 0.002). The entrapment hits concentrate at charges 5+, 1+
  and 4+ and in precursors with a single PIN row (FDP_c 0.031 / 0.022 against 0.002–0.006 with
  four or more rows).

### Ten more peak-group columns (not merged)

Columns added to the predicted-MS2 evidence of each row: spread and offset of the fragment mass
errors, the M+1 isotope of the fragments, apex offset, width and contrast of the median elution
profile, correlation of the MS1 monoisotopic and M+1 traces with that profile, presence of the
three most intense predicted fragments, and apex signal-to-noise.

| variant | slice 01 (paired, seeds higher) | slice 02 (paired, seeds higher) | FDP_c at q ≤ 0.01 |
|---|---|---|---|
| all ten, Percolator 3.7.1 (reference 5669 / 5798) | 5643 (−26, 1/5) | 5794 (−4, 3/5) | 0.0112 / 0.0068 |
| all ten, native (reference 5670 / 5825) | 5652 (−18, 2/5) | 5806 (−18, 1/5) | 0.0112 / 0.0059 |
| five leave-one-group-out sets, native | −24 … +26 | −31 … +12 | 0.0107–0.0118 / 0.0055–0.0078 |
| five single-group sets, native | −11 … +38 | −30 … +21 | 0.0103–0.0114 / 0.0060–0.0074 |

- No set moves N@FDP_c 1 % beyond seed noise, and every gain on one slice reverses on the other.
  Twelve sets were tried (all, none, five leave-one-out, five single); no other subsets.
- Alone, the columns separate accepted targets from decoys (AUC up to 0.95 for the profile
  contrast), and equally from entrapments, so there is no leakage. Percolator gives them weights
  near zero (|w| ≤ 0.06 except the apex signal-to-noise at −0.14, against +0.34 for `pf_corr`):
  they restate what `pf_corr`, `proj_angle`, the summed-apex and the MS1 columns already carry.
- Cost: +4 % search CPU, no memory.

### Pooling the other per-spectrum score terms over charges (not merged)

A DIA scan has no recorded charge, so the search runs one unit per charge and each unit computes
its per-spectrum quantities over its own candidates. After `RawScoreCal` (kept, above) four more
remain: the listwise score gap, the rank entropy, the delta score and the tailor denominator. The
first two enter RawScore.

| variant (native; reference is the kept build) | slice 01 | slice 02 | held-out 03 | held-out 5-min |
|---|---|---|---|---|
| reference | 5687 (5645–5760) | 5832 (5766–5881) | 5778 (5722–5819) | 13,602 (13,537–13,698) |
| gap and entropy pooled | +46 (4/5) | +65 (5/5) | +17 (3/5) | +16 (2/5) |
| all four pooled | +72 (4/5) | +44 (4/5) | – | – |
| delta and tailor pooled | −7 (1/5) | −24 (1/5) | – | – |

Pooling gap and entropy looked consistent on the tuning slices but stayed inside seed noise, and
the held-out slices do not confirm it; on the 5-min slice FDP_c at q ≤ 0.01 rises from 0.0115 to
0.0126. Dropped.

### Narrowing the fragment window to the measured spread

With the window recentred (above) and narrowed to four measured spreads (5.2 / 5.3 ppm instead of
7 ppm), native: slice 01 +141 (5/5), slice 02 −120 (0/5), FDP_c at q ≤ 0.01 0.0097 / 0.0118.
Dropped; the mode was removed from the code.

## Limits of this benchmark

- A slice is 70 s, about 42 cycles. 24 % of the PIN rows have fewer than five cycles on one side,
  so the ±5-cycle columns lose context, and retention-time features have little range (the RT
  refinement moves the median error from 0.19 to 0.16 min). Changes that act through retention
  time are under-tested.
- The five Percolator seeds share one PIN. Their range covers the training noise but not the
  noise from which entrapment precursors score near the threshold: a slice has 23–31 entrapments
  above q ≤ 0.01 (78 on the 5-min slice).
  A bootstrap over precursors gives the paired difference of two builds a standard deviation
  of 60–90 precursors on a 15-min slice and 120–140 on the 5-min slice, so one slice resolves a
  change of 2–3 % and the four together about 1.2 %. As a check, #402 → #410 comes out at
  +518 [+375, +703] and +334 [+183, +503], and #417 → main at +11 [−90, +72] and −23 [−75, +38].
- All four runs are one instrument, one sample type and one acquisition scheme.
- Quantification, diaPASEF and wide-window DIA are not covered.

## Reproduce

```bash
# BENCH_ROOT holds data/, tools/, bin/ and runs/; ANDES_SRC is the andes checkout.
P=https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/02/PXD070049
for f in LFQ_Astral_DIA_15min_50ng_Human_0{1,2,3} LFQ_Astral_DIA_5min_50ng_Human_01; do
  curl -fL -C - -o "$BENCH_ROOT/data/raw/$f.raw" "$P/$f.raw"
  ThermoRawFileParser -i "$BENCH_ROOT/data/raw/$f.raw" -o "$BENCH_ROOT/data/mzml" -f 1 -l 2
done
curl -fL -o human_reviewed.fasta.gz \
  'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29'
python3 -I build_library.py human_reviewed.fasta.gz "$BENCH_ROOT/data/frozen/pairlib_public.fasta" \
  "$BENCH_ROOT/data/frozen/pairlib_public_lib.tsv" --n 473461
for r in 01 02 03; do
  python3 -I slice_mzml.py "$BENCH_ROOT/data/mzml/LFQ_Astral_DIA_15min_50ng_Human_$r.mzML" \
    "$BENCH_ROOT/data/frozen/slice${r}_950_1020.mzML" 950 1020
done
python3 -I slice_mzml.py "$BENCH_ROOT/data/mzml/LFQ_Astral_DIA_5min_50ng_Human_01.mzML" \
  "$BENCH_ROOT/data/frozen_heldout/slice5min01_545_615.mzML" 545 615
sha256sum "$BENCH_ROOT"/data/frozen*/*    # compare with SHA256SUMS

./build_variant.sh main "$ANDES_SRC"                 # -> $BENCH_ROOT/bin/andes-main
./bench.sh base-main "$BENCH_ROOT/bin/andes-main" full
BENCH_FROZEN=$PWD/frozen_heldout.env ./bench.sh --raws "03 5m" held-main "$BENCH_ROOT/bin/andes-main" full
python3 -I compare.py base-main <tag>
```

Percolator 3.7.1 is expected at `$BENCH_ROOT/tools/env-perc/bin/percolator` (bioconda
`percolator=3.7.1`), or pass `--perc`. The UniProt download is the current release; the library
is only byte-identical to the one used here from release 2026_04 (sha256 of the gzipped FASTA
97e06262…5b48d).

## Open

- Recentring on the full runs (above).
- The search computes its expensive PSM features for all 50 rows of a scan and keeps 10 after
  the re-rank. Computing them for the kept rows only would save about 9 % of the CPU time
  (estimated from a build without them: 416 against 473 s), and predicting RT for the kept rows
  only about 4 %; both need the re-rank moved before those steps. Freeing the 40 dropped rows per
  scan after the re-rank is the largest memory item that scales with scans, but the run
  statistics and the TSV output read every row.
- The two in-andes Percolator 3.7.1 passes take longer than the search; the native engine gives
  the same counts in a few seconds.
