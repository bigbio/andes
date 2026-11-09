# andes benchmarks

andes is benchmarked on public PRIDE datasets for identification (high-res, TMT-labelled and
low-res LFQ, phospho, glyco) and for quantification (label-free with match between runs, and
TMT). Identifications are counted after Percolator at `q ≤ 0.01`. Where the database carries an
entrapment component, the page also gives the **true** false discovery proportion (FDP), which
is what tells you whether the nominal 1% holds. Other tools are named with the version compared.
Their numbers come either from runs on the same host with the same rescoring, or from a public
leaderboard. Section 2 gives the date and the andes commit behind each row.

## 1. Summary

Identification rows: `2051616b` (origin/main, 2026-10-23), default settings plus the dataset's mods,
Percolator 3.7.1 `-Y`, seeds 42, 1 and 2. Quantification rows were measured again on `09e1c00d`
(origin/main, 2026-11-05): the label-free rows (ProteoBench, truth designs) on `09e1c00d`, the
TMT rows on `0cfa03cc` (origin/main after #320 and #321); the one exception is the 500-run study,
which still carries its `2051616b` numbers. The later build changes only ion-trap HCD
identification and the Percolator PEP handling, neither of which enters the label-free rows
(Orbitrap and Astral data). Comparison rows keep the date they were measured.

### At a glance: where andes is strong and where it needs work

Verdicts compare andes with the best other tool measured on the same data: **well ahead**,
**ahead**, **comparable**, **behind**, **well behind**. 🟢 marks a strength, 🔴 a gap we are working on.

🟢 **Strong**
- Identification depth: well ahead of Comet on every standard set (Astral +66%, TMT low-res +29%,
  phospho +14% at a lower true FDP, UPS1 +10%), at about half Comet's wall time.
- Quantification depth: well ahead on 7 of 8 truth sets (+19% to +51% precursors over the best
  other tool); +54% PSMs on TMT PXD007683.
- ProteoBench QExactive: lower median error than the best leaderboard entry (0.161 against 0.190).
- DE true positives (common-tested family, sign-aware): ahead on QE (+22% over MaxQuant, +11% over FragPipe, with more false calls than both); tie on the UPS1 and Astral sets against MaxQuant as reported. MaxQuant 2.7 does not remove its match-between-runs rows; with its matches filtered at Match q-value ≤ 0.01 (83–94% of its transfers removed) andes is ahead on true positives on all four sets (+5% to +15%) and filtered MaxQuant has the fewest false calls on all four (andes 1.9–2.5x).
- TMT match between runs: the transferred values are as accurate as identified ones (excess −0.05
  points), against +1.53 (MaxQuant) and +0.40 (SIMSI-Transfer); `--tmt-unmix on` removes hrMS2
  interference.
- Protein CV: ahead on 3 of 6 truth sets.
- Precision of unchanged proteins: ahead on PXD009740, PXD009815 and QE (within-condition spread of MaxQuant 1.15–1.55x andes', FragPipe 1.25–1.40x).

🔴 **Needs work**
- DE false calls without a fold-change cut: well behind MaxQuant on PXD009740 and QE (2.4x and 3.0x), behind on PXD009815 and Astral (1.7x and 1.4x). Most are small shifts of the "unchanged" species that andes' tighter protein values make significant: above |log2 FC| ≥ 0.5 andes is ahead on QE and comparable on PXD009815; on PXD009740 296 of 320 come from comparisons against the highest UPS1 level, without which andes is level with MaxQuant (see *Why andes makes more false calls*).
- Gross error of bright identified cells: well behind (1.2–7.2% against MaxQuant 0.3–1.6%).
- False transfers into absent samples at the default threshold: behind (0.20% of transferred rows
  against FragPipe 0.09%; comparable at `quant_q ≤ 0.01`; on the classified reading behind on
  enriched rows, 0.052% against 0.010%; counting FragPipe's rows on precursors andes does not
  report gives 0.052% against 0.072%, but that reading uses andes' own table as the reference and is
  not neutral).
- TMT ratio accuracy (PXD007683): behind MaxQuant (protein MAE 0.213 against 0.159) and hrMS2
  interference behind (0.192 against 0.156).
- ProteoBench Astral: CV behind most submissions (0.112 against 0.070–0.095); median error behind
  the submissions that report fewer precursors (MaxQuant 0.160 at 20,283, Sage 0.167 at 30,319).
- Glyco recall of reference identifications: behind (79–89% of pGlyco2's).
- Error control: on a peptide-level shuffled entrapment UPS1 reads 0.92–0.94% at nominal 1% and
  timsTOF 0.79–0.91% (mean of six shuffle draws per run; single draws 0.52–1.12%), on target. The natural entrapments read
  higher (UPS1 1.95–2.07%, timsTOF up to 1.54%); about half of those hits fall on proteins with a
  close homolog in the target proteome (see the UPS1 and timsTOF paragraphs).

**Identification** (PSMs at Percolator `q ≤ 0.01`, 8 threads):

Legend: 🟢 andes better · 🟡 similar · 🔴 andes worse · ⚪ no comparison. The last column says whether a gap can be closed downstream of andes (filtering on its quality columns, a fold-change cut, normalisation or roll-up) or needs a change in the engine, and how much the gap matters: **minor** (usable as is), **moderate** (usable with the care stated), **major** (default results are hurt badly; use the fix). No result is unusable, and no gap is major. "Best other tool" is the strongest result another tool reported on the same data for that metric; percentages are andes relative to it.

| benchmark (data) | metric | andes | best other tool | result | impact · can it be fixed downstream? |
|---|---|---:|---|---|---|
| Astral HCD, PXD070049 (1 file) | PSMs at q ≤ 0.01 | 52,169 | Comet 2025.01: 31,435 | 🟢 +66% | — |
| Astral HCD | wall time (8 threads) | 114–117 s | Comet 2025.01: 217 s | 🟢 about 2x faster | — |
| Astral HCD | true FDP at nominal 1% | 0.89–1.08% | — | 🟢 on target | — |
| TMT low-res, PXD007683 a05058 (ion-trap CID) | PSMs | 13,516 | Comet 2025.01: 10,504 | 🟢 +29% | — |
| TMT high-res, PXD017803 TKO hrMS2 | PSMs | 6,628–6,686 | not measured | ⚪ | — |
| Low-res LFQ, UPS1 PXD001819 (ion-trap CID) | PSMs (`2ee76dd6`, 2 missed cleavages; was 16,225 at `2051616b`) | 16,632–16,675 | Java MS-GF+: 15,904 | 🟡 +5% | — |
| Low-res LFQ, UPS1 | true FDP, peptide-level shuffled entrapment (primary) | 0.92–0.94% | — | 🟢 within nominal 1% | — |
| Low-res LFQ, UPS1 | true FDP, natural *E. coli* entrapment (context, not a valid null here) | 1.95–2.07% (1.07–1.11% without the conserved proteins, T/E recomputed) | — | ⚪ context: not a calibration verdict; the shuffle above is the primary reading (0.92–0.94%, within nominal) | none. 45 of 94 hits (seed 42) sit on 8 *E. coli* proteins with a yeast homolog at ≥ 40% identity, 32 of them on two EF-Tu peptides. These peptides do not occur in any yeast protein, so they are either real *E. coli* protein in the recombinant UPS1 standard or yeast spectra matched to a near-identical *E. coli* peptide; neither has been shown |
| Phospho, PXD007653 (Q Exactive) | PSMs | 38,623 | Comet 2025.01: 33,984 (at 1.77% true FDP) | 🟢 +14%, at 0.80–0.99% FDP | — |
| timsTOF DDA-PASEF, PXD053748 and PXD064202 (4 runs) | PSMs at q ≤ 0.01 (shuffled-entrapment database, mean of 6 draws) | 15.6k–58.6k | Comet 2025.01: 14.1k–50.8k | 🟢 +11% to +17% | — |
| timsTOF DDA-PASEF | PSMs at a true FDP of 1% | 15.9k–58.8k | Comet 2025.01: 14.0k–49.5k | 🟢 +14% to +21% | — |
| timsTOF DDA-PASEF | wall time (8 threads) | 54–62 s | Comet 2025.01: 159–222 s | 🟢 about 3x faster | — |
| timsTOF DDA-PASEF | true FDP, peptide-level shuffled entrapment (primary) | 0.79–0.91% (mean of six draws per run; single draws 0.52–1.12%) | Comet 2025.01: 1.14–1.31% | 🟢 on target; Comet above nominal on all four runs | — |
| timsTOF DDA-PASEF | true FDP, natural yeast entrapment | 0.86–1.54% (0.72–1.02% without hits on conserved proteins) | — | 🟡 above nominal on one astrocyte run (1.44–1.54%) | Same conserved-homolog effect as UPS1; see the timsTOF paragraph |
| Glyco quick tier, PXD005553 | reference glycoPSMs confirmed | 88–90% | pGlyco2, MSFragger-Glyco (references) | 🔴 10–12% of reference IDs missed | **moderate**, usable. No: engine (glyco scoring) |
| Glyco deep tier, PXD005553 | reference glycoPSMs confirmed | 79% (pGlyco2), 88% (MSFragger-Glyco) | as above | 🔴 21% of pGlyco2's missed | **moderate**, usable. No: engine (glyco scoring) |

† On a cluster node (8 threads); the other walls are on the benchmark VM (8 cores), the host the
Comet times come from. Wall times compare only within one host.

**Quantification** (andes defaults, 2026-11-05: label-free rows at `09e1c00d`, TMT rows at `0cfa03cc`; the 500-run row 2026-10-23 at `2051616b`):

Legend: 🟢 andes better · 🟡 similar · 🔴 andes worse · ⚪ no comparison. The last column says whether a gap can be closed downstream of andes (filtering on its quality columns, a fold-change cut, normalisation or roll-up) or needs a change in the engine, and how much the gap matters: **minor** (usable as is), **moderate** (usable with the care stated), **major** (default results are hurt badly; use the fix). No result is unusable, and no gap is major. "Best other tool" is the strongest result another tool reported on the same data for that metric; percentages are andes relative to it.

| benchmark (data) | metric | andes | best other tool | result | impact · can it be fixed downstream? |
|---|---|---:|---|---|---|
| ProteoBench QExactive, PXD028735 (6 runs) | median \|ε\| (lower is better) | 0.161 | FragPipe 23.0: 0.190 | 🟢 15% lower error | — |
| ProteoBench QExactive | precursors in ≥ 3 runs | 63,796 | i2MassChroQ 1.0.16: 74,697 (\|ε\| 0.222); FragPipe 23.0: 63,701 | 🟡 similar to FragPipe; fewer than submissions with higher error | — |
| ProteoBench QExactive | CV | 0.117 | i2MassChroQ 1.2.0: 0.120 | 🟢 lowest | — |
| ProteoBench Astral, PXD070049 (6 runs) | median \|ε\| at similar depth (36–39k) | 0.186 | i2MassChroQ 1.2.6: 0.197; FragPipe 23.0: 0.212 | 🟢 lowest at this depth | — |
| ProteoBench Astral | median \|ε\|, any depth | 0.186 | MaxQuant 2.6.7.0: 0.160 (20,283 precursors); Sage 0.14.6: 0.167 (30,319) | 🔴 higher than tools reporting fewer precursors | **minor**, usable. Partly: a stricter transfer q trades depth for error, but identified rows alone also read 0.186, so the rest is in the engine (peak fit) |
| ProteoBench Astral | precursors in ≥ 3 runs | 36,675 | FragPipe 23.0: 39,057 | 🟡 similar | — |
| ProteoBench Astral | CV | 0.112 | Sage 0.14.6: 0.070 | 🔴 higher | **minor**, usable. Partly: identified rows alone give 0.100; normalisation changes little (0.113 raw) |
| Truth designs (8 sets) | precursors quantified | 10.7k–93.6k | MaxQuant 2.7 / FragPipe 21.1 per set | 🟢 +19% to +51% on 7 sets, 🟡 +3% on PXD001819 | — |
| Truth designs, PXD009740 / PXD009815 / QE / Astral | DE true positives (MBR on; common-tested family, sign-aware) | 1,300 / 844 / 1,047 / 685 | MaxQuant 1,304 / 844 / 857 / 650; FragPipe 1,237 / 849 / 946 / 689 | 🟡 tie / 🟡 tie / 🟢 +22% over MaxQuant, +11% over FragPipe / 🟡 tie | — |
| Truth designs, PXD009740 / PXD009815 / QE / Astral | DE false calls, no fold-change cut (common-tested family; wrong-sign calls included) | 305 / 153 / 110 / 296 | MaxQuant 129 / 92 / 37 / 210; FragPipe 302 / 301 / 43 / 389 | 🔴 more than MaxQuant on all four (2.4x / 1.7x / 3.0x / 1.4x) | **moderate**: apply a fold-change cut. With \|log2 FC\| ≥ 0.5 (seed 42, common family) andes is ahead on QE (11 against 18, with 1,046 against 852 true positives), comparable on PXD009815 (39 against 36) and behind on PXD009740 (167 against 74, all but a few against its highest UPS1 level) |
| Truth designs, PXD009740 / PXD009815 / QE / Astral | DE true positives / false calls against MaxQuant with its match rows filtered at Match q-value ≤ 0.01 (common-tested family of andes, filtered MaxQuant and FragPipe; sign-aware) | 1,156 / 272; 763 / 151; 983 / 105; 568 / 297 | MaxQuant, matches q ≤ 0.01: 1,077 / 131; 725 / 62; 851 / 42; 518 / 154. FragPipe: 1,089 / 274; 768 / 265; 905 / 41; 569 / 388 | 🟢 TP +7% / +5% / +15% / +10% over filtered MaxQuant; 🔴 false calls 2.1x / 2.4x / 2.5x / 1.9x of it | **moderate**: the filtered MaxQuant arm is the cleanest; andes' extra true positives come with more false calls |
| Truth designs, all sets | bright identified cells off by > 2x | 1.2–7.2% | MaxQuant 0.3–1.6% | 🔴 | **moderate**, few cells, usable. Partly: quality columns let a downstream weighting down-weight them; the fix is in the engine (template misfit on fronting or tailing peaks) |
| Truth designs, 6 sets | protein CV | 4.4–10.1% | FragPipe 4.9–13.5%, MaxQuant 7.0–11.6% | 🟢 on 3 sets, 🔴 on Astral and PXD001819 | **minor**, usable. Partly: dropping transfers or weighting by quality columns |
| UPS1 lowest three levels, PXD009740 (170/382/859 amol) / PXD009815 (10/50/100 amol) | precursors within 2x of truth, as reported | 502 / 24 | MaxQuant 417 / 23; FragPipe 458 / 50 | 🟢 ahead of both on PXD009740 and 🔴 behind FragPipe on PXD009815 as reported; at equal rows FragPipe leads on both (458 against 437; 50 against 22) and andes leads MaxQuant on PXD009740 (437 against 341) | **minor**, lowest amounts only. No: engine (faint-transfer detection) |
| Known-absent yeast, PXD014415 | false-transfer rate, all yeast rows (classified: enriched only 0.052%; with no-reference rows 0.052%) | 0.20% (0.08% at `quant_q ≤ 0.01`) | FragPipe 21.1: 0.09% (enriched only 0.010%; with no-reference rows 0.072%) | 🔴 at default, 🟡 at `quant_q ≤ 0.01`; classified: 🔴 on enriched rows (the no-reference reading uses andes' table as reference, not neutral) | **minor**, usable. Yes: filter `quant_q ≤ 0.01` (0.08%, at fewer transfers) |
| Blank run, PXD016921 | rows reported in the blank | 150 (85 identified, 65 transferred) | MaxQuant 105 (20 / 85, unfiltered matches; 29 (20 / 9) at Match q-value ≤ 0.01); FragPipe 21.1: 8 (8 / 0) | 🔴 | **moderate** for low-input or blank-adjacent runs. No: `quant_q ≤ 0.01` still leaves 144; mostly carry-over (real ions) |
| TMT knockout, PXD017803 SPS / RTS / hrMS2 | interference (lower is better) | 0.062 / 0.027 / 0.192 | MaxQuant 2.7: 0.054 / 0.038 / 0.156 | 🟡 / 🟢 / 🔴 (🟢 0.000 with `--tmt-unmix on`) | **moderate** on hrMS2, usable with the option. Yes for hrMS2: `--tmt-unmix on` (0.000) or a purity ≥ 0.9 filter (0.135) |
| TMT knockout, hrMS2 | PSMs with reporter ions evaluated (scans whose range records the reporters) | 3,278 | FragPipe 21.1: 3,038; MaxQuant 2.7: 2,241 | 🟢 | — |
| TMT known ratios, PXD007683 | PSMs | 149,121 | MaxQuant 2.7: 97,107 | 🟢 +54% | — |
| TMT known ratios, PXD007683 | protein ratio error (MAE A/C) | 0.213 | MaxQuant 2.7: 0.159 | 🔴 | **moderate**, usable. Partly: a stricter purity filter; the rest is extraction and roll-up |
| TMT match between runs, PXD014750 hrMS2 | excess wrong-value rate of transfers | −0.05 points | SIMSI-Transfer 0.8: +0.40; MaxQuant isobaric MBR: +1.53 | 🟢 | — |
| TMT match between runs | cells transferred | 3,636 | SIMSI-Transfer 0.8: 16,511 | 🔴 fewer | **minor** (opt-in mode). Yes: a looser `quant_q` keeps all 7,841 candidates at +0.40 points, the same as SIMSI-Transfer |
| Large study, PXD042233 (500 runs) | precursors in ≥ 3 runs | 243,895 | not measured | ⚪ | — |

Settings, error control and cost of every row: section 2.

**Opt-in modes** (Astral, PSMs / distinct peptides, Percolator seeds 42, 1 and 2; measured
2026-10-05 to 2026-10-13 on earlier builds and not re-measured at `2051616b`):

| mode | andes | without the mode |
|---|---|---|
| `--chimeric` (second and later peptides per spectrum), at an entrapment true FDP of 1% | 66,377–66,424 / 29,384–29,400 | 43,076–43,145 / 26,723–26,746 |
| `--refine` (PTM discovery, high-res only) | 48,312 PSMs, 28,728 peptides at `q ≤ 0.01`; its second-pass PSMs at ~1.2–1.8% FDP | 47,080 / 28,932 |

## 2. The benchmarks

**Where and how.** Every andes identification number on this page was measured on 2026-10-23 at
`2051616b`, except the UPS1 and timsTOF entrapment rows, which were re-measured on 2026-11-07 at
`2ee76dd6` (see their paragraphs), and every quantification number on 2026-11-05 at `09e1c00d` (label-free) or
`0cfa03cc` (TMT; the 500-run study excepted), unless a row says otherwise. PSM counts, wall time, CPU time and peak memory of the four standard
sets (Astral, TMT a05058, UPS1, phospho) come from the benchmark VM (8 cores, 8 threads,
`/usr/bin/time`, two runs except phospho), on the plain databases the Comet rows used. The true
FDP of every set and all other rows come from a cluster node (8 threads), with each target
database extended by an entrapment component: three peptide-level shuffled realisations
(`build_shuffled_entrap.py` seeds 1, 2 and 3) for Astral, phospho and glyco, one realisation for
UPS1, TMT, TKO and timsTOF, plus a foreign (reviewed *Arabidopsis*) or natural entrapment where
one was used before. Percolator 3.7.1 `-Y --num-threads 8`, seeds 42, 1 and 2. "At a true FDP
of 1%" is the deepest score-sorted list whose entrapment-estimated FDP stays at or below 1%.

**Standard, high-res (Astral).** One 15-minute DDA run of the ProteoBench human/yeast/E. coli
sample against the 31,889-entry HYE FASTA. 52,169–52,437 PSMs and 31,147–31,238 peptides at
`q ≤ 0.01`. With the peptide-level entrapment (T/E 1.075) the nominal 1% reads 0.89–1.08% true
FDP over nine database × seed combinations, and andes keeps 47,267–47,977 PSMs and
28,141–28,857 peptides at a true FDP of 1%; with the *Arabidopsis* entrapment (T/E 2.228),
0.99–1.05% and 49,965–50,185 PSMs. 114–117 s wall, 569–636 CPU-s, 4.8 GB. Comet 2025.01 was run
on the same VM with matched parameters and the same Percolator rescoring on 2026-10-05.

**TMT-labelled, low-res.** One ion-trap CID run of a TMT 10-plex sample against human + yeast
reviewed proteins (26,483 sequences): 13,443–13,516 PSMs and 11,852–11,922 peptides at
`q ≤ 0.01`; 40–47 s, 200–216 CPU-s, 1.7 GB. With a peptide-level entrapment (T/E 1.073) the
nominal 1% is 0.67–0.76% true FDP (13,325–13,371 PSMs at a true FDP of 1%); with *Arabidopsis*
(T/E 2.021), 0.83–1.05% (13,161–13,388).

**TMT-labelled, high-res (TKO hrMS2).** `Fig2_m23142_hrMS2` of PXD017803 (TMTpro, Orbitrap
MS2), `--protocol TMT`, against reviewed yeast + cRAP with a peptide-level entrapment (T/E 1.039):
6,628–6,686 PSMs at `q ≤ 0.01`, 0.98–1.31% true FDP, 6,607–6,662 PSMs at a true FDP of 1%.

**Low-res LFQ (UPS1).** One run of the UPS1 standard in yeast (`UPS1_5000amol_R1`), measured on
2026-11-07 at `2ee76dd6` with the documented defaults (no extra flags; 2 missed cleavages, the
default since #300), Percolator 3.7.1 `-Y` seeds 42, 1 and 7, entrapment evaluated with the
all-proteins rule. Two entrapment designs, each searched as its own database:

| design | T/E | PSMs at `q ≤ 0.01` | peptides | entrapment hits | true FDP at nominal 1% | PSMs at a true FDP of 1% |
|---|---:|---:|---:|---:|---:|---:|
| peptide-level shuffle of the yeast targets (primary) | 1.033 | 16,567–16,623 | 4,176–4,200 | 75–77 | **0.92–0.94%** | 16,617–16,733 |
| natural *E. coli* (4,531 `ENTRAP_` proteins) | 2.451 | 16,632–16,675 | 4,116–4,227 | 94–100 | 1.95–2.07% | 14,786–15,045 |
| natural, hits on conserved proteins removed | 2.451 (2.630 with those proteins removed from E) | same | same | 49–51 | 1.01–1.06% (1.07–1.11%) | 16,505–16,607 |

Cost on 8 threads: andes 13–15 s wall, 71–79 CPU-s, 1.4 GB; Percolator (three seeds) 38–40 s
each, about 190–200 CPU-s, 1.2 GB. Verdict on the shuffle: within nominal (0.92–0.94%). Verdict on
the natural design: above nominal (1.95–2.07%), and that reading is not a clean calibration:

- **Conserved proteins.** An `ENTRAP_` protein is called conserved when its best local alignment
  (BLOSUM62, gap open 11, extend 1, against the yeast target proteins preselected by shared
  4-mers) covers at least 60% of it at 40% identity or more; 269 of the 4,531 *E. coli*
  proteins qualify (719 at 30% / 50%; 23 at 60% / 60%). The 40% / 60% cut is a conventional
  ortholog-level threshold fixed in advance, not tuned; hit counts were not re-run at the other
  two cuts, only T/E without those proteins (2.872 at 30% / 50%, 2.456 at 60% / 60%).
- **Hit breakdown (seed 42, 94 hits on 41 protein groups).** 8 conserved groups hold 45 hits.
  The largest is EF-Tu (`EFTU1`/`EFTU2`, 65% identical to yeast EF-Tu): 32 PSMs from only two
  distinct peptides. Next are ornithine carbamoyltransferase (7 PSMs, one peptide) and enolase
  (5 PSMs). Elsewhere a protein typically has one hit. Repeated PSMs on the same few peptides of
  abundant, conserved proteins look like real material, not random errors. UPS1 is a recombinant
  standard expressed in *E. coli*, so host proteins such as EF-Tu can be present in the sample
  (the same pattern appeared in an earlier analysis, where EF-Tu and DnaK made up 40 of 100
  hits). This is consistent with, not proof of, real *E. coli* material; the data cannot separate
  it from a mis-assigned yeast peptide.
- **Shared peptides.** Counting a PSM as entrapment when any of its proteins is entrapment (the
  rule used before the shared-peptide fix) reads 3.09–3.23% on the same runs, because peptides
  shared with a yeast protein become false hits. The docs and the evaluator count a hit only
  when all proteins are entrapment and the I/L-folded peptide occurs in no target protein.
- **T/E.** 2.451 is computed from the actual database with the andes digest (2 missed
  cleavages, 6–50 residues, I/L folded, no proline rule), after moving the 2,486 entrapment
  peptides that occur in a yeast protein to the target half. At 1 missed cleavage (the digest of
  the earlier measurement) it is 2.351; with a proline rule 2.370 (1 missed) or 2.476 (2); at 3
  missed cleavages 2.550. A digest giving 2.524 was not reproduced. Every FDP here uses 2.451.
- **Old default digest.** With `--max-missed-cleavages 1` (the 2026-10-23 setting), same build:
  16,238–16,367 PSMs, natural FDP 1.94–1.99% (T/E 2.351), shuffle 0.88–0.96% (T/E 1.039).

Java MS-GF+ (v20240326, earlier session, not re-run) reported 15,904 PSMs on this file. Earlier
measurement (`2051616b`, 1 missed cleavage, seeds 42, 1 and 2): 16,049–16,242 PSMs, natural
2.07–2.19%, shuffle 0.91–0.98%.

**Phospho-enriched.** One single-shot EasyPhos mouse liver file
(`20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw`, as mzML), not in any model's
training data, 80 Da index slices. Against mouse reviewed proteins plus the whole-protein shuffled
entrapment that the Comet row used: 38,596–38,629 PSMs (27,693–27,724 phospho-bearing) and
20,422–20,461 peptides, 337 s, 2,160 CPU-s, 7.9 GB. Comet 2025.01 reported 33,984 at 1.77% true
FDP on 2026-10-06 (on the whole-protein shuffle, which can read low). andes's true FDP on three
peptide-level shuffles (T/E 1.061) is 0.80–0.99%, with 38,604–38,843 PSMs at a true FDP of 1%;
on *Arabidopsis* (T/E 1.403), 1.22–1.36% and 38,298–38,412. Phosphosite localisation is
**not** scored. Against the depositors' MaxQuant identifications (PEP ≤ 0.01, 23,563 scans)
andes covered 82.4% with 95.3% same-scan backbone agreement on 2026-10-07 (`40774aca`); not
re-measured.

**timsTOF DDA-PASEF.** Native `.d` input, 20 ppm precursor tolerance, 2 missed cleavages,
`--top-n 1`, 80 Da index slices, Carbamidomethyl C fixed, Oxidation M and protein N-terminal
Acetyl variable, against reviewed human proteins with a reviewed yeast entrapment proteome
(natural, T/E 3.784) or a peptide-level shuffle of the human targets (T/E 1.050). Measured on
2026-11-07 at `2ee76dd6`, Percolator 3.7.1 `-Y` seeds 42, 1 and 7, all-proteins rule, T/E at
2 missed cleavages (3.845 at 1). "Without conserved" removes hits on yeast proteins that align to
a human protein at 40% identity or more over 60% of their length (907 of 6,733), with the same T/E.

| run | shuffle: PSMs at `q ≤ 0.01` | shuffle: true FDP | shuffle: PSMs at a true FDP of 1% | natural yeast: PSMs / true FDP / without conserved / at 1% |
|---|---:|---:|---:|---|
| PXD053748 `..._4363` (plasma, 30 min) | 15,411–15,474 | 0.52–0.60% | 15,938–15,952 | 16,062–16,108 / 0.86–1.01% / 0.80–0.95% / 16,076–16,219 |
| PXD053748 `..._4367` (plasma, 50 min) | 18,900–18,957 | 0.82–0.94% | 19,031–19,152 | 19,456–19,507 / 0.91–0.98% / 0.79–0.86% / 19,528–19,576 |
| PXD064202 `..._2442` (astrocytes) | 36,052–36,125 | 0.94–1.01% | 36,047–36,178 | 36,873–36,982 / 1.44–1.54% / 0.93–1.02% / 36,040–36,071 |
| PXD064202 `..._2403` (astrocytes) | 58,814–58,873 | 1.12% | 58,542–58,637 | 59,739–59,761 / 1.17–1.24% / 0.72–0.78% / 59,068–59,260 |

Verdict on the shuffle: within nominal. The table is one shuffle draw, which the three Percolator
seeds re-score; on that draw `..._2403` reads 1.12%. Five further draws (same rules, T/E 1.050,
Percolator seeds 42 and 7, 2026-11-09) read 0.82–0.95% on `..._2403` and 0.70–0.96% on `..._2442`;
over six draws the means are 0.91% (SD 0.11) and 0.88% (SD 0.10), a spread that matches counting
noise on 250–320 entrapment hits. Per-run q-values (#333) do not change these numbers: each timsTOF
search holds one run. On the natural design `..._2442` is above nominal (1.44–1.54%); removing hits
on conserved proteins brings it to 0.93–1.02%. The conserved-protein pattern is the same as on
UPS1: at seed 42 on `..._2442`, 14 conserved groups hold 39 of 111 hits, among them yeast Ssb1/2
(15 PSMs from two distinct peptides, 64% identical to human HSP70), GAPDH and EF-Tu. A human
sample cannot contain yeast proteins, so these hits are best read as near-identical human
peptides assigned to the yeast homolog, not independent errors. The any-protein rule would
read 9.0–9.7% on the astrocyte runs. andes 40–52 s per run, 255–321 CPU-s, 2.1–2.5 GB on the
natural database (49–58 s, 307–371 CPU-s, 2.4–2.8 GB on the shuffle); Percolator is separate.
Earlier measurement (`2051616b`, 1 missed cleavage on the same inputs): natural yeast
0.78–1.55%, shuffle 0.62–1.03%.

**Against Comet 2025.01 (timsTOF).** The four runs were also searched with Comet 2025.01 on mzML
converted with tdf2mzml (default options): 20 ppm precursor, isotope error −1 to 3 (Comet has no
−1 to 2 setting), 2 missed cleavages, the same modifications, `fragment_bin_tol 0.02` with offset 0
(the convention of the other Comet rows), internal reversed decoys. Both tools searched the same
shuffled-entrapment database for six draws (seeds 1 and 11–15), and both PINs were rescored with
Percolator 3.7.1 `-Y` seed 42 (seed 7 moves the true FDP by at most 0.2 pp); andes build `2ee76dd6`,
8 threads each, on 2026-11-09. Means over the six draws:

| run | tool | PSMs at `q ≤ 0.01` | true FDP | PSMs at a true FDP of 1% | peptide sequences | wall | memory |
|---|---|---:|---:|---:|---:|---:|---:|
| `..._4363` | andes `.d` | 15,642 | 0.79% | 15,880 | 3,910 | 55 s | 2.4 GB |
| | Comet | 14,143 | 1.14% | 13,964 | 3,493 | 159 s | 6.2 GB |
| `..._4367` | andes `.d` | 18,916 | 0.90% | 19,133 | 3,865 | 56 s | 2.5 GB |
| | Comet | 16,139 | 1.15% | 15,865 | 3,396 | 176 s | 6.4 GB |
| `..._2442` | andes `.d` | 36,092 | 0.88% | 36,279 | 14,682 | 54 s | 2.3 GB |
| | Comet | 30,848 | 1.31% | 29,957 | 13,248 | 186 s | 6.0 GB |
| `..._2403` | andes `.d` | 58,578 | 0.91% | 58,828 | 24,146 | 62 s | 2.7 GB |
| | Comet | 50,814 | 1.31% | 49,517 | 22,007 | 222 s | 6.1 GB |

Verdict: ahead on every run. andes reports 10.6–17.2% more PSMs at `q ≤ 0.01` and 9.7–13.8% more
peptide sequences; its true FDP is within nominal (0.79–0.91%) while Comet's is above it on every run
(1.14–1.31%; 23 of 24 draws above 1%), so at a true FDP of 1% the gap is 13.7–21.1%. Draw-to-draw
SD is 0.08–0.21 pp for the FDP and about 1% of the PSM count. On the same mzML andes still reports
8.1–15.5% more PSMs; native `.d` adds 1.2–2.6%. andes is 2.9–3.6x faster and uses 2.2–2.6x less
memory. The tools build their decoys differently (andes permuted, Comet reversed), so part of
Comet's excess over nominal may come from the decoys; the depth at a true FDP of 1% does not depend
on them. Comet's settings were not tuned for PASEF data. Wall times exclude conversion (about 5–7
min per run for tdf2mzml) and Percolator.

**Glyco.** Mouse liver glycoproteome data (PXD005553), against mouse proteins plus an entrapment,
with `--glyco --glyco-species mouse --decoy-strategy sequon`. Two published reference
identification sets for the same spectra ship in [`glyco/truth/`](glyco/truth/): pGlyco2 (the
depositors') and MSFragger-Glyco (a later reanalysis, PXD031032). "Confirmed" is the share of a
reference's spectra that andes identifies with the same peptide.
- *Quick tier:* the `MouseLiver-Z-T-1` fraction as converted mzML, three peptide-level shuffled
  databases × Percolator seeds 42, 1 and 2: 6,992–7,072 glycoPSMs at 1%, true FDP 0.61–1.35%,
  pGlyco2 88.2–88.4% and MSFragger-Glyco 89.3–89.7% confirmed; 5.8–6.6 min and 5.4–5.6 GB on 8
  threads. Each run has only 20–45 entrapment hits, so the spread is mostly counting noise.
  Reading the `.raw` natively, as the commands below do, can give somewhat different counts.
- *Deep tier:* all five fractions, pooled before Percolator, five seeds. Measured 2026-09-05 at
  `14818d3e` on an older, smaller glycan list and converted mzML; not re-measured since.
- The open problem is selection, not generation: almost every reference spectrum produces a
  candidate row, but in the deep tier about a fifth of pGlyco2's spectra lose the per-scan
  choice.

### Quantification: label-free, ProteoBench

**Label-free quantification with match between runs (ProteoBench).**
- *Data and measurement:* the two ProteoBench "Quant LFQ DDA ion" modules, six runs each, where
  the human, yeast and E. coli ratios between conditions A and B are known: QExactive
  (`LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}`, PXD028735, Q Exactive HF-X) and
  Astral (`LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}`, PXD070049). Measured
  2026-11-05 at `09e1c00d` on a cluster node (8 threads per run), andes defaults throughout
  (§3, `--predict rt`), Percolator 3.7.1 through `--rescore` (seed 42). Scored with the official
  `proteobench` package 0.18.6 ("Custom" format, converted by
  [`proteobench/andes_to_proteobench.py`](proteobench/)).
- *Which quantity is scored:* `quantity_normalized` of the MBR table, which removes each run's
  loading level and its drift along the gradient, unless the table says raw. The raw `quantity`
  is what andes hands to quantms.
- *Results* at a minimum of 3 runs; precursors / median |ε| / E. coli |ε| / CV:

  | module | transfer q ≤ 0.01, normalised | raw quantity | transfer q ≤ 0.05, normalised | identified only, normalised |
  |---|---|---|---|---|
  | QExactive | **63,796 / 0.161 / 0.263 / 0.117** | 0.222 / 0.309 / 0.245 | 71,299 / 0.175 / 0.441 / 0.132 | 52,166 / 0.154 / 0.203 / 0.105 |
  | Astral | **36,675 / 0.186 / 0.252 / 0.112** | 0.204 / 0.254 / 0.113 | 43,286 / 0.205 / 0.349 / 0.127 | 31,739 / 0.186 / 0.235 / 0.100 |

  At least 1 / 6 runs (transfer q 0.01): QExactive 79,492 / 41,370, Astral 49,012 / 24,393.
  The last column is the same table with the transfers left out (match between runs off,
  `--evidence identified`); it is scored with the module's MBR flag off. Before this refresh
  the same measurement at `2051616b` gave 65,012 / 0.163 / 0.275 / 0.118 on QExactive and
  37,352 / 0.193 / 0.258 / 0.113 on Astral.
  The opt-in `--predict` settings `rt,ms2` and `ms2,peptdeep-rt,deeplc-rt` (AlphaPeptDeep and
  DeepLC models, `--predictor-models`, ONNX Runtime backend) were measured on 2026-10-21 at
  `a27a72c`, before the extraction change of #240: on QExactive 67,906 and 68,077 precursors
  against 65,391 for the default on that build, on Astral 34,413 and 33,987 against 33,876.
- *By abundance* (transfer q 0.01, normalised; quartiles of each module's own mean intensity,
  |ε| / E. coli |ε| / CV): QExactive Q1 0.279 / 0.461 / 0.195, Q2 0.169 / 0.263 / 0.130,
  Q3 0.136 / 0.205 / 0.105, Q4 0.122 / 0.183 / 0.077; Astral Q1 0.319 / 0.406 / 0.194, Q2 0.208 /
  0.303 / 0.138, Q3 0.161 / 0.212 / 0.103, Q4 0.122 / 0.197 / 0.060.
- *Cost* (six runs, CPU seconds, one run each on cluster nodes): QExactive search with prediction
  and Percolator 5.8k, library 4, extract 0.2k, finalize 0.2k per transfer threshold; Astral 5.8k,
  4, 30 and 0.1k. Each search takes 3–4 minutes of wall time, each extract about 6 seconds
  (QExactive) or 2 seconds (Astral), finalize 34 and 19 seconds.
- *Error control* is measured on the QExactive runs with the HYE database plus a peptide-level
  shuffled entrapment (T/E 1.075), andes defaults, Percolator 3.7.1 `-Y` seeds 42, 1 and 2
  outside andes (re-measured 2026-11-05 at `09e1c00d`). PSM FDP at q 0.01 is 0.84–0.85%;
  precursors quantified in at least 3 runs at transfer q 0.01 are 0.20–0.21% entrapment-only
  (0.26–0.27% at 0.05). Of the accepted transfers, the entrapment measures a false-donor share
  of 0.22–0.26% at transfer q 0.01 and 0.37–0.40% at 0.05; the identity-decoy estimate is
  0.28–0.29% and 0.42–0.47%. The estimate is above the measured share at both thresholds in all
  three seeds here (at `2051616b` it was below in one of six cases), so it tracks the
  false-donor share; three seeds do not make it an upper bound. The wrong-peak part has no
  entrapment check.

### Quantification: ProteoBench leaderboard comparison

**Comparison with the ProteoBench leaderboard.** Public results of both modules
(`Results_quant_ion_DDA` and `Results_quant_ion_DDA_Astral`, both last updated 2026-07-16; read
2026-10-21; scored by `proteobench` 0.17.0). For each tool the submission with the most
precursors and, when different, the one with the lowest median |ε|, at a minimum of 3 runs.
The last column says what limits the error of transferred (match-between-runs) values, taken
from the submission's parameters, comment and the tool's own output; "not checked" means the
submission's parameter file could not be retrieved (the ProteoBench data server was unreachable
on 2026-10-21).

*QExactive:*

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | error control on transferred values |
|---|---|---:|---:|---:|---:|---|
| andes, default (2026-11-05) | yes | 63,796 | 0.161 | 0.263 | 0.117 | PSM q ≤ 0.01; transfer q ≤ 0.01 (decoy-estimated wrong peaks and wrong donors) |
| andes, default, transfer q ≤ 0.05 | yes | 71,299 | 0.175 | 0.441 | 0.132 | transfer q ≤ 0.05 |
| andes, default, raw quantity | yes | 63,796 | 0.222 | 0.309 | 0.245 | as the default |
| andes, default, identified rows only | no | 52,166 | 0.154 | 0.203 | 0.105 | no transfers |
| FragPipe 23.0, most precursors | yes | 63,701 | 0.200 | 0.306 | 0.156 | identifications 1%; transfers at IonQuant MBR ion FDR 0.01 (submitted workflow) |
| FragPipe 23.0, lowest \|ε\| | yes | 54,374 | 0.190 | 0.257 | 0.154 | as above |
| MaxQuant 1.5.3.30, most precursors | yes | 44,299 | 0.210 | 0.386 | 0.219 | PSM and protein 1%; transfers not checked |
| MaxQuant 2.4.13.0, lowest \|ε\| | no | 29,737 | 0.192 | 0.263 | 0.199 | no transfers |
| i2MassChroQ 1.0.16, most precursors | yes | 74,697 | 0.222 | 0.460 | 0.132 | identifications about 1%; no transfer error setting in the submitted parameters |
| i2MassChroQ 1.2.0, lowest \|ε\| | no | 49,433 | 0.198 | 0.263 | 0.120 | no transfers |
| AlphaPept 0.5.0, most precursors | yes | 37,760 | 0.219 | 0.281 | 0.165 | PSM and protein 1%; transfers not checked |
| AlphaPept 0.5.3, lowest \|ε\| | no | 36,176 | 0.214 | 0.255 | 0.165 | no transfers |
| Sage 0.15.0-beta.2 | yes | 81,933 | 0.388 | 0.877 | 0.393 | peptide q ≤ 0.01 on identifications only; every charge state 2–4 of those peptides is integrated in every run and its peak q-value is not applied |
| MSAngel 2.2.10 | yes | 60,206 | 0.216 | 0.514 | 0.227 | PSM 1%; transfers not checked |
| ProlineStudio 2.3.0-SNAPSHOT | yes | 55,639 | 0.210 | 0.508 | 0.226 | PSM 1%; transfers not checked |
| WOMBAT 0.9.11, most precursors | yes | 55,943 | 0.211 | 0.692 | 0.235 | identifications 1%; transfers not checked |
| WOMBAT 0.9.11, lowest \|ε\| | yes | 41,016 | 0.195 | 0.394 | 0.213 | as above |
| PEAKS 13 | no | 72,741 | 0.289 | 0.741 | 0.216 | no transfers |
| quantms v1.3.1dev (OpenMS ProteomicsLFQ) | no | 33,040 | 0.209 | 0.387 | 0.224 | no transfers |
| quantms (OpenMS ProteomicsLFQ), our run (2026-10-17, its default) | no | 34,677 | 0.210 | 0.401 | 0.225 | no transfers |
| quantms (OpenMS ProteomicsLFQ), our run, its feature-linking MBR | yes | 39,257 | 0.221 | 0.563 | 0.234 | feature-quality score; no transfer error rate |

andes's rows use its own normalisation (`quantity_normalized`: median |ε| 0.161, CV 0.117). Its
raw quantity reads |ε| 0.222 and CV 0.245, less accurate than most submissions; normalisation
differs between tools, so the |ε| column does not compare extraction alone. Among the
submissions, FragPipe 23.0 reaches median |ε| 0.190 at 54,374 precursors, and i2MassChroQ 1.2.0
a CV of 0.120 without match between runs. At least six submissions report more precursors than
andes (Sage, the largest i2MassChroQ submissions and PEAKS 13 among them), at a median |ε| of
0.214–0.388. The two quantms runs marked "our run" are ours (2026-10-17, its development
version, raw intensities).

*On the precursors andes shares* (default at `09e1c00d`, transfer q 0.01; the submission's values | andes raw
| andes normalised, median |ε| / E. coli |ε| / CV):

| submission | shared | submission | andes raw | andes normalised | only in the submission |
|---|---:|---|---|---|---|
| FragPipe 23.0, most precursors | 51,722 | 0.189 / 0.274 / 0.155 | 0.207 / 0.276 / 0.241 | 0.147 / 0.225 / 0.108 | 11,979 at 0.272 / 0.708 |
| FragPipe 23.0, lowest \|ε\| | 48,938 | 0.187 / 0.243 / 0.155 | 0.203 / 0.268 / 0.241 | 0.144 / 0.207 / 0.108 | 5,436 at 0.220 / 0.510 |
| i2MassChroQ 1.0.16, most precursors | 50,910 | 0.207 / 0.348 / 0.125 | 0.209 / 0.291 / 0.243 | 0.152 / 0.241 / 0.112 | 23,787 at 0.269 / 0.948 |
| i2MassChroQ 1.2.4, 57,787 precursors | 49,853 | 0.202 / 0.348 / 0.123 | 0.209 / 0.291 / 0.242 | 0.151 / 0.242 / 0.112 | 7,934 at 0.279 / 0.835 |
| Sage 0.15.0-beta.2, our reproduction | 54,633 | 0.309 / 0.613 / 0.319 | 0.215 / 0.291 / 0.245 | 0.156 / 0.246 / 0.114 | 25,724 at 0.683 / 1.364 |

Sage's per-submission table was not retrievable, so its row uses a reproduction (2026-10-21): the
same Sage version on the same six mzML files with the submission's stated parameters and
per-charge output gives 80,357 precursors at |ε| 0.388, E. coli 0.878 and CV 0.394 (submission:
81,933, 0.388, 0.877, 0.393). Most of its precursors that andes does not report are charge
states that Sage never identified or are above its own peak q-value of 0.01, which it reports but
does not apply.

*Astral:*

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | error control on transferred values |
|---|---|---:|---:|---:|---:|---|
| andes, default (2026-11-05) | yes | 36,675 | 0.186 | 0.252 | 0.112 | PSM q ≤ 0.01; transfer q ≤ 0.01 (decoy-estimated) |
| andes, default, transfer q ≤ 0.05 | yes | 43,286 | 0.205 | 0.349 | 0.127 | transfer q ≤ 0.05 |
| andes, default, raw quantity | yes | 36,675 | 0.204 | 0.254 | 0.113 | as the default |
| andes, default, identified rows only | no | 31,739 | 0.186 | 0.235 | 0.100 | no transfers |
| FragPipe 23.0, most precursors | yes | 39,057 | 0.212 | 0.250 | 0.127 | identifications 1%; transfers not checked |
| FragPipe 22.1-build02, lowest \|ε\| | yes | 29,105 | 0.181 | 0.405 | 0.078 | as above |
| MaxQuant 2.6.7.0 | yes | 20,283 | 0.160 | 0.158 | 0.092 | PSM and protein 1%; transfers within a condition only (submission comment), not checked |
| i2MassChroQ 1.2.6, most precursors | yes | 36,469 | 0.197 | 0.417 | 0.095 | identifications about 1%; transfers not checked |
| i2MassChroQ 1.2.9, lowest \|ε\| | no | 24,803 | 0.178 | 0.245 | 0.077 | no transfers |
| Sage 0.14.6 | yes | 30,319 | 0.167 | 0.338 | 0.070 | peptide q ≤ 0.01 on identifications only; peak q-value not applied (as in QExactive) |
| DIA-NN 2.6.0, most precursors | yes | 93,237 | 0.277 | 0.691 | 0.159 | run and global precursor q ≤ 0.01 (submission comment) |
| DIA-NN 2.6.0, lowest \|ε\| | yes | 16,363 | 0.123 | 0.168 | 0.079 | as above, plus MS1-level q ≤ 0.01 |

On Astral, MaxQuant 2.6.7.0 (20,283 precursors, median |ε| 0.160) and Sage 0.14.6 (30,319, 0.167,
CV 0.070) report fewer precursors than andes (36,675, 0.186), i2MassChroQ 1.2.6 (36,469, 0.197,
CV 0.095) about as many, and FragPipe 23.0 more (39,057, 0.212); DIA-NN 2.6.0 reports 93,237 at
0.277, and its 16,363-precursor submission has the lowest median |ε| (0.123). andes's CV (0.112)
is above most submissions. Normalisation changes little here (raw |ε| 0.204, CV 0.113).
Shared-precursor comparisons need the per-submission tables, which could not be retrieved; they
are not given.

### Quantification: truth designs, match between runs on and off

**Quantification against truth designs, match between runs on and off.** Eight sets with a known
answer, quantified by andes, MaxQuant and FragPipe and read through one evaluation
([`quant-eval`](quant-eval/README.md)). andes was run on 2026-11-05 at `09e1c00d`, as were the
MaxQuant runs with match between runs off; the MaxQuant runs with it on and the FragPipe outputs
are from 2026-10 and were read again unchanged.
- *Sets and truth:*

  | set | what it is | runs | truth used |
  |---|---|---:|---|
  | PXD009740 | UPS1 standard at ten levels (170 amol to 250 fmol) in a constant yeast background | 30 | UPS1 amount per run; yeast constant |
  | PXD009815 | UPS1 standard at ten levels (10 amol to 50 fmol) in yeast | 40 | as above |
  | PXD001819 | UPS1 standard at nine levels in yeast; level and run order are confounded | 27 | as above (counts and accuracy; its differential expression is not a decision set, every tool reports a thousand or more yeast false positives) |
  | QE (ProteoBench) | human, yeast and *E. coli* at known ratios, Q Exactive HF-X (PXD028735) | 6 | yeast 2:1, *E. coli* 1:4, human 1:1 between conditions A and B |
  | Astral (ProteoBench) | the same mixture, Astral (PXD070049) | 6 | as above |
  | PXD030780 | human, yeast and *E. coli* at five mixing levels A–E | 25 | species amounts per group |
  | PXD014415 | human and human + yeast runs | 40 | yeast absent from 20 runs (false yeast rows) |
  | PXD016921 | a blank run among cell-input runs | 6 | no analyte in the blank run |

  Further truth sets are held back for later evaluations and were not run for this page.
- *andes:* origin/main `09e1c00d` with its defaults (trace-area share 0.1, the shared-feature
  downgrade, the quality columns); `quant library` (`--mbr-library-fdr 0.01`), `quant extract` per
  run and `quant finalize` (`--mbr-fdr 0.01`) over the search-time tables of the 2026-10-31
  searches (`--predict rt,ms2 --chimeric --refine --max-missed-cleavages 2`, 8 threads, one
  Percolator 3.7.1 `-Y` run per set, `--quant-dir`); the identifications were not searched again at
  `09e1c00d`. Quantity: the raw `quantity` column. *MBR on* is the default view of the finalize
  table (identified rows and transfers with `transfer_q ≤ 0.01`); *MBR off* is the identified rows
  of the same table. The finalize ran with the quant-model seeds 42, 1 and 7 (25 runs of
  PXD030780: 42 only); the tables show seed 42 and the text gives the seed range where it
  matters.
- *MaxQuant:* 2.7.0.0, Standard DDA with LFQ on, the set's own FASTA without contaminants, one
  experiment per raw file, 1% PSM and protein FDR, 2 missed cleavages, Carbamidomethyl C fixed,
  Oxidation M and protein N-terminal acetylation variable, run on the vendor files. The MBR-on
  runs have match between runs on (with its FDR option for PXD009740, PXD009815, QE and Astral);
  the MBR-off run of each set is the same parameter file with only `matchBetweenRuns` (and the MBR
  FDR option) switched off, re-run on the same files. The MBR-on runs of PXD001819 and of the blank
  run had the MBR FDR option off. Re-running them with the option on changed no number: MaxQuant
  2.7 does not remove match rows with that option, it only fills the `Match q-value` column of the
  evidence table. The MaxQuant rows here are therefore unfiltered transfers on every set, while
  andes' transfers are filtered at `transfer_q ≤ 0.01`. Filtering MaxQuant's matches at
  `Match q-value ≤ 0.01` cuts PXD001819 transfers from 43,373 to 2,550 and blank-run transfers
  from 85 to 9; a filtered MaxQuant arm is being added.
  No PXD014415 and PXD030780 MaxQuant run exists (PXD014415 has no MaxQuant column in the absence
  results below). The
  quantity is `Intensity` of the evidence table summed per precursor and run; the LFQ columns are
  not used.
- *FragPipe:* the 21.1 workflow (MSFragger, MSBooster, Percolator, Philosopher) with IonQuant
  1.10.12 run with match between runs (ion FDR 0.01, MaxLFQ) and **normalisation off**
  (`--normalization 0`), so that its ion intensities are raw like the others; the existing runs
  were used, none with match between runs off.
- *Search settings are not matched.* andes searched with chimeric spectra (`--chimeric`), the secondary-modification pass (`--refine`) and
  with PSM-level FDR, while MaxQuant and FragPipe applied a protein-level 1% FDR, so depth counts
  (precursors, proteins, rows) are not settings-matched. With settings matched to the tools
  (±10 ppm, isotope errors 0 to 2, no `--chimeric`, no `--refine`; seed 42) andes' precursors fall
  by 31% on PXD009740 and PXD009815 (to 27,226 and 29,600, level with MaxQuant and FragPipe), 12%
  on QE (81,987) and 6% on Astral (51,798). Its DE true positives do not change beyond seed spread
  (QE 1,049, Astral 673, PXD009815 831, PXD009740 1,293). So most of the UPS-set depth lead comes
  from the search settings, while the DE results come from the quantification. A protein-level 1%
  filter on andes removes 0.4-1.6% of precursors.
- *Evaluation rules.* Every tool's table is reduced to raw precursor-by-run intensities and goes
  through the same downstream steps: per-run median normalisation over the precursors quantified
  in all runs, Tukey median polish per protein over precursors and runs, proteins with at least two
  precursors tested with a moderated t, Benjamini–Hochberg 5% per condition pair. A true positive
  is a significant protein of a species whose amount changes, with the expected sign; a false
  call is a significant protein of the species whose amount is constant (yeast in the UPS1
  sets, human in the mixtures) or a significant protein of a changing species with the opposite
  sign (wrong-sign calls are counted as false and given in brackets); counts are given without and
  with a fold-change cut of |log2 FC| ≥ 0.5. Differential expression is scored on the
  common-tested family: the (protein, condition pair) rows that andes, MaxQuant and FragPipe
  (MBR on; PXD030780: andes and FragPipe) all test, with the moderated t, its variance prior and
  the Benjamini–Hochberg correction run on that set, so every arm is counted on identical rows;
  an arm with match between runs off is scored on the part of the family it can test. (An earlier
  version of this page tested each arm on its own proteins, which gave andes more rows than the
  others, 4,699 against MaxQuant's 4,312 on QE, and counted wrong-sign calls as true positives;
  both are corrected here.) The columns for fold-change cut and for true positives at pooled FDP
  10% are on the previous counting and were not recomputed on the common family. andes
  normalisation, protein roll-up and weighting are not part of the comparison: the same code
  treats all tools. Identified-cell accuracy is scored on the cells (precursor in run) that every
  tool in the table identifies, against the replicate-based expectation of the design; the
  bright decile is the top tenth of the predicted intensity and its gross error the share of
  cells off by more than 2x. UPS1 detection counts a precursor at a level when one of its values
  is within 2x of its own high-level reference. Absence is the yeast share of the transferred
  rows in the 20 human-only runs of PXD014415 (false-transfer rate) and the rows reported in the
  blank run. The masked wrong-peak rate removes 5% of the identified precursor-run pairs from the
  library, quantifies them again as transfers and counts those whose apex lies further from the
  original identified apex than that fit's peak width, at `quant_q ≤ 0.01` (andes only; no other
  tool exposes the quantity).
- *Depth* (precursors counted once if quantified in any run; cells are precursor-by-run values;
  PXD014415 and the blank run are absence designs):

  | set | tool | precursors | proteins (≥ 2 precursors) | identified cells | transferred cells |
  |---|---|---:|---:|---:|---:|
  | PXD009740 | andes, MBR on | 39,595 | 3,031 | 573,349 | 332,626 |
  |  | andes, MBR off | 39,550 | 3,027 | 573,349 | 0 |
  |  | MaxQuant, MBR on | 29,293 | 2,339 | 464,069 | 175,450 |
  |  | MaxQuant, MBR off | 29,292 | 2,339 | 464,069 | 0 |
  |  | FragPipe, MBR on | 25,888 | 2,259 | 414,258 | 141,367 |
  | PXD009815 | andes, MBR on | 42,686 | 3,362 | 798,525 | 463,945 |
  |  | andes, MBR off | 42,672 | 3,362 | 798,525 | 0 |
  |  | MaxQuant, MBR on | 27,459 | 2,453 | 586,274 | 230,204 |
  |  | MaxQuant, MBR off | 27,459 | 2,453 | 586,274 | 0 |
  |  | FragPipe, MBR on | 28,286 | 2,691 | 587,049 | 193,057 |
  | PXD001819 | andes, MBR on | 10,696 | 1,004 | 146,103 | 65,796 |
  |  | andes, MBR off | 10,696 | 1,004 | 146,103 | 0 |
  |  | MaxQuant, MBR on | 9,314 | 881 | 116,299 | 43,373 |
  |  | MaxQuant, MBR off | 9,314 | 881 | 116,299 | 0 |
  |  | FragPipe, MBR on | 10,406 | 937 | 118,961 | 63,595 |
  | QE (ProteoBench) | andes, MBR on | 92,828 | 6,797 | 350,579 | 71,126 |
  |  | andes, MBR off | 92,827 | 6,797 | 350,579 | 0 |
  |  | MaxQuant, MBR on | 51,609 | 4,810 | 186,178 | 25,465 |
  |  | MaxQuant, MBR off | 51,608 | 4,810 | 186,178 | 0 |
  |  | FragPipe, MBR on | 71,377 | 5,938 | 249,172 | 12,401 |
  | Astral (ProteoBench) | andes, MBR on | 55,246 | 4,901 | 212,633 | 28,812 |
  |  | andes, MBR off | 55,246 | 4,901 | 212,633 | 0 |
  |  | MaxQuant, MBR on | 24,745 | 2,855 | 93,669 | 29,072 |
  |  | MaxQuant, MBR off | 24,745 | 2,855 | 93,669 | 0 |
  |  | FragPipe, MBR on | 46,347 | 4,334 | 168,481 | 46,343 |
  | PXD030780 | andes, MBR on | 93,566 | 8,602 | 1,047,714 | 560,171 |
  |  | andes, MBR off | 93,563 | 8,601 | 1,047,714 | 0 |
  |  | FragPipe, MBR on | 65,586 | 6,885 | 824,966 | 351,518 |
  | PXD014415 | andes, MBR on | 76,099 | 5,470 | 1,494,110 | 801,869 |
  |  | andes, MBR off | 76,047 | 5,468 | 1,494,110 | 0 |
  |  | FragPipe, MBR on | 61,771 | 5,030 | 1,175,139 | 550,247 |
  | blank run (PXD016921) | andes, MBR on | 17,625 | 1,728 | 29,755 | 1,641 |
  |  | andes, MBR off | 17,625 | 1,728 | 29,755 | 0 |
  |  | MaxQuant, MBR on | 11,010 | 1,229 | 17,749 | 1,305 |
  |  | MaxQuant, MBR off | 11,010 | 1,229 | 17,749 | 0 |
  |  | FragPipe, MBR on | 14,153 | 1,403 | 24,638 | 4,015 |
  🟢 **Depth: well ahead** on 7 of 8 sets (+19% to +51% precursors over the best other tool);
  comparable on PXD001819 (+3%).
- *Accuracy and precision* (cells identified by every tool listed in the set; the andes cells
  are identified cells, so the andes MBR-off accuracy equals the MBR-on one and is not repeated;
  protein CV is the median replicate CV of matched features):

  | set | tool | \|err\| > 1 (%) | within 2x (%) | bright-decile gross (%) | precursor SD, top decile | protein CV (%) |
  |---|---|---:|---:|---:|---:|---:|
  | PXD009740 | andes, MBR on | 2.63 | 97.37 | 3.64 | 0.132 | 8.7 |
  |  | andes, MBR off | - | - | - | - | 10.2 |
  |  | MaxQuant, MBR on | 1.55 | 98.45 | 1.13 | 0.116 | 10.5 |
  |  | MaxQuant, MBR off | 1.42 | 98.58 | 1.14 | 0.116 | 10.5 |
  |  | FragPipe, MBR on | 4.02 | 95.98 | 2.37 | 0.126 | 13.5 |
  | PXD009815 | andes, MBR on | 1.42 | 98.58 | 5.97 | 0.056 | 4.4 |
  |  | andes, MBR off | - | - | - | - | 5.6 |
  |  | MaxQuant, MBR on | 0.35 | 99.65 | 0.55 | 0.066 | 7.0 |
  |  | MaxQuant, MBR off | 0.40 | 99.60 | 0.89 | 0.066 | 7.5 |
  |  | FragPipe, MBR on | 0.65 | 99.35 | 0.68 | 0.051 | 6.1 |
  | PXD001819 | andes, MBR on | 1.81 | 98.19 | 1.79 | 0.044 | 5.7 |
  |  | andes, MBR off | - | - | - | - | 6.8 |
  |  | MaxQuant, MBR on | 0.86 | 99.14 | 0.46 | 0.060 | 8.3 |
  |  | MaxQuant, MBR off | 0.89 | 99.11 | 0.74 | 0.060 | 7.9 |
  |  | FragPipe, MBR on | 1.80 | 98.20 | 3.76 | 0.038 | 5.2 |
  | QE (ProteoBench) | andes, MBR on | 3.75 | 96.25 | 6.44 | 0.208 | 9.2 |
  |  | andes, MBR off | - | - | - | - | 9.6 |
  |  | MaxQuant, MBR on | 2.92 | 97.08 | 1.62 | 0.184 | 11.6 |
  |  | MaxQuant, MBR off | 2.38 | 97.62 | 1.65 | 0.181 | 10.3 |
  |  | FragPipe, MBR on | 5.84 | 94.16 | 2.64 | 0.206 | 11.5 |
  | Astral (ProteoBench) | andes, MBR on | 1.62 | 98.38 | 4.39 | 0.057 | 5.8 |
  |  | andes, MBR off | - | - | - | - | 6.0 |
  |  | MaxQuant, MBR on | 0.69 | 99.31 | 0.26 | 0.072 | 7.0 |
  |  | MaxQuant, MBR off | 0.64 | 99.36 | 0.44 | 0.072 | 6.9 |
  |  | FragPipe, MBR on | 1.13 | 98.87 | 0.66 | 0.041 | 4.9 |
  | PXD030780 | andes, MBR on | 2.63 | 97.37 | 7.19 | 0.112 | 10.1 |
  |  | andes, MBR off | - | - | - | - | 12.1 |
  |  | FragPipe, MBR on | 1.44 | 98.56 | 0.30 | 0.096 | 10.2 |
  🔴 **Identified-cell accuracy: behind to well behind** on every set (bright-decile gross error
  1.2–7.2% against MaxQuant 0.3–1.6%). 🟢 **Protein CV: ahead** on PXD009740, PXD009815 and QE,
  comparable on PXD030780, behind FragPipe on Astral and PXD001819.
- *Differential expression* (BH 5% per pair; TP / FP is the common-tested family with the
  expected sign required, wrong-sign calls counted in FP and given in brackets; the seed range of
  andes with MBR on over seeds 42, 1 and 7 is 1,288–1,302 / 305–350 TP / FP on PXD009740,
  841–847 / 153–188 on PXD009815, 1,047–1,050 / 110–120 on QE and 681–685 / 296–301 on Astral,
  with 0–1, 7–13, 0 and 1 wrong-sign calls):

  | set | tool | TP / FP, common family, sign-aware (wrong sign) | TP / FP, \|log2 FC\| ≥ 0.5 (previous counting) | TP at pooled FDP 10% (previous counting) | ≥ 0.5 | TP at pooled FDP 10% |
  |---|---|---:|---:|---:|
  | PXD009740 | andes, MBR on | 1,300 / 305 (1) | 1,376 / 256 | 1,414 |
  |  | andes, MBR off | 1,068 / 180 (0) | 1,075 / 115 | 1,123 |
  |  | MaxQuant, MBR on | 1,304 / 129 (0) | 1,361 / 101 | 1,442 |
  |  | MaxQuant, MBR off | 1,047 / 188 (0) | 1,054 / 111 | 1,076 |
  |  | FragPipe, MBR on | 1,237 / 302 (0) | 1,262 / 284 | 1,272 |
  | PXD009815 | andes, MBR on | 844 / 153 (7) | 997 / 25 | 1,010 |
  |  | andes, MBR off | 640 / 76 (0) | 644 / 19 | 657 |
  |  | MaxQuant, MBR on | 844 / 92 (2) | 1,010 / 46 | 1,030 |
  |  | MaxQuant, MBR off | 541 / 54 (2) | 562 / 29 | 575 |
  |  | FragPipe, MBR on | 849 / 301 (2) | 937 / 107 | 915 |
  | PXD001819 | andes, MBR on | 375 / 3,296 (0) | 479 / 477 | 419 |
  |  | andes, MBR off | 266 / 1,833 (0) | 272 / 449 | 220 |
  |  | MaxQuant, MBR on | 358 / 939 (0) | 441 / 419 | 361 |
  |  | MaxQuant, MBR off | 240 / 861 (0) | 242 / 314 | 201 |
  |  | FragPipe, MBR on | 359 / 1,977 (0) | 441 / 247 | 346 |
  | QE (ProteoBench) | andes, MBR on | 1,047 / 110 (0) | 1,308 / 16 | 1,331 |
  |  | andes, MBR off | 1,007 / 62 (0) | 1,113 / 19 | 1,155 |
  |  | MaxQuant, MBR on | 857 / 37 (0) | 900 / 19 | 987 |
  |  | MaxQuant, MBR off | 824 / 43 (0) | 821 / 16 | 871 |
  |  | FragPipe, MBR on | 946 / 43 (1) | 984 / 23 | 1,016 |
  | Astral (ProteoBench) | andes, MBR on | 685 / 296 (1) | 752 / 26 | 705 |
  |  | andes, MBR off | 640 / 251 (0) | 669 / 31 | 618 |
  |  | MaxQuant, MBR on | 650 / 210 (2) | 643 / 39 | 608 |
  |  | MaxQuant, MBR off | 465 / 141 (0) | 463 / 33 | 439 |
  |  | FragPipe, MBR on | 689 / 389 (1) | 749 / 17 | 718 |
  | PXD030780 | andes, MBR on | 17,692 / 1,374 (23) | 13,122 / 167 | 20,063 |
  |  | andes, MBR off | 14,416 / 915 (29) | 10,263 / 212 | 16,133 |
  |  | FragPipe, MBR on | 18,227 / 1,751 (22) | 13,096 / 205 | 19,110 |
  🟢 **True positives: ahead** on QE (1,047 against MaxQuant 857, +22%, and FragPipe 946, +11%;
  seeds 1,047–1,050), with more false calls than both (110 against 37 and 43). **Tie** on PXD009740
  (1,300 against 1,304 and 1,237), PXD009815 (844 against 844 and 849), Astral (685 against 650
  and 689) and PXD030780 (17,692 against FragPipe 18,227), with more false calls than MaxQuant on the
  first three (305, 153 and 296 against 129, 92 and 210). 🔴 **False calls without a fold-change cut:
  well behind** MaxQuant on PXD009740 and QE (2.4x and 3.0x), behind on PXD009815 and Astral
  (1.7x and 1.4x); FragPipe has as many on PXD009740 (302), more on PXD009815 (301) and Astral
  (389), fewer on QE (43). Wrong-sign calls are few (0–13 per arm and set; andes 7 on PXD009815
  against 2 for MaxQuant and FragPipe). With |log2 FC| ≥ 0.5 (seed 42) false calls are 11 against
  MaxQuant 18 on QE, 39 against 36 on PXD009815 and 167 against 74 on PXD009740 (Astral not
  recomputed).
- *Why andes makes more false calls* (seed 42, PXD009740, PXD009815 and QE; analysis of the tables
  above, no re-run). The protein values of the unchanged species are tighter in andes than in the
  other tools: over every unchanged protein and pair, the within-condition spread of MaxQuant is
  1.15–1.55x andes' and that of FragPipe 1.25–1.40x. A small shift that the other tools blur
  becomes significant in andes. The calls are not made by transferred values: recomputed from
  identified cells only, the called proteins keep 80–95% of their fold change, and in 78–82% of
  them (PXD009740, PXD009815) most precursors of the protein move together across the replicates.
  On QE 60% of the calls follow a retention-time-local intensity shift that identified precursors of
  other unchanged proteins show as well. MaxQuant and FragPipe see the same direction on most called
  proteins, with 0.56–0.89 of andes' fold change and 1.4–2.4x its spread, so they rarely reach
  significance. On PXD009740, 296 of andes' 320 false calls are comparisons against the highest
  UPS1 level (250,844 amol), where the unchanged yeast shifts in every tool's data; without that
  level andes has 955 true positives and 24 false calls against MaxQuant's 953 and 31 (filtered
  MaxQuant 780 and 13, FragPipe 896 and 61). Calls significant and above a minimal fold change:

  | set | tool | all | \|log2 FC\| ≥ 0.3 | ≥ 0.5 | ≥ 1.0 |
  |---|---|---|---|---|---|
  | QE | andes | 1,049 / 113 | 1,049 / 83 | 1,046 / 11 | 476 / 0 |
  |  | MaxQuant | 857 / 37 | 857 / 35 | 852 / 18 | 303 / 2 |
  |  | FragPipe | 946 / 43 | 946 / 37 | 945 / 20 | 680 / 8 |
  | PXD009740 | andes | 1,287 / 320 | 1,287 / 310 | 1,287 / 167 | 1,271 / 27 |
  |  | MaxQuant | 1,292 / 129 | 1,292 / 128 | 1,292 / 74 | 1,238 / 25 |
  |  | FragPipe | 1,224 / 302 | 1,224 / 302 | 1,224 / 259 | 1,220 / 86 |
  | PXD009815 | andes | 853 / 183 | 852 / 60 | 846 / 39 | 800 / 21 |
  |  | MaxQuant | 842 / 92 | 842 / 66 | 838 / 36 | 778 / 6 |
  |  | FragPipe | 859 / 301 | 859 / 184 | 858 / 102 | 844 / 22 |

  Entries are true positives / false calls. The truth labels assume the unchanged species does not
  change at all; at andes' precision that assumption no longer holds for shifts below about 1.4-fold (|log2 FC| 0.5).
  Still behind: PXD009815 at large fold changes (21 against MaxQuant 6 at ≥ 1.0) and its wrong-sign
  calls (18 against 2).
- *UPS1 detection at the lowest levels* (UPS1 precursors within 2x of their own reference, summed
  over the three lowest levels; PXD009740 170, 382 and 859 amol, PXD009815 10, 50 and 100 amol,
  PXD001819 50, 125 and 250 amol). Rows are the UPS1 rows at those levels, within 2x is the share of
  rows, the median is the pooled median |log2 error| against the design reference, and the last
  column is the expected number of precursors within 2x if the arm reported only as many rows as
  FragPipe, drawn at random (the neutral comparison; the best-first ordering differs per tool):

  | set | tool | precursors within 2x, three lowest levels | rows | rows within 2x (%) | median \|log2 err\| | within 2x at FragPipe's row count, random rows |
  |---|---|---:|---:|---:|---:|---:|
  | PXD009740 | andes, MBR on | 502 | 1,600 | 60.3 | 0.84 | 436.7 |
  |  | andes, MBR off | 206 | 546 | 62.6 | 0.84 | 206.0 |
  |  | MaxQuant, MBR on | 417 | 1,642 | 38.9 | 1.20 | 340.7 |
  |  | MaxQuant, MBR off | 111 | 360 | 48.3 | 1.01 | 111.0 |
  |  | FragPipe, MBR on | 458 | 1,167 | 71.1 | 0.64 | 458.0 |
  | PXD009815 | andes, MBR on | 24 | 565 | 7.6 | 4.75 | 22.1 |
  |  | andes, MBR off | 3 | 153 | 2.0 | 4.27 | 3.0 |
  |  | MaxQuant, MBR on | 23 | 603 | 6.3 | 4.63 | 19.1 |
  |  | MaxQuant, MBR off | 0 | 108 | 0.0 | 4.92 | 0.0 |
  |  | FragPipe, MBR on | 50 | 443 | 15.6 | 3.73 | 50.0 |
  | PXD001819 | andes, MBR on | 0 | 53 | 0.0 | 6.00 | 0.0 |
  |  | andes, MBR off | 0 | 1 | 0.0 | 1.04 | 0.0 |
  |  | MaxQuant, MBR on | 1 | 68 | 1.5 | 4.25 | 0.5 |
  |  | MaxQuant, MBR off | 0 | 0 | - | - | 0.0 |
  |  | FragPipe, MBR on | 3 | 29 | 10.3 | 3.13 | 3.0 |

  Per level, MBR on (rows; precursors within 2x; share of rows within 2x; median |log2 error|):

  | set | level (amol) | andes | MaxQuant | FragPipe |
  |---|---:|---|---|---|
  | PXD009740 | 170 | 295; 66; 33.6%; 1.29 | 323; 38; 14.2%; 1.97 | 189; 74; 52.9%; 0.96 |
  |  | 382 | 467; 153; 59.3%; 0.85 | 532; 127; 35.2%; 1.26 | 345; 133; 67.0%; 0.63 |
  |  | 859 | 838; 283; 70.3%; 0.72 | 787; 252; 51.6%; 0.99 | 633; 251; 78.8%; 0.58 |
  | PXD009815 | 10 | 174; 0; 0.0%; 7.23 | 195; 0; 0.0%; 7.16 | 131; 1; 0.8%; 6.53 |
  |  | 50 | 194; 3; 2.1%; 4.71 | 189; 4; 2.1%; 4.39 | 140; 19; 17.1%; 3.97 |
  |  | 100 | 197; 21; 19.8%; 3.03 | 219; 19; 15.5%; 2.94 | 172; 30; 25.6%; 2.26 |

  As reported, andes is 🟢 **ahead** on PXD009740 (502 against MaxQuant 417 and FragPipe 458) and
  🔴 **behind** FragPipe on PXD009815 (24 against 50), but andes reports 1,600 rows on PXD009740
  against FragPipe's 1,167 and its share of rows within 2x is lower (60.3% against 71.1%). At equal
  rows FragPipe is **ahead** on both sets (PXD009740 458 against andes 437; PXD009815 50 against
  22) and andes is **ahead** of MaxQuant on PXD009740 (437 against 341) and **tied** on PXD009815
  (22 against 19). At 10 and 50 amol every tool over-reads: the median error is +6.5 to +7.4 log2
  at 10 amol and +4.0 to +4.7 at 50 amol (a factor of 2^4 to 2^7).
- *Absence* (PXD014415; yeast rows in the 20 human-only runs). Raw yeast share: andes with MBR on
  reports 451 identified and 769 transferred yeast rows (0.20% of its 389,423 transferred rows in
  those runs; seeds 0.19–0.20%), with `quant_q ≤ 0.01` 134 transferred rows (0.08%, 162,219
  transferred rows); with MBR off 451 identified yeast rows. FragPipe reports 208 identified and 243
  transferred (0.09% of 266,898). Many of these rows are not yeast-specific on their intensity
  profile, so the rows are also classed. Classes are defined from one reference
  table for every tool, andes's table of every row (2,892,797 rows, below-threshold rows included):
  for a precursor, d is the median over the runs that contain yeast minus the median over the
  human-only runs of the log2 intensity, normalised on the human proteins; the class is *enriched*
  when d is above a boundary learned from human control precursors by the number of human-only runs
  in which the precursor is seen (6.67, 5.62, 3.94 and 1.07 log2 for 1, 2–4, 5–9 and 10 or more
  runs), *not specific* when it is not, *only absent* when the precursor is seen in human-only runs
  only, and *no reference* when the reference table does not hold the precursor (a tool's
  precursors that andes does not report cannot be classed). Two readings are given: enriched rows
  only (a lower bound of the species-specific errors) and enriched plus no-reference rows. Transfers
  (enriched / not specific / only absent / no reference): andes 202 / 567 / 0 / 0, FragPipe
  27 / 52 / 0 / 164. Rates: enriched only, andes 0.052% against FragPipe 0.010%; with no-reference
  rows, andes 0.052% against FragPipe 0.072%. On identified rows (andes 38 / 412 / 1 / 0, FragPipe
  8 / 56 / 5 / 139) the rates are 0.0053% against 0.0014% (enriched only) and 0.0055% against
  0.0270% (with no-reference rows). The reading flips with the no-reference rows, which the profile
  cannot judge, so both are reported. *Blank run* (PXD016921): rows reported in the blank run, MBR
  off / on: andes 85 / 150 (144 at `quant_q ≤ 0.01`), MaxQuant 20 / 105 (the MBR-on run had the MBR
  FDR option off), FragPipe 8 with MBR; split into identified and transferred rows, MBR on: andes
  85 / 65, MaxQuant 20 / 85, FragPipe 8 / 0.
  🔴 **False transfers: behind** FragPipe at the default threshold on the raw share (0.20% against
  0.09%) and on enriched rows (0.052% against 0.010%), **ahead** when no-reference rows are counted
  (0.052% against 0.072%); comparable at `quant_q ≤ 0.01` on the raw share (0.08%). **Blank run:
  behind** (150 rows against FragPipe 8 and MaxQuant 105).
- *Masked wrong-peak rate* at `quant_q ≤ 0.01` (andes, pooled over runs; seeds 42 / 1 / 7; in
  brackets the rate of the MS1-only transfers, which are the less certain ones; the transfers
  with MS2 evidence read 0.18–0.65% in seed 42): PXD009740 0.50 / 0.50 / 0.49% (0.82), PXD009815
  0.55 / 0.55 / 0.60% (0.92), PXD001819 0.46 / 0.41 / 0.43% (1.01), QE 0.77 / 0.80 / 0.78%
  (1.30), Astral 0.66 / 0.66 / 0.62% (0.66), PXD014415 0.49 / 0.50 / 0.48% (1.09); the blank-run
  set has 126–130 masked pairs and none was wrong. Of the runs with at least 200 masked pairs, 0
  to 5 per set (of 6 to 40) are above 1%.
- *Matched true error.* On PXD014415 the false-transfer rate is the measured truth, so depth can
  be read at a tool's rate: at FragPipe's 0.09% andes keeps 487,849 transferred rows
  (`transfer_q ≤ 0.003`; seeds 0.88–0.89 × FragPipe's 550,247 rows), and 428,888 rows through
  `quant_q ≤ 0.015`; at its own `quant_q ≤ 0.01` it has 317,272 rows at 0.08%. MaxQuant has no
  PXD014415 run, and no tool reports a masked wrong-peak rate to compare.
- *Reading.* Match between runs raises the number of cells and the differential-expression
  true positives of every tool that offers it (common family, sign-aware): andes 1,068 to 1,300 on
  PXD009740, 1,007 to 1,047 on QE, 640 to 685 on Astral and 14,416 to 17,692 on PXD030780;
  MaxQuant 1,047 to 1,304, 824 to 857 and 465 to 650 on PXD009740, QE and Astral. The false-call
  counts move with them in different ways: andes's rise (180 to 305 on PXD009740, 62 to 110 on QE)
  while MaxQuant's fall on PXD009740 and QE (188 to 129, 43 to 37) and rise on Astral and PXD009815.
  On the UPS1 sets MaxQuant with match between runs has the lower false-call counts at the same
  true positives as andes (1,304 / 129 against 1,300 / 305 on PXD009740; 844 / 92 against 844 / 153
  on PXD009815). On the identified cells MaxQuant has the lower gross-error shares (bright decile
  0.3–1.6% against 1.2–7.2% for andes and 0.3–3.8% for FragPipe, MBR on) and andes the lower
  protein CV (4.4–10.1% against 7.0–11.6% for MaxQuant and 4.9–13.5% for FragPipe). andes
  quantifies more precursors and cells than the other tools on every set (see the search settings
  above: these depth counts are not settings-matched). With match between runs on, the UPS1
  precursors detected within 2x at the three lowest levels are 502 against 417 (MaxQuant) and 458
  (FragPipe) on PXD009740 and 24 against 23 and 50 on PXD009815, which at equal rows becomes 437
  against 341 and 458, and 22 against 19 and 50; with it off, 206 against 111 on PXD009740.
  FragPipe reports the fewest rows in the blank run and a lower false-transfer rate on PXD014415
  than andes on the raw share and on enriched rows (0.09% against 0.20%; 0.010% against 0.052%).

### Quantification: TMT and TMT match between runs

**TMT quantification.**
- *Data and measurement:* PXD017803 (yeast triple knockout, TMTpro 16-plex; two files each of
  hrMS2, SPS-MS3 and RTS-MS3) and PXD007683 (human + yeast TMT11 with known yeast ratios, 11
  SPS-MS3 fractions, ion-trap CID MS2), as mzML, with the kit impurity sheet. andes defaults,
  including the quantification gate (precursor purity ≥ 0.7 and Percolator PEP ≤ 0.05;
  §10 of DOCS.md), `--rescore`, 16 threads; measured 2026-11-05 at `0cfa03cc`. MaxQuant 2.7.0.0
  is the run of 2026-10-18 on the vendor files, with the same FASTA and impurity sheet, reporter
  extraction at MS2 or MS3, 1% PSM and protein FDR; its reporter intensities are used as
  reported (`Reporter intensity corrected`). Both tools' PSM tables go through one downstream
  (per-channel median normalisation over the set's PSMs, then the metrics below).
- *Extraction:* reporter extraction and impurity correction match OpenMS IsobaricAnalyzer
  exactly on the same scans.
- *Knockout interference* (median ratio of knockout to non-knockout channel intensity over the
  PSMs of the three knocked-out proteins, lower is better; hrMS2 only scans whose window starts
  below 126 Th):

  | acquisition | andes PSMs at 1% | PSMs quantified (gate passed) | andes interference | without the gate | purity ≥ 0.9 | MaxQuant 2.7.0.0: PSMs, interference |
  |---|---:|---:|---:|---:|---:|---|
  | SPS-MS3 | 9,806 | 8,635 | 0.062 | 0.117 | 0.023 | 6,883; 0.054 |
  | RTS-MS3 | 14,270 | 10,518 | 0.027 | 0.033 | 0.015 | 8,930; 0.038 |
  | hrMS2 | 13,474 | 3,481 | 0.192 | 0.202 | 0.135 | 9,677; 0.156 |

  Few knockout PSMs stand behind each value (34–91 for andes, 26–65 for MaxQuant), so read them
  as a direction; the two tools do not quantify the same PSMs, and MaxQuant reports fewer PSMs at
  its own 1% FDP on all three sets. The MaxQuant value is the lower one on SPS-MS3 and hrMS2, the
  andes value on RTS-MS3. FragPipe 21.1 reads 0.139, 0.055 and 0.214, and quantms
  (OpenMS IsobaricWorkflow, development version 96f263c) 0.072, 0.037 and 0.180 on the same three
  acquisitions. On hrMS2 only 4,502 of the 13,474 PSMs have a reporter scan whose range covers the
  reporter ions (3,481 of them pass the gate; 388 fail the purity and 633 the PEP criterion).
  `--tmt-unmix on` on hrMS2 sets the median knockout ratio to 0.000 at a background CV of 0.123
  (0.079 without).
  The MaxQuant PSM column counts all its PSMs; on hrMS2 most scans do not record the reporter
  ions, so on the PSMs that carry reporters the evaluated counts are andes 3,278, FragPipe 3,038
  and MaxQuant 2,241, and no quality gate can quantify more than the 4,502 reporter-bearing andes
  PSMs.
  Verdict: comparable to MaxQuant on SPS-MS3, 🟢 ahead on RTS-MS3, 🔴 behind on hrMS2 interference
  by default (🟢 well ahead with `--tmt-unmix on`); 🟢 most reporter-bearing PSMs quantified on hrMS2.
- *Known ratios (PXD007683):* 149,121 rank-1 target PSMs at 1% over the 11 fractions, 128,265
  quantified; yeast log2 ratios against the human-normalised expectation:

  | | A/B (expected 1.078) | A/C (1.688; nominal 3:1 = 1.585) | B/C (0.610) | protein-level MAE A/B, A/C, B/C | PSMs |
  |---|---:|---:|---:|---|---:|
  | andes, default gate | 1.191 | 1.589 | 0.397 | 0.230, 0.213, 0.213 | 149,121 |
  | andes, without the gate | 1.109 | 1.473 | 0.353 | 0.238, 0.266, 0.258 | 149,121 |
  | MaxQuant 2.7.0.0 | 1.263 | 1.700 | 0.424 | 0.222, 0.159, 0.186 | 97,107 |

  MaxQuant's yeast ratios are closer to the expectation on A/C and B/C (andes's on A/B) and its
  protein-level errors are lower on all three comparisons; it evaluates 8,427 proteins against
  10,797 for andes. 9.2 minutes and 17.8 GB for andes in one process (16 threads).
  Verdict: 🟢 **well ahead on depth** (+54% PSMs); 🔴 **behind on ratio accuracy** (protein MAE
  0.213 against 0.159; A/C 1.589 against 1.700, expected 1.688).
- *Match between runs for reporter ions (PXD014750, hrMS2).* Opt-in `--tmt --mbr` (hrMS2 only;
  refused for SPS-MS3 and ion-trap MS2). Data: the hrMS2 arm of PXD014750 (three TMT11 plexes of
  human and yeast peptides at changing mixing, each acquired twice; six files), searched against
  human + yeast + entrapment, Percolator `q ≤ 0.01`, `--rescore`, 16 threads. Truth: peptides
  of human origin only must be flat; peptides of yeast origin only must follow the plex's yeast
  profile when their identified PSMs in the other runs do. A value is wrong when its
  factor-normalised, impurity-corrected profile deviates by more than 0.5 log2 (RMS) from that
  expectation (an evaluation threshold, not a cut-off of any tool); the excess is the wrong rate
  of transferred values minus that of the same tool's identified values at matched species and
  signal. Only cells with a truth are scored, so the tools' denominators differ.

  | tool | match between runs | identified rows | transferred (cells or rows) | wrong, transferred | wrong, identified (matched) | excess, points (95% CI) |
  |---|---|---:|---:|---:|---:|---|
  | andes `--tmt` | off | 199,730 | 0 | — | 1.96% (all scored rows) | — |
  | andes `--tmt --mbr`, all candidates | on | 199,730 | 7,841 (6,260 scored) | 3.90% | 3.50% | +0.40 (−0.04 to +0.89) |
  | andes `--tmt --mbr`, `quant_q ≤ 0.01` | on | 199,730 | 3,636 (3,267 scored) | 2.79% | 2.84% | −0.05 (−0.56 to +0.51) |
  | MaxQuant 2.7.0.0 | on (isobaric MBR) | 137,539 evidence rows | 7,499 cells | 3.77% | 2.24% | +1.53 (+1.11 to +1.98) |
  | MaxQuant 2.7.0.0 | off | 137,539 evidence rows | 0 | — | — | — |
  | SIMSI-Transfer 0.8.0 (MaRaCluster 1.04.1), stringency p10 | on | MaxQuant scans | 16,511 cells | 2.25% | 1.85% | +0.40 (+0.18 to +0.64) |

  The andes rows come from this refresh (build `0cfa03cc`; `--mbr` adds 3.9% rows, 1.8% at `quant_q ≤ 0.01`; PSMs through the quantification gate
  131,161 without it, 134,797 with it; cost
  with the flag off 3.4k CPU-s, 8.0 minutes, 14.6 GB, with it 3.7k CPU-s, 8.4 minutes, 14.6 GB,
  one run each on one node type). The
  MaxQuant and SIMSI-Transfer rows are the earlier analysis of the same PXD014750 files with the
  same truth and wrong threshold (isobaric MBR as in MaxQuant ≥ 1.6.12.0: rows of type
  `MULTI-MATCH-MSMS` carry reporters; the "off" run is the same parameter file with the option
  off). The MaxQuant row counts all evidence rows of the types that carry an MS/MS scan, before
  reverse and contaminant hits are removed; andes's are PSMs at `q ≤ 0.01`. With `quant_q ≤ 0.01`
  andes accepts 3,636 transfers at an excess whose interval spans −0.6 to +0.5 points; its
  unfiltered candidates and both tools carry an excess of 0.4–1.5 points, and SIMSI-Transfer
  transfers the most cells. The entrapment share of the accepted andes transfers is 1.6% (1.7%
  of all candidates).

  Verdict: 🟢 **well ahead on accuracy** (excess −0.05 points against +1.53 and +0.40); behind
  SIMSI-Transfer on the number of cells transferred.

### Quantification: large study (500 runs)

**Large study (PXD042233, 500 runs).** Label-free DDA from Q Exactive HF, HF-X and Exploris 480,
`.raw` read natively, human reviewed proteins + contaminants, 5 ppm, 2 missed cleavages,
`--top-n 1`, Carbamidomethyl C fixed, Oxidation M and protein N-terminal Acetyl variable. Each
run is searched as its own job (8 threads, `--quant-dir`); the per-run PINs share one header
(it depends only on the command line) and are concatenated into one Percolator run (16 threads,
`--subset-max-train 1000000`); then `quant library`, `quant extract` per run (8 threads) and
`quant finalize` (16 threads). Measured 2026-10-23 at `2051616b` on cluster nodes.

| stage | wall | CPU | peak memory |
|---|---|---:|---:|
| search, 500 jobs | 0.7 min per run (at most 1.0) | 19.5 h | 3.6 GB |
| PIN concatenation (53.1 million rows) | 17 s | — | — |
| Percolator | 12.9 min | 0.7 h | 26.0 GB |
| `quant library` | 6.6 min | 0.1 h | 5.6 GB |
| `quant extract`, 500 jobs | 2.1 min per run (two runs 45 and 72 min) | 121.5 h | 12.7 GB (median 4.7) |
| `quant finalize`, transfer q 0.01, with the QPX views and protein inference | 25.8 min | 4.1 h | 45.2 GB |
| `quant finalize`, transfer q 0.05, without QPX views | 18.0 min | 3.9 h | 11.7 GB |

Outputs: 30.6 million PSMs at `q ≤ 0.01`. At transfer q 0.01, 194,370 accepted transfers
(175,087 `transfer_ms1`, 19,283 `transfer_ms2`) and 294,362 / 243,895 / 127,315 precursors
quantified in at least 1 / 3 / 50 runs (identified only: 294,362 / 243,730 / 126,741); 9,208
protein groups at protein q ≤ 0.01. At transfer q 0.05, 4,622,956 accepted transfers and
294,362 / 250,847 / 141,734 precursors. Extract is 83% of the CPU, and most of its transfer fits
are rejected (121 million `below_threshold` rows at q 0.01).

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
T/E that `build_databases.sh` printed), and peptides carrying `79.96`. Expect about 6 minutes on 8 threads
(337 s on the benchmark VM). The slice width is pinned because it sets peak memory (7.9 GB at 80 Da).

### Glyco (commands; scripts in [`glyco/`](glyco/))

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

For the quick tier, search `MouseLiver-Z-T-1` alone, run Percolator on its `.glyco.pin` without
pooling, and pass `--run MouseLiver-Z-T-1` to `score_vs_truth.py` (a single-file run writes
SpecIds without a file name). About 6 minutes and 6 GB on 8 threads. Pool the deep tier's
fractions: one fraction alone has too few glyco decoys for a stable 1% threshold.

### Label-free quantification, ProteoBench (commands; converter tested)

Inputs: the six runs of a module and the HYE FASTA linked from the ProteoBench module
documentation (`ProteoBenchFASTA_MixedSpecies_HYE.fasta`, 31,889 entries).
- *QExactive:* `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.raw` of PXD028735 (about 21 GB).
- *Astral:* `LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw` of PXD070049 (about 15 GB).

The measured runs read mzML peak-picked from these files. `R` is each run's stem. Every flag not
shown is the andes default (20 ppm precursor tolerance, `--max-missed-cleavages 1`, Carbamidomethyl C
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

### Truth-design quantification suite (commands; harness tested on synthetic data)

The evaluation code is [`quant-eval/`](quant-eval/README.md) (Python 3 with numpy, pandas, scipy
and pyarrow; `python3 tests/test_smoke.py` runs its tests). It reads finished tables and runs no
search. Per set, produce the inputs and write one manifest (the format is in
`quant_eval/suite.py`) naming the design files, the FASTA, the finalize table of each andes
quant-model seed, the masked-run directories and each tool's output:

```bash
# andes (one search per run, Percolator over the set, then the three quant stages)
andes --spectrum $R.mzML --database db.fasta --output-pin $R.pin --quant-dir q/ --threads 8 \
      --predict rt,ms2 --chimeric --refine --max-missed-cleavages 2
andes quant library  --quant-dir q/ --psms target.psms.tsv --decoy-psms decoy.psms.tsv
andes quant extract  --quant-dir q/ --run $R
andes quant finalize --quant-dir q/ -o set --mbr-fdr 0.01 --seed 42      # also --seed 1 and 7
# masked copy of the library for the wrong-peak measurement: the same three stages on a
# quant directory whose library lacks 5% of the identified pairs (mask_pairs.parquet)
# MaxQuant 2.7.0.0: one parameter file per set; the MBR-off run changes only
#   <matchBetweenRuns>False</matchBetweenRuns> (and <matchBetweenRunsFdr>False</matchBetweenRunsFdr>)
# FragPipe: IonQuant with --normalization 0
python3 -m quant_eval.suite --manifest manifest.json --out results/ --seeds 42,1,7
python3 -m quant_eval.board --results results/ > scoreboard.md
```

The `X_I` arm of the board is andes with match between runs off, `MQ_noMBR` is a second MaxQuant
entry of the manifest. ProteoBench submission files for the Q Exactive and Astral modules come from
`python3 -m quant_eval.proteobench` (raw intensities of the identified rows and the transfers
that pass the chosen rule) or from `proteobench/andes_to_proteobench.py`; neither uploads
anything.

### TMT quantification and match between runs (commands, not a tested script)

For the match-between-runs comparison on PXD014750 hrMS2 (six files, TMT11 impurity sheet, human +
yeast + entrapment database), the two runs differ by one flag:

```bash
andes --spectrum f1.mzML ... --spectrum f6.mzML --database hy_entrap.fasta --mods mods_tmt11.txt \
      --protocol TMT --rescore --fdr 0.01 --tmt tmt11 --tmt-level 2 --tmt-correction lot.txt \
      --threads 16 --output-tsv run.tsv [--mbr]
```

The truth evaluation reads `run.tmt.parquet` (columns `evidence_type`, `quant_q`, `quant_gate_passed`
and the channel columns) as described in the TMT paragraph of section 2.

### TMT quantification, knockout (commands, not a tested script)

The PXD017803 files `Fig2_m23142_hrMS2.raw` and `Fig2_m23145_SPSms3.raw`, reviewed yeast
proteins, and a mods file with Carbamidomethyl C and TMTpro (304.20715) on K and the peptide
N-terminus as fixed modifications, plus Oxidation M as variable. The measured runs read mzML
converted from these files (both files of each acquisition in one process), with 16 threads:

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
- **Entrapment design.** Shuffled entrapment is built peptide by peptide
  ([`build_shuffled_entrap.py`](build_shuffled_entrap.py)): each tryptic peptide of each target
  protein is shuffled with its C-terminal residue fixed, as in the peptide-level entrapment
  construction of Wen et al., Nat. Methods 2025. Every entrapment peptide then has a target of
  the same composition and mass. The earlier builder shuffled whole proteins: 89.6% of its
  entrapment peptides on the HYE database had no same-composition target peptide, and on the
  same Astral search it read 0.81% where the peptide-level databases read 0.92–1.02% and a
  foreign proteome (reviewed *Arabidopsis*, T/E 2.228) 0.99% (#225). On phospho it agreed with
  the peptide-level databases (1.05% against 0.95–1.04%). FDPs marked "whole-protein shuffle"
  were measured that way and may read low.
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

- **Label-free quantification rests on two datasets and one Percolator seed** (three for the
  error control), and its best scores use andes's own normalisation (`quantity_normalized`; the
  raw `quantity` is the quantms contract). On Astral andes's CV (0.113) is above most
  submissions, and its match-between-runs entrapment error is measured on QExactive only. The
  leaderboard was scored by `proteobench` 0.17.0 and andes by 0.18.6; several engines' transfer
  error control could not be checked.
- **Entrapment designs can disagree.** On UPS1 the natural *E. coli* entrapment reads the
  nominal 1% as 1.95–2.07% and a peptide-level shuffle as 0.92–0.94%; on the timsTOF astrocyte
  runs natural yeast reads 1.17–1.54% and the shuffle 0.94–1.12% on one draw (0.70–0.96% on five
  further draws). The shuffle is the primary
  calibration reading: the natural null is contaminated by conserved homologs (and, for UPS1, by
  possible real *E. coli* material), as the next paragraphs show. The 500-run study has no
  entrapment in this measurement.
- **Not re-measured at `2ee76dd6` or `2051616b`:** the glyco deep tier, the compared phospho FDP (Comet), the
  `--chimeric` and `--refine` opt-in table, the non-default `--predict` settings and the
  phospho agreement with MaxQuant.
- **Single files.** Phospho is one file, peptide-level only; `--refine` is measured on Astral
  only.
- **Glyco composition on fucose-rich tissue.** On heart and lung the backbone agrees with
  pGlyco2 94–99% of the time but the composition only 68–78% (97% on liver).
- **Not benchmarked:** iTRAQ, phosphosite localisation, timsTOF quantification. Java MS-GF+ was
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
