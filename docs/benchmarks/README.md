# andes benchmarks

andes is benchmarked on public PRIDE datasets for identification (high-res, TMT-labelled and
low-res LFQ, phospho, glyco) and for quantification (label-free with match between runs, and
TMT). Identifications are counted after Percolator at `q ≤ 0.01`. Where the database carries an
entrapment component, the page also gives the **true** false discovery proportion (FDP), which
is what tells you whether the nominal 1% holds. Other tools are named with the version compared.
Their numbers come either from runs on the same host with the same rescoring, or from a public
leaderboard. Section 2 gives the date and the andes commit behind each row.

## 1. Summary

Builds: identification `2051616b` (2026-10-23; timsTOF and entrapment rows `2ee76dd6`, 2026-11-07/09),
label-free quantification `09e1c00d` (2026-11-05), TMT `0cfa03cc`; Percolator 3.7.1 `-Y`.

One row per metric. 🟢 andes better · 🟡 similar · 🔴 andes worse, with how much: **minor** (usable
as is), **moderate** (usable with care), **major** (results hurt). ⚪ no other tool measured.
Numbers are andes against each tool on the same data; details and settings in section 2.

**Identification** (PSMs at Percolator `q ≤ 0.01`, 8 threads)

| data | metric | andes | other tool | result |
|---|---|---:|---|---|
| Astral HCD, PXD070049 | PSMs | 52,169 | Comet 31,435 | 🟢 +66% |
| Astral HCD | wall time | 114–117 s | Comet 217 s | 🟢 2x faster |
| Astral HCD | true FDP at 1% | 0.89–1.08% | — | 🟢 on target |
| timsTOF DDA-PASEF, 4 runs | PSMs | 15.6k–58.6k | Comet 14.1k–50.8k | 🟢 +11% to +17% |
| timsTOF DDA-PASEF | PSMs at a true FDP of 1% | 15.9k–58.8k | Comet 14.0k–49.5k | 🟢 +14% to +21% |
| timsTOF DDA-PASEF | true FDP at 1% | 0.79–0.91% | Comet 1.14–1.31% | 🟢 on target, Comet above |
| timsTOF DDA-PASEF | wall time | 54–62 s | Comet 159–222 s | 🟢 3x faster |
| TMT low-res, PXD007683 | PSMs | 13,516 | Comet 10,504 | 🟢 +29% |
| Phospho, PXD007653 | PSMs | 38,596–38,629 | Comet 33,984 (at 1.77% true FDP) | 🟢 +14%, at 0.80–0.99% FDP |
| Low-res UPS1, PXD001819 | PSMs | 16,567–16,675 | Java MS-GF+ 15,904 | 🟡 +4% to +5% |
| Low-res UPS1 | true FDP at 1% | 0.92–0.94% | — | 🟢 on target |
| TMT high-res, PXD017803 | PSMs | 6,628–6,686 | — | ⚪ |
| Glyco, PXD005553 | reference glycoPSMs confirmed | 79–90% | pGlyco2, MSFragger-Glyco (references) | 🔴 10–21% missed (moderate) |

**Quantification** (label-free at `09e1c00d`, TMT at `0cfa03cc`; DE counts are the common-tested
family, sign-aware, MBR on)

| data | metric | andes | MaxQuant 2.7 | FragPipe | result |
|---|---|---:|---:|---:|---|
| ProteoBench QE | median \|ε\| (lower is better) | 0.161 | — | 0.190 (54k precursors; 0.200 at 63.7k) | 🟢 15% lower |
| ProteoBench QE | precursors in ≥ 3 runs | 63,796 | — | 63,701 | 🟡 same |
| ProteoBench QE | CV | 0.117 | — | — (best entry 0.120) | 🟢 lowest |
| ProteoBench Astral | median \|ε\| | 0.186 | 0.160 (20k precursors) | 0.212 (39k) | 🟢 vs FragPipe; 🔴 vs MaxQuant at half the depth (minor) |
| ProteoBench Astral | CV | 0.112 | — | — (Sage 0.070) | 🔴 1.6x (minor) |
| Truth designs, 8 sets | precursors quantified | 10.7k–93.6k | | | 🟢 +19% to +51% on 7 sets, 🟡 on 1 |
| PXD009740 / PXD009815 / QE / Astral | DE true positives | 1,300 / 844 / 1,047 / 685 | 1,304 / 844 / 857 / 650 | 1,237 / 849 / 946 / 689 | 🟡 / 🟡 / 🟢 +22% / 🟡 |
| PXD009740 / PXD009815 / QE / Astral | DE false calls, no fold-change cut | 305 / 153 / 110 / 296 | 129 / 92 / 37 / 210 | 302 / 301 / 43 / 389 | 🔴 1.4–3.0x MaxQuant (moderate); mostly real small shifts, see note |
| PXD009740 / PXD009815 / QE | DE false calls at \|log2 FC\| ≥ 0.5 (seed 42) | 167 / 39 / 11 | 74 / 36 / 18 | 259 / 102 / 20 | 🔴 2.3x / 🟡 / 🟢 |
| PXD009740 / PXD009815 / QE / Astral | DE true positives / false calls, MaxQuant matches filtered at q ≤ 0.01 | 1,156 / 763 / 983 / 568; 272 / 151 / 105 / 297 | 1,077 / 725 / 851 / 518; 131 / 62 / 42 / 154 | | 🟢 +5% to +15% true positives; 🔴 1.9–2.5x false calls (moderate) |
| PXD009740 / PXD009815 / QE / Astral | precision of unchanged proteins (within-condition spread, tool ÷ andes) | 1 | 1.15 / 1.55 / 1.26 / — | 1.40 / 1.30 / 1.25 / — | 🟢 tighter than both |
| UPS1 lowest 3 levels, PXD009740 / PXD009815 | precursors within 2x of truth | 502 / 24 | 417 / 23 | 458 / 50 | 🟢 / 🔴 2x fewer than FragPipe (minor) |
| Truth designs, all sets | bright identified cells off by > 2x | 1.8–7.2% | 0.3–1.6% | 0.3–3.8% | 🔴 3–17x MaxQuant per set (moderate) |
| Truth designs, 6 sets | protein CV | 4.4–10.1% | 7.0–11.6% | 4.9–13.5% | 🟢 on 3 sets, 🔴 on 2 (minor) |
| Known-absent yeast, PXD014415 | false-transfer rate | 0.20% (0.08% at `quant_q ≤ 0.01`) | | 0.09% | 🔴 2x at default (minor) |
| Blank run, PXD016921 | rows reported in the blank | 150 (0 after #333) | 105 (29 with matches filtered) | 8 | 🔴 at `09e1c00d`; 🟢 on current main |
| TMT knockout, PXD017803 hrMS2 | interference (lower is better) | 0.192 (0.000 with `--tmt-unmix on`) | 0.156 | | 🔴 default, 🟢 with the option |
| TMT known ratios, PXD007683 | PSMs | 149,121 | 97,107 | | 🟢 +54% |
| TMT known ratios, PXD007683 | protein ratio error | 0.213 | 0.159 | | 🔴 1.3x (moderate) |
| TMT match between runs, PXD014750 | excess wrong-value rate of transfers | −0.05 points | +1.53 | | 🟢 (SIMSI-Transfer +0.40) |
| TMT match between runs | cells transferred | 3,636 | | | 🔴 4.5x fewer than SIMSI-Transfer (minor, opt-in) |
| Large study, PXD042233 (500 runs) | precursors in ≥ 3 runs | 243,895 | | | ⚪ |

Note on DE false calls: andes' protein values are tighter than the other tools', so small real
shifts of the "unchanged" species reach significance. On PXD009740, 296 of 320 false calls are
against its highest UPS1 level; without that level andes has 955 true positives and 24 false calls
against MaxQuant's 953 and 31. Analysis in section 2.3, *Differential expression*.

Settings, error control and cost of every row: section 2.

**Opt-in modes** (Astral, PSMs / distinct peptides, Percolator seeds 42, 1 and 2; measured
2026-10-05 to 2026-10-13 on earlier builds and not re-measured at `2051616b`):

| mode | andes | without the mode |
|---|---|---|
| `--chimeric` (second and later peptides per spectrum), at an entrapment true FDP of 1% | 66,377–66,424 / 29,384–29,400 | 43,076–43,145 / 26,723–26,746 |
| `--refine` (PTM discovery, high-res only) | 48,312 PSMs, 28,728 peptides at `q ≤ 0.01`; its second-pass PSMs at ~1.2–1.8% FDP | 47,080 / 28,932 |

## 2. The benchmarks

Icons and size words are as in section 1. Identification rows were measured on 2026-10-23 at
`2051616b`, except UPS1 and timsTOF (2026-11-07 and 2026-11-09 at `2ee76dd6`). Label-free
quantification rows were measured on 2026-11-05 at `09e1c00d`, TMT rows at `0cfa03cc`. A row says
so when it differs. PSM counts, wall time and memory of Astral, TMT, UPS1 and phospho come from
the benchmark VM (8 threads) on the plain databases the Comet rows used. True FDP and all other
rows come from a cluster node (8 threads) with an entrapment component added to the database
(rules in section 4). Percolator 3.7.1 `-Y --num-threads 8`, seeds 42, 1 and 2 unless a setup
line says otherwise. "At a true FDP of 1%" means the deepest score-sorted list whose
entrapment-estimated FDP stays at or below 1%.

### 2.1 Identification

#### Astral HCD (PXD070049)

Setup: one 15-minute DDA run of the ProteoBench human/yeast/E. coli sample, HYE FASTA (31,889
entries); andes `2051616b`, 2026-10-23; Comet 2025.01, same VM, matched parameters, same
Percolator, 2026-10-05.

| metric | andes | Comet 2025.01 | result |
|---|---:|---:|---|
| PSMs at `q ≤ 0.01` | 52,169–52,437 | 31,435 | 🟢 +66% |
| peptides at `q ≤ 0.01` | 31,147–31,238 | — | ⚪ |
| true FDP at nominal 1%, peptide-level entrapment (T/E 1.075; 9 database × seed runs) | 0.89–1.08% | — | 🟢 on target |
| PSMs / peptides at a true FDP of 1% | 47,267–47,977 / 28,141–28,857 | — | ⚪ |
| true FDP, *Arabidopsis* entrapment (T/E 2.228) | 0.99–1.05% (49,965–50,185 PSMs) | — | 🟢 on target |
| wall time; CPU time; memory | 114–117 s; 569–636 CPU-s; 4.8 GB | 217 s | 🟢 2x faster |

#### TMT-labelled, low-res (PXD007683, run a05058)

Setup: one ion-trap CID run of a TMT 10-plex sample; human + yeast reviewed proteins (26,483
sequences); andes `2051616b`; Comet 2025.01 as above.

| metric | andes | Comet 2025.01 | result |
|---|---:|---:|---|
| PSMs at `q ≤ 0.01` | 13,443–13,516 | 10,504 | 🟢 +29% |
| peptides at `q ≤ 0.01` | 11,852–11,922 | — | ⚪ |
| true FDP, peptide-level entrapment (T/E 1.073) | 0.67–0.76% (13,325–13,371 PSMs at 1%) | — | 🟢 on target |
| true FDP, *Arabidopsis* entrapment (T/E 2.021) | 0.83–1.05% (13,161–13,388 PSMs at 1%) | — | 🟢 on target |
| wall time; CPU time; memory | 40–47 s; 200–216 CPU-s; 1.7 GB | — | ⚪ |

#### TMT-labelled, high-res (PXD017803, TKO hrMS2)

Setup: `Fig2_m23142_hrMS2` (TMTpro, Orbitrap MS2), `--protocol TMT`, reviewed yeast + cRAP with
peptide-level entrapment (T/E 1.039); andes `2051616b`.

| metric | andes | other tool | result |
|---|---:|---|---|
| PSMs at `q ≤ 0.01` | 6,628–6,686 | — | ⚪ |
| true FDP at nominal 1% | 0.98–1.31% | — | ⚪ up to 0.3 points above nominal |
| PSMs at a true FDP of 1% | 6,607–6,662 | — | ⚪ |

#### Low-res LFQ (PXD001819, UPS1_5000amol_R1)

Setup: one run of the UPS1 standard in yeast; andes `2ee76dd6`, 2026-11-07, documented defaults
(2 missed cleavages, default since #300); Percolator seeds 42, 1 and 7; all-proteins entrapment
rule; two entrapment designs, each its own database. Java MS-GF+ v20240326 is from an earlier
session and was not re-run.

| metric | andes | Java MS-GF+ | result |
|---|---:|---:|---|
| PSMs at `q ≤ 0.01` (shuffle / natural database) | 16,567–16,623 / 16,632–16,675 | 15,904 | 🟡 +4% to +5% |
| peptides (shuffle / natural) | 4,176–4,200 / 4,116–4,227 | — | ⚪ |
| true FDP, peptide-level shuffle (T/E 1.033; 75–77 hits) | 0.92–0.94% | — | 🟢 on target |
| PSMs at a true FDP of 1%, shuffle | 16,617–16,733 | — | ⚪ |
| true FDP, natural *E. coli* (4,531 `ENTRAP_` proteins, T/E 2.451; 94–100 hits) | 1.95–2.07% | — | ⚪ context, not a calibration verdict (see first bullet) |
| PSMs at a true FDP of 1%, natural | 14,786–15,045 | — | ⚪ |
| natural, hits on conserved proteins removed (49–51 hits) | 1.01–1.06% (1.07–1.11% with T/E 2.630); 16,505–16,607 PSMs at 1% | — | 🟡 near nominal |
| andes wall; CPU; memory | 13–15 s; 71–79 CPU-s; 1.4 GB | — | ⚪ |
| Percolator, each of 3 seeds | 38–40 s; about 190–200 CPU-s; 1.2 GB | — | ⚪ |

- The natural reading is not a clean calibration. Of the 94 hits at seed 42, 45 sit on 8
  conserved groups (EF-Tu 32 PSMs from two peptides, ornithine carbamoyltransferase 7, enolase 5).
  UPS1 is made in *E. coli*, so host proteins may really be in the sample; the data cannot
  separate that from a mis-assigned yeast peptide. The shuffle (0.92–0.94%) is the primary reading.
- The any-protein rule (a PSM is entrapment if any of its proteins is) reads 3.09–3.23% on the
  same runs because shared peptides become false hits; the all-proteins rule is used here.
- T/E 2.451 is for the andes digest (2 missed cleavages). With `--max-missed-cleavages 1` (T/E
  2.351) the same build gives 16,238–16,367 PSMs, natural 1.94–1.99% and shuffle 0.88–0.96%
  (shuffle T/E 1.039); other digests give T/E 2.370–2.550.

#### Phospho-enriched (mouse liver, EasyPhos, single file)

Setup: `20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw` as mzML, not in any model's training
data; 80 Da index slices; mouse reviewed proteins plus the whole-protein shuffled entrapment the
Comet row used; andes `2051616b`; Comet 2025.01, 2026-10-06.

| metric | andes | Comet 2025.01 | result |
|---|---:|---:|---|
| PSMs at `q ≤ 0.01` | 38,596–38,629 (27,693–27,724 phospho-bearing) | 33,984 | 🟢 +14% |
| peptides at `q ≤ 0.01` | 20,422–20,461 | — | ⚪ |
| true FDP, 3 peptide-level shuffles (T/E 1.061) | 0.80–0.99% (38,604–38,843 PSMs at 1%) | 1.77% (whole-protein shuffle) | 🟢 on target; Comet above |
| true FDP, *Arabidopsis* entrapment (T/E 1.403) | 1.22–1.36% (38,298–38,412 PSMs at 1%) | — | 🔴 minor: above nominal |
| wall time; CPU time; memory | 337 s; 2,160 CPU-s; 7.9 GB | — | ⚪ |
| agreement with depositors' MaxQuant IDs (PEP ≤ 0.01, 23,563 scans) | 82.4% covered, 95.3% same-scan backbone | — | ⚪ |

- Phosphosite localisation is not scored.
- Comet's 1.77% is on the whole-protein shuffle, which can read low.
- The MaxQuant agreement was measured 2026-10-07 at `40774aca` and not repeated.

#### timsTOF DDA-PASEF (4 runs)

Setup: native `.d`, 20 ppm precursor, 2 missed cleavages, `--top-n 1`, 80 Da slices, Carbamidomethyl
C fixed, Oxidation M and protein N-terminal Acetyl variable; reviewed human proteins with a shuffle of
the human targets (T/E 1.050) or reviewed yeast (natural, T/E 3.784) as entrapment; andes
`2ee76dd6`, Percolator seeds 42, 1 and 7 for the single-draw tables, 2026-11-07 to 2026-11-09. Comet 2025.01 on tdf2mzml mzML, 20 ppm, isotope error −1 to 3, `fragment_bin_tol 0.02`,
internal reversed decoys; both tools on the same six shuffle draws (seeds 1, 11–15), Percolator seed 42, 2026-11-09; means over draws.

| metric (4 runs) | andes | Comet 2025.01 | result |
|---|---:|---:|---|
| PSMs at `q ≤ 0.01` | 15.6k–58.6k | 14.1k–50.8k | 🟢 +11% to +17% |
| peptide sequences | 3,865–24,146 | 3,396–22,007 | 🟢 +10% to +14% |
| PSMs at a true FDP of 1% | 15.9k–58.8k | 14.0k–49.5k | 🟢 +14% to +21% |
| true FDP at nominal 1% | 0.79–0.91% | 1.14–1.31% (23 of 24 draws above 1%) | 🟢 on target; Comet above |
| PSMs on the same mzML | — | — | 🟢 +8% to +16% (native `.d` adds 1.2–2.6%) |
| wall time | 54–62 s | 159–222 s | 🟢 3x faster |
| memory | 2.3–2.7 GB | 6.0–6.4 GB | 🟢 2.2–2.6x less |

Per run (mean over six draws):

| run | tool | PSMs | true FDP | PSMs at 1% | peptides | wall | memory |
|---|---|---:|---:|---:|---:|---:|---:|
| PXD053748 `..._4363` (plasma, 30 min) | andes | 15,642 | 0.79% | 15,880 | 3,910 | 55 s | 2.4 GB |
| | Comet | 14,143 | 1.14% | 13,964 | 3,493 | 159 s | 6.2 GB |
| PXD053748 `..._4367` (plasma, 50 min) | andes | 18,916 | 0.90% | 19,133 | 3,865 | 56 s | 2.5 GB |
| | Comet | 16,139 | 1.15% | 15,865 | 3,396 | 176 s | 6.4 GB |
| PXD064202 `..._2442` (astrocytes) | andes | 36,092 | 0.88% | 36,279 | 14,682 | 54 s | 2.3 GB |
| | Comet | 30,848 | 1.31% | 29,957 | 13,248 | 186 s | 6.0 GB |
| PXD064202 `..._2403` (astrocytes) | andes | 58,578 | 0.91% | 58,828 | 24,146 | 62 s | 2.7 GB |
| | Comet | 50,814 | 1.31% | 49,517 | 22,007 | 222 s | 6.1 GB |

andes alone, one shuffle draw and the natural yeast design (PSMs / true FDP; "no conserved" removes
hits on yeast proteins with a human homolog at 40% identity over 60% length, 907 of 6,733):

| run | shuffle PSMs | shuffle true FDP | shuffle PSMs at 1% | natural PSMs / true FDP / no conserved / PSMs at 1% |
|---|---:|---:|---:|---|
| `..._4363` | 15,411–15,474 | 0.52–0.60% | 15,938–15,952 | 16,062–16,108 / 0.86–1.01% / 0.80–0.95% / 16,076–16,219 |
| `..._4367` | 18,900–18,957 | 0.82–0.94% | 19,031–19,152 | 19,456–19,507 / 0.91–0.98% / 0.79–0.86% / 19,528–19,576 |
| `..._2442` | 36,052–36,125 | 0.94–1.01% | 36,047–36,178 | 36,873–36,982 / 1.44–1.54% / 0.93–1.02% / 36,040–36,071 |
| `..._2403` | 58,814–58,873 | 1.12% | 58,542–58,637 | 59,739–59,761 / 1.17–1.24% / 0.72–0.78% / 59,068–59,260 |

- One shuffle draw moves the true FDP: `..._2403` reads 1.12% on the draw above and 0.82–0.95% on
  five more; `..._2442` reads 0.70–0.96% on those. Over six draws the means are 0.91% (SD 0.11) and
  0.88% (SD 0.10), which matches counting noise on 250–320 hits. Draw-to-draw SD for Comet is
  0.08–0.21 points. The natural design on `..._2442` reads 1.44–1.54% (0.93–1.02% without
  conserved hits): 14 conserved groups hold 39 of 111 hits (yeast Ssb1/2: 15 PSMs from two
  peptides). A human sample cannot contain yeast proteins, so these are near-identical human
  peptides assigned to the homolog. The any-protein rule would read 9.0–9.7% on the astrocyte runs.
- The tools build decoys differently (andes permuted, Comet reversed), so part of Comet's excess
  over nominal may come from its decoys; the depth at a true FDP of 1% does not depend on them.
  Comet was not tuned for PASEF data. Wall times exclude conversion (5–7 min per run) and Percolator.
- andes cost per run: 40–52 s, 255–321 CPU-s, 2.1–2.5 GB on the natural database (49–58 s,
  307–371 CPU-s, 2.4–2.8 GB on the shuffle). Per-run q-values (#333) do not change these numbers
  because each timsTOF search holds one run.

#### Glyco (mouse liver, PXD005553)

Setup: `--glyco --glyco-species mouse --decoy-strategy sequon`, mouse proteins plus entrapment.
References for the same spectra are in [`glyco/truth/`](glyco/truth/): pGlyco2 (depositors) and
MSFragger-Glyco (reanalysis, PXD031032). "Confirmed" is the share of a reference's spectra that
andes identifies with the same peptide. Quick tier: fraction `MouseLiver-Z-T-1` as mzML, three
peptide-level shuffled databases × Percolator seeds 42, 1 and 2. Deep tier: all five fractions
pooled before Percolator, five seeds, measured 2026-09-05 at `14818d3e` (older, smaller glycan list,
mzML), not re-measured.

| metric | andes | reference | result |
|---|---:|---|---|
| glycoPSMs at 1% (quick tier) | 6,992–7,072 | — | ⚪ |
| true FDP (quick tier) | 0.61–1.35% | — | ⚪ 20–45 hits per run, mostly counting noise |
| spectra confirmed, quick tier | 88.2–88.4% | pGlyco2 | 🔴 11.6–11.8% missed (moderate) |
| spectra confirmed, quick tier | 89.3–89.7% | MSFragger-Glyco | 🔴 10.3–10.7% missed (moderate) |
| pGlyco2 spectra that lose the per-scan choice, deep tier | about 20% | pGlyco2 | 🔴 moderate |
| time; memory (quick tier, 8 threads) | 5.8–6.6 min; 5.4–5.6 GB | — | ⚪ |

- The open problem is selection, not generation: almost every reference spectrum yields a
  candidate row.
- Reading the `.raw` natively, as the commands in section 3 do, can give somewhat different counts.

### 2.2 Label-free quantification, ProteoBench

#### andes on the two ProteoBench modules

Setup: "Quant LFQ DDA ion" modules, six runs each: QExactive (PXD028735, Q Exactive HF-X) and
Astral (PXD070049, 15 min, 50 ng); human/yeast/E. coli at known A:B ratios. andes `09e1c00d`,
2026-11-05, defaults (`--predict rt`), Percolator through `--rescore` (seed 42), scored with
`proteobench` 0.18.6 ("Custom" format via [`proteobench/andes_to_proteobench.py`](proteobench/)).
Quantity: `quantity_normalized` (run loading and gradient drift removed) unless a row says raw.
Cells read precursors ≥ 3 runs / median |ε| / E. coli |ε| / CV.

| view | QExactive | Astral |
|---|---|---|
| default: transfer q ≤ 0.01, normalised | 63,796 / 0.161 / 0.263 / 0.117 | 36,675 / 0.186 / 0.252 / 0.112 |
| raw `quantity` (what quantms receives) | 63,796 / 0.222 / 0.309 / 0.245 | 36,675 / 0.204 / 0.254 / 0.113 |
| transfer q ≤ 0.05 | 71,299 / 0.175 / 0.441 / 0.132 | 43,286 / 0.205 / 0.349 / 0.127 |
| identified rows only (no match between runs) | 52,166 / 0.154 / 0.203 / 0.105 | 31,739 / 0.186 / 0.235 / 0.100 |
| precursors in ≥ 1 run / all 6 runs (default) | 79,492 / 41,370 | 49,012 / 24,393 |
| by abundance, quartile Q1 (lowest): \|ε\| / E. coli \|ε\| / CV | 0.279 / 0.461 / 0.195 | 0.319 / 0.406 / 0.194 |
| Q2 | 0.169 / 0.263 / 0.130 | 0.208 / 0.303 / 0.138 |
| Q3 | 0.136 / 0.205 / 0.105 | 0.161 / 0.212 / 0.103 |
| Q4 (highest) | 0.122 / 0.183 / 0.077 | 0.122 / 0.197 / 0.060 |
| opt-in `--predict rt,ms2` / `ms2,peptdeep-rt,deeplc-rt`: precursors (default on that build) | 67,906 / 68,077 (65,391) | 34,413 / 33,987 (33,876) |
| CPU-s, six runs: search + Percolator / library / extract / finalize per threshold | 5.8k / 4 / 0.2k / 0.2k | 5.8k / 4 / 30 / 0.1k |
| wall: search per run / extract per run / finalize | 3–4 min / 6 s / 34 s | 3–4 min / 2 s / 19 s |
| PSM FDP at `q ≤ 0.01` (entrapment, T/E 1.075) | 0.84–0.85% | — |
| entrapment-only share of precursors ≥ 3 runs, transfer q 0.01 / 0.05 | 0.20–0.21% / 0.26–0.27% | — |
| false-donor share of transfers, entrapment / identity-decoy estimate, q 0.01 | 0.22–0.26% / 0.28–0.29% | — |
| same, transfer q 0.05 | 0.37–0.40% / 0.42–0.47% | — |

- The error-control rows were measured on QExactive only, seeds 42, 1 and 2, at `09e1c00d`. The
  decoy estimate was above the measured share in all three seeds at both thresholds; three seeds do
  not make it an upper bound. The wrong-peak part has no entrapment check.
- The `--predict` rows are from 2026-10-21 at `a27a72c`, before the extraction change of #240.
- Raw `quantity` is less accurate than most leaderboard submissions; normalisation differs between
  tools, so the |ε| column does not compare extraction alone.

#### Against the ProteoBench leaderboard

Setup: public results of both modules (last updated 2026-07-16, read 2026-10-21, scored by
`proteobench` 0.17.0; andes by 0.18.6). Per tool the submission with the most precursors and, when
different, the one with the lowest median |ε|, at ≥ 3 runs. "Transfer control" is what limits the
error of match-between-runs values, from the submission's parameters, comment and output; "not
checked" means the parameter file could not be retrieved.

QExactive:

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | transfer control |
|---|---|---:|---:|---:|---:|---|
| andes, default | yes | 63,796 | 0.161 | 0.263 | 0.117 | transfer q ≤ 0.01 |
| andes, transfer q ≤ 0.05 | yes | 71,299 | 0.175 | 0.441 | 0.132 | transfer q ≤ 0.05 |
| andes, raw quantity | yes | 63,796 | 0.222 | 0.309 | 0.245 | as default |
| andes, identified rows only | no | 52,166 | 0.154 | 0.203 | 0.105 | none |
| FragPipe 23.0, most precursors | yes | 63,701 | 0.200 | 0.306 | 0.156 | IonQuant MBR ion FDR 0.01 |
| FragPipe 23.0, lowest \|ε\| | yes | 54,374 | 0.190 | 0.257 | 0.154 | as above |
| MaxQuant 1.5.3.30, most precursors | yes | 44,299 | 0.210 | 0.386 | 0.219 | not checked |
| MaxQuant 2.4.13.0, lowest \|ε\| | no | 29,737 | 0.192 | 0.263 | 0.199 | none |
| i2MassChroQ 1.0.16, most precursors | yes | 74,697 | 0.222 | 0.460 | 0.132 | no transfer setting |
| i2MassChroQ 1.2.0, lowest \|ε\| | no | 49,433 | 0.198 | 0.263 | 0.120 | none |
| AlphaPept 0.5.0, most precursors | yes | 37,760 | 0.219 | 0.281 | 0.165 | not checked |
| AlphaPept 0.5.3, lowest \|ε\| | no | 36,176 | 0.214 | 0.255 | 0.165 | none |
| Sage 0.15.0-beta.2 | yes | 81,933 | 0.388 | 0.877 | 0.393 | peptide q on IDs only; peak q not applied |
| MSAngel 2.2.10 | yes | 60,206 | 0.216 | 0.514 | 0.227 | not checked |
| ProlineStudio 2.3.0-SNAPSHOT | yes | 55,639 | 0.210 | 0.508 | 0.226 | not checked |
| WOMBAT 0.9.11, most precursors | yes | 55,943 | 0.211 | 0.692 | 0.235 | not checked |
| WOMBAT 0.9.11, lowest \|ε\| | yes | 41,016 | 0.195 | 0.394 | 0.213 | not checked |
| PEAKS 13 | no | 72,741 | 0.289 | 0.741 | 0.216 | none |
| quantms v1.3.1dev (OpenMS ProteomicsLFQ) | no | 33,040 | 0.209 | 0.387 | 0.224 | none |
| quantms, our run 2026-10-17, default | no | 34,677 | 0.210 | 0.401 | 0.225 | none |
| quantms, our run, feature-linking MBR | yes | 39,257 | 0.221 | 0.563 | 0.234 | feature quality score |

Astral:

| tool (submission) | MBR | precursors | median \|ε\| | E. coli \|ε\| | CV | transfer control |
|---|---|---:|---:|---:|---:|---|
| andes, default | yes | 36,675 | 0.186 | 0.252 | 0.112 | transfer q ≤ 0.01 |
| andes, transfer q ≤ 0.05 | yes | 43,286 | 0.205 | 0.349 | 0.127 | transfer q ≤ 0.05 |
| andes, raw quantity | yes | 36,675 | 0.204 | 0.254 | 0.113 | as default |
| andes, identified rows only | no | 31,739 | 0.186 | 0.235 | 0.100 | none |
| FragPipe 23.0, most precursors | yes | 39,057 | 0.212 | 0.250 | 0.127 | not checked |
| FragPipe 22.1-build02, lowest \|ε\| | yes | 29,105 | 0.181 | 0.405 | 0.078 | not checked |
| MaxQuant 2.6.7.0 | yes | 20,283 | 0.160 | 0.158 | 0.092 | within a condition only; not checked |
| i2MassChroQ 1.2.6, most precursors | yes | 36,469 | 0.197 | 0.417 | 0.095 | not checked |
| i2MassChroQ 1.2.9, lowest \|ε\| | no | 24,803 | 0.178 | 0.245 | 0.077 | none |
| Sage 0.14.6 | yes | 30,319 | 0.167 | 0.338 | 0.070 | peptide q on IDs only; peak q not applied |
| DIA-NN 2.6.0, most precursors | yes | 93,237 | 0.277 | 0.691 | 0.159 | run and global precursor q ≤ 0.01 |
| DIA-NN 2.6.0, lowest \|ε\| | yes | 16,363 | 0.123 | 0.168 | 0.079 | as above plus MS1-level q ≤ 0.01 |

On the precursors andes shares with a submission (default, transfer q 0.01; median |ε| / E. coli |ε|
/ CV; QExactive only, Astral per-submission tables could not be retrieved):

| submission | shared | submission | andes raw | andes normalised | only in the submission |
|---|---:|---|---|---|---|
| FragPipe 23.0, most precursors | 51,722 | 0.189 / 0.274 / 0.155 | 0.207 / 0.276 / 0.241 | 0.147 / 0.225 / 0.108 | 11,979 at 0.272 / 0.708 |
| FragPipe 23.0, lowest \|ε\| | 48,938 | 0.187 / 0.243 / 0.155 | 0.203 / 0.268 / 0.241 | 0.144 / 0.207 / 0.108 | 5,436 at 0.220 / 0.510 |
| i2MassChroQ 1.0.16, most precursors | 50,910 | 0.207 / 0.348 / 0.125 | 0.209 / 0.291 / 0.243 | 0.152 / 0.241 / 0.112 | 23,787 at 0.269 / 0.948 |
| i2MassChroQ 1.2.4, 57,787 precursors | 49,853 | 0.202 / 0.348 / 0.123 | 0.209 / 0.291 / 0.242 | 0.151 / 0.242 / 0.112 | 7,934 at 0.279 / 0.835 |
| Sage 0.15.0-beta.2, our reproduction | 54,633 | 0.309 / 0.613 / 0.319 | 0.215 / 0.291 / 0.245 | 0.156 / 0.246 / 0.114 | 25,724 at 0.683 / 1.364 |

- QExactive: at least six submissions report more precursors than andes (Sage and the largest
  i2MassChroQ and PEAKS submissions among them), at a median |ε| of 0.214–0.388. FragPipe 23.0
  reaches 0.190 at 54,374 precursors; i2MassChroQ 1.2.0 has a CV of 0.120 without match between runs.
- Astral: MaxQuant 2.6.7.0 and Sage 0.14.6 report fewer precursors than andes at lower median
  |ε| (0.160, 0.167); i2MassChroQ 1.2.6 about as many; FragPipe 23.0 more (39,057, 0.212); DIA-NN
  2.6.0 reports 93,237 at 0.277, and its 16,363-precursor submission has the lowest median |ε|
  (0.123). andes's CV (0.112) is above most submissions. Normalisation changes little here.
- Sage's per-submission table was not retrievable; the reproduction (2026-10-21, same version, same
  six mzML files, submission parameters, per-charge output) gives 80,357 precursors, |ε| 0.388,
  E. coli 0.878, CV 0.394 (submission: 81,933, 0.388, 0.877, 0.393). Most Sage precursors that
  andes lacks are charge states Sage never identified or that fail its own peak q-value of 0.01,
  which it reports but does not apply.

### 2.3 Label-free quantification, truth designs

Setup: eight sets with a known answer, quantified by andes, MaxQuant and FragPipe and read through
one evaluation ([`quant-eval`](quant-eval/README.md); rules in section 4). andes `09e1c00d`
defaults, run 2026-11-05 over the 2026-10-31 searches (`--predict rt,ms2 --chimeric --refine
--max-missed-cleavages 2`); `quant library` (`--mbr-library-fdr 0.01`), `quant extract`, `quant
finalize` (`--mbr-fdr 0.01`), quant-model seeds 42, 1 and 7 (PXD030780: 42); tables show seed 42.
Raw `quantity` column. MBR on = identified rows plus transfers with `transfer_q ≤ 0.01`; MBR off =
identified rows. MaxQuant 2.7.0.0 (Standard DDA, LFQ on, 1% PSM and protein FDR, 2 missed
cleavages, `Intensity` of the evidence table per precursor and run; MBR-off run re-run 2026-11-05,
MBR-on run from 2026-10). FragPipe 21.1 with IonQuant 1.10.12, MBR on (ion FDR 0.01), normalisation
off (`--normalization 0`), 2026-10 outputs, no MBR-off run.

| set | what it is | runs | truth |
|---|---|---:|---|
| PXD009740 | UPS1 at ten levels (170 amol–250 fmol) in constant yeast | 30 | UPS1 amount per run; yeast constant |
| PXD009815 | UPS1 at ten levels (10 amol–50 fmol) in yeast | 40 | as above |
| PXD001819 | UPS1 at nine levels in yeast; level and run order confounded | 27 | as above; counts and accuracy only (every tool reports 1,000+ yeast false positives in DE) |
| QE (ProteoBench) | human, yeast, *E. coli*; Q Exactive HF-X (PXD028735) | 6 | yeast 2:1, *E. coli* 1:4, human 1:1 (A vs B) |
| Astral (ProteoBench) | same mixture, Astral (PXD070049) | 6 | as above |
| PXD030780 | human, yeast, *E. coli*, five mixing levels A–E | 25 | species amounts per group |
| PXD014415 | human and human + yeast runs | 40 | yeast absent from 20 runs |
| PXD016921 | blank run among cell-input runs | 6 | no analyte in the blank |

- Settings are not matched: andes used `--chimeric`, `--refine` and PSM-level FDR; MaxQuant and
  FragPipe used protein-level 1% FDR. With matched settings (±10 ppm, isotope errors 0 to 2, no
  `--chimeric`, no `--refine`, seed 42) andes's precursors fall 31% on PXD009740 and PXD009815 (to
  27,226 and 29,600, level with the other tools), 12% on QE (81,987) and 6% on Astral (51,798),
  and DE true positives stay within seed spread (PXD009740 1,293, PXD009815 831, QE 1,049, Astral
  673). The depth lead on the UPS sets comes mostly from the search settings; the DE results come
  from the quantification. A protein-level 1% filter on andes removes 0.4–1.6% of precursors.
- MaxQuant 2.7 does not remove match rows with its MBR FDR option (it only fills the `Match
  q-value` column), so MaxQuant's MBR-on rows are unfiltered, while andes's transfers are filtered
  at `transfer_q ≤ 0.01`. Filtering MaxQuant at `Match q-value ≤ 0.01` cuts PXD001819 transfers from
  43,373 to 2,550 and blank-run transfers from 85 to 9. No MaxQuant run exists for PXD014415 or
  PXD030780.

#### Depth

Cells are precursor-by-run values. Each entry reads precursors / proteins (≥ 2 precursors) /
identified cells + transferred cells.

| set | andes, MBR on | MaxQuant, MBR on | FragPipe, MBR on | result |
|---|---|---|---|---|
| PXD009740 | 39,595 / 3,031 / 573,349 + 332,626 | 29,293 / 2,339 / 464,069 + 175,450 | 25,888 / 2,259 / 414,258 + 141,367 | 🟢 +35% |
| PXD009815 | 42,686 / 3,362 / 798,525 + 463,945 | 27,459 / 2,453 / 586,274 + 230,204 | 28,286 / 2,691 / 587,049 + 193,057 | 🟢 +51% |
| PXD001819 | 10,696 / 1,004 / 146,103 + 65,796 | 9,314 / 881 / 116,299 + 43,373 | 10,406 / 937 / 118,961 + 63,595 | 🟡 +3% |
| QE | 92,828 / 6,797 / 350,579 + 71,126 | 51,609 / 4,810 / 186,178 + 25,465 | 71,377 / 5,938 / 249,172 + 12,401 | 🟢 +30% |
| Astral | 55,246 / 4,901 / 212,633 + 28,812 | 24,745 / 2,855 / 93,669 + 29,072 | 46,347 / 4,334 / 168,481 + 46,343 | 🟢 +19% |
| PXD030780 | 93,566 / 8,602 / 1,047,714 + 560,171 | — | 65,586 / 6,885 / 824,966 + 351,518 | 🟢 +43% |
| PXD014415 | 76,099 / 5,470 / 1,494,110 + 801,869 | — | 61,771 / 5,030 / 1,175,139 + 550,247 | 🟢 +23% |
| blank run | 17,625 / 1,728 / 29,755 + 1,641 | 11,010 / 1,229 / 17,749 + 1,305 | 14,153 / 1,403 / 24,638 + 4,015 | 🟢 +25% |

Result is the precursor gain over the best other tool; settings are not matched (first bullet
above). With MBR off, precursor counts change by less than 0.2% for every tool and the transferred
cells are 0.

#### Accuracy and precision

Cells identified by every tool listed for the set. Entries read andes / MaxQuant / FragPipe (MBR on;
andes MBR-on accuracy equals MBR off because the cells are identified cells). Bright decile = top
tenth of predicted intensity; gross error = share of cells off by more than 2x. Protein CV is the
median replicate CV of matched features.

| set | cells off by > 2x (%) | bright-decile gross error (%) | precursor SD, top decile | protein CV (%) | protein CV, MBR off (%): andes / MaxQuant | result |
|---|---|---|---|---|---|---|
| PXD009740 | 2.63 / 1.55 / 4.02 | 3.64 / 1.13 / 2.37 | 0.132 / 0.116 / 0.126 | 8.7 / 10.5 / 13.5 | 10.2 / 10.5 | 🔴 accuracy; 🟢 CV |
| PXD009815 | 1.42 / 0.35 / 0.65 | 5.97 / 0.55 / 0.68 | 0.056 / 0.066 / 0.051 | 4.4 / 7.0 / 6.1 | 5.6 / 7.5 | 🔴 accuracy; 🟢 CV |
| PXD001819 | 1.81 / 0.86 / 1.80 | 1.79 / 0.46 / 3.76 | 0.044 / 0.060 / 0.038 | 5.7 / 8.3 / 5.2 | 6.8 / 7.9 | 🔴 vs MaxQuant; 🔴 CV vs FragPipe (minor) |
| QE | 3.75 / 2.92 / 5.84 | 6.44 / 1.62 / 2.64 | 0.208 / 0.184 / 0.206 | 9.2 / 11.6 / 11.5 | 9.6 / 10.3 | 🔴 accuracy; 🟢 CV |
| Astral | 1.62 / 0.69 / 1.13 | 4.39 / 0.26 / 0.66 | 0.057 / 0.072 / 0.041 | 5.8 / 7.0 / 4.9 | 6.0 / 6.9 | 🔴 accuracy; 🔴 CV vs FragPipe (minor) |
| PXD030780 | 2.63 / — / 1.44 | 7.19 / — / 0.30 | 0.112 / — / 0.096 | 10.1 / — / 10.2 | 12.1 / — | 🔴 accuracy; 🟡 CV |

- Bright-decile gross error is 1.8–7.2% for andes against 0.3–1.6% for MaxQuant (MBR on) and
  0.3–3.8% for FragPipe; accuracy is behind on every set. Protein CV is 4.4–10.1% for andes,
  7.0–11.6% for MaxQuant and 4.9–13.5% for FragPipe. MaxQuant MBR-off values are close to its MBR-on
  ones (gross error up to 1.65%).

#### Differential expression (BH 5% per pair)

Entries are true positives / false calls on the common-tested family, sign-aware (wrong-sign
calls, counted as false, in brackets). Seed 42. Across seeds 42, 1 and 7 andes MBR-on gives
1,288–1,302 / 305–350 on PXD009740, 841–847 / 153–188 on PXD009815, 1,047–1,050 / 110–120 on QE and
681–685 / 296–301 on Astral.

| set | andes, MBR on | andes, MBR off | MaxQuant, MBR on | MaxQuant, MBR off | FragPipe, MBR on | result |
|---|---|---|---|---|---|---|
| PXD009740 | 1,300 / 305 (1) | 1,068 / 180 (0) | 1,304 / 129 (0) | 1,047 / 188 (0) | 1,237 / 302 (0) | 🟡 true positives; 🔴 2.4x MaxQuant false calls (moderate) |
| PXD009815 | 844 / 153 (7) | 640 / 76 (0) | 844 / 92 (2) | 541 / 54 (2) | 849 / 301 (2) | 🟡 true positives; 🔴 1.7x MaxQuant false calls (moderate) |
| PXD001819 | 375 / 3,296 (0) | 266 / 1,833 (0) | 358 / 939 (0) | 240 / 861 (0) | 359 / 1,977 (0) | 🟡; not a decision set |
| QE | 1,047 / 110 (0) | 1,007 / 62 (0) | 857 / 37 (0) | 824 / 43 (0) | 946 / 43 (1) | 🟢 +22% over MaxQuant, +11% over FragPipe; 🔴 3.0x false calls (moderate) |
| Astral | 685 / 296 (1) | 640 / 251 (0) | 650 / 210 (2) | 465 / 141 (0) | 689 / 389 (1) | 🟡 true positives; 🔴 1.4x MaxQuant false calls (minor) |
| PXD030780 | 17,692 / 1,374 (23) | 14,416 / 915 (29) | — | — | 18,227 / 1,751 (22) | 🟡 |

MaxQuant matches filtered at `Match q-value ≤ 0.01` (true positives / false calls; PXD009740,
PXD009815, QE, Astral):

| tool | PXD009740 | PXD009815 | QE | Astral | result |
|---|---|---|---|---|---|
| andes, MBR on | 1,156 / 272 | 763 / 151 | 983 / 105 | 568 / 297 | 🟢 +5% to +15% true positives; 🔴 1.9–2.5x false calls (moderate) |
| MaxQuant, matches filtered | 1,077 / 131 | 725 / 62 | 851 / 42 | 518 / 154 | |

Spread of unchanged proteins within a condition (tool ÷ andes; PXD009740 / PXD009815 / QE): MaxQuant
1.15 / 1.55 / 1.26, FragPipe 1.40 / 1.30 / 1.25. True positives / false calls by minimal fold change
(seed 42 on main after #333, 2026-11-09; the table above is the mean over seeds at `09e1c00d`, e.g.
PXD009740 reads 1,287 / 320 here and 1,300 / 305 above):

| set | tool | all | \|log2 FC\| ≥ 0.3 | ≥ 0.5 | ≥ 1.0 |
|---|---|---|---|---|---|
| QE | andes | 1,049 / 113 | 1,049 / 83 | 1,046 / 11 | 476 / 0 |
| | MaxQuant | 857 / 37 | 857 / 35 | 852 / 18 | 303 / 2 |
| | FragPipe | 946 / 43 | 946 / 37 | 945 / 20 | 680 / 8 |
| PXD009740 | andes | 1,287 / 320 | 1,287 / 310 | 1,287 / 167 | 1,271 / 27 |
| | MaxQuant | 1,292 / 129 | 1,292 / 128 | 1,292 / 74 | 1,238 / 25 |
| | FragPipe | 1,224 / 302 | 1,224 / 302 | 1,224 / 259 | 1,220 / 86 |
| PXD009815 | andes | 853 / 183 | 852 / 60 | 846 / 39 | 800 / 21 |
| | MaxQuant | 842 / 92 | 842 / 66 | 838 / 36 | 778 / 6 |
| | FragPipe | 859 / 301 | 859 / 184 | 858 / 102 | 844 / 22 |

- At `|log2 FC| ≥ 0.5`, andes has 167 false calls against MaxQuant's 74 on PXD009740 (🔴 2.3x), 39
  against 36 on PXD009815 and 11 against 18 on QE; Astral was not recomputed. Still behind on
  PXD009815 at large fold changes (21 against 6 at ≥ 1.0) and in wrong-sign calls (18 against 2 at
  seed 42 after #333; 7 against 2 in the table above).
- Why: andes's protein values are tighter, so small real shifts of "unchanged" proteins reach
  significance; the truth labels assume no change at all, which stops holding below about 1.4-fold
  (|log2 FC| 0.5). The calls are not made by transferred values: from identified cells alone the
  called proteins keep 80–95% of their fold change. On PXD009740, 296 of 320 false calls are against
  the highest UPS1 level; without that level andes has 955 true positives and 24 false calls against
  MaxQuant's 953 and 31.
- Match between runs raises true positives for every tool that offers it (andes 1,068 to 1,300 on
  PXD009740, 1,007 to 1,047 on QE, 640 to 685 on Astral, 14,416 to 17,692 on PXD030780; MaxQuant
  1,047 to 1,304, 824 to 857, 465 to 650). andes's false calls rise with it (180 to 305, 62 to 110);
  MaxQuant's fall on PXD009740 and QE (188 to 129, 43 to 37) and rise on Astral and PXD009815.

#### UPS1 detection at the lowest three levels

Precursors within 2x of their own high-level reference, summed over the three lowest levels
(PXD009740 170, 382, 859 amol; PXD009815 10, 50, 100 amol; PXD001819 50, 125, 250 amol), MBR on.
"At FragPipe's row count" is the expected number within 2x if the tool reported only as many rows as
FragPipe, drawn at random.

| set | tool | within 2x | rows | rows within 2x (%) | median \|log2 err\| | at FragPipe's row count | result |
|---|---|---:|---:|---:|---:|---:|---|
| PXD009740 | andes | 502 | 1,600 | 60.3 | 0.84 | 436.7 | 🟢 vs MaxQuant; 🔴 vs FragPipe at equal rows (minor) |
| | MaxQuant | 417 | 1,642 | 38.9 | 1.20 | 340.7 | |
| | FragPipe | 458 | 1,167 | 71.1 | 0.64 | 458.0 | |
| PXD009815 | andes | 24 | 565 | 7.6 | 4.75 | 22.1 | 🔴 2x fewer than FragPipe (minor); 🟡 vs MaxQuant |
| | MaxQuant | 23 | 603 | 6.3 | 4.63 | 19.1 | |
| | FragPipe | 50 | 443 | 15.6 | 3.73 | 50.0 | |
| PXD001819 | andes | 0 | 53 | 0.0 | 6.00 | 0.0 | 🟡 none detects |
| | MaxQuant | 1 | 68 | 1.5 | 4.25 | 0.5 | |
| | FragPipe | 3 | 29 | 10.3 | 3.13 | 3.0 | |

Per level, MBR on (rows; within 2x; share within 2x; median |log2 error|):

| set | level (amol) | andes | MaxQuant | FragPipe |
|---|---:|---|---|---|
| PXD009740 | 170 | 295; 66; 33.6%; 1.29 | 323; 38; 14.2%; 1.97 | 189; 74; 52.9%; 0.96 |
| | 382 | 467; 153; 59.3%; 0.85 | 532; 127; 35.2%; 1.26 | 345; 133; 67.0%; 0.63 |
| | 859 | 838; 283; 70.3%; 0.72 | 787; 252; 51.6%; 0.99 | 633; 251; 78.8%; 0.58 |
| PXD009815 | 10 | 174; 0; 0.0%; 7.23 | 195; 0; 0.0%; 7.16 | 131; 1; 0.8%; 6.53 |
| | 50 | 194; 3; 2.1%; 4.71 | 189; 4; 2.1%; 4.39 | 140; 19; 17.1%; 3.97 |
| | 100 | 197; 21; 19.8%; 3.03 | 219; 19; 15.5%; 2.94 | 172; 30; 25.6%; 2.26 |

- With MBR off: PXD009740 andes 206 (546 rows, 62.6% within 2x, median error 0.84) against MaxQuant
  111 (360 rows, 48.3%, 1.01); PXD009815 andes 3 against MaxQuant 0; PXD001819 andes 0, MaxQuant 0.
- At equal rows FragPipe is ahead on both UPS1 sets (458 against andes 437; 50 against 22) and
  andes is ahead of MaxQuant on PXD009740 (437 against 341) and level on PXD009815 (22 against 19).
- At 10 and 50 amol every tool over-reads: median error +6.5 to +7.4 log2 at 10 amol and +4.0 to
  +4.7 at 50 amol (a factor of 2^4 to 2^7).

#### Absence: false transfers and the blank run

PXD014415: yeast rows among the transferred rows of the 20 human-only runs. PXD016921: rows
reported in the blank run. Classes (definition in section 4) are enriched / not specific / only
absent / no reference.

| metric | andes | MaxQuant 2.7 | FragPipe | result |
|---|---|---|---|---|
| yeast rows, identified (MBR off or on) | 451 | — | 208 | ⚪ |
| yeast rows, transferred; share of transferred rows | 769; 0.20% of 389,423 (seeds 0.19–0.20%) | — | 243; 0.09% of 266,898 | 🔴 2x at default (minor) |
| same at `quant_q ≤ 0.01` | 134; 0.08% of 162,219 | — | as above | 🟡 level |
| transferred classes | 202 / 567 / 0 / 0 | — | 27 / 52 / 0 / 164 | ⚪ |
| rate, enriched rows only | 0.052% | — | 0.010% | 🔴 (moderate) |
| rate, enriched plus no-reference rows | 0.052% | — | 0.072% | 🟢 |
| identified classes | 38 / 412 / 1 / 0 | — | 8 / 56 / 5 / 139 | ⚪ |
| identified rate, enriched only / with no-reference | 0.0053% / 0.0055% | — | 0.0014% / 0.0270% | 🔴 / 🟢 |
| rows in the blank run, MBR off / on | 85 / 150 (144 at `quant_q ≤ 0.01`; 0 on current main after #333) | 20 / 105 (29 with matches filtered) | 8 with MBR | 🔴 at `09e1c00d` |
| blank run, identified / transferred, MBR on | 85 / 65 | 20 / 85 | 8 / 0 | ⚪ |

- The reading flips with the no-reference rows, which the intensity profile cannot judge, so both
  are given. A tool's precursors that andes does not report cannot be classed.
- MaxQuant's MBR-on blank run had its MBR FDR option off; see the filtered-matches note above.

#### Masked wrong-peak rate (andes only)

At `quant_q ≤ 0.01`, pooled over runs; 5% of identified precursor-run pairs are masked from the
library, re-quantified as transfers, and counted when the apex lies further from the identified apex
than the fit's peak width. No other tool exposes this quantity.

| set | seeds 42 / 1 / 7 | MS1-only transfers, seed 42 |
|---|---|---:|
| PXD009740 | 0.50 / 0.50 / 0.49% | 0.82% |
| PXD009815 | 0.55 / 0.55 / 0.60% | 0.92% |
| PXD001819 | 0.46 / 0.41 / 0.43% | 1.01% |
| QE | 0.77 / 0.80 / 0.78% | 1.30% |
| Astral | 0.66 / 0.66 / 0.62% | 0.66% |
| PXD014415 | 0.49 / 0.50 / 0.48% | 1.09% |

- Transfers with MS2 evidence read 0.18–0.65% (seed 42). The blank-run set has 126–130 masked pairs
  and none was wrong. Among runs with at least 200 masked pairs, 0 to 5 per set (of 6 to 40) are
  above 1%.

#### Depth at matched true error (PXD014415)

| metric | andes | FragPipe | result |
|---|---:|---:|---|
| transferred rows at FragPipe's 0.09% (`transfer_q ≤ 0.003`) | 487,849 (0.88–0.89 × FragPipe's, over seeds) | 550,247 | 🔴 0.88–0.89x of FragPipe (minor) |
| transferred rows through `quant_q ≤ 0.015` | 428,888 | — | ⚪ |
| transferred rows at its own `quant_q ≤ 0.01` (0.08%) | 317,272 | — | ⚪ |

MaxQuant has no PXD014415 run. No tool reports a masked wrong-peak rate to compare.

### 2.4 Label-free quantification, large study

#### PXD042233 (500 runs)

Setup: label-free DDA from Q Exactive HF, HF-X and Exploris 480, `.raw` read natively, human
reviewed proteins + contaminants, 5 ppm, 2 missed cleavages, `--top-n 1`; one search job per run
(8 threads, `--quant-dir`), per-run PINs concatenated into one Percolator run (16 threads,
`--subset-max-train 1000000`), then `quant library`, `quant extract` per run, `quant finalize` (16
threads); andes `2051616b`, 2026-10-23, cluster nodes. No entrapment in this measurement.

| stage | wall | CPU | peak memory | output |
|---|---|---:|---:|---|
| search, 500 jobs | 0.7 min per run (at most 1.0) | 19.5 h | 3.6 GB | |
| PIN concatenation | 17 s | — | — | 53.1 million rows |
| Percolator | 12.9 min | 0.7 h | 26.0 GB | 30.6 million PSMs at `q ≤ 0.01` |
| `quant library` | 6.6 min | 0.1 h | 5.6 GB | |
| `quant extract`, 500 jobs | 2.1 min per run (two runs 45 and 72 min) | 121.5 h | 12.7 GB (median 4.7) | 121 million `below_threshold` rows at q 0.01 |
| `quant finalize`, transfer q 0.01, QPX views, protein inference | 25.8 min | 4.1 h | 45.2 GB | 194,370 transfers (175,087 `transfer_ms1`, 19,283 `transfer_ms2`); precursors in ≥ 1 / 3 / 50 runs 294,362 / 243,895 / 127,315 (identified only 294,362 / 243,730 / 126,741); 9,208 protein groups at protein q ≤ 0.01 |
| `quant finalize`, transfer q 0.05, no QPX views | 18.0 min | 3.9 h | 11.7 GB | 4,622,956 transfers; precursors in ≥ 1 / 3 / 50 runs 294,362 / 250,847 / 141,734 |

Extract is 83% of the CPU and most of its transfer fits are rejected. No other tool was run on this
set (⚪).

### 2.5 TMT quantification

#### Knockout interference and extraction (PXD017803)

Setup: yeast triple knockout, TMTpro 16-plex; two files each of hrMS2, SPS-MS3 and RTS-MS3, as
mzML, with the kit impurity sheet. andes defaults including the quantification gate (precursor purity
≥ 0.7, Percolator PEP ≤ 0.05), `--rescore`, 16 threads, `0cfa03cc`, 2026-11-05. MaxQuant 2.7.0.0
run 2026-10-18 on vendor files (reporter intensity corrected, 1% PSM and protein FDR). FragPipe 21.1;
quantms (OpenMS IsobaricWorkflow, `96f263c`). Interference is the median ratio of knockout to
non-knockout channel intensity over PSMs of the three knocked-out proteins (lower is better; hrMS2
only scans whose window starts below 126 Th). Reporter extraction and impurity correction match OpenMS
IsobaricAnalyzer exactly on the same scans.

| acquisition | andes PSMs at 1% (quantified) | andes interference | no gate / purity ≥ 0.9 | MaxQuant PSMs; interference | FragPipe | quantms | result vs MaxQuant |
|---|---:|---:|---|---|---:|---:|---|
| SPS-MS3 | 9,806 (8,635) | 0.062 | 0.117 / 0.023 | 6,883; 0.054 | 0.139 | 0.072 | 🟡 |
| RTS-MS3 | 14,270 (10,518) | 0.027 | 0.033 / 0.015 | 8,930; 0.038 | 0.055 | 0.037 | 🟢 |
| hrMS2 | 13,474 (3,481) | 0.192 | 0.202 / 0.135 | 9,677; 0.156 | 0.214 | 0.180 | 🔴 default (minor); 🟢 with `--tmt-unmix on` |

- `--tmt-unmix on` on hrMS2 gives a median knockout ratio of 0.000 at a background CV of 0.123 (0.079
  without it).
- Few knockout PSMs stand behind each value (34–91 andes, 26–65 MaxQuant); read them as a direction.
  The tools do not quantify the same PSMs, and MaxQuant reports fewer PSMs at its own 1% FDP on all
  three sets.
- On hrMS2 only 4,502 of 13,474 andes PSMs have a reporter scan covering the reporter ions (3,481
  pass the gate; 388 fail purity, 633 fail PEP). On those reporter-bearing PSMs the evaluated counts
  are andes 3,278, FragPipe 3,038 and MaxQuant 2,241: 🟢 most PSMs quantified.

#### Known ratios (PXD007683)

Setup: human + yeast TMT11 with known yeast ratios, 11 SPS-MS3 fractions, ion-trap CID MS2; andes
default gate, `--rescore`, 16 threads, `0cfa03cc`, 2026-11-05; MaxQuant 2.7.0.0 as above. Yeast log2
ratios against the human-normalised expectation (A/B 1.078; A/C 1.688, nominal 3:1 = 1.585; B/C 0.610).

| metric | andes | MaxQuant 2.7.0.0 | result |
|---|---:|---:|---|
| rank-1 target PSMs at 1% (quantified) | 149,121 (128,265) | 97,107 | 🟢 +54% |
| yeast log2 ratio A/B / A/C / B/C | 1.191 / 1.589 / 0.397 | 1.263 / 1.700 / 0.424 | 🟡 MaxQuant closer on A/C and B/C, andes on A/B |
| yeast ratios without the gate | 1.109 / 1.473 / 0.353 | — | ⚪ |
| protein-level MAE A/B / A/C / B/C | 0.230 / 0.213 / 0.213 | 0.222 / 0.159 / 0.186 | 🔴 1.3x on A/C (moderate) |
| same, without the gate | 0.238 / 0.266 / 0.258 | — | ⚪ |
| proteins evaluated | 10,797 | 8,427 | 🟢 |
| cost, one process, 16 threads | 9.2 min; 17.8 GB | — | ⚪ |

#### Match between runs for reporter ions (PXD014750, hrMS2)

Setup: opt-in `--tmt --mbr` (hrMS2 only; refused for SPS-MS3 and ion-trap MS2). Six files (three
TMT11 plexes of human and yeast peptides at changing mixing, each acquired twice), human + yeast +
entrapment database, `q ≤ 0.01`, `--rescore`, 16 threads; andes `0cfa03cc`. Truth: human-only
peptides must be flat; yeast-only peptides must follow the plex's yeast profile. A value is wrong
when its factor-normalised, impurity-corrected profile deviates by more than 0.5 log2 (RMS); the
excess is the wrong rate of transferred values minus that of the tool's identified values at matched
species and signal. Only cells with a truth are scored, so denominators differ. MaxQuant and
SIMSI-Transfer rows are an earlier analysis of the same files with the same truth and threshold.

| tool | MBR | identified rows | transferred | wrong, transferred | wrong, identified (matched) | excess, points (95% CI) | result |
|---|---|---:|---:|---:|---:|---|---|
| andes `--tmt` | off | 199,730 | 0 | — | 1.96% (all scored rows) | — | ⚪ |
| andes `--tmt --mbr`, all candidates | on | 199,730 | 7,841 (6,260 scored) | 3.90% | 3.50% | +0.40 (−0.04 to +0.89) | 🟡 |
| andes `--tmt --mbr`, `quant_q ≤ 0.01` | on | 199,730 | 3,636 (3,267 scored) | 2.79% | 2.84% | −0.05 (−0.56 to +0.51) | 🟢 well ahead on accuracy |
| MaxQuant 2.7.0.0 | on | 137,539 evidence rows | 7,499 cells | 3.77% | 2.24% | +1.53 (+1.11 to +1.98) | 🔴 |
| MaxQuant 2.7.0.0 | off | 137,539 evidence rows | 0 | — | — | — | ⚪ |
| SIMSI-Transfer 0.8.0 (MaRaCluster 1.04.1), p10 | on | MaxQuant scans | 16,511 cells | 2.25% | 1.85% | +0.40 (+0.18 to +0.64) | 🟡 |

- Cells transferred: andes 3,636 at `quant_q ≤ 0.01` against SIMSI-Transfer 16,511, 🔴 4.5x fewer
  (minor, opt-in). The entrapment share of accepted andes transfers is 1.6% (1.7% of all candidates).
- `--mbr` adds 3.9% rows (1.8% at `quant_q ≤ 0.01`). PSMs through the quantification gate: 131,161
  without it, 134,797 with it. Cost with the flag off 3.4k CPU-s, 8.0 min, 14.6 GB; with it 3.7k
  CPU-s, 8.4 min, 14.6 GB (one run each, one node type).
- The MaxQuant row counts all evidence rows of the types that carry an MS/MS scan, before reverse and
  contaminant hits are removed (isobaric MBR as in MaxQuant ≥ 1.6.12.0: `MULTI-MATCH-MSMS` rows); andes
  counts PSMs at `q ≤ 0.01`.

## 3. How to reproduce

You need `andes` on `PATH`, Docker (for the pinned Percolator image), `curl` and `python3`. Build andes
with `--features thermo` and install the .NET 8 runtime to read `.raw` directly (preferred over
converting; see [`reproduce/README.md`](reproduce/README.md)).

### Standard sets: Astral, TMT, UPS1 (tested scripts)

```bash
DATA=~/andes-bench                            # ~10 GB for these three
cd docs/benchmarks/reproduce
./build_databases.sh "$DATA"                  # UniProt download, ~1 min; prints the UPS1 T/E factor
./fetch_spectra.sh   "$DATA" astral tmt ups1  # 4.8 GB from PRIDE
./run.sh             "$DATA" astral tmt ups1  # search + Percolator; prints wall, PSMs, entrapment hits
```

`run.sh` uses 8 threads and Percolator 3.7.1 `--seed 42 -Y` with the `configs/` mods files; UniProt is
versioned and the original HYE FASTA is gone, so counts differ slightly (the scripts print a sha256
for each input).

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

Count rows at `q-value ≤ 0.01` with `entrap_common.py` (FDP = hits / PSMs × (1 + T/E), T/E from
`build_databases.sh`); it takes about 6 minutes on 8 threads, and the slice width is pinned because it
sets peak memory.

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

For the quick tier search `MouseLiver-Z-T-1` alone, skip pooling and pass `--run MouseLiver-Z-T-1` to
`score_vs_truth.py`; pool all fractions for the deep tier, because one fraction has too few decoys.

### Label-free quantification, ProteoBench (commands; converter tested)

Inputs: the six runs of a module (QExactive `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.raw`
of PXD028735, about 21 GB; Astral `LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw` of
PXD070049, about 15 GB) and `ProteoBenchFASTA_MixedSpecies_HYE.fasta` (31,889 entries). The measured
runs read peak-picked mzML; `R` is a run's stem; unlisted flags are the defaults.

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

6. Score with `proteobench` 0.18.6 (`DDAQuantIonModuleQExactive` or `DDAQuantIonAstralModule`
   `.benchmarking(custom_input.tsv, "Custom", ...)`).

The converter keeps `identified`, `transfer_ms2` and `transfer_ms1` rows; `--column quantity` writes
the raw quantity and `--evidence identified` the identification-only table (tests:
`python3 -I tests/test_proteobench_convert.py`). A one-process variant, `andes --spectrum ... --rescore
--lfq --mbr --quant-dir q/`, was not the configuration measured.

### Truth-design quantification suite (commands; harness tested on synthetic data)

[`quant-eval/`](quant-eval/README.md) needs Python 3 with numpy, pandas, scipy and pyarrow, reads
finished tables and runs no search; write one manifest per set (format in `quant_eval/suite.py`).

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

`X_I` on the board is andes with MBR off and `MQ_noMBR` a second MaxQuant manifest entry; ProteoBench
files come from `python3 -m quant_eval.proteobench` or `proteobench/andes_to_proteobench.py` and
nothing is uploaded.

### TMT match between runs (commands, not a tested script)

PXD014750 hrMS2 (six files, TMT11 impurity sheet, human + yeast + entrapment database); the two runs
differ by one flag, and the truth evaluation reads `run.tmt.parquet` as described in section 4.

```bash
andes --spectrum f1.mzML ... --spectrum f6.mzML --database hy_entrap.fasta --mods mods_tmt11.txt \
      --protocol TMT --rescore --fdr 0.01 --tmt tmt11 --tmt-level 2 --tmt-correction lot.txt \
      --threads 16 --output-tsv run.tsv [--mbr]
```

### TMT knockout (commands, not a tested script)

PXD017803 `Fig2_m23142_hrMS2.raw` and `Fig2_m23145_SPSms3.raw`, reviewed yeast proteins, and a mods
file with Carbamidomethyl C and TMTpro (304.20715) on K and the peptide N-terminus fixed, plus
Oxidation M variable; the measured runs read mzML (both files of each acquisition in one process, 16
threads).

```bash
andes --spectrum Fig2_m23142_hrMS2.raw  --database yeast.fasta --mods mods_tmtpro.txt \
      --precursor-tol 20ppm --protocol TMT --rescore --fdr 0.01 \
      --tmt tmt16 --tmt-level 2 --tmt-correction <kit impurity sheet> --output-tsv hr.tsv
andes --spectrum Fig2_m23145_SPSms3.raw --database yeast.fasta --mods mods_tmtpro.txt \
      --precursor-tol 20ppm --protocol TMT --rescore --fdr 0.01 \
      --tmt tmt16 --tmt-level 3 --tmt-correction <kit impurity sheet> --output-tsv sps.tsv
```

Channel layout and scoring are in section 4.

## 4. Method notes

- **True FDP by entrapment.** FDP = (entrapment hits / accepted) × (1 + T/E), with T/E the ratio of
  searchable target to entrapment peptides, measured for each database. A PSM is an entrapment hit
  only when all its proteins are entrapment and its peptide, with I folded to L, occurs in no target
  protein. `entrap_common.py` holds this rule for `run.sh`, `build_databases.sh` and the glyco
  evaluators. A natural entrapment (a foreign proteome) can be contaminated by conserved homologs;
  the shuffle is the primary reading.
- **Entrapment designs.** Peptide-level shuffles ([`build_shuffled_entrap.py`](build_shuffled_entrap.py);
  each tryptic peptide shuffled with its C-terminal residue fixed, after Wen et al., Nat. Methods
  2025) give every entrapment peptide a target of the same composition and mass. Three realisations
  (seeds 1, 2, 3) for Astral, phospho and glyco; one for UPS1, TMT, TKO and timsTOF. The earlier
  whole-protein shuffle left 89.6% of entrapment peptides on the HYE database without a
  same-composition target, and read 0.81% on Astral where peptide-level databases read 0.92–1.02% and
  *Arabidopsis* (T/E 2.228) 0.99% (#225); on phospho it agreed (1.05% against 0.95–1.04%). FDPs marked
  "whole-protein shuffle" may read low.
- **Conserved-protein filter.** An `ENTRAP_` protein is conserved when its best local alignment
  (BLOSUM62, gap open 11, extend 1, against target proteins preselected by shared 4-mers) covers at
  least 60% of it at 40% identity or more. For UPS1, 269 of 4,531 *E. coli* proteins qualify (719 at
  30% / 50%; 23 at 60% / 60%); the cut was fixed in advance. Hit counts were not re-run at the other
  cuts; T/E without those proteins is 2.872 at 30% / 50% and 2.456 at 60% / 60%.
- **Percolator.** Version 3.7.1, `--seed 42 -Y`, on plain FASTA with andes `XXX_` decoys, for every
  engine on the same host. The q-value has a floor of 1 / (number of top targets), so small or sparse
  runs return few identifications at 1%; pool glyco fractions for this reason.
- **ProteoBench terms.** "Min runs" is the number of the six runs in which a precursor
  (peptidoform and charge) must be quantified. Epsilon is the measured log2 A/B ratio minus the
  expected one; the headline is its median absolute value. CV is the median within-condition
  coefficient of variation, averaged over the two conditions. A transfer q-value scores the peak a
  precursor is matched to in a run where it was not identified; it is not an identification FDR.
- **Truth-design evaluation.** Every tool's table is reduced to raw precursor-by-run intensities, then:
  per-run median normalisation over precursors quantified in all runs, Tukey median polish per
  protein, proteins with ≥ 2 precursors tested with a moderated t, Benjamini–Hochberg 5% per
  condition pair. A true positive is a significant protein of a changing species with the expected
  sign; a false call is a significant protein of a constant species (yeast in UPS1 sets, human in the
  mixtures) or one of a changing species with the wrong sign. DE is scored on the common-tested family
  (rows that andes, MaxQuant and FragPipe all test; PXD030780: andes and FragPipe). Earlier versions
  tested each arm on its own proteins (4,699 against MaxQuant's 4,312 on QE) and counted wrong-sign
  calls as true positives; both are corrected. normalisation and roll-up are the same code for all
  tools. Accuracy is scored on cells every listed tool identifies, against the replicate-based
  expectation. UPS1 detection counts a precursor at a level when one value is within 2x of its own
  high-level reference.
- **Absence classes.** For a precursor, d is the median log2 intensity over runs that contain yeast
  minus that over human-only runs, normalised on human proteins. It is *enriched* when d exceeds a
  boundary learned from human control precursors by the number of human-only runs the precursor is
  seen in (6.67, 5.62, 3.94, 1.07 log2 for 1, 2–4, 5–9, 10 or more runs), *not specific* otherwise,
  *only absent* when seen in human-only runs only, and *no reference* when andes's table of every row
  (2,892,797 rows) lacks it. Enriched-only is a lower bound of species-specific errors.
- **Masked wrong-peak rate.** 5% of identified precursor-run pairs are removed from the library,
  quantified again as transfers, and counted when the apex lies further from the original apex than
  that fit's peak width, at `quant_q ≤ 0.01`.
- **TMT knockout layout.** MET6 (P05694) is knocked out in 126–128N, HIS4 (P00815) in 128C–130N and
  URA2 (P07259) in 130C–132N; 132C–134N is the parental strain. For each PSM of those proteins, divide
  the mean of its knockout channels by the mean of the other channels in `<stem>.tmt.parquet` (rows
  with `quant_gate_passed`), report the median over scans whose reporter window starts below 126 Th,
  and repeat on PSMs with `purity ≥ 0.9`. The MBR evaluation reads `run.tmt.parquet` (`evidence_type`,
  `quant_q`, `quant_gate_passed`, channel columns).
- **Comparisons.** Identification comparisons are runs on the same host with matched parameters and the
  same rescoring; ProteoBench comparisons are public leaderboard submissions. Wall times compare only
  within one host and session; counts are portable.

## 5. Known gaps

| gap | size | status / where to fix |
|---|---|---|
| Label-free quantification rests on two ProteoBench datasets and one Percolator seed (three for error control) | best scores use `quantity_normalized`; raw is `quantity`; Astral CV 0.113 above most submissions; MBR entrapment error measured on QExactive only; leaderboard scored by `proteobench` 0.17.0, andes by 0.18.6; several engines' transfer control unchecked | open |
| Entrapment designs disagree | UPS1: natural 1.95–2.07% against shuffle 0.92–0.94%; timsTOF astrocytes: natural 1.17–1.54% against shuffle 0.94–1.12% on one draw (0.70–0.96% on five more) | the shuffle is the primary reading (section 4); the 500-run study has no entrapment |
| Not re-measured at `2ee76dd6` or `2051616b` | glyco deep tier, Comet phospho FDP, `--chimeric` / `--refine` table, non-default `--predict` settings, phospho agreement with MaxQuant | open |
| Single files | phospho is one file, peptide-level only; `--refine` measured on Astral only | open |
| Glyco composition on fucose-rich tissue | heart and lung: backbone agrees with pGlyco2 94–99%, composition 68–78% (97% on liver) | open |
| Missing tool arms | no MaxQuant run for PXD014415 or PXD030780; FragPipe has no MBR-off run; MaxQuant MBR-on rows unfiltered except in the filtered-matches table | open |
| Not benchmarked | iTRAQ, phosphosite localisation, timsTOF quantification; Java MS-GF+ not re-run under current defaults (counts historical) | open |

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
