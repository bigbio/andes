# Identification benchmarks

¹ single run: one search with one Percolator seed, or one tool output. andes ranges span seeds,
draws or runs. Comet on the timsTOF runs is a mean over six draws; the run count behind the Java
MS-GF+ UPS1 count is not recorded.

Engines on the benchmark VM (8 threads, `/usr/bin/time`, plain databases, Percolator seed 42): andes
`2051616b` (2026-10-23), Comet 2025.01 (2026-10-05; phospho 2026-10-06). **‡** UPS1 andes is the
`2ee76dd6` run below (cluster node, entrapment database), so its cost does not compare with Comet's.

<!-- BEGIN generated: engines-cost (docs/benchmarks/render.py) -->
| dataset | engine | PSMs @ q≤0.01 | wall | CPU time | peak memory |
|---|---|---:|---:|---:|---:|
| Astral | andes | 52,169 | 114–117 s | 569–636 s | 4.8 GB |
| | Comet 2025.01 | 31,435¹ | 217 s¹ | 1,594 s¹ | 8.1 GB¹ |
| TMT a05058 | andes | 13,516 | 40–47 s | 200–216 s | 1.7 GB |
| | Comet 2025.01 | 10,504¹ | 77 s¹ | 568 s¹ | 2.9 GB¹ |
| UPS1 | andes | 16,567–16,675 ‡ | 13–15 s ‡ | 71–79 s ‡ | 1.4 GB ‡ |
| | Comet 2025.01 | 14,734¹ | 42 s¹ | 309 s¹ | 2.9 GB¹ |
| Phospho (PXD007653) | andes | 38,623 | 337 s | 2,160 s | 7.9 GB |
| | Comet 2025.01 | 33,984¹ | 1,096 s¹ | 8,399 s¹ | 6.0 GB¹ |
<!-- END generated: engines-cost -->

![PSMs at q ≤ 0.01 per engine](https://andes.quantms.org/figures/psms.svg)

![Wall time, CPU time and peak memory, andes as a percentage of Comet](https://andes.quantms.org/figures/cost_relative.svg)

![Entrapment-estimated true FDP at nominal 1%](https://andes.quantms.org/figures/true_fdp.svg)

Figures: drawn by the site from `website/benchmarks.json`, which `render.py` writes from
[`results/`](results/).

## Astral HCD (PXD070049)

Setup: one 15-minute DDA run of the ProteoBench human/yeast/E. coli sample, HYE FASTA (31,889
entries); andes `2051616b`, 2026-10-23; Comet 2025.01, same VM, matched parameters, same
Percolator, 2026-10-05.

| metric | andes | Comet 2025.01 | verdict | note |
|---|---:|---:|---|---|
| PSMs at `q ≤ 0.01` | 52,169–52,437 | 31,435¹ | 🟢 well ahead | +66% |
| peptides at `q ≤ 0.01` | 31,147–31,238 | — | ⚪ not compared | |
| true FDP at nominal 1%, peptide-level entrapment (T/E 1.075; 9 database × seed runs) | 0.89–1.08% | — | 🟢 on target | |
| PSMs / peptides at a true FDP of 1% | 47,267–47,977 / 28,141–28,857 | — | ⚪ not compared | |
| true FDP, *Arabidopsis* entrapment (T/E 2.228) | 0.99–1.05% (49,965–50,185 PSMs) | — | 🟢 on target | |
| wall time; CPU time; memory | 114–117 s; 569–636 CPU-s; 4.8 GB | 217 s¹ | 🟢 well ahead | wall 1.9x faster |

## TMT-labelled, low-res (PXD007683, run a05058)

Setup: one ion-trap CID run of a TMT 10-plex sample; human + yeast reviewed proteins (26,483
sequences); andes `2051616b`; Comet 2025.01 as above.

| metric | andes | Comet 2025.01 | verdict | note |
|---|---:|---:|---|---|
| PSMs at `q ≤ 0.01` | 13,443–13,516 | 10,504¹ | 🟢 ahead | +29% |
| peptides at `q ≤ 0.01` | 11,852–11,922 | — | ⚪ not compared | |
| true FDP, peptide-level entrapment (T/E 1.073) | 0.67–0.76% (13,325–13,371 PSMs at 1%) | — | 🟢 on target | |
| true FDP, *Arabidopsis* entrapment (T/E 2.021) | 0.83–1.05% (13,161–13,388 PSMs at 1%) | — | 🟢 on target | |
| wall time; CPU time; memory | 40–47 s; 200–216 CPU-s; 1.7 GB | 77 s; 568 CPU-s; 2.9 GB¹ | 🟢 well ahead | wall 1.6–1.9x faster |

## TMT-labelled, high-res (PXD017803, TKO hrMS2)

Setup: `Fig2_m23142_hrMS2` (TMTpro, Orbitrap MS2), `--protocol TMT`, reviewed yeast + cRAP with
peptide-level entrapment (T/E 1.039); andes `2051616b`.

| metric | andes | other tool | verdict | note |
|---|---:|---|---|---|
| PSMs at `q ≤ 0.01` | 6,628–6,686 | — | ⚪ not compared | |
| true FDP at nominal 1% | 0.98–1.31% | — | 🟡 near nominal | up to 0.3 points above nominal |
| PSMs at a true FDP of 1% | 6,607–6,662 | — | ⚪ not compared | |

## Low-res LFQ (PXD001819, UPS1_5000amol_R1)

Setup: one run of the UPS1 standard in yeast; andes `2ee76dd6`, 2026-11-07, documented defaults
(2 missed cleavages, default since #300); Percolator seeds 42, 1 and 7; all-proteins entrapment
rule; two entrapment designs, each its own database. Java MS-GF+ v20240326 is from an earlier
session and was not re-run.

| metric | andes | Java MS-GF+ | verdict | note |
|---|---:|---:|---|---|
| PSMs at `q ≤ 0.01` (shuffle / natural database) | 16,567–16,623 / 16,632–16,675 | 15,904 | 🟡 on par | +4% to +5% |
| peptides (shuffle / natural) | 4,176–4,200 / 4,116–4,227 | — | ⚪ not compared | |
| true FDP, peptide-level shuffle (T/E 1.033; 75–77 hits) | 0.92–0.94% | — | 🟢 on target | |
| PSMs at a true FDP of 1%, shuffle | 16,617–16,733 | — | ⚪ not compared | |
| true FDP, natural *E. coli* (4,531 `ENTRAP_` proteins, T/E 2.451; 94–100 hits) | 1.95–2.07% | — | ⚪ not compared | context, not a calibration reading (see first bullet) |
| PSMs at a true FDP of 1%, natural | 14,786–15,045 | — | ⚪ not compared | |
| natural, hits on conserved proteins removed (49–51 hits) | 1.01–1.06% (1.07–1.11% with T/E 2.630); 16,505–16,607 PSMs at 1% | — | 🟡 near nominal | crosses 1.1% with T/E 2.630 |
| andes wall; CPU; memory | 13–15 s; 71–79 CPU-s; 1.4 GB | — | ⚪ not compared | |
| Percolator, each of 3 seeds | 38–40 s; about 190–200 CPU-s; 1.2 GB | — | ⚪ not compared | |

- The natural reading is not a clean calibration. Of the 94 hits at seed 42, 45 sit on 8
  conserved groups (EF-Tu 32 PSMs from two peptides, ornithine carbamoyltransferase 7, enolase 5).
  UPS1 is made in *E. coli*, so host proteins may really be in the sample; the data cannot
  separate that from a mis-assigned yeast peptide. The shuffle (0.92–0.94%) is the primary reading.
- The any-protein rule (a PSM is entrapment if any of its proteins is) reads 3.09–3.23% on the
  same runs because shared peptides become false hits; the all-proteins rule is used here.
- T/E 2.451 is for the andes digest (2 missed cleavages). With `--max-missed-cleavages 1` (T/E
  2.351) the same build gives 16,238–16,367 PSMs, natural 1.94–1.99% and shuffle 0.88–0.96%
  (shuffle T/E 1.039); other digests give T/E 2.370–2.550.

## Phospho-enriched (mouse liver, EasyPhos, single file)

Setup: `20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw` as mzML, not in any model's training
data; 80 Da index slices; mouse reviewed proteins plus the whole-protein shuffled entrapment the
Comet row used; andes `2051616b`; Comet 2025.01, 2026-10-06.

| metric | andes | Comet 2025.01 | verdict | note |
|---|---:|---:|---|---|
| PSMs at `q ≤ 0.01` | 38,596–38,629 (27,693–27,724 phospho-bearing) | 33,984¹ | 🟢 ahead | +14% |
| peptides at `q ≤ 0.01` | 20,422–20,461 | — | ⚪ not compared | |
| true FDP, 3 peptide-level shuffles (T/E 1.061) | 0.80–0.99% (38,604–38,843 PSMs at 1%) | 1.77%¹ (whole-protein shuffle) | 🟢 on target | Comet above nominal, on a whole-protein shuffle that can read low |
| true FDP, *Arabidopsis* entrapment (T/E 1.403) | 1.22–1.36% (38,298–38,412 PSMs at 1%) | — | 🔴 above nominal | |
| wall time; CPU time; memory | 337 s; 2,160 CPU-s; 7.9 GB | — | ⚪ not compared | |
| agreement with depositors' MaxQuant IDs (PEP ≤ 0.01, 23,563 scans) | 82.4% covered, 95.3% same-scan backbone | — | ⚪ not compared | agreement with the depositors' identifications |

- Phosphosite localisation is not scored.
- Comet's 1.77% is on the whole-protein shuffle, which can read low.
- The MaxQuant agreement was measured 2026-10-07 at `40774aca` and not repeated.

## timsTOF DDA-PASEF (4 runs)

Setup: native `.d`, 20 ppm precursor, 2 missed cleavages, `--top-n 1`, 80 Da slices, Carbamidomethyl
C fixed, Oxidation M and protein N-terminal Acetyl variable; reviewed human proteins with a shuffle of
the human targets (T/E 1.050) or reviewed yeast (natural, T/E 3.784) as entrapment; andes
`2ee76dd6`, Percolator seeds 42, 1 and 7 for the single-draw tables, 2026-11-07 to 2026-11-09. Comet 2025.01 on tdf2mzml mzML, 20 ppm, isotope error −1 to 3, `fragment_bin_tol 0.02`,
internal reversed decoys; both tools on the same six shuffle draws (seeds 1, 11–15), Percolator seed 42, 2026-11-09; means over draws.

| metric (4 runs) | andes | Comet 2025.01 | verdict | note |
|---|---:|---:|---|---|
| PSMs at `q ≤ 0.01` | 15.6k–58.6k | 14.1k–50.8k | 🟢 ahead | +11% to +17% |
| peptide sequences | 3,865–24,146 | 3,396–22,007 | 🟢 ahead | +10% to +14% |
| PSMs at a true FDP of 1% | 15.9k–58.8k | 14.0k–49.5k | 🟢 ahead | +14% to +21% |
| true FDP at nominal 1% | 0.79–0.91% | 1.14–1.31% (23 of 24 draws above 1%) | 🟢 on target | Comet above nominal |
| PSMs on the same mzML | — | — | 🟡 on par | +8% to +16% (native `.d` adds 1.2–2.6%) |
| wall time | 54–62 s | 159–222 s | 🟢 well ahead | 3x faster |
| memory | 2.3–2.7 GB | 6.0–6.4 GB | 🟢 well ahead | 2.2–2.6x less |

Per-run tables (both tools; andes on one draw and on the natural design): [ADDITIONAL_BENCHMARKS.md](ADDITIONAL_BENCHMARKS.md#timstof-per-run).

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

## Glyco (mouse liver, PXD005553)

Setup: `--glyco --glyco-species mouse --decoy-strategy sequon`, mouse proteins plus entrapment.
References for the same spectra are in [`glyco/truth/`](glyco/truth/): pGlyco2 (depositors) and
MSFragger-Glyco (reanalysis, PXD031032). "Confirmed" is the share of a reference's spectra that
andes identifies with the same peptide. Quick tier: fraction `MouseLiver-Z-T-1` as mzML, three
peptide-level shuffled databases × Percolator seeds 42, 1 and 2. Deep tier: all five fractions
pooled before Percolator, five seeds, measured 2026-09-05 at `14818d3e` (older, smaller glycan list,
mzML), not re-measured.

| metric | andes | reference | verdict | note |
|---|---:|---|---|---|
| glycoPSMs at 1% (quick tier) | 6,992–7,072 | — | ⚪ not compared | |
| true FDP (quick tier) | 0.61–1.35% | — | 🟡 near nominal | 20–45 hits per run, mostly counting noise |
| spectra confirmed, quick tier | 88.2–88.4% | pGlyco2 | 🔴 behind | 11.6–11.8% missed |
| spectra confirmed, quick tier | 89.3–89.7% | MSFragger-Glyco | 🔴 behind | 10.3–10.7% missed |
| pGlyco2 spectra that lose the per-scan choice, deep tier | about 20% | pGlyco2 | 🔴 behind | about 20% missed |
| time; memory (quick tier, 8 threads) | 5.8–6.6 min; 5.4–5.6 GB | — | ⚪ not compared | |

- The open problem is selection, not generation: almost every reference spectrum yields a
  candidate row.
- Reading the `.raw` natively, as the [reproduce commands](reproduce.md#glyco-commands-scripts-in-glyco) do, can give somewhat different counts.
