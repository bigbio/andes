# andes benchmarks

andes is benchmarked on public PRIDE datasets for identification (high-res, TMT-labelled and
low-res LFQ, phospho, glyco) and for quantification (label-free with match between runs, and
TMT). Identifications are counted after Percolator at `q ≤ 0.01`. Where the database carries an
entrapment component, the page also gives the **true** false discovery proportion (FDP), which
is what tells you whether the nominal 1% holds. Other tools are named with the version compared.
Their numbers come either from runs on the same host with the same rescoring, or from a public
leaderboard. Each page gives the date and the andes commit behind each row.

| page | what it holds |
|---|---|
| this page | summary of every benchmark, opt-in modes |
| [Identification](identification.md) | Astral, TMT (low- and high-res), UPS1, phospho, timsTOF, glyco PXD005553; cost on one host |
| [Glycoproteomics, PXD065680](glyco/pxd065680.md) | five published glyco tools on one file |
| [Label-free quantification](lfq.md) | ProteoBench, eight truth designs, a 500-run study |
| [Glycopeptide LFQ, PXD005553](glyco/quant_mouse_liver/README.md) | label-free quantification and match between runs of glycopeptides |
| [TMT quantification](tmt.md) | knockout interference, known ratios, match between runs for reporter ions |
| [Methods](methods.md) | measurement conditions, entrapment rules, evaluation rules |
| [Known gaps](gaps.md) | what is not measured or not re-measured |
| [How to reproduce](reproduce.md) | scripts and commands |
| [Supplementary tables](ADDITIONAL_BENCHMARKS.md) | per-run, per-level and full leaderboard tables |

## Summary

Builds: identification `2051616b` (2026-10-23; timsTOF and entrapment rows `2ee76dd6`, 2026-11-07/09),
label-free quantification `09e1c00d` (2026-11-05), glycopeptide LFQ `0f036477` (2026-11-11), TMT
`0cfa03cc`; Percolator 3.7.1 `-Y`.

One row per metric, andes against each tool on the same data. Verdicts: 🟢 well ahead / ahead ·
🟡 on par · 🔴 behind / well behind · ⚪ not compared; for error control 🟢 on target · 🟡 near
nominal · 🔴 above nominal. The rule (by the ratio between the two values: under 1.1x on par, 1.5x
or more "well") is in [Methods](methods.md#verdicts); the note column gives the size. Details and
settings are on the page of each domain. ¹ single run: one search with one Percolator seed, or one
tool output; andes ranges span seeds, draws or runs. Unmarked values of other tools are means over
several draws (Comet on timsTOF) or have no recorded run count (Java MS-GF+ on UPS1, the TMT
match-between-runs rows).

**Identification** (PSMs at Percolator `q ≤ 0.01`, 8 threads)

<!-- BEGIN generated: identification-summary (docs/benchmarks/render.py) -->
| data | metric | andes | other tool | verdict | note |
|---|---|---:|---|---|---|
| Astral HCD, PXD070049 | PSMs | 52,169 | Comet 31,435¹ | 🟢 well ahead | +66% |
| Astral HCD | wall time | 114–117 s | Comet 217 s¹ | 🟢 well ahead | 1.9x faster |
| Astral HCD | true FDP at 1% | 0.89–1.08% | — | 🟢 on target | |
| timsTOF DDA-PASEF, 4 runs | PSMs | 15.6k–58.6k | Comet 14.1k–50.8k | 🟢 ahead | +11% to +17% per run |
| timsTOF DDA-PASEF | PSMs at a true FDP of 1% | 15.9k–58.8k | Comet 14.0k–49.5k | 🟢 ahead | +14% to +21% per run |
| timsTOF DDA-PASEF | true FDP at 1% | 0.79–0.91% | Comet 1.14–1.31% | 🟢 on target | Comet above nominal |
| timsTOF DDA-PASEF | wall time | 54–62 s | Comet 159–222 s | 🟢 well ahead | 3x faster |
| TMT low-res, PXD007683 | PSMs | 13,516 | Comet 10,504¹ | 🟢 ahead | +29% |
| Phospho, PXD007653 | PSMs | 38,596–38,629 | Comet 33,984¹ (at 1.77% true FDP) | 🟢 ahead | +14%; andes at 0.80–0.99% true FDP |
| Low-res UPS1, PXD001819 | PSMs | 16,567–16,675 | Java MS-GF+ 15,904 | 🟡 on par | +4% to +5% |
| Low-res UPS1 | true FDP at 1% | 0.92–0.94% | — | 🟢 on target | |
| TMT high-res, PXD017803 | PSMs | 6,628–6,686 | — | ⚪ not compared | no other tool run |
| Glyco, PXD005553 | reference glycoPSMs confirmed | 79–90% | pGlyco2, MSFragger-Glyco (references) | 🔴 behind | 10–21% of reference spectra missed |
<!-- END generated: identification-summary -->

**Quantification** (label-free at `09e1c00d`, TMT at `0cfa03cc`; DE counts are the common-tested
family, sign-aware, MBR on)

<!-- BEGIN generated: quantification-summary (docs/benchmarks/render.py) -->
| data | metric | andes | MaxQuant | FragPipe | verdict | note |
|---|---|---:|---:|---:|---|---|
| ProteoBench QE | median \|ε\| (lower is better) | 0.161 | 0.210¹ | 0.200¹ | 🟢 ahead of both | 23% / 20% lower |
| ProteoBench QE | precursors in ≥ 3 runs | 63,796 | 44,299¹ | 63,701¹ | 🟢 ahead of MaxQuant · 🟡 on par with FragPipe | +44% / same |
| ProteoBench QE | CV | 0.117 | 0.219¹ | 0.156¹ | 🟢 well ahead of MaxQuant · 🟢 ahead of FragPipe | 47% / 25% lower |
| ProteoBench Astral | median \|ε\| | 0.186 | 0.160¹ | 0.212¹ | 🔴 behind MaxQuant · 🟢 ahead of FragPipe | 1.16x MaxQuant at 55% of andes' depth / 12% lower |
| ProteoBench Astral | precursors in ≥ 3 runs | 36,675 | 20,283¹ | 39,057¹ | 🟢 well ahead of MaxQuant · 🟡 on par with FragPipe | +81% / 6% fewer |
| ProteoBench Astral | CV | 0.112 | 0.092¹ | 0.127¹ | 🔴 behind MaxQuant · 🟢 ahead of FragPipe | 1.2x MaxQuant at 55% of andes' depth / 12% lower |
| Truth designs, 8 sets | precursors quantified, each tool at its own search settings | 10.7k–93.6k | | | 🟢 well ahead on 1 set · 🟢 ahead on 6 · 🟡 on par on 1 | +19% to +51% over the best other tool |
| PXD009740 / PXD009815 / QE / Astral | precursors quantified, andes at matched search settings (no `--chimeric`, no `--refine`) | 27,226 / 29,600 / 81,987 / 51,798 | 29,293 / 27,459 / 51,609 / 24,745¹ | 25,888 / 28,286 / 71,377 / 46,347¹ | 🟡 on par / 🟡 on par / 🟢 ahead / 🟢 ahead | −7% / +5% / +15% / +12% over the best other tool |
| PXD009740 / PXD009815 / QE / Astral | DE true positives | 1,300 / 844 / 1,047 / 685 | 1,304 / 844 / 857 / 650¹ | 1,237 / 849 / 946 / 689¹ | 🟡 on par / 🟡 on par / 🟢 ahead of both / 🟡 on par | QE +22% / +11% |
| PXD009740 / PXD009815 / QE / Astral | DE false calls, no fold-change cut | 305 / 153 / 110 / 296 | 129 / 92 / 37 / 210¹ | 302 / 301 / 43 / 389¹ | MaxQuant: 🔴 well behind / 🔴 well behind / 🔴 well behind / 🔴 behind · FragPipe: 🟡 on par / 🟢 well ahead / 🔴 well behind / 🟢 ahead | 1.4–3.0x MaxQuant; mostly real small shifts, see note |
| PXD009740 / PXD009815 / QE | DE false calls at \|log2 FC\| ≥ 0.5 (seed 42) | 167 / 39 / 11 | 74 / 36 / 18¹ | 259 / 102 / 20¹ | MaxQuant: 🔴 well behind / 🟡 on par / 🟢 well ahead · FragPipe: 🟢 well ahead on all three | one andes seed |
| PXD009740 / PXD009815 / QE / Astral | DE true positives / false calls, MaxQuant matches filtered at q ≤ 0.01 | 1,156 / 763 / 983 / 568; 272 / 151 / 105 / 297 | 1,077 / 725 / 851 / 518; 131 / 62 / 42 / 154¹ | | true positives: 🟡 on par / 🟡 on par / 🟢 ahead / 🟡 on par · false calls: 🔴 well behind on all four | +5% to +15% true positives; 1.9–2.5x false calls |
| PXD009740 / PXD009815 / QE / Astral | precision of unchanged proteins (within-condition spread, tool ÷ andes) | 1 | 1.15 / 1.55 / 1.26¹ / — | 1.40 / 1.30 / 1.25¹ / — | MaxQuant: 🟢 ahead / 🟢 well ahead / 🟢 ahead · FragPipe: 🟢 ahead on all three | andes tighter than both |
| UPS1 lowest 3 levels, PXD009740 / PXD009815 | precursors within 2x of truth | 502 / 24 | 417 / 23¹ | 458 / 50¹ | PXD009740: 🟢 ahead of MaxQuant · 🟡 on par with FragPipe; PXD009815: 🟡 on par with MaxQuant · 🔴 well behind FragPipe | PXD009740 at FragPipe's row count: 437 against 458 |
| Truth designs, all sets | bright identified cells off by > 2x | 1.8–7.2% | 0.3–1.6%¹ | 0.3–3.8%¹ | MaxQuant: 🔴 well behind on all five sets · FragPipe: 🔴 well behind on five, 🟢 well ahead on PXD001819 | 3–17x MaxQuant per set |
| Truth designs, 6 sets | protein CV | 4.4–10.1% | 7.0–11.6%¹ | 4.9–13.5%¹ | MaxQuant: 🟢 ahead on four, 🟢 well ahead on PXD009815 · FragPipe: 🟢 well ahead on PXD009740, 🟢 ahead on two, 🟡 on par on two, 🔴 behind on Astral | |
| Known-absent yeast, PXD014415 | false-transfer rate | 0.20% (0.08% at `quant_q ≤ 0.01`) | | 0.09%¹ | 🔴 well behind FragPipe at default · 🟡 on par at `quant_q ≤ 0.01` | 2x at default |
| Blank run, PXD016921 | rows identified in the blank | 0 on main of 2026-11-09, after #333 `e841aae0` (85 at `09e1c00d`) | 20¹ | 8¹ | 🟢 well ahead of MaxQuant · 🟡 on par with FragPipe | counts below 10 are on par; see the blank-run notes |
| Blank run, PXD016921 | rows transferred into the blank | 0, not attempted | 85¹ unfiltered (9 at `Match q-value ≤ 0.01`) | 0, not scored | ⚪ not tested | andes and FragPipe make no transfer decision |
| Glycopeptide LFQ, PXD005553 (5 fractions, `0f036477`) | matrix fill, identified-only → +MBR | 58.6% → 85.1% (held-out 10%: 83.7% → 84.7%; the MBR scorer is trained on these fractions) | | | ⚪ not compared | no other tool run |
| TMT knockout, PXD017803 hrMS2 | interference (lower is better) | 0.192 (0.000 with `--tmt-unmix on`) | 0.156¹ | | 🔴 behind MaxQuant at default · 🟢 well ahead with `--tmt-unmix on` | 1.23x MaxQuant |
| TMT known ratios, PXD007683 | PSMs | 149,121 | 97,107¹ | | 🟢 well ahead | +54% |
| TMT known ratios, PXD007683 | protein ratio error | 0.213 | 0.159¹ | | 🔴 behind | 1.3x on A/C |
| TMT match between runs, PXD014750 | excess wrong-value rate of transfers | −0.05 points | +1.53 | | 🟢 ahead of MaxQuant · 🟡 on par with SIMSI-Transfer | 95% intervals: andes −0.56 to +0.51, MaxQuant +1.11 to +1.98, SIMSI-Transfer +0.40 (+0.18 to +0.64) |
| TMT match between runs | cells transferred | 3,636 | 7,499 | | 🔴 well behind MaxQuant and SIMSI-Transfer | SIMSI-Transfer 16,511; opt-in mode |
| Large study, PXD042233 (500 runs) | precursors in ≥ 3 runs | 243,895 | | | ⚪ not compared | no other tool run |
<!-- END generated: quantification-summary -->

ProteoBench rows use the MaxQuant and FragPipe leaderboard submissions with the most precursors
(MaxQuant 1.5.3.30 on QE, 2.6.7.0 on Astral; FragPipe 23.0); every other submission is in [Label-free quantification](lfq.md#proteobench).

Note on DE false calls: andes' protein values are tighter than the other tools', so small real
shifts of the "unchanged" species reach significance. On PXD009740, 296 of 320 false calls (seed
42 on main of 2026-11-09, after #333 `e841aae0`; 305 in the table, mean over seeds at `09e1c00d`) are against its highest UPS1
level; without that level andes has 955 true positives and 24 false calls
against MaxQuant's 953 and 31. Analysis in [Label-free quantification, *Differential expression*](lfq.md#differential-expression-bh-5-per-pair).

Settings, error control and cost of every row: [Identification](identification.md), [Label-free quantification](lfq.md), [TMT quantification](tmt.md) and [Methods](methods.md).

**Opt-in modes** (Astral, PSMs / distinct peptides, Percolator seeds 42, 1 and 2; measured
2026-10-05 to 2026-10-13 on earlier builds and not re-measured at `2051616b`):

<!-- BEGIN generated: opt-in (docs/benchmarks/render.py) -->
| mode | andes | without the mode |
|---|---|---|
| `--chimeric` (second and later peptides per spectrum), at an entrapment true FDP of 1% | 66,377–66,424 / 29,384–29,400 | 43,076–43,145 / 26,723–26,746 |
| `--refine` (PTM discovery, high-res only) | 48,312 PSMs, 28,728 peptides at `q ≤ 0.01`; its second-pass PSMs at ~1.2–1.8% FDP | 47,080 / 28,932 |
<!-- END generated: opt-in -->

![Chimeric mode on Astral, PSMs at a true FDP of 1%](https://andes.quantms.org/figures/chimeric.svg)

Figures: drawn by the site from `website/benchmarks.json`, which `render.py` writes from
[`results/`](results/).

## Layout

    docs/benchmarks/
      README.md      this page: summary and contents
      identification.md, lfq.md, tmt.md  per-domain results
      methods.md, gaps.md, reproduce.md  conditions and rules, known gaps, commands
      results/       headline numbers (JSON); render.py writes them into the pages and website/benchmarks.json
      check_links.py link and anchor check for these pages
      ADDITIONAL_BENCHMARKS.md  supplementary tables (per-run, per-level, full leaderboards, andes-only diagnostics)
      reproduce/     database, download and run scripts (no hardcoded paths)
      glyco/         glyco harness: pooling, yield, entrapment, reference scoring
      glyco/pxd065680.md  PXD065680 five-tool comparison (source of the site glyco page)
      glyco/quant_mouse_liver/  glycopeptide LFQ + MBR benchmark, PXD005553 (5 fractions)
      configs/       per-engine parameter files
      proteobench/   andes MBR output -> ProteoBench "Custom" table (python3 + pyarrow)
      entrap_common.py  shared entrapment rules (ALL-ENTRAP_, I/L rescue, T/E)
      tests/         tests: python3 -I tests/test_entrap.py (stdlib), tests/test_proteobench_convert.py

Bulk spectra and third-party engine binaries stay outside git.
