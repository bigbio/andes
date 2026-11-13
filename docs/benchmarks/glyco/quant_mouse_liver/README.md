# Glycopeptide LFQ + match-between-runs benchmark — mouse liver, 5 fractions

This benchmark exercises the intact-glycopeptide quantification path
(`andes --glyco --lfq --mbr`) on the five mouse-liver HCD fractions of
PXD005553 (deep tier: all fractions searched together, one pooled PIN). It
answers two questions:

1. **Identification** — does the search still produce a controlled set of
   glycopeptides on this data? (the existing harness already checks this; repeated
   here so the quantification has a fixed denominator)
2. **Quantification** — once identification is controlled, do the LFQ areas and
   the MBR transfers behave like measurements: complete matrices, consistent
   replicates, stable glycoform shares, stable protein roll-ups?

No external quantitative reference exists for these files: the pGlyco2 and
MSFragger-Glyco tables in `../truth/` carry no intensities, and pGlycoQuant
requires pGlyco3 result files that are not available here. Quantitative quality
is therefore assessed by internal consistency, not by agreement with another
engine's areas.

## Result in brief

| | value |
|---|---|
| glyco-PSMs at 1% FDR (pooled Percolator) | **33,823** |
| unique glycopeptides at 1% (topology-resolved) | **6,509** |
| entrapment FDP at 1% (sequon-corrected) | **0.93%** (95% CI 0.79–1.10) |
| of those PSMs quantified by `--lfq` | **97.5%** |
| ID-gated quantified features (identified + transferred) | **40,416** (32,967 + 7,449) |
| glycopeptides at composition / topology resolution | **4,661 / 6,392** |
| matrix fill, identified-only → +MBR | **58.6% → 84.6%** |
| glycopeptides in ≥3 runs, identified → +MBR | **2,466 → 4,150** |
| cross-run Pearson (log2 area), identified → +MBR | **0.889 → 0.827** |
| glycoform-share mean \|Δ\|, identified → +MBR | **0.075 → 0.063** |
| protein Top3 cross-run r, identified → +MBR | **0.940 → 0.911** |

**The most important result is a measurement-semantics one.** The `Accepted`
flag in the emitted tables is a *feature/measurement* q-value (XIC target vs its
+11 Th decoy twin), not an identification FDR. Used as shipped, the tables
report 254,240 quantified features whose entrapment share is **42.7%**
(sequon-corrected FDP ≈ **91.1%**). Joining each feature to the Percolator
q-value of its PSM (the `.lfq.tsv` `SpecId` matches the PIN byte for byte; a
transferred row carries `{donor SpecId}_mbr_{run}`) leaves 40,416 features
(15.9%) with a 0.85% entrapment share (FDP 1.82%; 0.88% for identified
features alone). Every quantification number below is computed **after** that
join. This is the documented workflow, but the emitted `Accepted` column alone
is easy to misread, and the gap between the two readings is large enough to
dominate any conclusion drawn from the tables.

## Data and search parameters

| | |
|---|---|
| spectra | `MouseLiver-Z-T-1.raw` … `MouseLiver-Z-T-5.raw` (PXD005553, HCD, Q Exactive) |
| database | `mouse_entrap.fasta`: 17,229 mouse + 17,229 same-length shuffled `ENTRAP_` proteins |
| decoys | `--decoy-strategy sequon` (shuffle, then restore the target sequon count; this replaced the older `sequon-reverse` name and recipe upstream) |
| glycan DB | `pGlyco-N-Mouse.gdb`, 6,660 entries with structure (topology) preserved |
| search | trypsin, ≤3 missed cleavages, 6–50 aa, 20 ppm, charge 2–5, oxidation M |
| quant | MS1, 10 ppm, ±60 s, ≥4 isotopes, cosine ≥ 0.7, feature FDR 0.05 |
| runtime | ≈14 min end-to-end (PIN phase 796 s); 112 threads; candidate index cached |
| andes | 0.4.0 (branch `glyco/lfq-mbr` @ 0866a63, on `origin/main` f25ca32) |

Outputs are written next to the PIN. Because `--output-pin` was given a name
already ending in `.glyco.pin`, andes appends its own suffix, so the files are
`mouseliver5.glyco.glyco.lfq.tsv`, `…glyco.pin`, and so on.

## 1. Identification layer (context)

Pooled Percolator 3.7.1 (`--seed 42 -Y --only-psms=false`) on the 211,392-row
PIN:

| gate | target PSMs | entrapment hits | sequon-corrected FDP |
|---|---:|---:|---:|
| q ≤ 0.01 | 33,823 | 148 | **0.93%** |
| q ≤ 0.05 | 35,747 | 690 | 4.11% |

The search-space correction is 2.13× (sample 50,425 N-X-S/T sequons vs
entrapment 44,593). Agreement with the two reference engines on the same
spectra:

| reference | spectra | recovered at andes 1% |
|---|---:|---:|
| pGlyco2 (depositors') | 17,855 | **89.1%** (15,912) |
| MSFragger-Glyco | 14,626 | **89.9%** (13,152) |

`eval_yield.py`: 33,823 glyco-PSMs, 6,509 unique glycopeptides, 2,306
compositions, 1,223 glycosites.

## 2. What the quantification tables contain

`Accepted` = feature q ≤ 0.05. Per run, before any identification gate:

| run | identified rows | transferred rows | accepted identified | accepted transferred |
|---|---:|---:|---:|---:|
| T-1 | 23,839 | 45,859 | 21,243 | 28,427 |
| T-2 | 21,801 | 47,586 | 19,387 | 29,765 |
| T-3 | 23,115 | 46,223 | 20,554 | 32,378 |
| T-4 | 23,293 | 46,314 | 20,532 | 31,215 |
| T-5 | 23,581 | 46,046 | 20,737 | 30,002 |
| **all** | **115,629** | **232,028** | **102,453** | **151,787** |

That is 254,240 accepted features over 60,790 (peptide, site, composition)
glycopeptides. The emitted tables agree with each other: `.glyco.matrix.tsv` has
60,790 rows and 225,633 filled cells, exactly the rows of `.glyco.quant.tsv`;
`site.tsv` glycoform shares sum to 1.0 within 1e-5 per (run, peptide, site) at
every one of 132,753 sites; and no protein has Top3 > summed area.

Their entrapment share and the effect of the identification gate:

![Identification gate](fig1_identification_gate.png)

| feature set | features | entrapment | sequon-corrected FDP |
|---|---:|---:|---:|
| accepted, no ID gate | 254,240 | 42.7% | 91.1% |
| accepted + ID q ≤ 0.01 | 40,416 | 0.85% | 1.82% |
| accepted + ID q > 0.05 / no q | 209,773 | 51.0% | 108.6% |

By evidence class, at the 1% gate: identified 136/32,967 entrap (0.41%),
transferred 209/7,449 (2.81%). Transfers are discussed in §6.

Gated feature counts per run (identified + transferred): 8,252 / 8,016 / 8,171 /
7,945 / 8,032 — flat across fractions, as replicates should be. Quantification
recovers 32,967 of the 33,823 accepted PSMs (97.5%); the remaining 2.5% did not
produce an accepted feature (XIC cosine/envelope gates or no usable MS1 scan).

## 3. Matrix completeness

4,661 glycopeptides × 5 runs. Without MBR, 13,664 cells are filled (58.6%);
MBR fills 6,060 more, for 19,724 (84.6%). The number of *identified*
glycopeptides does not change (4,661); MBR only adds measurements of existing
identifications.

![Matrix completeness](fig2_matrix_completeness.png)

| presence | identified | identified + MBR |
|---|---:|---:|
| ≥1 run | 4,661 | 4,661 |
| ≥2 runs | 3,134 | 4,454 |
| ≥3 runs | 2,466 | 4,150 |
| ≥4 runs | 1,962 | 3,646 |
| ≥5 runs | 1,441 | 2,813 |

## 4. Cross-run consistency

Pearson r on log2 area, pairwise over glycopeptides measured in both runs:

| | identified | identified + MBR |
|---|---:|---:|
| median pairwise r | **0.889** | **0.827** |
| range | 0.853–0.921 | 0.797–0.861 |

![Cross-run consistency](fig3_cross_run.png)

The drop is the expected MBR trade-off and is confined to transferred cells.
Splitting each pair by how the two measurements were obtained:

| class | median r | range | n per pair |
|---|---:|---|---:|
| both identified | **0.896** | 0.862–0.926 | 1,838–2,106 |
| one transferred | **0.555** | 0.499–0.660 | 880–1,091 |
| both transferred | **0.620** | 0.551–0.702 | 437–636 |

Glycoform shares agree better after MBR, because MBR adds the missing glycoforms
of a site in runs that identified only part of the family: mean |Δshare| 0.0752
→ **0.0631** (19,716 → 34,944 site-glycan pairs).

## 5. Protein roll-up

657 protein groups have at least one ID-gated glycopeptide. Summed glycopeptide
area vs Top3 (three most intense glycopeptides) is internally consistent, and
Top3 reproduces across runs:

- median cross-run Top3 Pearson r: **0.940** (identified) → **0.911** (+MBR).

Only 657 proteins is expected at this FDR with the 17,229-protein mouse
database; the ungated tables inflate this to 13,498, another symptom of
quantifying unfiltered rank-1 matches.

## 6. Transfers and their quality

- 232,028 transfers attempted across all runs; 175,379 measured; 151,787 at
  feature q ≤ 0.05; **7,449** also pass the source-identification gate.
- Transfer feature q: median 0.0174, p90 0.0337, max 0.0499 — the transfer
  competition uses its own target/decoy class, as designed.
- 7,449 transfers come from **3,346 distinct donor PSMs** (mean 2.23 target runs
  per donor).
- The donor subset that is actually transferable is **enriched in entrapment:
  2.18% of donor PSMs (73/3,346) vs 0.41% of identified features (136/32,967)**.
  Every feature-q bin shows the same effect (3.5% vs 0.41% in 0.010–0.015;
  1.5% vs 0.49% in 0.015–0.020). Transferability selects for abundant
  precursors, and false rank-1 matches are enriched among abundant features;
  the source q-value alone therefore does not make an MBR-only cell as
  trustworthy as a directly identified one.
- Consequence in the numbers: cross-run r for "one transferred" is ~0.56 vs
  ~0.90 for "both identified", and the gated transferred features carry 2.8%
  entrapment. MBR is a coverage tool, not an identification tool; the transfer
  q already separates the two classes, and downstream analysis should keep them
  distinguishable.

## 7. Glycan topology

The glycan database retains 6,660 tree topologies. Quantification is inherently
composition-level (two topologies of one composition are one intact precursor
mass), while the search, scoring and output keep the tree:

| | compositions | topologies | compositions with >1 topology |
|---|---:|---:|---:|
| all accepted (no gate) | 60,790 | 64,070 | 2,439 |
| ID-gated | 4,661 | 6,392 | 1,131 |

At the identification gate, **24% of quantified compositions (1,131/4,661)
retain more than one topology** — the isomer ambiguity is carried through to the
output rather than collapsed away, which is the point of the topology-retention
change. `eval_yield.py` counts 6,509 topology-resolved glycopeptides at 1%;
6,392 survive quantification. The two counters key the peptide slightly
differently (the PIN counter folds modified residues and reads the full tag), so
read this as agreement in scale, not an exact identity.

## 8. Limitations

- **No external quantitative reference.** Every quantitative claim here is
  internal (replicate consistency, completeness, share stability). It shows the
  path is self-consistent, not that it matches pGlycoQuant's areas.
- **Two q-value layers.** The emitted `Accepted` is a measurement q; the
  identification q must be joined from Percolator. The tables do not encode
  both. A future CLI pass that filters `.glyco.lfq.tsv` by an optional
  `--psms` argument would remove this trap; until then the join is mandatory.
- **Entrapment estimator is linear** (raw share × (1 + S/E)); at the extreme
  ungated values it extrapolates past 100%. It is accurate near the 1% gate,
  where it is used for the verdicts.
- **FQR is not implemented** (pGlycoQuant's two-component GMM over the +10 Da
  decoy). The feature q here is a picked target/decoy q-value against the
  +11 Th twin, and the transfer q is a separate class.
- **MBR RT model** is piecewise-linear per run with one knot per 300 anchors
  (14,560 anchors over 5 maps) and a narrow XIC window; it is designed for
  replicate/fraction matching, not for large RT shifts.
- **Numbers move with the search baseline.** This run is on `origin/main`
  (f25ca32) with `--decoy-strategy sequon`; an earlier version of this file
  was measured on the `glyco/topology-wurcs` tree with the retired
  `sequon-reverse` decoy recipe, which produces different decoy sequences.
  Compare runs within one baseline, not across decoy strategies.

## 9. Reproduction

The analysis scripts need `numpy`, `pandas`, `matplotlib` and `scipy`.

```bash
# 0. Thermo .raw input needs the `thermo` feature
cargo build --release -p andes --bin andes --features thermo

# 1. search + quantify (5 fractions -> one PIN + lfq/matrix tables)
DOTNET_ROOT=/home/daichengxin/.dotnet \
/home/daichengxin/gitrepo/andes/target/release/andes \
  --spectrum /Data/daicx/glyco-mouse/MouseLiver-Z-T-{1,2,3,4,5}.raw \
  --database /home/daichengxin/andes-bench/databases/mouse_entrap.fasta \
  --decoy-strategy sequon \
  --glyco --glyco-glycan-gdb /home/daichengxin/gitrepo/pro/pGlyco-N-Mouse.gdb \
  --fragmentation HCD --fragment-tol-ppm 20 \
  --lfq --mbr --threads 112 \
  --output-pin /tmp/glyco_quant_bench/mouseliver5.glyco.pin

# 2. pooled Percolator (OUT is the directory holding the PIN)
OUT=/tmp/glyco_quant_bench
docker run --rm -v "$OUT":/r quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 \
  percolator --seed 42 -Y --only-psms=false \
  --results-psms /r/out.psms --decoy-results-psms /r/out.dpsms \
  /r/mouseliver5.glyco.glyco.pin

# 3. identification checks (run from $OUT; H = this directory)
H=/home/daichengxin/gitrepo/andes/docs/benchmarks/glyco
cd "$OUT"
python3 "$H/eval_yield.py"  mouseliver5.glyco.glyco.pin out.psms
python3 "$H/eval_entrap.py" mouseliver5.glyco.glyco.pin out.psms 0.01 \
  /home/daichengxin/andes-bench/databases/mouse_entrap.fasta
python3 "$H/score_vs_truth.py" "$H/truth/pglyco2_mouse_liver.tsv.gz" \
  mouseliver5.glyco.glyco.pin out.psms 0.01
python3 "$H/score_vs_truth.py" "$H/truth/msfragger_glyco_mouse_liver.tsv.gz" \
  mouseliver5.glyco.glyco.pin out.psms 0.01

# 4. quantification analysis + figures
python3 "$H/quant_mouse_liver/glyco_quant_bench.py" --dir "$OUT" \
  --stem mouseliver5.glyco --out "$OUT/report" \
  --psms "$OUT/out.psms" \
  --fasta /home/daichengxin/andes-bench/databases/mouse_entrap.fasta
python3 "$H/quant_mouse_liver/glyco_quant_figs.py" \
  --json "$OUT/report/glyco_quant_numbers.json" --out "$OUT/report"
```

## Files

| file | content |
|---|---|
| `run_glyco_quant_bench.sh` | the andes invocation used here |
| `glyco_quant_bench.py` | analysis: gate join, matrix, consistency, proteins, entrapment |
| `glyco_quant_figs.py` | the three figures |
| `glyco_quant_report.md` | full generated tables (every run pair, every evidence split) |
| `glyco_quant_numbers.json` | machine-readable results |
| `fig1_identification_gate.png` | what the tables report vs the identification gate |
| `fig2_matrix_completeness.png` | glycopeptides reaching ≥k runs, identified vs +MBR |
| `fig3_cross_run.png` | 5×5 Pearson log2-area heatmaps |
