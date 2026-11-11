# Glycopeptide LFQ + match-between-runs benchmark — mouse liver, 5 fractions

This benchmark exercises the intact-glycopeptide quantification path
(`andes --glyco --lfq --glyco-mbr`) on the five mouse-liver HCD fractions of
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
| glyco-PSMs at 1% FDR (pooled Percolator) | **33,625** |
| unique glycopeptides at 1% (topology-resolved) | **6,476** |
| entrapment FDP at 1% (sequon-corrected) | **0.92%** (95% CI 0.78–1.08) |
| of those PSMs quantified by `--lfq` | **97.4%** |
| ID-gated quantified features (identified + transferred) | **40,430** (32,765 + 7,665) |
| glycopeptides at composition / topology resolution | **4,622 / 6,359** |
| matrix fill, identified-only → +MBR | **58.6% → 85.5%** |
| glycopeptides in ≥3 runs, identified → +MBR | **2,457 → 4,121** |
| cross-run Pearson (log2 area), identified → +MBR | **0.888 → 0.830** |
| glycoform-share mean \|Δ\|, identified → +MBR | **0.076 → 0.063** |
| protein Top3 cross-run r, identified → +MBR | **0.942 → 0.918** |

**The most important result is a measurement-semantics one.** The `Accepted`
flag in the emitted tables is a *feature/measurement* q-value (XIC target vs its
+11 Th decoy twin), not an identification FDR. Used as shipped, the tables
report 250,804 quantified features whose entrapment share is **42.7%**
(sequon-corrected FDP ≈ **90.9%**). Joining each feature to the Percolator
q-value of its PSM (the `.lfq.tsv` `SpecId` matches the PIN byte for byte; a
transferred row carries `{donor SpecId}_mbr_{run}`) leaves 40,430 features
(16.1%) with a 0.97% entrapment share (FDP 2.08%; 0.88% for identified
features alone). Every quantification number below is computed **after** that
join. This is the documented workflow, but the emitted `Accepted` column alone
is easy to misread, and the gap between the two readings is large enough to
dominate any conclusion drawn from the tables.

## Data and search parameters

| | |
|---|---|
| spectra | `MouseLiver-Z-T-1.raw` … `MouseLiver-Z-T-5.raw` (PXD005553, HCD, Q Exactive) |
| database | `mouse_entrap.fasta`: 17,229 mouse + 17,229 same-length shuffled `ENTRAP_` proteins |
| decoys | `--decoy-strategy sequon-reverse` (sequon density preserved) |
| glycan DB | `pGlyco-N-Mouse.gdb`, 6,660 entries with structure (topology) preserved |
| search | trypsin, ≤3 missed cleavages, 6–50 aa, 20 ppm, charge 2–5, oxidation M |
| quant | MS1, 10 ppm, ±60 s, ≥4 isotopes, cosine ≥ 0.7, feature FDR 0.05 |
| runtime | ≈11 min end-to-end (PIN phase 666 s); 112 threads; candidate index cached |
| andes | 0.4.0 (working tree `glyco/topology-wurcs`) |

Outputs are written next to the PIN. Because `--output-pin` was given a name
already ending in `.glyco.pin`, andes appends its own suffix, so the files are
`mouseliver5.glyco.glyco.lfq.tsv`, `…glyco.pin`, and so on.

## 1. Identification layer (context)

Pooled Percolator 3.7.1 (`--seed 42 -Y --only-psms=false`) on the 211,994-row
PIN:

| gate | target PSMs | entrapment hits | sequon-corrected FDP |
|---|---:|---:|---:|
| q ≤ 0.01 | 33,625 | 145 | **0.92%** |
| q ≤ 0.05 | 35,455 | 681 | 4.09% |

The search-space correction is 2.13× (sample 50,425 N-X-S/T sequons vs
entrapment 44,593). Agreement with the two reference engines on the same
spectra:

| reference | spectra | recovered at andes 1% |
|---|---:|---:|
| pGlyco2 (depositors') | 17,855 | **88.9%** (15,873) |
| MSFragger-Glyco | 14,626 | **89.7%** (13,114) |

`eval_yield.py`: 33,625 glyco-PSMs, 6,476 unique glycopeptides, 2,329
compositions, 1,209 glycosites.

## 2. What the quantification tables contain

`Accepted` = feature q ≤ 0.05. Per run, before any identification gate:

| run | identified rows | transferred rows | accepted identified | accepted transferred |
|---|---:|---:|---:|---:|
| T-1 | 23,472 | 44,756 | 20,869 | 27,874 |
| T-2 | 21,395 | 46,515 | 19,046 | 29,923 |
| T-3 | 22,673 | 45,157 | 20,094 | 31,928 |
| T-4 | 22,775 | 45,361 | 19,989 | 30,536 |
| T-5 | 23,232 | 44,946 | 20,364 | 30,181 |
| **all** | **113,547** | **226,735** | **100,362** | **150,442** |

That is 250,804 accepted features over 59,337 (peptide, site, composition)
glycopeptides. The emitted tables agree with each other: `.glyco.matrix.tsv` has
59,337 rows and 222,242 filled cells, exactly the rows of `.glyco.quant.tsv`;
`site.tsv` glycoform shares sum to 1.0 per (run, peptide, site) at every one of
130,946 sites; and no protein has Top3 > summed area.

Their entrapment share and the effect of the identification gate:

![Identification gate](fig1_identification_gate.png)

| feature set | features | entrapment | sequon-corrected FDP |
|---|---:|---:|---:|
| accepted, no ID gate | 250,804 | 42.7% | 90.9% |
| accepted + ID q ≤ 0.01 | 40,430 | 0.97% | 2.08% |
| accepted + ID q > 0.05 / no q | 206,478 | 51.0% | 108.6% |

By evidence class, at the 1% gate: identified 135/32,765 entrap (0.41%),
transferred 259/7,665 (3.38%). Transfers are discussed in §6.

Gated feature counts per run (identified + transferred): 8,191 / 8,130 / 8,122 /
7,939 / 8,048 — flat across fractions, as replicates should be. Quantification
recovers 32,765 of the 33,625 accepted PSMs (97.4%); the remaining 2.6% did not
produce an accepted feature (XIC cosine/envelope gates or no usable MS1 scan).

## 3. Matrix completeness

4,622 glycopeptides × 5 runs. Without MBR, 13,551 cells are filled (58.6%);
MBR fills 6,197 more, for 19,748 (85.5%). The number of *identified*
glycopeptides does not change (4,622); MBR only adds measurements of existing
identifications.

![Matrix completeness](fig2_matrix_completeness.png)

| presence | identified | identified + MBR |
|---|---:|---:|
| ≥1 run | 4,622 | 4,622 |
| ≥2 runs | 3,097 | 4,434 |
| ≥3 runs | 2,457 | 4,121 |
| ≥4 runs | 1,945 | 3,695 |
| ≥5 runs | 1,430 | 2,876 |

## 4. Cross-run consistency

Pearson r on log2 area, pairwise over glycopeptides measured in both runs:

| | identified | identified + MBR |
|---|---:|---:|
| median pairwise r | **0.888** | **0.830** |
| range | 0.850–0.923 | 0.795–0.859 |

![Cross-run consistency](fig3_cross_run.png)

The drop is the expected MBR trade-off and is confined to transferred cells.
Splitting each pair by how the two measurements were obtained:

| class | median r | range | n per pair |
|---|---:|---|---:|
| both identified | **0.896** | 0.858–0.928 | 1,830–2,086 |
| one transferred | **0.552** | 0.496–0.654 | 893–1,119 |
| both transferred | **0.603** | 0.532–0.667 | 445–709 |

Glycoform shares agree better after MBR, because MBR adds the missing glycoforms
of a site in runs that identified only part of the family: mean |Δshare| 0.0758
→ **0.0630** (19,566 → 35,265 site-glycan pairs).

## 5. Protein roll-up

647 protein groups have at least one ID-gated glycopeptide. Summed glycopeptide
area vs Top3 (three most intense glycopeptides) is internally consistent, and
Top3 reproduces across runs:

- median cross-run Top3 Pearson r: **0.942** (identified) → **0.918** (+MBR).

Only 647 proteins is expected at this FDR with the 17,229-protein mouse
database; the ungated tables inflate this to 13,382, another symptom of
quantifying unfiltered rank-1 matches.

## 6. Transfers and their quality

- 226,735 transfers attempted across all runs; 173,542 measured; 150,442 at
  feature q ≤ 0.05; **7,665** also pass the source-identification gate.
- Transfer feature q: median 0.0167, p90 0.0321, max 0.0481 — the transfer
  competition uses its own target/decoy class, as designed.
- 7,665 transfers come from **3,402 distinct donor PSMs** (mean 2.25 target runs
  per donor).
- The donor subset that is actually transferable is **enriched in entrapment:
  2.70% of donor PSMs (92/3,402) vs 0.41% of identified features (135/32,765)**.
  Every q bin shows the same effect (in the 0.008–0.01 bin, 15.3% of transfers
  vs 6.6% of identified features are entrapment). Transferability selects for
  abundant precursors, and false rank-1 matches are enriched among abundant
  features; the source q-value alone therefore does not make an MBR-only cell
  as trustworthy as a directly identified one.
- Consequence in the numbers: cross-run r for "one transferred" is ~0.55 vs
  ~0.90 for "both identified", and the gated transferred features carry 3.4%
  entrapment. MBR is a coverage tool, not an identification tool; the transfer
  q already separates the two classes, and downstream analysis should keep them
  distinguishable.

## 7. Glycan topology

The glycan database retains 6,660 tree topologies. Quantification is inherently
composition-level (two topologies of one composition are one intact precursor
mass), while the search, scoring and output keep the tree:

| | compositions | topologies | compositions with >1 topology |
|---|---:|---:|---:|
| all accepted (no gate) | 59,337 | 62,618 | 2,436 |
| ID-gated | 4,622 | 6,359 | 1,124 |

At the identification gate, **24% of quantified compositions (1,124/4,622)
retain more than one topology** — the isomer ambiguity is carried through to the
output rather than collapsed away, which is the point of the topology-retention
change. `eval_yield.py` counts 6,476 topology-resolved glycopeptides at 1%;
6,359 survive quantification. The two counters key the peptide slightly
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
  (14,264 anchors over 5 maps) and a narrow XIC window; it is designed for
  replicate/fraction matching, not for large RT shifts.

## 9. Reproduction

```bash
# 1. search + quantify (5 fractions -> one PIN + lfq/matrix tables)
DOTNET_ROOT=/home/daichengxin/.dotnet \
/home/daichengxin/gitrepo/andes/target/release/andes \
  --spectrum /Data/daicx/glyco-mouse/MouseLiver-Z-T-{1,2,3,4,5}.raw \
  --database /home/daichengxin/andes-bench/databases/mouse_entrap.fasta \
  --decoy-strategy sequon-reverse \
  --glyco --glyco-glycan-gdb /home/daichengxin/gitrepo/pro/pGlyco-N-Mouse.gdb \
  --fragmentation HCD --fragment-tol-ppm 20 \
  --lfq --glyco-mbr --threads 112 \
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
