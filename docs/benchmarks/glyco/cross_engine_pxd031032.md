# Cross-engine glycopeptide comparison — MouseLiver-Z-T-1

Date: 2026-10-21 · dataset: PXD031032 (reanalysis of PXD005553)

*Re-run 2026-10-21 after andes's Y-ion generation was aligned with GlycoDecipher's
single-glycosidic-bond-cleavage convention and the glycan tree topology was retained
(`pGlyco-N-Mouse.gdb`, 6,660 topology trees instead of 1,833 compositions). The five
reference-engine inputs are unchanged from the 2026-09-11 run, so every engine-to-engine
Jaccard below is byte-identical to that run — only the `andes` column moved.*

andes is compared against five published engines on the same raw file
`MouseLiver-Z-T-1.raw` (mouse liver N-glycoproteomics). Every engine is
normalised to a canonical glycan composition `(HexNAc, Hex, Fuc, NeuAc, NeuGc)`
and a bare peptide sequence (flanking residues, modification masses and the
glycan tag stripped). A "glycopeptide" is the distinct `(peptide, composition)`
pair; PSM rows are the accepted spectra behind those pairs.

Reproduce with `docs/benchmarks/glyco/compare_engines.py` (`--help` lists the inputs). The
result files are the ones deposited in PXD031032 and PXD005553.

## Engines and filters

| Engine | Result file | Filter |
|---|---|---|
| andes | `<stem>.glyco.pin` + native-rescore `.psms` | q ≤ 0.01, seed 42, single-cleavage Y-ions, topology-retained `pGlyco-N-Mouse.gdb` (6,660 trees) |
| Byonic | result workbook (`.xlsx`), PXD031032 | Score ≥ 300 |
| Glyco-Decipher | GPSM text export, PXD031032 | PeptideFDR ≤ 0.01 |
| StrucGP | result workbook (`.xlsx`), PXD031032 | all rows |
| pGlyco2 | `MouseLiver-Z-T-1-FDR.txt` (PXD005553) | TotalFDR ≤ 0.01, target only |
| MSFragger-Glyco | glyco `psm.tsv`, PXD031032 | Expectation ≤ 0.01, Delta Mass → composition |

## Identifications

| Engine | glycopeptides | PSM rows | proteins |
|---|---|---|---|
| **andes** | **2842** | 6874 | 449 (474 incl. entrapment) |
| Byonic | 2509 | 5997 | 401 |
| MSFragger-Glyco | 2507 | 6229 | 422 |
| Glyco-Decipher | 2301 | 5947 | 388 |
| StrucGP | 1789 | 3803 | 430 |
| pGlyco2 | 1761 | 3995 | 340 |

⚠ **The andes headline excludes entrapment hits.** 31 of its 6,905 accepted rows map to
entrapment proteins, which cannot be present in the sample and are therefore known-false.
No other engine here searched an entrapment database, so andes is the only row carrying
them. `compare_engines.py` excludes entrapment from the headline count (2,842) and prints
the inclusive figure (2,869) on a second line; the ranking against Byonic (2,509) does
not change either way.

31 entrapment rows in 6,905 is, separately, a reasonable error floor: it is the measured
cost of the 1% threshold, not a red flag.

## Pairwise overlap on (peptide, glycan composition)

| | Glyco-Decipher | StrucGP | Byonic | pGlyco2 | MSFragger-Glyco | andes |
|---|---|---|---|---|---|---|
| Glyco-Decipher | — | 0.372 | 0.538 | 0.378 | 0.527 | **0.538** |
| StrucGP | | — | 0.335 | 0.272 | 0.313 | 0.323 |
| Byonic | | | — | 0.352 | 0.551 | **0.541** |
| pGlyco2 | | | | — | 0.314 | 0.320 |
| MSFragger-Glyco | | | | | — | **0.550** |

Shared counts: andes–MSFragger-Glyco 1898, andes–Byonic 1878, andes–Glyco-Decipher
1799, andes–StrucGP 1130, andes–pGlyco2 1117. Protein-level: andes shares
375 / 357 / 362 / 273 / 324 with MSFragger-Glyco / Byonic / Glyco-Decipher / StrucGP / pGlyco2 respectively.

**Read:** andes identifies the most glycopeptides (2842), and its highest
agreement is with MSFragger-Glyco (0.550) and Byonic (0.541), followed closely
by Glyco-Decipher (0.538) — it covers 1898 of MSFragger-Glyco's 2507 (76%), 1878 of
Byonic's 2509 (75%) and 1799 of Glyco-Decipher's 2301 (78%). pGlyco2 and StrucGP
are the outliers at the low end on every row (0.27–0.38), consistent
with their smaller glycan databases and older scoring, not an andes-specific
divergence.

## Peptidoform agreement — pooled 5 fractions

The set-level table above is single-fraction `MouseLiver-Z-T-1`, because only that
fraction's result files for Glyco-Decipher/StrucGP/Byonic are on disk. MSFragger-Glyco and pGlyco2 have full
5-fraction references committed under `truth/`, so against those two the comparison is
pooled (all five fractions → one Percolator run), and it is scan-matched and
peptidoform-level: a "glycopeptidoform" is the distinct `(backbone, glycan composition)`
pair, and the agreement figure is the share of spectra *both* engines identified where they
call the **same backbone and the same glycan composition**.

| (pooled, 1% FDR) | pGlyco2 | MSFragger-Glyco |
|---|---:|---:|
| reference glycopeptidoforms | 2,286 | 1,785 |
| andes glycopeptidoforms | 4,723 | 4,723 |
| co-identified scans | 15,937 | 13,122 |
| **same-scan peptidoform** | **98.0%** | **96.1%** |
| same-scan backbone | 99.6% | 99.9% |
| peptide agreement (backbones found anywhere) | 94.0% | 95.3% |

The 2.0% (pGlyco2) / 3.9% (MSFragger-Glyco) of co-identified scans that disagree on the
peptidoform are backbone-correct but glycan-composition differences — the same defect class
the 2026-09-11 note flagged. Peptide-level agreement (94–95% of the reference's backbones
also found by andes) is high, and where both engines identify a scan the backbone agrees
99.6–99.9% of the time; the glycan composition is the residual disagreement.

The corresponding stage-decomposition on the same pooled run (`score_vs_truth.py`, bare
backbone match, no glycan) attributes every reference spectrum to one bucket: **pGlyco2**
CONFIRMED 88.9% / WRONG_TARGET_WON 4.9% / DECOY_WON 6.0% / NOT_EMITTED 0.1% · **MSFragger-Glyco**
CONFIRMED 89.7% / WRONG_TARGET_WON 4.6% / DECOY_WON 5.7% / NOT_EMITTED 0.0%. Selection
loss (wrong target + decoy) is 10.9% (pGlyco2) / 10.3% (MSFragger-Glyco).

## Timing — A/B against `main`

`MouseLiver-Z-T-1` (45,905 MS2 spectra; the `.mgf` and native `.raw` read the same
45,905 scans), 112 threads, `sequon-reverse` decoy, `mouse_entrap.fasta` (the
whole-protein shuffled entrapment used before #225). Three arms, all on the same host and thread count, so
the comparison is apples-to-apples:

1. **`main`** — the pre-PR code (`bb1ccdf`), built-in composition enumerator,
   default peptide-first path.
2. **branch default** — this branch's default peptide-first path
   (`--glyco-glycan-gdb` with the bundled mouse `.gdb`, no `--glyco-full-glycan-db`).
3. **`--glyco-full-glycan-db`** — the mass-driven path (skips the b/y postings
   index).

| metric | `main` (bb1ccdf) | branch default | `--glyco-full-glycan-db` |
|---|---:|---:|---:|
| total wall-clock | 3,388 s | 3,471 s | **135 s** |
| glyco PSM rows | 42,106 | 38,992 | 43,452 |
| targets @ 1% FDR (2× rule) | 4,997 | 5,485 | 5,325 |
| RawScore AUC (target > decoy) | 0.6533 | 0.6629 | 0.6465 |

**The honest speed baseline is `main`, not the branch default.** `main` finishes
the file in 3,388 s; `--glyco-full-glycan-db` finishes it in 135 s — **~25×
faster**. The previous table quoted "117×/77×" against the branch's own default
path, which is the wrong baseline: on this host the two peptide-first paths are
within ~2% of each other (3,388 s vs 3,471 s), so the "slower while using four
times the cores" gap the review cited was a cross-machine, cross-thread artifact,
not a property of the branch.

**Rows are not identifications.** `--glyco-full-glycan-db` emits 43,452 candidate
rows — +3.2% over `main`'s 42,106, +11.4% over the branch default's 38,992 — yet
identifies 5,325 targets at 1% FDR, which is *above* `main`'s 4,997 (+6.6%), not
below it. The extra rows are candidate emissions, not accepted PSMs; RawScore
separation moves 0.6533 → 0.6465. The old table's "−1.1% targets" compared
against the branch default (5,485), not against `main`.

**"targets @ 1% FDR (2× rule)"** is decoy counting, not entrapment: sort PSMs by
RawScore descending, `FDR = 2D/(D+T)` on the `sequon-reverse` `XXX_` decoys. The
factor 2 is exact because target:decoy is 1:1 by construction (see
`reproduce/README.md`); it is not an assumed entrapment scaling factor. Computed
by `compare_preperc.py`, which is committed alongside the other benchmark scripts.

**Single-fraction Percolator caveat.** Percolator's 3-fold cross-validation does
not converge on one fraction of the full-glycan-db output (43 k rows; one fold has
no separable training direction), so the single-fraction metric above is the
decoy-counting 2× rule. The pooled 5-fraction Percolator run is the comparable
q-value measurement (`eval_yield.py`), per the "pool before Percolator" rule.

## Golden fixtures

The two in-repo golden PINs move against `main` in one direction: `glyco.pin` keeps 120 of
121 rows and gains none, `glyco_highres.pin` keeps 117 of 121 and gains none. These are tiny
fixtures built for path coverage rather than yield, and the full-file numbers above go the
other way, but the direction has been consistent across every revision of this work and is
recorded here rather than left for someone to rediscover.

## Notes and caveats

- Filter thresholds are each engine's own convention (Byonic Score, Glyco-Decipher
  PeptideFDR, pGlyco2 TotalFDR, MSFragger-Glyco Expectation, andes native q-value); they
  are not directly calibrated to one another.
- MSFragger-Glyco's export carries the glycan only as a "Delta Mass" (residue-sum
  mass, no water subtracted), so the script maps it back to a composition via
  the bundled mouse `.gdb` within 0.02 Da (with ^13C M+1/M+2 precursor recovery).
  The peptide-level Expectation ≤ 0.01 is the FDR proxy — the file has no glycan
  q-value and no decoys. The Hex+NeuAc ↔ Fuc+NeuGc isobar (0.27 mDa) is resolved
  to the nearer mass; 3.1% of MSFragger-Glyco rows sit on that isobar.
- andes uses the topology-retained `pGlyco-N-Mouse.gdb` (6,660 distinct trees, 2,325
  compositions observed at 1% FDR). Byonic searched a wider glycan space (its 698
  distinct compositions include antenna variants), so exact composition matches are the
  conservative subset.
- pGlyco2 encodes the glycosylated Asn as `J` in its peptide column
  (`VSQVLHEGGHJVTK` → `…HNVTK`); the script maps `J→N`.
- Byonic peptide cells carry flanking residues **and** decimal masses
  (`K.HLLEN[+1864.634]ATASVSEAER.K`); the script strips flanks by first/last dot,
  not a full split, because the mass annotation contains a decimal point.
- Protein accessions are normalised to the canonical UniProt form (isoform suffix
  `-N` stripped), so `Q9Z2G6-2` and `Q9Z2G6` count as one protein.
