# Cross-engine glycopeptide comparison — MouseLiver-Z-T-1

Date: 2026-09-11 · dataset: PXD031032 (reanalysis of PXD005553)

andes is compared against five published engines on the same raw file
`MouseLiver-Z-T-1.raw` (mouse liver N-glycoproteomics). Every engine is
normalised to a canonical glycan composition `(HexNAc, Hex, Fuc, NeuAc, NeuGc)`
and a bare peptide sequence (flanking residues, modification masses and the
glycan tag stripped). A "glycopeptide" is the distinct `(peptide, composition)`
pair; PSM rows are the accepted spectra behind those pairs.

Reproduce with `docs/benchmarks/glyco/compare_engines.py`.

## Engines and filters

| Engine | Result file | Filter |
|---|---|---|
| andes | `<stem>.glyco.pin` + native-rescore `.psms` | q ≤ 0.01, seed 42, `--glyco-full-glycan-db`, `pGlyco-N-Mouse.gdb` (1833 comps) |
| Byonic | `MouseLiver-Z-T-1.raw_20210216_Byonic.xlsx` | Score ≥ 300 |
| Glyco-Decipher | `MouseLiver-Z-T-1_Glyco_Decipher_GPSM.txt` | PeptideFDR ≤ 0.01 |
| StrucGP | `MouseLiver-Z-T-1_result_StrucGP.xlsx` | all rows |
| pGlyco 2.0 | `MouseLiver-Z-T-1-FDR.txt` (PXD005553) | TotalFDR ≤ 0.01, target only |
| MSFragger-Glyco | `Mouse_MSFragger-Glyco_GlyTouCan_psm.tsv` | Expectation ≤ 0.01, Delta Mass → composition |

## Identifications

| Engine | glycopeptides | PSM rows | proteins |
|---|---|---|---|
| **andes** | **2868** | 6913 | 479 (459 excl. entrapment) |
| Byonic | 2509 | 5997 | 401 |
| MSFragger-Glyco | 2507 | 6229 | 422 |
| Glyco-Decipher | 2301 | 5947 | 388 |
| StrucGP | 1789 | 3803 | 430 |
| pGlyco 2.0 | 1761 | 3995 | 340 |

andes reports 24 entrapment-mapped rows at 1% FDR (of 459 non-entrapment
proteins) — a minor, expected FDR floor from the entrapment database, not a
red flag.

## Pairwise overlap on (peptide, glycan composition)

| | Glyco-Decipher | StrucGP | Byonic | pGlyco2 | MSFragger | andes |
|---|---|---|---|---|---|---|
| Glyco-Decipher | — | 0.372 | 0.538 | 0.378 | 0.527 | **0.545** |
| StrucGP | | — | 0.335 | 0.272 | 0.313 | 0.326 |
| Byonic | | | — | 0.352 | 0.551 | **0.551** |
| pGlyco2 | | | | — | 0.314 | 0.322 |
| MSFragger | | | | | — | **0.558** |

Shared counts: andes–MSFragger 1924, andes–Byonic 1911, andes–Glyco-Decipher
1823, andes–StrucGP 1145, andes–pGlyco2 1128. Protein-level: andes shares
371 / 366 / 357 / 327 / 274 with MSFragger / Glyco-Decipher / Byonic / pGlyco2 /
StrucGP respectively.

**Read:** andes identifies the most glycopeptides (2868), and its highest
agreement is with MSFragger-Glyco (0.558) and Byonic (0.551), followed closely
by Glyco-Decipher (0.545) — it covers 1924 of MSFragger's 2507 (77%), 1911 of
Byonic's 2509 (76%) and 1823 of Glyco-Decipher's 2301 (79%). pGlyco 2.0 (2017)
and StrucGP are the outliers at the low end on every row (0.27–0.38), consistent
with their smaller glycan databases and older scoring, not an andes-specific
divergence.

## Timing — A/B against `main`

`MouseLiver-Z-T-1` (45,905 MS2 spectra; the `.mgf` and native `.raw` read the same
45,905 scans), 112 threads, `sequon-reverse` decoy, `mouse_entrap.fasta` (1:1
shuffled-self entrapment). Three arms, all on the same host and thread count, so
the comparison is apples-to-apples:

1. **`main`** — the pre-PR code (`bb1ccdf`), built-in composition enumerator,
   default peptide-first path.
2. **branch default** — this branch's default peptide-first path
   (`--glyco-glycan-gdb pGlyco-N-Mouse.gdb`, no `--glyco-full-glycan-db`).
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

## Notes and caveats

- Filter thresholds are each engine's own convention (Byonic Score, Glyco-Decipher
  PeptideFDR, pGlyco2 TotalFDR, MSFragger Expectation, andes native q-value); they
  are not directly calibrated to one another.
- MSFragger-Glyco's export carries the glycan only as a "Delta Mass" (residue-sum
  mass, no water subtracted), so the script maps it back to a composition via
  `pGlyco-N-Mouse.gdb` within 0.02 Da (with ^13C M+1/M+2 precursor recovery).
  The peptide-level Expectation ≤ 0.01 is the FDR proxy — the file has no glycan
  q-value and no decoys. The Hex+NeuAc ↔ Fuc+NeuGc isobar (0.27 mDa) is resolved
  to the nearer mass; 3.1% of MSFragger rows sit on that isobar.
- andes uses `pGlyco-N-Mouse.gdb` (1833 compositions). Byonic searched a wider
  glycan space (its 698 distinct compositions include antenna variants), so exact
  composition matches are the conservative subset.
- pGlyco 2.0 encodes the glycosylated Asn as `J` in its peptide column
  (`VSQVLHEGGHJVTK` → `…HNVTK`); the script maps `J→N`.
- Byonic peptide cells carry flanking residues **and** decimal masses
  (`K.HLLEN[+1864.634]ATASVSEAER.K`); the script strips flanks by first/last dot,
  not a full split, because the mass annotation contains a decimal point.
- Protein accessions are normalised to the canonical UniProt form (isoform suffix
  `-N` stripped), so `Q9Z2G6-2` and `Q9Z2G6` count as one protein.
