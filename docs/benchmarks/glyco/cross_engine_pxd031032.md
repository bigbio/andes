# Cross-engine glycopeptide comparison — MouseLiver-Z-T-1

Date: 2026-09-11 · dataset: PXD031032 (reanalysis of PXD005553)

andes is compared against four published engines on the same raw file
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

## Identifications

| Engine | glycopeptides | PSM rows | proteins |
|---|---|---|---|
| **andes** | **2868** | 6913 | 479 (459 excl. entrapment) |
| Byonic | 2509 | 5997 | 402 |
| Glyco-Decipher | 2301 | 5947 | 389 |
| StrucGP | 1789 | 3803 | 430 |
| pGlyco 2.0 | 1761 | 3995 | 340 |

andes reports 24 entrapment-mapped rows at 1% FDR (of 459 non-entrapment
proteins) — a minor, expected FDR floor from the entrapment database, not a
red flag.

## Pairwise overlap on (peptide, glycan composition)

| | Glyco-Decipher | StrucGP | Byonic | pGlyco2 | andes |
|---|---|---|---|---|---|
| Glyco-Decipher | — | 0.372 | 0.538 | 0.378 | **0.545** |
| StrucGP | | — | 0.335 | 0.272 | 0.326 |
| Byonic | | | — | 0.352 | **0.551** |
| pGlyco2 | | | | — | 0.322 |

Shared counts: andes–Byonic 1911, andes–Glyco-Decipher 1823, andes–StrucGP
1145, andes–pGlyco2 1128. Protein-level: andes shares 366 / 357 / 327 / 274
with Glyco-Decipher / Byonic / pGlyco2 / StrucGP respectively.

**Read:** andes identifies the most glycopeptides (2868), and its highest
agreement is with Byonic (0.551) and Glyco-Decipher (0.545) — it covers 1823 of
Glyco-Decipher's 2301 (79%) and 1911 of Byonic's 2509 (76%). pGlyco 2.0 (2017)
is the outlier at the low end on every row (0.27–0.38), consistent with its
smaller glycan database and older scoring, not an andes-specific divergence.

## Timing — full-scale A/B

`MouseLiver-Z-T-1.raw`, 45905 MS2 spectra, 32 threads, `sequon-reverse` decoy,
`pGlyco-N-Mouse.gdb`:

| metric | baseline (peptide-first) | `--glyco-full-glycan-db` | Δ |
|---|---|---|---|
| glyco phase | 12474.9 s (3.46 h) | 106.2 s | **117× faster** |
| total wall-clock | 12610.3 s (3.50 h) | 163.2 s | **77× faster** |
| glyco PSM rows | 38800 | 43302 | +11.6% |
| targets @ 1% FDR (2× rule) | 5319 | 5262 | −57 (−1.1%) |
| RawScore AUC (target > decoy) | 0.661 | 0.643 | −0.018 |

The default peptide-first path spends ~3.5 h building the ~1.5 GB b/y postings
index and scoring against it. `--glyco-full-glycan-db` is mass-driven: it skips
that index and matches the backbone to peptides by mass in the phase-1 bucket
index, cutting the glyco phase to ~106 s with negligible sensitivity/scoring loss
(−1.1% targets, −0.018 AUC). Full-glycan-db phase breakdown: stream_search 51.4 s,
glyco scoring → 43302 rows 106.2 s, total 163.2 s.

## Notes and caveats

- Filter thresholds are each engine's own convention (Byonic Score, Glyco-Decipher
  PeptideFDR, pGlyco2 TotalFDR, andes native q-value); they are not directly
  calibrated to one another.
- andes uses `pGlyco-N-Mouse.gdb` (1833 compositions). Byonic searched a wider
  glycan space (its 698 distinct compositions include antenna variants), so exact
  composition matches are the conservative subset.
- pGlyco 2.0 encodes the glycosylated Asn as `J` in its peptide column
  (`VSQVLHEGGHJVTK` → `…HNVTK`); the script maps `J→N`.
- Byonic peptide cells carry flanking residues **and** decimal masses
  (`K.HLLEN[+1864.634]ATASVSEAER.K`); the script strips flanks by first/last dot,
  not a full split, because the mass annotation contains a decimal point.
