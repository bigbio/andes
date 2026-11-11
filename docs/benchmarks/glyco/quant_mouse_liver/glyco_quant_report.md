# Glycopeptide LFQ + MBR benchmark — mouse liver, 5 fractions (deep tier)

Runs: MouseLiver-Z-T-1, MouseLiver-Z-T-2, MouseLiver-Z-T-3, MouseLiver-Z-T-4, MouseLiver-Z-T-5.  Feature acceptance in the emitted tables: feature q ≤ 0.05 (`Accepted`).  Identification gate applied here: Percolator q ≤ 0.01 on the PSM each feature came from (for transferred rows, the source PSM).

## 0. Identification layer (pooled Percolator on the .glyco.pin)

Percolator target PSMs: 115,629; at q ≤ 0.01: 33,823.

| Percolator gate | target PSMs | entrapment | corrected FDP |
|---|---:|---:|---:|
| q ≤ 0.01 | 33,823 | 148 (0.44%) | 0.93% |
| q ≤ 0.05 | 35,747 | 690 (1.93%) | 4.11% |

## 1. What the quant tables contain

| run | identified rows | transferred rows | accepted identified | accepted transferred | accepted glycopeptides |
|---|---:|---:|---:|---:|---:|
| MouseLiver-Z-T-1 | 23,839 | 45,859 | 21,243 | 28,427 | 43,774 |
| MouseLiver-Z-T-2 | 21,801 | 47,586 | 19,387 | 29,765 | 43,452 |
| MouseLiver-Z-T-3 | 23,115 | 46,223 | 20,554 | 32,378 | 47,266 |
| MouseLiver-Z-T-4 | 23,293 | 46,314 | 20,532 | 31,215 | 46,006 |
| MouseLiver-Z-T-5 | 23,581 | 46,046 | 20,737 | 30,002 | 45,135 |

Accepted: **254,240 features**, **60,790 glycopeptides** (peptide, site, composition), **64,070** at topology resolution; 2,439 compositions carry more than one retained topology.

## 2. Identification control: feature q is not identification q

`Accepted` is a measurement/XIC q-value.  Joining each feature to the Percolator q-value of its PSM shows what it lets through:

| feature set | features | entrapment raw | sequon-corrected FDP |
|---|---:|---:|---:|
| accepted, no ID gate | 254,240 | 42.7% | 91.1% |
| accepted + ID q ≤ 0.01 | 40,416 | 0.9% | 1.8% |
| accepted + ID q > 0.05 / no q | 209,773 | 51.0% | 108.6% |

Entrapment by evidence class:

| evidence | gate | features | entrapment |
|---|---|---:|---:|
| Identified | no ID gate | 102,453 | 34,802 (33.97%) |
| Identified | q ≤ 0.01 | 32,967 | 136 (0.41%) |
| Transferred | no ID gate | 151,787 | 73,876 (48.67%) |
| Transferred | q ≤ 0.01 | 7,449 | 209 (2.81%) |

Only **40,416** of 254,240 accepted features (15.9%) pass the 1% identification gate; 84% of what the tables report as quantified glycopeptides come from PSMs that are mostly random rank-1 matches.  All sections after this one are computed **after** the identification gate.

Gated features: 40,416 (32,967 identified + 7,449 transferred).

| run | identified | transferred | unique glycopeptides (composition) | at topology resolution |
|---|---:|---:|---:|---:|
| MouseLiver-Z-T-1 | 6,967 | 1,285 | 3,990 | 4,585 |
| MouseLiver-Z-T-2 | 6,399 | 1,617 | 3,812 | 4,376 |
| MouseLiver-Z-T-3 | 6,667 | 1,504 | 4,074 | 4,617 |
| MouseLiver-Z-T-4 | 6,306 | 1,639 | 3,887 | 4,452 |
| MouseLiver-Z-T-5 | 6,628 | 1,404 | 3,961 | 4,526 |

Gated: **4,661** glycopeptides at composition resolution, **6,392** at topology resolution; 1,131 compositions carry more than one retained topology.

## 3. Transfers (MBR)

Accepted transfers (feature q ≤ 0.05): 7,449 after the ID gate / 151,787 before it.
Transfer feature q: median 0.0174, p90 0.0337, max 0.0499.

## 4. Quantification matrix, cross-run consistency, proteins

### All accepted (as emitted, no ID gate) — 254,240 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 60,790 |
| >=2 runs | 55,300 |
| >=3 runs | 47,908 |
| >=4 runs | 37,837 |
| >=5 runs | 23,798 |
| **unique glycopeptides** | **60,790** |
| matrix cells filled | 225,633 (74.2% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 33,853 | 0.663 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 36,179 | 0.667 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 34,988 | 0.668 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 33,879 | 0.634 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 36,850 | 0.676 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 35,448 | 0.679 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 35,092 | 0.657 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 39,384 | 0.721 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 37,686 | 0.706 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 36,460 | 0.699 |

Median pairwise r = 0.672.

Same comparison split by whether a measurement was identified or transferred in:

| run pair | class | n | Pearson |
|---|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both identified | 3,239 | 0.852 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | one transferred | 15,427 | 0.604 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both transferred | 15,187 | 0.616 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both identified | 3,460 | 0.861 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | one transferred | 16,398 | 0.605 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both transferred | 16,321 | 0.632 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both identified | 3,273 | 0.856 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | one transferred | 16,045 | 0.614 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both transferred | 15,670 | 0.623 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both identified | 3,306 | 0.866 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | one transferred | 15,157 | 0.575 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both transferred | 15,416 | 0.582 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both identified | 3,184 | 0.831 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | one transferred | 16,818 | 0.629 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both transferred | 16,848 | 0.642 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both identified | 3,161 | 0.842 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | one transferred | 15,966 | 0.629 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both transferred | 16,321 | 0.638 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both identified | 3,184 | 0.834 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | one transferred | 15,970 | 0.608 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both transferred | 15,938 | 0.613 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both identified | 3,266 | 0.850 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | one transferred | 18,833 | 0.683 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both transferred | 17,285 | 0.692 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both identified | 3,409 | 0.873 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | one transferred | 17,266 | 0.668 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both transferred | 17,011 | 0.671 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both identified | 3,319 | 0.866 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | one transferred | 16,652 | 0.673 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both transferred | 16,489 | 0.643 |

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0991 (359,819 pairs).

Proteins: 13,498; median cross-run Top3 r = 0.750.

Entrapment: 108,678 / 254,240 accepted features (42.75% raw); sequon-corrected FDP ≈ 91.08%.

### Identified only, ID-gated — 32,967 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 4,661 |
| >=2 runs | 3,134 |
| >=3 runs | 2,466 |
| >=4 runs | 1,962 |
| >=5 runs | 1,441 |
| **unique glycopeptides** | **4,661** |
| matrix cells filled | 13,664 (58.6% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 1,953 | 0.892 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 2,106 | 0.881 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 1,975 | 0.899 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 2,064 | 0.899 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 1,922 | 0.853 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 1,838 | 0.882 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 1,885 | 0.873 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 1,949 | 0.886 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 2,072 | 0.921 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 1,952 | 0.910 |

Median pairwise r = 0.889.

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0752 (19,716 pairs).

Proteins: 657; median cross-run Top3 r = 0.940.

Entrapment: 136 / 32,967 accepted features (0.41% raw); sequon-corrected FDP ≈ 0.88%.

### Identified + MBR, ID-gated — 40,416 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 4,661 |
| >=2 runs | 4,454 |
| >=3 runs | 4,150 |
| >=4 runs | 3,646 |
| >=5 runs | 2,813 |
| **unique glycopeptides** | **4,661** |
| matrix cells filled | 19,724 (84.6% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 3,416 | 0.824 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 3,595 | 0.819 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 3,458 | 0.834 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 3,485 | 0.830 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 3,522 | 0.797 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 3,354 | 0.821 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 3,431 | 0.810 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 3,603 | 0.841 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 3,620 | 0.861 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 3,460 | 0.846 |

Median pairwise r = 0.827.

Same comparison split by whether a measurement was identified or transferred in:

| run pair | class | n | Pearson |
|---|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both identified | 1,953 | 0.900 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | one transferred | 958 | 0.517 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both transferred | 505 | 0.636 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both identified | 2,106 | 0.891 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | one transferred | 1,027 | 0.526 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both transferred | 462 | 0.624 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both identified | 1,975 | 0.904 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | one transferred | 995 | 0.578 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both transferred | 488 | 0.600 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both identified | 2,064 | 0.906 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | one transferred | 984 | 0.544 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both transferred | 437 | 0.616 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both identified | 1,922 | 0.862 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | one transferred | 1,024 | 0.499 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both transferred | 576 | 0.631 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both identified | 1,838 | 0.890 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | one transferred | 880 | 0.550 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both transferred | 636 | 0.610 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both identified | 1,885 | 0.880 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | one transferred | 1,000 | 0.560 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both transferred | 546 | 0.551 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both identified | 1,949 | 0.891 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | one transferred | 1,091 | 0.600 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both transferred | 563 | 0.702 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both identified | 2,072 | 0.926 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | one transferred | 1,051 | 0.660 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both transferred | 497 | 0.643 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both identified | 1,952 | 0.916 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | one transferred | 952 | 0.638 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both transferred | 556 | 0.602 |

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0631 (34,944 pairs).

Proteins: 657; median cross-run Top3 r = 0.911.

Entrapment: 345 / 40,416 accepted features (0.85% raw); sequon-corrected FDP ≈ 1.82%.

## 5. Run log (quant / MBR lines)

- `Quantification: label-free MS1 (tol 10 ppm, rt window ±60 s, 4 isotopes, min cosine 0.7, feature FDR 0.05, smoothing Savitzky–Golay 5 scans)`
- `[glyco-quant] 95763 emitted rows are decoy-backbone hits and are not quantified (the PIN carries them for Percolator; a decoy precursor is not a glycopeptide)`
- `[glyco-quant] identified MouseLiver-Z-T-1: 23019 features in 23839 targets, 1.49s`
- `[glyco-quant] identified MouseLiver-Z-T-2: 21051 features in 21801 targets, 2.42s`
- `[glyco-quant] identified MouseLiver-Z-T-3: 22361 features in 23115 targets, 1.40s`
- `[glyco-quant] identified MouseLiver-Z-T-4: 22499 features in 23293 targets, 1.20s`
- `[glyco-quant] identified MouseLiver-Z-T-5: 22781 features in 23581 targets, 1.18s`
- `[glyco-quant] transferred MouseLiver-Z-T-1: 33376 features in 45859 targets, 1.51s`
- `[glyco-quant] transferred MouseLiver-Z-T-2: 34810 features in 47586 targets, 1.56s`
- `[glyco-quant] transferred MouseLiver-Z-T-3: 36511 features in 46223 targets, 1.32s`
- `[glyco-quant] transferred MouseLiver-Z-T-4: 35932 features in 46314 targets, 2.05s`
- `[glyco-quant] transferred MouseLiver-Z-T-5: 34750 features in 46046 targets, 1.97s`
- `[glyco-mbr] 5 run maps (14560 anchors); 232028 transfers attempted, 175379 measured, 151787 accepted`
- `[glyco-quant] quantified 287090 features from 347657 targets; 254240 at feature q <= 0.05 (accepted); 225633 glycopeptides, 132753 glycoform-site spectra across 5 runs`

