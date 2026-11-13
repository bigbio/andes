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
| MouseLiver-Z-T-1 | 23,839 | 46,608 | 21,519 | 24,569 | 40,085 |
| MouseLiver-Z-T-2 | 21,801 | 48,369 | 19,624 | 24,855 | 38,664 |
| MouseLiver-Z-T-3 | 23,115 | 46,960 | 20,835 | 25,832 | 40,977 |
| MouseLiver-Z-T-4 | 23,293 | 47,079 | 20,796 | 25,420 | 40,442 |
| MouseLiver-Z-T-5 | 23,581 | 46,776 | 21,030 | 25,548 | 40,870 |

Accepted: **230,028 features**, **61,724 glycopeptides** (peptide, site, composition), **65,036** at topology resolution; 2,465 compositions carry more than one retained topology.

## 2. Identification control: feature q is not identification q

`Accepted` is a measurement/XIC q-value.  Joining each feature to the Percolator q-value of its PSM shows what it lets through:

| feature set | features | entrapment raw | sequon-corrected FDP |
|---|---:|---:|---:|
| accepted, no ID gate | 230,028 | 41.8% | 89.0% |
| accepted + ID q ≤ 0.01 | 40,969 | 0.9% | 1.8% |
| accepted + ID q > 0.05 / no q | 184,826 | 51.0% | 108.7% |

Entrapment by evidence class:

| evidence | gate | features | entrapment |
|---|---|---:|---:|
| Identified | no ID gate | 103,804 | 35,362 (34.07%) |
| Identified | q ≤ 0.01 | 33,262 | 138 (0.41%) |
| Transferred | no ID gate | 126,224 | 60,741 (48.12%) |
| Transferred | q ≤ 0.01 | 7,707 | 212 (2.75%) |

Only **40,969** of 230,028 accepted features (17.8%) pass the 1% identification gate; 82% of what the tables report as quantified glycopeptides come from PSMs that are mostly random rank-1 matches.  All sections after this one are computed **after** the identification gate.

Gated features: 40,969 (33,262 identified + 7,707 transferred).

| run | identified | transferred | unique glycopeptides (composition) | at topology resolution |
|---|---:|---:|---:|---:|
| MouseLiver-Z-T-1 | 7,028 | 1,371 | 4,078 | 4,695 |
| MouseLiver-Z-T-2 | 6,446 | 1,668 | 3,844 | 4,428 |
| MouseLiver-Z-T-3 | 6,735 | 1,488 | 4,070 | 4,635 |
| MouseLiver-Z-T-4 | 6,369 | 1,680 | 3,936 | 4,523 |
| MouseLiver-Z-T-5 | 6,684 | 1,500 | 4,031 | 4,621 |

Gated: **4,692** glycopeptides at composition resolution, **6,436** at topology resolution; 1,141 compositions carry more than one retained topology.

## 3. Transfers (MBR)

Accepted transfers (feature q ≤ 0.05): 7,707 after the ID gate / 126,224 before it.
Transfer feature q: median 0.0259, p90 0.0313, max 0.0399.

## 4. Quantification matrix, cross-run consistency, proteins

### All accepted (as emitted, no ID gate) — 230,028 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 61,724 |
| >=2 runs | 44,276 |
| >=3 runs | 37,458 |
| >=4 runs | 32,003 |
| >=5 runs | 25,577 |
| **unique glycopeptides** | **61,724** |
| matrix cells filled | 201,038 (65.1% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 30,917 | 0.857 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 32,409 | 0.874 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 31,333 | 0.876 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 31,530 | 0.874 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 31,364 | 0.854 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 30,960 | 0.865 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 31,052 | 0.859 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 32,494 | 0.882 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 33,140 | 0.905 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 32,310 | 0.895 |

Median pairwise r = 0.874.

Same comparison split by whether a measurement was identified or transferred in:

| run pair | class | n | Pearson |
|---|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both identified | 3,266 | 0.855 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | one transferred | 13,561 | 0.848 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both transferred | 14,090 | 0.850 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both identified | 3,485 | 0.865 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | one transferred | 14,782 | 0.850 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both transferred | 14,142 | 0.889 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both identified | 3,300 | 0.856 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | one transferred | 13,794 | 0.873 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both transferred | 14,239 | 0.873 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both identified | 3,339 | 0.868 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | one transferred | 13,803 | 0.860 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both transferred | 14,388 | 0.878 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both identified | 3,204 | 0.829 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | one transferred | 13,449 | 0.829 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both transferred | 14,711 | 0.871 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both identified | 3,190 | 0.848 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | one transferred | 13,034 | 0.865 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both transferred | 14,736 | 0.857 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both identified | 3,212 | 0.835 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | one transferred | 13,156 | 0.852 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both transferred | 14,684 | 0.860 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both identified | 3,289 | 0.848 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | one transferred | 14,558 | 0.873 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both transferred | 14,647 | 0.891 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both identified | 3,432 | 0.865 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | one transferred | 15,141 | 0.896 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both transferred | 14,567 | 0.919 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both identified | 3,345 | 0.870 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | one transferred | 14,266 | 0.892 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both transferred | 14,699 | 0.895 |

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0634 (317,509 pairs).

Proteins: 13,552; median cross-run Top3 r = 0.892.

Entrapment: 96,103 / 230,028 accepted features (41.78% raw); sequon-corrected FDP ≈ 89.02%.

### Identified only, ID-gated — 33,262 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 4,692 |
| >=2 runs | 3,160 |
| >=3 runs | 2,483 |
| >=4 runs | 1,975 |
| >=5 runs | 1,448 |
| **unique glycopeptides** | **4,692** |
| matrix cells filled | 13,758 (58.6% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 1,968 | 0.891 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 2,121 | 0.880 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 1,991 | 0.896 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 2,082 | 0.896 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 1,928 | 0.856 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 1,850 | 0.883 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 1,895 | 0.869 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 1,960 | 0.883 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 2,083 | 0.918 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 1,965 | 0.909 |

Median pairwise r = 0.887.

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0754 (19,843 pairs).

Proteins: 659; median cross-run Top3 r = 0.936.

Entrapment: 138 / 33,262 accepted features (0.41% raw); sequon-corrected FDP ≈ 0.88%.

### Identified + MBR, ID-gated — 40,969 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 4,692 |
| >=2 runs | 4,288 |
| >=3 runs | 4,012 |
| >=4 runs | 3,718 |
| >=5 runs | 3,249 |
| **unique glycopeptides** | **4,692** |
| matrix cells filled | 19,959 (85.1% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 3,597 | 0.897 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 3,740 | 0.894 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 3,643 | 0.911 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 3,702 | 0.905 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 3,583 | 0.869 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 3,534 | 0.900 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 3,593 | 0.882 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 3,672 | 0.903 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 3,744 | 0.925 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 3,654 | 0.919 |

Median pairwise r = 0.902.

Same comparison split by whether a measurement was identified or transferred in:

| run pair | class | n | Pearson |
|---|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both identified | 1,968 | 0.903 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | one transferred | 1,054 | 0.774 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both transferred | 575 | 0.848 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both identified | 2,121 | 0.892 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | one transferred | 1,130 | 0.796 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both transferred | 489 | 0.886 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both identified | 1,991 | 0.906 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | one transferred | 1,069 | 0.858 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both transferred | 583 | 0.886 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both identified | 2,082 | 0.909 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | one transferred | 1,108 | 0.811 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both transferred | 512 | 0.887 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both identified | 1,928 | 0.869 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | one transferred | 997 | 0.737 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both transferred | 658 | 0.838 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both identified | 1,850 | 0.895 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | one transferred | 935 | 0.822 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both transferred | 749 | 0.863 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both identified | 1,895 | 0.880 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | one transferred | 1,045 | 0.773 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both transferred | 653 | 0.840 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both identified | 1,960 | 0.894 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | one transferred | 1,090 | 0.836 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both transferred | 622 | 0.900 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both identified | 2,083 | 0.926 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | one transferred | 1,107 | 0.858 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both transferred | 554 | 0.891 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both identified | 1,965 | 0.916 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | one transferred | 1,054 | 0.861 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both transferred | 635 | 0.894 |

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0516 (36,462 pairs).

Proteins: 659; median cross-run Top3 r = 0.964.

Entrapment: 350 / 40,969 accepted features (0.85% raw); sequon-corrected FDP ≈ 1.82%.

## 5. Run log (quant / MBR lines)

- `Quantification: intact glycopeptides from MS1 (tol 10 ppm, rt window ±120 s, 4-14 isotopes, min cosine 0.7, feature FDR 0.05, monoisotopic intensity)`
- `[glyco-quant] 95763 emitted rows are decoy-backbone hits and are not quantified (the PIN carries them for Percolator; a decoy precursor is not a glycopeptide)`
- `[glyco-quant] identified MouseLiver-Z-T-1: 23020 features in 23839 targets, 1.63s`
- `[glyco-quant] identified MouseLiver-Z-T-2: 21051 features in 21801 targets, 2.07s`
- `[glyco-quant] identified MouseLiver-Z-T-3: 22361 features in 23115 targets, 1.24s`
- `[glyco-quant] identified MouseLiver-Z-T-4: 22500 features in 23293 targets, 1.21s`
- `[glyco-quant] identified MouseLiver-Z-T-5: 22781 features in 23581 targets, 1.20s`
- `[glyco-quant] transferred MouseLiver-Z-T-1: 43969 features in 46608 targets, 3.36s`
- `[glyco-quant] transferred MouseLiver-Z-T-2: 46110 features in 48369 targets, 3.37s`
- `[glyco-quant] transferred MouseLiver-Z-T-3: 44978 features in 46960 targets, 3.01s`
- `[glyco-quant] transferred MouseLiver-Z-T-4: 45064 features in 47079 targets, 3.03s`
- `[glyco-quant] transferred MouseLiver-Z-T-5: 44783 features in 46776 targets, 3.04s`
- `[glyco-mbr] 5 run maps (14675 anchors); 235792 transfers attempted, 224904 measured, 126224 accepted`
- `[glyco-quant] quantified 336617 features from 351421 targets; 230028 accepted (103804 identified, quality-gated only; 126224 transferred at feature q <= 0.05); 201038 glycopeptides, 120257 glycoform-site spectra across 5 runs`

