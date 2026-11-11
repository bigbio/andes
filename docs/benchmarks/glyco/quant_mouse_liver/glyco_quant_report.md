# Glycopeptide LFQ + MBR benchmark — mouse liver, 5 fractions (deep tier)

Runs: MouseLiver-Z-T-1, MouseLiver-Z-T-2, MouseLiver-Z-T-3, MouseLiver-Z-T-4, MouseLiver-Z-T-5.  Feature acceptance in the emitted tables: feature q ≤ 0.05 (`Accepted`).  Identification gate applied here: Percolator q ≤ 0.01 on the PSM each feature came from (for transferred rows, the source PSM).

## 0. Identification layer (pooled Percolator on the .glyco.pin)

Percolator target PSMs: 113,534; at q ≤ 0.01: 33,625.

| Percolator gate | target PSMs | entrapment | corrected FDP |
|---|---:|---:|---:|
| q ≤ 0.01 | 33,625 | 145 (0.43%) | 0.92% |
| q ≤ 0.05 | 35,455 | 681 (1.92%) | 4.09% |

## 1. What the quant tables contain

| run | identified rows | transferred rows | accepted identified | accepted transferred | accepted glycopeptides |
|---|---:|---:|---:|---:|---:|
| MouseLiver-Z-T-1 | 23,472 | 44,756 | 20,869 | 27,874 | 42,880 |
| MouseLiver-Z-T-2 | 21,395 | 46,515 | 19,046 | 29,923 | 43,167 |
| MouseLiver-Z-T-3 | 22,673 | 45,157 | 20,094 | 31,928 | 46,402 |
| MouseLiver-Z-T-4 | 22,775 | 45,361 | 19,989 | 30,536 | 44,852 |
| MouseLiver-Z-T-5 | 23,232 | 44,946 | 20,364 | 30,181 | 44,941 |

Accepted: **250,804 features**, **59,337 glycopeptides** (peptide, site, composition), **62,618** at topology resolution; 2,436 compositions carry more than one retained topology.

## 2. Identification control: feature q is not identification q

`Accepted` is a measurement/XIC q-value.  Joining each feature to the Percolator q-value of its PSM shows what it lets through:

| feature set | features | entrapment raw | sequon-corrected FDP |
|---|---:|---:|---:|
| accepted, no ID gate | 250,804 | 42.7% | 90.9% |
| accepted + ID q ≤ 0.01 | 40,430 | 1.0% | 2.1% |
| accepted + ID q > 0.05 / no q | 206,478 | 51.0% | 108.6% |

Entrapment by evidence class:

| evidence | gate | features | entrapment |
|---|---|---:|---:|
| Identified | no ID gate | 100,362 | 33,850 (33.73%) |
| Identified | q ≤ 0.01 | 32,765 | 135 (0.41%) |
| Transferred | no ID gate | 150,442 | 73,132 (48.61%) |
| Transferred | q ≤ 0.01 | 7,665 | 259 (3.38%) |

Only **40,430** of 250,804 accepted features (16.1%) pass the 1% identification gate; 84% of what the tables report as quantified glycopeptides come from PSMs that are mostly random rank-1 matches.  All sections after this one are computed **after** the identification gate.

Gated features: 40,430 (32,765 identified + 7,665 transferred).

| run | identified | transferred | unique glycopeptides (composition) | at topology resolution |
|---|---:|---:|---:|---:|
| MouseLiver-Z-T-1 | 6,898 | 1,293 | 3,959 | 4,544 |
| MouseLiver-Z-T-2 | 6,384 | 1,746 | 3,867 | 4,449 |
| MouseLiver-Z-T-3 | 6,625 | 1,497 | 4,058 | 4,587 |
| MouseLiver-Z-T-4 | 6,253 | 1,686 | 3,908 | 4,474 |
| MouseLiver-Z-T-5 | 6,605 | 1,443 | 3,956 | 4,529 |

Gated: **4,622** glycopeptides at composition resolution, **6,359** at topology resolution; 1,124 compositions carry more than one retained topology.

## 3. Transfers (MBR)

Accepted transfers (feature q ≤ 0.05): 7,665 after the ID gate / 150,442 before it.
Transfer feature q: median 0.0167, p90 0.0321, max 0.0481.

## 4. Quantification matrix, cross-run consistency, proteins

### All accepted (as emitted, no ID gate) — 250,804 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 59,337 |
| >=2 runs | 54,152 |
| >=3 runs | 46,944 |
| >=4 runs | 37,635 |
| >=5 runs | 24,174 |
| **unique glycopeptides** | **59,337** |
| matrix cells filled | 222,242 (74.9% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 34,141 | 0.677 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 35,508 | 0.659 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 34,268 | 0.668 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 33,359 | 0.625 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 36,730 | 0.679 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 35,320 | 0.685 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 34,878 | 0.658 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 37,928 | 0.697 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 38,394 | 0.722 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 37,115 | 0.715 |

Median pairwise r = 0.678.

Same comparison split by whether a measurement was identified or transferred in:

| run pair | class | n | Pearson |
|---|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both identified | 3,192 | 0.852 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | one transferred | 15,727 | 0.611 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both transferred | 15,222 | 0.640 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both identified | 3,330 | 0.854 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | one transferred | 16,122 | 0.606 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both transferred | 16,056 | 0.616 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both identified | 3,223 | 0.860 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | one transferred | 15,576 | 0.614 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both transferred | 15,469 | 0.620 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both identified | 3,267 | 0.851 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | one transferred | 14,565 | 0.572 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both transferred | 15,527 | 0.569 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both identified | 3,120 | 0.825 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | one transferred | 16,658 | 0.626 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both transferred | 16,952 | 0.651 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both identified | 3,052 | 0.853 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | one transferred | 15,775 | 0.635 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both transferred | 16,493 | 0.644 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both identified | 3,066 | 0.839 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | one transferred | 15,454 | 0.607 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both transferred | 16,358 | 0.618 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both identified | 3,174 | 0.850 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | one transferred | 17,554 | 0.656 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both transferred | 17,200 | 0.662 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both identified | 3,320 | 0.874 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | one transferred | 17,984 | 0.694 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both transferred | 17,090 | 0.685 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both identified | 3,237 | 0.873 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | one transferred | 17,402 | 0.678 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both transferred | 16,476 | 0.677 |

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0980 (357,641 pairs).

Proteins: 13,382; median cross-run Top3 r = 0.757.

Entrapment: 106,982 / 250,804 accepted features (42.66% raw); sequon-corrected FDP ≈ 90.89%.

### Identified only, ID-gated — 32,765 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 4,622 |
| >=2 runs | 3,097 |
| >=3 runs | 2,457 |
| >=4 runs | 1,945 |
| >=5 runs | 1,430 |
| **unique glycopeptides** | **4,622** |
| matrix cells filled | 13,551 (58.6% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 1,947 | 0.892 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 2,086 | 0.883 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 1,963 | 0.900 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 2,042 | 0.896 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 1,908 | 0.850 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 1,830 | 0.881 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 1,881 | 0.874 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 1,921 | 0.883 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 2,052 | 0.923 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 1,936 | 0.912 |

Median pairwise r = 0.888.

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0758 (19,566 pairs).

Proteins: 647; median cross-run Top3 r = 0.942.

Entrapment: 135 / 32,765 accepted features (0.41% raw); sequon-corrected FDP ≈ 0.88%.

### Identified + MBR, ID-gated — 40,430 features

| presence | glycopeptides |
|---|---:|
| >=1 runs | 4,622 |
| >=2 runs | 4,434 |
| >=3 runs | 4,121 |
| >=4 runs | 3,695 |
| >=5 runs | 2,876 |
| **unique glycopeptides** | **4,622** |
| matrix cells filled | 19,748 (85.5% of grid) |

| run pair | n | Pearson log2 area |
|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | 3,488 | 0.841 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | 3,571 | 0.809 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | 3,464 | 0.837 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | 3,465 | 0.816 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | 3,565 | 0.795 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | 3,432 | 0.829 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | 3,462 | 0.807 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | 3,604 | 0.831 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | 3,667 | 0.859 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | 3,547 | 0.852 |

Median pairwise r = 0.830.

Same comparison split by whether a measurement was identified or transferred in:

| run pair | class | n | Pearson |
|---|---|---:|---:|
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both identified | 1,947 | 0.901 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | one transferred | 986 | 0.551 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-2 | both transferred | 555 | 0.667 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both identified | 2,086 | 0.890 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | one transferred | 1,032 | 0.506 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-3 | both transferred | 453 | 0.563 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both identified | 1,963 | 0.907 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | one transferred | 1,001 | 0.552 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-4 | both transferred | 500 | 0.645 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both identified | 2,042 | 0.904 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | one transferred | 978 | 0.496 |
| MouseLiver-Z-T-1 vs MouseLiver-Z-T-5 | both transferred | 445 | 0.571 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both identified | 1,908 | 0.858 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | one transferred | 1,063 | 0.533 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-3 | both transferred | 594 | 0.532 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both identified | 1,830 | 0.889 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | one transferred | 893 | 0.566 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-4 | both transferred | 709 | 0.594 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both identified | 1,881 | 0.880 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | one transferred | 993 | 0.522 |
| MouseLiver-Z-T-2 vs MouseLiver-Z-T-5 | both transferred | 588 | 0.557 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both identified | 1,921 | 0.889 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | one transferred | 1,096 | 0.587 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-4 | both transferred | 587 | 0.632 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both identified | 2,052 | 0.928 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | one transferred | 1,119 | 0.654 |
| MouseLiver-Z-T-3 vs MouseLiver-Z-T-5 | both transferred | 496 | 0.612 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both identified | 1,936 | 0.919 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | one transferred | 1,059 | 0.622 |
| MouseLiver-Z-T-4 vs MouseLiver-Z-T-5 | both transferred | 552 | 0.631 |

Glycoform share agreement at multi-run sites: mean |Δshare| 0.0630 (35,265 pairs).

Proteins: 647; median cross-run Top3 r = 0.918.

Entrapment: 394 / 40,430 accepted features (0.97% raw); sequon-corrected FDP ≈ 2.08%.

## 5. Run log (quant / MBR lines)

- `Quantification: label-free MS1 (tol 10 ppm, rt window ±60 s, 4 isotopes, min cosine 0.7, feature FDR 0.05, smoothing Savitzky–Golay 5 scans)`
- `[glyco-quant] 98447 emitted rows are decoy-backbone hits and are not quantified (the PIN carries them for Percolator; a decoy precursor is not a glycopeptide)`
- `[glyco-quant] identified MouseLiver-Z-T-1: 22688 features in 23472 targets, 1.42s`
- `[glyco-quant] identified MouseLiver-Z-T-2: 20679 features in 21395 targets, 1.43s`
- `[glyco-quant] identified MouseLiver-Z-T-3: 21922 features in 22673 targets, 1.69s`
- `[glyco-quant] identified MouseLiver-Z-T-4: 22007 features in 22775 targets, 1.66s`
- `[glyco-quant] identified MouseLiver-Z-T-5: 22416 features in 23232 targets, 1.18s`
- `[glyco-quant] transferred MouseLiver-Z-T-1: 32757 features in 44756 targets, 1.50s`
- `[glyco-quant] transferred MouseLiver-Z-T-2: 34953 features in 46515 targets, 1.53s`
- `[glyco-quant] transferred MouseLiver-Z-T-3: 35934 features in 45157 targets, 1.30s`
- `[glyco-quant] transferred MouseLiver-Z-T-4: 35261 features in 45361 targets, 1.28s`
- `[glyco-quant] transferred MouseLiver-Z-T-5: 34637 features in 44946 targets, 1.24s`
- `[glyco-mbr] 5 run maps (14264 anchors); 226735 transfers attempted, 173542 measured, 150442 accepted`
- `[glyco-quant] quantified 283254 features from 340282 targets; 250804 at feature q <= 0.05 (accepted); 222242 glycopeptides, 130946 glycoform-site spectra across 5 runs`

