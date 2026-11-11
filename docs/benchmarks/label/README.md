# `--label` benchmarks

Scripts behind the real-data numbers of DOCS.md §10 *Labelled pairs*.

| script | what it does |
|---|---|
| `openms_ffm.py MZML LABELS MAX_LABELLED TSV [TOL_PPM] [INTENSITY_CUTOFF]` | runs OpenMS FeatureFinderMultiplex (pyOpenMS) on one run and writes its multiplets (`rt`, `charge`, `map:mz:intensity` per channel) |
| `compare_triplex.py LABEL_PARQUET FFM_TSV` | matches andes label groups to OpenMS multiplets (10 ppm, ±30 s) and compares M/L and H/L per lysine / arginine class |
| `compare_duplex.py LABEL_PARQUET RUN FFM_TSV DESIGN_LOG2` | the same for a duplex run with a known log2 H/L |
| `match_by_q.py LABEL_PARQUET` | per run, the log-ratio of channels matched between channels by quality q-value bin, against groups with both channels identified |

Data: PXD003327 `Chris_Ecoli_1-2-4.mzML` (PRIDE); PXD022659 `Ecoli_LH_1_1_66_DDA1`,
`Ecoli_LH_1_10_66_DDA1` mzML (PanoramaPublic, `RawFiles/mzml`); PXD000714
`20120524_Renvoise_Triplex_rep1F3.mzXML` (PRIDE, converted to mzML with pyOpenMS).
