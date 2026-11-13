# `--label` benchmarks

Scripts behind the real-data numbers of DOCS.md §10 *Labelled pairs*.

| script | what it does |
|---|---|
| `openms_ffm.py MZML LABELS MAX_LABELLED TSV [TOL_PPM] [INTENSITY_CUTOFF]` | runs OpenMS FeatureFinderMultiplex (pyOpenMS) on one run and writes its multiplets (`rt`, `charge`, `map:mz:intensity` per channel) |
| `compare_triplex.py LABEL_PARQUET FFM_TSV` | matches andes label groups to OpenMS multiplets (10 ppm, ±30 s) and compares M/L and H/L per lysine / arginine class |
| `compare_duplex.py LABEL_PARQUET RUN FFM_TSV DESIGN_LOG2` | the same for a duplex run with a known log2 H/L |
| `mzml_to_mzxml.py MZML MZXML [--ms2-from CENTROIDED_MZML]` | writes the indexed mzXML MaxQuant 2.0.3 can read (see below) |
| `make_mqpar.py TEMPLATE OUT FASTA LABELS FILE...` | fills a MaxQuant parameter file for a SILAC run |
| `compare_maxquant.py LABEL_PARQUET EVIDENCE_TXT RUN DESIGN...` | matches andes channel ratios to MaxQuant evidence ratios on (sequence, charge); DESIGN is the log2 expected ratio per channel, or `class` to split lysine and arginine peptides |
| `andes_to_idxml.py Q_TSV PIN RUN MZML OUT` | writes one run's andes PSMs (PEP scores, labels as named modifications) as idXML for OpenMS MS1LabeledWorkflow |
| `compare_ms1labeled.py QPX_FEATURE_PARQUET LABEL_PARQUET RUN=DESIGN...` | matches MS1LabeledWorkflow ratios (its QPX feature table) to andes channel ratios against a known design |
| `match_by_q.py LABEL_PARQUET` | per run, the log-ratio of channels matched between channels by quality q-value bin, against groups with both channels identified |

PXD022659's HeLa ratio mixtures are Thermo `.raw` (PanoramaPublic `RawFiles/raw`, listed with
`?method=JSON`), read directly by MaxQuant and converted with ThermoRawFileParser for andes and
OpenMS. MS1LabeledWorkflow ran on andes's identifications after PeptideIndexer (target-decoy database
from DecoyDatabase, `-unmatched_action remove`), in the 3.6.0 release and built from OpenMS
`develop` (commit `fbb25191d7`) for its targeted quantification.

Data: PXD003327 `Chris_Ecoli_1-2-4.mzML` (PRIDE); PXD022659 `Ecoli_LH_1_1_66_DDA1`,
`Ecoli_LH_1_10_66_DDA1` mzML (PanoramaPublic, `RawFiles/mzml`); PXD000714
`20120524_Renvoise_Triplex_rep1F3.mzXML` (PRIDE, converted to mzML with pyOpenMS).

MaxQuant is 2.0.3.0 from Bioconda (`maxquant`, run under Mono), with the
template from `maxquant --create`. Its mzML reader keeps only MS2 isolation windows
that recur once per cycle (a DIA layout), so on DDA mzML it drops every MS2
spectrum. Under Mono it also fails to inflate zlib arrays. The vendor files of these runs are not public
(PXD003327 deposited mzML only; the PanoramaPublic `RawFiles/raw` folder needs a
login), so MaxQuant reads mzXML written by `mzml_to_mzxml.py`. That reader takes
the analyzer from a Thermo-style `filterLine` and misreads 64-bit peaks, so the
converter writes 32-bit peaks and one filter line per scan, and it drops empty
spectra, on which MaxQuant's MS1 peak detection and MS2 top-N filter crash. For PXD003327, MS1 stays profile
and MS2 comes centroided from OpenMS PeakPickerHiRes (`--ms2-from`), the same
input andes is compared on. For the TripleTOF runs of PXD022659, whose MS1 is centroided,
MaxQuant computes NaN masses from centroided MS1 and crashes in feature detection.
`--ms1-profile-resolution 30000` rewrites each centroid as a 5-point profile peak of that
FWHM, which gets past it, but no setting gave a usable run (DOCS.md §10). The FTMS MS/MS
tolerance was raised to 40 ppm, the TOF default, since the reader labels every scan FTMS.
