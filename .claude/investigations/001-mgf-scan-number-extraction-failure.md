# Investigation 001: MGF Scan Number Extraction Failure

**Status:** OPEN
**Date observed:** 2026-04-15
**Severity:** Medium — functional (spectra still searched, but scan numbers missing in output)
**Branch:** `feature/streaming-mzml-parser`

## What Was Observed

When running the baseline benchmark against MGF files, MS-GF+ emits repeated warnings:

```
Unable to extract the scan number from the title: id=PXD002047;TCGA-AA-A02O-01A-23_W_VU_20130205_A0218_10A_R_FR05.mzML;controllerType=0
Expected format is DatasetName.ScanStart.ScanEnd.Charge
```

The warning appeared for every spectrum in the MGF file (`test.mgf`), suggesting
the entire file uses a TITLE format that the parser cannot handle.

## Where It Was Observed

- **Run:** Baseline benchmark (`baseline/MSGFPlus.jar`, v2026.03.25)
- **Input:** `test.mgf` — MGF file with TITLE lines in PRIDE/ProteomeXchange format
- **Database:** `human-uniprot-contaminants.revCat.fasta`

## Relevant Code

### `MgfSpectrumParser.extractScanRangeFromTitle()` — the parser

```
src/main/java/edu/ucsd/msjava/parser/MgfSpectrumParser.java:278-316
```

The method splits the title on `.` and expects:
- `token.length > 3` → `DatasetName.ScanStart.ScanEnd.Charge`
- `token.length == 3 && title.endsWith(".")` → `DatasetName.ScanStart.ScanEnd.`

The PRIDE-format title `id=PXD002047;TCGA-AA-A02O-01A-23_W_VU_20130205_A0218_10A_R_FR05.mzML;controllerType=0`
splits to `["id=PXD002047;TCGA-AA-A02O-01A-23_W_VU_20130205_A0218_10A_R_FR05", "mzML;controllerType=0"]`
(only 2 tokens), so it falls through to the `else` branch and emits the warning.

### `MgfSpectrumParser.warnScanNotFoundInTitle()` — the warning

```
src/main/java/edu/ucsd/msjava/parser/MgfSpectrumParser.java:384-392
```

Capped at `MAX_SCAN_MISSING_WARNINGS` prints, then silently counts the rest.
Final total printed by `SpecKey.java:139`.

## Hypotheses

1. **Title format mismatch (most likely):** The MGF file uses a PRIDE/ProteomeXchange
   `TITLE` format that encodes the source file reference and controller info with
   semicolons, not the `Dataset.Start.End.Charge` convention. The parser has no
   fallback for alternative formats.

2. **Possible alternative scan encodings in TITLE:** Some MGF generators embed scan
   numbers as `scan=NNNN` or `scans=NNNN` within the TITLE string. The parser
   doesn't attempt to extract these.

3. **`index=` fallback:** When scan extraction fails, the spectrum gets assigned
   `index=N` as its ID (from `specIndexMap`). This means the mzIdentML output
   will reference spectra by index rather than native scan number, which may
   affect downstream tools that expect scan-based references.

## Impact

- **Search results:** Not affected — MS-GF+ still searches the spectra correctly.
- **Output traceability:** Degraded — mzIdentML references use index instead of
  native scan IDs, making it harder to trace PSMs back to the raw data.
- **Benchmark:** May cause metric discrepancies if downstream scripts parse scan
  numbers from the mzIdentML output.

## Potential Fixes

1. Add regex-based fallback in `extractScanRangeFromTitle()` to detect patterns like:
   - `scan=(\d+)` or `scans=(\d+)`
   - `spectrum=(\d+)`
   - `index=(\d+)`
2. Support PRIDE USI-style TITLE parsing: extract scan from
   `controllerType=0 controllerNumber=1 scan=NNNN` if present.
3. Allow users to specify a scan number extraction regex via CLI parameter.

## Next Steps

- [ ] Examine the actual MGF file to see the full TITLE line format
- [ ] Check if `scan=` or similar key-value pairs are embedded in the TITLE
- [ ] Review how other tools (MaxQuant, Comet, X!Tandem) handle non-standard TITLE formats
- [ ] Decide on backward-compatible fix approach
- [ ] Add unit test covering PRIDE-format TITLE strings
