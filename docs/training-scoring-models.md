# Training MS-GF+ scoring models

MS-GF+ ships with a set of pre-trained scoring models (`.param` files in
`src/main/resources/`) that cover the common combinations of activation
method, instrument type, enzyme, and protocol. The bundled set includes
HCD/QExactive/Tryp, HCD/HighRes/Tryp/TMT, and so on. If your data does
not match any of the bundled combinations, or if you want a model
specifically tuned for your instrument, you can train your own.

This page describes the recovered training entry point on this fork and
the end-to-end workflow.

## When to train a new model

You need a custom model only when:

- Your activation/instrument/enzyme/protocol combination has no bundled
  `.param` file. Run a search with `-inst HighRes -m HCD -e Tryp` and
  watch the startup log: if MS-GF+ falls back to a generic model, that's
  the signal.
- Your instrument's fragmentation pattern differs materially from the
  bundled training data (e.g. a new generation Astral run vs. a Q
  Exactive of 2014).
- You want to compare an in-house trained model against the bundled one
  to quantify the gain.

For most users on standard tryptic HCD runs, the bundled model is fine
and Phase B's calibrated precursor-window tightening (the
`-precursorCal auto` flag) is the bigger lever.

## What you need

1. **Spectra**: one or more mzML or MGF files of MS/MS data from the
   instrument and acquisition mode you want to train for.
2. **A protein database**: a FASTA covering the species in your
   training data.
3. **A modifications file**: standard MS-GF+ `Mods.txt` format.
4. **A target/decoy MS-GF+ search of (1) against (2) at standard 1%
   FDR**: this provides the annotated PSMs the trainer learns from.
   The trainer needs **a few hundred** confidently identified PSMs to
   produce a stable model; thousands is better.

The mzID input that upstream MS-GF+ supported was removed from this
fork. The trainer now accepts only TSV PSM lists. The standard MS-GF+
TSV writer (`DirectTSVWriter`, the one you get with the default search
output) produces TSVs that already match the trainer's expected format.

## Workflow

### Step 1 — Search your training data

Run an MS-GF+ search the same way you would for any project, with
`-tda 1` (target-decoy on) so the output has a `QValue` column the
trainer can filter on. Use a precision tolerance appropriate for the
instrument; for high-resolution data, `-precursorCal auto` is
recommended so Phase B's calibration tightens the window before the
search.

Example:

```sh
java -Xmx16G -jar MSGFPlus.jar \
  -s training.mzML \
  -d training.fasta \
  -mod Mods.txt \
  -t 10ppm \
  -tda 1 \
  -inst HighRes \
  -m HCD \
  -e Tryp \
  -precursorCal auto \
  -addFeatures 1 \
  -o training.tsv
```

The output `training.tsv` is the input to the trainer.

### Step 2 — Train the scoring model

Invoke `ScoringParamGen` with the same activation method, instrument
type, enzyme, and protocol you used in step 1. The trainer pulls the
high-confidence PSMs (default `QValue ≤ 0.01`), looks up each PSM's
spectrum in the directory you supply, and writes a `.param` file in
the current working directory.

```sh
java -Xmx4G -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen \
  -i training.tsv \
  -d /path/to/spectra-directory \
  -m HCD \
  -inst HighRes \
  -e Tryp \
  -protocol Standard
```

The output filename is derived from the data type: e.g. for the
arguments above, `HCD_HighRes_Tryp_Standard.param`.

### Step 3 — Use the model

Drop the `.param` file into `src/main/resources/` (rebuilding the JAR)
or supply it via the appropriate parameter file. MS-GF+ will pick it up
when the search's `(activationMethod, instrumentType, enzyme, protocol)`
tuple matches the filename.

## TSV input format

`AnnotatedSpectra` (the trainer's TSV reader) requires these columns,
identified by header name (case-insensitive):

| Column | Required | Description |
|---|---|---|
| `#SpecFile` | yes | Spectrum file name (matched by basename against `-d`). |
| `SpecID` | yes | Index or scan identifier; passed to `SpectraAccessor.getSpectrumById`. |
| `Peptide` | yes | Peptide sequence. May include `K.PEPTIDE.K` flanking residues; flankers are stripped. |
| `Charge` | yes | Integer charge state. |
| `FDR` / `EFDR` / `QValue` / `SpecQValue` | yes (any one) | Used to filter rows; default threshold is 0.01. |

Extra columns are ignored. The `Peptide` field accepts standard MS-GF+
modification syntax (e.g. `K.PEP+57.021M+15.995IDE.K`). The trainer
matches each PSM's `(SpecID, Charge)` against the spectrum file and
verifies that the peptide's theoretical mass is within 5 Da of the
spectrum's precursor mass; mismatches are reported and abort the run
(unless `-dropErrors 1` is supplied).

## CLI reference

```
java -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen [options]

Required:
  -i  <tsv1[,tsv2,...]>  Training result TSV files (mzID input not supported in this build)
  -d  <specDir>          Directory holding the spectrum files referenced by the TSVs
  -m  <activation>       Activation method (CID, ETD, HCD, UVPD, etc.)
  -inst <instrument>     Instrument type (LowRes, HighRes, QExactive, etc.)
  -e  <enzyme>           Enzyme name (Tryp, Chymotryp, LysC, AspN, etc.)

Optional:
  -protocol <name>       Protocol (default: NoProtocol/automatic)
  -thread <int>          Worker threads for parsing PSMs (default: 1)
  -dropErrors 0|1        Drop datasets with errors instead of failing (default: 0)
  -mgf 0|1               Also emit aggregated <dataType>.mgf (default: 0)
```

## Notes and limitations

- **mzID input was removed**: the upstream MS-GF+ `ui/ScoringParamGen`
  CLI accepted `.mzid` files via `MzIDParser`. That class was deleted
  from this fork in commit `9bf01c8`. If your existing pipeline produces
  mzID, pre-convert to TSV (e.g. with the upstream MS-GF+ JAR's
  `MzIDToTsv`) before passing to the trainer.
- **No `params/ParamManager`**: the upstream CLI used the (now-removed)
  `ParamManager` framework. The recovered entry point parses arguments
  manually; option semantics match the upstream CLI, but the help text
  is shorter.
- **Output goes to the current directory**: the `.param` file lands
  wherever you launched the JVM. There is no `-o` option; supply the
  intended output directory via `cd` before invoking, or copy the file
  afterwards.
- **Minimum training data**: the trainer applies an internal dedup
  (`(peptide, charge)` keyed, capped at 3 spectra per pair) before the
  partition step. Empirically (see `TestScoringParamGenSmoke`), the
  partition step refuses to fit under ~360 dedup-survived spectra at a
  single charge, silently emitting an empty partition set and aborting
  the model write. To stay above the floor, plan on **≥ 200 unique
  peptide identifications across the dominant charge state**, and
  preferably ≥ 500 unique peptides for a model with usable rank
  distributions.

## See also

- The recovered code lives at:
  - `src/main/java/edu/ucsd/msjava/ui/ScoringParamGen.java`
  - `src/main/java/edu/ucsd/msjava/msscorer/ScoringParameterGeneratorWithErrors.java`
  - `src/main/java/edu/ucsd/msjava/msscorer/ScoringParameterGenerator.java`
  - `src/main/java/edu/ucsd/msjava/msutil/AnnotatedSpectra.java`
  - `src/main/java/edu/ucsd/msjava/misc/TrainScoringParameters.java`
- The original upstream documentation page is still live at
  <https://msgfplus.github.io/msgfplus/ScoringParamGen.html>; the
  command-line semantics here match it apart from the mzID input
  caveat above.
