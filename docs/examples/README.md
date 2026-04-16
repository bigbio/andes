# Example configuration and sample outputs

Broader manuals (MS-GF+ CLI, BuildSA, MzID→TSV, changelog, …) are Markdown pages one level up: [`../README.md`](../README.md).

This directory holds **small text examples** shipped with the repository:

| File | Purpose |
|------|---------|
| `MSGFPlus_Params.txt` | Annotated MS-GF+ configuration / parameter reference (linked from the README and CLI help). |
| `enzymes.txt` | Example custom enzyme definitions (see CLI help for `params/enzymes.txt`). |
| `activationMethods.txt` | Example custom activation methods for `params/activationMethods.txt`. |
| `protocols.txt` | Example custom protocols for `params/protocols.txt`. |
| `Mods.txt` | Example modification file. |
| `test.mzid` | Tiny mzIdentML sample used in static documentation. |
| `test.tsv`, `test_Unrolled.tsv` | Example TSV exports for documentation. |

**Not stored here:** tutorial spreadsheets, plots, or bundled FASTA/index files. Those are either removed as non-essential bloat or live under `src/test/resources/` for automated tests.

To build a suffix-array index for your own database, run `BuildSA` (see `java … edu.ucsd.msjava.msdbsearch.BuildSA -h` after building the JAR).
