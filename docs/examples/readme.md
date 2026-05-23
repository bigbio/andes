# Example configuration and sample outputs

Broader manuals (MS-GF+ CLI, BuildSA, changelog, …) are Markdown pages one level up: [`../readme.md`](../readme.md).

This directory holds **small text examples** shipped with the repository:

| File | Purpose |
|------|---------|
| `MSGFPlus_Params.txt` | Annotated MS-GF+ configuration / parameter reference (linked from the README and CLI help). |
| `enzymes.txt` | Example custom enzyme definitions (see CLI help for `params/enzymes.txt`). |
| `activationMethods.txt` | Example custom activation methods for `params/activationMethods.txt`. |
| `protocols.txt` | Example custom protocols for `params/protocols.txt`. |
| `Mods.txt` | Example modification file. |
| `pxd001819_example.pin` | **Sample Percolator `.pin` output** from a PXD001819 (yeast + UPS1 on LTQ Orbitrap Velos) search. Header + 20 target PSMs + 10 decoy PSMs, chosen for peptide-sequence diversity so every column is represented. Use this to inspect the `.pin` schema without running a full search. Full column reference in [`../output.md`](../output.md). |
| `test.tsv`, `test_Unrolled.tsv` | Example TSV exports for documentation. |

**Not stored here:** tutorial spreadsheets, plots, or bundled FASTA/index files. Those are either removed as non-essential bloat or live under `src/test/resources/` for automated tests.

To build a suffix-array index for your own database, run `BuildSA` (see `java … edu.ucsd.msjava.msdbsearch.BuildSA -h` after building the JAR).
