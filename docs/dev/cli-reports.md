# Command-line reports

The end-of-run summary is the part of the output most users read. It is logged with
`log::info!(target: "summary", ...)`, so each line carries the standard log prefix (see
`clean-logs`), and the same content without the prefix goes to `<run>.statistics.log`
(`crates/output/src/statistics.rs`). The layout below is the message part of those lines.

## Layout

- A header line with andes version, git commit, and the command's key inputs.
- Short sections with a one-word title, each a left-aligned label column and a
  right-aligned number column:

```
andes 0.4.0 (2051616b)  search  run1.mzML  human.fasta

Input      spectra             48,213
           MS2 searched        47,901
           skipped                312   (no peaks 290, charge 0 22)
Model      hcd_highres_tryp    (instrument OrbitrapAstral, detected)
Results    PSMs  q<=0.01       52,169
           peptides            31,004
           proteins             5,812
Resources  wall                 114 s    CPU 569 s    peak RSS 4.8 GB
Output     run1.pin
```

- Width at most 100 columns. No box-drawing borders; plain spaces align columns and stay
  readable in log files and when piped.
- Every number has a label and a unit. Report rates as a percentage with its numerator.
- Always include wall time, CPU time and peak memory (see `performance`).
- Say what was decided automatically and why (model choice, calibration result, memory
  budget source), because that is what a user needs to reproduce or debug the run.

## Content rules

- Report counts at the FDR threshold the user chose, and name the threshold.
- Do not print raw-score statistics a user cannot act on; keep those in statistics.log.
- A mode (open search, DIA, glyco, xl, quant) adds its own section; it does not change the
  shape of the shared sections.
- Machine-readable output belongs in files, never parsed from the log. If a script needs a
  number, write it to a file. A new output or intermediate file is versioned Parquet, not
  TSV; PIN and Percolator files keep their formats.

## --help

Each flag: one sentence on what it does, its default, its unit. Group flags by stage. A new
flag also goes into DOCS.md in the same commit.
