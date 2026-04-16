# MS-GF+ Project — Claude Context

## Overview

MS-GF+ is a mass spectrometry database search tool for peptide identification.
The codebase is Java (Maven build). Benchmark harness scripts are local-only (not committed).

## Branch

Primary integration branch: `dev`

## Key Directories

- `src/main/java/edu/ucsd/msjava/` — core Java source
  - `msdbsearch/` — database search engine (DBScanner, ScoredSpectraMap)
  - `msutil/` — spectrum utilities (SpecKey, SpecKeyResult, SpectrumMetadata)
  - `mzid/` — mzIdentML output generation
  - `mzml/` — mzML parser (StaxMzMLParser — streaming rewrite)
  - `parser/` — input file parsers (MgfSpectrumParser, etc.)
  - `ui/` — CLI entry points (MSGFPlus, MSGFDB)
- Local benchmark harness/scripts are intentionally out-of-tree and not committed as `benchmark/`
- `src/test/` — unit tests

## Build

```bash
mvn -B verify
```

## Conventions

- Java 17+
- Maven for dependency management
- mzIdentML (`.mzid`) as primary output format
- TSV export via DirectTSVWriter
