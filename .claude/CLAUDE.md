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
  - `mzid/` — `DirectPinWriter` + `DirectTSVWriter` (only writers retained; all mzIdentML classes + consumers deleted)
  - `mzml/` — mzML parser (StaxMzMLParser — streaming rewrite)
  - `parser/` — input file parsers (MgfSpectrumParser, etc.)
  - `ui/` — CLI entry points (MSGFPlus, MSGFDB)
- Local benchmark harness/scripts are intentionally out-of-tree and not committed as `benchmark/`
- `src/test/` — unit tests

## Build

```bash
mvn -B verify
```

**Do NOT run full `mvn test` without scoping.** The suite includes `TestPrecursorCalIntegration` which runs 4 full MS-GF+ searches on the 82 MB `human-uniprot-contaminants.fasta` fixture and takes ≥ 90 min on an idle machine. For iteration, scope to relevant classes:

```bash
mvn -B -o test -Dtest='TestDirectPinWriter,TestMassCalibrator,TestPrecursorCalScaffolding'
```

## Conventions

- Java 17+
- Maven for dependency management
- Percolator `.pin` as the default output format (mzIdentML output removed; feed downstream via Percolator)
- TSV export via DirectTSVWriter
- Percolator `.pin` export via DirectPinWriter (PR #20 + PR #22)

## Performance-sensitive invariants (learned empirically)

- **Never wrap hot-path collections in `Map.copyOf` / `ImmutableCollections`.** Observed 2.2× Astral regression — likely a bad interaction between `Partition.hashCode` clustering and ImmutableCollections' open-addressing.
- **Any optional scoring-path feature behind a flag must be bit-identical to baseline when disabled.** Implement via `if (mode == OFF) return input_unchanged;` at the top of the entry point — do NOT rely on "multiply by zero" or "flag-dependent branch deep in the loop"; both reorder float ops.
- **Pre-passes (calibrators, samplers) must not mutate shared state.** MS-GF+'s `Spectrum` objects are shared across the pre-pass and main pass; mutating them in the pre-pass (e.g. via `scorer.getScoredSpectrum(spec)`) causes silent PSM-count drift when the main pass re-reads the mutated state.

## Benchmark harness

Local-only, gitignored (`benchmark/*` with `!benchmark/README.md` / `!benchmark/ci/` carve-outs). Three 3-arm scripts per dataset:

- `benchmark/run_pxd001819_3arm.sh` / `run_astral_3arm.sh` / `run_tmt_3arm.sh` — each runs baseline JAR / branch off / branch auto and produces `.pin` files
- `benchmark/compare_*_3arm_percolator.sh` — runs Percolator via Docker (biocontainers 3.7.1) on each pin; prints 1% / 5% FDR target counts
- See `~/.claude/projects/-Users-yperez-work-msgfplus/memory/reference_benchmark_infra.md` for full details (conda env, Docker image, dataset locations)

## Next planned work

**Chimeric two-pass cascade — DRAFT PR #42 (`feat/chimeric-dda-plus`, opt-in `--chimeric`).**
Beats Java on PSMs (Astral +101%, PXD +11%, all entrapment-FDP validated) and on speed
on all 3 datasets; **blocked from merge only by TMT PSMs (−5%)** per the beat-Java-on-both
gate. Reviewed (5-agent + 2 adversarial rounds), dead-code-cleaned, and GF/SpecE-parity
audited. See `docs/parity-analysis/notes/2026-05-31-cascade-optimized-multidataset-summary.md`.

TMT closing options (a future iteration, not the cascade): an **additive** Percolator
feature (e.g. DeltaRawScore) or a per-ion CID node-scoring trace. A native rescoring
pipeline (in-process Percolator + ms2pip/deeplc) is brainstormed in
`docs/superpowers/specs/2026-05-31-native-rescoring-pipeline-design.md` (4 open decisions).

**Abandoned:** the fragment-index candidate generator ("speed v2") was built and refuted
this session — Approach A degenerates (in-loop vote-all-touched), Approach B hits an
irreducible recall/speed tension (a top-K fragment prefilter drops exactly the secondary
co-isolated peptides that are the chimeric gain). Do not revisit fragment indexing for
chimeric speed; the cascade already beats Java on speed without it.
