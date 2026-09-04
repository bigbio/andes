# MS-GF+ Rust port — incremental JNI-bridge alternative

This document is a planning artefact, not a commitment. It is the
**incremental, JNI-bridged alternative** to the full Rust rewrite
described in [rust-full-rewrite-roadmap.md](./rust-full-rewrite-roadmap.md).
Read the full-rewrite roadmap first for the comprehensive plan; this
document exists for the case where we decide a full rewrite is too
ambitious and we want a smaller-bet variant.

| Aspect | Full rewrite (sister doc) | Incremental JNI port (this doc) |
|---|---|---|
| Surface ported to Rust | All ~31 K LOC: domain model, FASTA, mzML, search, GF/DP, model loader, CLI, output writers | ~5 K LOC: SA walk + scorer + grid only |
| Java retained | Reference oracle during shadow phase only | Long-term: parsing, calibrator, CLI, Percolator output, all I/O |
| Crates | 8 (`msgf-core` … `msgf-cli`, `msgf-diff`) | 1 (`msgfplus-inner`, `cdylib` for JNI/FFM) |
| Boundary | Pure Rust binary, Java only as oracle | JNI/FFM call once per task; no per-spectrum boundary crossings |
| Validation | Differential `--engine diff` mode across full pipeline | Bit-identical Percolator-pin output via shadow mode `-engine compare` |
| Effort | Multi-quarter, multi-engineer | Multi-month, single engineer |
| Risk surface | High (mzML quirks, GF/DP semantics, modifications, decoys all need parity) | Lower — everything outside the inner loop stays Java; Rust faces only the hot, arithmetic-heavy code |
| Win | Single static binary, JVM-free distribution, full ownership of perf knobs | Hot-path perf wins, smaller maintenance footprint; JNI overhead is the cap |

The default position remains: **don't rewrite the whole thing**.
Java MS-GF+ is mature, the test fixtures cover real edge cases, and we
have a working performance roadmap (Phase B has already shipped −10 % on
Astral; Experiment 2 explored further). A Rust rewrite is justified only
if there's a specific outcome we can't reach by incremental Java work.

If a full rewrite is the right call, follow
[rust-full-rewrite-roadmap.md](./rust-full-rewrite-roadmap.md). If the
goal is just to claw back hot-path performance with a contained
maintenance footprint, the rest of this document describes that
alternative.

## When a rewrite is worth it

| Signal | Status |
|---|---|
| Java GC stalls dominate the wall on the hot path | No. PhaseE telemetry showed GC was a small fraction of Astral wall. Heap pressure exists but is bounded. |
| JNI/native dependencies already in the build | Yes — DuckDB pulled in on the trainer build. So we already ship a native lib. |
| We can't get further perf in Java without exotic tricks | Partially. Experiment 2 was the obvious next Java lever and didn't clear the gate. The catalog/fragment-index path was cleaner-on-paper but failed in practice. |
| Cold start matters | Yes, mildly. JVM startup + class load is ~2 s; a Rust binary would be ~50 ms. Not the headline but real for short runs and CI. |
| We want a single static binary for distribution | Yes if we ship cluster/cloud workloads. The Java JAR + JRE bootstrap is friction. |
| Memory footprint matters for very large databases | Maybe. The compact suffix array structures already fit in JVM heaps without trouble for the proteome scale we benchmark. Multi-billion-residue databases would change this. |

**Verdict:** the case for a rewrite is "real but not urgent." A focused
port of the hot path, keeping the I/O and orchestration in Java, is the
shape that makes sense — full-rewrite ROI is not there.

## Scope (incremental, JNI-bridged)

Port only the inner loop. Java keeps:

- `StaxMzMLParser` and the MGF parser — mzML I/O is solved.
- `BuildSA` and the on-disk suffix-array layout — already compact, used
  by other tools, format stability matters.
- `DirectPinWriter` and the Percolator pipeline — output formatting.
- `MassCalibrator` (Phase B) — pre-pass logic, runs once per file, not
  on the hot path.
- The CLI (`MSGFPlus.java`, `MSGFPlusOptions.java`) — argument parsing,
  orchestration.

Rust gets:

- `DBScanner.dbSearch` — the SA-walk that enumerates peptide candidates
  and dispatches scoring.
- `ScoredSpectraMap` lookup — the precursor-mass-window search and
  per-spectrum scoring buffer.
- `NewRankScorer.score(...)` — the inner score computation that reads
  the bundled `.param` rank distributions.
- `CandidatePeptideGrid` — the per-residue mass + variant tracking
  during the SA walk.

The new Rust crate exposes a single `run_inner_search(...)` C ABI
function. Java calls it via JNI per file; the boundary crosses ~once
per task, not per-spectrum.

## Architecture

```
┌─────────────── Java (bigbio/msgfplus) ──────────────┐
│                                                     │
│   StaxMzMLParser  ──► spectrum batch                │
│   BuildSA         ──► CompactSuffixArray on disk    │
│   MassCalibrator  ──► precursor shift + tightening  │
│                                                     │
│   ConcurrentMSGFPlus                                │
│        │                                            │
│        ▼                                            │
│   ┌─── JNI ──── msgfplus_inner ────────────────┐    │
│   │                                            │    │
│   │   load_param_file(.param)                  │    │
│   │   build_spec_map(spectra_blob)             │    │
│   │   walk_and_score(sa_blob, spec_map, opts)  │    │
│   │   ──► [(spec_idx, peptide, score), ...]    │    │
│   │                                            │    │
│   └────────────────────────────────────────────┘    │
│        │                                            │
│        ▼                                            │
│   DirectPinWriter  ──► .pin                         │
│                                                     │
└─────────────────────────────────────────────────────┘
```

Boundary contract: input/output buffers passed as `ByteBuffer`s sized
upstream. No per-spectrum or per-peptide JNI calls. The inner search
runs to completion and returns a flat result array.

## Validation contract

The port must produce bit-identical Percolator-pin output vs the Java
engine on every fixture in `src/test/resources/`. Specifically:

- `TestPrecursorCalIntegration.precursorCalOffMatchesBaseline` — the
  hardest gate; row-for-row pin equality on `test.mgf` +
  `human-uniprot-contaminants.fasta`.
- The 3-dataset Phase B validation set (Astral / TMT / PXD001819) —
  same target/decoy counts at OFF and AUTO.
- The ProteoBench HYE fasta + Astral mzML pair — same Percolator 1 %
  FDR target counts.

Strategy:

- Java float math is `strictfp`-equivalent in modern JVMs; Rust uses
  IEEE 754 by default. Use `f32` and `f64` consistently with Java's
  types, avoid SIMD reordering for the score sums until equality is
  proven, then opt in to faster paths behind a flag.
- Add a "shadow mode" CLI flag (`-engine compare`) that runs both Java
  and Rust paths and aborts on the first mismatch. Used in CI on the
  full fixture matrix.
- Keep `-engine java` as the production fall-back for at least one
  release after Rust becomes default.

## Milestones

Each is 1–3 weeks of focused work for one engineer. Skipping a milestone
means accepting tighter coupling and a less measurable rollout.

### M0 — Plan + skeleton

- Decide JNI vs `jextract`-style binding (Java 22+ has a path that may
  be cleaner). Pick one.
- Cargo workspace at `rust/` with one crate (`msgfplus-inner`),
  `cdylib` target, no functional code yet.
- Build hookup in `pom.xml` so `mvn package` builds the Rust lib for
  the host platform and shades it into the JAR under
  `META-INF/native/`. Mirror the DuckDB pattern.
- A minimal JNI smoke test that calls a `version()` function across
  the boundary and asserts the returned string. Same shape as
  `TestDuckDbParquetSmoke`.

### M1 — Port `NewRankScorer`

- Rust loader for the bundled `.param` format.
- Rust `score(prm_grid, nominal_grid, num_mods, charge)` — pure
  function, no allocation in the hot path.
- Equality test: drive the Java and Rust scorers on the same
  hand-crafted `(grid, charge, num_mods)` tuples sampled from real
  search runs; assert byte-exact equality of every score float.

### M2 — Port `CompactSuffixArray` reader + walker

- Read the on-disk SA layout into Rust slices. No re-sort, no rebuild.
- Walker that emits `(start, length, peptide_bytes)` for each peptide
  the Java SA walk currently visits, in the same order.
- Equality test: dump the visit-order trace from Java and Rust on
  ecoli.fasta and assert identical sequences.

### M3 — Port `CandidatePeptideGrid`

- Per-residue mass + variant tracking. Stateful but bounded —
  comfortably fits the Rust ownership model.
- Equality test: drive both grids through the same residue+mod
  sequence and compare per-position outputs.

### M4 — Inner-search integration

- Combine M1–M3 into `walk_and_score(...)`. Wire through JNI.
- Run the full Phase B 3-dataset validation under `-engine compare`.
  The bar: zero pin-row drift, walls strictly faster than Java baseline.

### M5 — Production rollout

- Default to `-engine rust` for the next release.
- Keep `-engine java` opt-in for one release cycle.
- Remove the Java inner-loop after a clean release passes.

## Cross-platform packaging

Mirror the DuckDB precedent already in the repo: native `.so` / `.dylib`
/ `.dll` files for `linux_amd64`, `linux_arm64`, `osx_universal`, and
`windows_amd64` in `META-INF/native/`. CI matrix builds one `.jar`
artefact per host, OR a single fat jar with all four artefacts (the
DuckDB jar is the example — it embeds all four in one 73 MB blob).

Cross-compilation strategy: GitHub Actions matrix with `ubuntu-latest`,
`ubuntu-24.04-arm`, `macos-14` (Apple silicon), `windows-latest`. Each
job builds the Rust crate for its target, the final assembly job
gathers the four `.so/.dylib/.dll` outputs and runs `mvn package`.

## Open decisions before starting

1. **Java version floor.** The repo is on Java 17. JEP 442 (Foreign
   Function & Memory API) is preview in 21 and standard in 22 — using
   it means a hard Java-22 floor. Decide if we keep JNI or jump to FFM.
2. **`.param` model on the Rust side.** Do we re-read the binary
   `.param` files in Rust, or convert them to a Rust-native format at
   build time? Re-reading is one-time; converting buys us a chance to
   simplify the file format (the Halobacterium pilot showed the
   binary serialisation has some fragility).
3. **MSRV.** Pin a minimum supported Rust version in
   `rust-toolchain.toml`. Suggest 1.75 (current stable as of writing).
4. **Build trigger.** Default `mvn package` builds the Rust crate or
   skips it? Skipping by default keeps Java-only contributors
   unblocked but hides drift; building always is heavier.
5. **Telemetry continuity.** Rust path should emit the same
   `[Phase B telemetry]` and `[Experiment 2 telemetry]` log lines (if
   we add Experiment 2 back) so existing benchmark scripts keep
   working. Decide log format up front.

## What this plan is *not*

- Not a commitment to start. The bench-time we'd spend on M0 is real;
  the Java engine is performant enough today that this is a
  multi-quarter project, not a quarter project.
- Not a claim of speedup. M4 has a "strictly faster" gate; if we
  can't beat the Java baseline by enough to justify the maintenance
  cost of two engines, we abandon and document the lesson.
- Not a commitment to drop Java. The CLI and orchestration stay even
  after M5.

## Reference points in the existing codebase

The hot path that would be ported lives in:

- [DBScanner.java](../src/main/java/edu/ucsd/msjava/msdbsearch/DBScanner.java)
  — `dbSearch` is the SA walk; this is the function the Rust crate
  ultimately replaces.
- [ScoredSpectraMap.java](../src/main/java/edu/ucsd/msjava/msdbsearch/ScoredSpectraMap.java)
  — `getPepMassSpecKeyMap`, `preProcessSpectra`, the `subMap` lookup
  used in pairing.
- [NewRankScorer.java](../src/main/java/edu/ucsd/msjava/msscorer/NewRankScorer.java)
  — `.param` reader + `score(...)`.
- [CompactSuffixArray.java](../src/main/java/edu/ucsd/msjava/msdbsearch/CompactSuffixArray.java)
  — on-disk SA format; the Rust port reads the same bytes.
- [CandidatePeptideGrid.java](../src/main/java/edu/ucsd/msjava/msdbsearch/CandidatePeptideGrid.java)
  — per-residue mass + variant tracker.

Read those before M0 to size the actual Java surface accurately.
