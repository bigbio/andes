# msgf-rust port — design spec

**Status:** approved by user 2026-05-03 (brainstorm session).
**Branch:** `rust-implement` (forked clean from `dev`).
**Output of:** brainstorming Q1-Q11 + 6 design sections; consolidates the
two predecessor planning documents into a single actionable spec.

## TL;DR

`msgf-rust` is a from-scratch Rust implementation of MS-GF+ that fully
replaces the Java engine for supported workflows. End-state shipping
shape is Sage-style: single static binary on GitHub Releases per
platform, no JVM dependency.

The work lives on `rust-implement` as a long-lived branch; the merge to
`dev` is one literal PR carrying the entire rewrite. The merge gate is
1 % FDR target count parity (±0.5 %) on Astral + TMT + PXD001819 vs the
Java oracle.

## Goals

1. Replace Java MS-GF+ end-to-end for supported workflows (HCD/Tryp at
   minimum; additional combinations as parity allows).
2. Preserve scientific behaviour: per-PSM `SpecEValue` parity within
   tolerance (`|java − rust| / java ≤ 1e-3`); 1 % FDR target counts
   match Java baseline within ±0.5 % on the three sign-off datasets.
3. Ship a single static binary per platform — no JVM, no Maven, no
   shaded JAR.
4. Read existing Java `.param` files directly so the bundled scoring
   models and the recently-shipped msnet trainer pipeline (PR #27) stay
   relevant without migration.
5. Drop-in CLI compatibility for existing scripts and integrators
   (quantms): `msgf-rust -s file.mzML -d db.fasta -o out.pin ...` works
   exactly as `java -jar MSGFPlus.jar -s ...` does today.

## Non-goals

- *Not* aiming for a specific speedup. Match Java's wall on the three
  sign-off datasets at minimum; perf wins are bonus.
- *Not* adding new search features. Behaviour-preserving port only.
- *Not* introducing a new `.param` format. Read Java's existing binary
  layout directly.
- *Not* maintaining the Maven distribution post-merge. Existing Maven
  consumers (quantms) migrate to the binary download.
- *Not* publishing the workspace crates to crates.io. The only
  user-facing artefact is the `msgf-rust` binary.
- *Not* porting mzID writers. mzID was excised from this fork in
  commit `9bf01c8`; that decision stands.

## Architecture

### Workspace layout

`rust/` directory at the repo root, parallel to the existing Java
`src/`. Polyglot repo, single source of truth, atomic commits when
Java + Rust touch the same area.

```
astral-speed/
├── pom.xml                    ← Java unchanged
├── src/                       ← Java unchanged
├── rust/
│   ├── Cargo.toml             ← workspace root
│   ├── rust-toolchain.toml    ← MSRV pin (1.80)
│   ├── crates/
│   │   ├── core/              ← AAs, masses, mods, peptides, enzymes, tolerances
│   │   ├── spectra/           ← MGF + mzML readers (uses mzdata)
│   │   ├── fasta/             ← FASTA load, decoy gen, suffix array
│   │   ├── search/            ← Candidate generation, precursor matching, ranking
│   │   ├── gf/                ← Generating function / DP (highest-risk)
│   │   ├── model/             ← .param binary loader
│   │   ├── cli/               ← The binary (`name = "msgf-rust"`)
│   │   └── diff/              ← Parity tool (`name = "msgf-diff"`)
│   └── tests/parity/          ← Workspace integration fixtures
├── docs/                      ← Rust-related docs land here
└── .claude/plans/             ← Predecessor roadmaps
```

### Crate naming

Bare names — `core`, `spectra`, `fasta`, etc., not `msgf-` prefixed.
The crates are workspace-internal only; the only user-facing artefact
is the binary `msgf-rust`. These names are too generic for crates.io
publication, which is fine — the chosen distribution model
(GitHub Releases, single binary) doesn't involve crates.io.

### Crate dependencies

```
core ──────────────────────┐
   │                       │
   ↓                       │
fasta            model     │
   │             │         │
   └────→ search ←─────────┘
              │
              ↓
            gf

spectra ──→ uses mzdata

cli ──→ depends on core, spectra, fasta, search, gf, model
diff ──→ standalone (only depends on .pin reading)
```

### External dependencies

| Crate | Role | Justification |
|---|---|---|
| `mzdata` | mzML + MGF + MzMLb readers | Active proteomics-Rust ecosystem; gated on parsing-perf comparison vs Java in M0. If >30 % slower than `StaxMzMLParser`, fall back to a hand-rolled subset reader. |
| `rayon` | Data-parallel scheduler | Sage's choice; matches Java `ForkJoinPool` semantics. |
| `dashmap` (rayon feature) | Concurrent hashmaps | Where shared mutable state is unavoidable in parallel loops. May not be needed if pre-processing is sequential; decide in M0. |
| `clap` | CLI parsing | Supports legacy-style flags + modern subcommands in one tree via `clap::Parser`. |
| `byteorder` | Big-endian Java `DataInputStream` parity | Java's `.param` writer uses `DataOutputStream` (big-endian). |
| `thiserror` | Error types | No panics in library code. |
| `tracing` | Structured logging | Replaces logback usage. |

### MSRV

Pinned to **1.80** in `rust-toolchain.toml`. Reason: stable
`LazyLock`/`LazyCell` (replaces `once_cell` patterns), `let`-`else`,
mature iterator updates that simplify GF/DP transcription. Bump in M0
if specific deps need newer.

## Compatibility surface

### Two binaries

- `msgf-rust` — main user-facing CLI (`crates/cli/`).
- `msgf-diff` — parity comparison tool used by CI + nightly bench cron
  (`crates/diff/`).

### `msgf-rust` CLI — dual: drop-in legacy + Rust-native subcommands

Legacy form (drop-in for existing scripts and quantms) — detected by
absence of a known subcommand:

```sh
msgf-rust -s file.mzML -d db.fasta -mod Mods.txt -o out.pin \
          -tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 \
          -ntt 2 -minLength 6 -maxLength 40 -addFeatures 1 \
          -thread 8 -precursorCal auto
```

Rust-native subcommand form (recommended for new users):

```sh
msgf-rust search \
   --spectra file.mzML --database db.fasta --modifications Mods.txt --output out.pin \
   --target-decoy --tolerance 10ppm --isotope-error -1,2 \
   --activation HCD --instrument QExactive --enzyme Tryp --protocol Standard \
   --num-tolerable-termini 2 --min-length 6 --max-length 40 \
   --add-features --threads 8 --precursor-cal auto
```

Implementation: a single `clap::Parser` with a top-level dispatch.
First positional matched against known subcommand names (`search`,
`info`, `version`); fall through to legacy parser otherwise. Single
binary, two parsing paths.

### Subcommand inventory (v1)

- `search` — main search; drop-in equivalent for `java -jar MSGFPlus.jar -s ...`.
- `info` — inspect a `.param`, `.mzML`, or `.fasta` (debug aid for dev + power users).
- `version` — print version + build metadata (Rust version, git sha, mzdata version, .param schema rev).

(`train` is *not* a v1 subcommand. The Java trainer in PR #27 retains
that responsibility for v1; future Rust trainer integration has a slot
reserved.)

### Output formats — bit-equivalent to Java

- `.pin` (Percolator input): same column ordering as `DirectPinWriter`;
  same `Label` values (1 / -1); same field types; same `%.6f`/`%.6e`
  precision so downstream Percolator wrappers see byte-identical text
  (modulo legitimate float-order drift, which is bounded by the
  per-PSM SpecEValue tolerance gate).
- `.tsv` (Direct TSV output): same columns as `DirectTSVWriter` —
  `#SpecFile`, `SpecID`, `ScanNum`, `Title` (when MGF), `FragMethod`,
  `Precursor`, `IsotopeError`, `PrecursorError(ppm|Da)`, `Charge`,
  `Peptide`, `Protein`, `DeNovoScore`, `MSGFScore`, `SpecEValue`,
  `EValue`, [`QValue`, `PepQValue` when `-tda 1`], plus the
  `addFeatures` columns.
- mzID writing is *not* supported (matches the bigbio fork's existing
  position).

### `.param` bundling and lookup

- The bundled set of 31 `.param` files in
  `src/main/resources/ionstat/` gets embedded into the `msgf-rust`
  binary via `include_bytes!()` in `crates/model/`. Same lookup pattern
  as Java's `NewScorerFactory`: filename =
  `<activation>_<instrument>_<enzyme>[_<protocol>].param`.
- No CLI flag for a custom `.param` path in v1 (Java doesn't have one
  either; users swap the bundled file at build time).
- Binary format read directly from Java's `DataOutputStream` layout via
  `byteorder` crate; `crates/model/` exposes
  `Param::load_from_bytes(&[u8]) -> Result<Param>` mirroring
  `NewRankScorer.read()`.

### Error / exit code surface

Match Java's behaviour: exit 0 on success, non-zero on error.
All errors print to stderr with a clear message; no stack traces in
normal output. `--verbose` enables `tracing` debug output.

## Concurrency model

### Parallelism unit

Spectra. Same partitioning shape as Java's `ConcurrentMSGFPlus` —
split the spectrum list into 3× num-threads chunks (the existing
heuristic), schedule chunks concurrently via Rayon, merge results at
the end.

```rust
use rayon::prelude::*;

let results: Vec<MatchedPsm> = spectra
    .par_chunks(chunk_size)
    .flat_map(|chunk| score_chunk(chunk, &sa, &model, &params))
    .collect();
```

`score_chunk` is pure-CPU, takes everything by reference, allocates
only its own per-chunk buffers. Rayon's work-stealing handles
imbalance across chunks (some peptides match many SpecKeys, some
match few — same imbalance Java's `ForkJoinPool` already absorbs).

### State sharing

- *Read-only / shared by reference across threads:* `CompactSuffixArray`,
  `Param`, `ScoredSpectraMap` (after preprocessing), `AminoAcidSet`.
  All loaded once, immutable thereafter, passed as `&Sa`, `&Param`, etc.
  Zero synchronization overhead.
- *Per-chunk thread-local:* result buffer (`Vec<MatchedPsm>`),
  `CandidatePeptideGrid`, scoring scratch space. Allocated per
  `score_chunk` call, freed when the chunk finishes. No `Mutex`,
  no `Arc<RwLock>`.
- *Concurrent shared (rare path):* `dashmap::DashMap` for any case
  where multiple threads need to write to the same map. Currently
  anticipated only for the spectrum-mass index during pre-processing,
  before the parallel scoring loop starts. If pre-processing turns out
  to be sequential, `DashMap` may not be needed; decide in M0.

### Determinism for the parity gate

Spectrum-major output. Per-chunk results returned in chunk-input order
via `collect()`; within a chunk, scoring is deterministic by spectrum
index. Final `.pin` output ordered by spectrum index, identical to
Java.

### Spectrum mutation invariant

Rust enforces what `.claude/CLAUDE.md` warns about for Java: the
parallel iterator yields `&Spectrum`, not `&mut Spectrum`. Java's
manual `setCharge` / `getScoredSpectrum` patterns become impossible
here without an explicit `unsafe` block, which we won't write.
Whole class of bugs (the calibration pre-pass mutation we hit in
Phase B) becomes uninvitable.

## Validation stack

### `msgf-diff` tool

The workhorse of every gate. Reads two `.pin` files, joins rows on
`(SpecID, ScanNum, Charge, Peptide)`, compares each field with a
per-field tolerance config, exits non-zero if any tolerance exceeded.

```sh
msgf-diff compare java.pin rust.pin \
    --tolerance "SpecEValue:1e-3,EValue:1e-3,MSGFScore:1.0" \
    --report drift.json
```

Output: histogram of drift per field, top-N largest mismatches,
summary counts (rows-only-in-A, rows-only-in-B, rows-with-drift,
rows-bit-identical).

### Three-tier validation

| Tier | When | Where | Inputs | Gate |
|---|---|---|---|---|
| Per-push | every commit to `rust-implement` | GitHub Actions | `test.mgf` (small) + `BSA.fasta` | Per-PSM `\|java.specEValue − rust.specEValue\| / java.specEValue ≤ 1e-3`. Wall ~1-2 min. |
| Nightly | cron at `pride-linux-vm.ebi.ac.uk` | bench machine | `iprg-2013/F13.mgf` + `ecoli.fasta`, plus the 3 sign-off datasets at `-precursorCal off` | Same per-PSM tolerance + raw target/decoy count parity. Posts results to `STATUS.md` on the branch. Wall ~3h. |
| Sign-off | manual, pre-merge to `dev` | bench machine | Astral + TMT + PXD001819, both `precursorCal off` and `auto` | 1 % FDR target count parity ±0.5 % after Percolator. The user-facing metric. |

### GF/DP transcribe-then-optimize loop

Phase 6 (the GF/DP port) runs as a two-stage development loop:

1. **Stage 1 — direct transcription.** Java `msgf/` package goes
   line-by-line into `crates/gf/`. Same variable names, same control
   flow. `HashMap<K, V>` becomes `std::collections::HashMap<K, V>`,
   not `BTreeMap`. `int[]` becomes `Vec<i32>`, not `&[i32]`. Per-PSM
   SpecEValue tolerance must be ≤ 1e-3 on the small fixture before
   this stage merges.
2. **Stage 2 — optimize incrementally with parity brake.** Each
   optimization PR (within the branch) targets a specific hotspot:
   - Replace a `HashMap` with a `Vec`-backed dense lookup
   - Inline a hot-loop helper
   - Reuse a buffer across iterations
   - Use SIMD where fragment-ion masses align

   Each such PR runs the per-push gate. If tolerance breaks, the PR
   is reverted and a different angle tried. The `1e-3` ceiling is the
   hard brake.

The roadmap's "preserve behavior before optimizing" principle is
enforced operationally by this loop.

### `mzdata` performance gate (M0)

Before committing to `mzdata`, run a comparison harness:

```
benchmark/parse-mzdata-vs-staxmzmlparser.sh
```

Parses the Astral mzML through both engines, compares:
- Wall time (must be within 30 %)
- Spectrum count (must match exactly)
- Per-spectrum metadata: precursor m/z, charge, scan number, retention
  time (must match exactly)

Pass → ship with `mzdata`. Fail → roll our own subset reader for the
`mzML 1.1` indexed-centroid path (covers all 3 sign-off datasets).

### Test fixtures

| Tier | Fixture | Spectra | Wall (target) |
|---|---|---:|---:|
| Per-push | `test.mgf` + `BSA.fasta` (existing) | ~few | <30s round-trip |
| Per-push extended | `iprg-2013/F13.mgf` + `ecoli.fasta` (existing) | 1,406 | <2 min round-trip |
| Nightly + sign-off | Astral + TMT + PXD001819 (existing remote bench setup) | 50K-150K each | ~10-30 min/dataset |

All fixtures already exist in this repo or on the bench machine — no
new data acquisition needed.

## Phasing on `rust-implement`

| Phase | What | Per-PSM gate fixture | Wall (Claude-assisted) |
|---|---|---|---:|
| **M0** | Workspace + `rust-toolchain.toml` + `Cargo.toml` + clap skeleton + `msgf-diff` v0 + Java reference-capture harness + **`mzdata` perf gate run** | n/a | 5-7 days |
| **1 — `core`** | AAs, masses, mods, peptides, enzymes, tolerances, protocols. No I/O. | unit-level parity vs Java on hardcoded value tables | 5-7 days |
| **2 — `model`** | `.param` binary loader. Reads existing Java format directly. | byte-exact parity loading every bundled `.param` | 3-5 days |
| **3 — `spectra`** | mzML + MGF readers via `mzdata` (or fallback per M0 gate). | per-spectrum metadata parity on `test.mgf` + Astral mzML | 2-3 weeks |
| **4 — `fasta` + `search`** | FASTA load, decoy gen, suffix array reader, candidate generation, precursor matching, ranking queues. | candidate-set parity on traced spectra; same raw target/decoy pre-GF candidate counts | 2-3 weeks |
| **5 — cheap scoring** | The scorer path that evaluates PSMs before GF. | per-spectrum top-N candidate identity match | 5-10 days |
| **6 — `gf`** | Generating function / score-distribution DP. Highest correctness risk. Stage 1 direct transcription, Stage 2 optimize-with-brake. | per-PSM SpecEValue tolerance ≤ 1e-3 across `test.mgf`, then `F13.mgf`, then nightly cron on Astral | 5-8 weeks |
| **7 — output writers** | `.pin` + `.tsv` formatters in `cli`. | byte-equivalent output on `test.mgf` (modulo float-text precision allowed by per-PSM tolerance) | 5 days |
| **8 — shadow validation** | Long tail: edge cases, vendor quirks, mzML variants, mod-expansion subtleties. | nightly cron must show zero tolerance breaches over 7 consecutive nights before this phase is done | 3-6 weeks |
| **9 — production cutover** | Merge `rust-implement → dev`. Java becomes `--engine java` deprecated fallback. | sign-off gate (1 % FDR target counts ±0.5 % on Astral + TMT + PXD001819) is the merge gate. | 1.5-3 weeks |

**Total: 14-26 weeks (~3.5-6 months) with Claude assistance.**

**Phase parallelism note:** phases 2 (`model`) and 3 (`spectra`) are
independent — both only depend on `core`. They can be developed in
parallel by either two engineers or two Claude-assisted sessions.
The ordering above reflects dependency, not required serialisation.
Phases 4 onward depend on phase 3 (`search` reads spectra) so the
critical path runs through `core → spectra → fasta+search → gf`.

## Java overlap window

- *During phases 1-8:* Java engine on `dev` continues normal evolution.
  `rust-implement` rebases weekly to absorb. Java is the oracle.
- *At phase 9 merge:* msgf-rust ships as the default engine. Java
  engine renamed in CLI from default to `--engine java`; bundled JAR
  continues to be built by `mvn package` for one full release cycle.
- *One release cycle later (post-merge):* separate cleanup PR removes
  the Java engine source tree, the JAR build, and the `--engine java`
  flag. msgf-rust becomes the only implementation.

The release cycle gap is intentional: if msgf-rust hits a regression
in production we haven't caught, users have a fallback path.

## Telemetry continuity

msgf-rust emits the same `[Phase B telemetry]` log lines (and any
future telemetry) so existing benchmark scripts at
`pride-linux-vm.ebi.ac.uk` keep parsing them. CLI `--verbose`
controls emission, same way `-Dmsgfplus.phaseBTelemetry=true`
controls it in Java today.

## Open decisions deferred to M0

- Exact `mzdata` vs roll-our-own decision (gated on the M0 perf
  comparison).
- Whether `dashmap` is needed at all (depends on whether pre-processing
  ends up being parallelisable).
- Whether to bump MSRV beyond 1.80 (depends on what specific deps
  pull in).

## References

- `.claude/plans/rust-full-rewrite-roadmap.md` — predecessor full
  roadmap. Phases and validation strategy that informed this spec.
- `.claude/plans/rust-incremental-jni-alternative.md` — predecessor
  alternative scoping (kept for posterity; rejected by user in favour
  of full rewrite).
- `.claude/CLAUDE.md` — project-level invariants the port must respect
  (no shared `Spectrum` mutation, hot-path collection bans, etc.).
- Sage (https://github.com/lazear/sage) — closest comparable; pure-Rust
  proteomics search engine, single-static-binary distribution, Rayon-based
  parallelism. Threading model and distribution shape adopted directly.

## Change log

- 2026-05-03 — initial spec, brainstormed in session.
