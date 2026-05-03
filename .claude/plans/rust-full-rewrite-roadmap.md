# Plan: full Rust rewrite of MS-GF+ search engine

**Status: proposed**
Branch: `feat/rust-rewrite-plan` (worktree at
`/Users/yperez/work/msgfplus-workspace/astral-speed`).

## Why this exists

We want a realistic plan for moving the current Java implementation to
Rust without pretending this is a normal refactor. For this codebase,
the main risk is not compilation effort. It is silent scientific drift:
different candidate sets, score distributions, generating-function
behavior, or target/decoy outcomes that still look superficially
plausible.

Current codebase shape, as of 2026-05-02:
- `156` Java source files under `src/main/java`
- about `31k` lines of Java in the main source tree
- `34` Java test files
- largest engine areas:
  - `msdbsearch` ~ `6.7k` lines
  - `msscorer` ~ `3.5k` lines
  - `msgf` ~ `3.5k` lines

That means a full rewrite is possible, but it should be managed as a
new implementation with long-running differential validation against the
existing Java engine.

## Goals

- Reimplement the production search engine in Rust.
- Preserve scientific behavior closely enough that the Rust engine can
  replace the Java engine for supported workflows.
- Keep a staged path so we can compare Java and Rust outputs at each
  boundary instead of discovering drift only at the end.
- End with a Rust-native CLI and model-loading path that can ship
  independently of the JVM.

## Non-goals

- No assumption that a full rewrite automatically delivers a large speed
  win.
- No attempt to preserve every historical CLI quirk on day 1.
- No requirement to port every ancillary tool before the core search
  path is trustworthy.
- No “big bang” cutover where Java is removed before differential
  testing is complete.

## Recommendation

If we commit to a full rewrite, we should still execute it in phases:
- first build a Rust engine with a narrow, testable scope
- then run it in shadow mode against Java
- only then promote Rust to the default implementation

This is a better fit for the codebase than trying to port everything in
one pass and hoping the final answers line up.

## Rewrite principles

1. Preserve behavior before optimizing.
2. Treat every interface boundary as a test boundary.
3. Prefer exact replay of Java semantics over “cleaner” Rust semantics
   until parity is proven.
4. Make unsupported features explicit instead of silently approximating.
5. Keep Java as the oracle until Rust clears promotion gates.

## Architecture target

### Rust ownership

The final Rust engine should own:
- FASTA loading and indexing
- spectrum parsing for supported formats
- peptide candidate generation
- modification expansion
- cheap scoring
- generating-function / score-distribution computation
- ranking and result assembly
- output writing for supported formats

### Java overlap period

During the migration, Java remains the reference implementation for:
- parity comparisons
- golden-output generation
- temporary fallback for unsupported workflows

### Suggested Rust crate split

- `msgf-core`
  - amino acids, peptides, masses, tolerances, enzymes, protocols
- `msgf-spectra`
  - MGF and mzML readers, spectrum data model
- `msgf-fasta`
  - FASTA reading, decoy handling, sequence/index utilities
- `msgf-search`
  - candidate generation, precursor matching, ranking, result objects
- `msgf-gf`
  - generating functions, score distributions, DP machinery
- `msgf-model`
  - `.param` loading, scorer model types, future trainer hooks
- `msgf-cli`
  - user-facing command
- `msgf-diff`
  - differential-test and trace-comparison tools against Java outputs

## Migration phases

### Phase 0: freeze reference behavior

Before porting, define what “correct” means.

Deliverables:
- fixed benchmark datasets for:
  - small smoke test
  - standard proteomics
  - Astral
  - TMT
- captured Java reference outputs for each dataset
- explicit comparison rules for:
  - exact equality
  - float tolerance
  - acceptable ranking tie behavior

Reference artifacts should include:
- raw target count
- raw decoy count
- 1 % FDR targets
- 5 % FDR targets
- top PSM identities
- `SpecEValue` distribution summary
- selected per-spectrum traces for hard edge cases

Success criteria:
- every future Rust milestone compares against the same frozen oracle

### Phase 1: port the domain model only

Implement Rust equivalents of the stable low-level concepts:
- amino acids
- modifications
- peptide representation
- mass/tolerance math
- enzyme rules
- protocol / activation / instrument enums

This phase should not perform search yet.

High-risk details to preserve:
- nominal mass vs accurate mass semantics
- asymmetric precursor tolerances
- isotope error conventions
- peptide formatting with inline modifications
- enzyme cleavage edge cases

Success criteria:
- unit parity on low-level calculations
- exact reproduction of representative peptide and mass examples

### Phase 2: port model loading

Implement `.param` reading in Rust using the current Java binary format.

Source compatibility target:
- files currently loaded by [NewScorerFactory.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msscorer/NewScorerFactory.java:82)
- binary format currently written by [NewRankScorer.java](/Users/yperez/work/msgfplus-workspace/astral-speed/src/main/java/edu/ucsd/msjava/msscorer/NewRankScorer.java:663)

Do not redesign the format first. Load the current format exactly.

Success criteria:
- Rust can load all shipped `.param` files
- parsed metadata matches Java
- sampled partition and frequency tables match Java decoding

### Phase 3: port spectrum parsing

Implement Rust readers for the supported search formats first:
- MGF
- mzML

Preserve these behaviors:
- spectrum ID handling
- scan numbering
- precursor charge/mz extraction
- activation-method handling
- centroid/profile filtering semantics

Recommended rule:
- MGF first, mzML second
- no vendor-quirk expansion until current test datasets pass

Success criteria:
- Rust and Java produce the same logical spectrum metadata on frozen
  datasets
- supported-format declaration is explicit

### Phase 4: port candidate generation and precursor matching

This is the first engine phase where correctness risk becomes high.

Implement:
- FASTA traversal / indexing
- decoy handling
- enzyme-aware peptide boundary logic
- variable-mod expansion
- precursor-window matching
- per-spectrum candidate collection and ranking queues

This phase should stop before GF / final `SpecEValue` if needed, so we
can compare candidate sets first.

Success criteria:
- same candidate peptide identities for traced spectra
- same charge assignment behavior
- same raw target/decoy pre-GF candidate counts on fixed samples

### Phase 5: port cheap scoring

Implement the scorer path that evaluates peptide-spectrum candidates
before GF.

Risk areas:
- rank-based scorer interpretation
- ion-type lookup semantics
- spectrum preprocessing assumptions
- tie handling in bounded queues

Success criteria:
- candidate ordering matches Java on traced spectra
- per-spectrum top-N candidate sets are equal or explainably tied

### Phase 6: port generating-function / DP engine

This is the most correctness-sensitive phase.

Implement:
- graph construction
- partition-dependent scoring
- score-distribution DP
- spectral probability / `SpecEValue` computation

This phase should be treated as a standalone parity project.

Validation should compare:
- nominal mass windows
- score distribution support
- per-candidate spectral probabilities
- final rank ordering

Success criteria:
- no systematic `SpecEValue` drift on reference datasets
- native target/decoy and FDR summaries line up with Java

### Phase 7: port outputs and FDR-facing behavior

Implement Rust output for the formats we actively support:
- `.tsv`
- `.pin`

Preserve:
- peptide formatting
- protein accession formatting
- rank numbering
- field naming expected by downstream tools

Recommendation:
- keep FDR comparison external at first if needed
- then port the in-process result labeling and q-value surface

Success criteria:
- output files are consumable by current downstream workflows
- schema-compatible with present Java outputs

### Phase 8: shadow mode

Run Rust beside Java on the same datasets and compare results.

Suggested modes:
- `--engine java`
- `--engine rust`
- `--engine diff`

`diff` mode should compare:
- result counts
- top hits
- `SpecEValue`
- output rows
- optional candidate traces for debug builds

Success criteria:
- Rust clears all parity gates on supported workflows
- remaining mismatches are explained and documented

### Phase 9: production cutover

Promote Rust only for workflows that have passed shadow validation.

Recommended launch order:
1. MGF + standard tryptic search
2. mzML + standard tryptic search
3. Astral-focused supported workflow
4. TMT and broader protocol variants

Java should remain available as a fallback until the Rust path has
cleared a full release cycle.

## Validation strategy

### Exact-equality targets

These should be exact whenever practical:
- peptide sequence formatting
- protein accession formatting
- target/decoy labeling
- raw target and decoy counts
- top-hit identity for traced spectra

### Tolerance-based targets

These may need explicit float tolerances:
- internal score components
- probability values
- error summaries

### Golden dataset set

Keep at least four:
- tiny deterministic fixture
- mainstream HCD tryptic fixture
- Astral fixture
- TMT fixture

### Differential tooling

Build comparison tools early, not at the end.

Needed comparators:
- spectrum metadata diff
- candidate-set diff
- top-N scorer diff
- final output diff
- summary-metric diff

## Main risk areas

### Highest scientific risk

- generating-function / DP behavior
- modification expansion semantics
- peptide-boundary and enzyme logic
- rank/tie handling in bounded queues
- precursor tolerance and isotope-error handling

### Highest engineering risk

- mzML compatibility
- memory layout for large datasets
- replicating Java object-graph behavior with value-oriented Rust types
- preserving all supported workflow combinations without a giant if-tree

### Highest product risk

- spending months to reach a Rust engine that is faster but not yet
  trustworthy
- underestimating the time needed for parity tooling
- letting unsupported workflows fail silently instead of clearly

## Subsystem difficulty map

### Easier first ports

- enums and domain types
- tolerance math
- `.param` binary reader
- TSV / PIN writing

### Medium difficulty

- FASTA ingestion
- MGF parsing
- candidate queue logic
- direct scorer path

### Hard

- mzML behavior parity
- modification expansion semantics
- protein mapping / accession reconstruction

### Hardest

- generating-function engine
- final search parity on real datasets

## Staffing assumptions

This is best executed with:
- one owner for core search semantics
- one owner for parser / IO / CLI
- one owner or recurring effort on validation infrastructure

Single-person implementation is possible, but the validation burden will
stretch the timeline significantly.

## Suggested milestone order

1. Freeze Java golden datasets and parity rules.
2. Port domain model.
3. Port `.param` loader.
4. Port MGF parser.
5. Port FASTA + candidate generation.
6. Compare candidate sets before scoring.
7. Port cheap scoring.
8. Port GF / DP engine.
9. Port mzML.
10. Add shadow `diff` mode.
11. Port output writers.
12. Promote Rust workflow by workflow.

## Recommendation on branch strategy

Do not do this as one long-lived feature branch only.

Recommended structure:
- planning branch: `feat/rust-rewrite-plan`
- implementation branches by milestone:
  - `feat/rust-domain-model`
  - `feat/rust-param-loader`
  - `feat/rust-mgf-parser`
  - `feat/rust-candidate-engine`
  - `feat/rust-gf-engine`
  - `feat/rust-shadow-diff`

That keeps reviewable scope and makes it easier to pause after any
parity failure.

## Final recommendation

A full Rust rewrite is feasible, but it should be approved as a
multi-stage reimplementation project with explicit parity gates, not as
an optimistic “performance refactor.”

If at any point we cannot keep parity around candidate generation or the
GF engine, we should be willing to stop the rewrite from becoming the
default engine until those gaps are closed.
