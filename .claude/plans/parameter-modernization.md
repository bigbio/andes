# Plan: modernize MS-GF+ parameter handling

**Status: proposed**
Branch: `perf/search-sync-cleanup` (worktree at
`/Users/yperez/work/msgfplus-workspace/search-sync-cleanup`).

## Why this exists

The current parameter stack under `edu.ucsd.msjava.params` is doing
several jobs at once:
- command-line parsing
- type conversion
- validation
- help/usage rendering
- config-file alias handling
- backward-compatibility shims

That works, but it spreads option behavior across many small classes
(`Parameter`, `NumberParameter`, `RangeParameter`, `ToleranceParameter`,
`FileParameter`, enum wrappers, and `ParamManager`). The result is more
code than we need for a solved problem and a higher risk of subtle
parsing drift when new flags are added.

## Goals

- Reduce the amount of custom CLI parsing code.
- Keep existing MS-GF+ command-line behavior stable where practical.
- Preserve current config-file semantics in the first migration step.
- Keep `SearchParams` as the internal domain model for search settings.
- Improve help/usage generation and validation error consistency.

## Non-goals

- No search algorithm changes.
- No performance claim for the search itself; parsing happens once at
  startup and is not a runtime hotspot.
- No forced removal of legacy config-file aliases in phase 1.
- No broad package cleanup bundled into this effort.

## Recommended direction

Adopt `picocli` for command-line parsing and help generation, while
keeping a thin MSGF+-specific compatibility layer for:
- legacy option names and aliases
- config-file parsing
- repeated modification/custom-AA entries
- conversion into `SearchParams`, `AminoAcidSet`, `Tolerance`, and
  related domain objects

## Proposed migration shape

### Phase 1: introduce a typed CLI model beside `ParamManager`

- Add a new options class for `MSGFPlus` under `edu.ucsd.msjava.cli`.
- Represent flags as typed fields with defaults, required markers,
  and descriptions.
- Add custom `picocli` converters for:
  - precursor mass tolerance
  - integer and float ranges
  - output format
  - precursor calibration mode
  - file/directory validation
- Keep `ParamManager` intact during this phase.
- Add an adapter that maps parsed CLI options into the current
  `SearchParams` inputs.

Success criteria:
- `MSGFPlus` can parse its current CLI arguments through the new path.
- Generated help text is complete and readable.
- Existing tests for parameter behavior still pass or are updated
  mechanically where output formatting differs.

### Phase 2: preserve config-file compatibility explicitly

- Keep `ParamParser` or replace it with a thinner reader that still
  accepts the current `key=value` format.
- Centralize legacy config-name alias resolution in one place instead
  of scattering it through `ParamNameEnum`.
- Support repeated config entries for:
  - `DynamicMod`
  - `StaticMod`
  - `CustomAA`
- Feed config values into the same typed options model used by CLI.

Success criteria:
- Existing example parameter files still load.
- Duplicate-entry behavior for mods/custom amino acids is preserved.
- Command-line values continue to override config-file values.

### Phase 3: move validation out of the custom parameter hierarchy

- Replace per-type `parse()` methods with:
  - `picocli` conversion
  - explicit validation methods on the typed options object
  - targeted domain-level validation while building `SearchParams`
- Collapse or remove custom classes that are no longer needed:
  - `Parameter`
  - `NumberParameter`
  - `RangeParameter`
  - `IntParameter`
  - `FloatParameter`
  - `DoubleParameter`
  - `IntRangeParameter`
  - `FloatRangeParameter`
  - enum parameter wrappers

Success criteria:
- No user-visible behavior regressions on required flags, defaults,
  range checks, or enum choices.
- Validation failures still produce actionable messages.

### Phase 4: reduce `ParamManager` to compatibility-only or retire it

- If any remaining tools still depend on `ParamManager`, keep it only as
  a compatibility facade over the new parser.
- Otherwise remove `ParamManager` from the active CLI path.
- Decide whether `MSGFDB` migrates in the same PR series or follows
  after `MSGFPlus` is stable.

## Main risks

- Help text and error messages may change in ways that break tests or
  documentation.
- Config-file behavior is more important than it looks; it includes
  legacy aliases and repeated entries that generic CLI libraries do not
  model by default.
- `MSGFDB` and `MSGFPlus` share parts of the current stack, so an
  incomplete migration could increase duplication before it decreases.

## Validation plan

- Add focused tests for:
  - required arguments
  - default values
  - bad range syntax
  - enum parsing
  - file existence checks
  - config-file override precedence
  - repeated modification/custom-AA entries
- Keep existing `SearchParams` tests green.
- Run at least one end-to-end `MSGFPlus` smoke test on a known fixture.
- Compare old vs new parser outcomes for a representative set of real
  command lines and config files.

## Suggested implementation order

1. Add `picocli` dependency.
2. Build a typed `MSGFPlusOptions` class and converters.
3. Parse CLI into the new options class without removing `ParamManager`.
4. Add an adapter into the current `SearchParams` build path.
5. Port config-file handling.
6. Remove unused custom parameter classes.
7. Migrate `MSGFDB` only after `MSGFPlus` is stable.

## Recommendation on branch strategy

Do this in a dedicated refactor branch, not as part of a performance
cleanup PR. The expected win is maintainability and correctness, not
search throughput, and the surface area touches the public CLI.
