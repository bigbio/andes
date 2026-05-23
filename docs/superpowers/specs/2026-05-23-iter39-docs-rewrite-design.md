# iter39 — docs rewrite + CLI rename for the post-cutover state

**Branch:** `iter39-docs-rewrite` (cut from `master` HEAD `c863dae1`)
**Date:** 2026-05-23
**Status:** design approved, plan pending

---

## Context

PR #29 landed the Rust port of MS-GF+ as the production engine. The repo was
de-forked from `MSGFPlus/msgfplus` and renamed `bigbio/msgfplus` →
`bigbio/msgf-rust`. The Rust workspace is now at the repo root
(`Cargo.toml`, `crates/`, `resources/`, `test-fixtures/`). The Rust port beats
Java MS-GF+ at 1% FDR on all three benchmark datasets (Astral +0.98%,
PXD001819 within 0.3% at 3.3× wall, TMT +9.3% at 14% faster wall).

The current `README.md` and `docs/` tree predate the cutover. They describe
the Java tool: `mvn` build, JAR distribution, Java CLI flags, Java parameter
file templates. Most of it is stale.

This iteration treats msgf-rust as a new application and writes documentation
from scratch to fit it. It also takes the opportunity to clean up two
Java-historical CLI quirks: numeric-index enum flags and the singular `--mod`
flag for a file path.

## Goals

1. New `README.md` that serves both quantms pipeline operators and mass-spec
   researchers running searches directly, in a single linear narrative.
2. New single-file `DOCS.md` reference at the repo root.
3. New `CLI_MIGRATION.md` mapping Java MS-GF+ flags and legacy numeric IDs
   to the new Rust-idiomatic flag names.
4. CLI rename: replace numeric-ID enum flags with named values; rename
   `--ntt` → `--enzyme-specificity`; rename `--mod` → `--mods` with hidden
   alias.
5. Backward compatibility at runtime: the binary still accepts the legacy
   numeric forms (`--fragmentation 3`, etc.) and the old `--mod` name, so
   existing quantms scripts keep working without modification.
6. Delete the stale `docs/` user-facing tree.

## Non-goals (deferred to later iterations)

- Dockerfile rewrite (it still builds a Java JAR).
- One-time `cargo fmt` cleanup (~11k cosmetic lines).
- Thread-determinism tie-breaker fix.
- mdBook / GitHub Pages site.
- Porting Java's `ScoringParamGen` to Rust (acknowledged in `DOCS.md` as
  roadmap work; tracked as an open issue).

## Deliverables

| Path | Action | Purpose |
|---|---|---|
| `README.md` | rewrite | Linear front-door doc serving both audiences. ~190 lines. |
| `DOCS.md` | create | Single-file reference for CLI, formats, training, migration. ~505 lines. |
| `CLI_MIGRATION.md` | create | Java MS-GF+ → msgf-rust mapping + numeric-legacy → named-value table + worked examples. ~100 lines. |
| `crates/msgf-rust/src/bin/msgf-rust.rs` | edit | Add 4 `ValueEnum`-derived types, rename flags, update existing tests. |
| `crates/msgf-rust/tests/cli_smoke.rs` | edit | Add one new test: legacy numeric form and new named form produce identical output. |
| `docs/` user-facing tree | delete | All files listed in "docs/ deletion list" below. |
| `docs/superpowers/specs/` | excluded from deletion | Engineering-planning artifacts; not user-facing. |

## README.md content + structure

Linear flow, top-to-bottom. Order chosen so a researcher sees the "why
switch?" benchmark proof early, and an operator can jump straight to
Quick Start and recipes.

| # | Section | Content |
|---|---|---|
| 1 | Title + tagline + badges | CI, release, license, citation. ~8 lines. |
| 2 | What is this? | One paragraph: Rust port of MS-GF+, mzML/MGF + FASTA in, Percolator-ready `.pin` out. Names UCSD original team. ~10 lines. |
| 3 | Why msgf-rust? | Benchmark table: Rust vs Java MS-GF+ at 1% FDR on Astral / PXD001819 / TMT, plus wall-clock comparison. ~25 lines. |
| 4 | Install | Three options: (a) download a platform archive from GitHub Releases, (b) `cargo install --git`, (c) build from source. ~25 lines. |
| 5 | Quick Start | Minimal command: `msgf-rust --spectrum bsa.mgf --database bsa.fasta --output-pin out.pin`. Brief explanation of the `.pin` row. ~20 lines. |
| 6 | Common workflows | Four recipes: (a) Trypsin DDA + Percolator, (b) TMT search with mods, (c) Direct TSV output, (d) quantms pipeline integration. ~35 lines. |
| 7 | CLI summary | Table of ~15 most-used flags with one-line descriptions; link to `DOCS.md` for full reference. ~25 lines. |
| 8 | Auto-detection | Short paragraph: activation method auto-detected from mzML; param file auto-selected from (fragmentation, instrument, protocol). ~10 lines. |
| 9 | Parity vs Java MS-GF+ | One paragraph summary of what's bit-exact, what differs; link to `DOCS.md` known-divergences section. ~12 lines. |
| 10 | Citation | Cite Kim & Pevzner MS-GF+ paper. ~8 lines. |
| 11 | License | UCSD-Noncommercial; see `LICENSE`, `NOTICE`. ~6 lines. |
| 12 | Acknowledgments | UCSD original team, bigbio maintainers, quantms team. ~6 lines. |

**Total:** ~190 lines.

**Not in README** (lives in `DOCS.md` only): full CLI flag reference,
mods.txt grammar, PIN column-by-column reference, building from source in
detail, training notes, Java → Rust migration table, known-divergences
detail.

## DOCS.md content + structure

Single file, top-to-bottom. Each section is its own anchor for
deep-linking.

| # | Section | Content | ~lines |
|---|---|---|---|
| 0 | Table of contents | Anchor links to each section below. | 15 |
| 1 | CLI reference | Every flag, with description / default / value format, grouped by: required, search params, modifications, scoring, runtime, output. | 130 |
| 2 | Mods.txt format | Grammar, per-field rules, location vocabulary, `NumMods=N` directive, comment handling, 3 worked examples (cam-C + ox-M; TMT 10-plex; phospho-STY). | 50 |
| 3 | Output formats | 3a. PIN columns table. 3b. TSV columns table. 3c. Choosing between them. | 90 |
| 4 | Auto-detection | Activation-method detection from mzML CV params; param-file resolution table showing `(fragmentation, instrument, protocol) → bundled file`; instrument-class detection. | 35 |
| 5 | Building from source | Requirements (Rust 1.85+), `cargo build --release`, `cargo test --release` with notes on the 7 known-skipped tests + reasons, where the binary lands. | 30 |
| 6 | Training new `.param` files | The Rust port reuses Java MS-GF+'s `.param` files as-is. ScoringParamGen is not yet ported; tracked as roadmap work. Two paths for now: use bundled `.param` files, or train on `java-legacy` branch and point Rust at the output with `--param-file`. | 25 |
| 7 | Isobaric labeling | TMT and iTRAQ workflows: `--protocol` value, `--mods` entries, which bundled `.param` file gets auto-selected. | 35 |
| 8 | Java MS-GF+ → msgf-rust migration | 8a. Flag rename table (Java `-s` → Rust `--spectrum`, etc.). 8b. Numeric-legacy values (still accepted: `--fragmentation 3` works alongside `--fragmentation HCD`). 8c. Behavior differences (no mzXML, no mzIdentML, etc.). 8d. Known parity divergences. | 80 |
| 9 | License + citation | Full LICENSE excerpt + how to cite. | 15 |

**Total:** ~505 lines.

## CLI rename details

### Flag rename table

| Old (Java-style, current) | New (Rust-idiomatic) | Default | Accepted legacy form |
|---|---|---|---|
| `--fragmentation <0..=4>` | `--fragmentation <auto\|CID\|ETD\|HCD\|UVPD>` | `auto` | numeric 0..=4 |
| `--instrument <0..=3>` | `--instrument <low-res\|high-res\|TOF\|QExactive>` | `low-res` | numeric 0..=3 |
| `--protocol <0..=5>` | `--protocol <auto\|phospho\|iTRAQ\|iTRAQ-phospho\|TMT\|standard>` | `auto` | numeric 0..=5 |
| `--ntt <0\|1\|2>` | `--enzyme-specificity <non-specific\|semi\|fully>` | `fully` | numeric 0..=2 AND `--ntt` alias |
| `--mod <file>` | `--mods <file>` | (none) | `--mod` alias (hidden) |

Named-value conventions:
- Acronyms uppercase (community standard): HCD, CID, ETD, UVPD, TMT, iTRAQ, TOF.
- Brand names preserve common-form casing: QExactive.
- Descriptive values lowercase kebab-case: `auto`, `low-res`, `high-res`,
  `phospho`, `standard`, `non-specific`, `semi`, `fully`.
- clap parsing is case-insensitive — `--fragmentation hcd` works the same
  as `--fragmentation HCD`.

### Implementation per enum flag

```rust
#[derive(Clone, Copy, Debug, ValueEnum)]
enum Fragmentation {
    #[clap(name = "auto")] Auto,
    #[clap(name = "CID")]  Cid,
    #[clap(name = "ETD")]  Etd,
    #[clap(name = "HCD")]  Hcd,
    #[clap(name = "UVPD")] Uvpd,
}

#[arg(long, default_value = "auto", value_parser = parse_fragmentation)]
fragmentation: Fragmentation,

fn parse_fragmentation(s: &str) -> Result<Fragmentation, String> {
    // Canonical named value first (case-insensitive).
    if let Ok(v) = <Fragmentation as ValueEnum>::from_str(s, true) {
        return Ok(v);
    }
    // Legacy numeric ID (Java MS-GF+ compat).
    match s.parse::<u8>() {
        Ok(0) => Ok(Fragmentation::Auto),
        Ok(1) => Ok(Fragmentation::Cid),
        Ok(2) => Ok(Fragmentation::Etd),
        Ok(3) => Ok(Fragmentation::Hcd),
        Ok(4) => Ok(Fragmentation::Uvpd),
        _ => Err(format!(
            "invalid fragmentation `{s}`: expected auto|CID|ETD|HCD|UVPD \
             (or legacy 0..=4)"
        )),
    }
}
```

Same shape for `Instrument`, `Protocol`, `EnzymeSpecificity`.

### `--mods` rename

```rust
#[arg(long = "mods", alias = "mod", value_name = "MODFILE")]
mods: Option<PathBuf>,
```

`alias` (not `visible_alias`) means `--mod` is still accepted but `--help`
only shows `--mods`.

### Quantms compat policy

For v0.1.0 (the cutover release) the numeric form is "Java legacy" rather
than "deprecated Rust v0". Accept silently — no deprecation warning to
stderr. Migration is documented in `DOCS.md` §8 and `CLI_MIGRATION.md`.
Working quantms scripts keep working with zero changes.

### Internal code changes

- Replace `Option<u8>` enum fields + numeric-positional calls
  (`resolve_bundled_param(Some(3), Some(3), Some(4))`) with strongly-typed
  enums (`resolve_bundled_param(Fragmentation::Hcd, Instrument::QExactive,
  Protocol::Tmt)`).
- Update the 15 `param_resolver_tests` (~30 line diff).
- The auto-detect path (`resolve_bundled_param_for_activation`) now
  constructs the enum variants directly instead of numeric IDs.

## CLI_MIGRATION.md content

~100 lines. Two tables + worked examples.

- **Table A — Java MS-GF+ flag → msgf-rust flag.** Full mapping: `-s` →
  `--spectrum`, `-d` → `--database`, `-o` → `--output-pin`, `-mod` →
  `--mods`, `-tda 1` → "not needed, decoys auto-generated", `-inst N` →
  `--instrument <named-value>`, etc.
- **Table B — Numeric legacy → named values.** The same content as the
  Implementation table above, formatted for users porting scripts.
- **3 worked examples.** A Java MS-GF+ command line rewritten as a
  msgf-rust command line, side-by-side, for: (a) plain Trypsin DDA + 20ppm,
  (b) TMT 10-plex search, (c) phospho-STY search.

## docs/ deletion list

Delete (all in this PR):

- `docs/msgfplus.md`
- `docs/msgfdb_modfile.md`
- `docs/buildsa.md`
- `docs/output.md`
- `docs/readme.md`
- `docs/troubleshooting.md`
- `docs/training-scoring-models.md`
- `docs/isobariclabeling.md`
- `docs/changelog.md`
- `docs/parameterfiles/` (15 `.txt` files)
- `docs/examples/` (`Mods.txt`, `enzymes.txt`, etc. — content migrates into
  `DOCS.md` as inline examples)
- `docs/benchmarks/` (3 PNG figures from the Java perf comparison; stale)

Keep (excluded from deletion):

- `docs/superpowers/specs/` — engineering-planning subdirectory, not
  user-facing docs. This document lives here.

Already gitignored, no action:

- `docs/parity-analysis/` — local-only iter notes from iter1-38 development.

## Testing

| File | Change |
|---|---|
| `crates/msgf-rust/src/bin/msgf-rust.rs` (`param_resolver_tests`, 15 tests) | Update each from `resolve_bundled_param(Some(3), Some(3), Some(4))` → `resolve_bundled_param(Fragmentation::Hcd, Instrument::QExactive, Protocol::Tmt)`. Mechanical. |
| `crates/msgf-rust/tests/cli_smoke.rs` (7 existing integration tests) | The tests use `--fragmentation 3 --instrument 3 --protocol 4` strings; these still work (legacy accepted), so no behavior change is required. |
| `crates/msgf-rust/tests/cli_smoke.rs` (new test) | `cli_accepts_both_named_and_numeric_param_values`: run a search with `--fragmentation 3 --protocol 4` (legacy) and again with `--fragmentation HCD --protocol TMT` (canonical); assert PIN outputs are byte-identical. Guards the back-compat path. |

CI workflow (`.github/workflows/ci.yml`) — no change. The 7 currently-skipped
tests remain skipped for the reasons documented inline.

## Commit plan

One PR (`iter39-docs-rewrite` → `dev`), five reviewable commits in order:

1. `feat(cli): rename param flags to Rust-idiomatic named values with legacy compat` — CLI rename, enum types, custom parsers, updated `param_resolver_tests`, new round-trip test.
2. `docs: write new README.md (post-cutover, dual audience, linear narrative)` — replace `README.md`.
3. `docs: add DOCS.md (single-file reference)` — new `DOCS.md`.
4. `docs: add CLI_MIGRATION.md (Java → Rust + numeric legacy mapping)` — new file.
5. `docs: delete docs/ tree (content migrated to DOCS.md)` — `git rm -r` everything from the deletion list above; `docs/superpowers/` is preserved.

PR title: `iter39: docs + CLI rename for the post-cutover state`

## Risks

- **Risk:** A quantms script uses `--fragmentation 3` and we silently break it. **Mitigation:** the new round-trip integration test in `cli_smoke.rs` ensures legacy numeric values resolve to the same enum variants as the named values, locked in CI.
- **Risk:** Hidden `--mod` alias is missed by a user trying to migrate. **Mitigation:** `CLI_MIGRATION.md` calls it out as a top-line "what's renamed" entry.
- **Risk:** The deletion of `docs/parameterfiles/*.txt` breaks external links from third-party tooling that bundled those templates. **Mitigation:** Low — these were Java `-conf` templates; no equivalent Rust mechanism exists. `CLI_MIGRATION.md` covers the closest Rust path (direct CLI flags + `--param-file`).
- **Risk:** README + DOCS.md diverge from the binary over time. **Mitigation:** acceptable — both files are short enough that any future iteration that touches CLI flags or output format updates them in the same PR.

## Out of scope (re-affirming)

- Dockerfile rewrite
- One-time `cargo fmt`
- Thread-determinism tie-breaker
- mdBook / Pages site
- Porting ScoringParamGen
