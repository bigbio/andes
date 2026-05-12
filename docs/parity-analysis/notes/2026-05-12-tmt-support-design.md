# TMT support in msgf-rust — design note

**Date:** 2026-05-12
**Branch:** `feat/tmt-modified-support`
**Worktree:** `astral-speed-tmt`

## Motivation

On PXD007683 file `a05058.mzML` (TMT 6-plex CID dataset), Java MS-GF+ run
with `-m 1 -inst 1 -e 1 -protocol 4 -mod mods.txt` finds **28,790 targets**
in 3m07s; `msgf-rust` finds **19,937** in 6m01s — slower **and** ~31% fewer
PSMs. Java's flags decode to:

| Java flag      | Value | Meaning                              |
| -------------- | ----- | ------------------------------------ |
| `-m 1`         | 1     | Fragmentation: CID                   |
| `-inst 1`      | 1     | Instrument: HighRes                  |
| `-e 1`         | 1     | Enzyme: Trypsin                      |
| `-protocol 4`  | 4     | Protocol: TMT                        |
| `-mod mods.txt`| -     | Fixed TMT6plex on K and peptide N-term; CAM on C; Oxidation on M |

The current Rust binary cannot express any of these except enzyme
(implicit). With no `--mod`, the candidate-peptide masses for TMT-labelled
peptides are off by `+229.16` per K and per N-term, so the precursor
window never matches the real PSMs. That's the entire ~31% miss.

The slowdown is a secondary problem (likely related to scoring with the
wrong .param table) and is not addressed here.

## What already exists (the good news)

The underlying machinery for `--mod` + per-instrument scoring is **almost
entirely present** in `msgf-rust`. Specifically:

### 1. mods.txt parser — already implemented
`rust/crates/model/src/modification.rs::Modification::from_mods_txt_line`
parses Java's 5-field comma-separated format (`mass,residue,fix|opt,location,name`).
Tests cover: Carbamidomethyl, Oxidation, wildcard N-term, ProtNTerm Acetyl,
negative mass deltas, case-insensitive location names.

### 2. Mods-file loader — already implemented
`rust/crates/model/src/aa_set.rs::AminoAcidSetBuilder::add_mods_from_file`
reads a Mods.txt file, skips `#`-comments + blank lines, dispatches each
mod to `fixed_mods` or `variable_mods`, and surfaces line-numbered parse
errors. Unit-tested (`add_mods_from_file_parses_real_format`,
`add_mods_from_file_reports_line_number`).

### 3. Variable-mod expansion — already implemented
`rust/crates/search/src/candidate_gen.rs::expand_mod_combinations`
expands every span into all combinations of variable-mod applications,
consulting Anywhere + NTerm/ProtNTerm at position 0 and CTerm/ProtCTerm
at position n−1. The `max_variable_mods_per_peptide` cap is enforced
inside the recursion. Fixed mods are applied via `aa_set.variants_for`,
which has them folded into the residue's `Anywhere` variant.

### 4. .param files for every (frag, inst, protocol) combo — already bundled
`src/main/resources/ionstat/` contains ~50 binary `.param` files
including the four we need for TMT:
- `CID_LowRes_Tryp.param` and `HCD_QExactive_Tryp_TMT.param`
- `HCD_HighRes_Tryp_TMT.param` (likely match for `-m 1 -inst 1` though
  named HCD — Java's `-m 1` is CID and `-inst 1` is HighRes; see "Open
  question on .param selection" below)
- `CID_HighRes_Tryp.param` exists too.

The Rust `Param::load_from_file` already parses any of these. The
`.param` file's embedded `SpecDataType { activation, instrument, enzyme,
protocol }` is the source of truth for which scoring-table partition is
loaded — `msgf-rust`'s search/scoring code does **not** branch on
protocol/instrument enums at runtime. So "switching to TMT scoring
tables" = "loading the right .param file." No new scoring code needed.

### 5. Fixed-mod handling at peptide-mass level — already works
The candidate's peptide mass is computed from per-residue masses inside
`AminoAcidSetBuilder::build` (the `with_mod` path adds `mass_delta` to
the residue's mass). So loading a fixed `K + 229.163` mod via
`add_fixed_mod` correctly bumps every K in every candidate's mass.

## What is missing

**Sub-problem A — pure CLI plumbing:** the binary doesn't accept
`--mod`, `--protocol`, `--instrument`, `--fragmentation`. None of these
need new search/scoring code; they just need:
1. Clap flags on `Cli`.
2. Route `--mod` through `AminoAcidSetBuilder::add_mods_from_file`.
3. Translate `(--fragmentation, --instrument, --protocol)` into a
   bundled `.param` filename if `--param-file` was not given.

**Sub-problem B — Java-mod-file feature parity:** the existing Rust
parser is a strict subset of Java's. The gaps observed against real
`mods.txt` files in `src/test/resources/`:

| Feature                       | Java | Rust | Needed for TMT? |
| ----------------------------- | ---- | ---- | ----------------|
| `mass` numeric (e.g. `229.16`) | Yes  | Yes  | Yes (works)     |
| `composition` (e.g. `C2H3NO`)  | Yes  | No   | No — TMT mods.txt files in the wild use numeric mass; we'll require that |
| `NumMods=N` header line       | Yes  | No   | Helpful (sets max variable mods per peptide); easy fix |
| Multi-residue spec (`STY`)    | Yes  | No   | No for TMT (single residue K + wildcard \*); yes for phospho |
| `custom` AA entries           | Yes  | No   | No                                    |
| Comment-stripping mid-line (`mass,...,name # comment`) | Yes | No (only full-line `#`) | Nice-to-have |

For **this session**, only `NumMods=N` is needed to be a TMT-compatible
mods-file consumer for the typical case. Composition strings, multi-
residue specs, and custom AAs can be deferred.

## What we'll ship this session

### Scope (in this session)
1. **Add CLI flags** to `msgf-rust.rs`:
   - `--mod <PATH>` (alias for the existing builder method)
   - `--protocol <0..5>` (0=Automatic, 1=Phosphorylation, 2=iTRAQ, 3=iTRAQPhospho, 4=TMT, 5=Standard) — matches Java's index order
   - `--instrument <0..3>` (0=LowRes, 1=HighRes, 2=TOF, 3=QExactive) — matches Java
   - `--fragmentation <0..4>` (0=Auto/CID, 1=CID, 2=ETD, 3=HCD, 4=UVPD) — matches Java's `-m`
2. **Bundled-param auto-resolve.** When `--param-file` is not given,
   compose a filename from `{fragmentation}_{instrument}_Tryp[_{protocol}].param`
   (with `_TMT`, `_iTRAQ`, `_iTRAQPhospho`, `_Phosphorylation` suffixes
   for the four non-standard protocols), and fall back to the existing
   `HCD_QExactive_Tryp.param` default. Fail loudly if the selected
   filename doesn't exist on disk.
3. **Extend `add_mods_from_file`** to ignore `NumMods=N` lines and
   capture the value into a returned struct; the caller wires it into
   `params.max_variable_mods_per_peptide`. (Just a few lines in
   `aa_set.rs`.)
4. **Unit tests:**
   - `aa_set::add_mods_from_file` parses a TMT-style mods.txt (TMT6plex
     fixed on K, fixed on \*/N-term, CAM on C, Oxidation on M, `NumMods=3`)
     and produces the right fixed/variable lists.
   - `param_resolver` test — given (frag=CID, inst=HighRes, protocol=TMT),
     returns `CID_HighRes_Tryp_TMT.param` (and asserts the file exists in
     the source tree).
5. **CLI smoke test** — invoke binary with `--fragmentation 3 --instrument 3 --protocol 4 --mod <fixture>` on the existing BSA mgf and check it exits 0. (Doesn't test PSM recall; just that the new flags parse and the param resolver finds a real file.)
6. **PXD001819 no-regression** — verify with the existing
   `gf_java_parity` integration test that default flags produce
   unchanged output.

### Out of scope (deferred)
- Composition strings (`C2H3N1O1`) — can be added when we hit a
  real mods.txt that uses them. Most TMT/phospho parameter sheets in
  the wild quote numeric masses; mods.txt files in the test fixtures
  also use both forms and will need conversion to numeric.
  **NOTE:** the existing `src/test/resources/benchmark/PXD001819/mods.txt`
  uses composition strings, so the no-regression run will not load
  that file — it loads zero mods by default in the smoke path.
- Multi-residue specs (`STY`) — needed for phospho but not TMT.
- `custom` AA entries — niche.
- The performance gap (~2× wall time on PXD007683) — separate task.
- A new `--mods-txt-format` accepting Java's composition syntax — see
  follow-up note.

### Rough LOC estimate
- `aa_set.rs` `NumMods=` handling + return type change: **~30 LOC**
- `msgf-rust.rs` new flags + param resolver + builder wiring: **~120 LOC**
- New unit tests (TMT mods.txt + param resolver): **~80 LOC**
- New cli_smoke test for flag parsing: **~30 LOC**
- Design note (this file): **~150 LOC**

**Total:** ~410 LOC, well under the 500-LOC decision-gate threshold.

## Open questions / follow-ups

### Q1: Composition string support for real-world mods.txt
Both the project's own `PXD001819/mods.txt` and `TestCandidatePeptideGrid.txt`
use composition strings (`C2H3N1O1`, `HO3P`) instead of numeric masses.
Without composition-string support, those files cannot be parsed by
`msgf-rust`, even though `add_mods_from_file` is wired through. Java's
`Composition::getMass` handles C/H/N/O/S/P + Br/Cl/Fe/Se. Porting that
to Rust is straightforward (~60 LOC + element-mass table) and is the
single biggest follow-up.

### Q2: The exact .param to use for `-m 1 -inst 1 -e 1 -protocol 4`
The user's Java invocation is `-m 1 -inst 1 -e 1 -protocol 4`, which
decodes to (CID, HighRes, Trypsin, TMT). The bundled file directory has
`HCD_HighRes_Tryp_TMT.param` and `HCD_QExactive_Tryp_TMT.param` but no
`CID_HighRes_Tryp_TMT.param`. Java's `NewScorerFactory.get(...)` falls
back through the (activation × instrument × enzyme × protocol) table to
the most-similar entry that exists. For this session we hardwire the
fallback as "if the exact `.param` doesn't exist, error out with a
helpful message" — port of the Java fallback table is a follow-up.

### Q3: Whether `--fragmentation`/`--instrument`/`--protocol` should
also influence anything else besides .param selection
In `msgf-rust`'s current scoring path, the answer is "no": the loaded
`Param` is the single source of truth, and none of the search/scoring
code branches on the enums. If we ever port Java's runtime protocol-
inference (`AUTOMATIC` → looks at the mods present and picks one), this
would need re-thinking; for now the user must specify the protocol
explicitly (matching Java's `-protocol N` flag).

### Q4: Performance gap (separate problem)
Even after `--mod` is wired up, the 2× wall-time gap on PXD007683
(Java 3:07 vs Rust 6:01 with no mods) is unexplained. Hypothesis:
either (a) the fragment-tolerance default of 0.5 Da is wrong for HighRes
data and the wider tolerance is doing extra scoring work, or (b)
something in the spectra-loading path is slow on this particular file.
Out of scope for this session.

## Acceptance criteria for this session

1. `cargo build --workspace --release` succeeds.
2. `cargo test --workspace --lib` — all existing tests pass plus the new
   ones below.
3. `cargo test -p search --test gf_java_parity` — unchanged behavior.
4. New unit test: `aa_set::tmt_style_mods_file_parses` (TMT6plex K + N-term + CAM C + Ox M + NumMods=3).
5. New CLI smoke test asserts default flags still produce a non-empty PIN
   on the BSA fixture, and `--fragmentation 3 --instrument 3 --protocol 4`
   resolves to a real .param file.
6. Commits land on `feat/tmt-modified-support` (already checked out in
   this worktree).
