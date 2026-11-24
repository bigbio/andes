# Structure checks

`crates/andes/tests/structure.rs` turns the contributor rules into tests. Run them with
`cargo test -p andes --release --test structure`. Helpers are in
`crates/andes/tests/support/`, baselines in `crates/andes/tests/ratchets/`.

Never raise a baseline, and never regenerate one, to make a test pass. A baseline only
records violations that already existed; fix the code instead, or discuss the exception in
the pull request.

| Test | Rule | Baseline | How to fix |
|---|---|---|---|
| `every_skill_points_to_its_docs_dev_page` | Each `.claude/skills/<name>/SKILL.md` names `docs/dev/<name>.md`, that page exists, and the skill has at most 12 lines. | none | Point the skill at its page; keep the rules in the page. |
| `agents_md_lists_every_docs_dev_page` | `AGENTS.md` mentions every page under `docs/dev/`. | none | Add a row to the table in `AGENTS.md`. |
| `contributing_points_to_the_rules` | `CONTRIBUTING.md` links `docs/dev/` and `AGENTS.md`. | none | Restore the links. |
| `no_new_source_file_over_the_budget` | A file under `crates/*/src` that is not on `file_lines.tsv` has at most 1,500 lines. The files on the list (over the budget when the checks were added) are not gated: CI prints them, with how far each moved, from the ignored test `report_listed_files_over_the_budget`. | `file_lines.tsv` | Split the file into modules. When a listed file is split under the budget, remove its line (the regeneration command does it); from then on it is gated. |
| `no_new_function_over_the_budget` | A function in a file under `crates/*/src` has at most 300 lines (from its `fn` line to its closing brace, test code excluded), unless it is on `fn_lines.tsv` as `path::name`. The listed functions are not gated: CI prints them from the ignored test `report_listed_functions_over_the_budget`. | `fn_lines.tsv` | Split the function into steps (a function per stage, its inputs in a struct when there are many). When a listed function is split under the budget, remove its line (the regeneration command does it); from then on it is gated. |
| `no_new_unwrap_or_expect_in_live_code` | No new `.unwrap()` or `.expect(` in live code (not tests, comments or strings). `unwrap_err` and `expect_err` are not counted. A documented invariant is not exempt yet. | `unwrap_expect.tsv` | Return an error, or restructure with `let ... else` or `match`. |
| `no_print_macros_in_live_code` | No `println!`, `eprintln!`, `print!`, `eprint!` or `dbg!` in live code, except in the files on `PRINT_ALLOWED` in `structure.rs`. | none (allowlist in code) | Log through `log` (see `clean-logs`). A file whose output is a command's result goes on `PRINT_ALLOWED` with the reason. |
| `every_new_module_has_a_doc` | A source file starts with a `//!` doc; inner attributes may come first. | `module_docs.tsv` | Add a `//!` doc of 2-4 lines. |
| `every_cited_key_is_in_the_registry` | Every citation marker (`[cite:` then a registry key then `]`) in source, `docs/`, `README.md`, `DOCS.md`, `TRAIN.md` and `ARCHITECTURE.md` names a key in `docs/citations.yaml`. A marker whose key is not `[a-z0-9_]+`, or that is not closed within 64 characters, is malformed and fails. Prose with whitespace after `[cite:` is ignored. | none | Add the registry entry first, or correct the key. |
| `registry_entries_are_complete_and_unique` | Each registry entry has a valid key, authors, title, a year in 1900-2100 and a DOI or URL. | none | Complete or de-duplicate the entry. |
| `no_site_identifiers_in_committed_text` | No cluster, host or user path (`/home/<name>/`, `/Users/<name>/`, `/Data/<name>/`, shared-storage prefixes, site names) in committed text, `structure.rs` included. Site names are matched by hash. | none | Use a placeholder or an argument. |
| `codeowners_paths_exist` | Every path in `.github/CODEOWNERS` exists. | none | Update the path. |
| `architecture_names_every_crate` | `ARCHITECTURE.md` names every directory under `crates/` as `` `crates/<name>` ``. | none | Add the crate to the document. |

Test modules kept in their own files go under a `tests/` directory next to the module they test (for example `crates/search/src/engine/tests/dedup.rs`, declared with `#[cfg(test)] #[path = "engine/tests/dedup.rs"] mod dedup_tests;`). Files under `tests/` are not live code, so the unwrap, module-doc, print-macro and size checks skip them.

## Module and file names

- **One word where possible:** `dia/llr.rs`, not `dia_llr.rs`. A family of files that share a
  prefix lives in a folder named after that prefix.
- **A module with children is `name.rs` plus a `name/` folder, never `name/mod.rs`.**
- **Name files after what they do.** Use field terms (`pin`, `mbr`, `tdc`, `fdr`, `xic`, `ipf`,
  `tmt`, `ms1`, `ms2`) rather than private abbreviations, and no format or role suffixes such as
  `_cmd`, `_out`, `_table` or `_parquet`.
- **Do not name a module after a crate its own crate depends on**, so that `glyco::` or `quant::`
  never becomes ambiguous. This is why `search` has `glycopeptides` and `crosslinks`, and `output`
  has `quantities`.
- **Test files are the exception:** they describe a behaviour (`pin_goldens.rs`) and may use
  several words.

Two tests enforce the first two rules: `no_mod_rs_files` and `no_prefix_family_beside_its_folder`.
The second rejects a file `<prefix>_<rest>.rs` that sits beside a `<prefix>.rs` or `<prefix>/` and
names the folder path it belongs in.

## Baselines

A baseline is a TSV of `path<TAB>count` lines. The tests compare counts per file:

- A count above its baseline fails.
- A file with no baseline line and a count above 0 fails. After a rename or a split, move
  the line, or carry the counts to the new files only if the total does not rise, and say so
  in the pull request.
- A count below its baseline fails until the baseline is lowered, and at 0 until the line is
  removed, so a baseline cannot go stale.
- A line naming a file that no longer exists fails until it is removed or moved.

After a cleanup, lower the counts with
`cargo test -p andes --release --test structure -- --ignored write_<rule>_baseline`, where
`<rule>` is `file_lines`, `fn_lines`, `unwrap_expect` or `module_docs`. The command refuses to raise any
count or add a line.

Because lowering is strict, two parallel pull requests that clean the same file can conflict
on its baseline line; the second one rebases and sets the count to the new value.
