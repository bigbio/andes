//! Every `ANDES_*` environment variable named in the code must be a test-harness variable
//! listed in [`TEST_HARNESS_VARS`], and only test code may name one; the policy is in the
//! `cli-flags-not-env-vars` skill (`.claude/skills/cli-flags-not-env-vars/SKILL.md`).
//!
//! The scan looks for `ANDES_*` names inside any string literal of live code (comments are
//! ignored), so a key held in a constant or read through an imported `var` is found too.

#[path = "../../andes/tests/support/source_scan.rs"]
mod source_scan;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use source_scan::{blank, is_test_path, rust_files, test_item_ranges, Keep};

/// The test-harness variables: each selects an optional fixture at `cargo test` time and is
/// never read by the shipped binary. Paired with the file that reads it.
const TEST_HARNESS_VARS: &[(&str, &str)] = &[
    ("ANDES_GBDT_DIR", "crates/scoring/tests/gbdt_bench.rs"),
    ("ANDES_GLYCO_GDB", "crates/glyco/src/glycan_first.rs"),
    (
        "ANDES_GLYCO_MGF",
        "crates/glyco/tests/glycan_first_mouse.rs",
    ),
    ("ANDES_TEST_D", "crates/input/tests/timstof_d_loads.rs"),
    (
        "ANDES_TEST_PERCOLATOR_BIN",
        "crates/output/tests/percolator_integration.rs",
    ),
    ("ANDES_TEST_RAW", "crates/input/tests/thermo_raw.rs"),
    (
        "ANDES_TRAIN_BENCH",
        "crates/model-train/tests/yield_nonregression.rs",
    ),
    (
        "ANDES_V1_STORE",
        "crates/model-train/tests/partition_parity.rs",
    ),
];

/// One place where live code names an `ANDES_*` variable.
struct Use {
    file: PathBuf,
    line: usize,
    in_test_code: bool,
}

/// `ANDES_*` names inside the string literals of `strings` (a [`Keep::Strings`] blanking),
/// with their byte offsets.
fn names_in_strings(strings: &str, code: &str) -> Vec<(String, usize)> {
    let mut out = Vec::new();
    let mut from = 0;
    while let Some(p) = strings[from..].find("ANDES_") {
        let at = from + p;
        let end = strings[at..]
            .find(|c: char| !(c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_'))
            .map_or(strings.len(), |e| at + e);
        let in_literal = code.as_bytes()[at] == b' ';
        if in_literal && end > at + "ANDES_".len() {
            out.push((strings[at..end].to_string(), at));
        }
        from = end;
    }
    out
}

fn uses(root: &Path) -> BTreeMap<String, Vec<Use>> {
    let mut out: BTreeMap<String, Vec<Use>> = BTreeMap::new();
    for file in rust_files(root) {
        if file.ends_with("search/tests/env_registry.rs") {
            continue;
        }
        let Ok(text) = std::fs::read_to_string(&file) else {
            continue;
        };
        let strings = blank(&text, Keep::Strings);
        let code = blank(&text, Keep::Code);
        let tests = test_item_ranges(&code);
        let test_path = is_test_path(file.strip_prefix(root).unwrap_or(&file));
        for (name, at) in names_in_strings(&strings, &code) {
            out.entry(name).or_default().push(Use {
                file: file.clone(),
                line: text[..at].matches('\n').count() + 1,
                in_test_code: test_path || tests.iter().any(|r| r.contains(&at)),
            });
        }
    }
    out
}

fn crates_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .to_path_buf()
}

#[test]
fn every_env_var_is_a_listed_test_harness_variable() {
    let found = uses(&crates_root());
    assert!(
        !found.is_empty(),
        "scan found no variables at all — the scanner is broken, not the code"
    );
    let unlisted: Vec<&String> = found
        .keys()
        .filter(|v| !TEST_HARNESS_VARS.iter().any(|(l, _)| l == v))
        .collect();
    assert!(
        unlisted.is_empty(),
        "these environment variables are named by the code but are not test-harness \
         variables in TEST_HARNESS_VARS: {unlisted:?}\n\
         A setting that changes a result must be a typed CLI flag (see the \
         cli-flags-not-env-vars skill). Only a variable that selects a test fixture may be \
         added to the list."
    );
}

#[test]
fn only_test_code_names_an_andes_variable() {
    let shipped: Vec<String> = uses(&crates_root())
        .into_iter()
        .flat_map(|(name, uses)| {
            uses.into_iter()
                .filter(|u| !u.in_test_code)
                .map(move |u| format!("{name} at {}:{}", u.file.display(), u.line))
        })
        .collect();
    assert!(
        shipped.is_empty(),
        "shipped code names ANDES_* environment variables; make them CLI flags:\n{}",
        shipped.join("\n")
    );
}

#[test]
fn every_listed_variable_is_still_read_where_registered() {
    let root = crates_root();
    let found = uses(&root);
    let mut stale = Vec::new();
    for (name, path) in TEST_HARNESS_VARS {
        let registered = root.join(path.trim_start_matches("crates/"));
        let read_there = found
            .get(*name)
            .is_some_and(|us| us.iter().any(|u| u.file == registered));
        if !read_there {
            stale.push(format!("{name} (registered at {path})"));
        }
    }
    assert!(
        stale.is_empty(),
        "TEST_HARNESS_VARS lists variables the live code no longer names where registered: \
         {stale:?} — fix the path or remove the entry"
    );
}

#[test]
fn a_commented_out_read_or_a_constant_key_is_seen_correctly() {
    let src = "// std::env::var(\"ANDES_OLD\")\nconst KEY: &str = \"ANDES_NEW\";\n";
    let names = names_in_strings(&blank(src, Keep::Strings), &blank(src, Keep::Code));
    let names: Vec<&str> = names.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names, ["ANDES_NEW"]);
}
