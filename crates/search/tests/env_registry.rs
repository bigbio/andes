//! Every `ANDES_*` environment variable the code reads must be a test-harness variable
//! listed in [`TEST_HARNESS_VARS`]; the policy is in the `cli-flags-not-env-vars` skill
//! (`.claude/skills/cli-flags-not-env-vars/SKILL.md`).

use std::collections::BTreeSet;

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

fn vars_in_sources() -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .to_path_buf();
    collect(&root, &mut out);
    out
}

fn collect(dir: &std::path::Path, out: &mut BTreeSet<String>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for e in entries.flatten() {
        let p = e.path();
        if p.is_dir() {
            if p.file_name().is_some_and(|n| n == "target") {
                continue;
            }
            collect(&p, out);
        } else if p.extension().is_some_and(|x| x == "rs") {
            let Ok(text) = std::fs::read_to_string(&p) else {
                continue;
            };
            if p.file_name().is_some_and(|n| n == "env_registry.rs") {
                continue;
            }
            for pat in ["env::var(", "env::var_os("] {
                for (i, _) in text.match_indices(pat) {
                    let rest = &text[i + pat.len()..];
                    let Some(open) = rest.find('"') else { continue };
                    if rest[..open].trim() != "" {
                        continue;
                    }
                    let after = &rest[open + 1..];
                    if let Some(end) = after.find('"') {
                        let name = &after[..end];
                        if name.to_ascii_uppercase().starts_with("ANDES_") {
                            out.insert(name.to_string());
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn every_env_var_is_a_listed_test_harness_variable() {
    let found = vars_in_sources();
    assert!(
        !found.is_empty(),
        "scan found no variables at all — the scanner is broken, not the code"
    );
    let listed: BTreeSet<&str> = TEST_HARNESS_VARS.iter().map(|(v, _)| *v).collect();
    let unlisted: Vec<&String> = found
        .iter()
        .filter(|v| !listed.contains(v.as_str()))
        .collect();
    assert!(
        unlisted.is_empty(),
        "these environment variables are read by the code but are not test-harness \
         variables in TEST_HARNESS_VARS: {unlisted:?}\n\
         A setting that changes a result must be a typed CLI flag (see the \
         cli-flags-not-env-vars skill). Only a variable that selects a test fixture may be \
         added to the list."
    );
}

#[test]
fn test_harness_list_has_no_stale_entries() {
    let found = vars_in_sources();
    let stale: Vec<&str> = TEST_HARNESS_VARS
        .iter()
        .map(|(v, _)| *v)
        .filter(|v| !found.contains(*v))
        .collect();
    assert!(
        stale.is_empty(),
        "TEST_HARNESS_VARS lists variables the code no longer reads: {stale:?} — remove them"
    );
}
