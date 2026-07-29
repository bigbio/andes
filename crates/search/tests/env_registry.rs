//! Every `ANDES_*` environment variable must appear in `docs/ENV_VARS.md`.
//!
//! The engine reads 40+ environment variables. Some are diagnostics, some are
//! experiments — and some are *corrections that were written, validated, and then left
//! disabled by default*. That last category is the problem this guards.
//!
//! Two concrete cases motivated it. The Y-ladder is an unnormalised sum whose
//! expectation grows with glycan size; the correction for that bias exists, is gated
//! behind `ANDES_GLYCO_LADDER_NORM`, ships off, and was measured afterwards at +19
//! identifications on the benchmark — meaning the selector weights were re-tuned around
//! a bias whose fix was sitting in the tree. Separately, a deconvolution-aware c/z charge
//! guard was written and never wired to a call site at all.
//!
//! An undocumented switch is invisible twice over: users cannot find it, and we forget
//! it exists and tune around it. Requiring a registry entry makes the cost of hiding a
//! correction higher than the cost of deciding what its default should be.
//!
//! The registry also records the gating form. `var_os` means presence-only, so
//! `VAR=0` *enables* the behaviour — the opposite of what anyone typing that intends.
//! That bug has been fixed twice in this codebase already.

use std::collections::BTreeSet;

/// Scan the workspace sources for `"ANDES_*"` string literals.
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
            // Skip build artefacts; they contain generated copies of the sources.
            if p.file_name().is_some_and(|n| n == "target") {
                continue;
            }
            collect(&p, out);
        } else if p.extension().is_some_and(|x| x == "rs") {
            let Ok(text) = std::fs::read_to_string(&p) else {
                continue;
            };
            // This file names variables only inside its own documentation, and the
            // registry test itself must not be a source of truth.
            if p.file_name().is_some_and(|n| n == "env_registry.rs") {
                continue;
            }
            for (i, _) in text.match_indices("\"ANDES_") {
                let rest = &text[i + 1..];
                if let Some(end) = rest.find('"') {
                    let name = &rest[..end];
                    if name
                        .chars()
                        .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_')
                        && name.len() > "ANDES_".len()
                    {
                        out.insert(name.to_string());
                    }
                }
            }
        }
    }
}

#[test]
fn every_env_var_is_documented() {
    let registry = include_str!("../../../docs/ENV_VARS.md");
    let found = vars_in_sources();
    assert!(
        found.len() > 20,
        "scan found only {} variables — the scanner is broken, not the code",
        found.len()
    );

    let undocumented: Vec<&String> = found
        .iter()
        .filter(|v| !registry.contains(&format!("`{v}`")))
        .collect();

    assert!(
        undocumented.is_empty(),
        "these environment variables are read by the engine but absent from \
         docs/ENV_VARS.md: {undocumented:?}\n\
         Add a row with its gating form, location and purpose. If the variable gates a \
         CORRECTION that is off by default, say so explicitly and say why it has not \
         been made the default — that decision is the whole point of the registry."
    );
}

#[test]
fn registry_has_no_stale_entries() {
    // A registry that outlives its variables is as misleading as a missing one: it
    // suggests switches that no longer do anything.
    let registry = include_str!("../../../docs/ENV_VARS.md");
    let found = vars_in_sources();
    let mut stale = Vec::new();
    for line in registry.lines() {
        // Table rows only; prose in the header also mentions variable names.
        if !line.starts_with("| `ANDES_") {
            continue;
        }
        if let Some(name) = line
            .split('`')
            .nth(1)
            .filter(|n| n.starts_with("ANDES_"))
        {
            if !found.contains(name) {
                stale.push(name.to_string());
            }
        }
    }
    assert!(
        stale.is_empty(),
        "docs/ENV_VARS.md lists variables the code no longer reads: {stale:?} — \
         remove them so the registry stays trustworthy"
    );
}
