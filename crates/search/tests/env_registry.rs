use std::collections::BTreeSet;

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
fn every_env_var_is_documented() {
    let registry = include_str!("../../../docs/ENV_VARS.md");
    let found = vars_in_sources();
    assert!(
        !found.is_empty(),
        "scan found no variables at all — the scanner is broken, not the code"
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
    let registry = include_str!("../../../docs/ENV_VARS.md");
    let found = vars_in_sources();
    let mut stale = Vec::new();
    for line in registry.lines() {
        if !line.starts_with("| `ANDES_") {
            continue;
        }
        if let Some(name) = line.split('`').nth(1).filter(|n| n.starts_with("ANDES_")) {
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
