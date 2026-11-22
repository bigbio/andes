//! Structure checks on the contributor rules: every skill points to its `docs/dev` page,
//! `AGENTS.md` lists every page, and `CONTRIBUTING.md` points to the rules.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

#[path = "support/source_scan.rs"]
mod source_scan;

#[path = "support/ratchet.rs"]
mod ratchet;

use ratchet::compare;
use source_scan::{line_count, live_code, rel_key};

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

#[test]
fn every_skill_points_to_its_docs_dev_page() {
    let root = workspace_root();
    let skills = root.join(".claude").join("skills");
    let mut checked = 0;
    for entry in fs::read_dir(&skills).expect("read .claude/skills") {
        let dir = entry.expect("skill entry").path();
        if !dir.is_dir() {
            continue;
        }
        let name = dir.file_name().unwrap().to_string_lossy().into_owned();
        let skill = read(&dir.join("SKILL.md"));
        let page = format!("docs/dev/{name}.md");
        assert!(
            skill.contains(&page),
            ".claude/skills/{name}/SKILL.md must contain `{page}`"
        );
        assert!(root.join(&page).is_file(), "{page} does not exist");
        let lines = skill.lines().count();
        assert!(
            lines <= 12,
            ".claude/skills/{name}/SKILL.md has {lines} lines; the rules live in {page}, keep the skill to at most 12"
        );
        checked += 1;
    }
    assert!(checked > 0, "no skills found under .claude/skills");
}

#[test]
fn agents_md_lists_every_docs_dev_page() {
    let root = workspace_root();
    let agents = read(&root.join("AGENTS.md"));
    let mut checked = 0;
    for entry in fs::read_dir(root.join("docs").join("dev")).expect("read docs/dev") {
        let path = entry.expect("docs/dev entry").path();
        if path.extension().and_then(|e| e.to_str()) != Some("md") {
            continue;
        }
        let file = path.file_name().unwrap().to_string_lossy().into_owned();
        assert!(
            agents.contains(&file),
            "AGENTS.md does not list docs/dev/{file}"
        );
        checked += 1;
    }
    assert!(checked > 0, "no pages found under docs/dev");
}

#[test]
fn contributing_points_to_the_rules() {
    let contributing = read(&workspace_root().join("CONTRIBUTING.md"));
    assert!(
        contributing.contains("docs/dev/"),
        "CONTRIBUTING.md must link docs/dev/"
    );
    assert!(
        contributing.contains("AGENTS.md"),
        "CONTRIBUTING.md must link AGENTS.md"
    );
}

fn map(entries: &[(&str, usize)]) -> BTreeMap<String, usize> {
    entries.iter().map(|&(k, v)| (k.to_string(), v)).collect()
}

#[test]
fn ratchet_passes_equal_counts() {
    let m = map(&[("a/x.rs", 3), ("a/y.rs", 1)]);
    assert!(compare("r", &m, &m).is_empty());
}

#[test]
fn ratchet_fails_growth_and_names_file_and_counts() {
    let messages = compare("r", &map(&[("a/x.rs", 5)]), &map(&[("a/x.rs", 3)]));
    assert_eq!(messages.len(), 1);
    let m = &messages[0];
    assert!(m.contains("a/x.rs") && m.contains('5') && m.contains('3'), "{m}");
}

#[test]
fn ratchet_asks_to_lower_a_shrunk_baseline() {
    let messages = compare("r", &map(&[("a/x.rs", 2)]), &map(&[("a/x.rs", 3)]));
    assert_eq!(messages.len(), 1);
    assert!(messages[0].contains("lower the baseline to 2"), "{}", messages[0]);
}

#[test]
fn ratchet_reports_a_moved_file_by_name() {
    let messages = compare("r", &map(&[("a/y.rs", 7)]), &map(&[("a/x.rs", 7)]));
    assert!(
        messages.iter().any(|m| m.contains("moved?") && m.contains("a/y.rs")),
        "{messages:?}"
    );
}

#[test]
fn ratchet_keys_use_forward_slashes() {
    let root = source_scan::workspace_root();
    let path = root.join("crates").join("a").join("b.rs");
    assert_eq!(rel_key(&root, &path), "crates/a/b.rs");
}

#[test]
fn line_count_ignores_crlf() {
    assert_eq!(line_count("a\r\nb\r\n"), 2);
    assert_eq!(line_count("a\nb"), 2);
}

#[test]
fn live_code_skips_strings_comments_and_test_items() {
    let src = "let s = \"x.unwrap()\"; // y.unwrap()\n#[cfg(test)] mod t { fn f() { z.unwrap(); } }\nfn g() { w.unwrap(); }";
    let code = live_code(src);
    assert_eq!(code.len(), src.len());
    assert_eq!(code.matches(".unwrap()").count(), 1);
}
