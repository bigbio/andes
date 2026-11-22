//! Structure checks on the contributor rules: every skill points to its `docs/dev` page,
//! `AGENTS.md` lists every page, and `CONTRIBUTING.md` points to the rules.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

#[path = "support/source_scan.rs"]
mod source_scan;

#[path = "support/ratchet.rs"]
mod ratchet;

use ratchet::{baseline, compare, write_baseline};
use source_scan::{line_count, live_code, macro_calls, rel_key, source_files};

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

const MAX_FILE_LINES: usize = 1500;

/// For every source file, `measure` applied to its text, keyed by workspace-relative path.
/// Files with a zero measure stay in the map so the ratchet can tell a cleaned file from a
/// deleted one.
fn per_file(measure: impl Fn(&str) -> usize) -> BTreeMap<String, usize> {
    let root = workspace_root();
    source_files()
        .iter()
        .map(|path| (rel_key(&root, path), measure(&read(path))))
        .collect()
}

fn oversized_file_lines() -> BTreeMap<String, usize> {
    per_file(|src| {
        let lines = line_count(src);
        if lines > MAX_FILE_LINES {
            lines
        } else {
            0
        }
    })
}

#[test]
fn no_source_file_grows_past_the_budget() {
    let messages = compare("file_lines", &oversized_file_lines(), &baseline("file_lines"));
    assert!(
        messages.is_empty(),
        "source files over {MAX_FILE_LINES} lines may not grow; split them instead:\n{}",
        messages.join("\n")
    );
}

#[test]
#[ignore = "writes crates/andes/tests/ratchets/file_lines.tsv"]
fn write_file_lines_baseline() {
    write_baseline("file_lines", &oversized_file_lines());
}

fn unwrap_expect_counts() -> BTreeMap<String, usize> {
    per_file(|src| {
        let code = live_code(src);
        code.matches(".unwrap()").count() + code.matches(".expect(").count()
    })
}

#[test]
fn no_new_unwrap_or_expect_in_live_code() {
    let messages = compare("unwrap_expect", &unwrap_expect_counts(), &baseline("unwrap_expect"));
    assert!(
        messages.is_empty(),
        "return an error or document the invariant instead of unwrap/expect:\n{}",
        messages.join("\n")
    );
}

#[test]
#[ignore = "writes crates/andes/tests/ratchets/unwrap_expect.tsv"]
fn write_unwrap_expect_baseline() {
    write_baseline("unwrap_expect", &unwrap_expect_counts());
}

/// Files that may use print macros, with the reason. Everything else logs through `log`.
const PRINT_ALLOWED: &[(&str, &str)] = &[(
    "crates/andes/src/bin/andes/logging.rs",
    "the logger itself: the one place that writes log lines to stderr",
)];

#[test]
fn no_print_macros_in_live_code() {
    let root = workspace_root();
    let mut found = Vec::new();
    for path in source_files() {
        let key = rel_key(&root, &path);
        if PRINT_ALLOWED.iter().any(|&(allowed, _)| allowed == key) {
            continue;
        }
        let code = live_code(&read(&path));
        for name in ["println", "eprintln", "print", "eprint", "dbg"] {
            let n = macro_calls(&code, name).len();
            if n > 0 {
                found.push(format!("{key}: {n} x {name}!"));
            }
        }
    }
    assert!(
        found.is_empty(),
        "use log:: macros (docs/dev/clean-logs.md), or add the file to PRINT_ALLOWED with the command whose output it is:\n{}",
        found.join("\n")
    );
}

/// 1 when the file's first non-blank line is not a `//!` doc; inner attributes before it are
/// allowed.
fn lacks_module_doc(src: &str) -> usize {
    let first = src
        .lines()
        .map(str::trim)
        .find(|l| !l.is_empty() && !l.starts_with("#!["));
    usize::from(!first.is_some_and(|l| l.starts_with("//!")))
}

fn module_doc_counts() -> BTreeMap<String, usize> {
    per_file(lacks_module_doc)
}

#[test]
fn every_new_module_has_a_doc() {
    let messages = compare("module_docs", &module_doc_counts(), &baseline("module_docs"));
    assert!(
        messages.is_empty(),
        "start each source file with a `//!` doc of 2-4 lines:\n{}",
        messages.join("\n")
    );
}

#[test]
#[ignore = "writes crates/andes/tests/ratchets/module_docs.tsv"]
fn write_module_docs_baseline() {
    write_baseline("module_docs", &module_doc_counts());
}

fn citation_registry() -> Vec<serde_yaml::Value> {
    let text = read(&workspace_root().join("docs/citations.yaml"));
    serde_yaml::from_str(&text).expect("docs/citations.yaml is a YAML list")
}

fn entry_text<'a>(entry: &'a serde_yaml::Value, field: &str) -> &'a str {
    entry
        .get(field)
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or("")
        .trim()
}

fn is_citation_key(key: &str) -> bool {
    !key.is_empty()
        && key
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_')
}

fn cited_keys(text: &str) -> Vec<String> {
    let mut keys = Vec::new();
    let mut rest = text;
    while let Some(at) = rest.find("[cite:") {
        rest = &rest[at + "[cite:".len()..];
        if let Some(end) = rest.find(']') {
            if is_citation_key(&rest[..end]) {
                keys.push(rest[..end].to_string());
            }
        }
    }
    keys
}

#[test]
fn every_cited_key_is_in_the_registry() {
    let known: Vec<String> = citation_registry()
        .iter()
        .map(|e| entry_text(e, "key").to_string())
        .collect();
    let mut files = source_files();
    files.extend(source_scan::text_files(&["docs"], &["md"]));
    let mut missing = Vec::new();
    for file in files {
        for key in cited_keys(&read(&file)) {
            if !known.contains(&key) {
                missing.push(format!("{}: [cite:{key}]", file.display()));
            }
        }
    }
    assert!(
        missing.is_empty(),
        "citation keys missing from docs/citations.yaml:\n{}",
        missing.join("\n")
    );
}

#[test]
fn registry_entries_are_complete_and_unique() {
    let mut seen = std::collections::BTreeSet::new();
    let mut problems = Vec::new();
    for entry in citation_registry() {
        let key = entry_text(&entry, "key").to_string();
        if !is_citation_key(&key) {
            problems.push(format!("key {key:?} is empty or not [a-z0-9_]+"));
        }
        if !seen.insert(key.clone()) {
            problems.push(format!("{key}: duplicate key"));
        }
        for field in ["authors", "title"] {
            if entry_text(&entry, field).is_empty() {
                problems.push(format!("{key}: empty {field}"));
            }
        }
        let year = entry
            .get("year")
            .and_then(serde_yaml::Value::as_u64)
            .unwrap_or(0);
        if !(1900..=2100).contains(&year) {
            problems.push(format!("{key}: year {year} outside 1900-2100"));
        }
        if entry_text(&entry, "doi").is_empty() && entry_text(&entry, "url").is_empty() {
            problems.push(format!("{key}: needs a doi or a url"));
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

const SITE_LITERALS: [&str; 6] = [
    "/homes/",
    "/hps/",
    "/nfs/",
    "pst_prd",
    "codon-slurm",
    "hl-codon",
];

/// The site-specific patterns found in `text`, one entry per occurrence.
///
/// A path under `/home/<name>/`, `/Users/<name>/` or `/Data/<name>/` is a hit (`/home/runner/` is the
/// GitHub CI home and is allowed), as is any of [`SITE_LITERALS`].
fn site_hits(text: &str) -> Vec<String> {
    let mut hits = Vec::new();
    for literal in SITE_LITERALS {
        for _ in text.matches(literal) {
            hits.push(literal.to_string());
        }
    }
    for prefix in ["/home/", "/Users/", "/Data/"] {
        for (at, _) in text.match_indices(prefix) {
            let rest = &text[at + prefix.len()..];
            let word = rest
                .bytes()
                .take_while(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-'))
                .count();
            if word == 0 || rest.as_bytes().get(word) != Some(&b'/') {
                continue;
            }
            if prefix == "/home/" && rest.get(..word) == Some("runner") {
                continue;
            }
            hits.push(format!("{prefix}<name>/"));
        }
    }
    hits
}

#[test]
fn site_patterns_ignore_public_urls() {
    assert!(site_hits(
        "ftp.pride.ebi.ac.uk/pride/data https://www.ebi.ac.uk/pride doi.org/10.1/x"
    )
    .is_empty());
    assert_eq!(site_hits("/home/alice/x").len(), 1);
    assert_eq!(site_hits("/Users/bob/y").len(), 1);
    assert_eq!(site_hits("/hps/nobackup/x").len(), 1);
    assert_eq!(site_hits("pst_prd").len(), 1);
    assert_eq!(site_hits("/Data/carol/raw").len(), 1);
    assert!(site_hits("/home/runner/work").is_empty());
}

/// Committed text names no cluster, host or user path: paths are placeholders or
/// arguments. This file is exempt because its tests spell the patterns out.
#[test]
fn no_site_identifiers_in_committed_text() {
    let root = source_scan::workspace_root();
    let this_file = "crates/andes/tests/structure.rs";
    let mut files = source_scan::text_files(
        &["crates", "docs", "scripts", "website", ".claude", ".github"],
        &["rs", "md", "py", "sh", "yml", "yaml", "toml", "tsv", "json"],
    );
    for name in [
        "README.md",
        "DOCS.md",
        "TRAIN.md",
        "AGENTS.md",
        "CLAUDE.md",
        "CONTRIBUTING.md",
        "ARCHITECTURE.md",
    ] {
        let path = root.join(name);
        if path.is_file() {
            files.push(path);
        }
    }
    let mut found = Vec::new();
    for file in files {
        let key = rel_key(&root, &file);
        if key == this_file {
            continue;
        }
        let Ok(text) = fs::read_to_string(&file) else {
            continue;
        };
        for (n, line) in text.lines().enumerate() {
            for hit in site_hits(line) {
                found.push(format!("{key}:{}: {hit}", n + 1));
            }
        }
    }
    assert!(
        found.is_empty(),
        "site identifiers in committed text; replace each with a placeholder or an argument:\n{}",
        found.join("\n")
    );
}
