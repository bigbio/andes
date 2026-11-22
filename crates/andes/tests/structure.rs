//! Structure checks on the contributor rules: every skill points to its `docs/dev` page,
//! `AGENTS.md` lists every page, and `CONTRIBUTING.md` points to the rules.

use std::fs;
use std::path::{Path, PathBuf};

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
