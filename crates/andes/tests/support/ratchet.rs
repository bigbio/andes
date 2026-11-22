//! Per-file ratchets: a rule counts violations per file, the committed baseline under
//! `crates/andes/tests/ratchets/<rule>.tsv` records today's counts, and a count may only go
//! down. Lines are `path<TAB>count`; `#` starts a comment.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::PathBuf;

use super::source_scan::workspace_root;

fn baseline_path(rule: &str) -> PathBuf {
    workspace_root()
        .join("crates")
        .join("andes")
        .join("tests")
        .join("ratchets")
        .join(format!("{rule}.tsv"))
}

/// The committed baseline of `rule`; a missing file is an empty map.
pub fn baseline(rule: &str) -> BTreeMap<String, usize> {
    let path = baseline_path(rule);
    let Ok(text) = std::fs::read_to_string(&path) else {
        return BTreeMap::new();
    };
    let mut out = BTreeMap::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parsed = line
            .split_once('\t')
            .and_then(|(file, count)| Some((file.trim(), count.trim().parse::<usize>().ok()?)));
        match parsed {
            Some((file, count)) => {
                out.insert(file.to_string(), count);
            }
            None => panic!(
                "{}:{}: expected `path<TAB>count`, got `{line}`",
                path.display(),
                n + 1
            ),
        }
    }
    out
}

/// One message per violation of `rule`: growth, a count without a baseline line, a baseline
/// that can be lowered, or a baseline line for a file that no longer exists.
pub fn compare(
    rule: &str,
    counts: &BTreeMap<String, usize>,
    base: &BTreeMap<String, usize>,
) -> Vec<String> {
    let mut out = Vec::new();
    for (file, &count) in counts {
        match base.get(file) {
            Some(&allowed) if count > allowed => out.push(format!(
                "{rule}: {file} has {count}, baseline is {allowed}; remove the new occurrences (do not raise the baseline)"
            )),
            Some(&allowed) if count == 0 && allowed > 0 => out.push(format!(
                "{rule}: {file} has 0, below its baseline {allowed}; remove its baseline line"
            )),
            Some(&allowed) if count < allowed => out.push(format!(
                "{rule}: {file} has {count}, below its baseline {allowed}; lower the baseline to {count}"
            )),
            None if count > 0 => out.push(format!(
                "{rule}: {file} has {count} and no baseline line (baseline 0); remove the occurrences. If this file was split from or renamed from another, move or split its baseline line, keeping the total from rising"
            )),
            _ => {}
        }
    }
    for (file, &allowed) in base {
        if counts.contains_key(file) {
            continue;
        }
        let moved = counts
            .iter()
            .find(|&(new, &c)| c == allowed && !base.contains_key(new))
            .map(|(new, _)| format!(" (moved? {new} has exactly {allowed})"))
            .unwrap_or_default();
        out.push(format!(
            "{rule}: baseline line `{file}` ({allowed}) names a file that no longer exists; remove or move its baseline line{moved}"
        ));
    }
    out
}

/// The files whose count in `counts` is above their baseline, or that have a count above 0
/// and no baseline line. A regeneration may only lower or remove baseline lines, so any entry
/// here blocks it.
pub fn raised(counts: &BTreeMap<String, usize>, base: &BTreeMap<String, usize>) -> Vec<String> {
    let mut out = Vec::new();
    for (file, &count) in counts {
        match base.get(file) {
            Some(&allowed) if count > allowed => {
                out.push(format!("{file}: {count} is above its baseline {allowed}"));
            }
            None if count > 0 => out.push(format!("{file}: {count} has no baseline line")),
            _ => {}
        }
    }
    out
}

/// Rewrites the baseline of `rule` from `counts`: entries above 0 only, sorted by path, under
/// the existing leading `#` header lines. Panics, listing the files, when any count would rise
/// or a new line would appear; a baseline can only be lowered or shortened.
pub fn write_baseline(rule: &str, counts: &BTreeMap<String, usize>) {
    let refused = raised(counts, &baseline(rule));
    assert!(
        refused.is_empty(),
        "refusing to raise the {rule} baseline; fix the code instead:\n{}",
        refused.join("\n")
    );
    write_entries(rule, counts);
}

/// Writes `counts` as the baseline of `rule` without the no-raise check: entries above 0 only,
/// sorted by path, under the existing leading `#` header lines.
pub fn write_entries(rule: &str, counts: &BTreeMap<String, usize>) {
    let path = baseline_path(rule);
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).expect("create ratchets dir");
    }
    let existing = std::fs::read_to_string(&path).unwrap_or_default();
    let mut text: String = existing
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    for (file, c) in counts.iter().filter(|&(_, &c)| c > 0) {
        text.push_str(&format!("{file}\t{c}\n"));
    }
    std::fs::write(&path, text).unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
}
