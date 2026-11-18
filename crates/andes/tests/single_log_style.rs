//! Every stderr line of andes goes through the `log` facade and the one logger of the
//! binary (`src/bin/andes/logging.rs`), so the output has a single style. This test fails
//! when non-test code under `crates/*/src` prints to stderr directly.

use std::path::{Path, PathBuf};

const FORBIDDEN: [&str; 2] = ["eprintln!(", "eprint!("];

/// Source with every `#[cfg(test)] mod name { ... }` block removed, so unit tests may
/// still print.
fn without_test_modules(src: &str) -> String {
    let mut out = String::with_capacity(src.len());
    let mut rest = src;
    while let Some(at) = rest.find("#[cfg(test)]") {
        out.push_str(&rest[..at]);
        let after = &rest[at..];
        let Some(open) = after.find('{') else {
            rest = "";
            break;
        };
        let head = &after[..open];
        if !head.contains("mod ") || head.contains(';') {
            out.push_str("#[cfg(test)]");
            rest = &after["#[cfg(test)]".len()..];
            continue;
        }
        let mut depth = 0usize;
        let mut end = after.len();
        for (i, c) in after[open..].char_indices() {
            match c {
                '{' => depth += 1,
                '}' => {
                    depth -= 1;
                    if depth == 0 {
                        end = open + i + 1;
                        break;
                    }
                }
                _ => {}
            }
        }
        rest = &after[end..];
    }
    out.push_str(rest);
    out
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(dir).expect("read_dir") {
        let path = entry.expect("dir entry").path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

#[test]
fn no_crate_prints_to_stderr_outside_the_logger() {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let mut files = Vec::new();
    for entry in std::fs::read_dir(&crates).expect("read crates/") {
        let src = entry.expect("crate dir").path().join("src");
        if src.is_dir() {
            rust_files(&src, &mut files);
        }
    }
    files.sort();
    assert!(
        !files.is_empty(),
        "no sources found under {}",
        crates.display()
    );

    let mut offenders = Vec::new();
    for file in &files {
        let text = std::fs::read_to_string(file).expect("read source");
        let code = without_test_modules(&text);
        for line in code.lines() {
            let trimmed = line.trim_start();
            if trimmed.starts_with("//") {
                continue;
            }
            if FORBIDDEN.iter().any(|m| trimmed.contains(m)) {
                offenders.push(format!("{}: {}", file.display(), trimmed));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "print through log::{{error,warn,info,debug}}! instead of eprintln!/eprint!:\n{}",
        offenders.join("\n")
    );
}

#[test]
fn test_modules_are_skipped_and_other_code_is_kept() {
    let src = "fn a() {}\n#[cfg(test)]\nmod tests {\n    fn t() { eprintln!(\"{}\", 1); }\n}\nfn b() {}\n";
    let code = without_test_modules(src);
    assert!(!code.contains("eprintln!"));
    assert!(code.contains("fn a()") && code.contains("fn b()"));
}
