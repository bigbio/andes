//! Physical constants are defined once, in `crates/model/src/mass.rs`, so every algorithm
//! uses the same proton, isotope spacing and small-molecule masses. This test fails when
//! non-test code elsewhere under `crates/*/src` types one of them as a literal.
//!
//! Unimod modification deltas (for example `-17.026549` for Gln->pyro-Glu) are format
//! constants and stay literal: only proton, ¹H and isotope-spacing literals anywhere, and
//! constants *defined* as water, ammonia or carbon monoxide, are rejected.

use std::path::{Path, PathBuf};

/// Literal prefixes (digits only, underscores removed) of values that must come from
/// `model::mass`.
const ANYWHERE: [(&str, &str); 3] = [
    ("1.00727", "model::mass::PROTON"),
    ("1.00782", "model::mass::H"),
    ("1.00335", "model::mass::ISOTOPE"),
];

/// Prefixes rejected when they are the value of a `const`.
const AS_CONST: [(&str, &str); 3] = [
    ("18.01", "model::mass::H2O"),
    ("17.02", "model::mass::NH3"),
    ("27.99", "model::mass::CO"),
];

/// Literal copies kept on purpose, each with the reason. Replacing one changes search
/// results, so it is done in its own change with benchmark evidence (entrapment FDP on the
/// reference datasets), not as a refactor.
const KNOWN_EXCEPTIONS: [(&str, &str, &str); 2] = [
    (
        "search/src/match_engine.rs",
        "const H2O_LOSS_DA: f64 = 18.0106;",
        "rounded water loss; model::mass::H2O moves RawScore on ~2.4% of low-res PSMs",
    ),
    (
        "search/src/match_engine.rs",
        "const NH3_LOSS_DA: f64 = 17.0265;",
        "rounded ammonia loss; model::mass::NH3 moves RawScore on ~2.4% of low-res PSMs",
    ),
];

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

/// The line with comments and string literals blanked, and digit separators removed.
fn code_of(line: &str) -> String {
    let mut out = String::with_capacity(line.len());
    let mut in_str = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if in_str {
            if c == '\\' {
                chars.next();
            } else if c == '"' {
                in_str = false;
            }
            continue;
        }
        match c {
            '"' => in_str = true,
            '/' if chars.peek() == Some(&'/') => break,
            '_' => {}
            _ => out.push(c),
        }
    }
    out
}

/// The numeric literals of a code line, each with whether it is the value of a `const`.
fn literals(code: &str) -> Vec<(String, bool)> {
    let is_const = code.trim_start().starts_with("const ")
        || code.trim_start().starts_with("pub const ")
        || code.trim_start().starts_with("pub(crate) const ");
    let bytes = code.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let starts = bytes[i].is_ascii_digit()
            && (i == 0 || !(bytes[i - 1].is_ascii_alphanumeric() || bytes[i - 1] == b'.'));
        if starts {
            let j = (i..bytes.len())
                .find(|&k| !(bytes[k].is_ascii_digit() || bytes[k] == b'.'))
                .unwrap_or(bytes.len());
            out.push((code[i..j].to_string(), is_const));
            i = j;
        } else {
            i += 1;
        }
    }
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
fn physical_constants_come_only_from_model_mass() {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let home = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/mass.rs");
    let mut files = Vec::new();
    for entry in std::fs::read_dir(&crates).expect("read crates/") {
        let src = entry.expect("crate dir").path().join("src");
        if src.is_dir() {
            rust_files(&src, &mut files);
        }
    }
    files.sort();
    assert!(!files.is_empty(), "no sources under {}", crates.display());

    let mut offenders = Vec::new();
    for file in &files {
        if file.ends_with("model/src/mass.rs") || file == &home {
            continue;
        }
        let text = std::fs::read_to_string(file).expect("read source");
        for line in without_test_modules(&text).lines() {
            let excepted = KNOWN_EXCEPTIONS.iter().any(|(path, text, _)| {
                file.to_string_lossy().replace('\\', "/").ends_with(path) && line.trim() == *text
            });
            if excepted {
                continue;
            }
            let code = code_of(line);
            for (lit, is_const) in literals(&code) {
                let hit = ANYWHERE
                    .iter()
                    .chain(AS_CONST.iter().filter(|_| is_const))
                    .find(|(prefix, _)| lit.starts_with(prefix));
                if let Some((_, use_instead)) = hit {
                    offenders.push(format!(
                        "{}: {} — use {use_instead}",
                        file.display(),
                        line.trim()
                    ));
                }
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "physical constants typed outside crates/model/src/mass.rs:\n{}",
        offenders.join("\n")
    );
}

#[test]
fn every_known_exception_still_exists() {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    for (path, text, why) in KNOWN_EXCEPTIONS {
        let src = std::fs::read_to_string(crates.join(path)).expect("read exception file");
        assert!(
            src.lines().any(|l| l.trim() == text),
            "{path} no longer has `{text}` ({why}); remove it from KNOWN_EXCEPTIONS"
        );
    }
}

#[test]
fn comments_strings_and_unimod_deltas_are_not_flagged() {
    let code = code_of(r#"let s = "1.00335/z"; // proton 1.007276"#);
    assert!(literals(&code).is_empty());
    let delta = code_of(r#"m("Gln->pyro-Glu", -17.026549, &["Q"]),"#);
    assert!(literals(&delta).iter().all(|(_, c)| !c));
    let defined = code_of("const NH3: f64 = 17.026_549;");
    assert_eq!(literals(&defined), vec![("17.026549".to_string(), true)]);
}
