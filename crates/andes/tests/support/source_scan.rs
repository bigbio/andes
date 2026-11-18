//! A small Rust source scanner shared by the repository guard tests
//! (`single_log_style.rs`, `crates/search/tests/env_registry.rs`).
//!
//! It blanks comments and literals with spaces, keeping every byte offset and newline, so
//! the guards match code tokens only and can count braces exactly.

#![allow(dead_code)]

use std::ops::Range;
use std::path::{Path, PathBuf};

/// What [`blank`] keeps.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Keep {
    /// Blank comments, string and char literals: only code tokens remain.
    Code,
    /// Blank comments only: string literals stay.
    Strings,
}

fn is_ident(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

/// The source with comments (and, for [`Keep::Code`], string and char literals) replaced
/// by spaces. Newlines are kept and the result has the same length as `src`.
pub fn blank(src: &str, keep: Keep) -> String {
    let b = src.as_bytes();
    let mut out = b.to_vec();
    let mut i = 0;
    let wipe = |out: &mut Vec<u8>, from: usize, to: usize| {
        for byte in &mut out[from..to] {
            if *byte != b'\n' {
                *byte = b' ';
            }
        }
    };
    while i < b.len() {
        let prev_ident = i > 0 && is_ident(b[i - 1]);
        if b[i] == b'/' && b.get(i + 1) == Some(&b'/') {
            let end = b[i..]
                .iter()
                .position(|&c| c == b'\n')
                .map_or(b.len(), |p| i + p);
            wipe(&mut out, i, end);
            i = end;
        } else if b[i] == b'/' && b.get(i + 1) == Some(&b'*') {
            let mut depth = 0usize;
            let mut j = i;
            while j < b.len() {
                if b[j] == b'/' && b.get(j + 1) == Some(&b'*') {
                    depth += 1;
                    j += 2;
                } else if b[j] == b'*' && b.get(j + 1) == Some(&b'/') {
                    depth -= 1;
                    j += 2;
                    if depth == 0 {
                        break;
                    }
                } else {
                    j += 1;
                }
            }
            wipe(&mut out, i, j);
            i = j;
        } else if !prev_ident
            && (b[i] == b'r' || (b[i] == b'b' && b.get(i + 1) == Some(&b'r')))
            && {
                let start = if b[i] == b'b' { i + 2 } else { i + 1 };
                let hashes = b[start..].iter().take_while(|&&c| c == b'#').count();
                b.get(start + hashes) == Some(&b'"')
            }
        {
            let start = if b[i] == b'b' { i + 2 } else { i + 1 };
            let hashes = b[start..].iter().take_while(|&&c| c == b'#').count();
            let body = start + hashes + 1;
            let closing: Vec<u8> = std::iter::once(b'"')
                .chain(std::iter::repeat_n(b'#', hashes))
                .collect();
            let end = b[body..]
                .windows(closing.len())
                .position(|w| w == closing.as_slice())
                .map_or(b.len(), |p| body + p + closing.len());
            if keep == Keep::Code {
                wipe(&mut out, i, end);
            }
            i = end;
        } else if b[i] == b'"' {
            let mut j = i + 1;
            while j < b.len() && b[j] != b'"' {
                j += if b[j] == b'\\' { 2 } else { 1 };
            }
            let end = (j + 1).min(b.len());
            if keep == Keep::Code {
                wipe(&mut out, i, end);
            }
            i = end;
        } else if b[i] == b'\'' {
            let end = if b.get(i + 1) == Some(&b'\\') {
                b.get(i + 3..)
                    .and_then(|rest| rest.iter().position(|&c| c == b'\''))
                    .map(|p| i + 3 + p + 1)
            } else {
                let ch_len = src[i + 1..].chars().next().map_or(1, char::len_utf8);
                (b.get(i + 1 + ch_len) == Some(&b'\'')).then_some(i + 1 + ch_len + 1)
            };
            match end {
                Some(end) => {
                    if keep == Keep::Code {
                        wipe(&mut out, i, end);
                    }
                    i = end;
                }
                None => i += 1,
            }
        } else {
            i += 1;
        }
    }
    String::from_utf8(out).expect("blanking keeps UTF-8")
}

/// Byte ranges of the items marked `#[cfg(test)]` (a `mod` or `fn` with a body), found on
/// a [`Keep::Code`] blanking of the source.
pub fn test_item_ranges(code: &str) -> Vec<Range<usize>> {
    let b = code.as_bytes();
    let mut out = Vec::new();
    let mut from = 0;
    while let Some(p) = code[from..].find("#[cfg(test)]") {
        let at = from + p;
        let after = at + "#[cfg(test)]".len();
        let open = b[after..].iter().position(|&c| c == b'{' || c == b';');
        match open.map(|o| after + o) {
            Some(o) if b[o] == b'{' => {
                let mut depth = 0usize;
                let mut end = b.len();
                for (k, &c) in b.iter().enumerate().skip(o) {
                    if c == b'{' {
                        depth += 1;
                    } else if c == b'}' {
                        depth -= 1;
                        if depth == 0 {
                            end = k + 1;
                            break;
                        }
                    }
                }
                out.push(at..end);
                from = end;
            }
            _ => from = after,
        }
    }
    out
}

/// Whether `path` is test-only by location (`tests/`, `examples/` or `benches/`).
pub fn is_test_path(path: &Path) -> bool {
    path.components().any(|c| {
        matches!(
            c.as_os_str().to_str(),
            Some("tests" | "examples" | "benches")
        )
    })
}

/// Every `.rs` file under `dir`, sorted.
pub fn rust_files(dir: &Path) -> Vec<PathBuf> {
    fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
        let Ok(entries) = std::fs::read_dir(dir) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                if path.file_name().is_some_and(|n| n == "target") {
                    continue;
                }
                walk(&path, out);
            } else if path.extension().is_some_and(|e| e == "rs") {
                out.push(path);
            }
        }
    }
    let mut out = Vec::new();
    walk(dir, &mut out);
    out.sort();
    out
}

/// Byte offsets where a macro named `name` is invoked with any delimiter.
pub fn macro_calls(code: &str, name: &str) -> Vec<usize> {
    let b = code.as_bytes();
    let needle = format!("{name}!");
    let mut out = Vec::new();
    let mut from = 0;
    while let Some(p) = code[from..].find(&needle) {
        let at = from + p;
        let after = at + needle.len();
        let next = b[after..].iter().find(|c| !c.is_ascii_whitespace());
        let ident_before = at > 0 && is_ident(b[at - 1]);
        if !ident_before && matches!(next, Some(b'(' | b'[' | b'{')) {
            out.push(at);
        }
        from = after;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comments_strings_and_chars_are_blanked_with_offsets_kept() {
        let src = "let a = \"x{\"; // }\nlet c = '{'; /* } */ let r = r#\"}\"#;";
        let code = blank(src, Keep::Code);
        assert_eq!(code.len(), src.len());
        assert!(!code.contains('{') && !code.contains('}'));
        assert!(code.contains("let a") && code.contains('\n'));
        let strings = blank(src, Keep::Strings);
        assert!(strings.contains("\"x{\"") && !strings.contains("// }"));
    }

    #[test]
    fn an_escaped_quote_char_ends_at_its_own_closing_quote() {
        let src = "let q = ['\\'', 'x']; eprintln!(\"{}\", q[0]);";
        let code = blank(src, Keep::Code);
        assert!(!code.contains('\'') && !code.contains('x'));
        assert_eq!(macro_calls(&code, "eprintln").len(), 1);
    }

    #[test]
    fn lifetimes_are_not_char_literals() {
        let src = "fn f<'a>(x: &'a str) -> char { 'z' }";
        let code = blank(src, Keep::Code);
        assert!(code.contains("<'a>") && !code.contains("'z'"));
    }

    #[test]
    fn a_test_module_ends_at_its_real_closing_brace() {
        let src = "#[cfg(test)]\nmod tests { fn t() { let _ = '}'; } }\nfn shipped() {}";
        let code = blank(src, Keep::Code);
        let ranges = test_item_ranges(&code);
        assert_eq!(ranges.len(), 1);
        assert!(src[ranges[0].end..].contains("fn shipped"));
    }

    #[test]
    fn a_commented_attribute_does_not_hide_code() {
        let src = "// #[cfg(test)] mod tests\nfn shipped() { eprintln!(\"x\"); }";
        let code = blank(src, Keep::Code);
        assert!(test_item_ranges(&code).is_empty());
        assert_eq!(macro_calls(&code, "eprintln").len(), 1);
    }

    #[test]
    fn every_macro_delimiter_is_found_and_names_in_strings_are_not() {
        let code = blank(
            "eprintln!(\"a\"); eprintln![\"b\"]; eprintln! {\"c\"}; let s = \"eprintln!(x)\"; my_eprintln!(1);",
            Keep::Code,
        );
        assert_eq!(macro_calls(&code, "eprintln").len(), 3);
    }
}
