//! Every stderr line of andes goes through the `log` facade and the one logger of the
//! binary (`src/bin/andes/logging.rs`), so the output has a single style. This test fails
//! when non-test code under `crates/*/src` invokes `eprintln!` or `eprint!`, with any
//! delimiter. Comments, string literals and `#[cfg(test)]` items are ignored.

#[path = "support/source_scan.rs"]
mod source_scan;

use std::path::Path;

use source_scan::{blank, macro_calls, rust_files, test_item_ranges, Keep};

#[test]
fn no_crate_prints_to_stderr_outside_the_logger() {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let mut files = Vec::new();
    for entry in std::fs::read_dir(&crates).expect("read crates/") {
        let src = entry.expect("crate dir").path().join("src");
        if src.is_dir() {
            files.extend(rust_files(&src));
        }
    }
    assert!(
        !files.is_empty(),
        "no sources found under {}",
        crates.display()
    );

    let mut offenders = Vec::new();
    for file in &files {
        let text = std::fs::read_to_string(file).expect("read source");
        let code = blank(&text, Keep::Code);
        let tests = test_item_ranges(&code);
        for name in ["eprintln", "eprint"] {
            for at in macro_calls(&code, name) {
                if tests.iter().any(|r| r.contains(&at)) {
                    continue;
                }
                let line = text[..at].matches('\n').count() + 1;
                offenders.push(format!("{}:{line}", file.display()));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "print through log::{{error,warn,info,debug}}! instead of eprintln!/eprint!:\n{}",
        offenders.join("\n")
    );
}
