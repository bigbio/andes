use std::process::Command;

#[test]
fn train_help_shows_source_and_gbdt() {
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(["train", "--help"])
        .output()
        .unwrap();
    let h = String::from_utf8_lossy(&out.stdout);
    assert!(
        h.contains("--source"),
        "train --help must document --source:\n{h}"
    );
    assert!(
        h.contains("--gbdt"),
        "train --help must document --gbdt:\n{h}"
    );
}

#[test]
fn train_from_search_subcommand_exists() {
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(["train-from-search", "--help"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "train-from-search subcommand must exist"
    );
}
