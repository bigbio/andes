//! Tests of the directory defaults (`--dir` or the output directory).

use super::*;

#[test]
fn the_flag_wins() {
    assert_eq!(
        dir_or_output_dir(Some(Path::new("/shared")), Path::new("out/run.pin")),
        PathBuf::from("/shared")
    );
}

#[test]
fn the_default_is_the_output_directory() {
    assert_eq!(
        dir_or_output_dir(None, Path::new("results/run.pin")),
        PathBuf::from("results")
    );
    assert_eq!(
        dir_or_output_dir(None, Path::new("/data/out/run")),
        PathBuf::from("/data/out")
    );
}

#[test]
fn a_bare_output_name_means_the_working_directory() {
    assert_eq!(
        dir_or_output_dir(None, Path::new("run.pin")),
        PathBuf::from(".")
    );
}
