//! Tests of `.d` archive handling.

use super::inputs::check_distinct_d_archive_runs;
use std::path::PathBuf;

#[test]
fn archives_of_one_run_name_are_refused_before_unpacking() {
    let ok: Vec<PathBuf> = ["a.d.zip", "b.d.tar", "c.mzML"].map(PathBuf::from).to_vec();
    assert!(check_distinct_d_archive_runs(&ok).is_ok());
    let dup: Vec<PathBuf> = ["x/run.d.zip", "y/run.d.tar"].map(PathBuf::from).to_vec();
    let err = check_distinct_d_archive_runs(&dup).unwrap_err();
    assert!(
        err.contains("x/run.d.zip") && err.contains("y/run.d.tar"),
        "{err}"
    );
}
