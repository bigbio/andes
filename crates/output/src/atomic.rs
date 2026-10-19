//! Atomic output files: a writer writes `<path>.partial` and renames it to
//! `path` once complete, so a failed or interrupted run never leaves a
//! truncated file under the final name.

use std::path::{Path, PathBuf};

/// The temporary path a writer of `path` writes first.
pub(crate) fn partial(path: &Path) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(".partial");
    PathBuf::from(p)
}

/// Move the complete `partial` file to `path`.
pub(crate) fn commit(partial: &Path, path: &Path) -> std::io::Result<()> {
    std::fs::rename(partial, path)
}
