use std::path::{Path, PathBuf};

pub(crate) use self::partial_path as partial;

/// The sibling file a table is written to before it is renamed into place.
pub fn partial_path(path: &Path) -> PathBuf {
    let mut partial = path.as_os_str().to_owned();
    partial.push(".partial");
    PathBuf::from(partial)
}

/// Like `write_atomic` for writers that report `String` errors. `write` receives the
/// partial path and must prefix its own errors; a failed rename is reported as
/// `{path}: {error}`. The partial file is removed on any failure.
pub(crate) fn write_atomic_str<T>(
    path: &Path,
    write: impl FnOnce(&Path) -> Result<T, String>,
) -> Result<T, String> {
    let partial = partial_path(path);
    let written = write(&partial).and_then(|v| {
        std::fs::rename(&partial, path)
            .map(|()| v)
            .map_err(|e| format!("{}: {e}", path.display()))
    });
    if written.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}

pub fn write_atomic<T>(
    path: &Path,
    write: impl FnOnce(&Path) -> std::io::Result<T>,
) -> std::io::Result<T> {
    let tmp = partial(path);
    match write(&tmp) {
        Ok(v) => {
            std::fs::rename(&tmp, path)?;
            Ok(v)
        }
        Err(e) => {
            let _ = std::fs::remove_file(&tmp);
            Err(e)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_string_write_renames_on_success_and_cleans_up_on_failure() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.parquet");
        let err = write_atomic_str(&path, |tmp| {
            std::fs::write(tmp, "half").map_err(|e| e.to_string())?;
            Err::<(), _>("stop".to_string())
        });
        assert_eq!(err.unwrap_err(), "stop");
        assert!(!path.exists() && !partial(&path).exists());
        write_atomic_str(&path, |tmp| {
            std::fs::write(tmp, "whole").map_err(|e| e.to_string())
        })
        .unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "whole");
        assert!(!partial(&path).exists());
        let missing = dir.path().join("none").join("out.parquet");
        let err = write_atomic_str(&missing, |_| Ok(())).unwrap_err();
        assert!(
            err.starts_with(&format!("{}: ", missing.display())),
            "{err}"
        );
    }

    #[test]
    fn a_failed_write_leaves_neither_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.tsv");
        let err = write_atomic(&path, |tmp| {
            std::fs::write(tmp, "half")?;
            Err::<(), _>(std::io::Error::other("stop"))
        });
        assert!(err.is_err());
        assert!(!path.exists() && !partial(&path).exists());
        write_atomic(&path, |tmp| std::fs::write(tmp, "whole")).unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "whole");
        assert!(!partial(&path).exists());
    }
}
