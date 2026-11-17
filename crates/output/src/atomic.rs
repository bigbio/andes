use std::path::Path;

pub(crate) use crate::ms1_table::partial_path as partial;

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
