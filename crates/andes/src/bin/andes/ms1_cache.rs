//! Private scratch storage for a run's MS1 scans, as an MS1 table
//! (`output::ms1_table`) in a temporary file under `--temp-dir`. No run stays resident while
//! Percolator trains; quantification loads one run at a time, and the
//! `--quant-dir` artifact `<stem>.ms1.parquet` is a copy of this file.
use std::path::Path;

use output::ms1_table::{partial_path, read_ms1, write_ms1};
use quant::ms1_index::Ms1RunIndex;
use tempfile::NamedTempFile;

pub(crate) struct Ms1Cache {
    file: NamedTempFile,
}

impl Ms1Cache {
    /// Writes `index` of run `run` to a new temporary MS1 table in `dir`
    /// (created if missing).
    pub(super) fn store(index: &Ms1RunIndex, run: &str, dir: &Path) -> Result<Self, String> {
        crate::cli::ensure_scratch_dir(dir)?;
        let file = tempfile::Builder::new()
            .prefix("andes-ms1-")
            .suffix(".ms1.parquet")
            .tempfile_in(dir)
            .map_err(|e| crate::cli::scratch_dir_error(dir, &e))?;
        write_ms1(file.as_file(), run, index)
            .map_err(|e| format!("{}: {e}", file.path().display()))?;
        Ok(Self { file })
    }

    pub(super) fn load(&self) -> Result<Ms1RunIndex, String> {
        let path = self.file.path();
        let file = self
            .file
            .reopen()
            .map_err(|e| format!("{}: {e}", path.display()))?;
        read_ms1(file, path)
    }

    /// Copies the table to `path`, through `<path>.partial`, so the file
    /// appears only once complete. Every error names `path`.
    pub(super) fn persist(&self, path: &Path) -> Result<(), String> {
        let partial = partial_path(path);
        let copied = std::fs::copy(self.file.path(), &partial)
            .and_then(|_| std::fs::rename(&partial, path))
            .map_err(|e| format!("{}: {e}", path.display()));
        if copied.is_err() {
            let _ = std::fs::remove_file(&partial);
        }
        copied
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::scan::Ms1Scan;
    use output::ms1_table::read_ms1_table;

    fn index() -> Ms1RunIndex {
        Ms1RunIndex::new(vec![
            Ms1Scan {
                rt: 3.123456789,
                peaks: vec![(600.123456789, 1234.5678), (700.0, -0.0)],
                id: "scan=2".into(),
                faims_cv: Some(-45.0),
            },
            Ms1Scan {
                rt: 1.0,
                peaks: vec![],
                ..Ms1Scan::default()
            },
        ])
    }

    fn assert_same(a: &Ms1RunIndex, b: &Ms1RunIndex) {
        assert_eq!(a.len(), b.len());
        for i in 0..a.len() {
            assert_eq!(a.rt(i).to_bits(), b.rt(i).to_bits());
            assert_eq!(a.peaks(i).len(), b.peaks(i).len());
            for (p, q) in a.peaks(i).iter().zip(b.peaks(i)) {
                assert_eq!(p.0.to_bits(), q.0.to_bits());
                assert_eq!(p.1.to_bits(), q.1.to_bits());
            }
            assert_eq!(a.scan(i).id, b.scan(i).id);
            assert_eq!(a.scan(i).faims_cv, b.scan(i).faims_cv);
        }
    }

    #[test]
    fn lossless_repeatable_round_trip_and_cleanup() {
        let idx = index();
        let scratch = tempfile::tempdir().unwrap();
        let cache = Ms1Cache::store(&idx, "run", scratch.path()).unwrap();
        let path = cache.file.path().to_owned();
        assert_eq!(path.parent(), Some(scratch.path()));
        for _ in 0..2 {
            assert_same(&cache.load().unwrap(), &idx);
        }
        drop(cache);
        assert!(!path.exists());
    }

    #[test]
    fn empty_and_truncated_caches() {
        let scratch = tempfile::tempdir().unwrap();
        let empty = Ms1Cache::store(&Ms1RunIndex::default(), "e", scratch.path()).unwrap();
        assert!(empty.load().unwrap().is_empty());
        let cache = Ms1Cache::store(&index(), "run", scratch.path()).unwrap();
        cache.file.as_file().set_len(25).unwrap();
        let e = cache.load().unwrap_err();
        assert!(e.contains(".ms1.parquet"), "{e}");
    }

    #[test]
    fn persisted_copy_is_the_same_table() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.ms1.parquet");
        let idx = index();
        let cache = Ms1Cache::store(&idx, "run", dir.path()).unwrap();
        cache.persist(&path).unwrap();
        assert!(!partial_path(&path).exists());
        assert_eq!(
            std::fs::read(&path).unwrap(),
            std::fs::read(cache.file.path()).unwrap()
        );
        assert_same(&read_ms1_table(&path).unwrap(), &idx);
        let bad = dir.path().join("missing").join("x.ms1.parquet");
        let e = cache.persist(&bad).unwrap_err();
        assert!(e.contains("x.ms1.parquet"), "{e}");
    }

    /// A scratch directory that cannot be created is an error naming the flag.
    #[test]
    fn unwritable_scratch_dir_is_an_error_naming_the_flag() {
        let dir = tempfile::tempdir().unwrap();
        let blocker = dir.path().join("blocker");
        std::fs::write(&blocker, b"x").unwrap();
        let e = Ms1Cache::store(&index(), "run", &blocker.join("sub"))
            .err()
            .expect("must fail");
        assert!(e.contains("--temp-dir"), "{e}");
    }
}
