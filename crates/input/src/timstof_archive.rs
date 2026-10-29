//! Archived Bruker `.d` run folders (`feature = "timstof"`).
//!
//! Repositories ship timsTOF runs as archives of the `.d` folder: `X.d.zip`,
//! `X.d.tar`, `X.d.tar.gz` / `X.d.tgz`, and `X.d.gz` (a gzipped tar). The
//! reader needs the TDF database and its binary blob as files it can open and
//! memory-map, so an archive is unpacked into a directory the caller chooses
//! (the run's temporary directory), its run folder is located (the folder
//! holding `analysis.tdf` and `analysis.tdf_bin`, at the archive root or
//! nested inside it), and that folder is moved to `<dir>/<X>.d` so the run
//! keeps the archive's name.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// Error unpacking an archived `.d`.
#[derive(Debug, thiserror::Error)]
pub enum DArchiveError {
    #[error("cannot read the archive '{path}': {message}")]
    Corrupt { path: String, message: String },
    #[error("cannot unpack '{path}' into '{dir}': {source}")]
    Io {
        path: String,
        dir: String,
        source: std::io::Error,
    },
    #[error(
        "the archive '{path}' holds no Bruker .d run folder (no folder with analysis.tdf \
         and analysis.tdf_bin)"
    )]
    NoRun { path: String },
    #[error("the archive '{path}' holds {count} Bruker .d run folders; archive one run per file")]
    SeveralRuns { path: String, count: usize },
}

/// Archive layouts of a `.d` folder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DArchiveKind {
    Zip,
    Tar,
    /// A gzipped tar (`.d.tar.gz`, `.d.tgz`, `.d.gz`).
    TarGz,
}

/// File-name suffixes of archived `.d` folders, longest first, and the run
/// name each leaves when stripped.
const SUFFIXES: [(&str, DArchiveKind); 5] = [
    (".d.tar.gz", DArchiveKind::TarGz),
    (".d.zip", DArchiveKind::Zip),
    (".d.tar", DArchiveKind::Tar),
    (".d.tgz", DArchiveKind::TarGz),
    (".d.gz", DArchiveKind::TarGz),
];

/// The archive kind and run name (`X` of `X.d.zip`) of an archived `.d`,
/// from its file name (case-insensitive), or `None` for any other path.
pub fn d_archive_kind(path: &Path) -> Option<(DArchiveKind, String)> {
    let name = path.file_name()?.to_str()?;
    let lower = name.to_ascii_lowercase();
    SUFFIXES.iter().find_map(|(suffix, kind)| {
        let run = &name[..name.len().checked_sub(suffix.len())?];
        (lower.ends_with(suffix) && !run.is_empty()).then(|| (*kind, run.to_string()))
    })
}

/// Unpack the archived `.d` at `archive` under `dir` and return the run
/// folder, `<dir>/<X>.d`. `dir` must exist; the archive is unpacked into a
/// fresh sub-folder of it first.
pub fn unpack_d_archive(archive: &Path, dir: &Path) -> Result<PathBuf, DArchiveError> {
    let corrupt = |message: String| DArchiveError::Corrupt {
        path: archive.display().to_string(),
        message,
    };
    let io = |target: &Path, source: std::io::Error| DArchiveError::Io {
        path: archive.display().to_string(),
        dir: target.display().to_string(),
        source,
    };
    let (kind, run) = d_archive_kind(archive).ok_or_else(|| {
        corrupt("not a .d.zip, .d.tar, .d.tar.gz, .d.tgz or .d.gz file name".into())
    })?;
    let staging = dir.join(format!("{run}.d.unpack"));
    std::fs::create_dir_all(&staging).map_err(|e| io(&staging, e))?;
    let file = File::open(archive).map_err(|e| io(&staging, e))?;
    match kind {
        DArchiveKind::Zip => {
            let mut zip = zip::ZipArchive::new(BufReader::new(file))
                .map_err(|e| corrupt(format!("not a readable zip archive ({e})")))?;
            zip.extract(&staging)
                .map_err(|e| corrupt(format!("zip entry could not be extracted ({e})")))?;
        }
        DArchiveKind::Tar => unpack_tar(BufReader::new(file), &staging).map_err(corrupt)?,
        DArchiveKind::TarGz => {
            let mut file = file;
            let mut magic = [0u8; 2];
            let gz = file.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b];
            if !gz {
                return Err(corrupt("not a gzip file".into()));
            }
            file.seek(SeekFrom::Start(0)).map_err(|e| io(&staging, e))?;
            let decoder = flate2::read::MultiGzDecoder::new(BufReader::new(file));
            unpack_tar(decoder, &staging).map_err(corrupt)?;
        }
    }
    let mut runs = Vec::new();
    find_run_folders(&staging, 0, &mut runs);
    let found = match runs.len() {
        0 => {
            return Err(DArchiveError::NoRun {
                path: archive.display().to_string(),
            })
        }
        1 => runs.pop().expect("one run"),
        count => {
            return Err(DArchiveError::SeveralRuns {
                path: archive.display().to_string(),
                count,
            })
        }
    };
    let target = dir.join(format!("{run}.d"));
    std::fs::rename(&found, &target).map_err(|e| io(&target, e))?;
    Ok(target)
}

/// Unpack a tar stream into `dir` (entries that would land outside `dir` are
/// refused by the tar reader).
fn unpack_tar<R: Read>(reader: R, dir: &Path) -> Result<(), String> {
    let mut archive = tar::Archive::new(reader);
    archive
        .unpack(dir)
        .map_err(|e| format!("not a readable tar archive ({e})"))
}

/// Deepest nesting searched for the run folder below the unpack root.
const MAX_DEPTH: usize = 4;

/// Whether `dir` is a TDF run folder.
fn is_run_folder(dir: &Path) -> bool {
    dir.join("analysis.tdf").is_file() && dir.join("analysis.tdf_bin").is_file()
}

/// Collect the run folders at or below `dir` (not below a run folder, and
/// skipping the `__MACOSX` resource-fork folders some zip tools add).
fn find_run_folders(dir: &Path, depth: usize, out: &mut Vec<PathBuf>) {
    if is_run_folder(dir) {
        out.push(dir.to_path_buf());
        return;
    }
    if depth >= MAX_DEPTH {
        return;
    }
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    let mut subdirs: Vec<PathBuf> = entries
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.is_dir() && p.file_name().is_some_and(|n| n != "__MACOSX"))
        .collect();
    subdirs.sort();
    for sub in subdirs {
        find_run_folders(&sub, depth + 1, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// The two files of a run folder, as (name, bytes).
    const RUN_FILES: [(&str, &[u8]); 2] =
        [("analysis.tdf", b"sqlite"), ("analysis.tdf_bin", b"frames")];

    fn tar_bytes(prefix: &str) -> Vec<u8> {
        let mut builder = tar::Builder::new(Vec::new());
        for (name, data) in RUN_FILES {
            let mut header = tar::Header::new_gnu();
            header.set_size(data.len() as u64);
            header.set_mode(0o644);
            header.set_cksum();
            builder
                .append_data(&mut header, format!("{prefix}{name}"), data)
                .unwrap();
        }
        builder.into_inner().unwrap()
    }

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap()
    }

    fn zip_bytes(prefix: &str) -> Vec<u8> {
        let mut w = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
        let opts = zip::write::SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated);
        for (name, data) in RUN_FILES {
            w.start_file(format!("{prefix}{name}"), opts).unwrap();
            w.write_all(data).unwrap();
        }
        w.finish().unwrap().into_inner()
    }

    fn unpack(name: &str, bytes: &[u8]) -> (tempfile::TempDir, Result<PathBuf, DArchiveError>) {
        let tmp = tempfile::tempdir().unwrap();
        let archive = tmp.path().join(name);
        std::fs::write(&archive, bytes).unwrap();
        let out = tmp.path().join("out");
        std::fs::create_dir(&out).unwrap();
        let r = unpack_d_archive(&archive, &out);
        (tmp, r)
    }

    fn assert_run(r: &Result<PathBuf, DArchiveError>, run: &str) {
        let dir = r.as_ref().expect("unpacked");
        assert_eq!(
            dir.file_name().unwrap().to_str().unwrap(),
            format!("{run}.d")
        );
        assert_eq!(
            std::fs::read(dir.join("analysis.tdf_bin")).unwrap(),
            b"frames"
        );
        assert!(is_run_folder(dir));
    }

    #[test]
    fn archive_kinds_come_from_the_file_name() {
        let k = |n: &str| d_archive_kind(Path::new(n));
        assert_eq!(k("/x/A.d.zip"), Some((DArchiveKind::Zip, "A".into())));
        assert_eq!(k("A.D.ZIP"), Some((DArchiveKind::Zip, "A".into())));
        assert_eq!(k("A.d.tar"), Some((DArchiveKind::Tar, "A".into())));
        assert_eq!(k("A.d.tar.gz"), Some((DArchiveKind::TarGz, "A".into())));
        assert_eq!(k("A.d.tgz"), Some((DArchiveKind::TarGz, "A".into())));
        assert_eq!(k("A.d.gz"), Some((DArchiveKind::TarGz, "A".into())));
        assert_eq!(k("A.d"), None);
        assert_eq!(k("A.mzML.gz"), None);
        assert_eq!(k("A.zip"), None);
        assert_eq!(k(".d.zip"), None);
    }

    #[test]
    fn a_zip_with_the_run_folder_inside_unpacks() {
        let (_t, r) = unpack("run1.d.zip", &zip_bytes("run1.d/"));
        assert_run(&r, "run1");
    }

    #[test]
    fn a_flat_zip_of_the_folder_contents_unpacks() {
        let (_t, r) = unpack("run1.d.zip", &zip_bytes(""));
        assert_run(&r, "run1");
    }

    #[test]
    fn a_run_folder_nested_deeper_is_found_and_takes_the_archive_name() {
        let (_t, r) = unpack("run1.d.zip", &zip_bytes("export/2020/other.d/"));
        assert_run(&r, "run1");
    }

    #[test]
    fn tar_tar_gz_tgz_and_gz_unpack() {
        let (_t, r) = unpack("a.d.tar", &tar_bytes("a.d/"));
        assert_run(&r, "a");
        let (_t, r) = unpack("b.d.tar.gz", &gzip(&tar_bytes("b.d/")));
        assert_run(&r, "b");
        let (_t, r) = unpack("c.d.tgz", &gzip(&tar_bytes("")));
        assert_run(&r, "c");
        let (_t, r) = unpack("d.d.gz", &gzip(&tar_bytes("d.d/d.d/")));
        assert_run(&r, "d");
    }

    #[test]
    fn corrupt_and_foreign_archives_are_errors() {
        let (_t, r) = unpack("bad.d.zip", b"this is not a zip");
        assert!(matches!(r, Err(DArchiveError::Corrupt { .. })), "{r:?}");
        let (_t, r) = unpack("bad.d.tar.gz", b"plain text");
        assert!(matches!(r, Err(DArchiveError::Corrupt { .. })), "{r:?}");
        let truncated = gzip(&tar_bytes("x.d/"));
        let (_t, r) = unpack("bad.d.tgz", &truncated[..truncated.len() / 2]);
        assert!(r.is_err(), "{r:?}");
        // A well-formed archive without a run folder.
        let mut w = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
        w.start_file("readme.txt", zip::write::SimpleFileOptions::default())
            .unwrap();
        w.write_all(b"hello").unwrap();
        let (_t, r) = unpack("empty.d.zip", &w.finish().unwrap().into_inner());
        assert!(matches!(r, Err(DArchiveError::NoRun { .. })), "{r:?}");
    }

    #[test]
    fn two_run_folders_are_an_error() {
        let mut bytes = tar::Builder::new(Vec::new());
        for prefix in ["a.d/", "b.d/"] {
            for (name, data) in RUN_FILES {
                let mut header = tar::Header::new_gnu();
                header.set_size(data.len() as u64);
                header.set_mode(0o644);
                header.set_cksum();
                bytes
                    .append_data(&mut header, format!("{prefix}{name}"), data)
                    .unwrap();
            }
        }
        let (_t, r) = unpack("two.d.tar", &bytes.into_inner().unwrap());
        assert!(
            matches!(r, Err(DArchiveError::SeveralRuns { count: 2, .. })),
            "{r:?}"
        );
    }

    #[test]
    fn macos_resource_forks_are_ignored() {
        let mut w = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
        let opts = zip::write::SimpleFileOptions::default();
        for prefix in ["r.d/", "__MACOSX/r.d/"] {
            for (name, data) in RUN_FILES {
                w.start_file(format!("{prefix}{name}"), opts).unwrap();
                w.write_all(data).unwrap();
            }
        }
        let (_t, r) = unpack("r.d.zip", &w.finish().unwrap().into_inner());
        assert_run(&r, "r");
    }
}
