//! What the search gathers for the run report while it runs: the wall time and peak memory of
//! each stage (`StageClock`), the SHA-256 of every input file, hashed on a thread of its own from
//! the start of the run (`InputHasher`), and the command line. `methods` writes the Methods
//! paragraph and its references; `build` turns the run's results into an `output::RunReport`.

pub(crate) mod build;
pub(crate) mod methods;

use std::path::PathBuf;
use std::time::Instant;

use output::{InputFile, StageRow};

pub(crate) use crate::memory::peak_rss_bytes;

/// The stages of a run in the order they ran, each closed with its wall time and the process
/// peak memory at its end.
pub(crate) struct StageClock {
    rows: Vec<StageRow>,
    current: Option<(String, Instant)>,
}

impl StageClock {
    pub(crate) fn new() -> Self {
        Self {
            rows: Vec::new(),
            current: None,
        }
    }

    /// Closes the open stage, if any, and opens `name`.
    pub(crate) fn stage(&mut self, name: &str) {
        self.close();
        self.current = Some((name.to_string(), Instant::now()));
    }

    fn close(&mut self) {
        if let Some((name, start)) = self.current.take() {
            self.rows.push(StageRow {
                name,
                seconds: start.elapsed().as_secs_f64(),
                peak_rss_bytes: peak_rss_bytes(),
            });
        }
    }

    /// Closes the open stage and returns every stage in order.
    pub(crate) fn finish(mut self) -> Vec<StageRow> {
        self.close();
        self.rows
    }
}

/// The SHA-256 of the input files, computed on a thread of its own so that it overlaps the
/// search; `join` waits for it.
pub(crate) struct InputHasher {
    files: Vec<(String, PathBuf)>,
    handle: std::thread::JoinHandle<Vec<InputFile>>,
}

impl InputHasher {
    /// Starts hashing `files`, given as `(role, path)`, in order.
    pub(crate) fn start(files: Vec<(String, PathBuf)>) -> Self {
        let todo = files.clone();
        let handle = std::thread::spawn(move || todo.iter().map(|(r, p)| hash_one(r, p)).collect());
        Self { files, handle }
    }

    /// The hashed files in the order given. A file that could not be read carries its error;
    /// when the thread failed, every file says so.
    pub(crate) fn join(self) -> Vec<InputFile> {
        match self.handle.join() {
            Ok(files) => files,
            Err(_) => self
                .files
                .into_iter()
                .map(|(role, path)| InputFile {
                    role,
                    path: path.display().to_string(),
                    bytes: None,
                    sha256: Err("hashing thread failed".to_string()),
                })
                .collect(),
        }
    }
}

fn hash_one(role: &str, path: &std::path::Path) -> InputFile {
    let meta = std::fs::metadata(path);
    let (bytes, sha256) = match &meta {
        Ok(m) if m.is_dir() => (None, Err("a directory; not hashed".to_string())),
        Ok(m) => (
            Some(m.len()),
            search::decoy::file_sha256(path).map_err(|e| e.to_string()),
        ),
        Err(e) => (None, Err(e.to_string())),
    };
    InputFile {
        role: role.to_string(),
        path: path.display().to_string(),
        bytes,
        sha256,
    }
}

/// The command line the process was started with, each argument quoted when it holds
/// whitespace or a quote.
pub(crate) fn command_line() -> String {
    command_line_of(std::env::args_os().map(|a| a.to_string_lossy().into_owned()))
}

pub(crate) fn command_line_of(args: impl Iterator<Item = String>) -> String {
    args.map(|a| {
        if a.is_empty() || a.chars().any(|c| c.is_whitespace() || c == '\'' || c == '"') {
            format!("'{}'", a.replace('\'', "'\\''"))
        } else {
            a
        }
    })
    .collect::<Vec<_>>()
    .join(" ")
}

#[cfg(test)]
#[path = "report/tests/clock.rs"]
mod clock_tests;

#[cfg(test)]
#[path = "report/tests/methods.rs"]
mod methods_tests;

#[cfg(test)]
#[path = "report/tests/build.rs"]
mod build_tests;
