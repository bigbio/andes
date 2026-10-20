//! Percolator runner: resolve a backend (explicit path → `$PATH` → docker),
//! invoke Percolator on a PIN, and parse the target-PSM result TSV back into a
//! `SpecId → PercolatorPsm` map so the QPX writer can join PEP/q-value onto the
//! per-PSM rows.
//!
//! andes does not reimplement Percolator; it orchestrates the real binary (or
//! the pinned biocontainers docker image) with
//! `--seed 42 --results-psms … --decoy-results-psms … --weights … --only-psms <PIN>`.
//! The result TSV is parsed by HEADER NAME (not column index) so it survives the
//! column-offset differences between engines (e.g. a a comparison search engine-style extra `FileName`
//! column).

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Pinned biocontainers Percolator image used for the docker fallback.
pub const DEFAULT_PERCOLATOR_IMAGE: &str = "quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2";

/// One parsed Percolator target-PSM result row.
#[derive(Debug, Clone, PartialEq)]
pub struct PercolatorPsm {
    /// `PSMId` column — equals the andes PIN `SpecId` byte-for-byte.
    pub psm_id: String,
    /// `q-value` column (lower is better).
    pub q_value: f64,
    /// `posterior_error_probability` column (PEP; lower is better).
    pub pep: f64,
    /// `score` column (higher is better); `None` when the file has no such
    /// column or the value is missing or not a finite number.
    pub score: Option<f64>,
    /// `peptide` column (Percolator's flanked peptide string).
    pub peptide: String,
    /// `proteinIds` column joined with `\t` collapsed to a single string (the
    /// rest-of-line after `proteinIds`, tab-separated proteins re-joined by `;`).
    pub proteins: String,
}

/// How Percolator will be invoked.
#[derive(Debug, Clone, PartialEq)]
pub enum PercolatorBackend {
    /// A native binary at this path (explicit `--percolator-bin` or found on `$PATH`).
    Binary(PathBuf),
    /// The docker fallback with the given image tag.
    Docker { image: String },
}

impl PercolatorBackend {
    /// Human-readable provenance string for the run log.
    pub fn describe(&self) -> String {
        match self {
            PercolatorBackend::Binary(p) => format!("binary {}", p.display()),
            PercolatorBackend::Docker { image } => format!("docker {image}"),
        }
    }
}

/// Errors from backend resolution, invocation, or result parsing.
#[derive(Debug)]
pub enum PercolatorError {
    /// No backend could be resolved (no `--percolator-bin`, none on `$PATH`, no docker).
    NoBackend,
    /// The requested explicit binary path does not exist / is not executable.
    BinaryNotFound(PathBuf),
    /// Docker was requested/needed but the `docker` CLI is unavailable.
    DockerUnavailable,
    /// The Percolator process failed (non-zero exit). Carries captured stderr/stdout tail.
    Run {
        backend: String,
        status: String,
        log_tail: String,
    },
    /// Spawning the process failed.
    Spawn(std::io::Error),
    /// I/O while reading the result file.
    Io(std::io::Error),
    /// The result TSV was missing a required column.
    MissingColumn {
        column: &'static str,
        header: String,
    },
    /// A result row had fewer columns than the header requires (truncated/corrupt).
    MalformedRow { line: String },
    /// A `q-value`/PEP field was not a finite number — propagating it as NaN
    /// would silently corrupt FDR filtering downstream.
    NonFiniteScore { column: &'static str, value: String },
    /// The same `PSMId` appeared in more than one result row. `PSMId` is the
    /// join key into the idparquet and must be unique; a duplicate means the
    /// last-write-wins map would silently drop a PSM.
    DuplicatePsmId { psm_id: String },
}

impl std::fmt::Display for PercolatorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PercolatorError::NoBackend => write!(
                f,
                "no Percolator backend: pass --percolator-bin <path>, put `percolator` on PATH, \
                 or enable the docker fallback with --percolator-docker (needs the `docker` CLI \
                 and image {DEFAULT_PERCOLATOR_IMAGE})"
            ),
            PercolatorError::BinaryNotFound(p) => {
                write!(
                    f,
                    "--percolator-bin {} not found or not executable",
                    p.display()
                )
            }
            PercolatorError::DockerUnavailable => {
                write!(
                    f,
                    "docker fallback requested but the `docker` CLI is not available"
                )
            }
            PercolatorError::Run {
                backend,
                status,
                log_tail,
            } => write!(
                f,
                "Percolator ({backend}) failed ({status}); last log lines:\n{log_tail}"
            ),
            PercolatorError::Spawn(e) => write!(f, "failed to spawn Percolator: {e}"),
            PercolatorError::Io(e) => write!(f, "reading Percolator results: {e}"),
            PercolatorError::MissingColumn { column, header } => {
                write!(
                    f,
                    "Percolator result missing `{column}` column; header was: {header}"
                )
            }
            PercolatorError::MalformedRow { line } => {
                write!(f, "Percolator result row has too few columns: {line}")
            }
            PercolatorError::NonFiniteScore { column, value } => {
                write!(
                    f,
                    "Percolator result `{column}` is not a finite number: {value:?}"
                )
            }
            PercolatorError::DuplicatePsmId { psm_id } => {
                write!(
                    f,
                    "Percolator result has duplicate PSMId `{psm_id}` (join key must be unique)"
                )
            }
        }
    }
}

impl std::error::Error for PercolatorError {}

impl From<std::io::Error> for PercolatorError {
    fn from(e: std::io::Error) -> Self {
        PercolatorError::Io(e)
    }
}

/// Resolve which Percolator backend to use.
///
/// Order: (1) `explicit_bin` (`--percolator-bin`) — hard error if it does not
/// exist; (2) `percolator` on `$PATH`; (3) docker fallback when `force_docker`
/// is set OR nothing else was found AND the `docker` CLI is available. Returns
/// [`PercolatorError::NoBackend`] when none applies.
pub fn resolve_backend(
    explicit_bin: Option<&Path>,
    force_docker: bool,
    image: &str,
) -> Result<PercolatorBackend, PercolatorError> {
    if let Some(bin) = explicit_bin {
        if bin.exists() {
            return Ok(PercolatorBackend::Binary(bin.to_path_buf()));
        }
        return Err(PercolatorError::BinaryNotFound(bin.to_path_buf()));
    }
    if force_docker {
        if docker_available() {
            return Ok(PercolatorBackend::Docker {
                image: image.to_string(),
            });
        }
        return Err(PercolatorError::DockerUnavailable);
    }
    if let Some(p) = which_on_path("percolator") {
        return Ok(PercolatorBackend::Binary(p));
    }
    if docker_available() {
        return Ok(PercolatorBackend::Docker {
            image: image.to_string(),
        });
    }
    Err(PercolatorError::NoBackend)
}

/// Probe `$PATH` for an executable named `bin`. Returns the first hit.
fn which_on_path(bin: &str) -> Option<PathBuf> {
    let paths = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&paths) {
        let cand = dir.join(bin);
        if cand.is_file() {
            return Some(cand);
        }
    }
    None
}

/// True if the `docker` CLI can be invoked (`docker --version` exits cleanly).
fn docker_available() -> bool {
    Command::new("docker")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Output paths produced by a Percolator run (all in the PIN's parent dir).
struct ResultPaths {
    target: PathBuf,
    decoy: PathBuf,
    weights: PathBuf,
}

impl ResultPaths {
    fn new(out_dir: &Path, stem: &str) -> Self {
        ResultPaths {
            target: out_dir.join(format!("{stem}.percolator.target.psms.txt")),
            decoy: out_dir.join(format!("{stem}.percolator.decoy.psms.txt")),
            weights: out_dir.join(format!("{stem}.percolator.weights.txt")),
        }
    }

    /// Delete any result files left over from a previous run, so that a silent
    /// non-write by Percolator can never be mis-read as this run's output (the
    /// success check in `run_percolator` trusts `target.exists()`). Absent files
    /// are fine; a real removal failure is an error (we must not leave a stale
    /// file we would then read).
    fn clear_stale(&self) -> std::io::Result<()> {
        for p in [&self.target, &self.decoy, &self.weights] {
            match std::fs::remove_file(p) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
}

/// Build the argument vector for a NATIVE Percolator invocation (binary path
/// case). Mirrors `run_percolator_docker.sh`: `--seed 42`, results/decoy/weights,
/// `--only-psms`, then any user `extra_args`, then the PIN path last.
fn native_args(pin: &Path, out: &ResultPaths, extra_args: &[String]) -> Vec<String> {
    let mut args = vec![
        "--seed".to_string(),
        "42".to_string(),
        "--results-psms".to_string(),
        out.target.to_string_lossy().into_owned(),
        "--decoy-results-psms".to_string(),
        out.decoy.to_string_lossy().into_owned(),
        "--weights".to_string(),
        out.weights.to_string_lossy().into_owned(),
        "--only-psms".to_string(),
    ];
    args.extend(extra_args.iter().cloned());
    args.push(pin.to_string_lossy().into_owned());
    args
}

/// Build the full `docker run …` argument vector. The PIN's directory is bind
/// mounted read-only at `/pin` and the output dir at `/out`, exactly as the
/// reference script does, so percolator reads/writes in place. Result paths are
/// rewritten to their `/out` equivalents inside the container. Both mount
/// sources are made absolute: docker reads a relative or empty `-v` source as
/// a named volume (or rejects it), not as a host directory.
fn docker_args(
    image: &str,
    pin: &Path,
    out_dir: &Path,
    out: &ResultPaths,
    extra_args: &[String],
) -> std::io::Result<Vec<String>> {
    let pin_dir = host_dir(pin.parent().unwrap_or_else(|| Path::new(".")))?;
    let out_dir = host_dir(out_dir)?;
    let pin_name = pin
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let in_out = |p: &Path| {
        let name = p
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        format!("/out/{name}")
    };
    let mut args = vec![
        "run".to_string(),
        "--rm".to_string(),
        "--platform".to_string(),
        "linux/amd64".to_string(),
        "-v".to_string(),
        format!("{}:/pin:ro", pin_dir.display()),
        "-v".to_string(),
        format!("{}:/out", out_dir.display()),
        image.to_string(),
        "percolator".to_string(),
        "--seed".to_string(),
        "42".to_string(),
        "--results-psms".to_string(),
        in_out(&out.target),
        "--decoy-results-psms".to_string(),
        in_out(&out.decoy),
        "--weights".to_string(),
        in_out(&out.weights),
        "--only-psms".to_string(),
    ];
    args.extend(extra_args.iter().cloned());
    args.push(format!("/pin/{pin_name}"));
    Ok(args)
}

/// `dir` as an absolute host path; the empty parent of a bare file name is
/// the current directory.
fn host_dir(dir: &Path) -> std::io::Result<PathBuf> {
    if dir.as_os_str().is_empty() {
        std::env::current_dir()
    } else {
        std::path::absolute(dir)
    }
}

/// Run Percolator on `pin` using `backend`, writing result files into the PIN's
/// parent directory, then parse the TARGET PSMs into a `SpecId → PercolatorPsm`
/// map. `extra_args` are appended verbatim (from `--percolator-args`).
///
/// Returns the parsed map plus the resolved result file paths (so the caller can
/// log/keep them). The PEP/q join uses the map keyed by `PSMId` (== PIN SpecId).
pub fn run_percolator(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
) -> Result<HashMap<String, PercolatorPsm>, PercolatorError> {
    run_percolator_with_decoys(backend, pin, extra_args).map(|(targets, _)| targets)
}

/// [`run_percolator`], also returning the `PSMId`s of the decoy result file.
/// With target-decoy competition (`-Y`) Percolator keeps one PSM per
/// scan/`ExpMass` group in either file; the decoy ids let the QPX join tell a
/// target that lost the competition from one whose join drifted.
pub fn run_percolator_with_decoys(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
) -> Result<(HashMap<String, PercolatorPsm>, HashSet<String>), PercolatorError> {
    let out_dir = pin.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
    let out = result_paths(pin);
    // Remove any stale outputs from a previous run BEFORE launching, so the
    // `out.target.exists()` success check below can only pass on a file this
    // run actually wrote (guards against reading a previous run's results).
    out.clear_stale()?;

    let output = match backend {
        PercolatorBackend::Binary(bin) => Command::new(bin)
            .args(native_args(pin, &out, extra_args))
            .output()
            .map_err(PercolatorError::Spawn)?,
        PercolatorBackend::Docker { image } => Command::new("docker")
            .args(docker_args(image, pin, &out_dir, &out, extra_args)?)
            .output()
            .map_err(PercolatorError::Spawn)?,
    };

    if !output.status.success() || !out.target.exists() {
        let mut log = String::from_utf8_lossy(&output.stderr).into_owned();
        if log.trim().is_empty() {
            log = String::from_utf8_lossy(&output.stdout).into_owned();
        }
        let log_tail = tail_lines(&log, 20);
        return Err(PercolatorError::Run {
            backend: backend.describe(),
            status: output.status.to_string(),
            log_tail,
        });
    }

    let text = std::fs::read_to_string(&out.target)?;
    let targets = parse_psm_results(&text)?;
    let decoys = match std::fs::read_to_string(&out.decoy) {
        Ok(text) => parse_psm_results(&text)?.into_keys().collect(),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => HashSet::new(),
        Err(e) => return Err(e.into()),
    };
    Ok((targets, decoys))
}

/// The result files `run_percolator` writes for `pin`, next to it.
fn result_paths(pin: &Path) -> ResultPaths {
    let out_dir = pin.parent().unwrap_or_else(|| Path::new("."));
    let stem = pin
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    ResultPaths::new(out_dir, &stem)
}

/// The decoy-PSM result file a [`run_percolator`] call on `pin` writes; parse
/// it with [`parse_psm_results`].
pub fn decoy_results_path(pin: &Path) -> PathBuf {
    result_paths(pin).decoy
}

/// Keep the last `n` lines of `s`.
fn tail_lines(s: &str, n: usize) -> String {
    let lines: Vec<&str> = s.lines().collect();
    let start = lines.len().saturating_sub(n);
    lines[start..].join("\n")
}

/// Parse a Percolator score field, rejecting anything that is not a finite
/// number (empty, non-numeric, NaN, or ±inf) — such a value must not flow into
/// FDR filtering as a silent NaN.
fn parse_finite(s: &str, column: &'static str) -> Result<f64, PercolatorError> {
    match s.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(v),
        _ => Err(PercolatorError::NonFiniteScore {
            column,
            value: s.to_string(),
        }),
    }
}

/// Parse a Percolator PSM result TSV into a `PSMId → PercolatorPsm` map.
///
/// HEADER-NAME-DRIVEN: the column positions of `PSMId`, `q-value`,
/// `posterior_error_probability`, `peptide`, and `proteinIds` are located by
/// matching the header row, so a leading/extra column (e.g. a comparison search engine's `FileName`)
/// does not shift the parse. `proteinIds` is the LAST logical column; every tab
/// after it is a separate protein id, so the rest of the line is captured and
/// re-joined with `;`.
pub fn parse_psm_results(text: &str) -> Result<HashMap<String, PercolatorPsm>, PercolatorError> {
    let mut reader = PsmResultReader::new(text.as_bytes())?;
    let mut map = HashMap::new();
    while let Some(row) = reader.next_row()? {
        let psm_id = row.psm_id.to_string();
        let psm = PercolatorPsm {
            psm_id: psm_id.clone(),
            q_value: row.q_value,
            pep: row.pep,
            score: row.score,
            peptide: row.peptide.to_string(),
            proteins: row.proteins(),
        };
        if map.insert(psm_id.clone(), psm).is_some() {
            return Err(PercolatorError::DuplicatePsmId { psm_id });
        }
    }
    Ok(map)
}

/// One row of a Percolator PSM result file, borrowed from its line.
#[derive(Debug)]
pub struct PsmResultRow<'a> {
    /// `PSMId` column.
    pub psm_id: &'a str,
    /// `q-value` column.
    pub q_value: f64,
    /// `posterior_error_prob` column.
    pub pep: f64,
    /// `score` column; `None` when absent, unparseable or not finite.
    pub score: Option<f64>,
    /// `peptide` column.
    pub peptide: &'a str,
    /// `proteinIds` and every field after it.
    protein_fields: Vec<&'a str>,
}

impl PsmResultRow<'_> {
    /// The protein fields joined by `;`.
    pub fn proteins(&self) -> String {
        self.protein_fields.join(";")
    }
}

/// Column positions of a Percolator PSM result header.
#[derive(Debug, Clone, Copy)]
struct ResultColumns {
    id: usize,
    q: usize,
    pep: usize,
    peptide: usize,
    proteins: usize,
    score: Option<usize>,
}

impl ResultColumns {
    fn from_header(header: &str) -> Result<Self, PercolatorError> {
        let cols: Vec<&str> = header.split('\t').collect();
        let find = |name: &'static str| -> Result<usize, PercolatorError> {
            cols.iter()
                .position(|c| *c == name)
                .ok_or(PercolatorError::MissingColumn {
                    column: name,
                    header: header.to_string(),
                })
        };
        // Match a column by name-prefix. Percolator 3.7.x emits the PEP column as
        // `posterior_error_prob`; other builds/tools use the full
        // `posterior_error_probability`. A prefix match accepts both.
        let find_prefix = |prefix: &'static str| -> Result<usize, PercolatorError> {
            cols.iter()
                .position(|c| c.starts_with(prefix))
                .ok_or(PercolatorError::MissingColumn {
                    column: prefix,
                    header: header.to_string(),
                })
        };
        Ok(Self {
            id: find("PSMId")?,
            q: find("q-value")?,
            pep: find_prefix("posterior_error_prob")?,
            peptide: find("peptide")?,
            proteins: find("proteinIds")?,
            score: cols.iter().position(|c| *c == "score"),
        })
    }

    fn row<'a>(&self, line: &'a str) -> Result<PsmResultRow<'a>, PercolatorError> {
        let fields: Vec<&str> = line.split('\t').collect();
        // Fail loud on a row shorter than the highest non-protein index we read
        // (truncated/corrupt) — silently skipping it would drop a PSM and bias FDR.
        let max_scalar = self.id.max(self.q).max(self.pep).max(self.peptide);
        if fields.len() <= max_scalar {
            return Err(PercolatorError::MalformedRow {
                line: line.to_string(),
            });
        }
        // A non-finite or unparseable q-value/PEP must be an error, never a
        // silent NaN — NaN compares false in every FDR threshold check.
        let q_value = parse_finite(fields[self.q], "q-value")?;
        let pep = parse_finite(fields[self.pep], "posterior_error_prob")?;
        // The score is optional: a missing, unparseable or non-finite value is
        // no score. Callers that rank by it reject such PSMs themselves.
        let score = self
            .score
            .and_then(|c| fields.get(c))
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite());
        // proteinIds is the trailing column: every field from it on.
        let protein_fields = fields.get(self.proteins..).unwrap_or(&[]).to_vec();
        Ok(PsmResultRow {
            psm_id: fields[self.id],
            q_value,
            pep,
            score,
            peptide: fields[self.peptide],
            protein_fields,
        })
    }
}

/// Reads a Percolator PSM result file one row at a time, so a file of any
/// size is parsed with one line in memory. Lines are split as
/// [`str::lines`] splits them and empty lines are skipped; the columns are
/// found by header name as in [`parse_psm_results`].
pub struct PsmResultReader<R> {
    reader: R,
    line: String,
    columns: ResultColumns,
}

impl<R: std::io::BufRead> PsmResultReader<R> {
    /// Reads the header. Errors: a missing required column, or I/O.
    pub fn new(mut reader: R) -> Result<Self, PercolatorError> {
        let mut line = String::new();
        read_line(&mut reader, &mut line)?;
        let columns = ResultColumns::from_header(&line)?;
        Ok(Self {
            reader,
            line,
            columns,
        })
    }

    /// The next row, `None` at the end of the file. Errors: a malformed row,
    /// a non-finite q-value or PEP, or I/O.
    pub fn next_row(&mut self) -> Result<Option<PsmResultRow<'_>>, PercolatorError> {
        loop {
            if !read_line(&mut self.reader, &mut self.line)? {
                return Ok(None);
            }
            if !self.line.is_empty() {
                break;
            }
        }
        self.columns.row(&self.line).map(Some)
    }
}

/// Reads one line into `line` without its `\n` (and then `\r`) ending;
/// false at the end of the input.
fn read_line<R: std::io::BufRead>(
    reader: &mut R,
    line: &mut String,
) -> Result<bool, PercolatorError> {
    line.clear();
    if reader.read_line(line).map_err(PercolatorError::Io)? == 0 {
        return Ok(false);
    }
    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_pin_header() {
        // Real Percolator 3.7.x header: the PEP column is the abbreviated
        // `posterior_error_prob` (verified against a live percolator:3.7.1 run).
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    title_100_1\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\tsp|P2\n\
                    title_101_1\t1.0\t0.02\t0.03\tK.ANOTHERK.R\tsp|P3\n";
        let m = parse_psm_results(text).unwrap();
        assert_eq!(m.len(), 2);
        let a = &m["title_100_1"];
        assert_eq!(a.q_value, 0.001);
        assert_eq!(a.pep, 0.0005);
        assert_eq!(a.peptide, "R.PEPTIDEK.A");
        // Two trailing protein columns re-joined with ';'.
        assert_eq!(a.proteins, "sp|P1;sp|P2");
        assert_eq!(m["title_101_1"].proteins, "sp|P3");
        assert_eq!(a.score, Some(2.5));
        assert_eq!(m["title_101_1"].score, Some(1.0));
    }

    #[test]
    fn reader_splits_lines_as_str_lines() {
        // CRLF endings, an empty line, and a last line without an ending.
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\r\n\
                    a_1_1\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\tsp|P2\r\n\
                    \r\n\
                    b_1_1\t1.0\t0.02\t0.03\tK.ANOTHERK.R";
        let mut r = PsmResultReader::new(text.as_bytes()).unwrap();
        let a = r.next_row().unwrap().unwrap();
        assert_eq!((a.psm_id, a.score, a.peptide), ("a_1_1", Some(2.5), "R.PEPTIDEK.A"));
        assert_eq!(a.proteins(), "sp|P1;sp|P2");
        let b = r.next_row().unwrap().unwrap();
        assert_eq!((b.psm_id, b.q_value, b.pep), ("b_1_1", 0.02, 0.03));
        assert_eq!(b.proteins(), "");
        assert!(r.next_row().unwrap().is_none());
        let m = parse_psm_results(text).unwrap();
        assert_eq!(m.len(), 2);
        assert_eq!(m["a_1_1"].proteins, "sp|P1;sp|P2");
        assert!(matches!(
            PsmResultReader::new(&b""[..]),
            Err(PercolatorError::MissingColumn { .. })
        ));
    }

    #[test]
    fn score_column_is_optional_and_lenient() {
        let text = "PSMId\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    a_1_1\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\n";
        assert_eq!(parse_psm_results(text).unwrap()["a_1_1"].score, None);
        // A non-finite or unparseable score is no score; the row still parses.
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    a_1_1\tnan\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\n\
                    b_1_1\tinf\t0.002\t0.0005\tR.PEPTIDEK.A\tsp|P1\n\
                    c_1_1\tx\t0.003\t0.0005\tR.PEPTIDEK.A\tsp|P1\n";
        let m = parse_psm_results(text).unwrap();
        assert_eq!(m.len(), 3);
        assert!(m.values().all(|p| p.score.is_none()));
        assert_eq!(m["a_1_1"].q_value, 0.001);
        // A score column after the scalar columns does not lengthen the rows
        // that must be present.
        let text = "PSMId\tq-value\tposterior_error_prob\tpeptide\tproteinIds\tscore\n\
                    a_1_1\t0.001\t0.0005\tR.PEPTIDEK.A\n";
        assert_eq!(parse_psm_results(text).unwrap()["a_1_1"].score, None);
    }

    #[test]
    fn decoy_results_sit_next_to_the_pin() {
        assert_eq!(
            decoy_results_path(Path::new("/x/run.pin")),
            Path::new("/x/run.percolator.decoy.psms.txt")
        );
    }

    #[test]
    fn parse_leading_filename_column() {
        // a comparison search engine adds a leading FileName column, shifting q-value to col 4. The
        // by-NAME lookup must still find every field correctly.
        let text =
            "PSMId\tFileName\tscore\tq-value\tposterior_error_probability\tpeptide\tproteinIds\n\
                    sc_5_1\trun.mzML\t9.1\t0.004\t0.002\tK.SAGTPEPK.L\tDECOY_x\n";
        let m = parse_psm_results(text).unwrap();
        let r = &m["sc_5_1"];
        assert_eq!(r.q_value, 0.004);
        assert_eq!(r.pep, 0.002);
        assert_eq!(r.peptide, "K.SAGTPEPK.L");
        assert_eq!(r.proteins, "DECOY_x");
    }

    #[test]
    fn clear_stale_removes_prior_outputs() {
        let dir = std::env::temp_dir().join(format!("andes_perc_stale_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = ResultPaths::new(&dir, "x");
        std::fs::write(&out.target, b"stale").unwrap();
        std::fs::write(&out.weights, b"stale").unwrap();
        assert!(out.target.exists());
        out.clear_stale().unwrap();
        assert!(!out.target.exists(), "stale target must be removed");
        assert!(!out.weights.exists(), "stale weights must be removed");
        // Idempotent: clearing again when files are absent is Ok.
        out.clear_stale().unwrap();
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn parse_errors_on_non_finite_qvalue() {
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    x_1_1\t2.0\tnan\t0.01\tK.PK.A\tsp|P1\n";
        let err = parse_psm_results(text).unwrap_err();
        assert!(
            matches!(err, PercolatorError::NonFiniteScore { column, .. } if column == "q-value"),
            "got {err:?}"
        );
    }

    #[test]
    fn parse_errors_on_unparseable_pep() {
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    x_1_1\t2.0\t0.01\tnotanumber\tK.PK.A\tsp|P1\n";
        let err = parse_psm_results(text).unwrap_err();
        assert!(
            matches!(err, PercolatorError::NonFiniteScore { .. }),
            "got {err:?}"
        );
    }

    #[test]
    fn parse_errors_on_short_row() {
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    x_1_1\t2.0\t0.01\n";
        let err = parse_psm_results(text).unwrap_err();
        assert!(
            matches!(err, PercolatorError::MalformedRow { .. }),
            "got {err:?}"
        );
    }

    #[test]
    fn parse_errors_on_duplicate_psmid() {
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    dup_1_1\t2.0\t0.01\t0.005\tK.PK.A\tsp|P1\n\
                    dup_1_1\t1.0\t0.02\t0.006\tK.QK.A\tsp|P2\n";
        let err = parse_psm_results(text).unwrap_err();
        assert!(
            matches!(err, PercolatorError::DuplicatePsmId { ref psm_id } if psm_id == "dup_1_1"),
            "got {err:?}"
        );
    }

    #[test]
    fn parse_missing_qvalue_column_errors() {
        let text = "PSMId\tscore\tpeptide\tproteinIds\nx_1_1\t2.0\tK.PK.A\tsp|P1\n";
        let err = parse_psm_results(text).unwrap_err();
        match err {
            PercolatorError::MissingColumn { column, .. } => assert_eq!(column, "q-value"),
            other => panic!("expected MissingColumn, got {other:?}"),
        }
    }

    #[test]
    fn resolve_explicit_missing_binary_errors() {
        let err = resolve_backend(
            Some(Path::new("/no/such/percolator")),
            false,
            DEFAULT_PERCOLATOR_IMAGE,
        )
        .unwrap_err();
        assert!(matches!(err, PercolatorError::BinaryNotFound(_)));
    }

    #[test]
    fn resolve_explicit_existing_binary() {
        // Any existing file is accepted as the binary (executability is checked
        // at spawn time). Use this test's own source file as a stand-in path.
        let me = std::env::current_exe().unwrap();
        let b = resolve_backend(Some(&me), false, DEFAULT_PERCOLATOR_IMAGE).unwrap();
        assert_eq!(b, PercolatorBackend::Binary(me));
    }

    #[test]
    fn native_args_shape() {
        let out = ResultPaths::new(Path::new("/tmp/d"), "x");
        let args = native_args(
            Path::new("/tmp/d/x.pin"),
            &out,
            &["--testFDR".into(), "0.05".into()],
        );
        assert_eq!(args[0], "--seed");
        assert_eq!(args[1], "42");
        assert!(args.contains(&"--only-psms".to_string()));
        // extra args precede the PIN (which is last).
        assert_eq!(args.last().unwrap(), "/tmp/d/x.pin");
        let fdr_pos = args.iter().position(|a| a == "--testFDR").unwrap();
        assert!(fdr_pos < args.len() - 1);
    }

    #[test]
    fn docker_args_mounts_and_image() {
        let out = ResultPaths::new(Path::new("/data"), "run");
        let args = docker_args(
            DEFAULT_PERCOLATOR_IMAGE,
            Path::new("/data/run.pin"),
            Path::new("/data"),
            &out,
            &[],
        )
        .unwrap();
        assert_eq!(args[0], "run");
        assert!(args.contains(&"--platform".to_string()));
        let pin_mount = format!(
            "{}:/pin:ro",
            std::path::absolute("/data").unwrap().display()
        );
        assert!(args.contains(&pin_mount));
        assert!(args.contains(&DEFAULT_PERCOLATOR_IMAGE.to_string()));
        assert_eq!(args.last().unwrap(), "/pin/run.pin");
        // result paths rewritten under /out
        assert!(args
            .iter()
            .any(|a| a == "/out/run.percolator.target.psms.txt"));
    }

    #[test]
    fn docker_args_relative_pin_mounts_absolute_dirs() {
        let cwd = std::env::current_dir().unwrap();
        for (pin, dir) in [("out.pin", cwd.clone()), ("res/out.pin", cwd.join("res"))] {
            let pin = Path::new(pin);
            let out_dir = pin.parent().unwrap();
            let out = result_paths(pin);
            let args = docker_args(DEFAULT_PERCOLATOR_IMAGE, pin, out_dir, &out, &[]).unwrap();
            let mounts: Vec<&String> = args
                .iter()
                .zip(args.iter().skip(1))
                .filter(|(flag, _)| *flag == "-v")
                .map(|(_, m)| m)
                .collect();
            assert_eq!(
                mounts,
                [
                    &format!("{}:/pin:ro", dir.display()),
                    &format!("{}:/out", dir.display())
                ]
            );
        }
    }
}
