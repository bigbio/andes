//! Percolator runner: resolve a backend (explicit path → `$PATH` → docker),
//! invoke Percolator on a PIN, and parse the target-PSM result TSV back into a
//! `SpecId → PercolatorPsm` map so the QPX writer can join PEP/q-value onto the
//! per-PSM rows.
//!
//! andes does not reimplement Percolator; it orchestrates the real binary (or
//! the pinned biocontainers docker image) with
//! `--seed 42 --results-psms … --decoy-results-psms … --weights … --only-psms <PIN>`.
//! The result TSV is parsed by HEADER NAME (not column index) so it survives the
//! column-offset differences between engines (e.g. a Sage-style extra `FileName`
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
    /// The q-value of the pooled model (`pooled-q-value` column) when `q_value`
    /// is a per-run q-value ([`apply_per_run_q`]); `None` when it is the pooled one.
    pub pooled_q_value: Option<f64>,
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

impl PercolatorPsm {
    /// Report order of PSM rows: q-value, then PEP, both ascending, then `psm_id`.
    /// `psm_id` is unique within a run (it is the key of the result map), so the
    /// order is total and does not depend on the iteration order of that map.
    pub fn report_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.q_value
            .total_cmp(&other.q_value)
            .then_with(|| self.pep.total_cmp(&other.pep))
            .then_with(|| self.psm_id.cmp(&other.psm_id))
    }
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

/// The user's extra Percolator arguments with `--num-threads <threads>` added
/// in front, unless they already set it. Percolator trains one SVM per
/// cross-validation fold and by default uses three threads; on more cores it
/// runs about twice as fast with the same result (measured on a 1M-row PIN).
pub fn percolator_args_with_threads(extra: &[String], threads: usize) -> Vec<String> {
    if extra
        .iter()
        .any(|a| a == "--num-threads" || a.starts_with("--num-threads="))
    {
        return extra.to_vec();
    }
    let mut args = vec!["--num-threads".to_string(), threads.max(1).to_string()];
    args.extend(extra.iter().cloned());
    args
}

/// Target-decoy competition (`-Y`): Percolator keeps the best PSM of each
/// spectrum (scan and precursor mass) and estimates q-values from the winners,
/// instead of choosing a mode from the data. Every benchmark and calibration
/// figure of andes was measured this way. Added unless the user's arguments
/// already ask for it.
fn tdc_flag(extra_args: &[String]) -> Option<String> {
    let set = extra_args
        .iter()
        .any(|a| a == "-Y" || a == "--post-processing-tdc");
    (!set).then(|| "-Y".to_string())
}

/// Build the argument vector for a NATIVE Percolator invocation (binary path
/// case). Mirrors `run_percolator_docker.sh`: `--seed 42`, results/decoy/weights,
/// `--only-psms`, `-Y`, then any user `extra_args`, then the PIN path last.
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
    args.extend(tdc_flag(extra_args));
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
    args.extend(tdc_flag(extra_args));
    args.extend(extra_args.iter().cloned());
    args.push(format!("/pin/{pin_name}"));
    Ok(args)
}

/// The q-value below which a PSM counts as accepted when judging whether
/// Percolator's PEP is usable: Percolator's own selection FDR (its default).
const ACCEPTED_Q: f64 = 0.01;

/// Whether the PEP column is unusable: at least one PSM is accepted at
/// [`ACCEPTED_Q`] and every accepted PSM has PEP = 1. A posterior error
/// probability of 1 for each of a set whose false discovery rate is at most 1 %
/// is a contradiction, so the PEP estimate failed. Percolator's `-Y`
/// (target-decoy competition) PEP fit does this when a decoy scores far above
/// every other decoy and near the best targets: the fit then returns 1 for all
/// PSMs while the q-values stay valid. Whether it happens depends on the
/// random cross-validation split (the seed), not on the PIN columns.
pub fn pep_is_degenerate(targets: &HashMap<String, PercolatorPsm>) -> bool {
    let mut accepted = 0usize;
    for p in targets.values().filter(|p| p.q_value <= ACCEPTED_Q) {
        if p.pep < 1.0 {
            return false;
        }
        accepted += 1;
    }
    accepted > 0
}

/// Arguments of a static Percolator pass: the model of the main run applied as
/// it is (`--static`, weights from that run), without training. `--override`
/// keeps Percolator from replacing those weights by a default score vector
/// when they separate no PSM at q < 0.01, as in a run with few real PSMs.
/// `tdc` selects
/// target-decoy competition (`-Y`); without it the q-values and the PSM-level
/// PEP are estimated by Percolator's mix-max method on all PSMs. `pin`, the
/// two result files and `weights` are written as the backend sees them.
fn static_pass_args(
    pin: String,
    target: String,
    decoy: String,
    weights: String,
    extra_args: &[String],
    tdc: bool,
) -> Vec<String> {
    let mut args = vec![
        "--seed".to_string(),
        "42".to_string(),
        "--results-psms".to_string(),
        target,
        "--decoy-results-psms".to_string(),
        decoy,
        "--only-psms".to_string(),
        "--static".to_string(),
        "--override".to_string(),
        "--init-weights".to_string(),
        weights,
    ];
    if tdc {
        args.extend(tdc_flag(extra_args));
        args.extend(extra_args.iter().cloned());
    } else {
        args.extend(
            extra_args
                .iter()
                .filter(|a| *a != "-Y" && *a != "--post-processing-tdc")
                .cloned(),
        );
    }
    args.push(pin);
    args
}

/// Run a static pass (see [`static_pass_args`]) of `backend` on `pin` with the
/// model in `weights`, writing the results to `pass`. The PIN and the result
/// files are in `out_dir` or `pin`'s directory.
fn run_static_pass(
    backend: &PercolatorBackend,
    pin: &Path,
    out_dir: &Path,
    pass: &ResultPaths,
    weights: &Path,
    extra_args: &[String],
    tdc: bool,
) -> Result<std::process::Output, PercolatorError> {
    match backend {
        PercolatorBackend::Binary(bin) => Command::new(bin)
            .args(static_pass_args(
                pin.to_string_lossy().into_owned(),
                pass.target.to_string_lossy().into_owned(),
                pass.decoy.to_string_lossy().into_owned(),
                weights.to_string_lossy().into_owned(),
                extra_args,
                tdc,
            ))
            .output()
            .map_err(PercolatorError::Spawn),
        PercolatorBackend::Docker { image } => {
            let pin_dir = host_dir(pin.parent().unwrap_or_else(|| Path::new(".")))?;
            let in_out = |p: &Path| {
                let name = p
                    .file_name()
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_default();
                format!("/out/{name}")
            };
            let pin_name = pin
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();
            let mut args = vec![
                "run".to_string(),
                "--rm".to_string(),
                "--platform".to_string(),
                "linux/amd64".to_string(),
                "-v".to_string(),
                format!("{}:/pin:ro", pin_dir.display()),
                "-v".to_string(),
                format!("{}:/out", host_dir(out_dir)?.display()),
                image.clone(),
                "percolator".to_string(),
            ];
            args.extend(static_pass_args(
                format!("/pin/{pin_name}"),
                in_out(&pass.target),
                in_out(&pass.decoy),
                in_out(weights),
                extra_args,
                tdc,
            ));
            Command::new("docker")
                .args(args)
                .output()
                .map_err(PercolatorError::Spawn)
        }
    }
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
///
/// When the PIN holds several runs (a `FileName` column with more than one
/// value) the q-values are Percolator's per-run q-values, see
/// [`apply_per_run_q`]; the result files are rewritten with them.
pub fn run_percolator_with_decoys(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
) -> Result<(HashMap<String, PercolatorPsm>, HashSet<String>), PercolatorError> {
    let mut targets = run_percolator_pooled(backend, pin, extra_args)?;
    if let Some(summary) = apply_per_run_q(backend, pin, extra_args, &mut targets)? {
        eprintln!("{summary}");
    }
    let out = result_paths(pin);
    let decoys = match std::fs::read_to_string(&out.decoy) {
        Ok(text) => parse_psm_results(&text)?.into_keys().collect(),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => HashSet::new(),
        Err(e) => return Err(e.into()),
    };
    Ok((targets, decoys))
}

/// The model trained on all the PIN's PSMs: one Percolator run, with the
/// PEP repair of [`pep_is_degenerate`] when needed. The result files and the
/// weights stay next to the PIN; the q-values are the pooled ones.
pub fn run_percolator_pooled(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
) -> Result<HashMap<String, PercolatorPsm>, PercolatorError> {
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
    let mut targets = parse_psm_results(&text)?;
    if pep_is_degenerate(&targets) {
        // The main run's results are complete without the repair, so a failed
        // repair pass is a warning, not the loss of the search.
        if let Err(e) = repair_pep(backend, pin, &out_dir, &out, extra_args, &mut targets) {
            eprintln!("WARN: Percolator's PEP is 1 for every accepted PSM and the mix-max re-estimate failed: {e}");
        }
    }
    Ok(targets)
}

/// Result-file column that holds the pooled q-value of a PSM once the
/// `q-value` column holds the per-run one.
const POOLED_Q_COLUMN: &str = "pooled-q-value";
/// PIN column that names the run (input file) of each PSM.
const RUN_COLUMN: &str = "FileName";
/// Runs whose slices of the PIN are written and scored per scan of the PIN.
const RUNS_PER_SCAN: usize = 64;

/// The distinct values of the PIN's [`RUN_COLUMN`] in order of first
/// appearance, with the column's index; `None` when the PIN has no such
/// column.
fn pin_runs(pin: &Path) -> std::io::Result<Option<(usize, Vec<String>)>> {
    use std::io::{BufRead, BufReader};
    let mut lines = BufReader::new(std::fs::File::open(pin)?).lines();
    let Some(header) = lines.next().transpose()? else {
        return Ok(None);
    };
    let Some(col) = header.split('\t').position(|c| c == RUN_COLUMN) else {
        return Ok(None);
    };
    let mut seen = HashSet::new();
    let mut runs = Vec::new();
    for line in lines {
        let line = line?;
        if let Some(run) = line.split('\t').nth(col) {
            if !seen.contains(run) {
                seen.insert(run.to_string());
                runs.push(run.to_string());
            }
        }
    }
    Ok(Some((col, runs)))
}

/// Write the rows of `pin` of each run in `chunk` to its own PIN (the header
/// and that run's rows), `<dir>/<stem>.run<first + i>.pin`; returns the paths
/// in the order of `chunk`.
fn write_run_pins(
    pin: &Path,
    col: usize,
    chunk: &[String],
    first: usize,
    dir: &Path,
    stem: &str,
) -> std::io::Result<Vec<PathBuf>> {
    use std::io::{BufRead, BufReader, BufWriter, Write};
    let paths: Vec<PathBuf> = (0..chunk.len())
        .map(|i| dir.join(format!("{stem}.run{}.pin", first + i)))
        .collect();
    let index: HashMap<&str, usize> = chunk
        .iter()
        .enumerate()
        .map(|(i, r)| (r.as_str(), i))
        .collect();
    let mut writers = Vec::with_capacity(paths.len());
    let mut lines = BufReader::new(std::fs::File::open(pin)?).lines();
    let header = lines.next().transpose()?.unwrap_or_default();
    for p in &paths {
        let mut w = BufWriter::new(std::fs::File::create(p)?);
        writeln!(w, "{header}")?;
        writers.push(w);
    }
    for line in lines {
        let line = line?;
        if let Some(&i) = line.split('\t').nth(col).and_then(|r| index.get(r)) {
            writeln!(writers[i], "{line}")?;
        }
    }
    for mut w in writers {
        w.flush()?;
    }
    Ok(paths)
}

/// The `SpecId`s (first column) of the PSM rows of the PIN at `path`.
fn pin_spec_ids(path: &Path) -> std::io::Result<Vec<String>> {
    use std::io::{BufRead, BufReader};
    let mut ids = Vec::new();
    for line in BufReader::new(std::fs::File::open(path)?).lines().skip(1) {
        let line = line?;
        ids.push(line.split('\t').next().unwrap_or("").to_string());
    }
    Ok(ids)
}

/// Replace the `q-value` column of the Percolator result file `path` by the
/// values of `q` (keyed by `PSMId`); rows without a value keep theirs. The
/// replaced values stay in a `pooled-q-value` column before `proteinIds`
/// (kept as it is when the file already has one). Returns the number of rows
/// replaced.
fn rewrite_q_values(path: &Path, q: &HashMap<String, f64>) -> Result<usize, PercolatorError> {
    let text = match std::fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(0),
        Err(e) => return Err(e.into()),
    };
    let mut out = String::with_capacity(text.len() + text.len() / 16);
    let mut replaced = 0usize;
    // PSMId, q-value and proteinIds columns, and whether the pooled column is added.
    let mut cols: Option<(usize, usize, usize, bool)> = None;
    for (n, line) in text.lines().enumerate() {
        if n == 0 {
            let mut header: Vec<&str> = line.split('\t').collect();
            let find = |name: &'static str| {
                header.iter().position(|c| *c == name).ok_or_else(|| {
                    PercolatorError::MissingColumn {
                        column: name,
                        header: line.to_string(),
                    }
                })
            };
            let (id_col, q_col, prot_col) = (find("PSMId")?, find("q-value")?, find("proteinIds")?);
            let add = !header.contains(&POOLED_Q_COLUMN);
            if add {
                header.insert(prot_col, POOLED_Q_COLUMN);
            }
            cols = Some((id_col, q_col, prot_col, add));
            out.push_str(&header.join("\t"));
            out.push('\n');
            continue;
        }
        let (id_col, q_col, prot_col, add) = cols.expect("the header is the first line");
        let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        let old = fields.get(q_col).cloned();
        let new = fields
            .get(id_col)
            .and_then(|id| q.get(id))
            .map(|v| format!("{v:e}"));
        if let Some(v) = new {
            if let Some(f) = fields.get_mut(q_col) {
                *f = v;
                replaced += 1;
            }
        }
        if add {
            fields.insert(prot_col.min(fields.len()), old.unwrap_or_default());
        }
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    std::fs::write(path, out)?;
    Ok(replaced)
}

/// Per-run q-values for a PIN of several runs (a [`RUN_COLUMN`] with more
/// than one value): the model of the pooled run is applied to each run's PSMs
/// on its own (`--static`, the pooled weights, `-Y`) and Percolator estimates
/// the q-values of that run's target-decoy competition. The scores and the
/// PEP stay those of the pooled run; the `q-value` of `targets` and of both
/// result files is replaced. A PSM that the pass of its run leaves out of both
/// result files gets a q-value of 1. A run whose pass fails keeps the pooled
/// q-values of its PSMs, with a warning. The run PINs are scratch, removed after use.
///
/// Returns a one-line summary, or `None` when the PIN has one run or no
/// [`RUN_COLUMN`] (nothing changes).
pub fn apply_per_run_q(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
    targets: &mut HashMap<String, PercolatorPsm>,
) -> Result<Option<String>, PercolatorError> {
    let Some((col, runs)) = pin_runs(pin)? else {
        return Ok(None);
    };
    if runs.len() < 2 {
        return Ok(None);
    }
    let out_dir = pin.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
    let main = result_paths(pin);
    let stem = pin
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let mut target_q: HashMap<String, f64> = HashMap::new();
    let mut decoy_q: HashMap<String, f64> = HashMap::new();
    let mut failed: Vec<String> = Vec::new();
    // `SpecId`s of the PSMs of the runs whose pass worked.
    let mut covered: HashSet<String> = HashSet::new();
    for (c, chunk) in runs.chunks(RUNS_PER_SCAN).enumerate() {
        let first = c * RUNS_PER_SCAN;
        let paths = write_run_pins(pin, col, chunk, first, &out_dir, &stem)?;
        for (i, (run, run_pin)) in chunk.iter().zip(&paths).enumerate() {
            let pass = ResultPaths::new(&out_dir, &format!("{stem}.run{}", first + i));
            pass.clear_stale()?;
            let result = run_static_pass(
                backend,
                run_pin,
                &out_dir,
                &pass,
                &main.weights,
                extra_args,
                true,
            )
            .and_then(|output| {
                if output.status.success() && pass.target.exists() {
                    Ok(())
                } else {
                    let mut log = String::from_utf8_lossy(&output.stderr).into_owned();
                    if log.trim().is_empty() {
                        log = String::from_utf8_lossy(&output.stdout).into_owned();
                    }
                    Err(PercolatorError::Run {
                        backend: backend.describe(),
                        status: output.status.to_string(),
                        log_tail: tail_lines(&log, 5),
                    })
                }
            })
            .and_then(|()| {
                let t = parse_psm_results(&std::fs::read_to_string(&pass.target)?)?;
                let d = match std::fs::read_to_string(&pass.decoy) {
                    Ok(text) => parse_psm_results(&text)?,
                    Err(e) if e.kind() == std::io::ErrorKind::NotFound => HashMap::new(),
                    Err(e) => return Err(e.into()),
                };
                Ok((t, d))
            });
            // The run PIN and the pass results are scratch, gone either way.
            pass.clear_stale()?;
            if result.is_ok() {
                covered.extend(pin_spec_ids(run_pin)?);
            }
            std::fs::remove_file(run_pin)?;
            match result {
                Ok((t, d)) => {
                    target_q.extend(t.into_iter().map(|(id, p)| (id, p.q_value)));
                    decoy_q.extend(d.into_iter().map(|(id, p)| (id, p.q_value)));
                }
                Err(e) => {
                    eprintln!(
                        "WARN: Percolator per-run q-values of run {run} failed ({e}); \
                         its PSMs keep the pooled q-values"
                    );
                    failed.push(run.clone());
                }
            }
        }
    }
    // A PSM of a run whose pass worked that is in neither of its result
    // files has no verdict from Percolator: it is not accepted.
    for id in &covered {
        if !target_q.contains_key(id) && !decoy_q.contains_key(id) {
            target_q.insert(id.clone(), 1.0);
            decoy_q.insert(id.clone(), 1.0);
        }
    }
    let mut replaced = 0usize;
    for (id, psm) in targets.iter_mut() {
        if let Some(&q) = target_q.get(id) {
            psm.pooled_q_value.get_or_insert(psm.q_value);
            psm.q_value = q;
            replaced += 1;
        }
    }
    rewrite_q_values(&main.target, &target_q)?;
    rewrite_q_values(&main.decoy, &decoy_q)?;
    Ok(Some(format!(
        "Percolator: per-run q-values (pooled model applied to each of {} runs): {replaced} of {}          target PSMs re-estimated{}",
        runs.len(),
        targets.len(),
        if failed.is_empty() {
            String::new()
        } else {
            format!(", {} runs kept the pooled q-values", failed.len())
        }
    )))
}

/// Replaces the PEP of `targets` when Percolator's `-Y` PEP fit failed (see
/// [`pep_is_degenerate`]): a second, static pass with the weights of the main
/// run estimates the PEP by the mix-max method. The q-values, scores and the
/// set of PSMs stay those of the main run. If that PEP is degenerate too, the
/// main run's values are kept and a warning is printed.
fn repair_pep(
    backend: &PercolatorBackend,
    pin: &Path,
    out_dir: &Path,
    main: &ResultPaths,
    extra_args: &[String],
    targets: &mut HashMap<String, PercolatorPsm>,
) -> Result<(), PercolatorError> {
    let stem = pin
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let pass = ResultPaths::new(out_dir, &format!("{stem}.pep"));
    pass.clear_stale()?;
    let output = run_static_pass(
        backend,
        pin,
        out_dir,
        &pass,
        &main.weights,
        extra_args,
        false,
    )?;
    let result = if output.status.success() && pass.target.exists() {
        std::fs::read_to_string(&pass.target)
            .map_err(PercolatorError::Io)
            .and_then(|t| parse_psm_results(&t))
    } else {
        let mut log = String::from_utf8_lossy(&output.stderr).into_owned();
        if log.trim().is_empty() {
            log = String::from_utf8_lossy(&output.stdout).into_owned();
        }
        Err(PercolatorError::Run {
            backend: backend.describe(),
            status: output.status.to_string(),
            log_tail: tail_lines(&log, 20),
        })
    };
    // The repair files are scratch: gone whether or not the pass worked.
    pass.clear_stale()?;
    let repaired = result?;
    if pep_is_degenerate(&repaired) {
        eprintln!(
            "WARN: Percolator's PEP is 1 for every accepted PSM, also after the mix-max \
             re-estimate; the PEP column is not usable"
        );
        return Ok(());
    }
    let mut replaced = 0usize;
    for (id, psm) in targets.iter_mut() {
        if let Some(r) = repaired.get(id) {
            psm.pep = r.pep;
            replaced += 1;
        }
    }
    eprintln!(
        "Percolator: the target-decoy-competition PEP was 1 for every accepted PSM; PEP of {replaced} \
         PSMs re-estimated by the mix-max method on the same model (q-values unchanged)"
    );
    Ok(())
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
/// matching the header row, so a leading/extra column (e.g. Sage's `FileName`)
/// does not shift the parse. `proteinIds` is the LAST logical column; every tab
/// after it is a separate protein id, so the rest of the line is captured and
/// re-joined with `;`.
pub fn parse_psm_results(text: &str) -> Result<HashMap<String, PercolatorPsm>, PercolatorError> {
    let mut reader = PsmResultReader::new(text.as_bytes())?;
    let mut map = HashMap::new();
    while let Some(row) = reader.next_row()? {
        let psm_id = row.psm_id.to_string();
        let psm = PercolatorPsm {
            pooled_q_value: row.pooled_q,
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
    /// `pooled-q-value` column; `None` when the file has none.
    pub pooled_q: Option<f64>,
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
    pooled_q: Option<usize>,
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
            pooled_q: cols.iter().position(|c| *c == POOLED_Q_COLUMN),
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
        let pooled_q = self
            .pooled_q
            .and_then(|c| fields.get(c))
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite());
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
            pooled_q,
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

    fn psm(id: &str, q: f64, pep: f64) -> PercolatorPsm {
        PercolatorPsm {
            pooled_q_value: None,
            psm_id: id.to_string(),
            q_value: q,
            pep,
            score: None,
            peptide: "K.AAK.R".to_string(),
            proteins: "P".to_string(),
        }
    }

    #[test]
    fn report_order_is_q_then_pep_then_id_whatever_the_map_order() {
        let rows: Vec<PercolatorPsm> = (0..300)
            .map(|i| {
                psm(
                    &format!("run/scan={i}"),
                    (i % 3) as f64 * 0.001,
                    (i % 5) as f64 * 0.01,
                )
            })
            .collect();
        let sorted_ids = |map: &HashMap<String, PercolatorPsm>| -> Vec<String> {
            let mut v: Vec<&PercolatorPsm> = map.values().collect();
            v.sort_by(|a, b| a.report_cmp(b));
            v.iter().map(|p| p.psm_id.clone()).collect()
        };
        let first: HashMap<String, PercolatorPsm> =
            rows.iter().map(|p| (p.psm_id.clone(), p.clone())).collect();
        let want = sorted_ids(&first);
        // Each map has its own random hasher, so each iterates in its own order.
        for _ in 0..8 {
            let m: HashMap<String, PercolatorPsm> = rows
                .iter()
                .rev()
                .map(|p| (p.psm_id.clone(), p.clone()))
                .collect();
            assert_eq!(sorted_ids(&m), want);
        }
        let by_key: Vec<_> = want.iter().map(|id| &first[id]).collect();
        assert!(by_key.windows(2).all(|w| {
            (w[0].q_value, w[0].pep, w[0].psm_id.as_str())
                < (w[1].q_value, w[1].pep, w[1].psm_id.as_str())
        }));
    }

    #[test]
    fn percolator_gets_the_thread_count_unless_the_user_set_it() {
        let user = vec!["--testFDR".to_string(), "0.05".to_string()];
        assert_eq!(
            percolator_args_with_threads(&user, 16),
            ["--num-threads", "16", "--testFDR", "0.05"]
        );
        let set = vec!["--num-threads".to_string(), "4".to_string()];
        assert_eq!(percolator_args_with_threads(&set, 16), set);
        let set_eq = vec!["--num-threads=4".to_string()];
        assert_eq!(percolator_args_with_threads(&set_eq, 16), set_eq);
        assert_eq!(percolator_args_with_threads(&[], 0), ["--num-threads", "1"]);
    }

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
        assert_eq!(
            (a.psm_id, a.score, a.peptide),
            ("a_1_1", Some(2.5), "R.PEPTIDEK.A")
        );
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
    fn pep_of_one_for_every_accepted_psm_is_degenerate() {
        let map = |rows: &[(f64, f64)]| -> HashMap<String, PercolatorPsm> {
            rows.iter()
                .enumerate()
                .map(|(i, (q, pep))| (format!("r/{i}"), psm(&format!("r/{i}"), *q, *pep)))
                .collect()
        };
        // All accepted PSMs at PEP 1: the fit failed.
        assert!(pep_is_degenerate(&map(&[
            (0.0002, 1.0),
            (0.004, 1.0),
            (0.2, 0.3)
        ])));
        // One accepted PSM below 1 is a working fit.
        assert!(!pep_is_degenerate(&map(&[(0.0002, 1e-20), (0.004, 1.0)])));
        // PEP 1 only outside the accepted set is the normal tail.
        assert!(!pep_is_degenerate(&map(&[(0.001, 0.0003), (0.3, 1.0)])));
        // Nothing accepted: nothing to judge.
        assert!(!pep_is_degenerate(&map(&[(0.3, 1.0)])));
        assert!(!pep_is_degenerate(&HashMap::new()));
    }

    #[test]
    fn pep_pass_is_static_mix_max_without_tdc() {
        let args = static_pass_args(
            "/p/x.pin".into(),
            "/p/t".into(),
            "/p/d".into(),
            "/p/w".into(),
            &[
                "--num-threads".into(),
                "4".into(),
                "-Y".into(),
                "--post-processing-tdc".into(),
            ],
            false,
        );
        assert!(args.windows(2).any(|w| w == ["--init-weights", "/p/w"]));
        assert!(args.contains(&"--static".to_string()));
        assert!(!args
            .iter()
            .any(|a| a == "-Y" || a == "--post-processing-tdc"));
        assert!(args.windows(2).any(|w| w == ["--num-threads", "4"]));
        assert_eq!(args.last().unwrap(), "/p/x.pin");
    }

    #[test]
    fn per_run_pass_is_static_with_tdc_and_keeps_the_weights() {
        let args = static_pass_args(
            "/p/x.pin".into(),
            "/p/t".into(),
            "/p/d".into(),
            "/p/w".into(),
            &["-Y".into(), "--num-threads".into(), "4".into()],
            true,
        );
        assert_eq!(args.iter().filter(|a| *a == "-Y").count(), 1);
        for flag in ["--static", "--override"] {
            assert!(args.contains(&flag.to_string()), "{flag}");
        }
        assert!(args.windows(2).any(|w| w == ["--init-weights", "/p/w"]));
        assert_eq!(args.last().unwrap(), "/p/x.pin");
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
        // Sage adds a leading FileName column, shifting q-value to col 4. The
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
        assert_eq!(args.iter().filter(|a| *a == "-Y").count(), 1);
        // A user -Y (either spelling) is not repeated.
        for user in ["-Y", "--post-processing-tdc"] {
            let args = native_args(Path::new("/tmp/d/x.pin"), &out, &[user.into()]);
            assert!(!args.iter().any(|a| a == "-Y" && user != "-Y"));
            assert_eq!(args.iter().filter(|a| *a == user).count(), 1);
        }
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
        assert!(args.contains(&"-Y".to_string()));
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
