use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::process::Command;

pub const DEFAULT_PERCOLATOR_IMAGE: &str = "quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2";

#[derive(Debug, Clone, PartialEq)]
pub struct PercolatorPsm {
    pub psm_id: String,
    pub q_value: f64,
    pub pooled_q_value: Option<f64>,
    pub pep: f64,
    pub score: Option<f64>,
    pub peptide: String,
    pub proteins: String,
}

impl PercolatorPsm {
    pub fn report_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.q_value
            .total_cmp(&other.q_value)
            .then_with(|| self.pep.total_cmp(&other.pep))
            .then_with(|| self.psm_id.cmp(&other.psm_id))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum PercolatorBackend {
    Binary(PathBuf),
    Docker { image: String },
}

impl PercolatorBackend {
    pub fn describe(&self) -> String {
        match self {
            PercolatorBackend::Binary(p) => format!("binary {}", p.display()),
            PercolatorBackend::Docker { image } => format!("docker {image}"),
        }
    }
}

#[derive(Debug)]
pub enum PercolatorError {
    NoBackend,
    BinaryNotFound(PathBuf),
    DockerUnavailable,
    Run {
        backend: String,
        status: String,
        log_tail: String,
    },
    Spawn(std::io::Error),
    Io(std::io::Error),
    MissingColumn {
        column: &'static str,
        header: String,
    },
    MalformedRow { line: String },
    NonFiniteScore { column: &'static str, value: String },
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

fn docker_available() -> bool {
    Command::new("docker")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

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

fn tdc_flag(extra_args: &[String]) -> Option<String> {
    let set = extra_args
        .iter()
        .any(|a| a == "-Y" || a == "--post-processing-tdc");
    (!set).then(|| "-Y".to_string())
}

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

const ACCEPTED_Q: f64 = 0.01;

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
    ];
    if tdc {
        args.push("--override".to_string());
    }
    args.push("--init-weights".to_string());
    args.push(weights);
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

fn host_dir(dir: &Path) -> std::io::Result<PathBuf> {
    if dir.as_os_str().is_empty() {
        std::env::current_dir()
    } else {
        std::path::absolute(dir)
    }
}

pub fn run_percolator(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
) -> Result<HashMap<String, PercolatorPsm>, PercolatorError> {
    run_percolator_with_decoys(backend, pin, extra_args).map(|(targets, _)| targets)
}

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
    let decoys = match std::fs::File::open(&out.decoy) {
        Ok(f) => decoy_ids(std::io::BufReader::new(f))?,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => HashSet::new(),
        Err(e) => return Err(e.into()),
    };
    Ok((targets, decoys))
}

pub fn run_percolator_pooled(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
) -> Result<HashMap<String, PercolatorPsm>, PercolatorError> {
    let out_dir = pin.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
    let out = result_paths(pin);
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
        if let Err(e) = repair_pep(backend, pin, &out_dir, &out, extra_args, &mut targets) {
            eprintln!("WARN: Percolator's PEP is 1 for every accepted PSM and the mix-max re-estimate failed: {e}");
        }
    }
    Ok(targets)
}

const POOLED_Q_COLUMN: &str = "pooled-q-value";
const RUN_COLUMN: &str = "FileName";
const RUNS_PER_SCAN: usize = 64;
const LABEL_COLUMN: &str = "Label";

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

struct RunPin {
    path: PathBuf,
    targets: usize,
    decoys: usize,
}

fn run_pin_path(dir: &Path, stem: &str, i: usize) -> PathBuf {
    dir.join(format!("{stem}.run{i}.pin"))
}

fn write_run_pins(
    pin: &Path,
    col: usize,
    chunk: &[String],
    first: usize,
    dir: &Path,
    stem: &str,
) -> std::io::Result<Vec<RunPin>> {
    use std::io::{BufRead, BufReader, BufWriter, Write};
    let mut runs: Vec<RunPin> = (0..chunk.len())
        .map(|i| RunPin {
            path: run_pin_path(dir, stem, first + i),
            targets: 0,
            decoys: 0,
        })
        .collect();
    let index: HashMap<&str, usize> = chunk
        .iter()
        .enumerate()
        .map(|(i, r)| (r.as_str(), i))
        .collect();
    let mut lines = BufReader::new(std::fs::File::open(pin)?).lines();
    let header = lines.next().transpose()?.unwrap_or_default();
    let label = header.split('\t').position(|c| c == LABEL_COLUMN);
    let mut writers = Vec::with_capacity(runs.len());
    for r in &runs {
        let mut w = BufWriter::new(std::fs::File::create(&r.path)?);
        writeln!(w, "{header}")?;
        writers.push(w);
    }
    for line in lines {
        let line = line?;
        if let Some(&i) = line.split('\t').nth(col).and_then(|r| index.get(r)) {
            writeln!(writers[i], "{line}")?;
            match label.and_then(|c| line.split('\t').nth(c)) {
                Some("-1") => runs[i].decoys += 1,
                _ => runs[i].targets += 1,
            }
        }
    }
    for mut w in writers {
        w.flush()?;
    }
    Ok(runs)
}

fn pin_spec_ids(path: &Path) -> std::io::Result<Vec<String>> {
    use std::io::{BufRead, BufReader};
    let mut ids = Vec::new();
    for line in BufReader::new(std::fs::File::open(path)?).lines().skip(1) {
        let line = line?;
        ids.push(line.split('\t').next().unwrap_or("").to_string());
    }
    Ok(ids)
}

fn render_q_values(
    path: &Path,
    q: &HashMap<String, f64>,
) -> Result<Option<(String, String)>, PercolatorError> {
    let text = match std::fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(e.into()),
    };
    let mut out = String::with_capacity(text.len() + text.len() / 16);
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
            }
        }
        if add {
            fields.insert(prot_col.min(fields.len()), old.unwrap_or_default());
        }
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    Ok(Some((text, out)))
}

fn rewrite_result_q_values(
    main: &ResultPaths,
    target_q: &HashMap<String, f64>,
    decoy_q: &HashMap<String, f64>,
) -> Result<(), PercolatorError> {
    let target = render_q_values(&main.target, target_q)?;
    let decoy = render_q_values(&main.decoy, decoy_q)?;
    if let Some((_, new)) = &target {
        crate::atomic::write_atomic(&main.target, |tmp| std::fs::write(tmp, new))?;
    }
    if let Some((_, new)) = &decoy {
        if let Err(e) = crate::atomic::write_atomic(&main.decoy, |tmp| std::fs::write(tmp, new)) {
            if let Some((old, _)) = &target {
                let _ = crate::atomic::write_atomic(&main.target, |tmp| std::fs::write(tmp, old));
            }
            return Err(e.into());
        }
    }
    Ok(())
}

struct RunQ {
    target: HashMap<String, f64>,
    decoy: HashMap<String, f64>,
    ids: Vec<String>,
}

fn score_run(
    backend: &PercolatorBackend,
    run_pin: &Path,
    out_dir: &Path,
    pass: &ResultPaths,
    weights: &Path,
    extra_args: &[String],
) -> Result<RunQ, PercolatorError> {
    let ids = pin_spec_ids(run_pin)?;
    pass.clear_stale()?;
    let result = run_static_pass(backend, run_pin, out_dir, pass, weights, extra_args, true)
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
    let cleared = pass.clear_stale();
    let (t, d) = result?;
    cleared?;
    Ok(RunQ {
        target: t.into_iter().map(|(id, p)| (id, p.q_value)).collect(),
        decoy: d.into_iter().map(|(id, p)| (id, p.q_value)).collect(),
        ids,
    })
}

pub fn apply_per_run_q(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
    targets: &mut HashMap<String, PercolatorPsm>,
) -> Result<Option<String>, PercolatorError> {
    let (col, runs) = match pin_runs(pin) {
        Ok(Some(x)) => x,
        Ok(None) => return Ok(None),
        Err(e) => {
            eprintln!("WARN: per-run q-values skipped, the PIN cannot be read ({e}); the q-values stay pooled");
            return Ok(None);
        }
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
    let failed = std::cell::Cell::new(0usize);
    let unscorable = std::cell::Cell::new(0usize);
    let warn_failed = |run: &str, why: &dyn std::fmt::Display| {
        eprintln!(
            "WARN: Percolator per-run q-values of run {run} failed ({why}); \
             its PSMs keep the pooled q-values"
        );
        failed.set(failed.get() + 1);
    };
    for (c, chunk) in runs.chunks(RUNS_PER_SCAN).enumerate() {
        let first = c * RUNS_PER_SCAN;
        let written = write_run_pins(pin, col, chunk, first, &out_dir, &stem);
        match written {
            Err(e) => {
                for run in chunk {
                    warn_failed(run, &e);
                }
            }
            Ok(run_pins) => {
                for (i, (run, rp)) in chunk.iter().zip(&run_pins).enumerate() {
                    let outcome = if rp.targets == 0 || rp.decoys == 0 {
                        pin_spec_ids(&rp.path)
                            .map(|ids| {
                                unscorable.set(unscorable.get() + 1);
                                RunQ {
                                    target: HashMap::new(),
                                    decoy: HashMap::new(),
                                    ids,
                                }
                            })
                            .map_err(PercolatorError::Io)
                    } else {
                        let pass = ResultPaths::new(&out_dir, &format!("{stem}.run{}", first + i));
                        score_run(
                            backend,
                            &rp.path,
                            &out_dir,
                            &pass,
                            &main.weights,
                            extra_args,
                        )
                    };
                    match outcome {
                        Ok(r) => {
                            for id in &r.ids {
                                target_q
                                    .entry(id.clone())
                                    .or_insert_with(|| r.target.get(id).copied().unwrap_or(1.0));
                                decoy_q
                                    .entry(id.clone())
                                    .or_insert_with(|| r.decoy.get(id).copied().unwrap_or(1.0));
                            }
                        }
                        Err(e) => warn_failed(run, &e),
                    }
                }
            }
        }
        for i in 0..chunk.len() {
            let _ = std::fs::remove_file(run_pin_path(&out_dir, &stem, first + i));
        }
    }
    let mut replaced = 0usize;
    let mut note = "";
    match rewrite_result_q_values(&main, &target_q, &decoy_q) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("WARN: the per-run q-values could not be written to the result files ({e}); the q-values stay pooled");
            note = ", result files left pooled";
            target_q.clear();
        }
    }
    for (id, psm) in targets.iter_mut() {
        if let Some(&q) = target_q.get(id) {
            psm.pooled_q_value.get_or_insert(psm.q_value);
            psm.q_value = q;
            replaced += 1;
        }
    }
    Ok(Some(format!(
        "Percolator: per-run q-values (pooled model applied to each of {} runs): {replaced} of {} target PSMs re-estimated{}{}{note}",
        runs.len(),
        targets.len(),
        if unscorable.get() == 0 {
            String::new()
        } else {
            format!(", {} runs without targets or decoys got q = 1", unscorable.get())
        },
        if failed.get() == 0 {
            String::new()
        } else {
            format!(", {} runs kept the pooled q-values", failed.get())
        }
    )))
}

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

pub type RowScores = (HashMap<String, f64>, HashMap<String, f64>);

pub fn score_all_rows(
    backend: &PercolatorBackend,
    pin: &Path,
    extra_args: &[String],
) -> Result<RowScores, PercolatorError> {
    let main = result_paths(pin);
    let out_dir = pin.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
    let stem = pin
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let pass = ResultPaths::new(&out_dir, &format!("{stem}.rows"));
    pass.clear_stale()?;
    let output = run_static_pass(backend, pin, &out_dir, &pass, &main.weights, extra_args, false)?;
    let scores = |path: &Path| -> Result<HashMap<String, f64>, PercolatorError> {
        let text = std::fs::read_to_string(path).map_err(PercolatorError::Io)?;
        Ok(parse_psm_results(&text)?
            .into_iter()
            .filter_map(|(id, p)| p.score.map(|s| (id, s)))
            .collect())
    };
    let result = if output.status.success() && pass.target.exists() && pass.decoy.exists() {
        scores(&pass.target).and_then(|t| Ok((t, scores(&pass.decoy)?)))
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
    pass.clear_stale()?;
    result
}

fn result_paths(pin: &Path) -> ResultPaths {
    let out_dir = pin.parent().unwrap_or_else(|| Path::new("."));
    let stem = pin
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    ResultPaths::new(out_dir, &stem)
}

pub fn decoy_results_path(pin: &Path) -> PathBuf {
    result_paths(pin).decoy
}

fn tail_lines(s: &str, n: usize) -> String {
    let lines: Vec<&str> = s.lines().collect();
    let start = lines.len().saturating_sub(n);
    lines[start..].join("\n")
}

fn parse_finite(s: &str, column: &'static str) -> Result<f64, PercolatorError> {
    match s.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(v),
        _ => Err(PercolatorError::NonFiniteScore {
            column,
            value: s.to_string(),
        }),
    }
}

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

#[derive(Debug)]
pub struct PsmResultRow<'a> {
    pub psm_id: &'a str,
    pub q_value: f64,
    pub pooled_q: Option<f64>,
    pub pep: f64,
    pub score: Option<f64>,
    pub peptide: &'a str,
    protein_fields: Vec<&'a str>,
}

impl PsmResultRow<'_> {
    pub fn proteins(&self) -> String {
        self.protein_fields.join(";")
    }
}

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
        let max_scalar = self.id.max(self.q).max(self.pep).max(self.peptide);
        if fields.len() <= max_scalar {
            return Err(PercolatorError::MalformedRow {
                line: line.to_string(),
            });
        }
        let q_value = parse_finite(fields[self.q], "q-value")?;
        let pep = parse_finite(fields[self.pep], "posterior_error_prob")?;
        let pooled_q = self
            .pooled_q
            .and_then(|c| fields.get(c))
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite());
        let score = self
            .score
            .and_then(|c| fields.get(c))
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite());
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

fn decoy_ids<R: std::io::BufRead>(reader: R) -> Result<HashSet<String>, PercolatorError> {
    let mut rows = PsmResultReader::new(reader)?;
    let mut ids = HashSet::new();
    while let Some(row) = rows.next_row()? {
        if !ids.insert(row.psm_id.to_string()) {
            return Err(PercolatorError::DuplicatePsmId {
                psm_id: row.psm_id.to_string(),
            });
        }
    }
    Ok(ids)
}

pub struct PsmResultReader<R> {
    reader: R,
    line: String,
    columns: ResultColumns,
}

impl<R: std::io::BufRead> PsmResultReader<R> {
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
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    title_100_1\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\tsp|P2\n\
                    title_101_1\t1.0\t0.02\t0.03\tK.ANOTHERK.R\tsp|P3\n";
        let m = parse_psm_results(text).unwrap();
        assert_eq!(m.len(), 2);
        let a = &m["title_100_1"];
        assert_eq!(a.q_value, 0.001);
        assert_eq!(a.pep, 0.0005);
        assert_eq!(a.peptide, "R.PEPTIDEK.A");
        assert_eq!(a.proteins, "sp|P1;sp|P2");
        assert_eq!(m["title_101_1"].proteins, "sp|P3");
        assert_eq!(a.score, Some(2.5));
        assert_eq!(m["title_101_1"].score, Some(1.0));
    }

    #[test]
    fn reader_splits_lines_as_str_lines() {
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
        let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                    a_1_1\tnan\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\n\
                    b_1_1\tinf\t0.002\t0.0005\tR.PEPTIDEK.A\tsp|P1\n\
                    c_1_1\tx\t0.003\t0.0005\tR.PEPTIDEK.A\tsp|P1\n";
        let m = parse_psm_results(text).unwrap();
        assert_eq!(m.len(), 3);
        assert!(m.values().all(|p| p.score.is_none()));
        assert_eq!(m["a_1_1"].q_value, 0.001);
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
        assert!(pep_is_degenerate(&map(&[
            (0.0002, 1.0),
            (0.004, 1.0),
            (0.2, 0.3)
        ])));
        assert!(!pep_is_degenerate(&map(&[(0.0002, 1e-20), (0.004, 1.0)])));
        assert!(!pep_is_degenerate(&map(&[(0.001, 0.0003), (0.3, 1.0)])));
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
        for user in ["-Y", "--post-processing-tdc"] {
            let args = native_args(Path::new("/tmp/d/x.pin"), &out, &[user.into()]);
            assert!(!args.iter().any(|a| a == "-Y" && user != "-Y"));
            assert_eq!(args.iter().filter(|a| *a == user).count(), 1);
        }
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

    #[test]
    fn unreadable_decoy_file_leaves_the_target_file_pooled() {
        let dir = tempfile::tempdir().unwrap();
        let paths = ResultPaths::new(dir.path(), "s");
        let target = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                      a\t1.0\t0.01\t0.001\tK.AAK.R\tP\n";
        std::fs::write(&paths.target, target).unwrap();
        std::fs::write(
            &paths.decoy,
            "PSMId\tscore\tpeptide\tproteinIds\nd\t0.1\tK.KAA.R\tXXX_P\n",
        )
        .unwrap();
        let target_q = HashMap::from([("a".to_string(), 0.5)]);
        let decoy_q = HashMap::from([("d".to_string(), 0.5)]);
        assert!(rewrite_result_q_values(&paths, &target_q, &decoy_q).is_err());
        assert_eq!(std::fs::read_to_string(&paths.target).unwrap(), target);

        std::fs::write(
            &paths.decoy,
            "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\nd\t0.1\t0.2\t0.9\tK.KAA.R\tXXX_P\n",
        )
        .unwrap();
        rewrite_result_q_values(&paths, &target_q, &decoy_q).unwrap();
        for (path, pooled) in [(&paths.target, "0.01"), (&paths.decoy, "0.2")] {
            let text = std::fs::read_to_string(path).unwrap();
            let mut lines = text.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let col = |name: &str| row[header.iter().position(|c| *c == name).unwrap()];
            assert_eq!(col("q-value").parse::<f64>().unwrap(), 0.5);
            assert_eq!(col(POOLED_Q_COLUMN), pooled);
        }
    }
}
