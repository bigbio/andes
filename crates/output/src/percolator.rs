pub(crate) mod enzyme;

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
    Native { memory: rescore::MemLimit },
}

impl PercolatorBackend {
    pub fn describe(&self) -> String {
        match self {
            PercolatorBackend::Binary(p) => format!("binary {}", p.display()),
            PercolatorBackend::Docker { image } => format!("docker {image}"),
            PercolatorBackend::Native { memory } => match memory.sizing {
                rescore::Sizing::Fixed => format!(
                    "native ({:.1} GB in memory, spill to {})",
                    memory.bytes as f64 / (1u64 << 30) as f64,
                    memory.spill_dir.display()
                ),
                rescore::Sizing::FitInput { cap, .. } => format!(
                    "native (memory sized to the PIN, up to {:.1} GB, spill to {})",
                    memory.bytes.min(cap) as f64 / (1u64 << 30) as f64,
                    memory.spill_dir.display()
                ),
            },
        }
    }
}

struct Invocation {
    success: bool,
    status: String,
    log: String,
}

impl Invocation {
    fn of_process(output: std::process::Output) -> Self {
        let mut log = String::from_utf8_lossy(&output.stderr).into_owned();
        if log.trim().is_empty() {
            log = String::from_utf8_lossy(&output.stdout).into_owned();
        }
        Invocation {
            success: output.status.success(),
            status: output.status.to_string(),
            log,
        }
    }

    fn of_native(args: &[String], memory: &rescore::MemLimit) -> Self {
        match rescore::run(args, memory.clone()) {
            Ok(()) => Invocation {
                success: true,
                status: "finished".to_string(),
                log: String::new(),
            },
            Err(e) => Invocation {
                success: false,
                status: "error".to_string(),
                log: e,
            },
        }
    }

    fn failure(&self, backend: &PercolatorBackend, tail: usize) -> PercolatorError {
        PercolatorError::Run {
            backend: backend.describe(),
            status: self.status.clone(),
            log_tail: tail_lines(&self.log, tail),
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
    let mut args = vec![
        "--num-threads".to_string(),
        threads.clamp(1, 128).to_string(),
    ];
    args.extend(extra.iter().cloned());
    args
}

fn tdc_flag(extra_args: &[String]) -> Option<String> {
    let set = extra_args
        .iter()
        .any(|a| a == "-Y" || a == "--post-processing-tdc");
    (!set).then(|| "-Y".to_string())
}

/// `--seed 42` unless `extra_args` set the seed.
fn default_seed(extra_args: &[String]) -> Vec<String> {
    if extra_args.iter().any(|a| a == "--seed" || a == "-S") {
        Vec::new()
    } else {
        vec!["--seed".to_string(), "42".to_string()]
    }
}

fn native_args(pin: &Path, out: &ResultPaths, extra_args: &[String]) -> Vec<String> {
    let mut args = default_seed(extra_args);
    args.extend([
        "--results-psms".to_string(),
        out.target.to_string_lossy().into_owned(),
        "--decoy-results-psms".to_string(),
        out.decoy.to_string_lossy().into_owned(),
        "--weights".to_string(),
        out.weights.to_string_lossy().into_owned(),
        "--only-psms".to_string(),
    ]);
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
    ];
    args.extend(default_seed(extra_args));
    args.extend([
        "--results-psms".to_string(),
        in_out(&out.target),
        "--decoy-results-psms".to_string(),
        in_out(&out.decoy),
        "--weights".to_string(),
        in_out(&out.weights),
        "--only-psms".to_string(),
    ]);
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
    let mut args = default_seed(extra_args);
    args.extend([
        "--results-psms".to_string(),
        target,
        "--decoy-results-psms".to_string(),
        decoy,
        "--only-psms".to_string(),
        "--static".to_string(),
    ]);
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
) -> Result<Invocation, PercolatorError> {
    let local_args = || {
        static_pass_args(
            pin.to_string_lossy().into_owned(),
            pass.target.to_string_lossy().into_owned(),
            pass.decoy.to_string_lossy().into_owned(),
            weights.to_string_lossy().into_owned(),
            extra_args,
            tdc,
        )
    };
    match backend {
        PercolatorBackend::Binary(bin) => Command::new(bin)
            .args(local_args())
            .output()
            .map(Invocation::of_process)
            .map_err(PercolatorError::Spawn),
        PercolatorBackend::Native { memory } => Ok(Invocation::of_native(&local_args(), memory)),
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
                .map(Invocation::of_process)
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
        log::info!("{summary}");
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

    let ran = match backend {
        PercolatorBackend::Binary(bin) => Command::new(bin)
            .args(native_args(pin, &out, extra_args))
            .output()
            .map(Invocation::of_process)
            .map_err(PercolatorError::Spawn)?,
        PercolatorBackend::Docker { image } => Command::new("docker")
            .args(docker_args(image, pin, &out_dir, &out, extra_args)?)
            .output()
            .map(Invocation::of_process)
            .map_err(PercolatorError::Spawn)?,
        PercolatorBackend::Native { memory } => {
            Invocation::of_native(&native_args(pin, &out, extra_args), memory)
        }
    };

    if !ran.success || !out.target.exists() {
        return Err(ran.failure(backend, 20));
    }

    let text = std::fs::read_to_string(&out.target)?;
    let mut targets = parse_psm_results(&text)?;
    if pep_is_degenerate(&targets) {
        if let Err(e) = repair_pep(backend, pin, &out_dir, &out, extra_args, &mut targets) {
            log::warn!("Percolator's PEP is 1 for every accepted PSM and the mix-max re-estimate failed: {e}");
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
        .and_then(|ran| {
            if ran.success && pass.target.exists() {
                Ok(())
            } else {
                Err(ran.failure(backend, 5))
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
            log::warn!("per-run q-values skipped, the PIN cannot be read ({e}); the q-values stay pooled");
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
        log::warn!(
            "Percolator per-run q-values of run {run} failed ({why}); \
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
            log::warn!("the per-run q-values could not be written to the result files ({e}); the q-values stay pooled");
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
    let ran = run_static_pass(
        backend,
        pin,
        out_dir,
        &pass,
        &main.weights,
        extra_args,
        false,
    )?;
    let result = if ran.success && pass.target.exists() {
        std::fs::read_to_string(&pass.target)
            .map_err(PercolatorError::Io)
            .and_then(|t| parse_psm_results(&t))
    } else {
        Err(ran.failure(backend, 20))
    };
    pass.clear_stale()?;
    let repaired = result?;
    if pep_is_degenerate(&repaired) {
        log::warn!(
            "Percolator's PEP is 1 for every accepted PSM, also after the mix-max \
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
    log::info!(
        target: "percolator",
        "the target-decoy-competition PEP was 1 for every accepted PSM; PEP of {replaced} \
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
    let ran = run_static_pass(
        backend,
        pin,
        &out_dir,
        &pass,
        &main.weights,
        extra_args,
        false,
    )?;
    let scores = |path: &Path| -> Result<HashMap<String, f64>, PercolatorError> {
        let text = std::fs::read_to_string(path).map_err(PercolatorError::Io)?;
        Ok(parse_psm_results(&text)?
            .into_iter()
            .filter_map(|(id, p)| p.score.map(|s| (id, s)))
            .collect())
    };
    let result = if ran.success && pass.target.exists() && pass.decoy.exists() {
        scores(&pass.target).and_then(|t| Ok((t, scores(&pass.decoy)?)))
    } else {
        Err(ran.failure(backend, 20))
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
    /// The run (spectrum file) column a pooled result carries, `FileName` or `filename`.
    pub run: Option<&'a str>,
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
    run: Option<usize>,
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
            run: cols.iter().position(|c| c.eq_ignore_ascii_case(RUN_COLUMN)),
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
            run: self.run.and_then(|c| fields.get(c)).copied(),
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
#[path = "percolator/tests/percolator.rs"]
mod tests;
