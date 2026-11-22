//! `andes cohort`: plan the analysis of several projects' SDRFs.

use std::error::Error;
use std::path::{Path, PathBuf};

use crate::cli::{CohortCommand, CohortPlanArgs, CohortRunArgs, CohortXlFdrArgs};
use crate::cohort::{self, Plan};
use crate::config::RunConfig;
use output::stats::thousands;

pub(crate) fn run(cmd: CohortCommand) -> Result<(), Box<dyn Error>> {
    match cmd {
        CohortCommand::Plan(args) => plan(args),
        CohortCommand::XlFdr(args) => xl_fdr(args),
        CohortCommand::Run(args) => run_cohort(args),
    }
}

/// The cohort's plan, or every problem logged and an error.
fn load_plan(config: &Path) -> Result<Plan, Box<dyn Error>> {
    let cfg = RunConfig::load(config)?;
    let base = config.parent().unwrap_or_else(|| Path::new("."));
    match cohort::plan(&cfg, base) {
        Ok(p) => {
            for w in &p.warnings {
                log::warn!(target: "cohort", "{w}");
            }
            Ok(p)
        }
        Err(problems) => {
            for w in &problems.warnings {
                log::warn!(target: "cohort", "{w}");
            }
            for e in &problems.errors {
                log::error!(target: "cohort", "{e}");
            }
            let n = problems.errors.len();
            Err(format!(
                "{n} error{} in {}; no plan written",
                if n == 1 { "" } else { "s" },
                config.display()
            )
            .into())
        }
    }
}

fn write_plan(plan: &Plan, out: &Path) -> Result<(), Box<dyn Error>> {
    if let Some(dir) = out.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)
            .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    }
    let json = serde_json::to_string_pretty(plan)?;
    std::fs::write(out, json + "\n").map_err(|e| format!("cannot write {}: {e}", out.display()))?;
    log::info!(target: "cohort", "wrote {}", out.display());
    Ok(())
}

/// Give the spectrum file its run's name inside the run folder, so that the PIN `SpecId`
/// prefix is the run name. A link left by an earlier plan that points elsewhere is replaced.
fn stage(file: &Path, at: &Path) -> std::io::Result<()> {
    if at.symlink_metadata().is_ok() {
        if std::fs::read_link(at).is_ok_and(|t| t == file) {
            return Ok(());
        }
        std::fs::remove_file(at)?;
    }
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(file, at)
    }
    #[cfg(not(unix))]
    {
        std::fs::hard_link(file, at)
    }
}

/// Marker a run's folder gets once its search has finished, so that a rerun skips it; it
/// holds the run's [`fingerprint`].
const DONE: &str = "search.done";

/// A file's size and modification time (ns since the epoch).
fn identity(meta: &std::fs::Metadata) -> (u64, u128) {
    let mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_nanos());
    (meta.len(), mtime)
}

/// Every string of a config value, depth first.
fn strings<'a>(v: &'a serde_yaml::Value, out: &mut Vec<&'a str>) {
    match v {
        serde_yaml::Value::String(s) => out.push(s),
        serde_yaml::Value::Sequence(seq) => seq.iter().for_each(|x| strings(x, out)),
        serde_yaml::Value::Mapping(m) => m.values().for_each(|x| strings(x, out)),
        _ => {}
    }
}

/// What a run's search depends on: the spectrum file, the config's sections other than
/// `cohort`, the run's search flags, and every file those name (the database, a model),
/// each file by path, size and modification time. A finished search is kept only while its
/// fingerprint is unchanged, so editing the cohort section to add a project does not search
/// the other projects again, while replacing the database at the same path does.
fn fingerprint(
    r: &cohort::PlannedRun,
    args: &[String],
    config: &Path,
) -> Result<String, Box<dyn Error>> {
    let meta = std::fs::metadata(&r.file)
        .map_err(|e| format!("{}: cannot read {}: {e}", r.run, r.file.display()))?;
    let (len, mtime) = identity(&meta);
    let text = std::fs::read_to_string(config)
        .map_err(|e| format!("cannot read {}: {e}", config.display()))?;
    let mut yaml: serde_yaml::Value = serde_yaml::from_str(&text)?;
    if let Some(m) = yaml.as_mapping_mut() {
        m.remove("cohort");
    }
    let base = config.parent().unwrap_or_else(|| Path::new("."));
    let mut named = Vec::new();
    strings(&yaml, &mut named);
    named.extend(args.iter().map(String::as_str));
    let mut inputs = std::collections::BTreeSet::new();
    for s in named {
        let p = Path::new(s);
        let p = if p.is_absolute() {
            p.to_path_buf()
        } else {
            base.join(p)
        };
        if let Some(m) = std::fs::metadata(&p).ok().filter(|m| m.is_file()) {
            let (l, t) = identity(&m);
            inputs.insert(format!("input\t{}\t{l}\t{t}", p.display()));
        }
    }
    Ok(format!(
        "file\t{}\t{len}\t{mtime}\nargs\t{}\n{}\nconfig\t{}",
        r.file.display(),
        args.join("\u{1f}"),
        inputs.into_iter().collect::<Vec<_>>().join("\n"),
        serde_yaml::to_string(&yaml)?
    ))
}

/// The outputs the plan's groups read from a run's folder.
fn planned_outputs<'a>(plan: &'a Plan, r: &cohort::PlannedRun) -> Vec<&'a PathBuf> {
    plan.groups
        .iter()
        .flat_map(|g| g.inputs.iter())
        .filter(|p| p.starts_with(&r.dir))
        .collect()
}

/// Each run's project, by run name.
type Projects = std::collections::BTreeMap<String, String>;

/// The run -> project map of a plan written by `andes cohort plan`.
fn projects_of_plan(path: &Path) -> Result<Projects, Box<dyn Error>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let v: serde_json::Value =
        serde_json::from_str(&text).map_err(|e| format!("{}: not a plan: {e}", path.display()))?;
    let runs = v
        .get("runs")
        .and_then(|r| r.as_array())
        .ok_or_else(|| format!("{}: no `runs` list; not a plan", path.display()))?;
    let mut out = Projects::new();
    for r in runs {
        match (
            r.get("run").and_then(|x| x.as_str()),
            r.get("project").and_then(|x| x.as_str()),
        ) {
            (Some(run), Some(project)) => {
                out.insert(run.to_string(), project.to_string());
            }
            _ => {
                return Err(format!("{}: a run without `run` or `project`", path.display()).into())
            }
        }
    }
    Ok(out)
}

/// The protein q-value cut that picks a project's protein space.
const PROTEIN_SPACE_Q: f64 = 0.01;

/// Marks a finished linear stage; holds the fingerprint of the runs it searched.
const SPACE_DONE: &str = "protein_space.done";

/// A run's search flags, with the database set to its project's protein space when it has
/// one.
fn search_args(
    r: &cohort::PlannedRun,
    spaces: &std::collections::BTreeMap<String, PathBuf>,
) -> Vec<String> {
    let mut a = r.args.clone();
    if r.kind == cohort::Kind::Xl {
        if let Some(fasta) = spaces.get(&r.project) {
            let fasta = fasta.display().to_string();
            match a.iter().position(|x| x == "--database") {
                Some(i) if i + 1 < a.len() => a[i + 1] = fasta,
                _ => a.extend(["--database".to_string(), fasta]),
            }
        }
    }
    a
}

/// The group's run outputs by project, in plan order.
fn inputs_by_project(plan: &Plan, inputs: &[PathBuf]) -> Vec<(String, Vec<PathBuf>)> {
    let mut out: Vec<(String, Vec<PathBuf>)> = Vec::new();
    for p in inputs {
        let project = plan
            .runs
            .iter()
            .find(|r| p.starts_with(&r.dir))
            .map_or_else(String::new, |r| r.project.clone());
        match out.iter_mut().find(|(q, _)| *q == project) {
            Some((_, v)) => v.push(p.clone()),
            None => out.push((project, vec![p.clone()])),
        }
    }
    out
}

/// The database a run searches: its `--database` flag, else the config's.
fn database_of(r: &cohort::PlannedRun, cfg: &RunConfig, config: &Path) -> Option<PathBuf> {
    let flag = r
        .args
        .iter()
        .position(|a| a == "--database")
        .and_then(|i| r.args.get(i + 1))
        .map(PathBuf::from);
    let p = flag.or_else(|| cfg.io.database.clone())?;
    Some(if p.is_absolute() {
        p
    } else {
        config.parent().unwrap_or_else(|| Path::new(".")).join(p)
    })
}

/// Each cross-link project's protein space: a linear and monolink search of its runs
/// (`--xl-linear-only`, rescored together, see `space_search`), whose protein groups at protein q <=
/// `PROTEIN_SPACE_Q` become the FASTA its cross-link searches look in. A project whose runs
/// differ in fragment resolution gets one such search per resolution, and its space is the
/// union of theirs. A project whose runs have different search flags, or whose database is
/// unknown, keeps the full database.
/// The search is kept while its runs, settings and database are unchanged.
fn protein_spaces(
    plan: &Plan,
    cfg: &RunConfig,
    config: &Path,
    exe: &Path,
    threads: Option<usize>,
) -> Result<std::collections::BTreeMap<String, PathBuf>, Box<dyn Error>> {
    let mut spaces = std::collections::BTreeMap::new();
    let mut projects: Vec<&str> = plan
        .runs
        .iter()
        .filter(|r| r.kind == cohort::Kind::Xl)
        .map(|r| r.project.as_str())
        .collect();
    projects.sort_unstable();
    projects.dedup();
    for project in projects {
        let runs: Vec<&cohort::PlannedRun> = plan
            .runs
            .iter()
            .filter(|r| r.kind == cohort::Kind::Xl && r.project == project)
            .collect();
        if runs.iter().any(|r| r.args != runs[0].args) {
            log::warn!(
                target: "cohort",
                "{project}: its runs have different search flags; its cross-link searches use \
                 the full database"
            );
            continue;
        }
        let Some(db) = database_of(runs[0], cfg, config) else {
            log::warn!(
                target: "cohort",
                "{project}: no database found; its cross-link searches use the full database"
            );
            continue;
        };
        let dir = plan.output_dir.join("protein_space").join(project);
        std::fs::create_dir_all(&dir)
            .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
        let fasta = dir.join(format!("{project}.fasta"));
        let mut print = String::new();
        for r in &runs {
            print.push_str(&fingerprint(r, &r.args, config)?);
            print.push('\n');
        }
        let done = dir.join(SPACE_DONE);
        let fresh = std::fs::read_to_string(&done).is_ok_and(|old| old == print) && fasta.exists();
        if !fresh {
            let _ = std::fs::remove_file(&done);
            let start = std::time::Instant::now();
            let high: Vec<Option<bool>> = runs
                .iter()
                .map(|r| {
                    crate::spectra::detect_acquisition_metadata(&r.file)
                        .1
                        .map(|i| i.is_high_resolution())
                })
                .collect();
            let groups = resolution_groups(&high);
            let mut accessions: Vec<String> = Vec::new();
            let mut cut = 0.0_f64;
            for (res, members) in &groups {
                let stem = match (groups.len(), res) {
                    (1, _) | (_, None) => project.to_string(),
                    (_, Some(true)) => format!("{project}.high-res"),
                    (_, Some(false)) => format!("{project}.low-res"),
                };
                let group: Vec<&cohort::PlannedRun> = members.iter().map(|&i| runs[i]).collect();
                let table = space_search(project, &stem, &group, &dir, config, &db, exe, threads)?;
                let (found, c) = space_accessions(&table)?;
                cut = cut.max(c);
                for a in found {
                    if !accessions.contains(&a) {
                        accessions.push(a);
                    }
                }
            }
            let n = write_space_fasta(&db, &accessions, &fasta)?;
            log::info!(
                target: "cohort",
                "{project}: protein space of {} proteins (protein q <= {cut:.3}) [{:.0} s]",
                thousands(n),
                start.elapsed().as_secs_f64()
            );
            std::fs::write(&done, &print)?;
        }
        spaces.insert(project.to_string(), fasta);
    }
    Ok(spaces)
}

/// The runs of a project split by fragment resolution (`high[i]`: run `i` has high-resolution
/// MS2, `None` when its file does not record it), in order of first appearance, because one
/// search takes its model, fragment tolerance and score features from its first input. A run
/// of unknown resolution joins the first group of known resolution.
fn resolution_groups(high: &[Option<bool>]) -> Vec<(Option<bool>, Vec<usize>)> {
    let mut groups: Vec<(Option<bool>, Vec<usize>)> = Vec::new();
    for (i, h) in high.iter().enumerate() {
        if let Some(h) = h {
            match groups.iter_mut().find(|g| g.0 == Some(*h)) {
                Some(g) => g.1.push(i),
                None => groups.push((Some(*h), vec![i])),
            }
        }
    }
    let unknown: Vec<usize> = (0..high.len()).filter(|&i| high[i].is_none()).collect();
    match groups.first_mut() {
        Some(g) => {
            g.1.extend(unknown);
            g.1.sort_unstable();
        }
        None if !unknown.is_empty() => groups.push((None, unknown)),
        None => {}
    }
    groups
}

/// The protein groups of `runs` written to `dir`, their files named after `stem`: one linear
/// and monolink search (`--xl-linear-only`) of each run on its own, kept while that run and its
/// settings are unchanged, then one native Percolator run over their pooled PINs and protein
/// inference over its PSMs. Searching run by run bounds the memory by one run, where one
/// search of all the runs holds every spectrum's candidates until it writes its PIN. Returns
/// the `protein_groups.tsv`.
#[allow(clippy::too_many_arguments)]
fn space_search(
    project: &str,
    stem: &str,
    runs: &[&cohort::PlannedRun],
    dir: &Path,
    config: &Path,
    db: &Path,
    exe: &Path,
    threads: Option<usize>,
) -> Result<PathBuf, Box<dyn Error>> {
    let run_dir = dir.join("runs");
    std::fs::create_dir_all(&run_dir)
        .map_err(|e| format!("cannot create {}: {e}", run_dir.display()))?;
    log::info!(
        target: "cohort",
        "{project}: linear search of {} runs for its protein space ({stem})",
        runs.len()
    );
    let mut pins = Vec::with_capacity(runs.len());
    for (i, r) in runs.iter().enumerate() {
        let name = Path::new(&r.stage_as)
            .file_stem()
            .map_or_else(|| r.run.clone(), |s| s.to_string_lossy().into_owned());
        let pin = run_dir.join(format!("{name}.pin"));
        let done = run_dir.join(format!("{name}.done"));
        let print = fingerprint(r, &r.args, config)?;
        if !(std::fs::read_to_string(&done).is_ok_and(|old| old == print) && pin.exists()) {
            let _ = std::fs::remove_file(&done);
            std::fs::create_dir_all(&r.dir)
                .map_err(|e| format!("cannot create {}: {e}", r.dir.display()))?;
            let staged = r.dir.join(&r.stage_as);
            stage(&r.file, &staged).map_err(|e| {
                format!(
                    "cannot stage {} as {}: {e}",
                    r.file.display(),
                    staged.display()
                )
            })?;
            let mut cmd = std::process::Command::new(exe);
            cmd.arg("--config")
                .arg(config)
                .arg("--spectrum")
                .arg(&staged)
                .arg("--output-pin")
                .arg(&pin)
                .args(&r.args)
                .arg("--xl-linear-only");
            if let Some(t) = threads {
                cmd.arg("--threads").arg(t.to_string());
            }
            let log_path = run_dir.join(format!("{name}.andes.log"));
            let log = std::fs::File::create(&log_path)
                .map_err(|e| format!("cannot write {}: {e}", log_path.display()))?;
            cmd.stdout(log.try_clone()?).stderr(log);
            let status = cmd
                .status()
                .map_err(|e| format!("{project}: cannot start {}: {e}", exe.display()))?;
            if !status.success() {
                return Err(format!(
                    "{project}: the linear search of run {} failed ({status}); see {}",
                    r.run,
                    log_path.display()
                )
                .into());
            }
            std::fs::write(&done, &print)?;
            log::info!(
                target: "cohort",
                "{project}: protein-space run {}/{} searched",
                i + 1,
                runs.len()
            );
        }
        pins.push(pin);
    }
    let pooled = dir.join(format!("{stem}.pin"));
    let rows = pool_pins(&pins, &pooled)?;
    let spill = dir.join("rescore-spill");
    crate::cli::ensure_scratch_dir(&spill)?;
    let backend = output::PercolatorBackend::Native {
        memory: crate::memlimit::rescore_memory(None, spill),
    };
    let n_threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
    });
    let extra = output::percolator_args_with_threads(&[], n_threads);
    log::info!(
        target: "cohort",
        "{project}: rescoring {} PSM rows of {} runs ({stem})",
        thousands(rows),
        runs.len()
    );
    let (targets, _) = output::run_percolator_with_decoys(&backend, &pooled, &extra)
        .map_err(|e| format!("{project}: rescoring the protein-space PSMs: {e}"))?;
    let decoy_path = output::decoy_results_path(&pooled);
    let text = std::fs::read_to_string(&decoy_path)
        .map_err(|e| format!("cannot read {}: {e}", decoy_path.display()))?;
    let decoys =
        output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", decoy_path.display()))?;
    let proteins = protein_args(&runs[0].args, config)?;
    let inf = crate::inference_cmd::infer_from_maps(db, &proteins, targets, &decoys)?;
    let table = dir.join(format!("{stem}.protein_groups.tsv"));
    output::write_protein_groups_tsv(&table, &inf)
        .map_err(|e| format!("cannot write {}: {e}", table.display()))?;
    Ok(table)
}

/// Concatenate the PINs of single-run searches into `out`, header once. Their SpecIds are
/// already distinct (each starts with its run's file name) and their FileName column keeps
/// the runs apart for Percolator. Returns the number of PSM rows.
fn pool_pins(pins: &[PathBuf], out: &Path) -> Result<usize, Box<dyn Error>> {
    use std::io::{BufRead, BufReader, Write};
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(out).map_err(|e| format!("cannot write {}: {e}", out.display()))?,
    );
    let mut head: Option<String> = None;
    let mut rows = 0;
    for pin in pins {
        let mut lines = BufReader::new(
            std::fs::File::open(pin).map_err(|e| format!("cannot read {}: {e}", pin.display()))?,
        )
        .lines();
        let h = lines
            .next()
            .transpose()?
            .ok_or_else(|| format!("{}: empty PIN", pin.display()))?;
        match &head {
            None => {
                writeln!(w, "{h}")?;
                head = Some(h);
            }
            Some(first) if *first != h => {
                return Err(format!(
                    "{}: its PIN columns differ from {}",
                    pin.display(),
                    pins[0].display()
                )
                .into());
            }
            Some(_) => {}
        }
        for line in lines {
            let line = line?;
            if !line.is_empty() {
                writeln!(w, "{line}")?;
                rows += 1;
            }
        }
    }
    w.flush()?;
    Ok(rows)
}

/// The flags of `andes inference` (`ProteinArgs`) parsed from a search's flags and its
/// `--config`: the database digestion, modifications and decoy affixes the search used. The
/// search's other flags are not inference flags and are left out.
fn protein_args(
    search: &[String],
    config: &Path,
) -> Result<crate::cli::ProteinArgs, Box<dyn Error>> {
    use clap::{CommandFactory, FromArgMatches};
    #[derive(clap::Parser)]
    struct Proteins {
        #[command(flatten)]
        proteins: crate::cli::ProteinArgs,
    }
    let known: Vec<String> = Proteins::command()
        .get_arguments()
        .filter_map(|a| a.get_long())
        .filter(|l| *l != "config")
        .map(|l| format!("--{l}"))
        .collect();
    let mut argv = vec![
        "andes".to_string(),
        "--config".to_string(),
        config.display().to_string(),
    ];
    let mut keep = false;
    for a in search {
        let is_flag =
            a.starts_with('-') && !a[1..].starts_with(|c: char| c.is_ascii_digit() || c == '.');
        if is_flag {
            let name = a.split('=').next().unwrap_or(a);
            keep = known.iter().any(|k| k == name);
        }
        if keep {
            argv.push(a.clone());
        }
    }
    let m = Proteins::command().try_get_matches_from(argv)?;
    let mut p = Proteins::from_arg_matches(&m)?.proteins;
    crate::config::apply_protein_config(&mut p, &m)?;
    Ok(p)
}

/// The accessions of the protein groups at protein q <= `PROTEIN_SPACE_Q` in a
/// `protein_groups.tsv`, and the cut used. When no group reaches it (a sample of few
/// proteins cannot: the decoy pseudo-count keeps every q above it), the groups at the
/// lowest q any group reaches, with a warning.
fn space_accessions(groups: &Path) -> Result<(Vec<String>, f64), Box<dyn Error>> {
    let text = std::fs::read_to_string(groups)
        .map_err(|e| format!("cannot read {}: {e}", groups.display()))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |c: &str| {
        header
            .iter()
            .position(|h| *h == c)
            .ok_or_else(|| format!("{}: no `{c}` column", groups.display()))
    };
    let (acc, q) = (col("accessions")?, col("protein_q")?);
    let rows: Vec<(f64, &str)> = lines
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Some((f.get(q)?.parse::<f64>().ok()?, *f.get(acc)?))
        })
        .collect();
    let best = rows.iter().map(|r| r.0).fold(f64::INFINITY, f64::min);
    if !best.is_finite() {
        return Err(format!("{}: no protein groups", groups.display()).into());
    }
    let cut = if best <= PROTEIN_SPACE_Q {
        PROTEIN_SPACE_Q
    } else {
        log::warn!(
            target: "cohort",
            "{}: no protein group reaches protein q <= {PROTEIN_SPACE_Q}; the protein space \
             takes the groups at the lowest q reached, {best:.3}",
            groups.display()
        );
        best
    };
    let mut out: Vec<String> = Vec::new();
    for (_, a) in rows.iter().filter(|r| r.0 <= cut) {
        for x in a.split(';').map(str::trim).filter(|x| !x.is_empty()) {
            if !out.iter().any(|y| y == x) {
                out.push(x.to_string());
            }
        }
    }
    Ok((out, cut))
}

/// Write the entries of `db` whose accession (the header's first word, or its UniProt
/// accession) is in `accessions` to `out`, in database order, unless `out` already holds
/// exactly them (so that its unchanged file keeps the searches' fingerprints). Returns the
/// number of entries.
fn write_space_fasta(
    db: &Path,
    accessions: &[String],
    out: &Path,
) -> Result<usize, Box<dyn Error>> {
    let wanted: std::collections::HashSet<&str> = accessions.iter().map(String::as_str).collect();
    let text =
        std::fs::read_to_string(db).map_err(|e| format!("cannot read {}: {e}", db.display()))?;
    let mut kept = String::new();
    let mut n = 0;
    let mut keep = false;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            let word = h.split_whitespace().next().unwrap_or("");
            let uniprot = word.split('|').nth(1).unwrap_or("");
            keep = wanted.contains(word) || (!uniprot.is_empty() && wanted.contains(uniprot));
            n += usize::from(keep);
        }
        if keep {
            kept.push_str(line);
            kept.push('\n');
        }
    }
    if n == 0 {
        return Err(format!(
            "{}: none of the {} protein-space accessions is in {}",
            out.display(),
            accessions.len(),
            db.display()
        )
        .into());
    }
    if std::fs::read_to_string(out).ok().as_deref() != Some(kept.as_str()) {
        std::fs::write(out, kept).map_err(|e| format!("cannot write {}: {e}", out.display()))?;
    }
    Ok(n)
}

fn run_cohort(args: CohortRunArgs) -> Result<(), Box<dyn Error>> {
    let config = std::fs::canonicalize(&args.config)
        .map_err(|e| format!("cannot read {}: {e}", args.config.display()))?;
    let plan = load_plan(&config)?;
    report(&plan);
    write_plan(&plan, &plan.output_dir.join("plan.json"))?;
    let exe = std::env::current_exe()?;
    let cfg = RunConfig::load(&config)?;
    let spaces = protein_spaces(&plan, &cfg, &config, &exe, args.threads)?;
    let n = plan.runs.len();
    let mut skipped = 0;
    let mut redone = 0;
    for (i, r) in plan.runs.iter().enumerate() {
        let run_args = search_args(r, &spaces);
        let print = fingerprint(r, &run_args, &config)?;
        let missing = planned_outputs(&plan, r).into_iter().find(|p| !p.exists());
        match std::fs::read_to_string(r.dir.join(DONE)) {
            Ok(old) if old == print && missing.is_none() => {
                skipped += 1;
                continue;
            }
            Ok(old) if old == print => {
                redone += 1;
                log::info!(
                    target: "cohort",
                    "{}: {} is missing; searching again",
                    r.run,
                    missing.map_or(String::new(), |p| p.display().to_string())
                );
                std::fs::remove_file(r.dir.join(DONE))?;
            }
            Ok(_) => {
                redone += 1;
                log::info!(
                    target: "cohort",
                    "{}: its file, config or search flags changed; searching again",
                    r.run
                );
                std::fs::remove_file(r.dir.join(DONE))?;
            }
            Err(_) => {}
        }
        std::fs::create_dir_all(&r.dir)
            .map_err(|e| format!("cannot create {}: {e}", r.dir.display()))?;
        let staged = r.dir.join(&r.stage_as);
        stage(&r.file, &staged).map_err(|e| {
            format!(
                "cannot stage {} as {}: {e}",
                r.file.display(),
                staged.display()
            )
        })?;
        let log_path = r.dir.join("andes.log");
        let log = std::fs::File::create(&log_path)
            .map_err(|e| format!("cannot write {}: {e}", log_path.display()))?;
        let mut cmd = std::process::Command::new(&exe);
        cmd.arg("--config")
            .arg(&config)
            .arg("--spectrum")
            .arg(&staged)
            .arg("--output-pin")
            .arg(r.dir.join(format!("{}.pin", r.run)));
        if let Some(t) = args.threads {
            cmd.arg("--threads").arg(t.to_string());
        }
        cmd.args(&run_args).stdout(log.try_clone()?).stderr(log);
        let start = std::time::Instant::now();
        log::info!(target: "cohort", "search {}/{n}: {}", i + 1, r.run);
        let status = cmd
            .status()
            .map_err(|e| format!("{}: cannot start {}: {e}", r.run, exe.display()))?;
        if !status.success() {
            return Err(format!(
                "{}: the search failed ({status}); see {}",
                r.run,
                log_path.display()
            )
            .into());
        }
        std::fs::write(r.dir.join(DONE), &print)?;
        log::info!(
            target: "cohort",
            "search {}/{n}: {} done in {:.0} s",
            i + 1,
            r.run,
            start.elapsed().as_secs_f64()
        );
    }
    if skipped > 0 {
        log::info!(
            target: "cohort",
            "{} of {} runs already searched with the same settings; kept",
            thousands(skipped),
            thousands(n)
        );
    }
    if redone > 0 {
        log::info!(
            target: "cohort",
            "{} of {} runs searched again after a change",
            thousands(redone),
            thousands(n)
        );
    }
    let projects: Projects = plan
        .runs
        .iter()
        .map(|r| (r.run.clone(), r.project.clone()))
        .collect();
    let mut not_run: Vec<&str> = Vec::new();
    for g in &plan.groups {
        for stage in &g.stages {
            match *stage {
                "xl-fdr" => {
                    std::fs::create_dir_all(&g.dir)
                        .map_err(|e| format!("cannot create {}: {e}", g.dir.display()))?;
                    for (project, inputs) in inputs_by_project(&plan, &g.inputs) {
                        log::info!(
                            target: "cohort",
                            "group {}, project {project}: crosslink FDR over {} runs",
                            g.name,
                            inputs.len()
                        );
                        let db = spaces.get(&project).cloned().or_else(|| {
                            plan.runs
                                .iter()
                                .find(|r| r.project == project)
                                .and_then(|r| database_of(r, &cfg, &config))
                        });
                        xl_fdr_tables(
                            &inputs,
                            &g.dir.join(format!("{}.{project}", g.name)),
                            &projects,
                            db.as_deref(),
                        )?;
                    }
                }
                other => {
                    if !not_run.contains(&other) {
                        not_run.push(other);
                    }
                }
            }
        }
    }
    if !not_run.is_empty() {
        log::warn!(
            target: "cohort",
            "`andes cohort run` does not yet run the {} stage(s); those groups' runs are searched \
             and their outputs are in {}",
            not_run.join(", "),
            plan.output_dir.join("runs").display()
        );
    }
    Ok(())
}

/// One CSM read from a run's table.
struct CsmRow {
    run: String,
    /// The index of its table among the tables read, and its scan in that table: the key of
    /// its row in the table's sibling PIN.
    table: usize,
    scan: String,
    spectrum: String,
    class: xl::decoy::PairClass,
    link: xl::rollup::LinkClass,
    score: f64,
    score_source: &'static str,
    residue_pair: String,
    protein_pair: String,
    /// The two peptides (bare, I read as L) and the proteins each maps to, for the
    /// near-copy check of inter-links; empty when the table has no such columns.
    peptides: [Vec<u8>; 2],
    proteins: [Vec<String>; 2],
}

fn pair_class(s: &str) -> Option<xl::decoy::PairClass> {
    use xl::decoy::PairClass::*;
    [TT, TD, DD].into_iter().find(|c| c.name() == s)
}

fn link_class(s: &str) -> Option<xl::rollup::LinkClass> {
    xl::rollup::LinkClass::ALL
        .into_iter()
        .find(|c| c.name() == s)
}

/// Read one run's CSM table into `out`. Rows with the score `NA` (spectra without a CSM)
/// are counted, not read; any other score that is not a number, or a scored row without
/// its residue pair or protein pair, is an error. Returns (scored, unscored).
fn read_csms(path: &Path, table: usize, out: &mut Vec<CsmRow>) -> Result<(usize, usize), String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let run = name.strip_suffix(".xl.tsv").unwrap_or(&name).to_string();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |c: &str| {
        header.iter().position(|h| *h == c).ok_or_else(|| {
            format!(
                "{}: no `{c}` column; the table must come from an andes --xl search of this \
                 version (with the residue_pair, protein_pair and score_source columns)",
                path.display()
            )
        })
    };
    let (sp, cl, li, sc) = (col("spectrum")?, col("class")?, col("link")?, col("score")?);
    let (rp, pp) = (col("residue_pair")?, col("protein_pair")?);
    let src = col("score_source")?;
    let scan = col("scan")?;
    let opt = |c: &str| header.iter().position(|h| *h == c);
    let sides = [
        (opt("alpha_peptide"), opt("alpha_proteins")),
        (opt("beta_peptide"), opt("beta_proteins")),
    ];
    let (mut n, mut unscored) = (0, 0);
    for (k, line) in lines.enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let bad = |what: &str| format!("{}: row {}: bad {what}", path.display(), k + 2);
        if get(sc) == "NA" {
            unscored += 1;
            continue;
        }
        if get(src) != "pair" {
            return Err(format!(
                "{}: row {}: score_source `{}`; the run was rescored on its own \
                 (--xl-rescore-intra/--xl-rescore-inter), so its scores are not on the scale of \
                 the other runs. Search the runs without those flags to pool them",
                path.display(),
                k + 2,
                get(src)
            ));
        }
        let score = get(sc)
            .parse::<f64>()
            .ok()
            .filter(|s| s.is_finite())
            .ok_or_else(|| bad("score"))?;
        let pair = |i: usize, what: &str| {
            Some(get(i))
                .filter(|v| !v.is_empty() && *v != "NA")
                .map(str::to_string)
                .ok_or_else(|| bad(what))
        };
        out.push(CsmRow {
            run: run.clone(),
            table,
            scan: get(scan).to_string(),
            spectrum: get(sp).to_string(),
            class: pair_class(get(cl)).ok_or_else(|| bad("class"))?,
            link: link_class(get(li)).ok_or_else(|| bad("link"))?,
            score,
            score_source: "pair",
            residue_pair: pair(rp, "residue_pair")?,
            protein_pair: pair(pp, "protein_pair")?,
            peptides: sides
                .map(|(p, _)| p.map_or_else(Vec::new, |i| crate::xl_ambiguity::folded(get(i)))),
            proteins: sides.map(|(_, a)| {
                a.map_or_else(Vec::new, |i| {
                    get(i)
                        .split(';')
                        .filter(|x| !x.is_empty())
                        .map(str::to_string)
                        .collect()
                })
            }),
        });
        n += 1;
    }
    Ok((n, unscored))
}

fn xl_fdr(args: CohortXlFdrArgs) -> Result<(), Box<dyn Error>> {
    let projects = match &args.plan {
        Some(p) => projects_of_plan(p)?,
        None => Projects::new(),
    };
    xl_fdr_tables(
        &args.tables,
        &args.output,
        &projects,
        args.database.as_deref(),
    )
}

/// The q cuts of the per-project error report.
const PROJECT_REPORT_CUTS: [f64; 2] = [0.01, 0.05];

/// The project of a run: from the plan when one is given, else the part of the run name
/// before its first `.` (the default run name is `<project>.<file stem>`).
fn project_of<'a>(run: &'a str, projects: &'a Projects) -> &'a str {
    projects
        .get(run)
        .map(String::as_str)
        .unwrap_or_else(|| run.split_once('.').map_or(run, |(p, _)| p))
}

/// The decoy estimate of a share's false target-target CSMs, as a fraction of its TT:
/// (TD - DD) / TT, and (TD + DD) / TT for homotypic links, whose decoys are only TD and DD.
fn share_fdr(c: &xl::fdr::ClassCounts, link: xl::rollup::LinkClass) -> f64 {
    if c.tt == 0 {
        return 0.0;
    }
    let false_tt = if link == xl::rollup::LinkClass::Homotypic {
        (c.td + c.dd) as f64
    } else {
        (c.td as f64 - c.dd as f64).max(0.0)
    };
    false_tt / c.tt as f64
}

/// Every project's share of the CSMs each pooled cut of `PROJECT_REPORT_CUTS` accepts, per link
/// class, for the projects with a scored CSM of that class.
fn project_error(
    csms: &[CsmRow],
    q: &[Option<f64>],
    projects: &Projects,
) -> Vec<output::xl_projects::XlProjectRow> {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<(&str, usize), [xl::fdr::ClassCounts; PROJECT_REPORT_CUTS.len()]> =
        BTreeMap::new();
    for (c, q) in csms.iter().zip(q) {
        let li = xl::rollup::LinkClass::ALL
            .iter()
            .position(|&l| l == c.link)
            .unwrap_or(0);
        let e = counts
            .entry((project_of(&c.run, projects), li))
            .or_default();
        for (k, cut) in PROJECT_REPORT_CUTS.iter().enumerate() {
            if q.is_some_and(|v| v <= *cut) {
                e[k].add(c.class);
            }
        }
    }
    let mut rows = Vec::new();
    for ((project, li), at) in counts {
        let link = xl::rollup::LinkClass::ALL[li];
        for (c, cut) in at.iter().zip(PROJECT_REPORT_CUTS) {
            rows.push(output::xl_projects::XlProjectRow {
                project: project.to_string(),
                link: link.name().to_string(),
                q_cut: cut,
                tt: c.tt,
                td: c.td,
                dd: c.dd,
                decoy_fdr: share_fdr(c, link),
            });
        }
    }
    rows
}

fn report_project_error(rows: &[output::xl_projects::XlProjectRow]) {
    let shown: Vec<&output::xl_projects::XlProjectRow> = rows.iter().filter(|r| r.tt > 0).collect();
    if shown.is_empty() {
        return;
    }
    let width = shown
        .iter()
        .map(|r| r.project.len())
        .max()
        .unwrap_or(7)
        .max(7);
    log::info!(
        target: "cohort",
        "each project's share of the pooled list (decoy estimate = (TD - DD) / TT; \
         homotypic (TD + DD) / TT):"
    );
    log::info!(
        target: "cohort",
        "{:<width$}  {:<9}  {:>5}  {:>7}  {:>5}  {:>5}  {:>8}",
        "project",
        "link",
        "q",
        "TT",
        "TD",
        "DD",
        "estimate"
    );
    for r in shown {
        log::info!(
            target: "cohort",
            "{:<width$}  {:<9}  {:>4}%  {:>7}  {:>5}  {:>5}  {:>7.1}%",
            r.project,
            r.link,
            r.q_cut * 100.0,
            thousands(r.tt),
            thousands(r.td),
            thousands(r.dd),
            r.decoy_fdr * 100.0
        );
    }
}

/// The sibling PIN of a run's CSM table: `<base>.xl.pin` beside `<base>.xl.tsv`.
fn pin_of(table: &Path) -> PathBuf {
    let name = table
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    table.with_file_name(format!(
        "{}.xl.pin",
        name.strip_suffix(".xl.tsv").unwrap_or(&name)
    ))
}

/// Replace the pair score of the CSMs of `classes` with the score of one Percolator model
/// trained on all CSMs of `tables` (their sibling PINs pooled, the `dropped` features left
/// out), so that their q-values use the evidence of the whole group. Does nothing, with a
/// warning, when a table has no PIN beside it.
fn rescore(
    tables: &[PathBuf],
    csms: &mut [CsmRow],
    classes: &[xl::rollup::LinkClass],
    dropped: &[&str],
    what: &str,
) -> Result<(), Box<dyn Error>> {
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};
    if !csms.iter().any(|c| classes.contains(&c.link)) {
        return Ok(());
    }
    let pins: Vec<PathBuf> = tables.iter().map(|t| pin_of(t)).collect();
    if let Some(missing) = pins.iter().find(|p| !p.exists()) {
        log::warn!(
            target: "cohort",
            "{} not found; {what} keep the pair score (rescoring pools every run's PIN)",
            missing.display()
        );
        return Ok(());
    }
    let dir = tempfile::Builder::new()
        .prefix("andes-xl-rescore")
        .tempdir()?;
    let pooled = dir.path().join("rescore.pin");
    let mut out = std::io::BufWriter::new(std::fs::File::create(&pooled)?);
    let mut key_of: HashMap<String, (usize, String)> = HashMap::new();
    let mut head: Option<Vec<String>> = None;
    let mut keep: Vec<usize> = Vec::new();
    let (mut si, mut sc) = (0, 0);
    for (k, pin) in pins.iter().enumerate() {
        let mut lines = BufReader::new(
            std::fs::File::open(pin).map_err(|e| format!("cannot read {}: {e}", pin.display()))?,
        )
        .lines();
        let h: Vec<String> = lines
            .next()
            .transpose()?
            .unwrap_or_default()
            .split('\t')
            .map(str::to_string)
            .collect();
        match &head {
            None => {
                let pos = |c: &str| {
                    h.iter()
                        .position(|x| x == c)
                        .ok_or_else(|| format!("{}: no `{c}` column", pin.display()))
                };
                si = pos("SpecId")?;
                sc = pos("ScanNr")?;
                keep = (0..h.len())
                    .filter(|&i| {
                        let c = h[i].as_str();
                        !["FileName", "ExpMass", "CalcMass"].contains(&c) && !dropped.contains(&c)
                    })
                    .collect();
                let cols: Vec<&str> = keep.iter().map(|&i| h[i].as_str()).collect();
                writeln!(out, "{}", cols.join("\t"))?;
                head = Some(h);
            }
            Some(first) if *first != h => {
                return Err(format!(
                    "{}: its PIN columns differ from {}; search all runs with the same andes version",
                    pin.display(),
                    pins[0].display()
                )
                .into());
            }
            Some(_) => {}
        }
        for line in lines {
            let line = line?;
            let mut f: Vec<String> = line.split('\t').map(str::to_string).collect();
            let scan = f.get(sc).cloned().unwrap_or_default();
            let n: u64 = scan
                .parse()
                .map_err(|_| format!("{}: ScanNr `{scan}` is not a number", pin.display()))?;
            let id = format!("{k}|{}", f[si]);
            f[si] = id.clone();
            f[sc] = (k as u64 * 1_000_000_000 + n).to_string();
            key_of.insert(id, (k, scan));
            let cols: Vec<&str> = keep
                .iter()
                .filter_map(|&i| f.get(i).map(String::as_str))
                .collect();
            writeln!(out, "{}", cols.join("\t"))?;
        }
    }
    out.flush()?;
    drop(out);
    let backend = output::PercolatorBackend::Native {
        memory: crate::memlimit::rescore_memory(None, dir.path().to_path_buf()),
    };
    let t = std::time::Instant::now();
    let mut scored = output::run_percolator_pooled(&backend, &pooled, &[])
        .map_err(|e| format!("rescoring {what}: {e}"))?;
    let decoy_text = std::fs::read_to_string(output::decoy_results_path(&pooled))?;
    scored.extend(output::parse_psm_results(&decoy_text).map_err(|e| format!("{e}"))?);
    let mut score: HashMap<(usize, &str), f64> = HashMap::new();
    for (id, psm) in &scored {
        if let (Some((k, scan)), Some(s)) = (key_of.get(id), psm.score) {
            score.insert((*k, scan.as_str()), s);
        }
    }
    let mut unscored = 0usize;
    let mut n = 0usize;
    for c in csms.iter_mut().filter(|c| classes.contains(&c.link)) {
        match score.get(&(c.table, c.scan.as_str())) {
            Some(&s) => {
                c.score = s;
                c.score_source = "percolator";
                n += 1;
            }
            None => {
                c.score = f64::NEG_INFINITY;
                c.score_source = "percolator";
                unscored += 1;
            }
        }
    }
    log::info!(
        target: "cohort",
        "{what} rescored by one Percolator model over {} runs: {} CSMs ({} without a \
         score, ranked last) [{:.1}s]",
        thousands(tables.len()),
        thousands(n),
        thousands(unscored),
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Proteins a peptide's near-copies are collected for, at most; a CSM lists fewer.
const NEAR_COPY_CAP: usize = 64;

/// Move to `LinkClass::InterAmbiguous` every inter-link CSM with a peptide that has a
/// near-identical copy (`xl_ambiguity`) in a protein of `db` the CSM does not list for it.
/// Returns the number moved.
fn flag_ambiguous(csms: &mut [CsmRow], db: &Path) -> Result<usize, Box<dyn Error>> {
    use std::collections::HashMap;
    let t = std::time::Instant::now();
    let mut ids: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut peptides: Vec<Vec<u8>> = Vec::new();
    for c in csms
        .iter()
        .filter(|c| c.link == xl::rollup::LinkClass::Inter)
    {
        for p in &c.peptides {
            if !p.is_empty() && !ids.contains_key(p) {
                ids.insert(p.clone(), peptides.len());
                peptides.push(p.clone());
            }
        }
    }
    let copies = crate::xl_ambiguity::near_copies(db, &peptides, NEAR_COPY_CAP)?;
    let mut n = 0;
    for c in csms
        .iter_mut()
        .filter(|c| c.link == xl::rollup::LinkClass::Inter)
    {
        let ambiguous = (0..2).any(|side| {
            ids.get(&c.peptides[side])
                .is_some_and(|&i| copies[i].iter().any(|a| !c.proteins[side].contains(a)))
        });
        if ambiguous {
            c.link = xl::rollup::LinkClass::InterAmbiguous;
            n += 1;
        }
    }
    log::info!(
        target: "cohort",
        "{} of the inter-link CSMs have a peptide with a near-identical copy in another protein \
         of {}; their FDR is estimated apart (inter_ambiguous) [{:.1}s]",
        thousands(n),
        db.display(),
        t.elapsed().as_secs_f64()
    );
    Ok(n)
}

/// Crosslink q-values over every CSM of `tables`, written to `<output>.xl.*`. With the
/// database the runs searched, inter-links whose peptide has a near-copy in another protein
/// are estimated apart and every inter-link is rescored; without it inter-links keep the
/// pair score.
fn xl_fdr_tables(
    tables: &[PathBuf],
    output: &Path,
    projects: &Projects,
    database: Option<&Path>,
) -> Result<(), Box<dyn Error>> {
    let mut csms = Vec::new();
    for (k, t) in tables.iter().enumerate() {
        let (n, unscored) = read_csms(t, k, &mut csms)?;
        log::info!(
            target: "cohort",
            "{}: {} scored CSMs, {} spectra without a CSM",
            t.display(),
            thousands(n),
            thousands(unscored)
        );
    }
    rescore(
        tables,
        &mut csms,
        &[xl::rollup::LinkClass::Intra],
        &[],
        "intra-protein links",
    )?;
    match database {
        Some(db) => {
            flag_ambiguous(&mut csms, db)?;
            rescore(
                tables,
                &mut csms,
                &[
                    xl::rollup::LinkClass::Inter,
                    xl::rollup::LinkClass::InterAmbiguous,
                ],
                &crate::xl_fdr::GAP_FEATURES,
                "inter-links",
            )?;
        }
        None => log::warn!(
            target: "cohort",
            "no database given; inter-links keep the pair score and are not checked for \
             near-copies"
        ),
    }
    let items: Vec<xl::rollup::Item> = csms
        .iter()
        .map(|c| xl::rollup::Item {
            key: c.residue_pair.clone(),
            score: c.score,
            class: c.class,
            link: c.link,
            n_csms: 1,
        })
        .collect();
    let ppi: Vec<String> = csms.iter().map(|c| c.protein_pair.clone()).collect();
    let rows: Vec<usize> = (0..csms.len()).collect();
    let levels = crate::xl_fdr::levels_of(items, ppi, rows);
    crate::xl_fdr::report_levels(&levels);
    let with = |ext: &str| {
        let mut s = output.to_path_buf().into_os_string();
        s.push(ext);
        PathBuf::from(s)
    };
    let csm_path = with(".xl.tsv");
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(&csm_path)
            .map_err(|e| format!("cannot write {}: {e}", csm_path.display()))?,
    );
    use std::io::Write;
    writeln!(
        w,
        "run\tspectrum\tclass\tlink\tscore\tscore_source\tq_value\tresidue_pair\tprotein_pair"
    )?;
    let mut q = vec![None; csms.len()];
    for (&row, v) in levels.csm_rows.iter().zip(&levels.csm_q) {
        q[row] = *v;
    }
    for (c, q) in csms.iter().zip(&q) {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
            c.run,
            c.spectrum,
            c.class.name(),
            c.link.name(),
            c.score,
            c.score_source,
            q.map_or("NA".to_string(), |v| format!("{v:.5}")),
            c.residue_pair,
            c.protein_pair
        )?;
    }
    w.flush()?;
    crate::xl_fdr::write_levels(
        &with(".xl.residue_pairs.tsv"),
        &with(".xl.ppi.tsv"),
        &levels,
        &|row| format!("{}:{}", csms[row].run, csms[row].spectrum),
    )?;
    let shares = project_error(&csms, &q, projects);
    report_project_error(&shares);
    output::xl_projects::write_xl_projects(&with(".xl.projects.parquet"), &shares)?;
    log::info!(
        target: "cohort",
        "wrote {} ({} CSMs from {} runs), the residue-pair and protein-pair tables and the \
         per-project error table",
        csm_path.display(),
        thousands(csms.len()),
        thousands(tables.len())
    );
    Ok(())
}

fn plan(args: CohortPlanArgs) -> Result<(), Box<dyn Error>> {
    let plan = load_plan(&args.config)?;
    report(&plan);
    let out = args
        .output
        .unwrap_or_else(|| plan.output_dir.join("plan.json"));
    write_plan(&plan, &out)
}

fn report(plan: &Plan) {
    let projects: std::collections::BTreeSet<&str> =
        plan.runs.iter().map(|r| r.project.as_str()).collect();
    log::info!(
        target: "summary",
        "cohort {}: {} projects, {} runs, {} group{}",
        plan.cohort,
        projects.len(),
        plan.runs.len(),
        plan.groups.len(),
        if plan.groups.len() == 1 { "" } else { "s" }
    );
    let width = plan
        .groups
        .iter()
        .map(|g| g.name.len())
        .max()
        .unwrap_or(5)
        .max(5);
    log::info!(target: "summary", "{:<width$}  kind   runs  stages  projects", "group");
    for g in &plan.groups {
        log::info!(
            target: "summary",
            "{:<width$}  {:<5}  {:>4}  {}  {}",
            g.name,
            g.kind.name(),
            g.runs.len(),
            g.stages.join(" > "),
            g.projects.join(" ")
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csm_rows_are_read_by_header_and_unscored_rows_skipped() {
        let t = tempfile::tempdir().unwrap();
        let p = t.path().join("A.r1.xl.tsv");
        std::fs::write(
            &p,
            "spectrum\tscan\tscore_source\tclass\tlink\tq_value\tscore\tresidue_pair\tprotein_pair\n\
             A.r1/s1\t1\tpair\tTT\tinter\t0.01\t12.5\tP1:5 -- P2:9\tP1 -- P2\n\
             A.r1/s2\t2\tpair\tTD\tintra\tNA\tNA\tP1:5 -- P1:30\tP1 -- P1\n\
             A.r1/s3\t3\tpair\tDD\thomotypic\t0.5\t3.0\tP3:7 -- P3:7\tP3 -- P3\n",
        )
        .unwrap();
        let mut rows = Vec::new();
        assert_eq!(read_csms(&p, 4, &mut rows).unwrap(), (2, 1));
        assert_eq!(rows[0].run, "A.r1");
        assert_eq!((rows[0].table, rows[0].scan.as_str()), (4, "1"));
        assert_eq!(rows[1].scan, "3");
        assert_eq!(rows[0].class, xl::decoy::PairClass::TT);
        assert_eq!(rows[1].link, xl::rollup::LinkClass::Homotypic);
        let old = t.path().join("B.xl.tsv");
        std::fs::write(&old, "spectrum\tclass\tlink\tscore\n").unwrap();
        assert!(read_csms(&old, 0, &mut rows)
            .unwrap_err()
            .contains("residue_pair"));
    }

    #[test]
    fn a_malformed_score_or_missing_pair_is_an_error() {
        let t = tempfile::tempdir().unwrap();
        let head = "spectrum\tscan\tscore_source\tclass\tlink\tscore\tresidue_pair\tprotein_pair\n";
        for (name, row, what) in [
            (
                "a.xl.tsv",
                "s1\t1\tpair\tTT\tinter\t12x\tP1:5 -- P2:9\tP1 -- P2",
                "row 2: bad score",
            ),
            (
                "b.xl.tsv",
                "s1\t1\tpair\tTT\tinter\tnan\tP1:5 -- P2:9\tP1 -- P2",
                "bad score",
            ),
            (
                "c.xl.tsv",
                "s1\t1\tpair\tTT\tinter\t12.5\t\tP1 -- P2",
                "bad residue_pair",
            ),
            (
                "e.xl.tsv",
                "s1\t1\tpercolator\tTT\tinter\t2.1\tP1:5 -- P2:9\tP1 -- P2",
                "rescored on its own",
            ),
            (
                "d.xl.tsv",
                "s1\t1\tpair\tTT\tinter\t12.5\tP1:5 -- P2:9",
                "bad protein_pair",
            ),
        ] {
            let p = t.path().join(name);
            std::fs::write(&p, format!("{head}{row}\n")).unwrap();
            let e = read_csms(&p, 0, &mut Vec::new()).unwrap_err();
            assert!(e.contains(what), "{name}: {e}");
        }
    }

    #[test]
    fn a_table_pin_is_its_sibling_and_without_one_intra_links_keep_the_pair_score() {
        let t = tempfile::tempdir().unwrap();
        let table = t.path().join("A.r1.xl.tsv");
        assert_eq!(pin_of(&table), t.path().join("A.r1.xl.pin"));
        let mut rows = vec![CsmRow {
            run: "A.r1".into(),
            table: 0,
            scan: "7".into(),
            spectrum: "A.r1/s7".into(),
            class: xl::decoy::PairClass::TT,
            link: xl::rollup::LinkClass::Intra,
            score: 12.5,
            score_source: "pair",
            residue_pair: "P1:5 -- P2:9".into(),
            protein_pair: "P1 -- P2".into(),
            peptides: [Vec::new(), Vec::new()],
            proteins: [Vec::new(), Vec::new()],
        }];
        rescore(
            &[table],
            &mut rows,
            &[xl::rollup::LinkClass::Intra],
            &[],
            "intra-protein links",
        )
        .unwrap();
        assert_eq!((rows[0].score, rows[0].score_source), (12.5, "pair"));
    }

    #[test]
    fn a_project_mixing_fragment_resolutions_is_searched_once_per_resolution() {
        assert_eq!(
            resolution_groups(&[Some(false), Some(true), None, Some(false)]),
            vec![(Some(false), vec![0, 2, 3]), (Some(true), vec![1])]
        );
        assert_eq!(
            resolution_groups(&[Some(true), Some(true)]),
            vec![(Some(true), vec![0, 1])]
        );
        assert_eq!(resolution_groups(&[None, None]), vec![(None, vec![0, 1])]);
        assert!(resolution_groups(&[]).is_empty());
    }

    #[test]
    fn single_run_pins_pool_under_one_header_and_differing_columns_are_an_error() {
        let t = tempfile::tempdir().unwrap();
        let pin = |name: &str, text: &str| {
            let p = t.path().join(name);
            std::fs::write(&p, text).unwrap();
            p
        };
        let a = pin(
            "a.pin",
            "SpecId\tLabel\tFileName\na/1_1\t1\ta\na/2_1\t-1\ta\n",
        );
        let b = pin("b.pin", "SpecId\tLabel\tFileName\nb/1_1\t1\tb\n");
        let out = t.path().join("pooled.pin");
        assert_eq!(pool_pins(&[a.clone(), b], &out).unwrap(), 3);
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            "SpecId\tLabel\tFileName\na/1_1\t1\ta\na/2_1\t-1\ta\nb/1_1\t1\tb\n"
        );
        let c = pin("c.pin", "SpecId\tLabel\tScanNr\nc/1_1\t1\t1\n");
        let err = pool_pins(&[a, c], &out).unwrap_err().to_string();
        assert!(err.contains("PIN columns differ"), "{err}");
    }

    #[test]
    fn inference_takes_the_database_flags_of_a_search_and_leaves_its_others_out() {
        let t = tempfile::tempdir().unwrap();
        let config = t.path().join("run.yaml");
        std::fs::write(&config, "search: {}\n").unwrap();
        let search: Vec<String> = [
            "--enzyme",
            "trypsin,lysc",
            "--fixed-mods",
            "Carbamidomethyl (C)",
            "--xl",
            "DSSO (K, S, T, Y, Protein N-term)",
            "--variable-mods",
            "Oxidation (M)",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let p = protein_args(&search, &config).unwrap();
        assert_eq!(p.db.enzyme, "trypsin,lysc");
        assert!(p.db.fixed_mods.is_some());
        assert!(p.db.variable_mods.is_some());
    }

    #[test]
    fn the_protein_space_takes_the_groups_at_one_percent_or_the_best_q_reached() {
        let t = tempfile::tempdir().unwrap();
        let g = t.path().join("P.protein_groups.tsv");
        let head = "protein_group\tlead\taccessions\tprotein_q\n";
        std::fs::write(
            &g,
            format!("{head}1\tA\tsp|A|X;sp|B|Y\t0.002\n2\tC\tsp|C|Z\t0.01\n3\tD\tsp|D|W\t0.2\n"),
        )
        .unwrap();
        let (acc, cut) = space_accessions(&g).unwrap();
        assert_eq!(cut, PROTEIN_SPACE_Q);
        assert_eq!(acc, ["sp|A|X", "sp|B|Y", "sp|C|Z"]);
        std::fs::write(
            &g,
            format!("{head}1\tA\tsp|A|X\t0.032\n2\tC\tsp|C|Z\t0.032\n3\tD\tsp|D|W\t0.06\n"),
        )
        .unwrap();
        let (acc, cut) = space_accessions(&g).unwrap();
        assert!(
            (cut - 0.032).abs() < 1e-12,
            "a sample of few proteins cannot reach 1 %"
        );
        assert_eq!(acc, ["sp|A|X", "sp|C|Z"]);
    }

    #[test]
    fn the_space_fasta_keeps_the_chosen_entries_and_an_unchanged_file() {
        let t = tempfile::tempdir().unwrap();
        let db = t.path().join("db.fasta");
        std::fs::write(
            &db,
            ">sp|A|X one\nPEPK\nTIDE\n>sp|B|Y two\nAAAK\n>C three\nKKK\n",
        )
        .unwrap();
        let out = t.path().join("P.fasta");
        let n = write_space_fasta(&db, &["A".into(), "C".into()], &out).unwrap();
        assert_eq!(n, 2);
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            ">sp|A|X one\nPEPK\nTIDE\n>C three\nKKK\n"
        );
        let before = std::fs::metadata(&out).unwrap().modified().unwrap();
        std::thread::sleep(std::time::Duration::from_millis(20));
        write_space_fasta(&db, &["A".into(), "C".into()], &out).unwrap();
        assert_eq!(std::fs::metadata(&out).unwrap().modified().unwrap(), before);
        assert!(write_space_fasta(&db, &["Q".into()], &out).is_err());
    }

    #[test]
    fn a_crosslink_run_searches_its_project_space_in_place_of_its_database() {
        let run = |args: &[&str], kind| cohort::PlannedRun {
            run: "A.r1".into(),
            project: "A".into(),
            sdrf: PathBuf::from("a.sdrf.tsv"),
            file: PathBuf::from("r1.raw"),
            dir: PathBuf::from("runs/A.r1"),
            stage_as: "A.r1.raw".into(),
            kind,
            group: "xl".into(),
            args: args.iter().map(|s| s.to_string()).collect(),
        };
        let mut spaces = std::collections::BTreeMap::new();
        spaces.insert("A".to_string(), PathBuf::from("space/A.fasta"));
        assert_eq!(
            search_args(&run(&["--xl", "DSSO"], cohort::Kind::Xl), &spaces),
            ["--xl", "DSSO", "--database", "space/A.fasta"]
        );
        assert_eq!(
            search_args(
                &run(
                    &["--database", "db.fasta", "--xl", "DSSO"],
                    cohort::Kind::Xl
                ),
                &spaces
            ),
            ["--database", "space/A.fasta", "--xl", "DSSO"]
        );
        assert_eq!(
            search_args(&run(&["--lfq"], cohort::Kind::Lfq), &spaces),
            ["--lfq"]
        );
    }

    #[test]
    fn an_inter_link_with_a_near_copy_in_an_unlisted_protein_is_ambiguous() {
        let t = tempfile::tempdir().unwrap();
        let db = t.path().join("db.fasta");
        std::fs::write(
            &db,
            ">P1\nMKTLDEAGHSKR\n>P2\nAKTIDEAGNSKR\n>P3\nMKAAAAWWWKR\n>P4\nMKGGGGYYYKR\n",
        )
        .unwrap();
        let row = |a: &str, pa: &str, b: &str, pb: &str| CsmRow {
            run: "A.r1".into(),
            table: 0,
            scan: "1".into(),
            spectrum: "s".into(),
            class: xl::decoy::PairClass::TT,
            link: xl::rollup::LinkClass::Inter,
            score: 1.0,
            score_source: "pair",
            residue_pair: String::new(),
            protein_pair: String::new(),
            peptides: [
                crate::xl_ambiguity::folded(a),
                crate::xl_ambiguity::folded(b),
            ],
            proteins: [vec![pa.to_string()], vec![pb.to_string()]],
        };
        let mut rows = vec![
            row("K.TLDEAGHSK.R", "P1", "K.AAAAWWWK.R", "P3"),
            row("K.GGGGYYYK.R", "P4", "K.AAAAWWWK.R", "P3"),
        ];
        assert_eq!(flag_ambiguous(&mut rows, &db).unwrap(), 1);
        assert_eq!(
            rows[0].link,
            xl::rollup::LinkClass::InterAmbiguous,
            "P2 holds a one-substitution copy"
        );
        assert_eq!(rows[1].link, xl::rollup::LinkClass::Inter);
    }

    #[test]
    fn staging_replaces_a_link_to_another_file() {
        let t = tempfile::tempdir().unwrap();
        let (a, b, at) = (
            t.path().join("a.raw"),
            t.path().join("b.raw"),
            t.path().join("r.raw"),
        );
        std::fs::write(&a, "a").unwrap();
        std::fs::write(&b, "b").unwrap();
        stage(&a, &at).unwrap();
        stage(&a, &at).unwrap();
        stage(&b, &at).unwrap();
        assert_eq!(std::fs::read_to_string(&at).unwrap(), "b");
    }

    #[test]
    fn the_fingerprint_follows_the_search_settings_but_not_the_cohort_section() {
        let t = tempfile::tempdir().unwrap();
        let file = t.path().join("r1.raw");
        std::fs::write(&file, "spectra").unwrap();
        let config = t.path().join("c.yaml");
        let run = |args: &[&str]| cohort::PlannedRun {
            run: "A.r1".into(),
            project: "A".into(),
            sdrf: t.path().join("a.sdrf.tsv"),
            file: file.clone(),
            dir: t.path().join("runs/A.r1"),
            stage_as: "A.r1.raw".into(),
            kind: cohort::Kind::Xl,
            group: "xl".into(),
            args: args.iter().map(|s| s.to_string()).collect(),
        };
        std::fs::write(&config, "search: {threads: 4}\ncohort: {name: one}\n").unwrap();
        let base = fingerprint(
            &run(&["--xl", "DSSO"]),
            &["--xl".into(), "DSSO".into()],
            &config,
        )
        .unwrap();
        std::fs::write(&config, "search: {threads: 4}\ncohort: {name: two}\n").unwrap();
        assert_eq!(
            fingerprint(
                &run(&["--xl", "DSSO"]),
                &["--xl".into(), "DSSO".into()],
                &config
            )
            .unwrap(),
            base
        );
        assert_ne!(
            fingerprint(
                &run(&["--xl", "DSBU"]),
                &["--xl".into(), "DSBU".into()],
                &config
            )
            .unwrap(),
            base
        );
        std::fs::write(&config, "search: {threads: 8}\ncohort: {name: two}\n").unwrap();
        assert_ne!(
            fingerprint(
                &run(&["--xl", "DSSO"]),
                &["--xl".into(), "DSSO".into()],
                &config
            )
            .unwrap(),
            base
        );
        std::fs::write(&config, "search: {threads: 4}\ncohort: {name: two}\n").unwrap();
        std::fs::write(&file, "other spectra").unwrap();
        assert_ne!(
            fingerprint(
                &run(&["--xl", "DSSO"]),
                &["--xl".into(), "DSSO".into()],
                &config
            )
            .unwrap(),
            base
        );

        std::fs::write(t.path().join("db.fasta"), ">P1\nPEPTIDEK\n").unwrap();
        std::fs::write(&config, "io: {database: db.fasta}\ncohort: {name: two}\n").unwrap();
        let with_db = fingerprint(
            &run(&["--xl", "DSSO"]),
            &["--xl".into(), "DSSO".into()],
            &config,
        )
        .unwrap();
        assert!(with_db.contains("db.fasta"), "{with_db}");
        std::fs::write(t.path().join("db.fasta"), ">P1\nPEPTIDEKR\n").unwrap();
        assert_ne!(
            fingerprint(
                &run(&["--xl", "DSSO"]),
                &["--xl".into(), "DSSO".into()],
                &config
            )
            .unwrap(),
            with_db,
            "a database replaced at the same path"
        );
    }

    #[test]
    fn projects_come_from_the_plan_when_one_is_given() {
        let t = tempfile::tempdir().unwrap();
        let p = t.path().join("plan.json");
        std::fs::write(
            &p,
            r#"{"runs": [{"run": "A.v1.r1", "project": "A.v1"}, {"run": "r2", "project": "B"}]}"#,
        )
        .unwrap();
        let projects = projects_of_plan(&p).unwrap();
        assert_eq!(project_of("A.v1.r1", &projects), "A.v1");
        assert_eq!(project_of("r2", &projects), "B");
        assert_eq!(project_of("C.r3", &projects), "C");
        std::fs::write(&p, r#"{"runs": [{"run": "x"}]}"#).unwrap();
        assert!(projects_of_plan(&p).is_err());
    }

    #[test]
    fn each_project_share_is_counted_at_the_pooled_cuts() {
        use xl::decoy::PairClass::{DD, TD, TT};
        use xl::rollup::LinkClass::{Homotypic, Intra};
        let row = |run: &str, class, link| CsmRow {
            run: run.to_string(),
            table: 0,
            scan: String::new(),
            spectrum: String::new(),
            class,
            link,
            score: 0.0,
            score_source: "pair",
            residue_pair: String::new(),
            protein_pair: String::new(),
            peptides: [Vec::new(), Vec::new()],
            proteins: [Vec::new(), Vec::new()],
        };
        let csms = vec![
            row("A.r1", TT, Intra),
            row("A.r1", TT, Intra),
            row("A.r2", TD, Intra),
            row("A.r2", TT, Intra),
            row("B.r1", TT, Intra),
            row("B.r1", DD, Intra),
            row("B.r1", TD, Homotypic),
            row("B.r1", TT, Homotypic),
        ];
        let q = vec![
            Some(0.001),
            Some(0.02),
            Some(0.03),
            None,
            Some(0.004),
            Some(0.04),
            Some(0.01),
            Some(0.01),
        ];
        let rows = project_error(&csms, &q, &Projects::new());
        let at = |p: &str, link: &str, cut: f64| {
            rows.iter()
                .find(|r| r.project == p && r.link == link && r.q_cut == cut)
                .map(|r| (r.tt, r.td, r.dd, r.decoy_fdr))
                .unwrap()
        };
        assert_eq!(at("A", "intra", 0.01), (1, 0, 0, 0.0));
        assert_eq!(at("A", "intra", 0.05), (2, 1, 0, 0.5));
        assert_eq!(at("B", "intra", 0.05), (1, 0, 1, 0.0));
        assert_eq!(at("B", "homotypic", 0.01), (1, 1, 0, 1.0));
        assert_eq!(rows.len(), 6);
        assert_eq!(project_of("PXD1.run.with.dots", &Projects::new()), "PXD1");
    }
}
