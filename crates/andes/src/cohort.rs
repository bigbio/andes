//! Cohorts: the SDRFs of several projects analysed together.
//!
//! `andes cohort plan` reads the `cohort` section of a `--config` file and every SDRF it
//! lists, resolves each run's search settings, groups the runs that can share one FDR, and
//! writes the plan as JSON. It searches nothing; a job array or the pipeline runs the plan.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::config::RunConfig;
use crate::sdrf::{self, Quant, SdrfSettings};

/// Run name template used when the cohort section gives none.
pub const DEFAULT_RUN_NAME: &str = "{project}.{stem}";

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CohortCfg {
    pub name: Option<String>,
    /// Folder holding the spectrum files, in any layout below it; unset: next to each SDRF.
    pub data_dir: Option<PathBuf>,
    /// Folder for every run's and group's outputs; unset: `<name>/` next to the config.
    pub output_dir: Option<PathBuf>,
    pub sdrfs: Vec<CohortSdrf>,
    pub run_name: Option<String>,
    pub overrides: Vec<Override>,
    pub exclude: Vec<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CohortSdrf {
    pub project: String,
    pub sdrf: PathBuf,
}

/// Settings for one project that replace what its SDRF gives (or lacks).
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Override {
    pub project: String,
    #[serde(default)]
    pub search: OverrideSearch,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OverrideSearch {
    pub enzyme: Option<String>,
    pub precursor_tol: Option<String>,
    pub fixed_mods: Option<String>,
    pub variable_mods: Option<String>,
    pub xl: Option<String>,
}

/// The kind of data a run holds; it decides the stages its group runs after the search.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    Dda,
    Lfq,
    Tmt,
    Label,
    Xl,
}

impl Kind {
    pub fn name(&self) -> &'static str {
        match self {
            Kind::Dda => "dda",
            Kind::Lfq => "lfq",
            Kind::Tmt => "tmt",
            Kind::Label => "label",
            Kind::Xl => "xl",
        }
    }

    /// The stages a group of this kind runs once all its runs are searched.
    pub fn stages(&self) -> Vec<&'static str> {
        match self {
            Kind::Dda => vec!["percolator", "inference"],
            Kind::Lfq | Kind::Label => vec!["percolator", "inference", "quant"],
            Kind::Tmt => vec!["percolator", "inference", "tmt-quant"],
            Kind::Xl => vec!["xl-fdr"],
        }
    }
}

#[derive(Debug, Serialize)]
pub struct PlannedRun {
    /// Unique name of the run in the cohort, from the run name template.
    pub run: String,
    pub project: String,
    pub sdrf: PathBuf,
    /// The spectrum file: the SDRF's name found under `data_dir` (by its real path, symbolic
    /// links to folders resolved), or next to the SDRF.
    pub file: PathBuf,
    /// The run's output folder; its search writes `<dir>/<run>.pin` and the tables beside it.
    pub dir: PathBuf,
    /// The name to give the spectrum file when it is staged, so that its stem (the PIN
    /// `SpecId` prefix) is the run name.
    pub stage_as: String,
    pub kind: Kind,
    pub group: String,
    /// The `andes` search flags that set this run's settings; the runner adds the spectrum,
    /// `--config` and the outputs.
    pub args: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct PlannedGroup {
    pub name: String,
    pub kind: Kind,
    pub stages: Vec<&'static str>,
    pub projects: Vec<String>,
    pub runs: Vec<String>,
    /// The group's output folder; its stages write `<dir>/<group>.*`.
    pub dir: PathBuf,
    /// The run outputs the group's first stage reads, one per run.
    pub inputs: Vec<PathBuf>,
}

#[derive(Debug, Serialize)]
pub struct Plan {
    pub cohort: String,
    pub output_dir: PathBuf,
    pub runs: Vec<PlannedRun>,
    pub groups: Vec<PlannedGroup>,
    pub warnings: Vec<String>,
}

/// Every problem found while planning, so that one report lists them all.
#[derive(Debug, Default)]
pub struct Problems {
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}

/// One settings group of one project: the rows of its SDRF that share every setting.
struct Part {
    project: String,
    sdrf: PathBuf,
    settings: SdrfSettings,
}

/// The settings a run is searched with, after the overrides and the config file.
#[derive(Debug, Clone, PartialEq)]
struct Resolved {
    kind: Kind,
    enzyme: Option<String>,
    /// The enzyme the run is searched with, wherever it comes from; it names the group.
    group_enzyme: String,
    precursor_tol: Option<String>,
    fixed_mods: Option<Vec<String>>,
    variable_mods: Option<Vec<String>>,
    quant_flag: Option<(String, Option<String>)>,
    xl: Vec<String>,
}

impl Resolved {
    fn args(&self) -> Vec<String> {
        let mut a = Vec::new();
        let mut push = |flag: &str, v: &str| {
            a.push(flag.to_string());
            a.push(v.to_string());
        };
        if let Some(e) = &self.enzyme {
            push("--enzyme", e);
        }
        if let Some(t) = &self.precursor_tol {
            push("--precursor-tol", t);
        }
        if let Some(m) = &self.fixed_mods {
            push("--fixed-mods", &list(m));
        }
        if let Some(m) = &self.variable_mods {
            push("--variable-mods", &list(m));
        }
        if !self.xl.is_empty() {
            push("--xl", &self.xl.join("; "));
        }
        if let Some((flag, value)) = &self.quant_flag {
            a.push(flag.clone());
            if let Some(v) = value {
                a.push(v.clone());
            }
        }
        a
    }

    /// The key of the FDR group: runs share a group only when they share the kind of data,
    /// the enzyme and the crosslinkers.
    fn group_key(&self) -> (Kind, String, String) {
        let mut xl = self.xl.clone();
        xl.sort();
        (self.kind.clone(), self.group_enzyme.clone(), xl.join("; "))
    }
}

/// A modification list as the `--fixed-mods` / `--variable-mods` value; an empty list is the
/// empty string, which those flags read as no modifications.
fn list(m: &[String]) -> String {
    m.join(", ")
}

/// The `--xl` spec of an SDRF `comment[cross-linker]` cell: the reagent name and its sites.
fn xl_spec(cell: &str) -> Result<String, String> {
    let x = xl::spec::from_sdrf(cell, model::unimod::bundled()).map_err(|e| e.to_string())?;
    Ok(format!("{} ({})", x.name, x.describe_sites()))
}

/// Check an override's values once per project, so that a bad value is reported once.
fn check_override(project: &str, o: &OverrideSearch, problems: &mut Problems) {
    let mut check = |key: &str, v: &Option<String>, f: &dyn Fn(&str) -> Result<(), String>| {
        if let Some(v) = v {
            if let Err(e) = f(v) {
                problems.errors.push(format!(
                    "cohort.overrides: {project}: search.{key} `{v}`: {e}"
                ));
            }
        }
    };
    check("enzyme", &o.enzyme, &|v| {
        crate::cli::parse_enzyme_list(v).map(|_| ())
    });
    check("precursor_tol", &o.precursor_tol, &|v| {
        crate::parse_precursor_tol(v).map(|_| ())
    });
    let none = |v: &str| {
        if v.trim().eq_ignore_ascii_case("none") {
            String::new()
        } else {
            v.to_string()
        }
    };
    check("fixed_mods", &o.fixed_mods, &|v| {
        crate::cli::parse_fixed_mods(&none(v)).map(|_| ())
    });
    check("variable_mods", &o.variable_mods, &|v| {
        crate::cli::parse_variable_mods(&none(v)).map(|_| ())
    });
    check("xl", &o.xl, &|v| crate::cli::parse_xl(v).map(|_| ()));
}

/// An override value, else nothing when the config file sets the key for every run, else
/// the SDRF's value.
fn choose(over: Option<&String>, in_config: bool, sdrf: Option<String>) -> Option<String> {
    match over {
        Some(v) => Some(v.clone()),
        None if in_config => None,
        None => sdrf,
    }
}

fn resolve(
    part: &Part,
    ov: Option<&OverrideSearch>,
    cfg: &RunConfig,
    problems: &mut Problems,
) -> Option<Resolved> {
    let s = &part.settings;
    let p = &part.project;
    let sc = &cfg.search;
    let enzyme = choose(
        ov.and_then(|o| o.enzyme.as_ref()),
        sc.enzyme.is_some(),
        s.enzyme.clone(),
    );
    let precursor_tol = choose(
        ov.and_then(|o| o.precursor_tol.as_ref()),
        sc.precursor_tol.is_some(),
        s.precursor_tol.clone(),
    );
    let fixed = choose(
        ov.and_then(|o| o.fixed_mods.as_ref()),
        sc.fixed_mods.is_some(),
        s.fixed_mods.as_ref().map(|m| list(m)),
    );
    let variable = choose(
        ov.and_then(|o| o.variable_mods.as_ref()),
        sc.variable_mods.is_some(),
        s.variable_mods.as_ref().map(|m| list(m)),
    );
    let xl_over = ov.and_then(|o| o.xl.as_ref());
    let mut xl_sdrf = Vec::new();
    let mut ok = true;
    if xl_over.is_none() {
        for cell in &s.xl {
            match xl_spec(cell) {
                Ok(spec) => xl_sdrf.push(spec),
                Err(e) => {
                    problems.errors.push(format!(
                        "{p}: {}: `comment[cross-linker]` {cell}: {e}",
                        part.sdrf.display()
                    ));
                    ok = false;
                }
            }
        }
    }
    let xl = choose(
        xl_over,
        sc.xl.is_some(),
        (!xl_sdrf.is_empty()).then(|| xl_sdrf.join("; ")),
    );
    if precursor_tol.is_none() && sc.precursor_tol.is_none() {
        problems.warnings.push(format!(
            "{p}: {} gives no precursor tolerance; its runs use the default {} \
             (set cohort.overrides search.precursor_tol to change it)",
            part.sdrf.display(),
            crate::cli::default_value::<crate::cli::SearchArgs>("precursor_tol")
        ));
    }
    let to_list = |v: Option<String>| {
        v.map(|v| {
            if v.trim().is_empty() || v.trim().eq_ignore_ascii_case("none") {
                Vec::new()
            } else {
                v.split(", ").map(str::to_string).collect()
            }
        })
    };
    let xl: Vec<String> = xl
        .map(|v| v.split("; ").map(str::to_string).collect())
        .unwrap_or_default();
    let (kind, quant_flag) = if !xl.is_empty() || sc.xl.is_some() {
        (Kind::Xl, None)
    } else {
        let q = &cfg.quant;
        if q.tmt.is_some() || q.lfq.is_some() || q.label.is_some() {
            let kind = if q.tmt.is_some() {
                Kind::Tmt
            } else if q.label.is_some() {
                Kind::Label
            } else if q.lfq == Some(true) {
                Kind::Lfq
            } else {
                Kind::Dda
            };
            let sdrf_kind = match &s.quant {
                Some(Quant::Lfq) => Kind::Lfq,
                Some(Quant::Isobaric(_)) => Kind::Tmt,
                Some(Quant::Label(_)) => Kind::Label,
                None => Kind::Dda,
            };
            if s.quant.is_some() && sdrf_kind != kind {
                problems.warnings.push(format!(
                    "{p}: the config's quant section makes its runs {}, not the {} the label \
                     column of {} gives",
                    kind.name(),
                    sdrf_kind.name(),
                    part.sdrf.display()
                ));
            }
            (kind, None)
        } else {
            match &s.quant {
                Some(Quant::Lfq) => (Kind::Lfq, Some(("--lfq".into(), None))),
                Some(Quant::Isobaric(plex)) => {
                    (Kind::Tmt, Some(("--tmt".into(), Some(plex.clone()))))
                }
                Some(Quant::Label(spec)) => {
                    (Kind::Label, Some(("--label".into(), Some(spec.clone()))))
                }
                None => (Kind::Dda, None),
            }
        }
    };
    let group_enzyme = enzyme
        .clone()
        .or_else(|| sc.enzyme.clone())
        .unwrap_or_else(|| crate::cli::default_value::<crate::cli::DatabaseArgs>("enzyme"));
    ok.then(|| Resolved {
        kind,
        enzyme,
        group_enzyme,
        precursor_tol,
        fixed_mods: to_list(fixed),
        variable_mods: to_list(variable),
        quant_flag,
        xl,
    })
}

/// The run name for one spectrum file: the template with `{project}` and `{stem}` filled.
pub fn run_name(template: &str, project: &str, file: &Path) -> String {
    let name = file
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let (stem, _) = split_name(&name);
    template
        .replace("{project}", project)
        .replace("{stem}", stem)
}

/// A file name split into its stem and its format suffix: the longest [`CONVERTED`] suffix
/// it ends with, ignoring case (`x.d.zip` -> `x` + `.d.zip`), else its last extension.
fn split_name(name: &str) -> (&str, &str) {
    let lower = name.to_ascii_lowercase();
    let known = CONVERTED
        .iter()
        .map(|e| format!(".{}", e.to_ascii_lowercase()))
        .filter(|e| e.len() < lower.len() && lower.ends_with(e.as_str()))
        .map(|e| e.len())
        .max();
    let at = match known {
        Some(n) => name.len() - n,
        None => name.rfind('.').filter(|&i| i > 0).unwrap_or(name.len()),
    };
    name.split_at(at)
}

fn check_template(t: &str) -> Result<(), String> {
    if !t.contains("{stem}") {
        return Err(format!("run_name `{t}` must contain {{stem}}"));
    }
    let rest = t.replace("{project}", "").replace("{stem}", "");
    if rest.contains('{') || rest.contains('}') {
        return Err(format!(
            "run_name `{t}`: only {{project}} and {{stem}} can be filled in"
        ));
    }
    if t.contains('/') || t.contains('\\') {
        return Err(format!("run_name `{t}` must not contain a path separator"));
    }
    Ok(())
}

/// Whether andes reads the file as spectra: mzML, Thermo `.raw`, Bruker `.d` or a `.d`
/// archive, or MGF, the first and last optionally gzipped. Any other name would be read as
/// MGF.
fn readable(file: &Path) -> bool {
    if crate::spectra::is_d_archive(file) {
        return true;
    }
    let name = file
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    [".mzml", ".raw", ".d", ".mgf"]
        .iter()
        .any(|e| name.ends_with(e))
}

/// The spectrum files below a folder by lower-cased file name. A Bruker `.d` folder is one
/// entry; its contents are not listed. Symbolic links to folders are followed, and each
/// folder is listed once, under its real path, however many links lead to it, so a link
/// loop ends and a file reached through a link keeps the project folder it is really in.
fn index_files(root: &Path) -> Result<BTreeMap<String, Vec<PathBuf>>, String> {
    let mut out: BTreeMap<String, Vec<PathBuf>> = BTreeMap::new();
    let mut seen: BTreeSet<PathBuf> = BTreeSet::new();
    let mut todo = vec![root.to_path_buf()];
    while let Some(dir) = todo.pop() {
        let real = std::fs::canonicalize(&dir)
            .map_err(|e| format!("cohort.data_dir: cannot read {}: {e}", dir.display()))?;
        if !seen.insert(real.clone()) {
            continue;
        }
        let dir = real;
        let entries = std::fs::read_dir(&dir)
            .map_err(|e| format!("cohort.data_dir: cannot read {}: {e}", dir.display()))?;
        for e in entries {
            let path = e
                .map_err(|e| format!("cohort.data_dir: cannot read {}: {e}", dir.display()))?
                .path();
            let name = path
                .file_name()
                .map(|n| n.to_string_lossy().to_ascii_lowercase())
                .unwrap_or_default();
            if path.is_dir() && !name.ends_with(".d") {
                todo.push(path);
            } else {
                out.entry(name).or_default().push(path);
            }
        }
    }
    for v in out.values_mut() {
        v.sort();
    }
    Ok(out)
}

/// Spectrum formats a file named in an SDRF may have been converted to, tried when the
/// named file itself is not under `data_dir`.
const CONVERTED: &[&str] = &[
    "mzML", "mzML.gz", "raw", "d", "d.zip", "d.tar", "d.tar.gz", "d.tgz", "mgf", "mgf.gz",
];

/// Find the file an SDRF names under `data_dir`: the name itself, else a converted copy
/// with the same stem. When `<data_dir>/<project>/` exists only that folder is searched;
/// otherwise the whole of `data_dir` is, except the folders of the other `projects`, so a
/// project never takes another project's file of the same name.
fn locate(
    index: &BTreeMap<String, Vec<PathBuf>>,
    root: &Path,
    project: &str,
    projects: &[&str],
    name: &str,
) -> Result<PathBuf, String> {
    let lower = name.to_ascii_lowercase();
    let (stem, _) = split_name(&lower);
    let real = |d: PathBuf| std::fs::canonicalize(&d).unwrap_or(d);
    let own = root.join(project);
    let (scope, others): (PathBuf, Vec<PathBuf>) = if own.is_dir() {
        (real(own), Vec::new())
    } else {
        (
            real(root.to_path_buf()),
            projects
                .iter()
                .filter(|p| **p != project)
                .map(|p| root.join(p))
                .filter(|d| d.is_dir())
                .map(real)
                .collect(),
        )
    };
    let mut names = vec![lower.clone()];
    names.extend(
        CONVERTED
            .iter()
            .map(|e| format!("{stem}.{}", e.to_ascii_lowercase()))
            .filter(|n| *n != lower),
    );
    for n in names {
        let found: Vec<&PathBuf> = index
            .get(&n)
            .map(|v| {
                v.iter()
                    .filter(|p| p.starts_with(&scope) && !others.iter().any(|o| p.starts_with(o)))
                    .collect()
            })
            .unwrap_or_default();
        match found.len() {
            0 => continue,
            1 => return Ok(found[0].clone()),
            k => {
                return Err(format!(
                    "{name}: {k} files under {} match ({}); keep one, or put each project's files \
                     in <data_dir>/<project>/",
                    scope.display(),
                    found
                        .iter()
                        .map(|p| p.display().to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                ))
            }
        }
    }
    Err(format!("{name}: not found under {}", scope.display()))
}

fn slug(s: &str) -> String {
    let mut out = String::new();
    for c in s.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
        } else if !out.ends_with('-') {
            out.push('-');
        }
    }
    out.trim_matches('-').to_string()
}

fn group_name(key: &(Kind, String, String), taken: &BTreeSet<String>) -> String {
    let linkers: Vec<&str> = key
        .2
        .split("; ")
        .filter(|s| !s.is_empty())
        .map(|s| s.split(" (").next().unwrap_or(s))
        .collect();
    let mut base = key.0.name().to_string();
    if !linkers.is_empty() {
        base = format!("{base}-{}", slug(&linkers.join("-")));
    }
    base = format!("{base}-{}", slug(&key.1));
    let mut name = base.clone();
    let mut n = 2;
    while taken.contains(&name) {
        name = format!("{base}-{n}");
        n += 1;
    }
    name
}

/// Resolve a cohort config into a plan. `base` is the directory the config's relative
/// SDRF paths are read from. Returns the plan, or every error found.
pub fn plan(cfg: &RunConfig, base: &Path) -> Result<Plan, Problems> {
    let mut problems = Problems::default();
    let Some(cohort) = &cfg.cohort else {
        problems
            .errors
            .push("the config has no `cohort` section".into());
        return Err(problems);
    };
    if cohort.sdrfs.is_empty() {
        problems.errors.push("cohort.sdrfs lists no SDRF".into());
    }
    let template = cohort.run_name.as_deref().unwrap_or(DEFAULT_RUN_NAME);
    if let Err(e) = check_template(template) {
        problems.errors.push(format!("cohort.{e}"));
    }

    let mut projects: Vec<&str> = Vec::new();
    for s in &cohort.sdrfs {
        if projects.contains(&s.project.as_str()) {
            problems.errors.push(format!(
                "cohort.sdrfs: project {} is listed twice",
                s.project
            ));
        }
        projects.push(&s.project);
    }
    for p in &cohort.exclude {
        if !projects.contains(&p.as_str()) {
            problems
                .errors
                .push(format!("cohort.exclude: {p} is not in cohort.sdrfs"));
        }
    }
    let mut overridden: Vec<&str> = Vec::new();
    for o in &cohort.overrides {
        if !projects.contains(&o.project.as_str()) {
            problems.errors.push(format!(
                "cohort.overrides: {} is not in cohort.sdrfs",
                o.project
            ));
        }
        if overridden.contains(&o.project.as_str()) {
            problems
                .errors
                .push(format!("cohort.overrides: {} has two entries", o.project));
        }
        overridden.push(&o.project);
        check_override(&o.project, &o.search, &mut problems);
    }

    let abs = |p: &Path| {
        if p.is_absolute() {
            p.to_path_buf()
        } else {
            base.join(p)
        }
    };
    let data_dir = cohort.data_dir.as_deref().map(abs);
    let index = match &data_dir {
        Some(d) => match index_files(d) {
            Ok(i) => Some(i),
            Err(e) => {
                problems.errors.push(e);
                None
            }
        },
        None => None,
    };
    let cohort_name = cohort.name.clone().unwrap_or_else(|| "cohort".into());
    let output_dir = cohort
        .output_dir
        .as_deref()
        .map(abs)
        .unwrap_or_else(|| base.join(&cohort_name));

    let mut parts: Vec<Part> = Vec::new();
    for s in cohort
        .sdrfs
        .iter()
        .filter(|s| !cohort.exclude.contains(&s.project))
    {
        let path = if s.sdrf.is_absolute() {
            s.sdrf.clone()
        } else {
            base.join(&s.sdrf)
        };
        let text = match std::fs::read_to_string(&path) {
            Ok(t) => t,
            Err(e) => {
                problems.errors.push(format!(
                    "{}: cannot read {}: {e}",
                    s.project,
                    path.display()
                ));
                continue;
            }
        };
        let dir = path.parent().unwrap_or_else(|| Path::new("."));
        match sdrf::split(&text, dir) {
            Ok(groups) => {
                if groups.len() > 1 {
                    problems.warnings.push(format!(
                        "{}: {} splits into {} settings groups; each is searched with its own \
                         settings",
                        s.project,
                        path.display(),
                        groups.len()
                    ));
                }
                for settings in groups {
                    parts.push(Part {
                        project: s.project.clone(),
                        sdrf: path.clone(),
                        settings,
                    });
                }
            }
            Err(e) => problems
                .errors
                .push(format!("{}: {}: {e}", s.project, path.display())),
        }
    }

    let mut runs: Vec<PlannedRun> = Vec::new();
    let mut keys: Vec<((Kind, String, String), String)> = Vec::new();
    let mut taken: BTreeSet<String> = BTreeSet::new();
    let mut names: BTreeSet<String> = BTreeSet::new();
    let mut owner: BTreeMap<PathBuf, String> = BTreeMap::new();
    for part in &parts {
        let ov = cohort
            .overrides
            .iter()
            .find(|o| o.project == part.project)
            .map(|o| &o.search);
        let Some(r) = resolve(part, ov, cfg, &mut problems) else {
            continue;
        };
        let key = r.group_key();
        let group = match keys.iter().find(|(k, _)| *k == key) {
            Some((_, g)) => g.clone(),
            None => {
                let g = group_name(&key, &taken);
                taken.insert(g.clone());
                keys.push((key, g.clone()));
                g
            }
        };
        let args = r.args();
        let mut files: Vec<PathBuf> = Vec::new();
        let mut missing: Vec<String> = Vec::new();
        for f in &part.settings.data_files {
            let name = f
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();
            match (&data_dir, &index) {
                (Some(root), Some(idx)) => match locate(idx, root, &part.project, &projects, &name)
                {
                    Ok(p) => files.push(p),
                    Err(e) => missing.push(e),
                },
                (Some(_), None) => {}
                (None, _) => {
                    if !f.exists() {
                        missing.push(format!("{name}: not found next to the SDRF"));
                    }
                    files.push(f.clone());
                }
            }
        }
        if let Some(first) = missing.first() {
            let msg = format!(
                "{}: {} of {} data files in {} are missing ({first}{})",
                part.project,
                missing.len(),
                part.settings.data_files.len(),
                part.sdrf.display(),
                if missing.len() > 1 { "; ..." } else { "" }
            );
            if data_dir.is_some() {
                problems.errors.push(msg);
            } else {
                problems.warnings.push(format!(
                    "{msg}; set cohort.data_dir to the folder that holds them"
                ));
            }
        }
        let unread: Vec<&PathBuf> = files.iter().filter(|f| !readable(f)).collect();
        if let Some(first) = unread.first() {
            problems.warnings.push(format!(
                "{}: {} of {} data files in {} are not spectrum files andes reads (e.g. {}); \
                 convert or extract each to mzML, Thermo .raw, Bruker .d or MGF with the same \
                 stem and plan again",
                part.project,
                unread.len(),
                files.len(),
                part.sdrf.display(),
                first.display()
            ));
        }
        for file in &files {
            let run = run_name(template, &part.project, file);
            let id = std::fs::canonicalize(file).unwrap_or_else(|_| file.clone());
            if let Some(first) = owner.get(&id) {
                problems.errors.push(format!(
                    "{} is the spectrum file of both {first} and {run}; give each project its own \
                     copy (<data_dir>/<project>/) or list the file in one SDRF",
                    file.display()
                ));
            } else {
                owner.insert(id, run.clone());
            }
            if !names.insert(run.clone()) {
                problems.errors.push(format!(
                    "run name {run} is given to two spectrum files (from {}); make cohort.run_name \
                     include {{project}} or rename the files",
                    part.sdrf.display()
                ));
            }
            let fname = file
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();
            let (_, ext) = split_name(&fname);
            runs.push(PlannedRun {
                stage_as: format!("{run}{ext}"),
                dir: output_dir.join("runs").join(&run),
                run,
                project: part.project.clone(),
                sdrf: part.sdrf.clone(),
                file: file.clone(),
                kind: r.kind.clone(),
                group: group.clone(),
                args: args.clone(),
            });
        }
    }

    if !problems.errors.is_empty() {
        return Err(problems);
    }
    let groups = keys
        .iter()
        .map(|(key, name)| {
            let members: Vec<&PlannedRun> = runs.iter().filter(|r| &r.group == name).collect();
            let mut projects: Vec<String> = Vec::new();
            for r in &members {
                if !projects.contains(&r.project) {
                    projects.push(r.project.clone());
                }
            }
            let input = |r: &PlannedRun| match key.0 {
                Kind::Xl => r.dir.join(format!("{}.xl.tsv", r.run)),
                _ => r.dir.join(format!("{}.pin", r.run)),
            };
            PlannedGroup {
                name: name.clone(),
                kind: key.0.clone(),
                stages: key.0.stages(),
                projects,
                runs: members.iter().map(|r| r.run.clone()).collect(),
                dir: output_dir.join("groups").join(name),
                inputs: members.iter().map(|r| input(r)).collect(),
            }
        })
        .collect();
    Ok(Plan {
        cohort: cohort_name,
        output_dir,
        runs,
        groups,
        warnings: problems.warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "source name\tcomment[data file]\tcomment[cleavage agent details]\t\
                        comment[modification parameters]\tcomment[precursor mass tolerance]\t\
                        comment[label]\tcomment[cross-linker]";

    fn sdrf(rows: &[(&str, &str, &str)]) -> String {
        let mut s = String::from(HEAD);
        for (file, tol, xl) in rows {
            s.push_str(&format!(
                "\ns\t{file}\tNT=Trypsin;AC=MS:1001251\tNT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed\t\
                 {tol}\tlabel free sample\t{xl}"
            ));
        }
        s
    }

    const DSSO: &str = "NT=DSSO;AC=XLMOD:02126;CL=yes;TA=K,S,T,Y,nterm;MH=54.01;ML=85.98";

    fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        p
    }

    fn cfg(yaml: &str) -> RunConfig {
        RunConfig::from_yaml_str(yaml).unwrap()
    }

    #[test]
    fn two_projects_share_a_group_and_get_distinct_run_names() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(
            d,
            "a.sdrf.tsv",
            &sdrf(&[("r1.raw", "10 ppm", DSSO), ("r2.raw", "10 ppm", DSSO)]),
        );
        std::fs::create_dir_all(d.join("b")).unwrap();
        write(d, "b/b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", DSSO)]));
        write(d, "r1.raw", "");
        write(d, "r2.raw", "");
        write(d, "b/r1.raw", "");
        let c = cfg("cohort:\n  name: t\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b/b.sdrf.tsv}\n");
        let p = plan(&c, d).unwrap();
        let names: Vec<&str> = p.runs.iter().map(|r| r.run.as_str()).collect();
        assert_eq!(names, ["A.r1", "A.r2", "B.r1"]);
        assert_eq!(p.groups.len(), 1);
        assert_eq!(p.groups[0].kind, Kind::Xl);
        assert_eq!(p.groups[0].stages, ["xl-fdr"]);
        assert_eq!(p.groups[0].projects, ["A", "B"]);
        assert_eq!(p.runs[0].stage_as, "A.r1.raw");
        assert!(p.warnings.is_empty(), "{:?}", p.warnings);
        let a = &p.runs[0].args;
        let value = |flag: &str| &a[a.iter().position(|x| x == flag).unwrap() + 1];
        crate::cli::parse_xl(value("--xl")).unwrap();
        crate::cli::parse_fixed_mods(value("--fixed-mods")).unwrap();
        assert_eq!(value("--variable-mods"), "", "no variable mods in the SDRF");
        crate::cli::parse_variable_mods(value("--variable-mods")).unwrap();
    }

    #[test]
    fn disagreeing_rows_split_into_settings_groups() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(
            d,
            "a.sdrf.tsv",
            &sdrf(&[("r1.raw", "10 ppm", ""), ("r2.raw", "20 ppm", "")]),
        );
        let c = cfg("cohort:\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n");
        let p = plan(&c, d).unwrap();
        assert_eq!(p.runs.len(), 2);
        let tol = |r: &PlannedRun| {
            let i = r.args.iter().position(|a| a == "--precursor-tol").unwrap();
            r.args[i + 1].clone()
        };
        assert_eq!(tol(&p.runs[0]), "10ppm");
        assert_eq!(tol(&p.runs[1]), "20ppm");
        assert_eq!(
            p.groups.len(),
            1,
            "precursor tolerance does not split the FDR group"
        );
        assert_eq!(p.groups[0].kind, Kind::Lfq);
        assert!(p.warnings.iter().any(|w| w.contains("2 settings groups")));
        assert!(sdrf::parse(&std::fs::read_to_string(d.join("a.sdrf.tsv")).unwrap(), d).is_err());
    }

    #[test]
    fn every_problem_is_reported_at_once() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(
            d,
            "b.sdrf.tsv",
            &sdrf(&[("r1.raw", "10 ppm", "NT=unknown reagent;AC=XLMOD:99999")]),
        );
        let c = cfg("cohort:\n  run_name: \"{stem}\"\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n    - {project: C, sdrf: missing.sdrf.tsv}\n  overrides:\n    - {project: Z, search: {enzyme: trypsin}}\n");
        let e = plan(&c, d).unwrap_err().errors;
        assert!(
            e.iter()
                .any(|m| m.starts_with("B:") && m.contains("cross-linker")),
            "{e:?}"
        );
        assert!(
            e.iter()
                .any(|m| m.starts_with("C:") && m.contains("cannot read")),
            "{e:?}"
        );
        assert!(e.iter().any(|m| m.contains("overrides: Z")), "{e:?}");
    }

    #[test]
    fn duplicate_run_names_are_an_error() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        let c = cfg("cohort:\n  run_name: \"{stem}\"\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n");
        let e = plan(&c, d).unwrap_err().errors;
        assert!(e.iter().any(|m| m.contains("run name r1")), "{e:?}");
    }

    #[test]
    fn override_replaces_the_sdrf_and_config_keys_win_over_the_sdrf() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "", DSSO)]));
        let c = cfg("search: {enzyme: lysc}\ncohort:\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n  overrides:\n    - {project: A, search: {precursor_tol: 5ppm, xl: DSBU}}\n");
        let p = plan(&c, d).unwrap();
        let a = &p.runs[0].args;
        assert!(
            !a.contains(&"--enzyme".to_string()),
            "config enzyme applies: {a:?}"
        );
        let i = a.iter().position(|x| x == "--precursor-tol").unwrap();
        assert_eq!(a[i + 1], "5ppm");
        let i = a.iter().position(|x| x == "--xl").unwrap();
        assert_eq!(a[i + 1], "DSBU");
        assert_eq!(p.groups[0].name, "xl-dsbu-lysc");
    }

    #[test]
    fn excluded_projects_are_not_planned() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "b.sdrf.tsv", &sdrf(&[("r9.raw", "10 ppm", "")]));
        let c = cfg("cohort:\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n  exclude: [B]\n");
        let p = plan(&c, d).unwrap();
        assert_eq!(p.runs.len(), 1);
        assert_eq!(p.runs[0].project, "A");
    }

    #[test]
    fn unreadable_data_files_are_a_warning() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(
            d,
            "a.sdrf.tsv",
            &sdrf(&[("r1.zip", "10 ppm", ""), ("r2.d.zip", "10 ppm", "")]),
        );
        let c = cfg("cohort:\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n");
        let w = plan(&c, d).unwrap().warnings;
        assert!(
            w.iter()
                .any(|m| m.contains("1 of 2 data files") && m.contains("r1.zip")),
            "{w:?}"
        );
        assert!(readable(Path::new("x.mzML.gz")) && readable(Path::new("x.d")));
    }

    #[test]
    fn data_files_are_found_under_data_dir_by_name() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", DSSO)]));
        write(
            d,
            "b.sdrf.tsv",
            &sdrf(&[("r1.raw", "10 ppm", DSSO), ("r2.d.zip", "10 ppm", DSSO)]),
        );
        for sub in ["raw/A", "raw/B", "raw/other"] {
            std::fs::create_dir_all(d.join(sub)).unwrap();
        }
        write(d, "raw/A/r1.mzML", "");
        write(d, "raw/B/r1.raw", "");
        write(d, "raw/other/r1.raw", "");
        write(d, "raw/B/r2.d.zip", "");
        let yaml = "cohort:\n  name: t\n  data_dir: raw\n  output_dir: out\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n";
        let p = plan(&cfg(yaml), d).unwrap();
        let file = |run: &str| p.runs.iter().find(|r| r.run == run).unwrap();
        assert_eq!(file("A.r1").file, d.join("raw/A/r1.mzML"), "converted copy");
        assert_eq!(file("A.r1").stage_as, "A.r1.mzML");
        assert_eq!(
            file("B.r1").file,
            d.join("raw/B/r1.raw"),
            "project folder wins"
        );
        assert_eq!(file("B.r2").stage_as, "B.r2.d.zip");
        assert_eq!(file("B.r2").dir, d.join("out/runs/B.r2"));
        assert_eq!(
            p.groups[0].dir,
            d.join("out/groups").join(&p.groups[0].name)
        );
        assert!(p.groups[0]
            .inputs
            .contains(&d.join("out/runs/A.r1/A.r1.xl.tsv")));

        std::fs::remove_file(d.join("raw/B/r1.raw")).unwrap();
        write(
            d,
            "b.sdrf.tsv",
            &sdrf(&[("r1.raw", "10 ppm", DSSO), ("r9.raw", "10 ppm", DSSO)]),
        );
        let e = plan(&cfg(yaml), d).unwrap_err().errors;
        assert!(
            e.iter()
                .any(|m| m.contains("B:") && m.contains("2 of 2 data files")),
            "{e:?}"
        );
    }

    #[test]
    fn two_projects_with_the_same_file_name_get_their_own_file_and_run() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        let yaml = "cohort:\n  data_dir: raw\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n";
        std::fs::create_dir_all(d.join("raw")).unwrap();
        write(d, "raw/r1.raw", "");
        let flat = t.path().join("raw/A");
        std::fs::create_dir_all(&flat).unwrap();
        std::fs::create_dir_all(d.join("raw/B")).unwrap();
        write(d, "raw/A/r1.raw", "");
        write(d, "raw/B/r1.raw", "");
        let p = plan(&cfg(yaml), d).unwrap();
        let a = p.runs.iter().find(|r| r.project == "A").unwrap();
        let b = p.runs.iter().find(|r| r.project == "B").unwrap();
        assert_eq!((a.run.as_str(), b.run.as_str()), ("A.r1", "B.r1"));
        assert_eq!(a.file, d.join("raw/A/r1.raw"));
        assert_eq!(b.file, d.join("raw/B/r1.raw"));
        assert_ne!(a.stage_as, b.stage_as);

        std::fs::remove_dir_all(d.join("raw/A")).unwrap();
        std::fs::remove_dir_all(d.join("raw/B")).unwrap();
        std::fs::create_dir_all(d.join("raw/x")).unwrap();
        write(d, "raw/x/r1.raw", "");
        let e = plan(&cfg(yaml), d).unwrap_err().errors;
        assert!(e.iter().any(|m| m.contains("2 files under")), "{e:?}");
    }

    #[test]
    fn one_spectrum_file_given_to_two_projects_is_an_error() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "r1.raw", "");
        let c = cfg("cohort:\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n");
        let e = plan(&c, d).unwrap_err().errors;
        assert!(
            e.iter()
                .any(|m| m.contains("spectrum file of both A.r1 and B.r1")),
            "{e:?}"
        );
    }

    #[test]
    fn a_project_without_a_folder_does_not_take_another_projects_file() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        std::fs::create_dir_all(d.join("raw/B")).unwrap();
        write(d, "raw/B/r1.raw", "");
        let yaml = "cohort:\n  data_dir: raw\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n";
        let e = plan(&cfg(yaml), d).unwrap_err().errors;
        assert!(
            e.iter()
                .any(|m| m.starts_with("A:") && m.contains("1 of 1 data files")),
            "{e:?}"
        );
        write(d, "raw/r1.raw", "");
        let p = plan(&cfg(yaml), d).unwrap();
        let file = |run: &str| p.runs.iter().find(|r| r.run == run).unwrap().file.clone();
        assert_eq!(file("A.r1"), d.join("raw/r1.raw"));
        assert_eq!(file("B.r1"), d.join("raw/B/r1.raw"));
    }

    #[test]
    fn the_config_quant_section_sets_the_kind_of_data() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "r1.raw", "");
        let sdrfs = "cohort:\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n";
        let p = plan(&cfg(sdrfs), d).unwrap();
        assert_eq!(p.groups[0].kind, Kind::Lfq, "the SDRF says label free");
        assert!(p.runs[0].args.contains(&"--lfq".to_string()));

        let p = plan(&cfg(&format!("quant: {{tmt: tmt10}}\n{sdrfs}")), d).unwrap();
        assert_eq!(p.groups[0].kind, Kind::Tmt);
        assert_eq!(p.groups[0].stages, ["percolator", "inference", "tmt-quant"]);
        assert!(
            !p.runs[0].args.iter().any(|a| a == "--lfq" || a == "--tmt"),
            "the runner's --config applies the config: {:?}",
            p.runs[0].args
        );
        assert!(
            p.warnings.iter().any(|w| w.contains("makes its runs tmt")),
            "{:?}",
            p.warnings
        );

        let p = plan(&cfg(&format!("quant: {{lfq: true}}\n{sdrfs}")), d).unwrap();
        assert_eq!(p.groups[0].kind, Kind::Lfq);
        assert!(
            p.warnings.is_empty(),
            "config and SDRF agree: {:?}",
            p.warnings
        );
        let p = plan(&cfg(&format!("quant: {{lfq: false}}\n{sdrfs}")), d).unwrap();
        assert_eq!(p.groups[0].kind, Kind::Dda);
    }

    #[cfg(unix)]
    #[test]
    fn linked_folders_are_followed_once_and_keep_their_project() {
        let t = tempfile::tempdir().unwrap();
        let d = &std::fs::canonicalize(t.path()).unwrap();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        std::fs::create_dir_all(d.join("store/b")).unwrap();
        std::fs::create_dir_all(d.join("raw/A")).unwrap();
        write(d, "raw/A/r1.raw", "");
        write(d, "store/b/r1.raw", "");
        std::os::unix::fs::symlink(d.join("store/b"), d.join("raw/B")).unwrap();
        std::os::unix::fs::symlink(d.join("store/b"), d.join("raw/A/alias")).unwrap();
        std::os::unix::fs::symlink(d.join("raw"), d.join("raw/loop")).unwrap();
        let yaml = "cohort:\n  data_dir: raw\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n";
        let p = plan(&cfg(yaml), d).unwrap();
        let file = |run: &str| p.runs.iter().find(|r| r.run == run).unwrap().file.clone();
        assert_eq!(file("A.r1"), d.join("raw/A/r1.raw"));
        assert_eq!(file("B.r1"), d.join("store/b/r1.raw"));
        let index = index_files(&d.join("raw")).unwrap();
        assert_eq!(index["r1.raw"].len(), 2, "{:?}", index["r1.raw"]);
    }

    #[test]
    fn defaults_come_from_the_flags() {
        assert_eq!(
            crate::cli::default_value::<crate::cli::DatabaseArgs>("enzyme"),
            "trypsin"
        );
        assert!(!crate::cli::default_value::<crate::cli::SearchArgs>("precursor_tol").is_empty());
    }

    #[test]
    fn file_names_split_at_the_format_suffix() {
        assert_eq!(split_name("x.d.zip"), ("x", ".d.zip"));
        assert_eq!(split_name("run.1.RAW"), ("run.1", ".RAW"));
        assert_eq!(split_name("a.mzML.gz"), ("a", ".mzML.gz"));
        assert_eq!(split_name("noext"), ("noext", ""));
    }

    #[test]
    fn run_name_template_is_checked() {
        assert!(check_template("{project}.{stem}").is_ok());
        assert!(check_template("{project}").is_err());
        assert!(check_template("{stem}.{lab}").is_err());
        assert!(check_template("{project}/{stem}").is_err());
    }
}
