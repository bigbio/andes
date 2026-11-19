//! Cohorts: the SDRFs of several projects analysed together.
//!
//! `andes cohort plan` reads the `cohort` section of a `--config` file and every SDRF it
//! lists, resolves each run's search settings, groups the runs that can share one FDR, and
//! writes the plan as JSON. It searches nothing; a job array or the pipeline runs the plan.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::config::RunConfig;
use crate::sdrf::{self, Quant, SdrfSettings};

/// Run name template used when the cohort section gives none.
pub const DEFAULT_RUN_NAME: &str = "{project}.{stem}";

/// The enzyme of a run whose SDRF and config give none: the `--enzyme` default.
const DEFAULT_ENZYME: &str = "trypsin";

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CohortCfg {
    pub name: Option<String>,
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
    /// The spectrum file as the SDRF names it, resolved next to the SDRF.
    pub file: PathBuf,
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
}

#[derive(Debug, Serialize)]
pub struct Plan {
    pub cohort: String,
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

fn list(m: &[String]) -> String {
    if m.is_empty() {
        "none".into()
    } else {
        m.join(", ")
    }
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
    check("fixed_mods", &o.fixed_mods, &|v| {
        crate::cli::parse_fixed_mods(v).map(|_| ())
    });
    check("variable_mods", &o.variable_mods, &|v| {
        crate::cli::parse_variable_mods(v).map(|_| ())
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
            "{p}: {} gives no precursor tolerance; its runs use the default 20 ppm \
             (set cohort.overrides search.precursor_tol to change it)",
            part.sdrf.display()
        ));
    }
    let to_list = |v: Option<String>| {
        v.map(|v| {
            if v == "none" {
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
        match &s.quant {
            Some(Quant::Lfq) => (Kind::Lfq, Some(("--lfq".into(), None))),
            Some(Quant::Isobaric(plex)) => (Kind::Tmt, Some(("--tmt".into(), Some(plex.clone())))),
            Some(Quant::Label(spec)) => (Kind::Label, Some(("--label".into(), Some(spec.clone())))),
            None => (Kind::Dda, None),
        }
    };
    let group_enzyme = enzyme
        .clone()
        .or_else(|| sc.enzyme.clone())
        .unwrap_or_else(|| DEFAULT_ENZYME.into());
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
    let stem = file
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    template
        .replace("{project}", project)
        .replace("{stem}", &stem)
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
        let unread: Vec<&PathBuf> = part
            .settings
            .data_files
            .iter()
            .filter(|f| !readable(f))
            .collect();
        if let Some(first) = unread.first() {
            problems.warnings.push(format!(
                "{}: {} of {} data files in {} are not spectrum files andes reads (e.g. {}); \
                 convert or extract each and stage the result under its run's stage_as name",
                part.project,
                unread.len(),
                part.settings.data_files.len(),
                part.sdrf.display(),
                first.display()
            ));
        }
        for file in &part.settings.data_files {
            let run = run_name(template, &part.project, file);
            if !names.insert(run.clone()) {
                problems.errors.push(format!(
                    "run name {run} is given to two spectrum files (from {}); make cohort.run_name \
                     include {{project}} or rename the files",
                    part.sdrf.display()
                ));
            }
            let ext = file
                .extension()
                .map(|e| format!(".{}", e.to_string_lossy()))
                .unwrap_or_default();
            runs.push(PlannedRun {
                stage_as: format!("{run}{ext}"),
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
            PlannedGroup {
                name: name.clone(),
                kind: key.0.clone(),
                stages: key.0.stages(),
                projects,
                runs: members.iter().map(|r| r.run.clone()).collect(),
            }
        })
        .collect();
    Ok(Plan {
        cohort: cohort.name.clone().unwrap_or_else(|| "cohort".into()),
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
        let d = t.path();
        write(
            d,
            "a.sdrf.tsv",
            &sdrf(&[("r1.raw", "10 ppm", DSSO), ("r2.raw", "10 ppm", DSSO)]),
        );
        write(d, "b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", DSSO)]));
        let c = cfg("cohort:\n  name: t\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n");
        let p = plan(&c, d).unwrap();
        let names: Vec<&str> = p.runs.iter().map(|r| r.run.as_str()).collect();
        assert_eq!(names, ["A.r1", "A.r2", "B.r1"]);
        assert_eq!(p.groups.len(), 1);
        assert_eq!(p.groups[0].kind, Kind::Xl);
        assert_eq!(p.groups[0].stages, ["xl-fdr"]);
        assert_eq!(p.groups[0].projects, ["A", "B"]);
        assert_eq!(p.runs[0].stage_as, "A.r1.raw");
        assert!(p.warnings.is_empty(), "{:?}", p.warnings);
        let xl = p.runs[0].args.iter().position(|a| a == "--xl").unwrap();
        crate::cli::parse_xl(&p.runs[0].args[xl + 1]).unwrap();
    }

    #[test]
    fn disagreeing_rows_split_into_settings_groups() {
        let t = tempfile::tempdir().unwrap();
        let d = t.path();
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
        let d = t.path();
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
        let d = t.path();
        write(d, "a.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        write(d, "b.sdrf.tsv", &sdrf(&[("r1.raw", "10 ppm", "")]));
        let c = cfg("cohort:\n  run_name: \"{stem}\"\n  sdrfs:\n    - {project: A, sdrf: a.sdrf.tsv}\n    - {project: B, sdrf: b.sdrf.tsv}\n");
        let e = plan(&c, d).unwrap_err().errors;
        assert!(e.iter().any(|m| m.contains("run name r1")), "{e:?}");
    }

    #[test]
    fn override_replaces_the_sdrf_and_config_keys_win_over_the_sdrf() {
        let t = tempfile::tempdir().unwrap();
        let d = t.path();
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
        let d = t.path();
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
        let d = t.path();
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
    fn run_name_template_is_checked() {
        assert!(check_template("{project}.{stem}").is_ok());
        assert!(check_template("{project}").is_err());
        assert!(check_template("{stem}.{lab}").is_err());
        assert!(check_template("{project}/{stem}").is_err());
    }
}
