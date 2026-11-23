//! Tests that `ParamSources` records where each search setting came from, with the value as
//! written, through the same merge order the search command uses.

use clap::{CommandFactory, FromArgMatches};
use serde_yaml::Value;

use crate::config::{apply, RunConfig};
use crate::param_sources::{from_command_line, ParamSources, Source};
use crate::params_render::effective_config;
use crate::sdrf::{self, SdrfSettings};

fn baseline_matches(extra: &[&str]) -> clap::ArgMatches {
    let base = [
        "andes",
        "--spectrum",
        "x",
        "--database",
        "y",
        "--output-pin",
        "o",
    ];
    crate::TopCli::command()
        .try_get_matches_from(base.iter().chain(extra))
        .expect("parse args")
}

/// The sources after the command line, the config text and an SDRF setting (when given) are
/// merged in the order of the search command.
fn merged(yaml: &str, extra: &[&str], sdrf: Option<SdrfSettings>) -> ParamSources {
    let m = baseline_matches(extra);
    let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
    let mut sources = from_command_line(&m);
    let cfg = RunConfig::from_yaml_str(yaml).unwrap();
    apply(cfg, &mut a, &m, &mut sources).unwrap();
    if let Some(s) = sdrf {
        sdrf::apply(&s, &mut a, &m, &mut sources).unwrap();
    }
    sources
}

fn text(v: &str) -> Value {
    Value::String(v.to_string())
}

#[test]
fn command_line_value_is_recorded_as_command_line() {
    let s = merged("{}", &["--precursor-tol", "10ppm"], None);
    assert_eq!(s.source("precursor_tol"), Source::CommandLine);
    assert_eq!(s.value("precursor_tol"), Some(&text("10ppm")));
}

#[test]
fn config_value_is_recorded_with_its_yaml_value() {
    let s = merged("search:\n  precursor_tol: 10ppm\n", &[], None);
    assert_eq!(s.source("precursor_tol"), Source::Config);
    assert_eq!(s.value("precursor_tol"), Some(&text("10ppm")));
}

#[test]
fn command_line_beats_config_in_sources() {
    let s = merged(
        "search:\n  precursor_tol: 10ppm\n",
        &["--precursor-tol", "5ppm"],
        None,
    );
    assert_eq!(s.source("precursor_tol"), Source::CommandLine);
    assert_eq!(s.value("precursor_tol"), Some(&text("5ppm")));
}

#[test]
fn config_beats_sdrf_in_sources() {
    let settings = SdrfSettings {
        precursor_tol: Some("30ppm".to_string()),
        enzyme: Some("trypsin".to_string()),
        ..SdrfSettings::default()
    };
    let s = merged("search:\n  precursor_tol: 10ppm\n", &[], Some(settings));
    assert_eq!(s.source("precursor_tol"), Source::Config);
    assert_eq!(s.value("precursor_tol"), Some(&text("10ppm")));
    assert_eq!(s.source("enzyme"), Source::Sdrf);
    assert_eq!(s.value("enzyme"), Some(&text("trypsin")));
}

#[test]
fn untouched_flag_is_default() {
    let s = merged("{}", &[], None);
    assert_eq!(s.source("precursor_tol"), Source::Default);
    assert_eq!(s.value("precursor_tol"), None);
}

/// The YAML text `<section>:\n  <key>: <value>` that would set `id` to its recorded value.
fn config_text(id: &str, value: &Value) -> Option<String> {
    let p = crate::params::spec(id)?;
    let mut inner = serde_yaml::Mapping::new();
    inner.insert(text(p.key), value.clone());
    let mut outer = serde_yaml::Mapping::new();
    outer.insert(text(p.section.yaml_name()?), Value::Mapping(inner));
    serde_yaml::to_string(&Value::Mapping(outer)).ok()
}

#[test]
fn recorded_values_parse_as_config() {
    let cases: [&[&str]; 4] = [
        &["--spectrum", "a.mzML"],
        &["--spectrum", "a.mzML", "--spectrum", "b.mzML"],
        &["--fixed-mods", "Carbamidomethyl (C)", "--lfq"],
        &["--min-peaks", "10", "--tmt-min-purity", "0.5", "--precursor-tol", "20ppm"],
    ];
    for args in cases {
        let all: Vec<&str> = ["andes", "--database", "y", "--output-pin", "o"]
            .into_iter()
            .chain(args.iter().copied())
            .collect();
        let m = crate::TopCli::command().try_get_matches_from(all).unwrap();
        let sources = from_command_line(&m);
        for p in crate::params::PARAMS {
            let Some(value) = sources.value(p.id) else {
                continue;
            };
            let Some(yaml) = config_text(p.id, value) else {
                continue;
            };
            let cfg = RunConfig::from_yaml_str(&yaml)
                .unwrap_or_else(|e| panic!("{} recorded as {yaml:?}: {e}", p.id));
            match p.id {
                "spectrum" => {
                    let want: Vec<std::path::PathBuf> = args
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i > 0 && args[i - 1] == "--spectrum")
                        .map(|(_, a)| std::path::PathBuf::from(a))
                        .collect();
                    assert_eq!(cfg.io.spectrum, Some(want));
                }
                "fixed_mods" => {
                    assert_eq!(cfg.search.fixed_mods, Some(vec![text("Carbamidomethyl (C)")]));
                }
                "lfq" => assert_eq!(cfg.quant.lfq, Some(true)),
                "min_peaks" => assert_eq!(value, &Value::from(10)),
                "tmt_min_purity" => assert_eq!(value, &Value::from(0.5)),
                "precursor_tol" => assert_eq!(value, &text("20ppm")),
                _ => {}
            }
        }
    }
}

/// Parameters whose sample cannot be given alone, with the reason; none today.
const EXEMPT: &[(&str, &str)] = &[];

/// Arguments a parameter's sample needs beside it (clap `requires`), in the sample run and in
/// the baseline it is replayed onto.
const COMPANIONS: &[(&str, &[&str])] = &[("xl_linear_only", &["--xl", "DSSO"])];

/// The arguments every sample run has, except the sampled flag itself.
const BASE: [(&str, &str); 3] = [
    ("--database", "y"),
    ("--output-pin", "o"),
    ("--spectrum", "x"),
];

fn flag_text(arg: &clap::Arg) -> Option<String> {
    arg.get_long().map(|l| format!("--{l}"))
}

/// Values worth trying for an argument, most specific first: the row's example, defaults, the
/// listed possible
/// values (`false` before `true`), then a sample of the unit.
fn candidates(arg: &clap::Arg, unit: crate::params::Unit, example: Option<&str>) -> Vec<String> {
    let mut out: Vec<String> = example.map(String::from).into_iter().collect();
    out.extend(
        arg.get_default_values()
            .iter()
            .map(|v| v.to_string_lossy().into_owned()),
    );
    let mut possible: Vec<String> = arg
        .get_possible_values()
        .iter()
        .map(|v| v.get_name().to_string())
        .collect();
    possible.sort_by_key(|v| v != "false");
    out.extend(possible);
    out.push(crate::params_table::unit_sample(unit).to_string());
    out.extend(["x", "1", "0.5", "20ppm"].map(String::from));
    out
}

fn takes_no_value(arg: &clap::Arg) -> bool {
    arg.get_num_args().is_some_and(|n| n.max_values() == 0)
}

/// The command line of a sample run without the sample: `andes`, the base arguments other than
/// `flag`, and the companions of `id`.
fn baseline_args(flag: &str, id: &str) -> Vec<String> {
    let mut all: Vec<String> = vec!["andes".to_string()];
    for &(f, v) in BASE.iter().filter(|b| b.0 != flag) {
        all.extend([f.to_string(), v.to_string()]);
    }
    if let Some((_, extra)) = COMPANIONS.iter().find(|c| c.0 == id) {
        all.extend(extra.iter().map(|e| e.to_string()));
    }
    all
}

/// The matches of the first sample of `p` that clap accepts, with the flag text. With
/// `changed_first` the candidates that differ from the clap default are tried first, so that a
/// replayed value has something to change.
fn sample(
    cmd: &clap::Command,
    p: &crate::params::ParamSpec,
    changed_first: bool,
) -> Result<(clap::ArgMatches, String), String> {
    let arg = cmd
        .get_arguments()
        .find(|a| a.get_id() == p.id)
        .ok_or_else(|| format!("{}: no clap argument", p.id))?;
    let flag = flag_text(arg).ok_or_else(|| format!("{}: no long flag", p.id))?;
    let delimiter = arg.get_value_delimiter();
    let mut tried: Vec<Vec<String>> = Vec::new();
    if takes_no_value(arg) {
        tried.push(vec![flag.clone()]);
    } else {
        let mut pieces = candidates(arg, p.unit, p.example);
        if changed_first {
            let defaults: Vec<String> = arg
                .get_default_values()
                .iter()
                .map(|v| v.to_string_lossy().into_owned())
                .collect();
            pieces.sort_by_key(|piece| defaults.contains(piece));
        }
        for (i, piece) in pieces.iter().enumerate() {
            let value = match delimiter {
                Some(d) => format!("{piece}{d}{}", pieces.get(i + 1).unwrap_or(piece)),
                None => piece.clone(),
            };
            tried.push(vec![flag.clone(), value]);
        }
    }
    let mut last_error = String::new();
    for t in &tried {
        let mut all = baseline_args(&flag, p.id);
        all.extend(t.iter().cloned());
        match crate::TopCli::command().try_get_matches_from(all) {
            Ok(m) => return Ok((m, flag)),
            Err(e) => last_error = e.to_string(),
        }
    }
    Err(format!("{}: no sample parses ({last_error})", p.id))
}

/// The rows the exhaustive tests cover: every config parameter not exempted.
fn sampled_rows() -> impl Iterator<Item = &'static crate::params::ParamSpec> {
    crate::params::PARAMS.iter().filter(|p| {
        p.section != crate::params::Section::CliOnly && !EXEMPT.iter().any(|(id, _)| *id == p.id)
    })
}

#[test]
fn every_cli_value_records_as_a_config_value() {
    let mut cmd = crate::TopCli::command();
    cmd.build();
    let mut failures = Vec::new();
    for p in sampled_rows() {
        let (m, _) = match sample(&cmd, p, false) {
            Ok(s) => s,
            Err(e) => {
                failures.push(e);
                continue;
            }
        };
        let delimited = cmd
            .get_arguments()
            .find(|a| a.get_id() == p.id)
            .is_some_and(|a| a.get_value_delimiter().is_some());
        let sources = from_command_line(&m);
        let Some(value) = sources.value(p.id) else {
            failures.push(format!("{}: nothing recorded", p.id));
            continue;
        };
        if delimited && !matches!(value, Value::Sequence(v) if v.len() == 2) {
            failures.push(format!("{}: a delimited list recorded as {value:?}", p.id));
        }
        let Some(yaml) = config_text(p.id, value) else {
            failures.push(format!("{}: no config text", p.id));
            continue;
        };
        if let Err(e) = RunConfig::from_yaml_str(&yaml) {
            failures.push(format!("{}: recorded {yaml:?}: {e}", p.id));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Parameters whose settings cannot be compared through `Debug` after a replay, with the
/// reason; none today.
const REPLAY_EXEMPT: &[(&str, &str)] = &[];

/// The first lines where two pretty `Debug` dumps differ, for a readable failure.
fn first_differences(a: &str, b: &str) -> String {
    let differing: Vec<String> = a
        .lines()
        .zip(b.lines())
        .filter(|(x, y)| x != y)
        .take(3)
        .map(|(x, y)| format!("    command line: {}\n    replayed:     {}", x.trim(), y.trim()))
        .collect();
    if differing.is_empty() {
        "    the dumps differ in length".to_string()
    } else {
        differing.join("\n")
    }
}

#[test]
fn every_cli_value_has_the_same_meaning_after_replay() {
    let mut cmd = crate::TopCli::command();
    cmd.build();
    let mut failures = Vec::new();
    for p in sampled_rows().filter(|p| !REPLAY_EXEMPT.iter().any(|(id, _)| *id == p.id)) {
        let (m, flag) = match sample(&cmd, p, true) {
            Ok(s) => s,
            Err(e) => {
                failures.push(e);
                continue;
            }
        };
        let direct = match crate::SearchArgs::from_arg_matches(&m) {
            Ok(a) => a,
            Err(e) => {
                failures.push(format!("{}: {e}", p.id));
                continue;
            }
        };
        let sources = from_command_line(&m);
        let Some(yaml) = sources.value(p.id).and_then(|v| config_text(p.id, v)) else {
            failures.push(format!("{}: nothing recorded", p.id));
            continue;
        };
        let base_args = baseline_args(&flag, p.id);
        let baseline = match crate::TopCli::command().try_get_matches_from(base_args) {
            Ok(m) => m,
            Err(e) => {
                failures.push(format!("{}: the baseline does not parse: {e}", p.id));
                continue;
            }
        };
        let replayed = crate::SearchArgs::from_arg_matches(&baseline)
            .map_err(|e| e.to_string())
            .and_then(|mut a| {
                let cfg = RunConfig::from_yaml_str(&yaml)?;
                apply(cfg, &mut a, &baseline, &mut ParamSources::default())?;
                Ok(a)
            });
        let replayed = match replayed {
            Ok(a) => a,
            Err(e) => {
                failures.push(format!("{}: replaying {yaml:?}: {e}", p.id));
                continue;
            }
        };
        let (a, b) = (format!("{direct:#?}"), format!("{replayed:#?}"));
        if a != b {
            failures.push(format!(
                "{}: {flag} and the replayed {yaml:?} differ:\n{}",
                p.id,
                first_differences(&a, &b)
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// An SDRF of one label-free run: protease, a fixed and a variable modification, the
/// tolerances and the dissociation method.
const REPLAY_SDRF: &str = "source name\tcomment[data file]\tcomment[label]\t\
comment[cleavage agent details]\tcomment[modification parameters]\t\
comment[modification parameters]\tcomment[precursor mass tolerance]\t\
comment[fragment mass tolerance]\tcomment[dissociation method]\n\
s1\ta.raw\tlabel free sample\tNT=Trypsin;AC=MS:1001251\t\
NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed;PP=Anywhere\t\
NT=Phospho;TA=S,T,Y;MT=Variable;PP=Anywhere\t10 ppm\t0.02 Da\tNT=HCD;AC=PRIDE:0000590\n";

/// The settings of a search merged from `m`, the config text and the SDRF fixture, with their
/// sources, in the order of the search command.
fn merged_with_sdrf(m: &clap::ArgMatches, yaml: &str) -> (crate::SearchArgs, ParamSources) {
    let settings =
        sdrf::parse(REPLAY_SDRF, std::path::Path::new("/d")).expect("the SDRF fixture parses");
    let mut a = crate::SearchArgs::from_arg_matches(m).expect("search args");
    let mut sources = from_command_line(m);
    let cfg = RunConfig::from_yaml_str(yaml).unwrap_or_else(|e| panic!("{e}\n{yaml}"));
    apply(cfg, &mut a, m, &mut sources).expect("the config applies");
    sdrf::apply(&settings, &mut a, m, &mut sources).expect("the SDRF applies");
    (a, sources)
}

#[test]
fn config_and_sdrf_replay_gives_the_same_settings() {
    let m = baseline_matches(&[
        "--sdrf",
        "design.sdrf.tsv",
        "--precursor-tol",
        "5ppm",
        "--top-n",
        "7",
        "--xl",
        "DSSO; DSBU",
    ]);
    let yaml = "search:\n  max_missed_cleavages: 1\n  enzyme: chymotrypsin\n\
                quant:\n  lfq_rt_window: 90\n";
    let (first, sources) = merged_with_sdrf(&m, yaml);
    let printed = effective_config(&sources).expect("printable");

    let again_m = crate::TopCli::command()
        .try_get_matches_from(["andes"])
        .expect("parse args");
    let (again, again_sources) = merged_with_sdrf(&again_m, &printed);
    assert_eq!(
        format!("{first:#?}"),
        format!("{again:#?}"),
        "the replay of this config differs:\n{printed}"
    );
    for id in ["precursor_tol", "enzyme", "fixed_mods", "lfq", "xl"] {
        assert_eq!(
            again_sources.source(id),
            Source::Config,
            "{id} is read from the printed config"
        );
    }
}

#[test]
fn several_crosslinkers_record_as_a_list() {
    let s = merged("{}", &["--xl", "DSSO; DSBU ;"], None);
    assert_eq!(
        s.value("xl"),
        Some(&Value::Sequence(vec![text("DSSO"), text("DSBU")]))
    );
}
