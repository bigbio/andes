//! Tests that `ParamSources` records where each search setting came from, with the value as
//! written, through the same merge order the search command uses.

use clap::{CommandFactory, FromArgMatches};
use serde_yaml::Value;

use crate::config::{apply, RunConfig};
use crate::param_sources::{from_command_line, ParamSources, Source};
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

fn flag_text(arg: &clap::Arg) -> Option<String> {
    arg.get_long().map(|l| format!("--{l}"))
}

/// Values worth trying for an argument, most specific first: the row's example, defaults, the listed possible
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

#[test]
fn every_cli_value_records_as_a_config_value() {
    let mut cmd = crate::TopCli::command();
    cmd.build();
    let base = [
        ("--database", "y"),
        ("--output-pin", "o"),
        ("--spectrum", "x"),
    ];
    let mut failures = Vec::new();
    for p in crate::params::PARAMS
        .iter()
        .filter(|p| p.section != crate::params::Section::CliOnly)
    {
        if EXEMPT.iter().any(|(id, _)| *id == p.id) {
            continue;
        }
        let Some(arg) = cmd.get_arguments().find(|a| a.get_id() == p.id) else {
            failures.push(format!("{}: no clap argument", p.id));
            continue;
        };
        let Some(flag) = flag_text(arg) else {
            failures.push(format!("{}: no long flag", p.id));
            continue;
        };
        let delimiter = arg.get_value_delimiter();
        let mut tried: Vec<Vec<String>> = Vec::new();
        if takes_no_value(arg) {
            tried.push(vec![flag.clone()]);
        } else {
            let pieces = candidates(arg, p.unit, p.example);
            for (i, piece) in pieces.iter().enumerate() {
                let value = match delimiter {
                    Some(d) => format!("{piece}{d}{}", pieces.get(i + 1).unwrap_or(piece)),
                    None => piece.clone(),
                };
                tried.push(vec![flag.clone(), value]);
            }
        }
        let mut last_error = String::new();
        let mut parsed = None;
        for t in &tried {
            let mut all: Vec<String> = vec!["andes".to_string()];
            for &(f, v) in base.iter().filter(|b| b.0 != flag) {
                all.extend([f.to_string(), v.to_string()]);
            }
            all.extend(t.iter().cloned());
            match crate::TopCli::command().try_get_matches_from(all) {
                Ok(m) => {
                    parsed = Some(m);
                    break;
                }
                Err(e) => last_error = e.to_string(),
            }
        }
        let Some(m) = parsed else {
            failures.push(format!("{}: no sample parses ({last_error})", p.id));
            continue;
        };
        let sources = from_command_line(&m);
        let Some(value) = sources.value(p.id) else {
            failures.push(format!("{}: nothing recorded", p.id));
            continue;
        };
        if delimiter.is_some() && !matches!(value, Value::Sequence(v) if v.len() == 2) {
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
