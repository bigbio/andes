//! Keeps the parameter table in `params.rs` consistent with the search command's clap
//! declaration and with the config file schema.

use std::collections::BTreeSet;

use clap::{Arg, ArgAction, CommandFactory};

use crate::config::RunConfig;
use crate::params::{spec, ParamSpec, Section, Unit, Unset, PARAMS};

fn search_command() -> clap::Command {
    let mut cmd = crate::TopCli::command();
    cmd.build();
    cmd
}

fn is_builtin(arg: &Arg) -> bool {
    matches!(arg.get_id().as_str(), "help" | "version")
}

fn cli_ids(cmd: &clap::Command) -> BTreeSet<String> {
    cmd.get_arguments()
        .filter(|a| !is_builtin(a))
        .map(|a| a.get_id().to_string())
        .collect()
}

fn yaml_scalar(v: &str) -> String {
    if v.parse::<f64>().is_ok() || v == "true" || v == "false" {
        v.to_string()
    } else {
        format!("\"{}\"", v.replace('"', "\\\""))
    }
}

pub(crate) fn unit_sample(unit: Unit) -> &'static str {
    match unit {
        Unit::PpmOrDa => "20ppm",
        Unit::Fraction => "0.5",
        Unit::Path | Unit::None => "x",
        Unit::Ppm
        | Unit::Da
        | Unit::Minutes
        | Unit::Seconds
        | Unit::Gb
        | Unit::Count
        | Unit::Number
        | Unit::Charge => "1",
    }
}

fn config_value(arg: &Arg, p: &ParamSpec) -> String {
    if matches!(p.id, "fixed_mods" | "variable_mods") {
        return "[\"Oxidation (M)\"]".to_string();
    }
    let action = arg.get_action();
    if matches!(action, ArgAction::SetTrue) {
        return "true".to_string();
    }
    let defaults: Vec<String> = arg
        .get_default_values()
        .iter()
        .map(|v| yaml_scalar(&v.to_string_lossy()))
        .collect();
    let scalar = defaults
        .first()
        .cloned()
        .unwrap_or_else(|| unit_sample(p.unit).to_string());
    if matches!(action, ArgAction::Append) {
        let items = if defaults.is_empty() {
            scalar
        } else {
            defaults.join(", ")
        };
        return format!("[{items}]");
    }
    scalar
}

#[test]
fn every_search_arg_has_one_param_row() {
    let cmd = search_command();
    let cli = cli_ids(&cmd);
    let mut rows = BTreeSet::new();
    let mut duplicated = Vec::new();
    for p in PARAMS {
        if !rows.insert(p.id.to_string()) {
            duplicated.push(p.id);
        }
    }
    let without_row: Vec<&String> = cli.difference(&rows).collect();
    let without_arg: Vec<&String> = rows.difference(&cli).collect();
    assert!(
        duplicated.is_empty() && without_row.is_empty() && without_arg.is_empty(),
        "params.rs and the CLI disagree\n  arguments with no row: {without_row:?}\n  \
         rows with no argument: {without_arg:?}\n  ids with two rows: {duplicated:?}"
    );
}

#[test]
fn every_config_key_parses() {
    let cmd = search_command();
    let mut failures = Vec::new();
    for p in PARAMS.iter().filter(|p| p.section != Section::CliOnly) {
        let Some(arg) = cmd.get_arguments().find(|a| a.get_id().as_str() == p.id) else {
            failures.push(format!("{}: no such argument", p.id));
            continue;
        };
        let Some(section) = p.section.yaml_name() else {
            continue;
        };
        let yaml = format!("{section}:\n  {}: {}\n", p.key, config_value(arg, p));
        if let Err(e) = RunConfig::from_yaml_str(&yaml) {
            failures.push(format!("{}: {e} in {yaml:?}", p.id));
        }
    }
    assert!(
        failures.is_empty(),
        "config keys that do not parse:\n{}",
        failures.join("\n")
    );
}

#[test]
fn every_search_arg_is_documented() {
    let cmd = search_command();
    let mut problems = Vec::new();
    for arg in cmd.get_arguments().filter(|a| !is_builtin(a)) {
        let id = arg.get_id().as_str();
        let help = arg.get_help().map(|h| h.to_string()).unwrap_or_default();
        let long = arg
            .get_long_help()
            .map(|h| h.to_string())
            .unwrap_or_default();
        if help.trim().is_empty() && long.trim().is_empty() {
            problems.push(format!("{id}: no help text"));
        }
        let has_default = !arg.get_default_values().is_empty();
        let unset = spec(id).map(|p| p.unset);
        if !has_default && unset == Some(Unset::ClapDefault) {
            problems.push(format!(
                "{id}: no clap default, so its row needs Auto, Off or Required"
            ));
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

#[test]
fn units_match_flag_names() {
    for p in PARAMS {
        let expected = if p.id.ends_with("_ppm") {
            Some(Unit::Ppm)
        } else if p.id.ends_with("_da") {
            Some(Unit::Da)
        } else if p.id.ends_with("_gb") {
            Some(Unit::Gb)
        } else {
            None
        };
        if let Some(unit) = expected {
            assert_eq!(p.unit, unit, "{}: the id names its unit", p.id);
        }
    }
}

#[test]
fn cli_only_matches_config_mapping() {
    let mut keys = BTreeSet::new();
    for p in PARAMS {
        match p.section.yaml_name() {
            None => assert_eq!(p.key, p.id, "{}: a CLI-only row keeps its id as key", p.id),
            Some(section) => assert!(
                keys.insert((section, p.key)),
                "{}: config key {section}.{} is used by two rows",
                p.id,
                p.key
            ),
        }
    }
    let cli_only: Vec<&str> = PARAMS
        .iter()
        .filter(|p| p.section == Section::CliOnly)
        .map(|p| p.id)
        .collect();
    assert!(
        cli_only.contains(&"config") && cli_only.contains(&"quiet"),
        "CliOnly rows: {cli_only:?}"
    );
}
