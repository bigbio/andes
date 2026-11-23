//! Records where each search setting came from: its default, the `--config` file, the SDRF or
//! the command line.
//!
//! The merge in `main` starts from [`from_command_line`] and lets `config::apply` and
//! `sdrf::apply` add what they set. A setting nobody recorded is a [`Source::Default`]. Each
//! record carries the value as a config value, in the shape the `RunConfig` key expects: the
//! command-line text typed by the argument, the YAML value of the config key, or the SDRF cell
//! turned into that shape. Later consumers (the run report, `--print-config`) read the
//! result through [`ParamSources::source`] and [`ParamSources::value`]. A value that has no
//! config form (text that is not UTF-8, a value YAML cannot hold) is noted in
//! [`ParamSources::unwritable`], so that `--print-config` refuses to write a config that would
//! not reproduce the run.

use std::collections::BTreeMap;

use clap::parser::ValueSource;
use clap::{ArgAction, CommandFactory};
use serde_yaml::Value;

use crate::params::{spec, Unit, PARAMS};

/// Where a setting's value came from. The command line beats the config file, which beats the
/// SDRF, which beats the default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Default,
    Config,
    Sdrf,
    CommandLine,
}

/// The source and written value of every search setting that was not left at its default,
/// keyed by clap id, and the settings whose value has no config form.
#[derive(Debug, Default, Clone)]
pub struct ParamSources {
    set: BTreeMap<&'static str, (Source, Option<Value>)>,
    unwritable: Vec<String>,
}

impl ParamSources {
    /// Notes that `id` was set from `source` to `value`. An id outside the parameter table is
    /// a programming error and is ignored in release builds.
    pub fn record(&mut self, id: &str, source: Source, value: Value) {
        match spec(id) {
            Some(p) => {
                self.set.insert(p.id, (source, Some(value)));
            }
            None => debug_assert!(false, "`{id}` is not a search parameter"),
        }
    }

    /// Like [`record`](Self::record), for a value that is not yet a YAML value (a list of
    /// strings, a path, a flag). A value YAML cannot hold records the source alone and is
    /// noted as unwritable.
    pub fn record_as<T: serde::Serialize + ?Sized>(&mut self, id: &str, source: Source, v: &T) {
        match serde_yaml::to_value(v) {
            Ok(value) => self.record(id, source, value),
            Err(e) => {
                self.record_without_value(id, source);
                self.unwritable.push(format!("`{id}`: {e}"));
            }
        }
    }

    /// Notes that `id` was set from `source` and that no config form of its value is known.
    /// An id outside the parameter table is a programming error and is ignored in release
    /// builds.
    pub fn record_without_value(&mut self, id: &str, source: Source) {
        match spec(id) {
            Some(p) => {
                self.set.insert(p.id, (source, None));
            }
            None => debug_assert!(false, "`{id}` is not a search parameter"),
        }
    }

    /// The settings whose value was given but has no config form, each with the reason.
    pub fn unwritable(&self) -> &[String] {
        &self.unwritable
    }

    /// Notes that `id` was set by the config file whose parsed text is `written`, with the value
    /// of its section and key in that file.
    pub fn record_config(&mut self, id: &str, written: &Value) {
        let value = spec(id)
            .and_then(|p| written.get(p.section.yaml_name()?)?.get(p.key))
            .cloned()
            .unwrap_or_default();
        self.record(id, Source::Config, value);
    }

    /// The source of `id`; [`Source::Default`] when nothing was recorded.
    pub fn source(&self, id: &str) -> Source {
        self.set.get(id).map_or(Source::Default, |(s, _)| *s)
    }

    /// The value recorded for `id` as a config value; `None` when the source set the setting but
    /// no config form of its value is known (an SDRF crosslinker cell).
    pub fn value(&self, id: &str) -> Option<&Value> {
        self.set.get(id).and_then(|(_, v)| v.as_ref())
    }
}

/// Every search setting given on the command line, as the config value that would set it: a
/// list for an argument that takes several values (and for the modification and crosslinker
/// lists, split as their parsers split them), a boolean for a switch, a number for a numeric
/// unit when the text is one, otherwise the raw text. Text that is not UTF-8 records the source
/// alone and is noted as unwritable.
pub fn from_command_line(m: &clap::ArgMatches) -> ParamSources {
    let mut cmd = crate::TopCli::command();
    cmd.build();
    let mut out = ParamSources::default();
    for p in PARAMS {
        let Some(arg) = cmd.get_arguments().find(|a| a.get_id() == p.id) else {
            debug_assert!(false, "`{}` has a row but no clap argument", p.id);
            continue;
        };
        let given = matches!(m.try_contains_id(p.id), Ok(true))
            && m.value_source(p.id) == Some(ValueSource::CommandLine);
        if !given {
            continue;
        }
        let value = if let Ok(Some(flag)) = m.try_get_one::<bool>(p.id) {
            Value::Bool(*flag)
        } else {
            let delimiter = arg.get_value_delimiter();
            let Some(raw) = utf8_values(m, p.id) else {
                out.record_without_value(p.id, Source::CommandLine);
                out.unwritable
                    .push(format!("`{}`: the command-line value is not UTF-8 text", p.id));
                continue;
            };
            let raw: Vec<String> = match delimiter {
                Some(d) => raw.iter().flat_map(|r| r.split(d)).map(str::to_string).collect(),
                None => raw,
            };
            let typed = |t: &String| typed_text(p.unit, t);
            if matches!(p.id, "fixed_mods" | "variable_mods") {
                let joined = raw.join(", ");
                let items = model::mod_spec::split_list(&joined)
                    .map(|v| v.into_iter().map(|i| Value::from(i.trim())).collect::<Vec<_>>());
                Value::Sequence(items.unwrap_or_else(|_| raw.iter().map(typed).collect()))
            } else if p.id == "xl" {
                let items = raw.iter().flat_map(|r| crate::cli::parse::split_xl(r));
                Value::Sequence(items.map(Value::from).collect())
            } else if takes_several(arg) {
                Value::Sequence(raw.iter().map(typed).collect())
            } else {
                raw.first().map_or(Value::Bool(true), typed)
            }
        };
        out.record(p.id, Source::CommandLine, value);
    }
    out
}

/// The raw command-line texts of `id`; `None` when one of them is not UTF-8.
fn utf8_values(m: &clap::ArgMatches, id: &str) -> Option<Vec<String>> {
    m.get_raw(id).map_or(Some(Vec::new()), |vals| {
        vals.map(|v| v.to_str().map(str::to_string)).collect()
    })
}

/// Whether clap collects more than one value for the argument.
pub(crate) fn takes_several(arg: &clap::Arg) -> bool {
    matches!(arg.get_action(), ArgAction::Append)
        || arg.get_value_delimiter().is_some()
        || arg.get_num_args().is_some_and(|n| n.max_values() > 1)
}

/// A YAML number when the unit is numeric (a charge range such as `2..5` stays text) and the
/// text is a number, otherwise the text. Integers are tried as `u64`, then `i64`, then the
/// text as a finite `f64`, so every integer the CLI accepts stays an integer.
pub(crate) fn typed_text(unit: Unit, text: &str) -> Value {
    let numeric = matches!(
        unit,
        Unit::Count
            | Unit::Number
            | Unit::Fraction
            | Unit::Gb
            | Unit::Seconds
            | Unit::Ppm
            | Unit::Da
            | Unit::Charge
    );
    let number = || {
        if let Ok(u) = text.parse::<u64>() {
            return Some(Value::from(u));
        }
        if let Ok(i) = text.parse::<i64>() {
            return Some(Value::from(i));
        }
        text.parse::<f64>().ok().filter(|f| f.is_finite()).map(Value::from)
    };
    if numeric {
        number().unwrap_or_else(|| Value::from(text))
    } else {
        Value::from(text)
    }
}
