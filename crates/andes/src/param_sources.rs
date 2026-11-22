//! Records where each search setting came from: its default, the `--config` file, the SDRF or
//! the command line.
//!
//! The merge in `main` starts from [`from_command_line`] and lets `config::apply` and
//! `sdrf::apply` add what they set. A setting nobody recorded is a [`Source::Default`]. Each
//! record carries the value as it was written: the raw command-line text, the YAML value of the
//! config key, or the SDRF cell. Later consumers (the run report, `--print-config`) read the
//! result through [`ParamSources::source`] and [`ParamSources::value`].

use std::collections::BTreeMap;

use clap::parser::ValueSource;
use serde_yaml::Value;

use crate::params::{spec, PARAMS};

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
/// keyed by clap id.
#[derive(Debug, Default, Clone)]
pub struct ParamSources(BTreeMap<&'static str, (Source, Option<Value>)>);

impl ParamSources {
    /// Notes that `id` was set from `source` to `value`. An id outside the parameter table is
    /// a programming error and is ignored in release builds.
    pub fn record(&mut self, id: &str, source: Source, value: Value) {
        match spec(id) {
            Some(p) => {
                self.0.insert(p.id, (source, Some(value)));
            }
            None => debug_assert!(false, "`{id}` is not a search parameter"),
        }
    }

    /// Like [`record`](Self::record), for a value that is not yet a YAML value (a list of
    /// strings, a path, a flag).
    pub fn record_as<T: serde::Serialize + ?Sized>(&mut self, id: &str, source: Source, v: &T) {
        self.record(id, source, serde_yaml::to_value(v).unwrap_or_default());
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
        self.0.get(id).map_or(Source::Default, |(s, _)| *s)
    }

    /// The value recorded for `id`, as written.
    pub fn value(&self, id: &str) -> Option<&Value> {
        self.0.get(id).and_then(|(_, v)| v.as_ref())
    }
}

/// Every search setting given on the command line, with its raw text: one string, or a list
/// when the flag was given more than one value. A flag without a value is `true`.
pub fn from_command_line(m: &clap::ArgMatches) -> ParamSources {
    let mut out = ParamSources::default();
    for p in PARAMS {
        let given = matches!(m.try_contains_id(p.id), Ok(true))
            && m.value_source(p.id) == Some(ValueSource::CommandLine);
        if !given {
            continue;
        }
        let mut raw: Vec<Value> = m
            .get_raw(p.id)
            .map(|vals| vals.map(|v| Value::from(v.to_string_lossy().into_owned())).collect())
            .unwrap_or_default();
        let value = match raw.len() {
            0 => Value::Bool(true),
            1 => raw.remove(0),
            _ => Value::Sequence(raw),
        };
        out.record(p.id, Source::CommandLine, value);
    }
    out
}
