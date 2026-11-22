//! Keeps the generated files in step with the declaration: `config.example.yaml`, the DOCS.md
//! flag tables between their markers, and the `--print-config` text.

use clap::CommandFactory;

use crate::config::RunConfig;
use crate::param_sources::{ParamSources, Source};
use crate::params::{Section, PARAMS};
use crate::params_render::{config_example, docs_reference, effective_config};

const START: &str = "<!-- generated:search-flags:start -->";
const END: &str = "<!-- generated:search-flags:end -->";

fn repo_file(name: &str) -> String {
    format!("{}/../../{name}", env!("CARGO_MANIFEST_DIR"))
}

/// The byte range of DOCS.md between the marker lines.
fn generated_span(docs: &str) -> (usize, usize) {
    let start = docs.find(START).expect("DOCS.md has the start marker") + START.len() + 1;
    let end = docs.find(END).expect("DOCS.md has the end marker");
    assert!(start <= end, "the start marker comes before the end marker");
    (start, end)
}

/// Whether `text` names `flag` as a whole flag (`--xl` is not in `--xl-oracle`).
fn mentions(text: &str, flag: &str) -> bool {
    text.match_indices(flag).any(|(i, _)| {
        text[i + flag.len()..]
            .chars()
            .next()
            .is_none_or(|c| !(c.is_ascii_alphanumeric() || c == '-' || c == '_'))
    })
}

#[test]
fn config_example_is_generated() {
    let committed = std::fs::read_to_string(repo_file("config.example.yaml"))
        .expect("read config.example.yaml");
    assert!(
        committed == config_example(),
        "config.example.yaml is out of date; run `cargo test -p andes write_config_example -- --ignored`"
    );
}

#[test]
fn docs_reference_is_generated() {
    let docs = std::fs::read_to_string(repo_file("DOCS.md")).expect("read DOCS.md");
    let (start, end) = generated_span(&docs);
    assert!(
        docs[start..end] == docs_reference(),
        "DOCS.md flag tables are out of date; run `cargo test -p andes write_docs_reference -- --ignored`"
    );
}

#[test]
fn generated_example_config_parses() {
    RunConfig::from_yaml_str(&config_example()).expect("the generated example parses");
}

#[test]
fn internal_flags_are_not_advertised() {
    let mut cmd = crate::TopCli::command();
    cmd.build();
    let example = config_example();
    let docs = docs_reference();
    let mut leaks = Vec::new();
    for arg in cmd.get_arguments().filter(|a| a.is_hide_set()) {
        let Some(long) = arg.get_long() else { continue };
        let flag = format!("--{long}");
        if mentions(&example, &flag) || mentions(&docs, &flag) {
            leaks.push(flag.clone());
        }
        let key = PARAMS
            .iter()
            .find(|p| arg.get_id() == p.id && p.section != Section::CliOnly)
            .map(|p| p.key);
        if let Some(key) = key {
            let written = example.lines().any(|l| {
                let l = l.trim_start().trim_start_matches('#').trim_start();
                l.starts_with(&format!("{key}:"))
            });
            if written {
                leaks.push(format!("key {key} of {flag}"));
            }
        }
    }
    assert!(leaks.is_empty(), "hidden flags advertised: {leaks:?}");
}

#[test]
fn effective_config_parses() {
    let text = effective_config(&ParamSources::default()).expect("printable");
    RunConfig::from_yaml_str(&text).expect("the default effective config parses");
}

#[test]
fn effective_config_names_sources() {
    let mut sources = ParamSources::default();
    sources.record_as("psm_fdr", Source::CommandLine, &0.02);
    sources.record_as("max_spectra", Source::Config, &5u64);
    let text = effective_config(&sources).expect("printable");
    let line = |key: &str| text.lines().find(|l| l.trim_start().starts_with(key));
    assert!(line("psm_fdr: 0.02").is_some_and(|l| l.contains("# command line")));
    assert!(line("max_spectra: 5").is_some_and(|l| l.contains("# config")));
    assert!(line("# min_length: 6").is_some_and(|l| l.contains("# default")));
    assert!(line("# xl:").is_some_and(|l| l.contains("unset: off")));
    RunConfig::from_yaml_str(&text).expect("an effective config with set values parses");
}

#[test]
fn effective_config_leaves_defaults_unset() {
    let mut sources = ParamSources::default();
    sources.record_as("top_n", Source::CommandLine, &7u64);
    let text = effective_config(&sources).expect("printable");
    let cfg = RunConfig::from_yaml_str(&text).expect("parses");
    assert_eq!(cfg.search.max_missed_cleavages, None);
    assert_eq!(cfg.search.top_n, Some(7));
}

#[test]
fn example_config_sets_only_the_paths() {
    let cfg = RunConfig::from_yaml_str(&config_example()).expect("parses");
    assert_eq!(cfg.search.max_missed_cleavages, None);
}

#[test]
fn effective_config_hides_internal_flags_until_set() {
    let text = effective_config(&ParamSources::default()).expect("printable");
    assert!(!text.contains("max_spectra"), "internal key shown by default");
}

#[test]
fn printed_strings_read_back_unchanged() {
    let texts = [
        "line\nbreak",
        "tab\there",
        "quote \" and \\ backslash",
        "yes",
        "- dash",
        "#hash",
        "colon: space",
        "ünïcode",
        "bell\u{7}",
        "123",
        "",
    ];
    for text in texts {
        let mut sources = ParamSources::default();
        sources.record_as("percolator_image", Source::CommandLine, text);
        let printed = effective_config(&sources).expect("printable");
        let cfg = RunConfig::from_yaml_str(&printed)
            .unwrap_or_else(|e| panic!("{text:?}: {e}\n{printed}"));
        assert_eq!(
            cfg.rescoring.percolator_image.as_deref(),
            Some(text),
            "{printed}"
        );
    }
}

#[test]
fn printed_config_closes_with_the_command_line_only_flags() {
    let mut sources = ParamSources::default();
    sources.record_as("xl_oracle", Source::CommandLine, &5u64);
    sources.record_as("ipf_debug", Source::CommandLine, &true);
    sources.record_as("quiet", Source::CommandLine, &true);
    let text = effective_config(&sources).expect("printable");
    let note = "# set on the command line, no config key (pass again to reproduce): ";
    assert!(text.contains(&format!("{note}--xl-oracle 5\n")), "{text}");
    assert!(text.contains(&format!("{note}--ipf-debug\n")), "{text}");
    assert!(!text.contains("--quiet"), "{text}");
    RunConfig::from_yaml_str(&text).expect("the closing comments keep the YAML valid");
}

/// A value whose serialisation always fails, as a value YAML cannot hold.
struct Unwritable;

impl serde::Serialize for Unwritable {
    fn serialize<S: serde::Serializer>(&self, _: S) -> Result<S::Ok, S::Error> {
        Err(serde::ser::Error::custom("no YAML form"))
    }
}

#[test]
fn a_value_without_config_form_refuses_to_print() {
    let mut sources = ParamSources::default();
    sources.record_as("percolator_image", Source::CommandLine, &Unwritable);
    assert_eq!(sources.source("percolator_image"), Source::CommandLine);
    assert_eq!(sources.unwritable().len(), 1);
    assert!(effective_config(&sources).is_err());
}

#[test]
#[ignore = "writes config.example.yaml; run on the build host and commit the file"]
fn write_config_example() {
    std::fs::write(repo_file("config.example.yaml"), config_example())
        .expect("write config.example.yaml");
}

#[test]
#[ignore = "writes DOCS.md between its markers; run on the build host and commit the file"]
fn write_docs_reference() {
    let path = repo_file("DOCS.md");
    let docs = std::fs::read_to_string(&path).expect("read DOCS.md");
    let (start, end) = generated_span(&docs);
    let updated = format!("{}{}{}", &docs[..start], docs_reference(), &docs[end..]);
    std::fs::write(&path, updated).expect("write DOCS.md");
}
