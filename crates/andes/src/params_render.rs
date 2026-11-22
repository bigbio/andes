//! Renders the search parameters from their one declaration: `config.example.yaml`, the DOCS.md
//! flag tables and the `--print-config` YAML.
//!
//! The declaration is the search command's clap arguments (flag, default, help, visibility,
//! heading) joined with the `PARAMS` table (section, key, unit, range, unset meaning, legacy
//! form, example). Nothing here holds a flag name or a default of its own. A hidden argument is
//! internal and is never advertised; an argument hidden from the short help is "(advanced)".

use clap::{Arg, ArgAction, CommandFactory};
use serde_yaml::Value;

use crate::param_sources::{takes_several, typed_text, ParamSources, Source};
use crate::params::{ParamSpec, Section, Unit, Unset, PARAMS};

const EXAMPLE_HEADER: &str = r##"# andes run configuration: `andes --config this-file.yaml`.
#
# Every key is optional; a key left out keeps its default. A flag on the command line wins over
# this file (command line > config > default). Unknown keys are an error.
# Values are written as on the command line ("20ppm", "2..5", auto). The values below are the
# defaults, except the example paths; a commented key without a default shows an example value.
# Each comment names the key's flag; "advanced" flags are listed by `andes --help` only,

"##;

const EXAMPLE_TAIL: &str = r##"# Read by `andes cohort plan` and `andes cohort run`; a search ignores this section. The sections above are the
# settings every run of the cohort shares; each run's own settings come from its SDRF row.
# cohort:
#   name: human-xl                 # the default output folder, <name>/ next to this file
#   data_dir: /data/raw            # spectrum files, any layout; <data_dir>/<project>/ is used when present
#   output_dir: /data/human-xl     # runs/<run>/, groups/<group>/ and plan.json
#   sdrfs:                         # one entry per project; paths relative to this file
#     - {project: PXD019110, sdrf: PXD019110.sdrf.tsv}
#     - {project: PXD059766, sdrf: PXD059766.sdrf.tsv}
#   run_name: "{project}.{stem}"   # unique run name; {project} and {stem} are filled in
#   overrides:                     # replace or fill in what a project's SDRF gives
#     - project: PXD019110
#       search: {precursor_tol: 10ppm}   # keys: enzyme, precursor_tol, fixed_mods, variable_mods, xl
#   exclude: []                    # projects listed above but left out of the plan
"##;

/// The section banners of `config.example.yaml`: the sections a banner covers and its title.
const BANNERS: &[(&[Section], &str)] = &[
    (&[Section::Io], "Input and output"),
    (&[Section::Search], "Search"),
    (&[Section::Scoring], "Calibration, model and runtime"),
    (&[Section::Decoys], "Decoys"),
    (&[Section::Chimeric, Section::Refine], "Modes"),
    (&[Section::Rescoring], "Rescoring and FDR"),
    (&[Section::Glyco], "Glycopeptides"),
    (&[Section::Quant], "Quantification"),
];

/// Widest key text that the comments of a section line up behind.
const ALIGN_CAP: usize = 40;

/// One written setting: the key text (with its `# ` when commented), the trailing comment and
/// the comment block under it.
struct Line {
    key: String,
    note: String,
    block: Option<&'static str>,
}

fn command() -> clap::Command {
    let mut cmd = crate::TopCli::command();
    cmd.build();
    cmd
}

fn is_builtin(arg: &Arg) -> bool {
    matches!(arg.get_id().as_str(), "help" | "version")
}

fn flag(arg: &Arg) -> String {
    arg.get_long()
        .map_or_else(|| arg.get_id().to_string(), |l| format!("--{l}"))
}

fn collapsed(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// The first sentence of `help` without its closing period; an `e.g.` or `i.e.` does not end it.
fn first_sentence(help: &str) -> String {
    let text = collapsed(help);
    let bytes = text.as_bytes();
    for (i, _) in text.match_indices('.') {
        let ends = !matches!(bytes.get(i + 1), Some(b) if *b != b' ');
        let abbreviation = text[..i].ends_with("e.g") || text[..i].ends_with("i.e");
        if ends && !abbreviation {
            return text[..i].to_string();
        }
    }
    text
}

fn help_of(arg: &Arg) -> String {
    arg.get_help().map(ToString::to_string).unwrap_or_default()
}

fn is_list(arg: &Arg, p: &ParamSpec) -> bool {
    takes_several(arg) || matches!(p.id, "fixed_mods" | "variable_mods")
}

/// The config value of the argument's own text or texts: a number for a numeric unit, a
/// boolean for `true` and `false`, a list when the argument takes several values.
fn value_of(arg: &Arg, p: &ParamSpec, texts: &[String]) -> Value {
    let one = |t: &String| match t.as_str() {
        "true" => Value::Bool(true),
        "false" => Value::Bool(false),
        _ => typed_text(p.unit, t),
    };
    if is_list(arg, p) {
        Value::Sequence(texts.iter().map(one).collect())
    } else {
        texts.first().map_or(Value::Null, one)
    }
}

/// The clap default as a config value; `None` when clap holds none or hides it (a default that
/// depends on the machine). A switch defaults to off.
fn default_value(arg: &Arg, p: &ParamSpec) -> Option<Value> {
    match arg.get_action() {
        ArgAction::SetTrue => return Some(Value::Bool(false)),
        ArgAction::SetFalse => return Some(Value::Bool(true)),
        _ => {}
    }
    if arg.is_hide_default_value_set() {
        return None;
    }
    let mut texts: Vec<String> = Vec::new();
    for d in arg.get_default_values() {
        let d = d.to_string_lossy();
        match arg.get_value_delimiter() {
            Some(sep) => texts.extend(d.split(sep).map(str::to_string)),
            None => texts.push(d.into_owned()),
        }
    }
    (!texts.is_empty()).then(|| value_of(arg, p, &texts))
}

fn example_value(arg: &Arg, p: &ParamSpec) -> Option<Value> {
    p.example.map(|e| value_of(arg, p, &[e.to_string()]))
}

fn unset_text(unset: Unset) -> String {
    match unset {
        Unset::ClapDefault | Unset::Off => "off".to_string(),
        Unset::Auto(text) => text.to_string(),
        Unset::Required => "required".to_string(),
    }
}

fn scalar(s: &str) -> String {
    let word = |c: char| c.is_ascii_alphanumeric() || "_./+-".contains(c);
    let first_ok = s
        .chars()
        .next()
        .is_some_and(|c| c.is_ascii_alphanumeric() || c == '_' || c == '/' || c == '.');
    let reserved = matches!(
        s.to_ascii_lowercase().as_str(),
        "true" | "false" | "null" | "yes" | "no" | "on" | "off" | "y" | "n" | "~"
    );
    if first_ok && s.chars().all(word) && !reserved && s.parse::<f64>().is_err() {
        s.to_string()
    } else {
        format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
    }
}

/// One-line YAML of a value: scalars plain where that is safe, lists and maps in flow style.
fn flow(v: &Value) -> String {
    match v {
        Value::Null => "null".to_string(),
        Value::Bool(b) => b.to_string(),
        Value::Number(n) => n.to_string(),
        Value::String(s) => scalar(s),
        Value::Sequence(items) => {
            let items: Vec<String> = items.iter().map(flow).collect();
            format!("[{}]", items.join(", "))
        }
        Value::Mapping(m) => {
            let pairs: Vec<String> = m
                .iter()
                .map(|(k, v)| format!("{}: {}", flow(k), flow(v)))
                .collect();
            format!("{{{}}}", pairs.join(", "))
        }
        Value::Tagged(t) => flow(&t.value),
    }
}

/// The lines of one config section, comments lined up and comment blocks under their rows.
fn write_section(out: &mut String, name: &str, lines: &[Line]) {
    out.push_str(name);
    out.push_str(":\n");
    let width = lines
        .iter()
        .map(|l| l.key.len())
        .filter(|w| *w <= ALIGN_CAP)
        .max()
        .unwrap_or(0);
    for l in lines {
        let pad = width.saturating_sub(l.key.len()) + 1;
        out.push_str(&format!("  {}{}# {}\n", l.key, " ".repeat(pad), l.note));
        for b in l.block.into_iter().flat_map(str::lines) {
            if b.is_empty() {
                out.push_str("  #\n");
            } else {
                out.push_str(&format!("  # {b}\n"));
            }
        }
    }
}

/// The banner-and-section body of a config file; `line_for` turns one row into its line, or
/// `None` to leave it out.
fn write_sections(
    cmd: &clap::Command,
    line_for: &dyn Fn(&Arg, &ParamSpec) -> Option<Line>,
) -> String {
    let mut out = String::new();
    for (i, (sections, title)) in BANNERS.iter().enumerate() {
        if i > 0 {
            out.push('\n');
        }
        out.push_str(&format!("# ── {title} "));
        let used = out.lines().last().map_or(0, |l| l.chars().count());
        out.push_str(&"─".repeat(78usize.saturating_sub(used)));
        out.push('\n');
        for (j, section) in sections.iter().enumerate() {
            if j > 0 {
                out.push('\n');
            }
            let lines: Vec<Line> = PARAMS
                .iter()
                .filter(|p| p.section == *section)
                .filter_map(|p| {
                    let arg = cmd.get_arguments().find(|a| a.get_id() == p.id)?;
                    line_for(arg, p)
                })
                .collect();
            write_section(&mut out, section.yaml_name().unwrap_or_default(), &lines);
        }
    }
    out
}

/// The whole `config.example.yaml`: every public parameter on one line, set to its default, or
/// commented out with an example value when it has none. Only the tests that keep the committed
/// file in step call it.
#[allow(dead_code)]
pub fn config_example() -> String {
    let cmd = command();
    let body = write_sections(&cmd, &|arg, p| {
        if arg.is_hide_set() {
            return None;
        }
        let default = default_value(arg, p);
        let required = p.unset == Unset::Required;
        let (value, example, commented) = match (default, example_value(arg, p)) {
            (Some(d), _) => (flow(&d), false, false),
            (None, Some(e)) => (flow(&e), true, !required),
            (None, None) => ("...".to_string(), false, true),
        };
        let advanced = if arg.is_hide_short_help_set() {
            " (advanced)"
        } else {
            ""
        };
        let marker = if example { " (example)" } else { "" };
        Some(Line {
            key: format!("{}{}: {value}", if commented { "# " } else { "" }, p.key),
            note: format!(
                "{}{advanced}: {}{marker}",
                flag(arg),
                first_sentence(&help_of(arg))
            ),
            block: p.example_block,
        })
    });
    format!("{EXAMPLE_HEADER}{body}\n{EXAMPLE_TAIL}")
}

fn unit_text(unit: Unit) -> &'static str {
    match unit {
        Unit::None => "",
        Unit::Ppm => "ppm",
        Unit::Da => "Da",
        Unit::PpmOrDa => "ppm or Da",
        Unit::Minutes => "min",
        Unit::Seconds => "s",
        Unit::Gb => "GB",
        Unit::Fraction => "fraction",
        Unit::Count => "count",
        Unit::Number => "number",
        Unit::Charge => "charge",
        Unit::Path => "path",
    }
}

fn cell(text: &str) -> String {
    text.replace('|', "\\|")
}

fn code_cell(text: Option<&str>) -> String {
    text.map_or_else(|| "—".to_string(), |t| format!("`{}`", cell(t)))
}

fn default_cell(arg: &Arg, p: &ParamSpec) -> String {
    if matches!(arg.get_action(), ArgAction::SetTrue | ArgAction::SetFalse) {
        return "*(off)*".to_string();
    }
    let defaults: Vec<String> = arg
        .get_default_values()
        .iter()
        .map(|d| d.to_string_lossy().into_owned())
        .collect();
    if defaults.is_empty() || arg.is_hide_default_value_set() {
        format!("*({})*", cell(&unset_text(p.unset)))
    } else {
        format!("`{}`", cell(&defaults.join(",")))
    }
}

fn description(arg: &Arg) -> String {
    let long = arg.get_long_help().map(ToString::to_string);
    cell(&collapsed(&long.unwrap_or_else(|| help_of(arg))))
}

/// The DOCS.md flag tables: one per help heading, in declaration order, public flags only. Only
/// the tests that keep DOCS.md in step call it.
#[allow(dead_code)]
pub fn docs_reference() -> String {
    let cmd = command();
    let mut tables: Vec<(String, Vec<String>)> = Vec::new();
    for arg in cmd.get_arguments().filter(|a| !is_builtin(a) && !a.is_hide_set()) {
        let Some(p) = PARAMS.iter().find(|p| arg.get_id() == p.id) else {
            continue;
        };
        let heading = arg.get_help_heading().unwrap_or("Other").to_string();
        let advanced = if arg.is_hide_short_help_set() {
            " (advanced)"
        } else {
            ""
        };
        let unit = unit_text(p.unit);
        let row = format!(
            "| `{}`{advanced} | {} | {} | {} | {} | {} |",
            flag(arg),
            if unit.is_empty() { "—" } else { unit },
            default_cell(arg, p),
            code_cell(p.range),
            description(arg),
            code_cell(p.legacy),
        );
        match tables.iter_mut().find(|(h, _)| *h == heading) {
            Some((_, rows)) => rows.push(row),
            None => tables.push((heading, vec![row])),
        }
    }
    let rendered: Vec<String> = tables
        .iter()
        .map(|(heading, rows)| {
            format!(
                "### {heading}\n\n| Flag | Unit | Default | Range | Description | Legacy form |\n\
                 |---|---|---|---|---|---|\n{}\n",
                rows.join("\n")
            )
        })
        .collect();
    rendered.join("\n")
}

fn source_text(source: Source) -> &'static str {
    match source {
        Source::Default => "default",
        Source::Config => "config",
        Source::Sdrf => "SDRF",
        Source::CommandLine => "command line",
    }
}

/// The `--print-config` YAML: every config parameter under its section with the value the run
/// uses (the recorded value, else the clap default) and a comment naming its source. A
/// parameter with neither is commented out with what unset means. Internal parameters appear
/// only when something set them. `--print-config` (andes P1 task 4) is its reader.
#[allow(dead_code)]
pub fn effective_config(sources: &ParamSources) -> String {
    let cmd = command();
    let body = write_sections(&cmd, &|arg, p| {
        let source = sources.source(p.id);
        if arg.is_hide_set() && source == Source::Default {
            return None;
        }
        let value = sources
            .value(p.id)
            .cloned()
            .or_else(|| default_value(arg, p))
            .filter(|_| source == Source::Default || sources.value(p.id).is_some());
        let (key, note) = match value {
            Some(v) => (
                format!("{}: {}", p.key, flow(&v)),
                source_text(source).to_string(),
            ),
            None if source == Source::Default => (
                format!("# {}:", p.key),
                format!("unset: {}", unset_text(p.unset)),
            ),
            None => (
                format!("# {}:", p.key),
                format!("{}; no config form of the value", source_text(source)),
            ),
        };
        Some(Line {
            key,
            note,
            block: None,
        })
    });
    format!(
        "# The settings this run uses (`--print-config`); each comment names where the value came from.\n\n{body}"
    )
}
