//! The run report: the inputs with their checksums, every parameter with the source of its value,
//! the time and peak memory of each stage, the counts at each FDR level, the software, models and
//! decoys, and a Methods paragraph with its references. The search appends its text to
//! `<stem>.statistics.log` and writes the same rows to `<stem>.statistics.parquet`.
//!
//! The Parquet twin (`andes.run_report.schema`) is one long table with the columns `section`,
//! `key`, `value`, `source` (nullable), `number` (float64, nullable) and `fdr` (float64,
//! nullable). Seconds, bytes, counts and FDR levels are numbers, so no reader parses text.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::statistics::thousands;
use crate::versioned::{write_table, TableVersion};

/// The major version of the run-report table; readers reject any other.
pub const RUN_REPORT_SCHEMA_MAJOR: u32 = 1;

const RUN_REPORT: TableVersion = TableVersion {
    key: "andes.run_report.schema",
    version: "1.0",
    major: RUN_REPORT_SCHEMA_MAJOR,
    what: "run report",
};

/// The widest line of the text sections; longer text is wrapped at word boundaries.
const LINE_WIDTH_CHARS: usize = 100;

/// Bytes in the unit the report prints memory in (GiB, printed as `GB` like the run log).
const BYTES_PER_GB: f64 = (1u64 << 30) as f64;

/// One input file: its role (`spectrum`, `database`, `sdrf`, `config`), its path, its size and
/// its SHA-256, or the error that kept it from being read.
#[derive(Debug, Clone, PartialEq)]
pub struct InputFile {
    pub role: String,
    pub path: String,
    pub bytes: Option<u64>,
    pub sha256: Result<String, String>,
}

/// One parameter: its id, its value as one line of YAML and where the value came from
/// (`default`, `config`, `SDRF` or `command line`). `number` is the value when it is a number.
#[derive(Debug, Clone, PartialEq)]
pub struct ParamRow {
    pub id: String,
    pub value: String,
    pub source: String,
    pub number: Option<f64>,
}

/// One stage of the run: its wall time and the peak memory of the process when it ended.
#[derive(Debug, Clone, PartialEq)]
pub struct StageRow {
    pub name: String,
    pub seconds: f64,
    pub peak_rss_bytes: Option<u64>,
}

/// One count: what is counted (`spectra`, `PSMs`, `peptides`, `protein groups`, ...), at which
/// q-value threshold when it has one, and how many.
#[derive(Debug, Clone, PartialEq)]
pub struct CountRow {
    pub what: String,
    pub fdr: Option<f64>,
    pub n: usize,
}

/// One reference of the Methods paragraph: the registry key and the text
/// `Authors (year). Title. Journal. doi:...`.
#[derive(Debug, Clone, PartialEq)]
pub struct Reference {
    pub key: String,
    pub text: String,
}

/// Everything the run report holds, as plain data; `render` and `write_parquet` turn it into
/// the text sections and the Parquet twin.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RunReport {
    pub andes_version: String,
    pub command_line: String,
    pub inputs: Vec<InputFile>,
    pub params: Vec<ParamRow>,
    pub stages: Vec<StageRow>,
    pub counts: Vec<CountRow>,
    /// `(role, model id)` of each model the run used.
    pub models: Vec<(String, String)>,
    /// How the decoys were made, such as `shuffled sequences (seed 42), prefix XXX_`.
    pub decoys: String,
    /// The decoy key of the searched database (`andes database`'s key of the same settings).
    pub decoy_key: String,
    pub methods: String,
    pub references: Vec<Reference>,
    /// The file name of the Parquet twin, named in the text where it lists all parameters.
    pub parquet_name: String,
}

/// A section header in the style of the run summary.
fn header(out: &mut String, title: &str) {
    out.push_str(&format!("──────── {title} ────────\n"));
}

/// `text` wrapped at word boundaries to lines of at most `LINE_WIDTH_CHARS`, each starting with
/// `indent`; a word longer than a line stays whole.
fn wrapped(out: &mut String, indent: &str, text: &str) {
    let mut line = String::new();
    for word in text.split_whitespace() {
        let width = indent.chars().count() + line.chars().count() + 1 + word.chars().count();
        if !line.is_empty() && width > LINE_WIDTH_CHARS {
            out.push_str(&format!("{indent}{line}\n"));
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    if !line.is_empty() {
        out.push_str(&format!("{indent}{line}\n"));
    }
}

fn gb(bytes: u64) -> String {
    format!("{:.2} GB", bytes as f64 / BYTES_PER_GB)
}

fn checksum_text(sha256: &Result<String, String>) -> String {
    match sha256 {
        Ok(hex) => format!("sha256 {hex}"),
        Err(e) => format!("checksum unavailable: {e}"),
    }
}

fn fdr_text(fdr: Option<f64>) -> String {
    fdr.map_or_else(String::new, |q| format!("q<={q}"))
}

impl RunReport {
    /// The text sections appended to `statistics.log`, in order: Inputs, Parameters (only those
    /// not at their default), Stages, Counts, Software, Methods, References.
    pub fn render(&self) -> String {
        let mut s = String::new();
        self.render_inputs(&mut s);
        self.render_params(&mut s);
        self.render_stages(&mut s);
        self.render_counts(&mut s);
        self.render_software(&mut s);
        header(&mut s, "Methods");
        wrapped(&mut s, "  ", &self.methods);
        header(&mut s, "References");
        for (i, r) in self.references.iter().enumerate() {
            wrapped(&mut s, "  ", &format!("[{}] {}", i + 1, r.text));
        }
        s
    }

    fn render_inputs(&self, s: &mut String) {
        header(s, "Inputs");
        let role_w = self.inputs.iter().map(|f| f.role.len()).max().unwrap_or(0);
        for f in &self.inputs {
            let size = f.bytes.map_or_else(
                || "size unknown".to_string(),
                |b| format!("{} bytes", thousands(b as usize)),
            );
            s.push_str(&format!("  {:<role_w$}  {}  ({size})\n", f.role, f.path));
            s.push_str(&format!("  {:<role_w$}  {}\n", "", checksum_text(&f.sha256)));
        }
    }

    fn render_params(&self, s: &mut String) {
        header(s, "Parameters");
        let set: Vec<&ParamRow> = self.params.iter().filter(|p| p.source != "default").collect();
        let id_w = set.iter().map(|p| p.id.len()).max().unwrap_or(0);
        for p in &set {
            s.push_str(&format!("  {:<id_w$}  {}  ({})\n", p.id, p.value, p.source));
        }
        s.push_str(&format!(
            "  ({} parameters at their default; all are in {})\n",
            self.params.len() - set.len(),
            self.parquet_name
        ));
    }

    fn render_stages(&self, s: &mut String) {
        header(s, "Stages");
        let name_w = self.stages.iter().map(|r| r.name.len()).max().unwrap_or(0).max("stage".len());
        s.push_str(&format!(
            "  {:<name_w$}  {:>10}  peak memory (process, so far)\n",
            "stage", "wall"
        ));
        for r in &self.stages {
            let peak = r.peak_rss_bytes.map_or_else(|| "-".to_string(), gb);
            s.push_str(&format!(
                "  {:<name_w$}  {:>8.2} s  {peak}\n",
                r.name, r.seconds
            ));
        }
    }

    fn render_counts(&self, s: &mut String) {
        header(s, "Counts");
        let what_w = self.counts.iter().map(|c| c.what.len()).max().unwrap_or(0);
        for c in &self.counts {
            s.push_str(&format!(
                "  {:<what_w$}  {:<9}  {:>12}\n",
                c.what,
                fdr_text(c.fdr),
                thousands(c.n)
            ));
        }
    }

    fn render_software(&self, s: &mut String) {
        header(s, "Software");
        s.push_str(&format!("  andes      {}\n", self.andes_version));
        wrapped(s, "  command    ", &self.command_line);
        for (role, id) in &self.models {
            s.push_str(&format!("  model      {role}: {id}\n"));
        }
        s.push_str(&format!("  decoys     {}\n", self.decoys));
        s.push_str(&format!("  decoy key  {}\n", self.decoy_key));
    }

    /// The rows of the Parquet twin, in the order of the text sections: every parameter,
    /// including those at their default.
    fn rows(&self) -> Vec<Row> {
        let mut rows = Vec::new();
        for f in &self.inputs {
            rows.push(Row::new("input", &f.path, &checksum_text(&f.sha256))
                .source(&f.role)
                .number(f.bytes.map(|b| b as f64)));
        }
        for p in &self.params {
            rows.push(Row::new("parameter", &p.id, &p.value).source(&p.source).number(p.number));
        }
        for r in &self.stages {
            rows.push(Row::new("stage_seconds", &r.name, &format!("{:.2} s", r.seconds))
                .number(Some(r.seconds)));
            let peak = r.peak_rss_bytes.map_or_else(|| "-".to_string(), gb);
            rows.push(Row::new("stage_peak_rss_bytes", &r.name, &peak)
                .number(r.peak_rss_bytes.map(|b| b as f64)));
        }
        for c in &self.counts {
            let mut row = Row::new("count", &c.what, &thousands(c.n)).number(Some(c.n as f64));
            row.fdr = c.fdr;
            rows.push(row);
        }
        rows.push(Row::new("software", "andes_version", &self.andes_version));
        rows.push(Row::new("software", "command_line", &self.command_line));
        for (role, id) in &self.models {
            rows.push(Row::new("model", role, id));
        }
        rows.push(Row::new("software", "decoys", &self.decoys));
        rows.push(Row::new("software", "decoy_key", &self.decoy_key));
        rows.push(Row::new("methods", "paragraph", &self.methods));
        for (i, r) in self.references.iter().enumerate() {
            rows.push(Row::new("reference", &r.key, &r.text).number(Some((i + 1) as f64)));
        }
        rows
    }

    /// Writes the Parquet twin to `path` through a partial file, so a reader never sees half a
    /// table.
    pub fn write_parquet(&self, path: &Path) -> std::io::Result<()> {
        let rows = self.rows();
        let text = |f: fn(&Row) -> &str| -> ArrayRef {
            Arc::new(StringArray::from_iter_values(rows.iter().map(f)))
        };
        let schema = Arc::new(Schema::new(vec![
            Field::new("section", DataType::Utf8, false),
            Field::new("key", DataType::Utf8, false),
            Field::new("value", DataType::Utf8, false),
            Field::new("source", DataType::Utf8, true),
            Field::new("number", DataType::Float64, true),
            Field::new("fdr", DataType::Float64, true),
        ]));
        let columns: Vec<ArrayRef> = vec![
            text(|r| r.section.as_str()),
            text(|r| r.key.as_str()),
            text(|r| r.value.as_str()),
            Arc::new(StringArray::from(
                rows.iter().map(|r| r.source.as_deref()).collect::<Vec<_>>(),
            )),
            Arc::new(Float64Array::from(rows.iter().map(|r| r.number).collect::<Vec<_>>())),
            Arc::new(Float64Array::from(rows.iter().map(|r| r.fdr).collect::<Vec<_>>())),
        ];
        let batch = RecordBatch::try_new(schema.clone(), columns).map_err(|e| e.to_string());
        write_table(path, &RUN_REPORT, RUN_REPORT.version, schema, [batch], &[])
            .map_err(std::io::Error::other)
    }
}

/// One row of the Parquet twin.
struct Row {
    section: String,
    key: String,
    value: String,
    source: Option<String>,
    number: Option<f64>,
    fdr: Option<f64>,
}

impl Row {
    fn new(section: &str, key: &str, value: &str) -> Self {
        Self {
            section: section.to_string(),
            key: key.to_string(),
            value: value.to_string(),
            source: None,
            number: None,
            fdr: None,
        }
    }

    fn source(mut self, source: &str) -> Self {
        self.source = Some(source.to_string());
        self
    }

    fn number(mut self, number: Option<f64>) -> Self {
        self.number = number;
        self
    }
}

#[cfg(test)]
#[path = "report/tests/report.rs"]
mod tests;
