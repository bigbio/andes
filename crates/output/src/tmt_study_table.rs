//! The study-level isobaric table of `andes quant finalize` on a TMT study,
//! `<base>.tmt.parquet` (schema key `andes.tmt_study.schema`).
//!
//! One row per (run, PSM) at the pooled PSM q threshold, with the reporter channels as
//! columns, as in the per-run `<stem>.tmt.parquet`:
//! - `filename`, `plex`, `fraction`, `instrument`: the run and its SDRF design (null without
//!   an SDRF);
//! - the PSM: `scannr`, `spec_id`, `quant_scan`, `quant_ms_level`, `rt`, `peptide`
//!   (peptidoform), `sequence`, `charge`, `proteins` (the search's, `;`-separated);
//! - `psm_q_value`, `psm_pep`, `psm_score`: the pooled rescoring of every run;
//! - the reporters: `<prefix>_<channel>` (corrected when the runs were corrected),
//!   `raw_<prefix>_<channel>` (with correction), `reporter_intensity_sum`, `reporter_sn_sum`,
//!   `sn_<prefix>_<channel>`, and with unmixing `interference` and `unmix_<prefix>_<channel>`;
//! - `purity`, `quant_gate_passed` (the run's gate: purity and summed S/N), `quant_weight`
//!   (purity squared times 1 - pooled PEP when the gate passes, else 0);
//! - the study's protein group of the peptide at the protein q threshold: `pg_accessions`
//!   (`;`-separated), `anchor_protein`, `unique`, `pg_q_value` (null when the peptide has no
//!   group at the threshold);
//! - `sample_<prefix>_<channel>` (SDRF source name) and `role_<prefix>_<channel>` (`sample`,
//!   `reference` or `blank`) of the row's run; null for a channel the SDRF does not name.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, Float32Builder, Float64Builder, Int32Builder, StringBuilder,
    UInt32Builder, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use quant::isobaric::Plex;

use crate::versioned_parquet::{read_table, ReadTable, TableVersion, TableWriter, BATCH_ROWS};

const TABLE: TableVersion = TableVersion {
    key: "andes.tmt_study.schema",
    version: "1.0",
    major: 1,
    what: "isobaric study table",
};
/// Metadata key of the schema version; a per-run table has `andes.tmt.schema` instead.
pub const TMT_STUDY_SCHEMA_KEY: &str = TABLE.key;
pub const TMT_STUDY_PSM_FDR_KEY: &str = "andes.tmt_study.psm_fdr";
pub const TMT_STUDY_PROTEIN_FDR_KEY: &str = "andes.tmt_study.protein_fdr";
/// The `--pep` threshold of the rows, or `none`.
pub const TMT_STUDY_PEP_KEY: &str = "andes.tmt_study.pep";
pub const TMT_STUDY_RUNS_KEY: &str = "andes.tmt_study.runs";
pub const TMT_STUDY_SDRF_KEY: &str = "andes.tmt_study.sdrf";

/// One (run, PSM) row.
#[derive(Debug, Clone, PartialEq)]
pub struct TmtStudyRow {
    pub filename: String,
    pub scannr: i32,
    pub spec_id: String,
    pub quant_scan: String,
    pub quant_ms_level: u8,
    pub rt: Option<f64>,
    pub peptide: String,
    pub sequence: String,
    pub charge: u8,
    pub proteins: String,
    pub psm_q_value: f64,
    pub psm_pep: f64,
    pub psm_score: Option<f64>,
    pub purity: Option<f32>,
    pub intensities: Vec<f32>,
    pub raw: Option<Vec<f32>>,
    pub intensity_sum: f32,
    pub sn_sum: Option<f32>,
    pub sn: Vec<Option<f32>>,
    pub interference: Option<f32>,
    pub unmixed: Vec<Option<f32>>,
    pub quant_gate_passed: bool,
    pub quant_weight: f32,
    pub group: Option<StudyGroup>,
}

/// The protein group of a row's peptide.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyGroup {
    pub accessions: String,
    pub anchor: String,
    pub unique: bool,
    pub q_value: f64,
}

/// The SDRF design of one run, shared by its rows.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RunDesignColumns {
    pub plex: Option<u32>,
    pub fraction: Option<String>,
    pub instrument: Option<String>,
    /// Per plex channel, the sample (SDRF source name) and its role.
    pub channels: Vec<Option<(String, &'static str)>>,
}

/// Layout of the reporter columns, the same for every run of a study.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TmtStudyLayout {
    pub plex: Plex,
    pub corrected: bool,
    pub unmixed: bool,
}

fn prefix(plex: Plex) -> &'static str {
    if matches!(plex, Plex::Itraq4 | Plex::Itraq8) {
        "itraq"
    } else {
        "tmt"
    }
}

fn schema(l: TmtStudyLayout) -> Schema {
    let p = prefix(l.plex);
    let names: Vec<&str> = l.plex.channels().iter().map(|c| c.name).collect();
    let block = |head: &str, nullable: bool| -> Vec<Field> {
        names
            .iter()
            .map(|c| Field::new(format!("{head}{p}_{c}"), DataType::Float32, nullable))
            .collect()
    };
    let mut f = vec![
        Field::new("filename", DataType::Utf8, false),
        Field::new("plex", DataType::UInt32, true),
        Field::new("fraction", DataType::Utf8, true),
        Field::new("instrument", DataType::Utf8, true),
        Field::new("scannr", DataType::Int32, false),
        Field::new("spec_id", DataType::Utf8, false),
        Field::new("quant_scan", DataType::Utf8, false),
        Field::new("quant_ms_level", DataType::UInt8, false),
        Field::new("rt", DataType::Float64, true),
        Field::new("peptide", DataType::Utf8, false),
        Field::new("sequence", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("proteins", DataType::Utf8, false),
        Field::new("psm_q_value", DataType::Float64, false),
        Field::new("psm_pep", DataType::Float64, false),
        Field::new("psm_score", DataType::Float64, true),
        Field::new("purity", DataType::Float32, true),
    ];
    f.extend(block("", false));
    if l.corrected {
        f.extend(block("raw_", false));
    }
    f.push(Field::new(
        "reporter_intensity_sum",
        DataType::Float32,
        false,
    ));
    f.push(Field::new("reporter_sn_sum", DataType::Float32, true));
    f.extend(block("sn_", true));
    if l.unmixed {
        f.push(Field::new("interference", DataType::Float32, true));
        f.extend(block("unmix_", true));
    }
    f.extend([
        Field::new("quant_gate_passed", DataType::Boolean, false),
        Field::new("quant_weight", DataType::Float32, false),
        Field::new("pg_accessions", DataType::Utf8, true),
        Field::new("anchor_protein", DataType::Utf8, true),
        Field::new("unique", DataType::Boolean, true),
        Field::new("pg_q_value", DataType::Float64, true),
    ]);
    for c in &names {
        f.push(Field::new(format!("sample_{p}_{c}"), DataType::Utf8, true));
    }
    for c in &names {
        f.push(Field::new(format!("role_{p}_{c}"), DataType::Utf8, true));
    }
    Schema::new(f)
}

/// Streams the study table run by run.
pub struct TmtStudyWriter {
    layout: TmtStudyLayout,
    schema: Arc<Schema>,
    writer: TableWriter,
    rows: usize,
}

impl TmtStudyWriter {
    pub fn create(
        path: &Path,
        layout: TmtStudyLayout,
        meta: &[(&str, String)],
    ) -> Result<Self, String> {
        let schema = Arc::new(schema(layout));
        let writer = TableWriter::create(
            path,
            &TABLE,
            TABLE.version,
            schema.clone(),
            meta,
            BATCH_ROWS,
        )?;
        Ok(Self {
            layout,
            schema,
            writer,
            rows: 0,
        })
    }

    /// Writes the rows of one run with the run's design.
    pub fn push_run(
        &mut self,
        rows: &[TmtStudyRow],
        design: &RunDesignColumns,
    ) -> Result<(), String> {
        for chunk in rows.chunks(BATCH_ROWS) {
            let b = self.batch(chunk, design)?;
            self.writer.write(&b)?;
            self.rows += chunk.len();
        }
        Ok(())
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn finish(self) -> Result<(), String> {
        self.writer.finish()
    }

    fn batch(&self, rows: &[TmtStudyRow], d: &RunDesignColumns) -> Result<RecordBatch, String> {
        let n_ch = self.layout.plex.n_channels();
        if d.channels.len() != n_ch {
            return Err(format!(
                "internal: {} channel samples for {n_ch} channels",
                d.channels.len()
            ));
        }
        let mut filename = StringBuilder::new();
        let mut plex = UInt32Builder::new();
        let mut fraction = StringBuilder::new();
        let mut instrument = StringBuilder::new();
        let mut scannr = Int32Builder::new();
        let mut spec_id = StringBuilder::new();
        let mut quant_scan = StringBuilder::new();
        let mut level = UInt8Builder::new();
        let mut rt = Float64Builder::new();
        let mut peptide = StringBuilder::new();
        let mut sequence = StringBuilder::new();
        let mut charge = UInt8Builder::new();
        let mut proteins = StringBuilder::new();
        let mut q = Float64Builder::new();
        let mut pep = Float64Builder::new();
        let mut score = Float64Builder::new();
        let mut purity = Float32Builder::new();
        let new_block = || (0..n_ch).map(|_| Float32Builder::new()).collect::<Vec<_>>();
        let (mut shown, mut raw, mut sn, mut unmix) =
            (new_block(), new_block(), new_block(), new_block());
        let mut sum = Float32Builder::new();
        let mut sn_sum = Float32Builder::new();
        let mut interference = Float32Builder::new();
        let mut gate = BooleanBuilder::new();
        let mut weight = Float32Builder::new();
        let mut pg = StringBuilder::new();
        let mut anchor = StringBuilder::new();
        let mut unique = BooleanBuilder::new();
        let mut pg_q = Float64Builder::new();
        let mut samples: Vec<StringBuilder> = (0..n_ch).map(|_| StringBuilder::new()).collect();
        let mut roles: Vec<StringBuilder> = (0..n_ch).map(|_| StringBuilder::new()).collect();
        for r in rows {
            let bad = r.intensities.len() != n_ch
                || r.sn.len() != n_ch
                || (self.layout.corrected && r.raw.as_ref().is_none_or(|v| v.len() != n_ch))
                || (self.layout.unmixed && r.unmixed.len() != n_ch);
            if bad {
                return Err(format!(
                    "{}: reporter values do not match the {n_ch} channels",
                    r.spec_id
                ));
            }
            filename.append_value(&r.filename);
            plex.append_option(d.plex);
            fraction.append_option(d.fraction.as_deref());
            instrument.append_option(d.instrument.as_deref());
            scannr.append_value(r.scannr);
            spec_id.append_value(&r.spec_id);
            quant_scan.append_value(&r.quant_scan);
            level.append_value(r.quant_ms_level);
            rt.append_option(r.rt);
            peptide.append_value(&r.peptide);
            sequence.append_value(&r.sequence);
            charge.append_value(r.charge);
            proteins.append_value(&r.proteins);
            q.append_value(r.psm_q_value);
            pep.append_value(r.psm_pep);
            score.append_option(r.psm_score);
            purity.append_option(r.purity);
            for (b, &v) in shown.iter_mut().zip(&r.intensities) {
                b.append_value(v);
            }
            if let Some(v) = &r.raw {
                for (b, &v) in raw.iter_mut().zip(v) {
                    b.append_value(v);
                }
            }
            sum.append_value(r.intensity_sum);
            sn_sum.append_option(r.sn_sum);
            for (b, &v) in sn.iter_mut().zip(&r.sn) {
                b.append_option(v);
            }
            interference.append_option(r.interference);
            for (b, &v) in unmix.iter_mut().zip(&r.unmixed) {
                b.append_option(v);
            }
            gate.append_value(r.quant_gate_passed);
            weight.append_value(r.quant_weight);
            let g = r.group.as_ref();
            pg.append_option(g.map(|g| g.accessions.as_str()));
            anchor.append_option(g.map(|g| g.anchor.as_str()));
            unique.append_option(g.map(|g| g.unique));
            pg_q.append_option(g.map(|g| g.q_value));
            for (i, c) in d.channels.iter().enumerate() {
                samples[i].append_option(c.as_ref().map(|c| c.0.as_str()));
                roles[i].append_option(c.as_ref().map(|c| c.1));
            }
        }
        let fin = |b: &mut Vec<Float32Builder>| -> Vec<ArrayRef> {
            b.iter_mut()
                .map(|b| Arc::new(b.finish()) as ArrayRef)
                .collect()
        };
        let mut cols: Vec<ArrayRef> = vec![
            Arc::new(filename.finish()),
            Arc::new(plex.finish()),
            Arc::new(fraction.finish()),
            Arc::new(instrument.finish()),
            Arc::new(scannr.finish()),
            Arc::new(spec_id.finish()),
            Arc::new(quant_scan.finish()),
            Arc::new(level.finish()),
            Arc::new(rt.finish()),
            Arc::new(peptide.finish()),
            Arc::new(sequence.finish()),
            Arc::new(charge.finish()),
            Arc::new(proteins.finish()),
            Arc::new(q.finish()),
            Arc::new(pep.finish()),
            Arc::new(score.finish()),
            Arc::new(purity.finish()),
        ];
        cols.extend(fin(&mut shown));
        if self.layout.corrected {
            cols.extend(fin(&mut raw));
        }
        cols.push(Arc::new(sum.finish()));
        cols.push(Arc::new(sn_sum.finish()));
        cols.extend(fin(&mut sn));
        if self.layout.unmixed {
            cols.push(Arc::new(interference.finish()));
            cols.extend(fin(&mut unmix));
        }
        cols.push(Arc::new(gate.finish()));
        cols.push(Arc::new(weight.finish()));
        cols.push(Arc::new(pg.finish()));
        cols.push(Arc::new(anchor.finish()));
        cols.push(Arc::new(unique.finish()));
        cols.push(Arc::new(pg_q.finish()));
        cols.extend(samples.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
        cols.extend(roles.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
        RecordBatch::try_new(self.schema.clone(), cols).map_err(|e| e.to_string())
    }
}

/// Reads a study table, checking its schema version.
pub fn read_tmt_study(path: &Path) -> Result<ReadTable, String> {
    read_table(path, &TABLE)
}

/// Whether the Parquet file's footer names the study table's schema key (a per-run
/// `<stem>.tmt.parquet` does not).
pub fn is_tmt_study_table(path: &Path) -> Result<bool, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let reader = parquet::file::reader::SerializedFileReader::new(file)
        .map_err(|e| format!("{}: {e}", path.display()))?;
    use parquet::file::reader::FileReader;
    Ok(reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .is_some_and(|kv| kv.iter().any(|k| k.key == TABLE.key)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Array, StringArray, UInt32Array};

    fn row(spec: &str, gate: bool) -> TmtStudyRow {
        TmtStudyRow {
            filename: "r1".into(),
            scannr: 7,
            spec_id: spec.into(),
            quant_scan: "ms3".into(),
            quant_ms_level: 3,
            rt: Some(60.0),
            peptide: "PEPTIDEK".into(),
            sequence: "PEPTIDEK".into(),
            charge: 2,
            proteins: "P1".into(),
            psm_q_value: 0.001,
            psm_pep: 0.01,
            psm_score: Some(3.0),
            purity: Some(0.9),
            intensities: vec![1.0; 6],
            raw: None,
            intensity_sum: 6.0,
            sn_sum: None,
            sn: vec![None; 6],
            interference: None,
            unmixed: Vec::new(),
            quant_gate_passed: gate,
            quant_weight: if gate { 0.81 } else { 0.0 },
            group: gate.then(|| StudyGroup {
                accessions: "P1".into(),
                anchor: "P1".into(),
                unique: true,
                q_value: 0.0,
            }),
        }
    }

    #[test]
    fn rows_carry_the_run_design_and_the_schema_key() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("s.tmt.parquet");
        let layout = TmtStudyLayout {
            plex: Plex::Tmt6,
            corrected: false,
            unmixed: false,
        };
        let mut w = TmtStudyWriter::create(&path, layout, &[]).unwrap();
        let mut channels = vec![None; 6];
        channels[0] = Some(("GIS".to_string(), "reference"));
        let design = RunDesignColumns {
            plex: Some(2),
            fraction: Some("3".into()),
            instrument: None,
            channels,
        };
        w.push_run(&[row("r1/a_7_1", true), row("r1/b_8_1", false)], &design)
            .unwrap();
        assert_eq!(w.rows(), 2);
        w.finish().unwrap();
        assert!(is_tmt_study_table(&path).unwrap());
        let t = read_tmt_study(&path).unwrap();
        let b = &t.batches[0];
        let col = |n: &str| b.column(b.schema().index_of(n).unwrap()).clone();
        let plex = col("plex");
        assert_eq!(
            plex.as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .value(1),
            2
        );
        let s = col("sample_tmt_126");
        let s = s.as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(s.value(0), "GIS");
        assert_eq!(col("role_tmt_127").null_count(), 2);
        assert_eq!(col("pg_accessions").null_count(), 1);
    }
}
