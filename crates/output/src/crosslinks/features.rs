//! The cross-link feature table (`<base>.xl.features.parquet`): one row per quantified
//! cross-linked precursor and run, with the evidence it rests on, its MS1 quantity and the
//! quality of its feature.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, BooleanArray, Float64Array, Int32Array, StringArray, UInt32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::versioned::{read_table, write_table, ReadTable, TableVersion, BATCH_ROWS};

const TABLE: TableVersion = TableVersion {
    key: "andes.xl_features.schema",
    version: "1.0",
    major: 1,
    what: "cross-link feature table",
};

/// One cross-linked precursor in one run.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct XlFeatureRow {
    pub run: String,
    /// The precursor key: linker, both peptides with their modifications and link sites, charge.
    pub precursor: String,
    pub residue_pair: String,
    pub protein_pair: String,
    pub link: String,
    pub crosslinker: String,
    /// The bare sequence, the peptide with its modifications, the 1-based link site and the
    /// kind of site (`nterm`, `cterm` or the linked residue) of each peptide.
    pub alpha_peptide: String,
    pub alpha_peptidoform: String,
    pub alpha_site: u32,
    pub alpha_site_kind: String,
    pub beta_peptide: String,
    pub beta_peptidoform: String,
    pub beta_site: u32,
    pub beta_site_kind: String,
    pub charge: u32,
    pub mono_mz: f64,
    /// `identified` (a cross-link spectrum match in this run) or `transferred` (measured at
    /// the retention time the alignment predicts, without one).
    pub evidence: String,
    /// Whether the value is reported as a quantity: an identification whose isotope pattern
    /// fits, or a transfer at or below the transfer q-value threshold.
    pub accepted: bool,
    /// The MS1 area over the whole isotope envelope; null without a feature.
    pub quantity: Option<f64>,
    /// The q-value of the feature against the feature at its shifted decoy mass.
    pub quant_q: Option<f64>,
    pub feature_score: Option<f64>,
    /// Where the feature was looked for: the retention time of the identification, or the
    /// predicted one of a transfer with its standard deviation, in seconds.
    pub anchor_rt_s: f64,
    pub anchor_sd_s: f64,
    pub apex_rt_s: Option<f64>,
    pub rt_start_s: Option<f64>,
    pub rt_stop_s: Option<f64>,
    pub apex_intensity: Option<f64>,
    pub cosine: Option<f64>,
    pub n_isotopes: u32,
    pub n_isotopes_found: u32,
    pub scan: i32,
    pub csm_score: Option<f64>,
    pub n_csms: u32,
    pub residue_pair_q: f64,
}

fn schema() -> Schema {
    let text = |name: &str| Field::new(name, DataType::Utf8, false);
    let count = |name: &str| Field::new(name, DataType::UInt32, false);
    let number = |name: &str| Field::new(name, DataType::Float64, false);
    let maybe = |name: &str| Field::new(name, DataType::Float64, true);
    Schema::new(vec![
        text("run"),
        text("precursor"),
        text("residue_pair"),
        text("protein_pair"),
        text("link"),
        text("crosslinker"),
        text("alpha_peptide"),
        text("alpha_peptidoform"),
        count("alpha_site"),
        text("alpha_site_kind"),
        text("beta_peptide"),
        text("beta_peptidoform"),
        count("beta_site"),
        text("beta_site_kind"),
        count("charge"),
        number("mono_mz"),
        text("evidence"),
        Field::new("accepted", DataType::Boolean, false),
        maybe("quantity"),
        maybe("quant_q"),
        maybe("feature_score"),
        number("anchor_rt_s"),
        number("anchor_sd_s"),
        maybe("apex_rt_s"),
        maybe("rt_start_s"),
        maybe("rt_stop_s"),
        maybe("apex_intensity"),
        maybe("cosine"),
        count("n_isotopes"),
        count("n_isotopes_found"),
        Field::new("scan", DataType::Int32, false),
        maybe("csm_score"),
        count("n_csms"),
        number("residue_pair_q"),
    ])
}

fn batch(schema: &Arc<Schema>, rows: &[XlFeatureRow]) -> Result<RecordBatch, String> {
    let text = |f: &dyn Fn(&XlFeatureRow) -> &str| -> ArrayRef {
        Arc::new(StringArray::from_iter_values(rows.iter().map(f)))
    };
    let count = |f: &dyn Fn(&XlFeatureRow) -> u32| -> ArrayRef {
        Arc::new(UInt32Array::from_iter_values(rows.iter().map(f)))
    };
    let number = |f: &dyn Fn(&XlFeatureRow) -> f64| -> ArrayRef {
        Arc::new(Float64Array::from_iter_values(rows.iter().map(f)))
    };
    let maybe = |f: &dyn Fn(&XlFeatureRow) -> Option<f64>| -> ArrayRef {
        Arc::new(rows.iter().map(f).collect::<Float64Array>())
    };
    let cols: Vec<ArrayRef> = vec![
        text(&|r| &r.run),
        text(&|r| &r.precursor),
        text(&|r| &r.residue_pair),
        text(&|r| &r.protein_pair),
        text(&|r| &r.link),
        text(&|r| &r.crosslinker),
        text(&|r| &r.alpha_peptide),
        text(&|r| &r.alpha_peptidoform),
        count(&|r| r.alpha_site),
        text(&|r| &r.alpha_site_kind),
        text(&|r| &r.beta_peptide),
        text(&|r| &r.beta_peptidoform),
        count(&|r| r.beta_site),
        text(&|r| &r.beta_site_kind),
        count(&|r| r.charge),
        number(&|r| r.mono_mz),
        text(&|r| &r.evidence),
        Arc::new(BooleanArray::from_iter(
            rows.iter().map(|r| Some(r.accepted)),
        )),
        maybe(&|r| r.quantity),
        maybe(&|r| r.quant_q),
        maybe(&|r| r.feature_score),
        number(&|r| r.anchor_rt_s),
        number(&|r| r.anchor_sd_s),
        maybe(&|r| r.apex_rt_s),
        maybe(&|r| r.rt_start_s),
        maybe(&|r| r.rt_stop_s),
        maybe(&|r| r.apex_intensity),
        maybe(&|r| r.cosine),
        count(&|r| r.n_isotopes),
        count(&|r| r.n_isotopes_found),
        Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.scan))),
        maybe(&|r| r.csm_score),
        count(&|r| r.n_csms),
        number(&|r| r.residue_pair_q),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

pub fn write_xl_features(path: &Path, rows: &[XlFeatureRow]) -> Result<(), String> {
    let schema = Arc::new(schema());
    let batches = rows.chunks(BATCH_ROWS).map(|c| batch(&schema, c));
    write_table(path, &TABLE, TABLE.version, schema.clone(), batches, &[])
}

pub fn read_xl_features(path: &Path) -> Result<ReadTable, String> {
    read_table(path, &TABLE)
}

#[cfg(test)]
mod tests {
    use arrow::array::Array;

    use super::*;
    use crate::versioned::column;

    #[test]
    fn rows_round_trip_with_nulls_and_the_schema_version() {
        let t = tempfile::tempdir().unwrap();
        let p = t.path().join("c.xl.features.parquet");
        let found = XlFeatureRow {
            run: "run1".into(),
            precursor: "BS3|AAK@3|CCK@3/4".into(),
            residue_pair: "P1:10 -- P1:40".into(),
            evidence: "identified".into(),
            accepted: true,
            quantity: Some(1.5e7),
            quant_q: Some(0.0),
            charge: 4,
            n_isotopes: 7,
            ..XlFeatureRow::default()
        };
        let missing = XlFeatureRow {
            run: "run2".into(),
            evidence: "transferred".into(),
            ..found.clone()
        };
        let missing = XlFeatureRow {
            accepted: false,
            quantity: None,
            quant_q: None,
            ..missing
        };
        write_xl_features(&p, &[found, missing]).unwrap();
        let back = read_xl_features(&p).unwrap();
        assert_eq!(back.meta.get(TABLE.key).map(String::as_str), Some("1.0"));
        let b = &back.batches[0];
        assert_eq!(b.num_rows(), 2);
        let evidence = column::<StringArray>(b, "evidence", &p).unwrap();
        let quantity = column::<Float64Array>(b, "quantity", &p).unwrap();
        let accepted = column::<BooleanArray>(b, "accepted", &p).unwrap();
        assert_eq!(
            (evidence.value(0), evidence.value(1)),
            ("identified", "transferred")
        );
        assert_eq!(quantity.value(0), 1.5e7);
        assert!(quantity.is_null(1));
        assert!(accepted.value(0) && !accepted.value(1));
    }
}
