//! Quantification tables `<stem>.tmt.parquet` and `<stem>.lfq.parquet`, and
//! the record types the QPX feature parquet is written from.
//!
//! `tmt.parquet` has one row per PSM with reporter ions and one column per
//! channel; `lfq.parquet` one row per quantified feature per run. Both carry
//! their schema version in the key-value metadata (`andes.tmt.schema`,
//! `andes.lfq.schema`; a reader rejects another major). Missing values are
//! null.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, Float32Builder, Float64Builder, Int32Builder, ListBuilder,
    StringBuilder, UInt32Builder, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::versioned_parquet::{read_table, write_table, ReadTable, TableVersion, BATCH_ROWS};

use quant::isobaric::Plex;
use quant::lfq::FeatureQuant;

/// One modification on a quantified peptidoform (for the parquet).
#[derive(Debug, Clone, PartialEq)]
pub struct ModRecord {
    pub name: String,
    pub accession: Option<String>,
    /// 1-based residue position.
    pub position: i32,
    pub amino_acid: char,
}

/// Identification side of a quantified row, shared by both modules.
#[derive(Debug, Clone)]
pub struct QuantId {
    /// Index of the input file (`--spectrum` order).
    pub file_idx: usize,
    /// File stem of the run (`run_file_name`).
    pub run: String,
    pub spec_id: String,
    pub scan: i32,
    pub rt_seconds: Option<f64>,
    pub sequence: String,
    pub peptidoform: String,
    pub modifications: Vec<ModRecord>,
    pub charge: u8,
    pub proteins: Vec<String>,
    pub is_decoy: bool,
    pub calculated_mz: f64,
    pub observed_mz: f64,
    pub q_value: Option<f64>,
    pub pep: Option<f64>,
    pub missed_cleavages: Option<i16>,
    /// 0-based index of the PSM row in `psms.parquet` (the feature's `psm_ids`).
    pub psm_index: i64,
}

/// One quantified PSM of an isobaric run.
#[derive(Debug, Clone)]
pub struct TmtRow {
    pub id: QuantId,
    /// Native id of the scan the reporters were read from (the MS2 itself, or
    /// its MS3).
    pub quant_scan_id: String,
    pub quant_ms_level: u8,
    pub purity: Option<f32>,
    pub raw: Vec<f32>,
    pub corrected: Vec<f32>,
    /// Per-channel signal-to-noise of the corrected intensities (the raw ones
    /// without a correction matrix), when the reporter scan records noise.
    pub sn: Option<Vec<f32>>,
    /// Suggested weight of this PSM in a protein or peptide roll-up (0 for a
    /// PSM the quantification gate removed).
    pub quant_weight: f32,
    /// The PSM passed the quantification gate.
    pub quant_gate_passed: bool,
    /// `--tmt-unmix`: the scan's estimated co-isolation interference
    /// fraction and the channel intensities with it removed (`None` when the
    /// scan could not be modelled).
    pub interference: Option<f32>,
    pub unmixed: Option<Vec<f32>>,
}

impl TmtRow {
    /// Summed raw reporter intensity.
    pub fn intensity_sum(&self) -> f32 {
        self.raw.iter().map(|v| v.max(0.0)).sum()
    }

    /// Summed reporter signal-to-noise, when noise is known.
    pub fn sn_sum(&self) -> Option<f32> {
        self.sn.as_ref().map(|sn| sn.iter().sum())
    }
}

/// One quantified precursor feature of a label-free run.
#[derive(Debug, Clone)]
pub struct LfqRow {
    pub id: QuantId,
    pub feature: FeatureQuant,
    pub score: f32,
    pub feature_q_value: Option<f64>,
    /// The decoy twin's score (None when the twin found no peak).
    pub decoy_score: Option<f32>,
}

const TMT_TABLE: TableVersion = TableVersion {
    key: "andes.tmt.schema",
    version: "1.1",
    major: 1,
    what: "isobaric quantification table",
};
const LFQ_TABLE: TableVersion = TableVersion {
    key: "andes.lfq.schema",
    version: "1.0",
    major: 1,
    what: "label-free quantification table",
};
/// `tmt.parquet` metadata: the plex name and whether a correction matrix was
/// applied (`true` / `false`).
pub const TMT_PLEX_KEY: &str = "andes.tmt.plex";
pub const TMT_CORRECTED_KEY: &str = "andes.tmt.corrected";
/// `tmt.parquet` metadata: whether the interference-unmixing block was
/// written (`true` / `false`).
pub const TMT_UNMIX_KEY: &str = "andes.tmt.unmix";
/// `lfq.parquet` metadata: the run stems in input order (newline-joined),
/// and the feature q-value and cosine cuts of the `accepted` column.
pub const LFQ_RUNS_KEY: &str = "andes.lfq.runs";
pub const LFQ_FEATURE_FDR_KEY: &str = "andes.lfq.feature_fdr";
pub const LFQ_MIN_COSINE_KEY: &str = "andes.lfq.min_cosine";

/// The column-name prefix of a plex's channels (`tmt` or `itraq`).
fn channel_prefix(plex: Plex) -> &'static str {
    if matches!(plex, Plex::Itraq4 | Plex::Itraq8) {
        "itraq"
    } else {
        "tmt"
    }
}

fn tmt_schema(plex: Plex, corrected: bool, unmixed: bool) -> Schema {
    let prefix = channel_prefix(plex);
    let mut fields = vec![
        Field::new("filename", DataType::Utf8, false),
        Field::new("scannr", DataType::Int32, false),
        Field::new("spec_id", DataType::Utf8, false),
        Field::new("quant_scan", DataType::Utf8, false),
        Field::new("quant_ms_level", DataType::UInt8, false),
        Field::new("rt", DataType::Float64, true),
        Field::new("peptide", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("proteins", DataType::Utf8, false),
        Field::new("is_decoy", DataType::Boolean, false),
        Field::new("psm_q_value", DataType::Float64, true),
        Field::new("psm_pep", DataType::Float64, true),
        Field::new("purity", DataType::Float32, true),
    ];
    for c in plex.channels() {
        fields.push(Field::new(
            format!("{prefix}_{}", c.name),
            DataType::Float32,
            false,
        ));
    }
    if corrected {
        for c in plex.channels() {
            fields.push(Field::new(
                format!("raw_{prefix}_{}", c.name),
                DataType::Float32,
                false,
            ));
        }
    }
    fields.extend([
        Field::new("reporter_intensity_sum", DataType::Float32, false),
        Field::new("reporter_sn_sum", DataType::Float32, true),
        Field::new("quant_weight", DataType::Float32, false),
        Field::new("quant_gate_passed", DataType::Boolean, false),
    ]);
    for c in plex.channels() {
        fields.push(Field::new(
            format!("sn_{prefix}_{}", c.name),
            DataType::Float32,
            true,
        ));
    }
    if unmixed {
        fields.push(Field::new("interference", DataType::Float32, true));
        for c in plex.channels() {
            fields.push(Field::new(
                format!("unmix_{prefix}_{}", c.name),
                DataType::Float32,
                true,
            ));
        }
    }
    Schema::new(fields)
}

fn tmt_batch(
    schema: &Arc<Schema>,
    n_channels: usize,
    corrected: bool,
    unmixed: bool,
    rows: &[TmtRow],
) -> Result<RecordBatch, String> {
    let mut filename = StringBuilder::new();
    let mut scannr = Int32Builder::new();
    let mut spec_id = StringBuilder::new();
    let mut quant_scan = StringBuilder::new();
    let mut level = UInt8Builder::new();
    let mut rt = Float64Builder::new();
    let mut peptide = StringBuilder::new();
    let mut charge = UInt8Builder::new();
    let mut proteins = StringBuilder::new();
    let mut decoy = BooleanBuilder::new();
    let mut q = Float64Builder::new();
    let mut pep = Float64Builder::new();
    let mut purity = Float32Builder::new();
    let mut shown: Vec<Float32Builder> = (0..n_channels).map(|_| Float32Builder::new()).collect();
    let mut raw: Vec<Float32Builder> = (0..n_channels).map(|_| Float32Builder::new()).collect();
    let mut sum = Float32Builder::new();
    let mut sn_sum = Float32Builder::new();
    let mut weight = Float32Builder::new();
    let mut passed = BooleanBuilder::new();
    let mut sn: Vec<Float32Builder> = (0..n_channels).map(|_| Float32Builder::new()).collect();
    let mut interference = Float32Builder::new();
    let mut unmix: Vec<Float32Builder> = (0..n_channels).map(|_| Float32Builder::new()).collect();
    for r in rows {
        let id = &r.id;
        if r.raw.len() != n_channels
            || r.corrected.len() != n_channels
            || r.sn.as_ref().is_some_and(|v| v.len() != n_channels)
            || r.unmixed.as_ref().is_some_and(|v| v.len() != n_channels)
        {
            return Err(format!(
                "{}: reporter values do not match the {n_channels} channels",
                id.spec_id
            ));
        }
        filename.append_value(&id.run);
        scannr.append_value(id.scan);
        spec_id.append_value(&id.spec_id);
        quant_scan.append_value(&r.quant_scan_id);
        level.append_value(r.quant_ms_level);
        rt.append_option(id.rt_seconds);
        peptide.append_value(&id.peptidoform);
        charge.append_value(id.charge);
        proteins.append_value(id.proteins.join(";"));
        decoy.append_value(id.is_decoy);
        q.append_option(id.q_value);
        pep.append_option(id.pep);
        purity.append_option(r.purity);
        let values = if corrected { &r.corrected } else { &r.raw };
        for (b, &v) in shown.iter_mut().zip(values) {
            b.append_value(v);
        }
        for (b, &v) in raw.iter_mut().zip(&r.raw) {
            b.append_value(v);
        }
        sum.append_value(r.intensity_sum());
        sn_sum.append_option(r.sn_sum());
        weight.append_value(r.quant_weight);
        passed.append_value(r.quant_gate_passed);
        for (i, b) in sn.iter_mut().enumerate() {
            b.append_option(r.sn.as_ref().map(|v| v[i]));
        }
        interference.append_option(r.interference);
        for (i, b) in unmix.iter_mut().enumerate() {
            b.append_option(r.unmixed.as_ref().map(|v| v[i]));
        }
    }
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(filename.finish()),
        Arc::new(scannr.finish()),
        Arc::new(spec_id.finish()),
        Arc::new(quant_scan.finish()),
        Arc::new(level.finish()),
        Arc::new(rt.finish()),
        Arc::new(peptide.finish()),
        Arc::new(charge.finish()),
        Arc::new(proteins.finish()),
        Arc::new(decoy.finish()),
        Arc::new(q.finish()),
        Arc::new(pep.finish()),
        Arc::new(purity.finish()),
    ];
    cols.extend(shown.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
    if corrected {
        cols.extend(raw.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
    }
    cols.push(Arc::new(sum.finish()));
    cols.push(Arc::new(sn_sum.finish()));
    cols.push(Arc::new(weight.finish()));
    cols.push(Arc::new(passed.finish()));
    cols.extend(sn.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
    if unmixed {
        cols.push(Arc::new(interference.finish()));
        cols.extend(unmix.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
    }
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

/// Write `tmt.parquet`: one row per PSM with reporter ions, including the
/// ones the quantification gate removed (`quant_gate_passed` false,
/// `quant_weight` 0). Columns: the PSM (`filename scannr spec_id quant_scan
/// quant_ms_level rt peptide charge proteins is_decoy psm_q_value psm_pep
/// purity`), the channel intensities in kit order (`tmt_126 … tmt_131`,
/// corrected when a correction matrix was applied), the raw values in a
/// `raw_` block when one was, then the summed raw reporter intensity, summed
/// signal-to-noise, roll-up weight and gate result, and a per-channel `sn_`
/// block (null when the reporter scan carries no noise). With `unmixed`
/// (`--tmt-unmix`), the estimated interference fraction and an `unmix_`
/// channel block follow (null for a scan that could not be modelled). `rt` is
/// in seconds, `proteins` is `;`-joined.
pub fn write_tmt_parquet(
    path: &Path,
    plex: Plex,
    rows: &[TmtRow],
    corrected: bool,
    unmixed: bool,
) -> Result<(), String> {
    let schema = Arc::new(tmt_schema(plex, corrected, unmixed));
    let n = plex.n_channels();
    let batches = rows
        .chunks(BATCH_ROWS)
        .map(|c| tmt_batch(&schema, n, corrected, unmixed, c));
    write_table(
        path,
        &TMT_TABLE,
        TMT_TABLE.version,
        schema.clone(),
        batches,
        &[
            (TMT_PLEX_KEY, plex.name().to_string()),
            (TMT_CORRECTED_KEY, corrected.to_string()),
            (TMT_UNMIX_KEY, unmixed.to_string()),
        ],
    )
}

/// Read `tmt.parquet`: its record batches and key-value metadata. Rejects
/// another major schema version; every error names the file.
pub fn read_tmt_parquet(path: &Path) -> Result<ReadTable, String> {
    read_table(path, &TMT_TABLE)
}

fn lfq_schema() -> Schema {
    Schema::new(vec![
        Field::new("filename", DataType::Utf8, false),
        Field::new("peptide", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("proteins", DataType::Utf8, false),
        Field::new("scannr", DataType::Int32, false),
        Field::new("spec_id", DataType::Utf8, false),
        Field::new("psm_q_value", DataType::Float64, true),
        Field::new("calc_mz", DataType::Float64, false),
        Field::new("rt_anchor", DataType::Float64, true),
        Field::new("rt_apex", DataType::Float64, false),
        Field::new("rt_start", DataType::Float64, false),
        Field::new("rt_stop", DataType::Float64, false),
        Field::new("n_scans", DataType::UInt32, false),
        Field::new("area", DataType::Float64, false),
        Field::new("apex_intensity", DataType::Float32, false),
        Field::new(
            "isotope_areas",
            DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
            false,
        ),
        Field::new("cosine", DataType::Float32, false),
        Field::new("n_isotopes", DataType::UInt32, false),
        Field::new("score", DataType::Float32, false),
        Field::new("decoy_score", DataType::Float32, true),
        Field::new("feature_q_value", DataType::Float64, true),
        Field::new("accepted", DataType::Boolean, false),
    ])
}

fn lfq_batch(
    schema: &Arc<Schema>,
    rows: &[LfqRow],
    accepted: &[bool],
) -> Result<RecordBatch, String> {
    let mut filename = StringBuilder::new();
    let mut peptide = StringBuilder::new();
    let mut charge = UInt8Builder::new();
    let mut proteins = StringBuilder::new();
    let mut scannr = Int32Builder::new();
    let mut spec_id = StringBuilder::new();
    let mut psm_q = Float64Builder::new();
    let mut calc_mz = Float64Builder::new();
    let mut rt_anchor = Float64Builder::new();
    let mut rt_apex = Float64Builder::new();
    let mut rt_start = Float64Builder::new();
    let mut rt_stop = Float64Builder::new();
    let mut n_scans = UInt32Builder::new();
    let mut area = Float64Builder::new();
    let mut apex_intensity = Float32Builder::new();
    let mut iso = ListBuilder::new(Float64Builder::new());
    let mut cosine = Float32Builder::new();
    let mut n_iso = UInt32Builder::new();
    let mut score = Float32Builder::new();
    let mut decoy_score = Float32Builder::new();
    let mut feature_q = Float64Builder::new();
    let mut acc = BooleanBuilder::new();
    for (r, &ok) in rows.iter().zip(accepted) {
        let f = &r.feature;
        filename.append_value(&r.id.run);
        peptide.append_value(&r.id.peptidoform);
        charge.append_value(r.id.charge);
        proteins.append_value(r.id.proteins.join(";"));
        scannr.append_value(r.id.scan);
        spec_id.append_value(&r.id.spec_id);
        psm_q.append_option(r.id.q_value);
        calc_mz.append_value(r.id.calculated_mz);
        rt_anchor.append_option(r.id.rt_seconds);
        rt_apex.append_value(f.apex_rt);
        rt_start.append_value(f.rt_start);
        rt_stop.append_value(f.rt_stop);
        n_scans.append_value(f.n_scans);
        area.append_value(f.area);
        apex_intensity.append_value(f.apex_intensity);
        iso.values().append_slice(&f.isotope_areas);
        iso.append(true);
        cosine.append_value(f.cosine);
        n_iso.append_value(f.n_isotopes_found);
        score.append_value(r.score);
        decoy_score.append_option(r.decoy_score);
        feature_q.append_option(r.feature_q_value);
        acc.append_value(ok);
    }
    let cols: Vec<ArrayRef> = vec![
        Arc::new(filename.finish()),
        Arc::new(peptide.finish()),
        Arc::new(charge.finish()),
        Arc::new(proteins.finish()),
        Arc::new(scannr.finish()),
        Arc::new(spec_id.finish()),
        Arc::new(psm_q.finish()),
        Arc::new(calc_mz.finish()),
        Arc::new(rt_anchor.finish()),
        Arc::new(rt_apex.finish()),
        Arc::new(rt_start.finish()),
        Arc::new(rt_stop.finish()),
        Arc::new(n_scans.finish()),
        Arc::new(area.finish()),
        Arc::new(apex_intensity.finish()),
        Arc::new(iso.finish()),
        Arc::new(cosine.finish()),
        Arc::new(n_iso.finish()),
        Arc::new(score.finish()),
        Arc::new(decoy_score.finish()),
        Arc::new(feature_q.finish()),
        Arc::new(acc.finish()),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

/// What `lfq.parquet` records beside its rows.
#[derive(Debug, Clone)]
pub struct LfqTableMeta<'a> {
    /// Run stems in input order, including runs without a feature.
    pub runs: &'a [String],
    /// Feature q-value cut of the `accepted` column.
    pub feature_fdr: f64,
    /// Envelope cosine cut of the `accepted` column.
    pub min_cosine: f32,
}

/// Write `lfq.parquet` (long): one row per quantified feature in one run
/// with its PSM anchor, apex, boundaries, areas (`isotope_areas` a list,
/// monoisotope first), the decoy-competition result and `accepted`
/// (`accepted[i]` for `rows[i]`: feature q-value and cosine within the run's
/// cuts). The per-run wide view is the `accepted` rows pivoted on
/// `(peptide, charge)` × `filename` over `area`, keeping the highest `score`
/// when a run has more than one.
pub fn write_lfq_parquet(
    path: &Path,
    rows: &[LfqRow],
    accepted: &[bool],
    meta: &LfqTableMeta<'_>,
) -> Result<(), String> {
    if accepted.len() != rows.len() {
        return Err(format!(
            "{}: {} accepted flags for {} rows",
            path.display(),
            accepted.len(),
            rows.len()
        ));
    }
    let schema = Arc::new(lfq_schema());
    let batches = rows
        .chunks(BATCH_ROWS)
        .zip(accepted.chunks(BATCH_ROWS))
        .map(|(r, a)| lfq_batch(&schema, r, a));
    write_table(
        path,
        &LFQ_TABLE,
        LFQ_TABLE.version,
        schema.clone(),
        batches,
        &[
            (LFQ_RUNS_KEY, meta.runs.join("\n")),
            (LFQ_FEATURE_FDR_KEY, meta.feature_fdr.to_string()),
            (LFQ_MIN_COSINE_KEY, meta.min_cosine.to_string()),
        ],
    )
}

/// Read `lfq.parquet`: its record batches and key-value metadata. Rejects
/// another major schema version; every error names the file.
pub fn read_lfq_parquet(path: &Path) -> Result<ReadTable, String> {
    read_table(path, &LFQ_TABLE)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::versioned_parquet::{column, opt};
    use arrow::array::{
        Array, BooleanArray, Float32Array, Float64Array, Int32Array, ListArray, StringArray,
        UInt32Array,
    };

    fn id(run: &str, file_idx: usize, pep: &str, z: u8) -> QuantId {
        QuantId {
            file_idx,
            run: run.to_string(),
            spec_id: format!("{run}.scan=5"),
            scan: 5,
            rt_seconds: Some(100.0),
            sequence: pep.to_string(),
            peptidoform: pep.to_string(),
            modifications: vec![],
            charge: z,
            proteins: vec!["P1".into()],
            is_decoy: false,
            calculated_mz: 500.0,
            observed_mz: 500.001,
            q_value: Some(0.001),
            pep: None,
            missed_cleavages: Some(0),
            psm_index: 0,
        }
    }

    fn feature(area: f64) -> FeatureQuant {
        FeatureQuant {
            apex_rt: 101.0,
            rt_start: 95.0,
            rt_stop: 107.0,
            area,
            apex_intensity: 1000.0,
            isotope_areas: vec![area * 0.6, area * 0.4],
            cosine: 0.98,
            n_isotopes_found: 2,
            rt_delta_s: 1.0,
            n_scans: 13,
        }
    }

    fn tmt_row() -> TmtRow {
        TmtRow {
            id: id("run1", 0, "PEPTIDEK", 2),
            quant_scan_id: "scan=5".into(),
            quant_ms_level: 2,
            purity: Some(0.9),
            raw: vec![1.0, 2.0, 3.0, 4.0],
            corrected: vec![1.1, 2.1, 3.1, 4.1],
            sn: Some(vec![11.0, 21.0, 31.0, 41.0]),
            quant_weight: 0.81,
            quant_gate_passed: true,
            interference: None,
            unmixed: None,
        }
    }

    fn names(t: &ReadTable) -> Vec<String> {
        t.batches[0]
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect()
    }

    fn f32_at(b: &RecordBatch, name: &str, path: &Path) -> Vec<Option<f32>> {
        let a = column::<Float32Array>(b, name, path).unwrap();
        (0..a.len()).map(|i| opt(a, i)).collect()
    }

    #[test]
    fn tmt_parquet_has_one_column_per_channel_raw_block_and_gated_rows() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.tmt.parquet");
        let mut gated = tmt_row();
        gated.sn = None;
        gated.purity = None;
        gated.id.q_value = None;
        gated.quant_weight = 0.0;
        gated.quant_gate_passed = false;
        write_tmt_parquet(&p, Plex::Itraq4, &[tmt_row(), gated], true, false).unwrap();
        let t = read_tmt_parquet(&p).unwrap();
        assert_eq!(
            t.meta.get(TMT_PLEX_KEY).map(String::as_str),
            Some("itraq4plex")
        );
        assert_eq!(
            t.meta.get(TMT_CORRECTED_KEY).map(String::as_str),
            Some("true")
        );
        let header = names(&t);
        assert_eq!(header.len(), 13 + 8 + 4 + 4);
        assert_eq!(header[13], "itraq_114");
        assert_eq!(header[17], "raw_itraq_114");
        assert_eq!(
            header[21..25],
            [
                "reporter_intensity_sum",
                "reporter_sn_sum",
                "quant_weight",
                "quant_gate_passed"
            ]
        );
        assert_eq!(header[25], "sn_itraq_114");
        let b = &t.batches[0];
        assert_eq!(b.num_rows(), 2);
        let pep = column::<StringArray>(b, "peptide", &p).unwrap();
        assert_eq!(pep.value(0), "PEPTIDEK");
        assert_eq!(f32_at(b, "purity", &p), [Some(0.9), None]);
        assert_eq!(f32_at(b, "itraq_114", &p), [Some(1.1), Some(1.1)]);
        assert_eq!(f32_at(b, "raw_itraq_114", &p), [Some(1.0), Some(1.0)]);
        assert_eq!(f32_at(b, "reporter_intensity_sum", &p), [Some(10.0); 2]);
        assert_eq!(f32_at(b, "reporter_sn_sum", &p), [Some(104.0), None]);
        assert_eq!(f32_at(b, "quant_weight", &p), [Some(0.81), Some(0.0)]);
        assert_eq!(f32_at(b, "sn_itraq_117", &p), [Some(41.0), None]);
        let q = column::<Float64Array>(b, "psm_q_value", &p).unwrap();
        assert_eq!((opt(q, 0), opt(q, 1)), (Some(0.001), None));
        let passed = column::<BooleanArray>(b, "quant_gate_passed", &p).unwrap();
        assert_eq!((passed.value(0), passed.value(1)), (true, false));
        let scan = column::<Int32Array>(b, "scannr", &p).unwrap();
        assert_eq!(scan.value(0), 5);
        // Without a correction matrix the channels hold the raw values and
        // there is no raw block.
        write_tmt_parquet(&p, Plex::Itraq4, &[tmt_row()], false, false).unwrap();
        let t = read_tmt_parquet(&p).unwrap();
        assert_eq!(names(&t).len(), 13 + 4 + 4 + 4);
        assert_eq!(f32_at(&t.batches[0], "itraq_114", &p), [Some(1.0)]);
        // A row whose channel count does not match the plex is refused.
        let mut short = tmt_row();
        short.raw.pop();
        assert!(write_tmt_parquet(&p, Plex::Itraq4, &[short], false, false).is_err());
        assert!(!dir.path().join("x.tmt.parquet.partial").exists());
    }

    #[test]
    fn tmt_parquet_appends_the_unmix_block_after_the_existing_columns() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("u.tmt.parquet");
        let mut modelled = tmt_row();
        modelled.interference = Some(0.25);
        modelled.unmixed = Some(vec![0.5, 1.5, 2.5, 3.5]);
        write_tmt_parquet(&p, Plex::Itraq4, &[modelled, tmt_row()], true, true).unwrap();
        let t = read_tmt_parquet(&p).unwrap();
        assert_eq!(t.meta.get(TMT_UNMIX_KEY).map(String::as_str), Some("true"));
        let header = names(&t);
        // The columns of a run without --tmt-unmix keep their positions.
        assert_eq!(header.len(), 13 + 8 + 4 + 4 + 1 + 4);
        assert_eq!(header[25], "sn_itraq_114");
        assert_eq!(header[29], "interference");
        assert_eq!(header[30], "unmix_itraq_114");
        let b = &t.batches[0];
        assert_eq!(f32_at(b, "interference", &p), [Some(0.25), None]);
        assert_eq!(f32_at(b, "unmix_itraq_117", &p), [Some(3.5), None]);
        assert_eq!(f32_at(b, "itraq_117", &p), [Some(4.1), Some(4.1)]);
        let mut bad = tmt_row();
        bad.unmixed = Some(vec![1.0]);
        assert!(write_tmt_parquet(&p, Plex::Itraq4, &[bad], true, true).is_err());
        write_tmt_parquet(&p, Plex::Itraq4, &[tmt_row()], true, false).unwrap();
        let t = read_tmt_parquet(&p).unwrap();
        assert_eq!(t.meta.get(TMT_UNMIX_KEY).map(String::as_str), Some("false"));
        assert_eq!(names(&t).len(), 13 + 8 + 4 + 4);
    }

    fn lfq_rows() -> Vec<LfqRow> {
        vec![
            LfqRow {
                id: id("a", 0, "PEPTIDEK", 2),
                feature: feature(1e6),
                score: 0.9,
                feature_q_value: Some(0.001),
                decoy_score: None,
            },
            LfqRow {
                id: id("b", 1, "PEPTIDEK", 2),
                feature: feature(2e6),
                score: 0.8,
                feature_q_value: Some(0.002),
                decoy_score: Some(0.1),
            },
            LfqRow {
                id: id("b", 1, "OTHERK", 3),
                feature: feature(5e5),
                score: 0.5,
                feature_q_value: None,
                decoy_score: Some(0.6),
            },
        ]
    }

    #[test]
    fn lfq_parquet_is_long_with_the_accepted_flag_and_runs() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.lfq.parquet");
        let runs = ["a".to_string(), "b".to_string(), "empty".to_string()];
        let meta = LfqTableMeta {
            runs: &runs,
            feature_fdr: 0.01,
            min_cosine: 0.7,
        };
        let rows = lfq_rows();
        write_lfq_parquet(&p, &rows, &[true, true, false], &meta).unwrap();
        let t = read_lfq_parquet(&p).unwrap();
        assert_eq!(
            t.meta.get(LFQ_RUNS_KEY).map(String::as_str),
            Some("a\nb\nempty")
        );
        assert_eq!(
            t.meta.get(LFQ_FEATURE_FDR_KEY).map(String::as_str),
            Some("0.01")
        );
        assert_eq!(
            t.meta.get(LFQ_MIN_COSINE_KEY).map(String::as_str),
            Some("0.7")
        );
        assert_eq!(
            names(&t),
            [
                "filename",
                "peptide",
                "charge",
                "proteins",
                "scannr",
                "spec_id",
                "psm_q_value",
                "calc_mz",
                "rt_anchor",
                "rt_apex",
                "rt_start",
                "rt_stop",
                "n_scans",
                "area",
                "apex_intensity",
                "isotope_areas",
                "cosine",
                "n_isotopes",
                "score",
                "decoy_score",
                "feature_q_value",
                "accepted"
            ]
        );
        let b = &t.batches[0];
        assert_eq!(b.num_rows(), 3);
        let area = column::<Float64Array>(b, "area", &p).unwrap();
        assert_eq!(area.value(1), 2e6);
        let fq = column::<Float64Array>(b, "feature_q_value", &p).unwrap();
        assert_eq!((opt(fq, 0), opt(fq, 2)), (Some(0.001), None));
        assert_eq!(f32_at(b, "decoy_score", &p), [None, Some(0.1), Some(0.6)]);
        let iso = column::<ListArray>(b, "isotope_areas", &p).unwrap();
        let first = iso.value(0);
        let first = first.as_any().downcast_ref::<Float64Array>().unwrap();
        assert_eq!(first.values().to_vec(), rows[0].feature.isotope_areas);
        let acc = column::<BooleanArray>(b, "accepted", &p).unwrap();
        assert_eq!(
            (0..3).map(|i| acc.value(i)).collect::<Vec<_>>(),
            [true, true, false]
        );
        let n = column::<UInt32Array>(b, "n_scans", &p).unwrap();
        assert_eq!(n.value(0), 13);
        assert!(write_lfq_parquet(&p, &rows, &[true], &meta).is_err());
    }

    #[test]
    fn quant_readers_reject_other_majors_and_each_others_files() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.parquet");
        let schema = Arc::new(tmt_schema(Plex::Itraq4, false, false));
        let b = tmt_batch(&schema, 4, false, false, &[tmt_row()]);
        write_table(&p, &TMT_TABLE, "2.1", schema, [b], &[]).unwrap();
        let e = read_tmt_parquet(&p).err().unwrap();
        assert!(e.contains("x.parquet") && e.contains("schema 2.1"), "{e}");
        let e = read_lfq_parquet(&p).err().unwrap();
        assert!(e.contains(LFQ_TABLE.key), "{e}");
    }
}
