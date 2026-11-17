use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanArray, BooleanBuilder, Float64Array, Float64Builder, StringArray,
    StringBuilder, UInt32Array, UInt32Builder, UInt64Array,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::versioned_parquet::{
    column, opt, read_table, write_table, TableVersion, TableWriter, BATCH_ROWS,
};

const MBR_TABLE: TableVersion = TableVersion {
    key: "andes.mbr.schema",
    version: "1.3",
    major: 1,
    what: "MBR table",
};
const VALIDATION_TABLE: TableVersion = TableVersion {
    key: "andes.mbr_validation.schema",
    version: "1.2",
    major: 1,
    what: "MBR validation table",
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceType {
    Identified,
    IdentifiedBelowThreshold,
    TransferMs2,
    TransferMs1,
    SharedFeature,
    BelowThreshold,
}

impl EvidenceType {
    pub fn is_transfer(self) -> bool {
        matches!(
            self,
            EvidenceType::TransferMs2
                | EvidenceType::TransferMs1
                | EvidenceType::SharedFeature
                | EvidenceType::BelowThreshold
        )
    }

    pub fn as_str(self) -> &'static str {
        match self {
            EvidenceType::Identified => "identified",
            EvidenceType::IdentifiedBelowThreshold => "identified_below_threshold",
            EvidenceType::TransferMs2 => "transfer_ms2",
            EvidenceType::TransferMs1 => "transfer_ms1",
            EvidenceType::SharedFeature => "shared_feature",
            EvidenceType::BelowThreshold => "below_threshold",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MbrRow<'a> {
    pub precursor_key: &'a str,
    pub run: &'a str,
    pub evidence_type: EvidenceType,
    pub quantity: Option<f64>,
    pub quantity_se: Option<f64>,
    pub quantity_normalized: Option<f64>,
    pub transfer_q: Option<f64>,
    pub quant_q: Option<f64>,
    pub apex_rt_s: Option<f64>,
    pub library_q: f64,
    pub claim: Option<MbrClaim<'a>>,
    pub measurement_q: Option<f64>,
    pub transfer_pep: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MbrClaim<'a> {
    pub claimed_by: Option<&'a str>,
    pub dt_s: Option<f64>,
    pub identified: bool,
    pub co_apex: bool,
    pub n_co_apex: u32,
    pub n_identified_co_apex: u32,
}

pub const MBR_COLUMNS: [&str; 18] = [
    "precursor_key",
    "run",
    "evidence_type",
    "quantity",
    "quantity_se",
    "quantity_normalized",
    "transfer_q",
    "quant_q",
    "apex_rt_s",
    "library_q",
    "claimed_by",
    "claim_dt_s",
    "claim_identified",
    "claim_co_apex",
    "claim_n_co_apex",
    "claim_n_identified_co_apex",
    "measurement_q",
    "transfer_pep",
];

pub const VALIDATION_Q: f64 = 0.01;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Accepted {
    pub pseudo: usize,
    pub pairs: usize,
    pub wrong_peak: usize,
    pub wrong_measurement: usize,
    pub pooled: usize,
    pub mass_only: usize,
    pub rt_only: usize,
    pub wrong_peak_only: usize,
    pub decoys: [usize; 3],
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MbrValidation {
    pub q_threshold: f64,
    pub pseudo_transfers: usize,
    pub pseudo_scored: usize,
    pub wrong_peak_pairs: usize,
    pub wrong_peak: usize,
    pub wrong_measurement: usize,
    pub transfers: usize,
    pub transfers_scored: usize,
    pub at_q01: Accepted,
    pub at_threshold: Accepted,
    pub quant_feature_threshold: Option<f64>,
    pub identified: usize,
    pub identified_scored: usize,
    pub identified_q01: usize,
    pub identified_below_threshold: usize,
    pub ms2_rows: [usize; 3],
    pub ms2_supported: [usize; 3],
    pub ms2_pseudo_rows: usize,
    pub ms2_pseudo_supported: usize,
    pub transfer_ms2: usize,
    pub target_entries: usize,
    pub identity_entries: usize,
    pub decoy_transfers: [usize; 3],
    pub shared_feature: usize,
}

fn add_counts<const N: usize>(a: &mut [usize; N], b: &[usize; N]) {
    for (x, y) in a.iter_mut().zip(b) {
        *x += y;
    }
}

impl Accepted {
    pub fn add(&mut self, o: &Accepted) {
        let Accepted {
            pseudo,
            pairs,
            wrong_peak,
            wrong_measurement,
            pooled,
            mass_only,
            rt_only,
            wrong_peak_only,
            decoys,
        } = o;
        self.pseudo += pseudo;
        self.pairs += pairs;
        self.wrong_peak += wrong_peak;
        self.wrong_measurement += wrong_measurement;
        self.pooled += pooled;
        self.mass_only += mass_only;
        self.rt_only += rt_only;
        self.wrong_peak_only += wrong_peak_only;
        add_counts(&mut self.decoys, decoys);
    }
}

impl MbrValidation {
    pub fn add_counts(&mut self, o: &MbrValidation) {
        let MbrValidation {
            q_threshold: _,
            pseudo_transfers,
            pseudo_scored,
            wrong_peak_pairs,
            wrong_peak,
            wrong_measurement,
            transfers,
            transfers_scored,
            at_q01,
            at_threshold,
            quant_feature_threshold: _,
            identified,
            identified_scored,
            identified_q01,
            identified_below_threshold,
            ms2_rows,
            ms2_supported,
            ms2_pseudo_rows,
            ms2_pseudo_supported,
            transfer_ms2,
            target_entries: _,
            identity_entries: _,
            decoy_transfers,
            shared_feature,
        } = o;
        self.pseudo_transfers += pseudo_transfers;
        self.pseudo_scored += pseudo_scored;
        self.wrong_peak_pairs += wrong_peak_pairs;
        self.wrong_peak += wrong_peak;
        self.wrong_measurement += wrong_measurement;
        self.transfers += transfers;
        self.transfers_scored += transfers_scored;
        self.at_q01.add(at_q01);
        self.at_threshold.add(at_threshold);
        self.identified += identified;
        self.identified_scored += identified_scored;
        self.identified_q01 += identified_q01;
        self.identified_below_threshold += identified_below_threshold;
        add_counts(&mut self.ms2_rows, ms2_rows);
        add_counts(&mut self.ms2_supported, ms2_supported);
        self.ms2_pseudo_rows += ms2_pseudo_rows;
        self.ms2_pseudo_supported += ms2_pseudo_supported;
        self.transfer_ms2 += transfer_ms2;
        add_counts(&mut self.decoy_transfers, decoy_transfers);
        self.shared_feature += shared_feature;
    }
}

const DECOY_ACCEPTED_COLUMNS: [&str; 4] = [
    "accepted_wrong_peak_only",
    "accepted_mass_decoys",
    "accepted_rt_decoys",
    "accepted_identity_decoys",
];

const WRONG_MEASUREMENT_COLUMNS: [&str; 2] = [
    "accepted_wrong_measurement",
    "accepted_wrong_measurement_rate",
];

const ACCEPTED_COLUMNS: [&str; 8] = [
    "pseudo_accepted",
    "pseudo_acceptance",
    "accepted_pairs",
    "accepted_wrong_peak",
    "accepted_wrong_peak_rate",
    "accepted_pooled",
    "accepted_mass_only",
    "accepted_rt_only",
];

pub fn validation_columns() -> Vec<String> {
    let mut c: Vec<String> = [
        "q_threshold",
        "pseudo_transfers",
        "pseudo_scored",
        "wrong_peak_pairs",
        "wrong_peak",
        "wrong_peak_rate",
        "transfers",
        "transfers_scored",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for suffix in ["_q01", "_at_threshold"] {
        c.extend(ACCEPTED_COLUMNS.iter().map(|n| format!("{n}{suffix}")));
    }
    c.extend(
        [
            "quant_feature_threshold",
            "identified",
            "identified_scored",
            "identified_accepted_q01",
            "identified_below_threshold",
            "ms2_target_rows",
            "ms2_target_supported",
            "ms2_mass_decoy_rows",
            "ms2_mass_decoy_supported",
            "ms2_rt_decoy_rows",
            "ms2_rt_decoy_supported",
            "ms2_pseudo_rows",
            "ms2_pseudo_supported",
            "transfer_ms2",
        ]
        .iter()
        .map(|s| s.to_string()),
    );
    for suffix in ["_q01", "_at_threshold"] {
        c.extend(
            DECOY_ACCEPTED_COLUMNS
                .iter()
                .map(|n| format!("{n}{suffix}")),
        );
    }
    c.extend(
        [
            "target_entries",
            "identity_entries",
            "mass_decoy_transfers",
            "rt_decoy_transfers",
            "identity_decoy_transfers",
            "shared_feature",
            "wrong_measurement",
            "wrong_measurement_rate",
        ]
        .iter()
        .map(|s| s.to_string()),
    );
    for suffix in ["_q01", "_at_threshold"] {
        c.extend(
            WRONG_MEASUREMENT_COLUMNS
                .iter()
                .map(|n| format!("{n}{suffix}")),
        );
    }
    c
}

fn ratio(num: usize, den: usize) -> Option<f64> {
    (den > 0).then(|| num as f64 / den as f64)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ValidationCell {
    Count(u64),
    Value(Option<f64>),
}

impl MbrValidation {
    pub fn pseudo_acceptance(&self, a: &Accepted) -> Option<f64> {
        ratio(a.pseudo, self.pseudo_transfers)
    }

    pub fn wrong_peak_rate(&self) -> Option<f64> {
        ratio(self.wrong_peak, self.wrong_peak_pairs)
    }

    pub fn wrong_measurement_rate(&self) -> Option<f64> {
        ratio(self.wrong_measurement, self.wrong_peak_pairs)
    }

    pub fn cells(&self) -> Vec<ValidationCell> {
        use ValidationCell::{Count, Value};
        let n = |v: usize| Count(v as u64);
        let mut c = vec![
            Value(Some(self.q_threshold)),
            n(self.pseudo_transfers),
            n(self.pseudo_scored),
            n(self.wrong_peak_pairs),
            n(self.wrong_peak),
            Value(self.wrong_peak_rate()),
            n(self.transfers),
            n(self.transfers_scored),
        ];
        for a in [&self.at_q01, &self.at_threshold] {
            c.extend([
                n(a.pseudo),
                Value(self.pseudo_acceptance(a)),
                n(a.pairs),
                n(a.wrong_peak),
                Value(a.wrong_peak_rate()),
                n(a.pooled),
                n(a.mass_only),
                n(a.rt_only),
            ]);
        }
        c.push(Value(self.quant_feature_threshold));
        c.extend([
            n(self.identified),
            n(self.identified_scored),
            n(self.identified_q01),
            n(self.identified_below_threshold),
        ]);
        for f in 0..3 {
            c.extend([n(self.ms2_rows[f]), n(self.ms2_supported[f])]);
        }
        c.extend([
            n(self.ms2_pseudo_rows),
            n(self.ms2_pseudo_supported),
            n(self.transfer_ms2),
        ]);
        for a in [&self.at_q01, &self.at_threshold] {
            c.extend([
                n(a.wrong_peak_only),
                n(a.decoys[0]),
                n(a.decoys[1]),
                n(a.decoys[2]),
            ]);
        }
        c.extend([
            n(self.target_entries),
            n(self.identity_entries),
            n(self.decoy_transfers[0]),
            n(self.decoy_transfers[1]),
            n(self.decoy_transfers[2]),
            n(self.shared_feature),
            n(self.wrong_measurement),
            Value(self.wrong_measurement_rate()),
        ]);
        for a in [&self.at_q01, &self.at_threshold] {
            c.extend([n(a.wrong_measurement), Value(a.wrong_measurement_rate())]);
        }
        c.iter()
            .map(|v| match *v {
                Value(x) => Value(x.filter(|x| x.is_finite())),
                count => count,
            })
            .collect()
    }
}

impl Accepted {
    pub fn wrong_peak_rate(&self) -> Option<f64> {
        ratio(self.wrong_peak, self.pairs)
    }

    pub fn wrong_measurement_rate(&self) -> Option<f64> {
        ratio(self.wrong_measurement, self.pairs)
    }
}

fn mbr_schema() -> Schema {
    let f = |name: &str| Field::new(name, DataType::Float64, true);
    Schema::new(vec![
        Field::new("precursor_key", DataType::Utf8, false),
        Field::new("run", DataType::Utf8, false),
        Field::new("evidence_type", DataType::Utf8, false),
        f("quantity"),
        f("quantity_se"),
        f("quantity_normalized"),
        f("transfer_q"),
        f("quant_q"),
        f("apex_rt_s"),
        Field::new("library_q", DataType::Float64, false),
        Field::new("claimed_by", DataType::Utf8, true),
        f("claim_dt_s"),
        Field::new("claim_identified", DataType::Boolean, true),
        Field::new("claim_co_apex", DataType::Boolean, true),
        Field::new("claim_n_co_apex", DataType::UInt32, true),
        Field::new("claim_n_identified_co_apex", DataType::UInt32, true),
        f("measurement_q"),
        f("transfer_pep"),
    ])
}

struct MbrColumnsBuilder {
    key: StringBuilder,
    run: StringBuilder,
    ev: StringBuilder,
    floats: Vec<Float64Builder>,
    library_q: Float64Builder,
    claimed_by: StringBuilder,
    claim_dt: Float64Builder,
    claim_identified: BooleanBuilder,
    claim_co_apex: BooleanBuilder,
    claim_n_co: UInt32Builder,
    claim_n_id_co: UInt32Builder,
    measurement_q: Float64Builder,
    transfer_pep: Float64Builder,
    rows: usize,
}

impl MbrColumnsBuilder {
    fn new() -> Self {
        Self {
            key: StringBuilder::new(),
            run: StringBuilder::new(),
            ev: StringBuilder::new(),
            floats: (0..6).map(|_| Float64Builder::new()).collect(),
            library_q: Float64Builder::new(),
            claimed_by: StringBuilder::new(),
            claim_dt: Float64Builder::new(),
            claim_identified: BooleanBuilder::new(),
            claim_co_apex: BooleanBuilder::new(),
            claim_n_co: UInt32Builder::new(),
            claim_n_id_co: UInt32Builder::new(),
            measurement_q: Float64Builder::new(),
            transfer_pep: Float64Builder::new(),
            rows: 0,
        }
    }

    fn push(&mut self, r: &MbrRow<'_>) {
        self.key.append_value(r.precursor_key);
        self.run.append_value(r.run);
        self.ev.append_value(r.evidence_type.as_str());
        let values = [
            r.quantity,
            r.quantity_se,
            r.quantity_normalized,
            r.transfer_q,
            r.quant_q,
            r.apex_rt_s,
        ];
        for (b, v) in self.floats.iter_mut().zip(values) {
            b.append_option(v.filter(|v| v.is_finite()));
        }
        self.library_q.append_value(r.library_q);
        let c = r.claim.as_ref();
        self.claimed_by.append_option(c.and_then(|c| c.claimed_by));
        self.claim_dt
            .append_option(c.and_then(|c| c.dt_s).filter(|v| v.is_finite()));
        self.claim_identified.append_option(c.map(|c| c.identified));
        self.claim_co_apex.append_option(c.map(|c| c.co_apex));
        self.claim_n_co.append_option(c.map(|c| c.n_co_apex));
        self.claim_n_id_co
            .append_option(c.map(|c| c.n_identified_co_apex));
        self.measurement_q
            .append_option(r.measurement_q.filter(|v| v.is_finite()));
        self.transfer_pep
            .append_option(r.transfer_pep.filter(|v| v.is_finite()));
        self.rows += 1;
    }

    fn batch(&mut self, schema: &Arc<Schema>) -> Result<RecordBatch, String> {
        self.rows = 0;
        let mut cols: Vec<ArrayRef> = vec![
            Arc::new(self.key.finish()),
            Arc::new(self.run.finish()),
            Arc::new(self.ev.finish()),
        ];
        cols.extend(
            self.floats
                .iter_mut()
                .map(|b| Arc::new(b.finish()) as ArrayRef),
        );
        cols.push(Arc::new(self.library_q.finish()));
        cols.push(Arc::new(self.claimed_by.finish()));
        cols.push(Arc::new(self.claim_dt.finish()));
        cols.push(Arc::new(self.claim_identified.finish()));
        cols.push(Arc::new(self.claim_co_apex.finish()));
        cols.push(Arc::new(self.claim_n_co.finish()));
        cols.push(Arc::new(self.claim_n_id_co.finish()));
        cols.push(Arc::new(self.measurement_q.finish()));
        cols.push(Arc::new(self.transfer_pep.finish()));
        RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
    }
}

fn mbr_batch(schema: &Arc<Schema>, rows: &[MbrRow<'_>]) -> Result<RecordBatch, String> {
    let mut b = MbrColumnsBuilder::new();
    for r in rows {
        b.push(r);
    }
    b.batch(schema)
}

pub struct MbrParquetWriter {
    schema: Arc<Schema>,
    columns: MbrColumnsBuilder,
    table: TableWriter,
    path: std::path::PathBuf,
}

impl MbrParquetWriter {
    pub fn create(path: &Path) -> Result<Self, String> {
        let schema = Arc::new(mbr_schema());
        let table = TableWriter::create(path, &MBR_TABLE, MBR_TABLE.version, schema.clone(), &[])?;
        Ok(Self {
            schema,
            columns: MbrColumnsBuilder::new(),
            table,
            path: path.to_path_buf(),
        })
    }

    pub fn push(&mut self, row: &MbrRow<'_>) -> Result<(), String> {
        self.columns.push(row);
        if self.columns.rows == BATCH_ROWS {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), String> {
        let b = self
            .columns
            .batch(&self.schema)
            .map_err(|e| format!("{}: {e}", self.path.display()))?;
        self.table.write(&b)
    }

    pub fn finish(mut self) -> Result<(), String> {
        if self.columns.rows > 0 {
            self.flush()?;
        }
        self.table.finish()
    }
}

pub fn write_mbr_parquet(path: &Path, rows: &[MbrRow<'_>]) -> Result<(), String> {
    let schema = Arc::new(mbr_schema());
    let batches = rows.chunks(BATCH_ROWS).map(|c| mbr_batch(&schema, c));
    write_table(
        path,
        &MBR_TABLE,
        MBR_TABLE.version,
        schema.clone(),
        batches,
        &[],
    )
}

#[derive(Debug, Clone, PartialEq)]
pub struct MbrRecord {
    pub precursor_key: String,
    pub run: String,
    pub evidence_type: String,
    pub quantity: Option<f64>,
    pub quantity_se: Option<f64>,
    pub quantity_normalized: Option<f64>,
    pub transfer_q: Option<f64>,
    pub quant_q: Option<f64>,
    pub apex_rt_s: Option<f64>,
    pub library_q: f64,
    pub claimed_by: Option<String>,
    pub claim_co_apex: Option<bool>,
    pub claim_n_co_apex: Option<u32>,
    pub measurement_q: Option<f64>,
}

pub fn read_mbr_parquet(path: &Path) -> Result<Vec<MbrRecord>, String> {
    let t = read_table(path, &MBR_TABLE)?;
    let mut out = Vec::new();
    for b in &t.batches {
        let text = |name: &str| column::<StringArray>(b, name, path);
        let float = |name: &str| column::<Float64Array>(b, name, path);
        let (key, run, ev) = (text("precursor_key")?, text("run")?, text("evidence_type")?);
        let [q, se, norm, tq, qq, rt] = [
            "quantity",
            "quantity_se",
            "quantity_normalized",
            "transfer_q",
            "quant_q",
            "apex_rt_s",
        ]
        .map(float);
        let (q, se, norm, tq, qq, rt) = (q?, se?, norm?, tq?, qq?, rt?);
        let lib = float("library_q")?;
        let claimed = b
            .column_by_name("claimed_by")
            .map(|_| text("claimed_by"))
            .transpose()?;
        let co_apex = b
            .column_by_name("claim_co_apex")
            .map(|_| column::<BooleanArray>(b, "claim_co_apex", path))
            .transpose()?;
        let n_co = b
            .column_by_name("claim_n_co_apex")
            .map(|_| column::<UInt32Array>(b, "claim_n_co_apex", path))
            .transpose()?;
        let mq = b
            .column_by_name("measurement_q")
            .map(|_| float("measurement_q"))
            .transpose()?;
        for i in 0..b.num_rows() {
            out.push(MbrRecord {
                precursor_key: key.value(i).to_string(),
                run: run.value(i).to_string(),
                evidence_type: ev.value(i).to_string(),
                quantity: opt(q, i),
                quantity_se: opt(se, i),
                quantity_normalized: opt(norm, i),
                transfer_q: opt(tq, i),
                quant_q: opt(qq, i),
                apex_rt_s: opt(rt, i),
                library_q: lib.value(i),
                claimed_by: claimed
                    .filter(|c| !c.is_null(i))
                    .map(|c| c.value(i).to_string()),
                claim_co_apex: co_apex.filter(|c| !c.is_null(i)).map(|c| c.value(i)),
                claim_n_co_apex: n_co.and_then(|c| opt(c, i)),
                measurement_q: mq.and_then(|c| opt(c, i)),
            });
        }
    }
    Ok(out)
}

pub fn write_mbr_validation(path: &Path, v: &MbrValidation) -> Result<(), String> {
    let names = validation_columns();
    let cells = v.cells();
    debug_assert_eq!(names.len(), cells.len());
    let mut fields = Vec::with_capacity(names.len());
    let mut cols: Vec<ArrayRef> = Vec::with_capacity(names.len());
    for (name, cell) in names.iter().zip(&cells) {
        match *cell {
            ValidationCell::Count(n) => {
                fields.push(Field::new(name, DataType::UInt64, false));
                cols.push(Arc::new(UInt64Array::from(vec![n])));
            }
            ValidationCell::Value(x) => {
                fields.push(Field::new(name, DataType::Float64, true));
                cols.push(Arc::new(Float64Array::from(vec![x])));
            }
        }
    }
    let schema = Arc::new(Schema::new(fields));
    let batch = RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string());
    write_table(
        path,
        &VALIDATION_TABLE,
        VALIDATION_TABLE.version,
        schema,
        [batch],
        &[],
    )
}

pub fn read_mbr_validation(path: &Path) -> Result<Vec<(String, ValidationCell)>, String> {
    let t = read_table(path, &VALIDATION_TABLE)?;
    let b = match t.batches.as_slice() {
        [b] if b.num_rows() == 1 => b,
        _ => return Err(format!("{}: expected one validation row", path.display())),
    };
    let schema = b.schema();
    let mut out = Vec::with_capacity(b.num_columns());
    for f in schema.fields() {
        let name = f.name();
        let cell = match f.data_type() {
            DataType::UInt64 => {
                ValidationCell::Count(column::<UInt64Array>(b, name, path)?.value(0))
            }
            DataType::Float64 => {
                ValidationCell::Value(opt(column::<Float64Array>(b, name, path)?, 0))
            }
            other => {
                return Err(format!(
                    "{}: column {name} has unexpected type {other}",
                    path.display()
                ))
            }
        };
        out.push((name.clone(), cell));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_by_row_writer_gives_the_same_file() {
        let dir = tempfile::tempdir().unwrap();
        let base = rows();
        let many: Vec<MbrRow<'static>> = (0..BATCH_ROWS + 1234)
            .map(|i| MbrRow {
                quantity: Some(i as f64),
                ..base[i % base.len()].clone()
            })
            .collect();
        for set in [&many[..0], &many[..3], &many[..BATCH_ROWS], &many[..]] {
            let (a, b) = (dir.path().join("a.parquet"), dir.path().join("b.parquet"));
            write_mbr_parquet(&a, set).unwrap();
            let mut w = MbrParquetWriter::create(&b).unwrap();
            for r in set {
                w.push(r).unwrap();
            }
            w.finish().unwrap();
            assert_eq!(
                std::fs::read(&a).unwrap(),
                std::fs::read(&b).unwrap(),
                "{}",
                set.len()
            );
            assert!(!dir.path().join("b.parquet.partial").exists());
        }
        let c = dir.path().join("c.parquet");
        drop(MbrParquetWriter::create(&c).unwrap());
        assert!(!c.exists() && !dir.path().join("c.parquet.partial").exists());
    }

    fn rows() -> Vec<MbrRow<'static>> {
        vec![
            MbrRow {
                precursor_key: "PEPTIDEK/2",
                run: "a",
                evidence_type: EvidenceType::Identified,
                quantity: Some(1.5e6),
                quantity_se: Some(2e4),
                quantity_normalized: Some(1.25e6),
                transfer_q: None,
                quant_q: Some(0.003),
                apex_rt_s: Some(612.25),
                library_q: 0.001,
                claim: None,
                measurement_q: None,
                transfer_pep: None,
            },
            MbrRow {
                precursor_key: "PEPTIDEK/2",
                run: "b",
                evidence_type: EvidenceType::BelowThreshold,
                quantity: Some(f64::NAN),
                quantity_se: None,
                quantity_normalized: None,
                transfer_q: Some(0.2),
                quant_q: Some(0.2),
                apex_rt_s: None,
                library_q: 0.001,
                claim: Some(MbrClaim {
                    claimed_by: Some("OTHERPEPK/2"),
                    dt_s: Some(-1.5),
                    identified: true,
                    co_apex: true,
                    n_co_apex: 2,
                    n_identified_co_apex: 1,
                }),
                measurement_q: None,
                transfer_pep: None,
            },
            MbrRow {
                precursor_key: "PEPTIDEK/2",
                run: "c",
                evidence_type: EvidenceType::TransferMs1,
                quantity: Some(9e5),
                quantity_se: Some(1e4),
                quantity_normalized: None,
                transfer_q: Some(0.004),
                quant_q: Some(0.006),
                apex_rt_s: Some(611.0),
                library_q: 0.001,
                claim: Some(MbrClaim {
                    claimed_by: None,
                    dt_s: None,
                    identified: false,
                    co_apex: false,
                    n_co_apex: 0,
                    n_identified_co_apex: 0,
                }),
                measurement_q: Some(0.006),
                transfer_pep: None,
            },
        ]
    }

    #[test]
    fn mbr_parquet_round_trips_with_nulls_for_missing_values() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr.parquet");
        write_mbr_parquet(&p, &rows()).unwrap();
        assert!(!dir.path().join("x.mbr.parquet.partial").exists());
        let back = read_mbr_parquet(&p).unwrap();
        let want: Vec<MbrRecord> = rows()
            .iter()
            .map(|r| MbrRecord {
                precursor_key: r.precursor_key.into(),
                run: r.run.into(),
                evidence_type: r.evidence_type.as_str().into(),
                quantity: r.quantity.filter(|v| v.is_finite()),
                quantity_se: r.quantity_se,
                quantity_normalized: r.quantity_normalized,
                transfer_q: r.transfer_q,
                quant_q: r.quant_q,
                apex_rt_s: r.apex_rt_s,
                library_q: r.library_q,
                claimed_by: r
                    .claim
                    .as_ref()
                    .and_then(|c| c.claimed_by)
                    .map(String::from),
                claim_co_apex: r.claim.as_ref().map(|c| c.co_apex),
                claim_n_co_apex: r.claim.as_ref().map(|c| c.n_co_apex),
                measurement_q: r.measurement_q,
            })
            .collect();
        assert_eq!(back, want);
        assert_eq!(back[2].measurement_q, Some(0.006));
        assert_eq!(back[1].quantity, None);
        let t = read_table(&p, &MBR_TABLE).unwrap();
        let names: Vec<String> = t.batches[0]
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect();
        assert_eq!(names, MBR_COLUMNS);
        assert_eq!(t.meta.get(MBR_TABLE.key).map(String::as_str), Some("1.3"));
        write_mbr_parquet(&p, &[]).unwrap();
        assert!(read_mbr_parquet(&p).unwrap().is_empty());
    }

    #[test]
    fn shared_feature_is_a_named_transfer_type() {
        assert_eq!(EvidenceType::SharedFeature.as_str(), "shared_feature");
        assert!(EvidenceType::SharedFeature.is_transfer());
        assert!(!EvidenceType::Identified.is_transfer());
    }

    #[test]
    fn mbr_readers_reject_other_majors_and_each_others_files() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr.parquet");
        let schema = Arc::new(mbr_schema());
        let b = mbr_batch(&schema, &rows());
        write_table(&p, &MBR_TABLE, "2.0", schema, [b], &[]).unwrap();
        let e = read_mbr_parquet(&p).unwrap_err();
        assert!(
            e.contains("x.mbr.parquet") && e.contains("schema 2.0"),
            "{e}"
        );
        write_mbr_parquet(&p, &rows()).unwrap();
        let e = read_mbr_validation(&p).unwrap_err();
        assert!(e.contains(VALIDATION_TABLE.key), "{e}");
    }

    #[test]
    fn validation_is_one_typed_row_with_rates_at_both_thresholds() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr_validation.parquet");
        let v = MbrValidation {
            q_threshold: 0.05,
            pseudo_transfers: 10,
            pseudo_scored: 9,
            wrong_peak_pairs: 8,
            wrong_peak: 1,
            wrong_measurement: 2,
            transfers: 100,
            transfers_scored: 90,
            at_q01: Accepted {
                pseudo: 8,
                pairs: 0,
                wrong_peak: 0,
                wrong_measurement: 0,
                pooled: 40,
                mass_only: 41,
                rt_only: 39,
                wrong_peak_only: 44,
                decoys: [1, 2, 3],
            },
            at_threshold: Accepted {
                pseudo: 9,
                pairs: 4,
                wrong_peak: 1,
                wrong_measurement: 3,
                pooled: 60,
                mass_only: 61,
                rt_only: 59,
                wrong_peak_only: 66,
                decoys: [4, 5, 6],
            },
            quant_feature_threshold: None,
            identified: 50,
            identified_scored: 48,
            identified_q01: 45,
            identified_below_threshold: 0,
            ms2_rows: [100, 90, 80],
            ms2_supported: [30, 1, 2],
            ms2_pseudo_rows: 10,
            ms2_pseudo_supported: 4,
            transfer_ms2: 25,
            target_entries: 70,
            identity_entries: 2,
            decoy_transfers: [95, 85, 7],
            shared_feature: 6,
        };
        write_mbr_validation(&p, &v).unwrap();
        let cells = read_mbr_validation(&p).unwrap();
        let names: Vec<&str> = cells.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, validation_columns());
        assert_eq!(cells.iter().map(|(_, c)| *c).collect::<Vec<_>>(), v.cells());
        let get = |name: &str| cells.iter().find(|(n, _)| n == name).unwrap().1;
        use ValidationCell::{Count, Value};
        assert_eq!(get("q_threshold"), Value(Some(0.05)));
        assert_eq!(get("wrong_peak_rate"), Value(Some(0.125)));
        assert_eq!(get("pseudo_acceptance_q01"), Value(Some(0.8)));
        assert_eq!(get("pseudo_acceptance_at_threshold"), Value(Some(0.9)));
        assert_eq!(get("accepted_wrong_peak_rate_q01"), Value(None));
        assert_eq!(
            get("accepted_wrong_peak_rate_at_threshold"),
            Value(Some(0.25))
        );
        assert_eq!(get("accepted_rt_only_q01"), Count(39));
        assert_eq!(get("accepted_pooled_at_threshold"), Count(60));
        assert_eq!(get("quant_feature_threshold"), Value(None));
        assert_eq!(get("identified_scored"), Count(48));
        assert_eq!(get("identified_accepted_q01"), Count(45));
        assert_eq!(get("identified_below_threshold"), Count(0));
        assert_eq!(get("ms2_target_rows"), Count(100));
        assert_eq!(get("ms2_mass_decoy_supported"), Count(1));
        assert_eq!(get("ms2_rt_decoy_rows"), Count(80));
        assert_eq!(get("ms2_pseudo_supported"), Count(4));
        assert_eq!(get("transfer_ms2"), Count(25));
        assert_eq!(get("accepted_wrong_peak_only_q01"), Count(44));
        assert_eq!(get("accepted_identity_decoys_q01"), Count(3));
        assert_eq!(get("accepted_mass_decoys_at_threshold"), Count(4));
        assert_eq!(get("identity_entries"), Count(2));
        assert_eq!(get("identity_decoy_transfers"), Count(7));
        assert_eq!(get("shared_feature"), Count(6));
        assert_eq!(get("wrong_measurement"), Count(2));
        assert_eq!(get("wrong_measurement_rate"), Value(Some(0.25)));
        assert_eq!(get("accepted_wrong_measurement_q01"), Count(0));
        assert_eq!(get("accepted_wrong_measurement_rate_q01"), Value(None));
        assert_eq!(get("accepted_wrong_measurement_at_threshold"), Count(3));
        assert_eq!(
            get("accepted_wrong_measurement_rate_at_threshold"),
            Value(Some(0.75))
        );
        assert_eq!(
            names.last(),
            Some(&"accepted_wrong_measurement_rate_at_threshold")
        );
    }
}
