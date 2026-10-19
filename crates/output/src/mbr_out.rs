//! Match-between-runs outputs.
//!
//! `<base>.mbr.parquet`: long format, one row per (precursor, run) target that
//! was identified in the run or transferred into it, every evidence type
//! included:
//!
//! | column | content |
//! |---|---|
//! | `precursor_key` | peptidoform and charge key |
//! | `run` | run stem |
//! | `evidence_type` | `identified`, `identified_below_threshold` (identified, quality q-value above the identified threshold or missing), `transfer_ms2` (transfer q-value at or below the threshold, and an MS2 scan of the run supports the transfer), `transfer_ms1` (transfer q-value at or below the threshold, without that MS2 support) or `below_threshold` |
//! | `quantity`, `quantity_se` | fitted MS1 quantity and its standard error |
//! | `quantity_normalized` | `quantity` with the run's normalisation offset at the row's apex RT removed (a level per run and an RT-local term fitted across runs); null when the run is not normalised. `quantity` stays the raw measurement |
//! | `transfer_q` | internal transfer quant-quality q-value (not an identification FDR); null for identified rows |
//! | `quant_q` | internal quality q-value of the row's MS1 feature against its decoys (not an identification FDR): the identified-row estimate for identified rows, `transfer_q` for transfers |
//! | `apex_rt_s` | fitted apex, seconds |
//! | `library_q` | internal library entry gate (not an identification FDR) |
//!
//! The three text columns are UTF-8, the others float64; missing or
//! non-finite values are null (`library_q` is never null). A
//! `below_threshold` or `identified_below_threshold` row keeps its quantity
//! so a downstream tool can treat it as censored rather than missing.
//!
//! `<base>.mbr_validation.parquet`: one row with the pseudo-transfer
//! validation (acceptance rate, wrong-peak rate) and the accepted
//! transfer counts of the pooled and per-family decoy estimates, always at
//! q <= 0.01 (columns ending `_q01`) and at the run's transfer q-value
//! threshold (columns ending `_at_threshold`, equal to the `_q01` ones when the
//! threshold is 0.01), then the identified quality threshold (null when none)
//! and the identified row counts: all, with a quality q-value, at
//! `quant_q <= 0.01`, and reported as `identified_below_threshold`; then the
//! MS2 evidence of the transfer rows: per hypothesis (target, mass-shift
//! decoy, RT-shift decoy) the rows with MS2 evidence and how many of them an
//! MS2 scan supports, the pseudo-transfers likewise, and the accepted
//! transfers reported as `transfer_ms2`.
//!
//! Counts are uint64, thresholds and rates float64 (null when undefined).
//! It is a separate file because its one row of run-level counts shares no
//! column with the per-(precursor, run) rows.
//!
//! Both files carry their schema version in the key-value metadata
//! (`andes.mbr.schema`, `andes.mbr_validation.schema`; a reader rejects
//! another major) and are written to `<path>.partial` and renamed into place.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, Float64Array, Float64Builder, StringArray, StringBuilder, UInt64Array,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::versioned_parquet::{column, opt, read_table, write_table, TableVersion, BATCH_ROWS};

const MBR_TABLE: TableVersion = TableVersion {
    key: "andes.mbr.schema",
    version: "1.0",
    major: 1,
    what: "MBR table",
};
const VALIDATION_TABLE: TableVersion = TableVersion {
    key: "andes.mbr_validation.schema",
    version: "1.0",
    major: 1,
    what: "MBR validation table",
};

/// Why a quantity is reported for a (precursor, run).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceType {
    Identified,
    /// Identified, with a quality q-value above the identified threshold (or
    /// none).
    IdentifiedBelowThreshold,
    /// An accepted transfer that an MS2 scan of the run supports.
    TransferMs2,
    TransferMs1,
    BelowThreshold,
}

impl EvidenceType {
    pub fn as_str(self) -> &'static str {
        match self {
            EvidenceType::Identified => "identified",
            EvidenceType::IdentifiedBelowThreshold => "identified_below_threshold",
            EvidenceType::TransferMs2 => "transfer_ms2",
            EvidenceType::TransferMs1 => "transfer_ms1",
            EvidenceType::BelowThreshold => "below_threshold",
        }
    }
}

/// One row of `<base>.mbr.parquet`.
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
}

pub const MBR_COLUMNS: [&str; 10] = [
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
];

/// The default transfer q-value threshold, at which the validation is
/// always reported.
pub const VALIDATION_Q: f64 = 0.01;

/// Accepted counts at one transfer q-value threshold.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Accepted {
    /// Pseudo-transfers at or below the threshold.
    pub pseudo: usize,
    /// Accepted pseudo-transfers among `MbrValidation::wrong_peak_pairs`.
    pub pairs: usize,
    /// Of those, at the wrong peak.
    pub wrong_peak: usize,
    /// Transfer targets at or below the threshold, both decoy families.
    pub pooled: usize,
    /// Mass-shift decoys only.
    pub mass_only: usize,
    /// RT-shift decoys only.
    pub rt_only: usize,
    /// Transfer targets at or below the threshold by the wrong-peak decoy
    /// families alone (without the identity decoys).
    pub wrong_peak_only: usize,
    /// Decoy transfer rows at or below the threshold: mass-shift, RT-shift,
    /// identity decoys.
    pub decoys: [usize; 3],
}

/// The pseudo-transfer validation and transfer counts of one finalize.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MbrValidation {
    /// The run's transfer q-value threshold (`at_threshold`).
    pub q_threshold: f64,
    /// Pseudo-transfer rows (identified pairs measured as transfers).
    pub pseudo_transfers: usize,
    /// Pseudo-transfers with a q-value.
    pub pseudo_scored: usize,
    /// Pseudo-transfers whose identified fit exists, both fitted.
    pub wrong_peak_pairs: usize,
    /// Of those, apex further than the identified width from the identified apex.
    pub wrong_peak: usize,
    /// Transfer targets.
    pub transfers: usize,
    /// Transfer targets with a q-value.
    pub transfers_scored: usize,
    /// Accepted at [`VALIDATION_Q`].
    pub at_q01: Accepted,
    /// Accepted at `q_threshold`.
    pub at_threshold: Accepted,
    /// The identified quality q-value threshold, when one is applied.
    pub quant_feature_threshold: Option<f64>,
    /// Identified target rows.
    pub identified: usize,
    /// Identified target rows with a quality q-value.
    pub identified_scored: usize,
    /// Identified target rows with a quality q-value at or below
    /// [`VALIDATION_Q`].
    pub identified_q01: usize,
    /// Identified target rows reported as `identified_below_threshold`.
    pub identified_below_threshold: usize,
    /// Transfer rows with MS2 evidence, per hypothesis: target, mass-shift
    /// decoy, RT-shift decoy.
    pub ms2_rows: [usize; 3],
    /// Of those, the rows an MS2 scan supports.
    pub ms2_supported: [usize; 3],
    /// Pseudo-transfers with MS2 evidence.
    pub ms2_pseudo_rows: usize,
    /// Of those, the ones an MS2 scan supports.
    pub ms2_pseudo_supported: usize,
    /// Accepted transfers reported as `transfer_ms2`.
    pub transfer_ms2: usize,
    /// Target and identity-decoy library entries.
    pub target_entries: usize,
    pub identity_entries: usize,
    /// Decoy transfer rows with a q-value: mass-shift, RT-shift, identity
    /// decoys.
    pub decoy_transfers: [usize; 3],
    /// Transfer targets not transferred because their best donor PSM is
    /// above the transfer donor q-value.
    pub donor_gated: usize,
}

/// Columns of the decoy acceptance, each followed by `_q01` and by
/// `_at_threshold`, after the MS2 columns.
const DECOY_ACCEPTED_COLUMNS: [&str; 4] = [
    "accepted_wrong_peak_only",
    "accepted_mass_decoys",
    "accepted_rt_decoys",
    "accepted_identity_decoys",
];

/// Columns of the accepted counts, each followed by `_q01` and by
/// `_at_threshold`.
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

/// The columns of `<base>.mbr_validation.parquet`, in order.
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
            "donor_gated",
        ]
        .iter()
        .map(|s| s.to_string()),
    );
    c
}

fn ratio(num: usize, den: usize) -> Option<f64> {
    (den > 0).then(|| num as f64 / den as f64)
}

/// One cell of the validation table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ValidationCell {
    /// A count (uint64).
    Count(u64),
    /// A threshold or rate (float64, null when undefined).
    Value(Option<f64>),
}

impl MbrValidation {
    /// Pseudo-transfers accepted in `a` over all pseudo-transfers. Optimistic
    /// as a transfer sensitivity: the alignment map of a pseudo-transfer's
    /// run was fitted with that pair as an anchor.
    pub fn pseudo_acceptance(&self, a: &Accepted) -> Option<f64> {
        ratio(a.pseudo, self.pseudo_transfers)
    }

    pub fn wrong_peak_rate(&self) -> Option<f64> {
        ratio(self.wrong_peak, self.wrong_peak_pairs)
    }

    /// The validation row, one cell per [`validation_columns`] entry.
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
            n(self.donor_gated),
        ]);
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
    ])
}

fn mbr_batch(schema: &Arc<Schema>, rows: &[MbrRow<'_>]) -> Result<RecordBatch, String> {
    let mut key = StringBuilder::new();
    let mut run = StringBuilder::new();
    let mut ev = StringBuilder::new();
    let mut floats: Vec<Float64Builder> = (0..6).map(|_| Float64Builder::new()).collect();
    let mut library_q = Float64Builder::new();
    for r in rows {
        key.append_value(r.precursor_key);
        run.append_value(r.run);
        ev.append_value(r.evidence_type.as_str());
        let values = [
            r.quantity,
            r.quantity_se,
            r.quantity_normalized,
            r.transfer_q,
            r.quant_q,
            r.apex_rt_s,
        ];
        for (b, v) in floats.iter_mut().zip(values) {
            b.append_option(v.filter(|v| v.is_finite()));
        }
        library_q.append_value(r.library_q);
    }
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(key.finish()),
        Arc::new(run.finish()),
        Arc::new(ev.finish()),
    ];
    cols.extend(floats.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
    cols.push(Arc::new(library_q.finish()));
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

/// Writes `<base>.mbr.parquet` rows in the given order.
pub fn write_mbr_parquet(path: &Path, rows: &[MbrRow<'_>]) -> Result<(), String> {
    let schema = Arc::new(mbr_schema());
    let batches = rows.chunks(BATCH_ROWS).map(|c| mbr_batch(&schema, c));
    write_table(path, &MBR_TABLE, MBR_TABLE.version, schema.clone(), batches, &[])
}

/// One row of `<base>.mbr.parquet` read back.
#[derive(Debug, Clone, PartialEq)]
pub struct MbrRecord {
    pub precursor_key: String,
    pub run: String,
    /// One of the [`EvidenceType::as_str`] names.
    pub evidence_type: String,
    pub quantity: Option<f64>,
    pub quantity_se: Option<f64>,
    pub quantity_normalized: Option<f64>,
    pub transfer_q: Option<f64>,
    pub quant_q: Option<f64>,
    pub apex_rt_s: Option<f64>,
    pub library_q: f64,
}

/// Reads `<base>.mbr.parquet`, in file order. Rejects another major schema
/// version; every error names the file.
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
            });
        }
    }
    Ok(out)
}

/// Writes the one-row `<base>.mbr_validation.parquet`.
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

/// Reads `<base>.mbr_validation.parquet`: every column with its value, in
/// file order. Rejects another major schema version; every error names the
/// file.
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
            DataType::Float64 => ValidationCell::Value(opt(column::<Float64Array>(b, name, path)?, 0)),
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
            })
            .collect();
        assert_eq!(back, want);
        assert_eq!(back[1].quantity, None);
        // The column names and order are the documented ones.
        let t = read_table(&p, &MBR_TABLE).unwrap();
        let names: Vec<String> = t.batches[0]
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect();
        assert_eq!(names, MBR_COLUMNS);
        assert_eq!(t.meta.get(MBR_TABLE.key).map(String::as_str), Some("1.0"));
        // An empty table is still a valid file.
        write_mbr_parquet(&p, &[]).unwrap();
        assert!(read_mbr_parquet(&p).unwrap().is_empty());
    }

    #[test]
    fn mbr_readers_reject_other_majors_and_each_others_files() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr.parquet");
        let schema = Arc::new(mbr_schema());
        let b = mbr_batch(&schema, &rows());
        write_table(&p, &MBR_TABLE, "2.0", schema, [b], &[]).unwrap();
        let e = read_mbr_parquet(&p).unwrap_err();
        assert!(e.contains("x.mbr.parquet") && e.contains("schema 2.0"), "{e}");
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
            transfers: 100,
            transfers_scored: 90,
            at_q01: Accepted {
                pseudo: 8,
                pairs: 0,
                wrong_peak: 0,
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
            donor_gated: 3,
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
        assert_eq!(get("accepted_wrong_peak_rate_at_threshold"), Value(Some(0.25)));
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
        assert_eq!(get("donor_gated"), Count(3));
    }
}
