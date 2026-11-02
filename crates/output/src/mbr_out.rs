//! Match-between-runs outputs.
//!
//! `<base>.mbr.parquet`: long format, one row per (precursor, run) target that
//! was identified in the run or transferred into it, every evidence type
//! included:
//!
//! | column | content |
//! |---|---|
//! | `precursor_key` | `<peptidoform>/<charge>` of the precursor, a positional-isomer class, in the site form of its best PSM |
//! | `run` | run stem |
//! | `evidence_type` | `identified`, `identified_below_threshold` (identified, quality q-value above the identified threshold or missing), `transfer_ms2` (transfer q-value at or below the threshold, and an MS2 scan of the run supports the transfer), `transfer_ms1` (transfer q-value at or below the threshold, without that MS2 support) or `below_threshold` |
//! | `quantity`, `quantity_se` | fitted MS1 quantity and its standard error |
//! | `quantity_normalized` | `quantity` with the run's normalisation offset at the row's apex RT removed (a level per run and an RT-local term fitted across runs); null when the run is not normalised. `quantity` stays the raw measurement |
//! | `transfer_q` | internal transfer quant-quality q-value (not an identification FDR); null for identified rows |
//! | `quant_q` | internal quality q-value of the row's MS1 feature against its decoys (not an identification FDR): the identified-row estimate for identified rows; for an accepted transfer the running mean, over the run's accepted transfers in ascending probability of a wrong measurement, of that probability; `transfer_q` for a transfer below the threshold |
//! | `apex_rt_s` | fitted apex, seconds |
//! | `library_q` | internal library entry gate (not an identification FDR) |
//! | `shared_peak_with` | for a transfer, the key of a look-alike precursor (m/z within the extract tolerance) identified in the run whose identified peak the transfer's apex lies on (within its identified width); null otherwise. Such a transfer may measure the look-alike rather than its own precursor (schema 1.1) |
//! | `quantity_fixed`, `quantity_fixed_se` | for a transfer, the quantity with the apex fixed at the predicted RT and the width at the prior's (the same width as `quantity`), not clamped at 0, and its standard error: near 0 on average when the precursor is absent, where `quantity` takes the largest bump near the prediction; null for identified rows and for an extract of schema 4.0 (schema 1.1) |
//!
//! The text columns are UTF-8, the others float64; missing or
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
//! transfers reported as `transfer_ms2`; then (schema 1.1) the transfer
//! targets whose apex lies on a look-alike's identified peak, how many of
//! them are accepted at the threshold and how many of those an MS2 scan
//! supports that tells the two apart, and the pseudo-transfers on a
//! look-alike's peak.
//!
//! `<base>.mbr_validation_runs.parquet` ([`write_mbr_runs`]): the transfer
//! counts per (run, prior stratum), with the stratum's range of the prior,
//! the K-decoy rate and the mean local FDR of each cell's accepted
//! transfers, at q <= 0.01 and at the threshold; the cells add up to the
//! counts of the validation table.
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

use crate::versioned_parquet::{
    column, opt, read_table, write_table, TableVersion, TableWriter, BATCH_ROWS,
};

const MBR_TABLE: TableVersion = TableVersion {
    key: "andes.mbr.schema",
    version: "1.1",
    major: 1,
    what: "MBR table",
};
const VALIDATION_TABLE: TableVersion = TableVersion {
    key: "andes.mbr_validation.schema",
    version: "1.1",
    major: 1,
    what: "MBR validation table",
};

const RUNS_TABLE: TableVersion = TableVersion {
    key: "andes.mbr_validation_runs.schema",
    version: "1.0",
    major: 1,
    what: "MBR per-run validation table",
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
    /// The key of the look-alike precursor whose identified peak this
    /// transfer lies on.
    pub shared_peak_with: Option<&'a str>,
    /// A transfer's quantity with the apex fixed at its predicted RT, not
    /// clamped, and its standard error.
    pub quantity_fixed: Option<f64>,
    pub quantity_fixed_se: Option<f64>,
}

pub const MBR_COLUMNS: [&str; 13] = [
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
    "shared_peak_with",
    "quantity_fixed",
    "quantity_fixed_se",
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
    /// Transfer targets whose apex lies on the identified peak of a
    /// look-alike precursor in the run.
    pub shared_peak_transfers: usize,
    /// Of those, accepted at `q_threshold`.
    pub shared_peak_accepted: usize,
    /// Of those, with an MS2 scan that supports the transfer and tells it
    /// apart from the look-alike (not an I/L swap).
    pub shared_peak_accepted_ms2: usize,
    /// Pseudo-transfers (identified, so present) on a look-alike's peak.
    pub shared_peak_pseudo: usize,
}

fn add_counts<const N: usize>(a: &mut [usize; N], b: &[usize; N]) {
    for (x, y) in a.iter_mut().zip(b) {
        *x += y;
    }
}

impl Accepted {
    /// Adds the counts of `o` (of disjoint rows).
    pub fn add(&mut self, o: &Accepted) {
        let Accepted {
            pseudo,
            pairs,
            wrong_peak,
            pooled,
            mass_only,
            rt_only,
            wrong_peak_only,
            decoys,
        } = o;
        self.pseudo += pseudo;
        self.pairs += pairs;
        self.wrong_peak += wrong_peak;
        self.pooled += pooled;
        self.mass_only += mass_only;
        self.rt_only += rt_only;
        self.wrong_peak_only += wrong_peak_only;
        add_counts(&mut self.decoys, decoys);
    }
}

impl MbrValidation {
    /// Adds the row counts of `o` (of disjoint rows); the thresholds and the
    /// library entry counts stay.
    pub fn add_counts(&mut self, o: &MbrValidation) {
        let MbrValidation {
            q_threshold: _,
            pseudo_transfers,
            pseudo_scored,
            wrong_peak_pairs,
            wrong_peak,
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
            donor_gated,
            shared_peak_transfers,
            shared_peak_accepted,
            shared_peak_accepted_ms2,
            shared_peak_pseudo,
        } = o;
        self.pseudo_transfers += pseudo_transfers;
        self.pseudo_scored += pseudo_scored;
        self.wrong_peak_pairs += wrong_peak_pairs;
        self.wrong_peak += wrong_peak;
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
        self.donor_gated += donor_gated;
        self.shared_peak_transfers += shared_peak_transfers;
        self.shared_peak_accepted += shared_peak_accepted;
        self.shared_peak_accepted_ms2 += shared_peak_accepted_ms2;
        self.shared_peak_pseudo += shared_peak_pseudo;
    }
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
            "shared_peak_transfers",
            "shared_peak_accepted_at_threshold",
            "shared_peak_accepted_ms2_at_threshold",
            "shared_peak_pseudo",
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
            n(self.shared_peak_transfers),
            n(self.shared_peak_accepted),
            n(self.shared_peak_accepted_ms2),
            n(self.shared_peak_pseudo),
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
        Field::new("shared_peak_with", DataType::Utf8, true),
        f("quantity_fixed"),
        f("quantity_fixed_se"),
    ])
}

/// The columns of `<base>.mbr.parquet` being filled row by row.
struct MbrColumnsBuilder {
    key: StringBuilder,
    run: StringBuilder,
    ev: StringBuilder,
    floats: Vec<Float64Builder>,
    library_q: Float64Builder,
    shared: StringBuilder,
    fixed: [Float64Builder; 2],
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
            shared: StringBuilder::new(),
            fixed: [Float64Builder::new(), Float64Builder::new()],
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
        self.shared.append_option(r.shared_peak_with);
        for (b, v) in self
            .fixed
            .iter_mut()
            .zip([r.quantity_fixed, r.quantity_fixed_se])
        {
            b.append_option(v.filter(|v| v.is_finite()));
        }
        self.rows += 1;
    }

    /// The batch of the rows pushed since the last call.
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
        cols.push(Arc::new(self.shared.finish()));
        cols.extend(
            self.fixed
                .iter_mut()
                .map(|b| Arc::new(b.finish()) as ArrayRef),
        );
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

/// `<base>.mbr.parquet` written row by row, in batches of the same size as
/// [`write_mbr_parquet`]'s (the same rows give the same file).
pub struct MbrParquetWriter {
    schema: Arc<Schema>,
    columns: MbrColumnsBuilder,
    table: TableWriter,
    path: std::path::PathBuf,
}

impl MbrParquetWriter {
    /// Creates `<path>.partial`; [`MbrParquetWriter::finish`] renames it into
    /// place.
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

/// Writes `<base>.mbr.parquet` rows in the given order.
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
    /// `None` also for a schema-1.0 file, which has no such column.
    pub shared_peak_with: Option<String>,
    /// `None` also for a schema-1.0 file.
    pub quantity_fixed: Option<f64>,
    pub quantity_fixed_se: Option<f64>,
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
        let shared = match b.schema().index_of("shared_peak_with") {
            Ok(_) => Some(text("shared_peak_with")?),
            Err(_) => None,
        };
        let optional_float = |name: &str| match b.schema().index_of(name) {
            Ok(_) => float(name).map(Some),
            Err(_) => Ok(None),
        };
        let (fq, fse) = (
            optional_float("quantity_fixed")?,
            optional_float("quantity_fixed_se")?,
        );
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
                shared_peak_with: shared
                    .filter(|c| !c.is_null(i))
                    .map(|c| c.value(i).to_string()),
                quantity_fixed: fq.and_then(|c| opt(c, i)),
                quantity_fixed_se: fse.and_then(|c| opt(c, i)),
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

/// The transfer counts of one (run, prior stratum) of
/// `<base>.mbr_validation_runs.parquet`. Index `[0]` of a pair is at
/// q <= [`VALIDATION_Q`], `[1]` at the run's transfer q-value threshold;
/// decoy families are mass shift, RT shift, identity decoy.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RunStratumCounts {
    /// The stratum's range of the prior of presence: above `prior_lo` (from
    /// it, for the first stratum) up to `prior_hi`.
    pub prior_lo: f64,
    pub prior_hi: f64,
    /// Transfer targets with a q-value.
    pub transfers: usize,
    /// Decoy transfer rows with a q-value, per family.
    pub decoys: [usize; 3],
    /// Accepted transfer targets.
    pub accepted: [usize; 2],
    /// Accepted decoy transfer rows, per family.
    pub accepted_decoys: [[usize; 3]; 2],
    /// Sum of the local FDR over the accepted transfer targets: their
    /// expected false transfers.
    pub expected_false: [f64; 2],
    /// Accepted transfer targets on a look-alike's identified peak (at the
    /// threshold).
    pub shared_peak_accepted: usize,
}

impl RunStratumCounts {
    /// Adds the counts of `o` (of disjoint rows); the prior range stays.
    pub fn add(&mut self, o: &RunStratumCounts) {
        self.transfers += o.transfers;
        add_counts(&mut self.decoys, &o.decoys);
        add_counts(&mut self.accepted, &o.accepted);
        for (a, b) in self.accepted_decoys.iter_mut().zip(&o.accepted_decoys) {
            add_counts(a, b);
        }
        for (a, b) in self.expected_false.iter_mut().zip(&o.expected_false) {
            *a += b;
        }
        self.shared_peak_accepted += o.shared_peak_accepted;
    }

    /// The wrong-peak decoys per target of this cell (a target and its
    /// decoys share their run and stratum).
    pub fn decoy_ratio(&self) -> Option<f64> {
        ratio(self.decoys[0] + self.decoys[1], self.transfers)
    }

    /// The K-decoy rate of the cell's accepted transfers at threshold `i`
    /// (0: q 0.01, 1: the run's threshold):
    /// `((D_mass + D_rt) / K + D_identity) / T`, `K` the cell's
    /// [`RunStratumCounts::decoy_ratio`]. Noisy for a small cell; the sum
    /// over cells of the numerators and of `T` gives a coarser one.
    pub fn decoy_rate(&self, i: usize) -> Option<f64> {
        let k = self.decoy_ratio().filter(|k| *k > 0.0)?;
        let d = &self.accepted_decoys[i];
        let t = self.accepted[i];
        (t > 0).then(|| ((d[0] + d[1]) as f64 / k + d[2] as f64) / t as f64)
    }

    /// The local-FDR rate of the cell's accepted transfers at threshold `i`.
    pub fn model_rate(&self, i: usize) -> Option<f64> {
        let t = self.accepted[i];
        (t > 0).then(|| self.expected_false[i] / t as f64)
    }
}

/// The columns of `<base>.mbr_validation_runs.parquet`, in order.
pub const RUNS_COLUMNS: [&str; 23] = [
    "run",
    "prior_stratum",
    "prior_lo",
    "prior_hi",
    "transfers",
    "mass_decoys",
    "rt_decoys",
    "identity_decoys",
    "decoy_ratio",
    "accepted_q01",
    "accepted_mass_decoys_q01",
    "accepted_rt_decoys_q01",
    "accepted_identity_decoys_q01",
    "decoy_rate_q01",
    "expected_false_q01",
    "model_rate_q01",
    "accepted_at_threshold",
    "accepted_mass_decoys_at_threshold",
    "accepted_rt_decoys_at_threshold",
    "accepted_identity_decoys_at_threshold",
    "decoy_rate_at_threshold",
    "expected_false_at_threshold",
    "model_rate_at_threshold",
];

/// The column after [`RUNS_COLUMNS`]: accepted shared-peak transfers.
const RUNS_SHARED: &str = "shared_peak_accepted_at_threshold";

/// Writes `<base>.mbr_validation_runs.parquet`: one row per (run, prior
/// stratum) in the order given, with [`RUNS_COLUMNS`] and the accepted
/// transfers on a look-alike's peak. Counts are uint64, the prior range,
/// ratios, expected false transfers and rates float64 (null when
/// undefined).
pub fn write_mbr_runs(path: &Path, rows: &[(&str, u8, RunStratumCounts)]) -> Result<(), String> {
    let mut fields = vec![
        Field::new("run", DataType::Utf8, false),
        Field::new("prior_stratum", DataType::UInt64, false),
    ];
    let float = |n: &str| Field::new(n, DataType::Float64, true);
    let count = |n: &str| Field::new(n, DataType::UInt64, false);
    for &n in &RUNS_COLUMNS[2..] {
        let is_float = n.starts_with("prior_")
            || n.ends_with("ratio")
            || n.contains("_rate_")
            || n.starts_with("expected_false");
        fields.push(if is_float { float(n) } else { count(n) });
    }
    fields.push(count(RUNS_SHARED));
    let schema = Arc::new(Schema::new(fields));
    let run: Vec<&str> = rows.iter().map(|r| r.0).collect();
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(run)),
        Arc::new(UInt64Array::from(
            rows.iter().map(|r| u64::from(r.1)).collect::<Vec<_>>(),
        )),
    ];
    let f = |g: &dyn Fn(&RunStratumCounts) -> Option<f64>| -> ArrayRef {
        Arc::new(Float64Array::from(
            rows.iter()
                .map(|r| g(&r.2).filter(|v| v.is_finite()))
                .collect::<Vec<_>>(),
        ))
    };
    let n = |g: &dyn Fn(&RunStratumCounts) -> usize| -> ArrayRef {
        Arc::new(UInt64Array::from(
            rows.iter().map(|r| g(&r.2) as u64).collect::<Vec<_>>(),
        ))
    };
    cols.push(f(&|c| Some(c.prior_lo)));
    cols.push(f(&|c| Some(c.prior_hi)));
    cols.push(n(&|c| c.transfers));
    for k in 0..3 {
        cols.push(n(&|c| c.decoys[k]));
    }
    cols.push(f(&|c| c.decoy_ratio()));
    for i in 0..2 {
        cols.push(n(&|c| c.accepted[i]));
        for k in 0..3 {
            cols.push(n(&|c| c.accepted_decoys[i][k]));
        }
        cols.push(f(&|c| c.decoy_rate(i)));
        cols.push(f(&|c| Some(c.expected_false[i])));
        cols.push(f(&|c| c.model_rate(i)));
    }
    cols.push(n(&|c| c.shared_peak_accepted));
    let batch = RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string());
    write_table(path, &RUNS_TABLE, RUNS_TABLE.version, schema, [batch], &[])
}

/// Reads `<base>.mbr_validation_runs.parquet` back: per row, the run, the
/// stratum and the counts. Rejects another major schema version; every
/// error names the file.
pub fn read_mbr_runs(path: &Path) -> Result<Vec<(String, u8, RunStratumCounts)>, String> {
    let t = read_table(path, &RUNS_TABLE)?;
    let mut out = Vec::new();
    for b in &t.batches {
        let run = column::<StringArray>(b, "run", path)?;
        let n = |name: &str| column::<UInt64Array>(b, name, path);
        let f = |name: &str| column::<Float64Array>(b, name, path);
        let stratum = n("prior_stratum")?;
        let (lo, hi) = (f("prior_lo")?, f("prior_hi")?);
        let transfers = n("transfers")?;
        let decoys = [n("mass_decoys")?, n("rt_decoys")?, n("identity_decoys")?];
        let at = |suffix: &str| -> Result<_, String> {
            Ok((
                n(&format!("accepted_{suffix}"))?,
                [
                    n(&format!("accepted_mass_decoys_{suffix}"))?,
                    n(&format!("accepted_rt_decoys_{suffix}"))?,
                    n(&format!("accepted_identity_decoys_{suffix}"))?,
                ],
                f(&format!("expected_false_{suffix}"))?,
            ))
        };
        let (a0, ad0, e0) = at("q01")?;
        let (a1, ad1, e1) = at("at_threshold")?;
        let shared = n(RUNS_SHARED)?;
        for i in 0..b.num_rows() {
            let v = |a: &UInt64Array| a.value(i) as usize;
            out.push((
                run.value(i).to_string(),
                stratum.value(i) as u8,
                RunStratumCounts {
                    prior_lo: opt(lo, i).unwrap_or(f64::NAN),
                    prior_hi: opt(hi, i).unwrap_or(f64::NAN),
                    transfers: v(transfers),
                    decoys: decoys.map(v),
                    accepted: [v(a0), v(a1)],
                    accepted_decoys: [ad0.map(v), ad1.map(v)],
                    expected_false: [
                        opt(e0, i).unwrap_or(f64::NAN),
                        opt(e1, i).unwrap_or(f64::NAN),
                    ],
                    shared_peak_accepted: v(shared),
                },
            ));
        }
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
        // More than one batch, with a partial last one.
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
        // Dropped unfinished, no file is left.
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
                shared_peak_with: None,
                quantity_fixed: None,
                quantity_fixed_se: None,
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
                shared_peak_with: Some("PEPTLDEK/2"),
                quantity_fixed: Some(-1.5e3),
                quantity_fixed_se: Some(2e3),
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
                shared_peak_with: r.shared_peak_with.map(String::from),
                quantity_fixed: r.quantity_fixed,
                quantity_fixed_se: r.quantity_fixed_se,
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
        assert_eq!(t.meta.get(MBR_TABLE.key).map(String::as_str), Some("1.1"));
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
            shared_peak_transfers: 12,
            shared_peak_accepted: 6,
            shared_peak_accepted_ms2: 1,
            shared_peak_pseudo: 2,
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
        assert_eq!(get("donor_gated"), Count(3));
        assert_eq!(get("shared_peak_transfers"), Count(12));
        assert_eq!(get("shared_peak_accepted_at_threshold"), Count(6));
        assert_eq!(get("shared_peak_accepted_ms2_at_threshold"), Count(1));
        assert_eq!(get("shared_peak_pseudo"), Count(2));
    }

    #[test]
    fn runs_table_round_trips_with_its_rates() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr_validation_runs.parquet");
        let a = RunStratumCounts {
            prior_lo: 0.0,
            prior_hi: 0.05,
            transfers: 100,
            decoys: [100, 100, 3],
            accepted: [40, 50],
            accepted_decoys: [[1, 1, 0], [2, 2, 1]],
            expected_false: [0.5, 1.25],
            shared_peak_accepted: 4,
        };
        let empty = RunStratumCounts {
            prior_lo: 0.05,
            prior_hi: 1.0,
            ..RunStratumCounts::default()
        };
        let rows = vec![("a", 0u8, a.clone()), ("a", 1u8, empty.clone())];
        write_mbr_runs(&p, &rows).unwrap();
        let back = read_mbr_runs(&p).unwrap();
        assert_eq!(
            back,
            vec![
                ("a".to_string(), 0u8, a.clone()),
                ("a".to_string(), 1u8, empty.clone())
            ]
        );
        // K = 2 decoys per target: (2 / 2 + 0) / 40 at q 0.01, (4 / 2 + 1) / 50.
        assert_eq!(a.decoy_ratio(), Some(2.0));
        assert_eq!(a.decoy_rate(0), Some(1.0 / 40.0));
        assert_eq!(a.decoy_rate(1), Some(3.0 / 50.0));
        assert_eq!(a.model_rate(1), Some(1.25 / 50.0));
        assert_eq!((empty.decoy_rate(0), empty.model_rate(0)), (None, None));
        let t = read_table(&p, &RUNS_TABLE).unwrap();
        let names: Vec<String> = t.batches[0]
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect();
        let mut want: Vec<&str> = RUNS_COLUMNS.to_vec();
        want.push(RUNS_SHARED);
        assert_eq!(names, want);
        let mut sum = a.clone();
        sum.add(&a);
        assert_eq!((sum.transfers, sum.accepted_decoys[1][2]), (200, 2));
        assert_eq!(sum.expected_false, [1.0, 2.5]);
    }

    #[test]
    fn a_schema_1_0_mbr_table_reads_without_shared_peaks() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("old.mbr.parquet");
        let new = mbr_schema();
        let old = Arc::new(Schema::new(
            new.fields()
                .iter()
                .filter(|f| {
                    !["shared_peak_with", "quantity_fixed", "quantity_fixed_se"]
                        .contains(&f.name().as_str())
                })
                .cloned()
                .collect::<Vec<_>>(),
        ));
        let full = mbr_batch(&Arc::new(new), &rows()).unwrap();
        let cols: Vec<ArrayRef> = (0..old.fields().len())
            .map(|i| full.column(i).clone())
            .collect();
        let b = RecordBatch::try_new(old.clone(), cols).map_err(|e| e.to_string());
        write_table(&p, &MBR_TABLE, "1.0", old, [b], &[]).unwrap();
        let back = read_mbr_parquet(&p).unwrap();
        assert_eq!(back.len(), 2);
        assert!(back
            .iter()
            .all(|r| r.shared_peak_with.is_none() && r.quantity_fixed.is_none()));
    }
}
