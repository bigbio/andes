//! Parquet reader: load [`Param`] models from the combined-schema Parquet store.
//!
//! Usage:
//! ```rust,ignore
//! let store = ModelStore::open(&path)?;
//! let param = store.load_param("my_model")?;
//! ```

use std::path::{Path, PathBuf};

use arrow::array::{
    Array, BinaryArray, BooleanArray, Float32Array, Int32Array, Int64Array, ListArray, StringArray,
    StructArray,
};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rustc_hash::FxHashMap;

use model::tolerance::Tolerance;
use scoring_crate::param_model::{
    FragmentOffsetFrequency, IonType, Param, Partition, PrecursorOffsetFrequency, SpecDataType,
};

use crate::counts::CountStats;
use crate::select::{parse_experiment_class, SelectionEntry};
use crate::store::write::SourceLedger;
use crate::TrainError;

// ── public API ───────────────────────────────────────────────────────────────

/// A handle to a Parquet model store.  Created by [`ModelStore::open`].
pub struct ModelStore {
    /// One or more Parquet part files that make up this store.
    ///
    /// - Single-file store: a one-element vector holding the file path.
    /// - Partitioned store (a directory of `protocol=<P>/models.parquet`
    ///   partitions): one entry per partition file, sorted for determinism.
    ///
    /// All read operations scan every part and union the results, so a
    /// partitioned store is indistinguishable from the equivalent single file.
    parts: Vec<PathBuf>,
    /// Full manifest rows read at open time (union across all parts).
    manifest: Vec<RawManifestEntry>,
}

/// Raw manifest entry as stored in the Parquet file (before type conversion).
#[derive(Debug, Clone)]
pub struct RawManifestEntry {
    /// Lowercased model identifier, e.g. `"hcd_qexactive_tryp"`.
    pub model_id: String,
    /// Activation method string as stored, e.g. `"HCD"`, `"CID"`.
    pub activation: String,
    /// Instrument string as stored, e.g. `"QExactive"`, `"LowRes"`.
    pub instrument: String,
    /// Enzyme string as stored, e.g. `"Trypsin"`, `"ArgC"`.
    pub enzyme: String,
    /// Protocol string as stored, e.g. `"Automatic"`, `"TMT"`, `"Phosphorylation"`.
    pub protocol: String,
}

impl ModelStore {
    /// Open the store and read the manifest.
    ///
    /// `path` may be either:
    /// - a single Parquet file (the historical single-file layout), or
    /// - a **directory** of per-protocol Parquet partitions written in the
    ///   Hive-style layout `protocol=<Protocol>/models.parquet` (or, more
    ///   generally, any tree of `*.parquet` files). Every partition is read and
    ///   the manifests/tables are unioned into one in-memory store, exactly as
    ///   if they had been written to a single file.
    pub fn open(path: &Path) -> Result<Self, TrainError> {
        let parts = resolve_parts(path)?;
        // Build a model_id → part index as we read manifests (finding 3.7). In a
        // partitioned store every model must live in exactly one part; a
        // duplicate `model_id` across parts (e.g. a stale temp partition left
        // behind) otherwise makes `load_param` path-order-dependent — it stops at
        // the first part that contains the id. Reject the ambiguity up front.
        let mut manifest: Vec<RawManifestEntry> = Vec::new();
        let mut owner: FxHashMap<String, PathBuf> = FxHashMap::default();
        for part in &parts {
            for entry in read_manifest(part)? {
                if let Some(prev) = owner.get(&entry.model_id) {
                    if prev != part {
                        return Err(TrainError::Other(format!(
                            "duplicate model_id '{}' found in two partitions ({} and {}); \
                             a partitioned store must hold each model in exactly one part \
                             (remove the stale/duplicate partition)",
                            entry.model_id,
                            prev.display(),
                            part.display(),
                        )));
                    }
                } else {
                    owner.insert(entry.model_id.clone(), part.clone());
                }
                manifest.push(entry);
            }
        }
        Ok(Self { parts, manifest })
    }

    /// List all model IDs in the store.
    pub fn model_ids(&self) -> Vec<String> {
        self.manifest.iter().map(|e| e.model_id.clone()).collect()
    }

    /// Return the raw manifest entries (one per model).
    pub fn manifest_entries(&self) -> &[RawManifestEntry] {
        &self.manifest
    }

    /// Convert each manifest entry to a [`SelectionEntry`] for use with
    /// [`crate::select::select`].
    ///
    /// The mapping from the parquet `protocol` column to `experiment_class` is:
    /// - `"Automatic"` or `"Standard"` → empty set
    /// - `"TMT"` → `{"tmt"}`
    /// - `"Phosphorylation"` → `{"phospho"}`
    /// - `"iTRAQ"` → `{"itraq"}`
    /// - `"iTRAQPhospho"` → `{"itraq", "phospho"}`
    ///
    /// Any unrecognised protocol falls back to empty set (same as "standard").
    pub fn selection_entries(&self) -> Vec<SelectionEntry> {
        self.manifest
            .iter()
            .map(|e| SelectionEntry {
                model_id: e.model_id.clone(),
                activation: e.activation.clone(),
                instrument: e.instrument.clone(),
                enzyme: e.enzyme.clone(),
                experiment_class: protocol_to_experiment_class(&e.protocol),
            })
            .collect()
    }

    /// Load and reconstruct a [`Param`] for the given `model_id`.
    ///
    /// In a partitioned store, each model lives in exactly one partition; the
    /// reconstructor scans every part and stops once it has assembled the model.
    pub fn load_param(&self, model_id: &str) -> Result<Param, TrainError> {
        reconstruct_param(&self.parts, model_id)
    }

    /// Return the source ledger entries for the given `model_id`.
    ///
    /// Returns an empty `Vec` if the store contains no `"source"` rows for
    /// this model (e.g. stores written by the legacy [`super::write_models`]).
    pub fn load_sources(&self, model_id: &str) -> Result<Vec<SourceLedger>, TrainError> {
        let mut ledgers: Vec<SourceLedger> = Vec::new();
        for part in &self.parts {
            ledgers.extend(read_sources(part, model_id)?);
        }
        Ok(ledgers)
    }

    /// Reconstruct the [`CountStats`] for `(model_id, source_id)` from the
    /// `"stat"` rows in the store.
    ///
    /// Returns [`TrainError::NoModel`] if no stat rows are found for the
    /// given `(model_id, source_id)` pair.
    pub fn load_source_stats(
        &self,
        model_id: &str,
        source_id: &str,
    ) -> Result<CountStats, TrainError> {
        // A model's stat rows live in exactly one partition. Scan parts and
        // return the first part that actually contains stats for this pair.
        // Only a `NoModel` miss is a reason to keep probing the next partition;
        // any other error (Parquet/schema) is real and must surface immediately.
        for part in &self.parts {
            match read_source_stats(part, model_id, source_id) {
                Ok(stats) => return Ok(stats),
                Err(TrainError::NoModel(_)) => continue,
                Err(e) => return Err(e),
            }
        }
        Err(TrainError::NoModel(format!(
            "source_stats({model_id}, {source_id})"
        )))
    }
}

// ── partition resolution ──────────────────────────────────────────────────────

/// Resolve `path` to the list of Parquet part files that constitute the store.
///
/// - A regular file ⇒ a single-element list (single-file store).
/// - A directory ⇒ every `*.parquet` file found beneath it (recursively),
///   sorted for deterministic iteration order. This covers the Hive-style
///   `protocol=<P>/models.parquet` partition layout.
fn resolve_parts(path: &Path) -> Result<Vec<PathBuf>, TrainError> {
    if path.is_dir() {
        let mut parts: Vec<PathBuf> = Vec::new();
        collect_parquet_files(path, &mut parts)?;
        if parts.is_empty() {
            return Err(TrainError::Other(format!(
                "no *.parquet partition files found under directory {}",
                path.display()
            )));
        }
        parts.sort();
        Ok(parts)
    } else {
        Ok(vec![path.to_owned()])
    }
}

/// Recursively collect `*.parquet` files under `dir` into `out`.
///
/// Temporary / in-flight artifacts are ignored (finding 3.7) so a partial write
/// left behind by an interrupted store update doesn't shadow or duplicate a real
/// partition: files whose name starts with `.` (hidden / atomic-write temp) or
/// ends with `.tmp.parquet`, and any `.crc`-style sidecars (non-`.parquet`
/// extension) are skipped.
fn collect_parquet_files(dir: &Path, out: &mut Vec<PathBuf>) -> Result<(), TrainError> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let p = entry.path();
        let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if p.is_dir() {
            // Skip hidden / temp staging dirs (e.g. `.tmp`, `_temporary`) before
            // descending so a stale partial write never contributes parquet files.
            if name.starts_with('.') || name.starts_with('_') {
                continue;
            }
            collect_parquet_files(&p, out)?;
            continue;
        }
        if p.extension().and_then(|e| e.to_str()) != Some("parquet") {
            continue;
        }
        // Skip hidden/temp partition files (e.g. `.models.parquet.tmp.parquet`,
        // `.~tmp.parquet`, or any dotfile) so a stale temp doesn't become a
        // phantom duplicate partition.
        if name.starts_with('.') || name.ends_with(".tmp.parquet") {
            continue;
        }
        out.push(p);
    }
    Ok(())
}

/// Convert the parquet `protocol` column value to a `BTreeSet<String>` experiment class.
///
/// Maps the stored `protocol` column naming to the lowercase slug-set used by
/// [`crate::select`].
///
/// `iTRAQPhospho` is intentionally mapped to the single slug `"itraqphospho"`
/// (not `{"itraq","phospho"}`) so that:
/// - exact-match (step 1) in [`crate::select::select`] finds combo models
///   like `hcd_qexactive_tryp_itraqphospho` when a `SelectionKey` with
///   `{"itraqphospho"}` is used, and
/// - when no iTRAQPhospho-specific model is bundled, select() falls through
///   to the empty-class fallback (the protocol-less model) rather than
///   spuriously matching a `{"phospho"}` subset entry.
pub fn protocol_to_experiment_class(protocol: &str) -> std::collections::BTreeSet<String> {
    match protocol {
        "Automatic" | "Standard" => std::collections::BTreeSet::new(),
        "TMT"                    => parse_experiment_class("tmt"),
        "Phosphorylation"        => parse_experiment_class("phospho"),
        "iTRAQ"                  => parse_experiment_class("itraq"),
        // Keep as a single opaque slug — do NOT split into {"itraq","phospho"}.
        "iTRAQPhospho"           => {
            let mut s = std::collections::BTreeSet::new();
            s.insert("itraqphospho".to_string());
            s
        }
        other                    => parse_experiment_class(other),
    }
}

// ── manifest reader ──────────────────────────────────────────────────────────

fn read_manifest(path: &Path) -> Result<Vec<RawManifestEntry>, TrainError> {
    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut entries: Vec<RawManifestEntry> = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;
        let record_kind = batch
            .column_by_name("record_kind")
            .ok_or_else(|| TrainError::Other("missing record_kind column".into()))?;
        let model_id_col = batch
            .column_by_name("model_id")
            .ok_or_else(|| TrainError::Other("missing model_id column".into()))?;

        let rk = record_kind.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
            TrainError::Other("record_kind not StringArray".into())
        })?;
        let mid = model_id_col.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
            TrainError::Other("model_id not StringArray".into())
        })?;

        for i in 0..batch.num_rows() {
            if !rk.is_null(i) && rk.value(i) == "manifest" {
                let activation = str_col(&batch, "activation")
                    .map(|c| if c.is_null(i) { String::new() } else { c.value(i).to_string() })
                    .unwrap_or_default();
                let instrument = str_col(&batch, "instrument")
                    .map(|c| if c.is_null(i) { String::new() } else { c.value(i).to_string() })
                    .unwrap_or_default();
                let enzyme = str_col(&batch, "enzyme")
                    .map(|c| if c.is_null(i) { String::new() } else { c.value(i).to_string() })
                    .unwrap_or_default();
                let protocol = str_col(&batch, "protocol")
                    .map(|c| if c.is_null(i) { String::new() } else { c.value(i).to_string() })
                    .unwrap_or_default();
                entries.push(RawManifestEntry {
                    model_id: mid.value(i).to_string(),
                    activation,
                    instrument,
                    enzyme,
                    protocol,
                });
            }
        }
    }
    Ok(entries)
}

// ── param reconstructor ──────────────────────────────────────────────────────

struct ManifestRow {
    activation: String,
    instrument: String,
    enzyme: Option<String>,
    protocol: String,
    version: i32,
    mme_val: f32,
    mme_is_ppm: bool,
    apply_deconv: bool,
    deconv_tol: f32,
    num_segments: i32,
    max_rank: i32,
    error_scaling_factor: i32,
    min_charge: i32,
    max_charge: i32,
    num_precursor_off: i32,
    charge_hist: Vec<(i32, i32)>,
    /// Serialized GBDT blob, if present. `None` for models written before this
    /// column was added (backward-compatible: column absent ⇒ `gbdt_peak_model = None`).
    gbdt_bytes: Option<Vec<u8>>,
    /// Serialized fragment-intensity GBDT blob, if present. `None` for models
    /// written before this column existed (backward-compatible: column absent ⇒
    /// `frag_intensity_model = None`).
    frag_intensity_bytes: Option<Vec<u8>>,
    /// Serialized rich-ion GBDT blob, if present. `None` for models
    /// written before this column existed (backward-compatible: column absent ⇒
    /// `rich_ion_model = None`).
    rich_ion_bytes: Option<Vec<u8>>,
}

fn reconstruct_param(parts: &[PathBuf], model_id: &str) -> Result<Param, TrainError> {
    // We do a full scan of each part file and collect relevant rows. For large
    // stores this is acceptable as a first implementation; predicate pushdown
    // can be added later. In a partitioned store every row for a given model
    // lives in exactly one partition, so accumulating across parts is a no-op
    // for parts that don't contain this model.

    let mut manifest_opt: Option<ManifestRow> = None;
    // partition list (in stored order)
    let mut partitions: Vec<Partition> = Vec::new();
    let mut precursor_off_map: FxHashMap<i32, Vec<PrecursorOffsetFrequency>> =
        FxHashMap::default();
    let mut frag_off_table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>> =
        FxHashMap::default();
    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
        FxHashMap::default();
    let mut ion_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
    let mut noise_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
    let mut ion_existence_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();

    for path in parts {
        // A model lives in exactly one partition; once we've assembled it from
        // a part (manifest found), no later part can contribute more rows.
        if manifest_opt.is_some() {
            break;
        }
        let file = std::fs::File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;

        let record_kind = str_col(&batch, "record_kind")?;
        let mid_col = str_col(&batch, "model_id")?;

        for i in 0..batch.num_rows() {
            if mid_col.value(i) != model_id {
                continue;
            }
            let rk = if record_kind.is_null(i) { "" } else { record_kind.value(i) };
            match rk {
                "manifest" => {
                    manifest_opt = Some(parse_manifest_row(&batch, i)?);
                }
                "table" => {
                    let table_kind = str_col(&batch, "table_kind")?;
                    if table_kind.is_null(i) { continue; }
                    let tk = table_kind.value(i);

                    let part_charge = i32_col(&batch, "part_charge")?;
                    let part_mass_bits = i32_col(&batch, "part_mass_bits")?;
                    let part_seg = i32_col(&batch, "part_seg")?;

                    match tk {
                        "partition" => {
                            let charge = part_charge.value(i);
                            let mass_bits = part_mass_bits.value(i) as u32;
                            let seg = part_seg.value(i);
                            partitions.push(Partition {
                                charge,
                                parent_mass: f32::from_bits(mass_bits),
                                seg_num: seg,
                            });
                        }
                        "precursor_off" => {
                            let charge = part_charge.value(i);
                            let prec_arr = batch
                                .column_by_name("precursor_offsets")
                                .ok_or_else(|| TrainError::Other("missing precursor_offsets".into()))?;
                            let list = prec_arr
                                .as_any()
                                .downcast_ref::<ListArray>()
                                .ok_or_else(|| TrainError::Other("precursor_offsets not ListArray".into()))?;
                            if list.is_null(i) { continue; }
                            let sa = list.value(i);
                            let structs = sa
                                .as_any()
                                .downcast_ref::<StructArray>()
                                .ok_or_else(|| TrainError::Other("precursor_offsets item not StructArray".into()))?;
                            let rc_col = structs.column_by_name("reduced_charge")
                                .ok_or_else(|| TrainError::Other("missing reduced_charge".into()))?;
                            let off_col = structs.column_by_name("offset")
                                .ok_or_else(|| TrainError::Other("missing offset".into()))?;
                            let tol_is_ppm_col = structs.column_by_name("tol_is_ppm")
                                .ok_or_else(|| TrainError::Other("missing tol_is_ppm".into()))?;
                            let tol_val_col = structs.column_by_name("tol_val")
                                .ok_or_else(|| TrainError::Other("missing tol_val".into()))?;
                            let freq_col = structs.column_by_name("frequency")
                                .ok_or_else(|| TrainError::Other("missing frequency".into()))?;
                            let rcs = rc_col.as_any().downcast_ref::<Int32Array>()
                                .ok_or_else(|| TrainError::Other("column reduced_charge has unexpected arrow type".into()))?;
                            let offs = off_col.as_any().downcast_ref::<Float32Array>()
                                .ok_or_else(|| TrainError::Other("column offset has unexpected arrow type".into()))?;
                            let tols = tol_is_ppm_col.as_any().downcast_ref::<BooleanArray>()
                                .ok_or_else(|| TrainError::Other("column tol_is_ppm has unexpected arrow type".into()))?;
                            let tolvs = tol_val_col.as_any().downcast_ref::<Float32Array>()
                                .ok_or_else(|| TrainError::Other("column tol_val has unexpected arrow type".into()))?;
                            let freqs = freq_col.as_any().downcast_ref::<Float32Array>()
                                .ok_or_else(|| TrainError::Other("column frequency has unexpected arrow type".into()))?;
                            let entries = precursor_off_map.entry(charge).or_default();
                            for j in 0..structs.len() {
                                let is_ppm = tols.value(j);
                                let tol_raw = tolvs.value(j) as f64;
                                let tolerance = if is_ppm {
                                    Tolerance::Ppm(tol_raw)
                                } else {
                                    Tolerance::Da(tol_raw)
                                };
                                entries.push(PrecursorOffsetFrequency {
                                    reduced_charge: rcs.value(j),
                                    offset: offs.value(j),
                                    tolerance,
                                    frequency: freqs.value(j),
                                });
                            }
                        }
                        "frag_off" => {
                            let charge = part_charge.value(i);
                            let mass = f32::from_bits(part_mass_bits.value(i) as u32);
                            let seg = part_seg.value(i);
                            let part = Partition { charge, parent_mass: mass, seg_num: seg };

                            let values_arr = list_col(&batch, "values")?;
                            if values_arr.is_null(i) { continue; }
                            let flat_arr = values_arr.value(i);
                            let flat = flat_arr.as_any().downcast_ref::<Float32Array>()
                                .ok_or_else(|| TrainError::Other("frag_off values not Float32Array".into()))?;
                            let len = flat.len();
                            if len % 4 != 0 {
                                return Err(TrainError::Other(
                                    format!("frag_off flat length {} not multiple of 4", len),
                                ));
                            }

                            // Read the parallel frag_off_loss_classes column if present.
                            // Old stores (without this column) yield all loss_class=0.
                            let loss_classes_vec: Vec<i32> = match batch.column_by_name("frag_off_loss_classes") {
                                Some(col) => {
                                    if let Some(list) = col.as_any().downcast_ref::<ListArray>() {
                                        if list.is_null(i) {
                                            Vec::new()
                                        } else {
                                            let item = list.value(i);
                                            if let Some(arr) = item.as_any().downcast_ref::<Int32Array>() {
                                                (0..arr.len()).map(|k| arr.value(k)).collect()
                                            } else {
                                                Vec::new()
                                            }
                                        }
                                    } else {
                                        Vec::new()
                                    }
                                }
                                None => Vec::new(),
                            };

                            let mut frags: Vec<FragmentOffsetFrequency> =
                                Vec::with_capacity(len / 4);
                            let mut j = 0;
                            let mut entry_idx = 0usize;
                            while j + 3 < len {
                                let is_prefix_f = flat.value(j);
                                let ion_charge = flat.value(j + 1) as i32;
                                let offset_bits_f = flat.value(j + 2);
                                // Use to_bits() to recover the exact u32 bit pattern;
                                // the writer stored f32::from_bits(offset_bits) so this
                                // round-trips without precision loss.
                                let offset_bits = offset_bits_f.to_bits();
                                let frequency = flat.value(j + 3);
                                // Read per-entry loss_class from parallel column; default 0.
                                let lc: u8 = loss_classes_vec
                                    .get(entry_idx)
                                    .map(|&v| v.clamp(0, 255) as u8)
                                    .unwrap_or(0);
                                let ion_type = if is_prefix_f > 0.5 {
                                    IonType::Prefix { charge: ion_charge, offset_bits, loss_class: lc }
                                } else if is_prefix_f < -0.5 {
                                    IonType::Noise
                                } else {
                                    IonType::Suffix { charge: ion_charge, offset_bits, loss_class: lc }
                                };
                                frags.push(FragmentOffsetFrequency { ion_type, frequency });
                                j += 4;
                                entry_idx += 1;
                            }
                            frag_off_table.insert(part, frags);
                        }
                        "rank_dist" => {
                            let charge = part_charge.value(i);
                            let mass = f32::from_bits(part_mass_bits.value(i) as u32);
                            let seg = part_seg.value(i);
                            let part = Partition { charge, parent_mass: mass, seg_num: seg };

                            let ik_col = str_col(&batch, "ion_kind")?;
                            let ic_col = i32_col(&batch, "ion_charge")?;
                            let iob_col = i32_col(&batch, "ion_offset_bits")?;

                            let kind = ik_col.value(i);
                            let ic = ic_col.value(i);
                            let iob = iob_col.value(i) as u32;
                            // Read ion_loss_class; default 0 for old stores (column absent).
                            let lc: u8 = match batch.column_by_name("ion_loss_class") {
                                Some(col) => match col.as_any().downcast_ref::<Int32Array>() {
                                    Some(arr) if !arr.is_null(i) => arr.value(i).clamp(0, 255) as u8,
                                    _ => 0,
                                },
                                None => 0,
                            };
                            let ion_type = decode_ion_type(kind, ic, iob, lc)?;

                            let values_arr = list_col(&batch, "values")?;
                            if values_arr.is_null(i) { continue; }
                            let v_arr = values_arr.value(i);
                            let v = v_arr.as_any().downcast_ref::<Float32Array>()
                                .ok_or_else(|| TrainError::Other("rank_dist values not Float32Array".into()))?;
                            let freqs: Vec<f32> = (0..v.len()).map(|k| v.value(k)).collect();

                            rank_dist_table
                                .entry(part)
                                .or_default()
                                .insert(ion_type, freqs);
                        }
                        "ion_err" | "noise_err" | "ion_existence" => {
                            let charge = part_charge.value(i);
                            let mass = f32::from_bits(part_mass_bits.value(i) as u32);
                            let seg = part_seg.value(i);
                            let part = Partition { charge, parent_mass: mass, seg_num: seg };

                            let values_arr = list_col(&batch, "values")?;
                            if values_arr.is_null(i) { continue; }
                            let v_arr = values_arr.value(i);
                            let v = v_arr.as_any().downcast_ref::<Float32Array>()
                                .ok_or_else(|| TrainError::Other("dist values not Float32Array".into()))?;
                            let vals: Vec<f32> = (0..v.len()).map(|k| v.value(k)).collect();
                            match tk {
                                "ion_err" => { ion_err_dist_table.insert(part, vals); }
                                "noise_err" => { noise_err_dist_table.insert(part, vals); }
                                "ion_existence" => { ion_existence_table.insert(part, vals); }
                                _ => {}
                            }
                        }
                        _ => {} // unknown table_kind: skip
                    }
                }
                _ => {} // unknown record_kind: skip
            }
        }
    }
    } // end per-part loop

    let manifest = manifest_opt
        .ok_or_else(|| TrainError::NoModel(model_id.to_string()))?;

    // Reconstruct SpecDataType.
    let activation = model::activation::ActivationMethod::from_name(&manifest.activation)
        .ok_or_else(|| TrainError::Other(format!("unknown activation: {}", manifest.activation)))?;
    let instrument = model::instrument::InstrumentType::from_name(&manifest.instrument)
        .ok_or_else(|| TrainError::Other(format!("unknown instrument: {}", manifest.instrument)))?;
    let enzyme = match &manifest.enzyme {
        Some(e) => Some(
            model::enzyme::Enzyme::from_name(e)
                .ok_or_else(|| TrainError::Other(format!("unknown enzyme: {e}")))?,
        ),
        None => None,
    };
    let protocol = model::protocol::Protocol::from_name(&manifest.protocol)
        .ok_or_else(|| TrainError::Other(format!("unknown protocol: {}", manifest.protocol)))?;

    let data_type = SpecDataType { activation, instrument, enzyme, protocol };

    let mme = if manifest.mme_is_ppm {
        Tolerance::Ppm(manifest.mme_val as f64)
    } else {
        Tolerance::Da(manifest.mme_val as f64)
    };

    // `partitions` was stored in the binary-reader sorted order (charge → seg_num → parent_mass).
    // The binary reader sorts defensively; the parquet store stores in sorted order too.
    // Sort to match the loader invariant.
    partitions.sort();

    let mut param = Param {
        version: manifest.version,
        data_type,
        mme,
        apply_deconvolution: manifest.apply_deconv,
        deconvolution_error_tolerance: manifest.deconv_tol,
        charge_hist: manifest.charge_hist,
        min_charge: manifest.min_charge,
        max_charge: manifest.max_charge,
        num_segments: manifest.num_segments,
        partitions,
        num_precursor_off: manifest.num_precursor_off,
        precursor_off_map,
        frag_off_table,
        max_rank: manifest.max_rank,
        rank_dist_table,
        error_scaling_factor: manifest.error_scaling_factor,
        ion_err_dist_table,
        noise_err_dist_table,
        ion_existence_table,
        partition_ion_types_cache: FxHashMap::default(),
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    };
    param.rebuild_cache();

    // Decode and attach the GBDT signal/noise model if a blob was stored.
    // Enforce the feature-contract at load (finding 3.1): a blob whose declared
    // feature count doesn't match the in-process per-peak extractor would
    // otherwise mis-score via NaN→default_left traversal rather than erroring.
    if let Some(bytes) = manifest.gbdt_bytes.as_ref() {
        let model = scoring_crate::gbdt_eval::GbdtPeakModel::from_bytes(bytes)
            .map_err(|e| TrainError::Other(format!("decode gbdt_model_bytes for '{model_id}': {e}")))?;
        model
            .validate_n_features(scoring_crate::peak_features::N_FEATURES, "peak")
            .map_err(|e| TrainError::Other(format!("gbdt_model_bytes for '{model_id}': {e}")))?;
        param.gbdt_peak_model = Some(model);
    }

    // Decode and attach the fragment-intensity regressor if a blob was stored.
    if let Some(bytes) = manifest.frag_intensity_bytes.as_ref() {
        let model = scoring_crate::gbdt_eval::GbdtPeakModel::from_bytes(bytes)
            .map_err(|e| TrainError::Other(
                format!("decode frag_intensity_model_bytes for '{model_id}': {e}")
            ))?;
        model
            .validate_n_features(scoring_crate::frag_features::N_FRAG_FEATURES, "frag-intensity")
            .map_err(|e| TrainError::Other(format!("frag_intensity_model_bytes for '{model_id}': {e}")))?;
        param.frag_intensity_model = Some(std::sync::Arc::new(model));
    }

    // Decode and attach the rich-ion classifier if a blob was stored.
    if let Some(bytes) = manifest.rich_ion_bytes.as_ref() {
        let model = scoring_crate::gbdt_eval::GbdtPeakModel::from_bytes(bytes)
            .map_err(|e| TrainError::Other(
                format!("decode rich_ion_model_bytes for '{model_id}': {e}")
            ))?;
        model
            .validate_n_features(scoring_crate::ion_features::N_ION_FEATURES, "rich-ion")
            .map_err(|e| TrainError::Other(format!("rich_ion_model_bytes for '{model_id}': {e}")))?;
        param.rich_ion_model = Some(std::sync::Arc::new(model));
    }

    Ok(param)
}

fn parse_manifest_row(
    batch: &arrow::record_batch::RecordBatch,
    i: usize,
) -> Result<ManifestRow, TrainError> {
    let activation = str_col(batch, "activation")?.value(i).to_string();
    let instrument = str_col(batch, "instrument")?.value(i).to_string();
    let enzyme_col = str_col(batch, "enzyme")?;
    let enzyme = if enzyme_col.is_null(i) {
        None
    } else {
        Some(enzyme_col.value(i).to_string())
    };
    let protocol = str_col(batch, "protocol")?.value(i).to_string();
    let version = i32_col(batch, "version")?.value(i);
    let mme_val = f32_col(batch, "mme_val")?.value(i);
    let mme_is_ppm = bool_col(batch, "mme_is_ppm")?.value(i);
    let apply_deconv = bool_col(batch, "apply_deconv")?.value(i);
    let deconv_tol = f32_col(batch, "deconv_tol")?.value(i);
    let num_segments = i32_col(batch, "num_segments")?.value(i);
    let max_rank = i32_col(batch, "max_rank")?.value(i);
    let error_scaling_factor = i32_col(batch, "error_scaling_factor")?.value(i);
    let min_charge = i32_col(batch, "min_charge")?.value(i);
    let max_charge = i32_col(batch, "max_charge")?.value(i);
    let num_precursor_off = i32_col(batch, "num_precursor_off")?.value(i);

    // Parse charge_hist.
    let ch_col = batch
        .column_by_name("charge_hist")
        .ok_or_else(|| TrainError::Other("missing charge_hist".into()))?;
    let ch_list = ch_col.as_any().downcast_ref::<ListArray>()
        .ok_or_else(|| TrainError::Other("charge_hist not ListArray".into()))?;
    let mut charge_hist: Vec<(i32, i32)> = Vec::new();
    if !ch_list.is_null(i) {
        let sa = ch_list.value(i);
        let structs = sa.as_any().downcast_ref::<StructArray>()
            .ok_or_else(|| TrainError::Other("charge_hist item not StructArray".into()))?;
        let charges = structs.column_by_name("charge")
            .ok_or_else(|| TrainError::Other("missing charge_hist.charge".into()))?;
        let counts = structs.column_by_name("count")
            .ok_or_else(|| TrainError::Other("missing charge_hist.count".into()))?;
        let ca = charges.as_any().downcast_ref::<Int32Array>()
            .ok_or_else(|| TrainError::Other("column charge_hist.charge has unexpected arrow type".into()))?;
        let co = counts.as_any().downcast_ref::<Int32Array>()
            .ok_or_else(|| TrainError::Other("column charge_hist.count has unexpected arrow type".into()))?;
        for j in 0..structs.len() {
            charge_hist.push((ca.value(j), co.value(j)));
        }
    }

    // Read the optional GBDT blob. Column absent (old store) → None (backward-compat).
    // Column present but wrong type → hard error (schema corruption).
    let gbdt_bytes: Option<Vec<u8>> = match batch.column_by_name("gbdt_model_bytes") {
        None => None, // backward-compat: old store without the column
        Some(col) => {
            let arr = col.as_any().downcast_ref::<BinaryArray>()
                .ok_or_else(|| TrainError::Other("gbdt_model_bytes column is not BinaryArray".into()))?;
            if arr.is_null(i) { None } else { Some(arr.value(i).to_vec()) }
        }
    };

    // Read the optional fragment-intensity model blob (backward-compatible: column absent → None).
    let frag_intensity_bytes: Option<Vec<u8>> = match batch.column_by_name("frag_intensity_model_bytes") {
        None => None, // backward-compat: old store without the column
        Some(col) => {
            let arr = col.as_any().downcast_ref::<BinaryArray>()
                .ok_or_else(|| TrainError::Other(
                    "frag_intensity_model_bytes column is not BinaryArray".into()
                ))?;
            if arr.is_null(i) { None } else { Some(arr.value(i).to_vec()) }
        }
    };

    // Read the optional rich-ion model blob (backward-compatible: column absent → None).
    let rich_ion_bytes: Option<Vec<u8>> = match batch.column_by_name("rich_ion_model_bytes") {
        None => None, // backward-compat: old store without the column
        Some(col) => {
            let arr = col.as_any().downcast_ref::<BinaryArray>()
                .ok_or_else(|| TrainError::Other(
                    "rich_ion_model_bytes column is not BinaryArray".into()
                ))?;
            if arr.is_null(i) { None } else { Some(arr.value(i).to_vec()) }
        }
    };

    Ok(ManifestRow {
        activation,
        instrument,
        enzyme,
        protocol,
        version,
        mme_val,
        mme_is_ppm,
        apply_deconv,
        deconv_tol,
        num_segments,
        max_rank,
        error_scaling_factor,
        min_charge,
        max_charge,
        num_precursor_off,
        charge_hist,
        gbdt_bytes,
        frag_intensity_bytes,
        rich_ion_bytes,
    })
}

/// Reconstruct an [`IonType`] from the stored kind/charge/offset_bits/loss_class fields.
///
/// `loss_class` defaults to `0` for stores written before the `ion_loss_class` column
/// was added (backward-compatible: old files simply omit the column).
fn decode_ion_type(kind: &str, charge: i32, offset_bits: u32, loss_class: u8) -> Result<IonType, TrainError> {
    match kind {
        "prefix" => Ok(IonType::Prefix { charge, offset_bits, loss_class }),
        "suffix" => Ok(IonType::Suffix { charge, offset_bits, loss_class }),
        "noise" => Ok(IonType::Noise),
        other => Err(TrainError::Other(format!("unknown ion_kind: {other}"))),
    }
}

// ── column accessors ─────────────────────────────────────────────────────────

fn str_col<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a StringArray, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not StringArray")))
}

fn i32_col<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a Int32Array, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<Int32Array>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not Int32Array")))
}

fn f32_col<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a Float32Array, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<Float32Array>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not Float32Array")))
}

fn bool_col<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a BooleanArray, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<BooleanArray>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not BooleanArray")))
}

fn list_col<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a ListArray, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<ListArray>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not ListArray")))
}

fn i64_col<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a Int64Array, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not Int64Array")))
}

// ── source ledger reader ──────────────────────────────────────────────────────

fn read_sources(path: &Path, model_id: &str) -> Result<Vec<SourceLedger>, TrainError> {
    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut ledgers: Vec<SourceLedger> = Vec::new();

    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;

        // Check if this batch has source/stat columns at all.
        if batch.column_by_name("source_id").is_none() {
            continue;
        }

        let record_kind = str_col(&batch, "record_kind")?;
        let mid_col = str_col(&batch, "model_id")?;

        for i in 0..batch.num_rows() {
            if mid_col.value(i) != model_id {
                continue;
            }
            if record_kind.is_null(i) || record_kind.value(i) != "source" {
                continue;
            }

            let source_id = opt_str_col(&batch, "source_id", i);
            let dataset = opt_str_col(&batch, "dataset", i);
            let date = opt_str_col(&batch, "date", i);
            let instrument = opt_str_col(&batch, "src_instrument", i);
            let experiment_class = opt_str_col(&batch, "src_experiment_class", i);

            let n_psms = if let Ok(col) = i64_col(&batch, "n_psms") {
                if col.is_null(i) { 0 } else { col.value(i) }
            } else { 0 };

            let weight = if let Ok(col) = f32_col(&batch, "weight") {
                if col.is_null(i) { 1.0 } else { col.value(i) }
            } else { 1.0 };

            let train_fdr = if let Ok(col) = f32_col(&batch, "train_fdr") {
                if col.is_null(i) { 0.01 } else { col.value(i) }
            } else { 0.01 };

            ledgers.push(SourceLedger {
                source_id,
                dataset,
                n_psms,
                date,
                weight,
                train_fdr,
                instrument,
                experiment_class,
            });
        }
    }

    Ok(ledgers)
}

fn opt_str_col(
    batch: &arrow::record_batch::RecordBatch,
    name: &str,
    i: usize,
) -> String {
    match batch.column_by_name(name) {
        Some(col) => match col.as_any().downcast_ref::<StringArray>() {
            Some(arr) => if arr.is_null(i) { String::new() } else { arr.value(i).to_string() },
            None => String::new(),
        },
        None => String::new(),
    }
}

// ── per-source stats reader ───────────────────────────────────────────────────

fn read_source_stats(path: &Path, model_id: &str, source_id: &str) -> Result<CountStats, TrainError> {
    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut stats = CountStats::new();
    let mut found_any = false;

    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;

        // Skip batches that don't have stat columns.
        if batch.column_by_name("counts").is_none() {
            continue;
        }
        if batch.column_by_name("source_id").is_none() {
            continue;
        }

        let record_kind = str_col(&batch, "record_kind")?;
        let mid_col = str_col(&batch, "model_id")?;
        let sid_col = str_col(&batch, "source_id")?;

        for i in 0..batch.num_rows() {
            if mid_col.value(i) != model_id {
                continue;
            }
            if record_kind.is_null(i) || record_kind.value(i) != "stat" {
                continue;
            }
            if sid_col.is_null(i) || sid_col.value(i) != source_id {
                continue;
            }

            found_any = true;

            let table_kind_col = str_col(&batch, "table_kind")?;
            if table_kind_col.is_null(i) { continue; }
            let tk = table_kind_col.value(i);

            let counts_list = list_col(&batch, "counts")?;

            match tk {
                "rank" => {
                    // Reconstruct (Partition, IonType) from columns.
                    let part = read_partition(&batch, i)?;
                    let ion = read_ion_type(&batch, i)?;

                    if counts_list.is_null(i) { continue; }
                    let arr = counts_list.value(i);
                    let int64_arr = arr.as_any().downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    let counts: Vec<u64> = (0..int64_arr.len()).map(|j| int64_arr.value(j) as u64).collect();
                    stats.rank.insert((part, ion), counts);
                }
                "error" => {
                    let part = read_partition(&batch, i)?;
                    if counts_list.is_null(i) { continue; }
                    let arr = counts_list.value(i);
                    let int64_arr = arr.as_any().downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    let counts: Vec<u64> = (0..int64_arr.len()).map(|j| int64_arr.value(j) as u64).collect();
                    stats.error.insert(part, counts);
                }
                "noise_error" => {
                    let part = read_partition(&batch, i)?;
                    if counts_list.is_null(i) { continue; }
                    let arr = counts_list.value(i);
                    let int64_arr = arr.as_any().downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    let counts: Vec<u64> = (0..int64_arr.len()).map(|j| int64_arr.value(j) as u64).collect();
                    stats.noise_error.insert(part, counts);
                }
                "existence" => {
                    // counts[idx] = count for existence key (partition, idx).
                    let part = read_partition(&batch, i)?;
                    if counts_list.is_null(i) { continue; }
                    let arr = counts_list.value(i);
                    let int64_arr = arr.as_any().downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    for j in 0..int64_arr.len() {
                        let count = int64_arr.value(j) as u64;
                        if count > 0 {
                            stats.existence.insert((part, j as u32), count);
                        }
                    }
                }
                "charge" => {
                    // counts parallel to charge_keys.
                    let charge_keys_col = list_col(&batch, "charge_keys")?;
                    if counts_list.is_null(i) || charge_keys_col.is_null(i) { continue; }

                    let counts_arr = counts_list.value(i);
                    let keys_arr = charge_keys_col.value(i);
                    let int64_arr = counts_arr.as_any().downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("charge counts not Int64Array".into()))?;
                    let int32_arr = keys_arr.as_any().downcast_ref::<Int32Array>()
                        .ok_or_else(|| TrainError::Other("charge_keys not Int32Array".into()))?;
                    for j in 0..int64_arr.len().min(int32_arr.len()) {
                        let charge_key = int32_arr.value(j);
                        let count = int64_arr.value(j) as u64;
                        if count > 0 {
                            stats.charge.insert(charge_key, count);
                        }
                    }
                }
                _ => {} // unknown stat table_kind: skip
            }
        }
    }

    if !found_any {
        return Err(TrainError::NoModel(format!("source_stats({model_id}, {source_id})")));
    }

    Ok(stats)
}

fn read_partition(
    batch: &arrow::record_batch::RecordBatch,
    i: usize,
) -> Result<Partition, TrainError> {
    let charge = i32_col(batch, "part_charge")?.value(i);
    let mass_bits = i32_col(batch, "part_mass_bits")?.value(i) as u32;
    let seg = i32_col(batch, "part_seg")?.value(i);
    Ok(Partition {
        charge,
        parent_mass: f32::from_bits(mass_bits),
        seg_num: seg,
    })
}

fn read_ion_type(
    batch: &arrow::record_batch::RecordBatch,
    i: usize,
) -> Result<IonType, TrainError> {
    let ik_col = str_col(batch, "ion_kind")?;
    let ic_col = i32_col(batch, "ion_charge")?;
    let iob_col = i32_col(batch, "ion_offset_bits")?;
    let kind = ik_col.value(i);
    let ic = ic_col.value(i);
    let iob = iob_col.value(i) as u32;
    // Read ion_loss_class if present; default to 0 for old stores that lack this column.
    let loss_class: u8 = match batch.column_by_name("ion_loss_class") {
        Some(col) => match col.as_any().downcast_ref::<Int32Array>() {
            Some(arr) if !arr.is_null(i) => arr.value(i).clamp(0, 255) as u8,
            _ => 0,
        },
        None => 0,
    };
    decode_ion_type(kind, ic, iob, loss_class)
}

// ── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Resolve the bundled partitioned store directory relative to this file's
    /// manifest dir (same convention the binary uses via `CARGO_MANIFEST_DIR`).
    fn bundled_store_path() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models")
    }

    /// One real partition file of the bundled store (the `Automatic` protocol,
    /// which carries `hcd_qexactive_tryp`), for tests that need a single
    /// Parquet file to copy around.
    fn bundled_partition_file() -> std::path::PathBuf {
        bundled_store_path().join("protocol=Automatic/models.parquet")
    }

    #[test]
    fn selection_entries_returns_bundled_count_with_hcd_qexactive_tryp() {
        let path = bundled_store_path();
        let store = ModelStore::open(&path)
            .expect("failed to open bundled model store");
        let entries = store.selection_entries();
        // The bundle ships 17 own-trained models, one selection entry each. (Two
        // seed-copy phospho slugs with no own training data were dropped; their
        // queries route to the nearest own model via select_nearest.)
        assert_eq!(
            entries.len(),
            17,
            "expected 17 selection entries, got {}",
            entries.len()
        );
        let found = entries
            .iter()
            .any(|e| e.model_id == "hcd_qexactive_tryp");
        assert!(
            found,
            "expected an entry with model_id == \"hcd_qexactive_tryp\""
        );
        let found_tmt = entries
            .iter()
            .any(|e| e.model_id == "cid_lowres_tryp_tmt");
        assert!(
            found_tmt,
            "expected an entry with model_id == \"cid_lowres_tryp_tmt\""
        );
    }

    #[test]
    fn duplicate_model_id_across_parts_is_rejected() {
        // Finding 3.7: two partitions that both declare the same model_id make
        // load order-dependent, so `open` must reject it.
        let src = bundled_partition_file();
        let dir = tempfile::tempdir().expect("tempdir");
        // Two real partitions, each a full copy of one bundled partition → every
        // model_id in it appears twice across parts.
        std::fs::copy(&src, dir.path().join("part_a.parquet")).expect("copy a");
        std::fs::copy(&src, dir.path().join("part_b.parquet")).expect("copy b");
        let msg = match ModelStore::open(dir.path()) {
            Ok(_) => panic!("duplicate model_id across parts must be rejected"),
            Err(e) => format!("{e}"),
        };
        assert!(msg.contains("duplicate model_id"), "unexpected error: {msg}");
    }

    #[test]
    fn temp_suffix_partition_is_ignored() {
        // Finding 3.7: a stale temp partition (dotfile / `.tmp.parquet`) must be
        // skipped so it neither shadows nor duplicates a real partition.
        let src = bundled_partition_file();
        let dir = tempfile::tempdir().expect("tempdir");
        std::fs::copy(&src, dir.path().join("models.parquet")).expect("copy real");
        // A leftover temp copy of the same data — would otherwise trip the
        // duplicate-id guard if it weren't ignored.
        std::fs::copy(&src, dir.path().join(".models.parquet.tmp.parquet")).expect("copy temp");
        std::fs::copy(&src, dir.path().join("models.tmp.parquet")).expect("copy temp2");
        let store = ModelStore::open(dir.path()).expect("temp partitions must be ignored");
        assert!(
            store.model_ids().iter().any(|id| id == "hcd_qexactive_tryp"),
            "real partition must still load"
        );
    }
}
