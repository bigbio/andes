//! Parquet reader: load [`Param`] models from the combined-schema Parquet store.
//!
//! Usage:
//! ```rust,ignore
//! let store = ModelStore::open(&path)?;
//! let param = store.load_param("my_model")?;
//! ```

use std::path::{Path, PathBuf};

use arrow::array::{
    Array, BooleanArray, Float32Array, Int32Array, ListArray, StringArray, StructArray,
};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rustc_hash::FxHashMap;

use model::tolerance::Tolerance;
use scoring_crate::param_model::{
    FragmentOffsetFrequency, IonType, Param, Partition, PrecursorOffsetFrequency, SpecDataType,
};

use crate::TrainError;

// ── public API ───────────────────────────────────────────────────────────────

/// A handle to a Parquet model store.  Created by [`ModelStore::open`].
pub struct ModelStore {
    path: PathBuf,
    /// IDs of models available in the file (parsed at open time).
    manifest: Vec<String>,
}

impl ModelStore {
    /// Open the store and read the manifest (one pass through the Parquet file).
    pub fn open(path: &Path) -> Result<Self, TrainError> {
        let manifest = read_manifest(path)?;
        Ok(Self { path: path.to_owned(), manifest })
    }

    /// List all model IDs in the store.
    pub fn model_ids(&self) -> &[String] {
        &self.manifest
    }

    /// Load and reconstruct a [`Param`] for the given `model_id`.
    pub fn load_param(&self, model_id: &str) -> Result<Param, TrainError> {
        reconstruct_param(&self.path, model_id)
    }
}

// ── manifest reader ──────────────────────────────────────────────────────────

fn read_manifest(path: &Path) -> Result<Vec<String>, TrainError> {
    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut ids: Vec<String> = Vec::new();
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
                ids.push(mid.value(i).to_string());
            }
        }
    }
    Ok(ids)
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
}

fn reconstruct_param(path: &Path, model_id: &str) -> Result<Param, TrainError> {
    // We do a full scan of the file and collect relevant rows. For large stores
    // this is acceptable as a first implementation; predicate pushdown can be
    // added later.

    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

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
                            let rcs = rc_col.as_any().downcast_ref::<Int32Array>().unwrap();
                            let offs = off_col.as_any().downcast_ref::<Float32Array>().unwrap();
                            let tols = tol_is_ppm_col.as_any().downcast_ref::<BooleanArray>().unwrap();
                            let tolvs = tol_val_col.as_any().downcast_ref::<Float32Array>().unwrap();
                            let freqs = freq_col.as_any().downcast_ref::<Float32Array>().unwrap();
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
                            let mut frags: Vec<FragmentOffsetFrequency> =
                                Vec::with_capacity(len / 4);
                            let mut j = 0;
                            while j + 3 < len {
                                let is_prefix_f = flat.value(j);
                                let ion_charge = flat.value(j + 1) as i32;
                                let offset_bits_f = flat.value(j + 2);
                                // Use to_bits() to recover the exact u32 bit pattern;
                                // the writer stored f32::from_bits(offset_bits) so this
                                // round-trips without precision loss.
                                let offset_bits = offset_bits_f.to_bits();
                                let frequency = flat.value(j + 3);
                                let ion_type = if is_prefix_f > 0.5 {
                                    IonType::Prefix { charge: ion_charge, offset_bits }
                                } else if is_prefix_f < -0.5 {
                                    IonType::Noise
                                } else {
                                    IonType::Suffix { charge: ion_charge, offset_bits }
                                };
                                frags.push(FragmentOffsetFrequency { ion_type, frequency });
                                j += 4;
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
                            let ion_type = decode_ion_type(kind, ic, iob)?;

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
    };
    param.rebuild_cache();
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
        let ca = charges.as_any().downcast_ref::<Int32Array>().unwrap();
        let co = counts.as_any().downcast_ref::<Int32Array>().unwrap();
        for j in 0..structs.len() {
            charge_hist.push((ca.value(j), co.value(j)));
        }
    }

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
    })
}

fn decode_ion_type(kind: &str, charge: i32, offset_bits: u32) -> Result<IonType, TrainError> {
    match kind {
        "prefix" => Ok(IonType::Prefix { charge, offset_bits }),
        "suffix" => Ok(IonType::Suffix { charge, offset_bits }),
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
