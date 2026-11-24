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
use scoring_crate::params::{
    FragmentOffsetFrequency, IonType, Param, Partition, PrecursorOffsetFrequency, SpecDataType,
};

use crate::counts::CountStats;
use crate::select::{parse_experiment_class, SelectionEntry};
use crate::store::write::SourceLedger;
use crate::TrainError;

pub struct ModelStore {
    parts: Vec<PathBuf>,
    manifest: Vec<RawManifestEntry>,
}

#[derive(Debug, Clone)]
pub struct RawManifestEntry {
    pub model_id: String,
    pub activation: String,
    pub instrument: String,
    pub enzyme: String,
    pub protocol: String,
}

impl ModelStore {
    pub fn open(path: &Path) -> Result<Self, TrainError> {
        let parts = resolve_parts(path)?;
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

    pub fn model_ids(&self) -> Vec<String> {
        self.manifest.iter().map(|e| e.model_id.clone()).collect()
    }

    pub fn manifest_entries(&self) -> &[RawManifestEntry] {
        &self.manifest
    }

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

    pub fn load_param(&self, model_id: &str) -> Result<Param, TrainError> {
        reconstruct_param(&self.parts, model_id)
    }

    pub fn load_sources(&self, model_id: &str) -> Result<Vec<SourceLedger>, TrainError> {
        let mut ledgers: Vec<SourceLedger> = Vec::new();
        for part in &self.parts {
            ledgers.extend(read_sources(part, model_id)?);
        }
        Ok(ledgers)
    }

    pub fn load_source_stats(
        &self,
        model_id: &str,
        source_id: &str,
    ) -> Result<CountStats, TrainError> {
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

fn collect_parquet_files(dir: &Path, out: &mut Vec<PathBuf>) -> Result<(), TrainError> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let p = entry.path();
        let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if p.is_dir() {
            if name.starts_with('.') || name.starts_with('_') {
                continue;
            }
            collect_parquet_files(&p, out)?;
            continue;
        }
        if p.extension().and_then(|e| e.to_str()) != Some("parquet") {
            continue;
        }
        if name.starts_with('.') || name.ends_with(".tmp.parquet") {
            continue;
        }
        out.push(p);
    }
    Ok(())
}

pub fn protocol_to_experiment_class(protocol: &str) -> std::collections::BTreeSet<String> {
    match protocol {
        "Automatic" | "Standard" => std::collections::BTreeSet::new(),
        "TMT" => parse_experiment_class("tmt"),
        "Phosphorylation" => parse_experiment_class("phospho"),
        "iTRAQ" => parse_experiment_class("itraq"),
        "iTRAQPhospho" => {
            let mut s = std::collections::BTreeSet::new();
            s.insert("itraqphospho".to_string());
            s
        }
        other => parse_experiment_class(other),
    }
}

fn read_manifest(path: &Path) -> Result<Vec<RawManifestEntry>, TrainError> {
    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder
        .build()
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut entries: Vec<RawManifestEntry> = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;
        let record_kind = batch
            .column_by_name("record_kind")
            .ok_or_else(|| TrainError::Other("missing record_kind column".into()))?;
        let model_id_col = batch
            .column_by_name("model_id")
            .ok_or_else(|| TrainError::Other("missing model_id column".into()))?;

        let rk = record_kind
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| TrainError::Other("record_kind not StringArray".into()))?;
        let mid = model_id_col
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| TrainError::Other("model_id not StringArray".into()))?;

        for i in 0..batch.num_rows() {
            if !rk.is_null(i) && rk.value(i) == "manifest" {
                let activation = str_col(&batch, "activation")
                    .map(|c| {
                        if c.is_null(i) {
                            String::new()
                        } else {
                            c.value(i).to_string()
                        }
                    })
                    .unwrap_or_default();
                let instrument = str_col(&batch, "instrument")
                    .map(|c| {
                        if c.is_null(i) {
                            String::new()
                        } else {
                            c.value(i).to_string()
                        }
                    })
                    .unwrap_or_default();
                let enzyme = str_col(&batch, "enzyme")
                    .map(|c| {
                        if c.is_null(i) {
                            String::new()
                        } else {
                            c.value(i).to_string()
                        }
                    })
                    .unwrap_or_default();
                let protocol = str_col(&batch, "protocol")
                    .map(|c| {
                        if c.is_null(i) {
                            String::new()
                        } else {
                            c.value(i).to_string()
                        }
                    })
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
    gbdt_bytes: Option<Vec<u8>>,
    frag_intensity_bytes: Option<Vec<u8>>,
    rich_ion_bytes: Option<Vec<u8>>,
}

fn reconstruct_param(parts: &[PathBuf], model_id: &str) -> Result<Param, TrainError> {
    let mut manifest_opt: Option<ManifestRow> = None;
    let mut partitions: Vec<Partition> = Vec::new();
    let mut precursor_off_map: FxHashMap<i32, Vec<PrecursorOffsetFrequency>> = FxHashMap::default();
    let mut frag_off_table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>> =
        FxHashMap::default();
    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
        FxHashMap::default();
    let mut ion_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
    let mut noise_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
    let mut ion_existence_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();

    for path in parts {
        if manifest_opt.is_some() {
            break;
        }
        let file = std::fs::File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        let reader = builder
            .build()
            .map_err(|e| TrainError::Parquet(e.to_string()))?;

        for batch in reader {
            let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;

            let record_kind = str_col(&batch, "record_kind")?;
            let mid_col = str_col(&batch, "model_id")?;

            for i in 0..batch.num_rows() {
                if mid_col.value(i) != model_id {
                    continue;
                }
                let rk = if record_kind.is_null(i) {
                    ""
                } else {
                    record_kind.value(i)
                };
                match rk {
                    "manifest" => {
                        manifest_opt = Some(parse_manifest_row(&batch, i)?);
                    }
                    "table" => {
                        let table_kind = str_col(&batch, "table_kind")?;
                        if table_kind.is_null(i) {
                            continue;
                        }
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
                                let prec_arr =
                                    batch.column_by_name("precursor_offsets").ok_or_else(|| {
                                        TrainError::Other("missing precursor_offsets".into())
                                    })?;
                                let list = prec_arr
                                    .as_any()
                                    .downcast_ref::<ListArray>()
                                    .ok_or_else(|| {
                                        TrainError::Other("precursor_offsets not ListArray".into())
                                    })?;
                                if list.is_null(i) {
                                    continue;
                                }
                                let sa = list.value(i);
                                let structs =
                                    sa.as_any().downcast_ref::<StructArray>().ok_or_else(|| {
                                        TrainError::Other(
                                            "precursor_offsets item not StructArray".into(),
                                        )
                                    })?;
                                let rc_col =
                                    structs.column_by_name("reduced_charge").ok_or_else(|| {
                                        TrainError::Other("missing reduced_charge".into())
                                    })?;
                                let off_col = structs
                                    .column_by_name("offset")
                                    .ok_or_else(|| TrainError::Other("missing offset".into()))?;
                                let tol_is_ppm_col =
                                    structs.column_by_name("tol_is_ppm").ok_or_else(|| {
                                        TrainError::Other("missing tol_is_ppm".into())
                                    })?;
                                let tol_val_col = structs
                                    .column_by_name("tol_val")
                                    .ok_or_else(|| TrainError::Other("missing tol_val".into()))?;
                                let freq_col = structs
                                    .column_by_name("frequency")
                                    .ok_or_else(|| TrainError::Other("missing frequency".into()))?;
                                let rcs = rc_col.as_any().downcast_ref::<Int32Array>().ok_or_else(
                                    || {
                                        TrainError::Other(
                                            "column reduced_charge has unexpected arrow type"
                                                .into(),
                                        )
                                    },
                                )?;
                                let offs = off_col
                                    .as_any()
                                    .downcast_ref::<Float32Array>()
                                    .ok_or_else(|| {
                                        TrainError::Other(
                                            "column offset has unexpected arrow type".into(),
                                        )
                                    })?;
                                let tols = tol_is_ppm_col
                                    .as_any()
                                    .downcast_ref::<BooleanArray>()
                                    .ok_or_else(|| {
                                    TrainError::Other(
                                        "column tol_is_ppm has unexpected arrow type".into(),
                                    )
                                })?;
                                let tolvs = tol_val_col
                                    .as_any()
                                    .downcast_ref::<Float32Array>()
                                    .ok_or_else(|| {
                                        TrainError::Other(
                                            "column tol_val has unexpected arrow type".into(),
                                        )
                                    })?;
                                let freqs = freq_col
                                    .as_any()
                                    .downcast_ref::<Float32Array>()
                                    .ok_or_else(|| {
                                        TrainError::Other(
                                            "column frequency has unexpected arrow type".into(),
                                        )
                                    })?;
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
                                let part = Partition {
                                    charge,
                                    parent_mass: mass,
                                    seg_num: seg,
                                };

                                let values_arr = list_col(&batch, "values")?;
                                if values_arr.is_null(i) {
                                    continue;
                                }
                                let flat_arr = values_arr.value(i);
                                let flat = flat_arr
                                    .as_any()
                                    .downcast_ref::<Float32Array>()
                                    .ok_or_else(|| {
                                        TrainError::Other("frag_off values not Float32Array".into())
                                    })?;
                                let len = flat.len();
                                if len % 4 != 0 {
                                    return Err(TrainError::Other(format!(
                                        "frag_off flat length {} not multiple of 4",
                                        len
                                    )));
                                }

                                let loss_classes_vec: Vec<i32> = match batch
                                    .column_by_name("frag_off_loss_classes")
                                {
                                    Some(col) => {
                                        if let Some(list) = col.as_any().downcast_ref::<ListArray>()
                                        {
                                            if list.is_null(i) {
                                                Vec::new()
                                            } else {
                                                let item = list.value(i);
                                                if let Some(arr) =
                                                    item.as_any().downcast_ref::<Int32Array>()
                                                {
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
                                    let offset_bits = offset_bits_f.to_bits();
                                    let frequency = flat.value(j + 3);
                                    let lc: u8 = loss_classes_vec
                                        .get(entry_idx)
                                        .map(|&v| v.clamp(0, 255) as u8)
                                        .unwrap_or(0);
                                    let ion_type = if is_prefix_f > 0.5 {
                                        IonType::Prefix {
                                            charge: ion_charge,
                                            offset_bits,
                                            loss_class: lc,
                                        }
                                    } else if is_prefix_f < -0.5 {
                                        IonType::Noise
                                    } else {
                                        IonType::Suffix {
                                            charge: ion_charge,
                                            offset_bits,
                                            loss_class: lc,
                                        }
                                    };
                                    frags.push(FragmentOffsetFrequency {
                                        ion_type,
                                        frequency,
                                    });
                                    j += 4;
                                    entry_idx += 1;
                                }
                                frag_off_table.insert(part, frags);
                            }
                            "rank_dist" => {
                                let charge = part_charge.value(i);
                                let mass = f32::from_bits(part_mass_bits.value(i) as u32);
                                let seg = part_seg.value(i);
                                let part = Partition {
                                    charge,
                                    parent_mass: mass,
                                    seg_num: seg,
                                };

                                let ik_col = str_col(&batch, "ion_kind")?;
                                let ic_col = i32_col(&batch, "ion_charge")?;
                                let iob_col = i32_col(&batch, "ion_offset_bits")?;

                                let kind = ik_col.value(i);
                                let ic = ic_col.value(i);
                                let iob = iob_col.value(i) as u32;
                                let lc: u8 = match batch.column_by_name("ion_loss_class") {
                                    Some(col) => match col.as_any().downcast_ref::<Int32Array>() {
                                        Some(arr) if !arr.is_null(i) => {
                                            arr.value(i).clamp(0, 255) as u8
                                        }
                                        _ => 0,
                                    },
                                    None => 0,
                                };
                                let ion_type = decode_ion_type(kind, ic, iob, lc)?;

                                let values_arr = list_col(&batch, "values")?;
                                if values_arr.is_null(i) {
                                    continue;
                                }
                                let v_arr = values_arr.value(i);
                                let v = v_arr.as_any().downcast_ref::<Float32Array>().ok_or_else(
                                    || {
                                        TrainError::Other(
                                            "rank_dist values not Float32Array".into(),
                                        )
                                    },
                                )?;
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
                                let part = Partition {
                                    charge,
                                    parent_mass: mass,
                                    seg_num: seg,
                                };

                                let values_arr = list_col(&batch, "values")?;
                                if values_arr.is_null(i) {
                                    continue;
                                }
                                let v_arr = values_arr.value(i);
                                let v = v_arr.as_any().downcast_ref::<Float32Array>().ok_or_else(
                                    || TrainError::Other("dist values not Float32Array".into()),
                                )?;
                                let vals: Vec<f32> = (0..v.len()).map(|k| v.value(k)).collect();
                                match tk {
                                    "ion_err" => {
                                        ion_err_dist_table.insert(part, vals);
                                    }
                                    "noise_err" => {
                                        noise_err_dist_table.insert(part, vals);
                                    }
                                    "ion_existence" => {
                                        ion_existence_table.insert(part, vals);
                                    }
                                    _ => {}
                                }
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    let manifest = manifest_opt.ok_or_else(|| TrainError::NoModel(model_id.to_string()))?;

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

    let data_type = SpecDataType {
        activation,
        instrument,
        enzyme,
        protocol,
    };

    let mme = if manifest.mme_is_ppm {
        Tolerance::Ppm(manifest.mme_val as f64)
    } else {
        Tolerance::Da(manifest.mme_val as f64)
    };

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

    if let Some(bytes) = manifest.gbdt_bytes.as_ref() {
        let model = scoring_crate::gbdt::GbdtPeakModel::from_bytes(bytes).map_err(|e| {
            TrainError::Other(format!("decode gbdt_model_bytes for '{model_id}': {e}"))
        })?;
        model
            .validate_features(&scoring_crate::features::peaks::FEATURE_NAMES, "peak")
            .map_err(|e| TrainError::Other(format!("gbdt_model_bytes for '{model_id}': {e}")))?;
        param.gbdt_peak_model = Some(model);
    }

    if let Some(bytes) = manifest.frag_intensity_bytes.as_ref() {
        let model = scoring_crate::gbdt::GbdtPeakModel::from_bytes(bytes).map_err(|e| {
            TrainError::Other(format!(
                "decode frag_intensity_model_bytes for '{model_id}': {e}"
            ))
        })?;
        model
            .validate_features(
                &scoring_crate::features::fragment::FRAG_FEATURE_NAMES,
                "frag-intensity",
            )
            .map_err(|e| {
                TrainError::Other(format!("frag_intensity_model_bytes for '{model_id}': {e}"))
            })?;
        param.frag_intensity_model = Some(std::sync::Arc::new(model));
    }

    if let Some(bytes) = manifest.rich_ion_bytes.as_ref() {
        let model = scoring_crate::gbdt::GbdtPeakModel::from_bytes(bytes).map_err(|e| {
            TrainError::Other(format!("decode rich_ion_model_bytes for '{model_id}': {e}"))
        })?;
        model
            .validate_features(&scoring_crate::features::ion::ION_FEATURE_NAMES, "rich-ion")
            .map_err(|e| {
                TrainError::Other(format!("rich_ion_model_bytes for '{model_id}': {e}"))
            })?;
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

    let ch_col = batch
        .column_by_name("charge_hist")
        .ok_or_else(|| TrainError::Other("missing charge_hist".into()))?;
    let ch_list = ch_col
        .as_any()
        .downcast_ref::<ListArray>()
        .ok_or_else(|| TrainError::Other("charge_hist not ListArray".into()))?;
    let mut charge_hist: Vec<(i32, i32)> = Vec::new();
    if !ch_list.is_null(i) {
        let sa = ch_list.value(i);
        let structs = sa
            .as_any()
            .downcast_ref::<StructArray>()
            .ok_or_else(|| TrainError::Other("charge_hist item not StructArray".into()))?;
        let charges = structs
            .column_by_name("charge")
            .ok_or_else(|| TrainError::Other("missing charge_hist.charge".into()))?;
        let counts = structs
            .column_by_name("count")
            .ok_or_else(|| TrainError::Other("missing charge_hist.count".into()))?;
        let ca = charges
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or_else(|| {
                TrainError::Other("column charge_hist.charge has unexpected arrow type".into())
            })?;
        let co = counts
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or_else(|| {
                TrainError::Other("column charge_hist.count has unexpected arrow type".into())
            })?;
        for j in 0..structs.len() {
            charge_hist.push((ca.value(j), co.value(j)));
        }
    }

    let gbdt_bytes: Option<Vec<u8>> = match batch.column_by_name("gbdt_model_bytes") {
        None => None,
        Some(col) => {
            let arr = col.as_any().downcast_ref::<BinaryArray>().ok_or_else(|| {
                TrainError::Other("gbdt_model_bytes column is not BinaryArray".into())
            })?;
            if arr.is_null(i) {
                None
            } else {
                Some(arr.value(i).to_vec())
            }
        }
    };

    let frag_intensity_bytes: Option<Vec<u8>> =
        match batch.column_by_name("frag_intensity_model_bytes") {
            None => None,
            Some(col) => {
                let arr = col.as_any().downcast_ref::<BinaryArray>().ok_or_else(|| {
                    TrainError::Other("frag_intensity_model_bytes column is not BinaryArray".into())
                })?;
                if arr.is_null(i) {
                    None
                } else {
                    Some(arr.value(i).to_vec())
                }
            }
        };

    let rich_ion_bytes: Option<Vec<u8>> = match batch.column_by_name("rich_ion_model_bytes") {
        None => None,
        Some(col) => {
            let arr = col.as_any().downcast_ref::<BinaryArray>().ok_or_else(|| {
                TrainError::Other("rich_ion_model_bytes column is not BinaryArray".into())
            })?;
            if arr.is_null(i) {
                None
            } else {
                Some(arr.value(i).to_vec())
            }
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

fn decode_ion_type(
    kind: &str,
    charge: i32,
    offset_bits: u32,
    loss_class: u8,
) -> Result<IonType, TrainError> {
    match kind {
        "prefix" => Ok(IonType::Prefix {
            charge,
            offset_bits,
            loss_class,
        }),
        "suffix" => Ok(IonType::Suffix {
            charge,
            offset_bits,
            loss_class,
        }),
        "noise" => Ok(IonType::Noise),
        other => Err(TrainError::Other(format!("unknown ion_kind: {other}"))),
    }
}

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

fn read_sources(path: &Path, model_id: &str) -> Result<Vec<SourceLedger>, TrainError> {
    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder
        .build()
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut ledgers: Vec<SourceLedger> = Vec::new();

    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;

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
                if col.is_null(i) {
                    0
                } else {
                    col.value(i)
                }
            } else {
                0
            };

            let weight = if let Ok(col) = f32_col(&batch, "weight") {
                if col.is_null(i) {
                    1.0
                } else {
                    col.value(i)
                }
            } else {
                1.0
            };

            let train_fdr = if let Ok(col) = f32_col(&batch, "train_fdr") {
                if col.is_null(i) {
                    0.01
                } else {
                    col.value(i)
                }
            } else {
                0.01
            };

            let decay_half_life_days = f32_col(&batch, "decay_half_life_days")
                .ok()
                .filter(|col| !col.is_null(i))
                .map(|col| col.value(i));

            ledgers.push(SourceLedger {
                source_id,
                dataset,
                n_psms,
                date,
                weight,
                train_fdr,
                decay_half_life_days,
                instrument,
                experiment_class,
            });
        }
    }

    Ok(ledgers)
}

fn opt_str_col(batch: &arrow::record_batch::RecordBatch, name: &str, i: usize) -> String {
    match batch.column_by_name(name) {
        Some(col) => match col.as_any().downcast_ref::<StringArray>() {
            Some(arr) => {
                if arr.is_null(i) {
                    String::new()
                } else {
                    arr.value(i).to_string()
                }
            }
            None => String::new(),
        },
        None => String::new(),
    }
}

fn read_source_stats(
    path: &Path,
    model_id: &str,
    source_id: &str,
) -> Result<CountStats, TrainError> {
    let file = std::fs::File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder
        .build()
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut stats = CountStats::new();
    let mut found_any = false;

    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;

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
            if table_kind_col.is_null(i) {
                continue;
            }
            let tk = table_kind_col.value(i);

            let counts_list = list_col(&batch, "counts")?;

            match tk {
                "rank" => {
                    let part = read_partition(&batch, i)?;
                    let ion = read_ion_type(&batch, i)?;

                    if counts_list.is_null(i) {
                        continue;
                    }
                    let arr = counts_list.value(i);
                    let int64_arr = arr
                        .as_any()
                        .downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    let counts: Vec<u64> = (0..int64_arr.len())
                        .map(|j| int64_arr.value(j) as u64)
                        .collect();
                    stats.rank.insert((part, ion), counts);
                }
                "error" => {
                    let part = read_partition(&batch, i)?;
                    if counts_list.is_null(i) {
                        continue;
                    }
                    let arr = counts_list.value(i);
                    let int64_arr = arr
                        .as_any()
                        .downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    let counts: Vec<u64> = (0..int64_arr.len())
                        .map(|j| int64_arr.value(j) as u64)
                        .collect();
                    stats.error.insert(part, counts);
                }
                "noise_error" => {
                    let part = read_partition(&batch, i)?;
                    if counts_list.is_null(i) {
                        continue;
                    }
                    let arr = counts_list.value(i);
                    let int64_arr = arr
                        .as_any()
                        .downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    let counts: Vec<u64> = (0..int64_arr.len())
                        .map(|j| int64_arr.value(j) as u64)
                        .collect();
                    stats.noise_error.insert(part, counts);
                }
                "existence" => {
                    let part = read_partition(&batch, i)?;
                    if counts_list.is_null(i) {
                        continue;
                    }
                    let arr = counts_list.value(i);
                    let int64_arr = arr
                        .as_any()
                        .downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("counts not Int64Array".into()))?;
                    for j in 0..int64_arr.len() {
                        let count = int64_arr.value(j) as u64;
                        if count > 0 {
                            stats.existence.insert((part, j as u32), count);
                        }
                    }
                }
                "charge" => {
                    let charge_keys_col = list_col(&batch, "charge_keys")?;
                    if counts_list.is_null(i) || charge_keys_col.is_null(i) {
                        continue;
                    }

                    let counts_arr = counts_list.value(i);
                    let keys_arr = charge_keys_col.value(i);
                    let int64_arr = counts_arr
                        .as_any()
                        .downcast_ref::<Int64Array>()
                        .ok_or_else(|| TrainError::Other("charge counts not Int64Array".into()))?;
                    let int32_arr = keys_arr
                        .as_any()
                        .downcast_ref::<Int32Array>()
                        .ok_or_else(|| TrainError::Other("charge_keys not Int32Array".into()))?;
                    for j in 0..int64_arr.len().min(int32_arr.len()) {
                        let charge_key = int32_arr.value(j);
                        let count = int64_arr.value(j) as u64;
                        if count > 0 {
                            stats.charge.insert(charge_key, count);
                        }
                    }
                }
                _ => {}
            }
        }
    }

    if !found_any {
        return Err(TrainError::NoModel(format!(
            "source_stats({model_id}, {source_id})"
        )));
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
    let loss_class: u8 = match batch.column_by_name("ion_loss_class") {
        Some(col) => match col.as_any().downcast_ref::<Int32Array>() {
            Some(arr) if !arr.is_null(i) => arr.value(i).clamp(0, 255) as u8,
            _ => 0,
        },
        None => 0,
    };
    decode_ion_type(kind, ic, iob, loss_class)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bundled_store_path() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models")
    }

    fn bundled_partition_file() -> std::path::PathBuf {
        bundled_store_path().join("protocol=Automatic/models.parquet")
    }

    #[test]
    fn selection_entries_returns_bundled_count_with_hcd_qexactive_tryp() {
        let path = bundled_store_path();
        let store = ModelStore::open(&path).expect("failed to open bundled model store");
        let entries = store.selection_entries();
        assert_eq!(
            entries.len(),
            19,
            "expected 19 selection entries, got {}",
            entries.len()
        );
        let found = entries.iter().any(|e| e.model_id == "hcd_qexactive_tryp");
        assert!(
            found,
            "expected an entry with model_id == \"hcd_qexactive_tryp\""
        );
        let found_tmt = entries.iter().any(|e| e.model_id == "cid_lowres_tryp_tmt");
        assert!(
            found_tmt,
            "expected an entry with model_id == \"cid_lowres_tryp_tmt\""
        );
    }

    #[test]
    fn duplicate_model_id_across_parts_is_rejected() {
        let src = bundled_partition_file();
        let dir = tempfile::tempdir().expect("tempdir");
        std::fs::copy(&src, dir.path().join("part_a.parquet")).expect("copy a");
        std::fs::copy(&src, dir.path().join("part_b.parquet")).expect("copy b");
        let msg = match ModelStore::open(dir.path()) {
            Ok(_) => panic!("duplicate model_id across parts must be rejected"),
            Err(e) => format!("{e}"),
        };
        assert!(
            msg.contains("duplicate model_id"),
            "unexpected error: {msg}"
        );
    }

    #[test]
    fn temp_suffix_partition_is_ignored() {
        let src = bundled_partition_file();
        let dir = tempfile::tempdir().expect("tempdir");
        std::fs::copy(&src, dir.path().join("models.parquet")).expect("copy real");
        std::fs::copy(&src, dir.path().join(".models.parquet.tmp.parquet")).expect("copy temp");
        std::fs::copy(&src, dir.path().join("models.tmp.parquet")).expect("copy temp2");
        let store = ModelStore::open(dir.path()).expect("temp partitions must be ignored");
        assert!(
            store
                .model_ids()
                .iter()
                .any(|id| id == "hcd_qexactive_tryp"),
            "real partition must still load"
        );
    }
}
