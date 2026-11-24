//! Reading MSnet PSM tables (Parquet) as training labels.

use std::fs::File;
use std::path::Path;

use model::{ModLocation, Modification, ResidueSpec, Spectrum};

pub(crate) struct MsnetPsm {
    pub(crate) spectrum: Spectrum,
    pub(crate) peptide: model::peptide::Peptide,
    pub(crate) charge: u8,
}

pub(crate) fn build_msnet_peptide(
    seq: &str,
    res_mod_pos: &[i32],
    res_mod_delta: &[f64],
    nterm_delta: f64,
    cterm_delta: f64,
) -> Result<model::peptide::Peptide, Box<dyn std::error::Error>> {
    if res_mod_pos.len() != res_mod_delta.len() {
        return Err(format!(
            "res_mod_pos ({}) and res_mod_delta ({}) length mismatch",
            res_mod_pos.len(),
            res_mod_delta.len()
        )
        .into());
    }
    let bytes = seq.as_bytes();
    if bytes.is_empty() {
        return Err("empty peptide sequence".into());
    }
    let n = bytes.len();

    let mut residue_delta = vec![0.0f64; n];
    let mut residue_modded = vec![false; n];
    for (&pos1, &delta) in res_mod_pos.iter().zip(res_mod_delta.iter()) {
        if pos1 < 1 || (pos1 as usize) > n {
            return Err(
                format!("res_mod_pos {pos1} out of range for sequence of length {n}").into(),
            );
        }
        let idx = (pos1 - 1) as usize;
        residue_delta[idx] += delta;
        residue_modded[idx] = true;
    }
    if nterm_delta != 0.0 {
        residue_delta[0] += nterm_delta;
        residue_modded[0] = true;
    }
    if cterm_delta != 0.0 {
        residue_delta[n - 1] += cterm_delta;
        residue_modded[n - 1] = true;
    }

    let mut residues = Vec::with_capacity(n);
    for (i, &r) in bytes.iter().enumerate() {
        let aa = model::AminoAcid::standard(r)
            .ok_or_else(|| format!("non-standard residue {:?} at position {}", r as char, i + 1))?;
        if residue_modded[i] {
            let m = Modification {
                name: "msnet".to_string(),
                mass_delta: residue_delta[i],
                residue: ResidueSpec::Specific(r),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            };
            residues.push(aa.with_mod(m));
        } else {
            residues.push(aa);
        }
    }

    Ok(model::peptide::Peptide::new(residues, b'_', b'-'))
}

pub(crate) fn read_msnet_parquet(path: &Path) -> Result<Vec<MsnetPsm>, Box<dyn std::error::Error>> {
    use arrow::array::{Array, Float32Array, Float64Array, Int32Array, ListArray, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = File::open(path).map_err(|e| format!("opening {}: {e}", path.display()))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| format!("parquet reader for {}: {e}", path.display()))?
        .build()
        .map_err(|e| format!("building parquet reader for {}: {e}", path.display()))?;

    let mut out = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| format!("reading batch from {}: {e}", path.display()))?;

        let col = |name: &str| -> Result<&std::sync::Arc<dyn Array>, Box<dyn std::error::Error>> {
            batch
                .column_by_name(name)
                .ok_or_else(|| format!("missing column '{name}' in {}", path.display()).into())
        };

        let seq = col("seq")?
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or("column 'seq' is not a STRING column")?;
        let charge = col("charge")?
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or("column 'charge' is not an INT32 column")?;
        let prec_mz = col("prec_mz")?
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or("column 'prec_mz' is not a DOUBLE column")?;
        let res_mod_pos = col("res_mod_pos")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_pos' is not a LIST column")?;
        let res_mod_delta = col("res_mod_delta")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_delta' is not a LIST column")?;
        let nterm = col("nterm_delta")?
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or("column 'nterm_delta' is not a DOUBLE column")?;
        let cterm = col("cterm_delta")?
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or("column 'cterm_delta' is not a DOUBLE column")?;
        let mz = col("mz")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'mz' is not a LIST column")?;
        let intensity = col("intensity")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'intensity' is not a LIST column")?;

        let list_i32 =
            |list: &ListArray, i: usize| -> Result<Vec<i32>, Box<dyn std::error::Error>> {
                if list.is_null(i) {
                    return Ok(Vec::new());
                }
                let v = list.value(i);
                let a = v
                    .as_any()
                    .downcast_ref::<Int32Array>()
                    .ok_or("list element is not INT32")?;
                Ok((0..a.len()).map(|j| a.value(j)).collect())
            };
        let list_f64 =
            |list: &ListArray, i: usize| -> Result<Vec<f64>, Box<dyn std::error::Error>> {
                if list.is_null(i) {
                    return Ok(Vec::new());
                }
                let v = list.value(i);
                let a = v
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .ok_or("list element is not DOUBLE")?;
                Ok((0..a.len()).map(|j| a.value(j)).collect())
            };
        let list_f32 =
            |list: &ListArray, i: usize| -> Result<Vec<f32>, Box<dyn std::error::Error>> {
                if list.is_null(i) {
                    return Ok(Vec::new());
                }
                let v = list.value(i);
                let a = v
                    .as_any()
                    .downcast_ref::<Float32Array>()
                    .ok_or("list element is not FLOAT")?;
                Ok((0..a.len()).map(|j| a.value(j)).collect())
            };

        for i in 0..batch.num_rows() {
            let seq_s = seq.value(i);
            let ch = charge.value(i);
            if !(1..=255).contains(&ch) {
                return Err(format!("invalid charge {ch} at row {i} of {}", path.display()).into());
            }
            let charge_u8 = ch as u8;

            let positions = list_i32(res_mod_pos, i)?;
            let deltas = list_f64(res_mod_delta, i)?;
            let peptide =
                build_msnet_peptide(seq_s, &positions, &deltas, nterm.value(i), cterm.value(i))?;

            let mzs = list_f32(mz, i)?;
            let ints = list_f32(intensity, i)?;
            if mzs.len() != ints.len() {
                return Err(format!(
                    "mz ({}) and intensity ({}) length mismatch at row {i} of {}",
                    mzs.len(),
                    ints.len(),
                    path.display()
                )
                .into());
            }
            let mut peaks: Vec<(f64, f32)> = mzs
                .iter()
                .zip(ints.iter())
                .map(|(&m, &it)| (m as f64, it))
                .collect();
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            let spectrum = Spectrum {
                title: format!("row{i}"),
                precursor_mz: prec_mz.value(i),
                precursor_intensity: None,
                precursor_charge: Some(charge_u8 as i32),
                rt_seconds: None,
                scan: None,
                peaks,
                activation_method: None,
                isolation_lower_offset: None,
                isolation_upper_offset: None,
                isolation_target_mz: None,
                mobility_window: None,
            };

            out.push(MsnetPsm {
                spectrum,
                peptide,
                charge: charge_u8,
            });
        }
    }
    Ok(out)
}
