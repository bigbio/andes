use crate::pin::RowMeta;
use std::path::Path;

pub fn fmt_g(v: f64, _precision: usize) -> String {
    v.to_string()
}

pub fn fmt_fixed(v: f64, _precision: usize) -> String {
    v.to_string()
}

pub struct Results<'a> {
    pub rows: &'a [usize],
    pub scores: &'a [f64],
    pub q: &'a [f64],
    pub pep: &'a [f64],
}

pub fn write_results(
    target: &Path,
    decoy: &Path,
    _meta: &[RowMeta],
    _file_names: &[String],
    _res: &Results,
) -> Result<(), String> {
    std::fs::write(target, "").map_err(|e| e.to_string())?;
    std::fs::write(decoy, "").map_err(|e| e.to_string())
}

pub fn write_weights(
    path: &Path,
    _names: &[String],
    _w_norm: &[Vec<f64>],
    _w_raw: &[Vec<f64>],
) -> Result<(), String> {
    std::fs::write(path, "").map_err(|e| e.to_string())
}
