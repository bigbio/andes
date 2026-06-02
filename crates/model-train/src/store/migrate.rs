//! Bulk migration of binary `.param` files into a single Parquet model store.
//!
//! Call [`migrate_dir`] with a directory containing `*.param` files and a
//! destination path for the output Parquet file.  Each file becomes one model
//! in the store; the model ID is the filename stem **lowercased** (e.g.
//! `HCD_QExactive_Tryp.param` → `hcd_qexactive_tryp`).

use std::path::{Path, PathBuf};

use scoring_crate::param_model::Param;

use crate::store::write::write_models;
use crate::TrainError;

/// Migrate all `*.param` files in `ionstat` into a single Parquet store at
/// `out`.  Returns a `Vec<(model_id, source_path)>` for every migrated file.
///
/// The files are processed in sorted (by filename) order for determinism.
pub fn migrate_dir(ionstat: &Path, out: &Path) -> Result<Vec<(String, PathBuf)>, TrainError> {
    // Collect all *.param files, sorted by file name for determinism.
    let mut entries: Vec<PathBuf> = std::fs::read_dir(ionstat)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("param"))
        .collect();
    entries.sort();

    if entries.is_empty() {
        return Err(TrainError::Other(format!(
            "no *.param files found in {}",
            ionstat.display()
        )));
    }

    // Load every param and derive a model ID from the filename stem.
    let mut params: Vec<Param> = Vec::with_capacity(entries.len());
    let mut id_and_path: Vec<(String, PathBuf)> = Vec::with_capacity(entries.len());

    for path in &entries {
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| {
                TrainError::Other(format!("cannot derive stem from {}", path.display()))
            })?;
        let model_id = stem.to_lowercase();
        let param = Param::load_from_file(path)
            .map_err(|e| TrainError::Other(format!("failed to load {}: {e}", path.display())))?;
        params.push(param);
        id_and_path.push((model_id, path.clone()));
    }

    // Build the slice expected by write_models: Vec<(String, &Param)>.
    let model_slice: Vec<(String, &Param)> = id_and_path
        .iter()
        .zip(params.iter())
        .map(|((id, _), p)| (id.clone(), p))
        .collect();

    write_models(out, &model_slice)?;

    Ok(id_and_path)
}
