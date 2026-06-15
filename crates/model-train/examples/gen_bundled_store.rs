//! Regenerate the bundled Parquet model store from binary `.param` files.
//!
//! NOTE: the 39 source `.param` files were migrated into
//! `resources/models.parquet` and **removed from the tree**. The originals
//! are archived under `internal-docs/model-archive/`. To regenerate the
//! store, restore the `.param` files into a directory and point this example
//! at it.
//!
//! Usage (from the workspace root):
//!   cargo run -p model-train --example gen_bundled_store -- [PARAM_DIR] [OUT_PARQUET]
//!
//! Defaults: PARAM_DIR = resources/legacy-params, OUT_PARQUET = resources/models.parquet
use std::path::{Path, PathBuf};

fn main() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let default_dir = Path::new(manifest).join("../../resources/legacy-params");
    let args: Vec<String> = std::env::args().collect();
    let param_dir: PathBuf = args.get(1).map(PathBuf::from).unwrap_or_else(|| default_dir.clone());
    let out: PathBuf = args
        .get(2)
        .map(PathBuf::from)
        .unwrap_or_else(|| default_dir.join("models.parquet"));

    if !param_dir.exists() {
        eprintln!(
            "param dir {} does not exist.\n\
             The bundled .param files were migrated into resources/models.parquet and \
             archived under internal-docs/model-archive/. Restore them into a directory \
             first, then pass it: cargo run -p model-train --example gen_bundled_store -- <PARAM_DIR> [OUT]",
            param_dir.display()
        );
        std::process::exit(1);
    }

    let ids = model_train::store::migrate_dir(&param_dir, &out)
        .expect("migrate .param files into the parquet store");
    println!("wrote {} models to {}", ids.len(), out.display());
}
