//! Regenerate the bundled Parquet model store from binary `.param` files.
//!
//! NOTE: the 39 source `.param` files were migrated into
//! `resources/ionstat/models.parquet` and **removed from the tree**. To
//! regenerate the store you must first restore the `.param` files into a
//! directory — e.g. from git history:
//!   git checkout <pre-migration-rev> -- resources/ionstat
//! and then point this example at that directory.
//!
//! Usage (from the workspace root):
//!   cargo run -p model-train --example gen_bundled_store -- [PARAM_DIR] [OUT_PARQUET]
//!
//! Defaults: PARAM_DIR = resources/ionstat, OUT_PARQUET = resources/ionstat/models.parquet
use std::path::{Path, PathBuf};

fn main() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let default_dir = Path::new(manifest).join("../../resources/ionstat");
    let args: Vec<String> = std::env::args().collect();
    let param_dir: PathBuf = args.get(1).map(PathBuf::from).unwrap_or_else(|| default_dir.clone());
    let out: PathBuf = args
        .get(2)
        .map(PathBuf::from)
        .unwrap_or_else(|| default_dir.join("models.parquet"));

    if !param_dir.exists() {
        eprintln!(
            "param dir {} does not exist.\n\
             The bundled .param files were removed after migration — restore them first \
             (e.g. `git checkout <pre-migration-rev> -- resources/ionstat`) and pass the \
             directory: cargo run -p model-train --example gen_bundled_store -- <PARAM_DIR> [OUT]",
            param_dir.display()
        );
        std::process::exit(1);
    }

    let ids = model_train::store::migrate_dir(&param_dir, &out)
        .expect("migrate .param files into the parquet store");
    println!("wrote {} models to {}", ids.len(), out.display());
}
