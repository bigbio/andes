//! Regenerate the bundled Parquet model store from the binary `.param` files.
//!
//! Run from the workspace root:
//!   cargo run -p model-train --example gen_bundled_store
//!
//! Reads every `resources/ionstat/*.param` and writes a single
//! `resources/ionstat/models.parquet` that replaces them.
use std::path::Path;

fn main() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let ionstat = Path::new(manifest).join("../../resources/ionstat");
    let out = ionstat.join("models.parquet");
    let ids = model_train::store::migrate_dir(&ionstat, &out)
        .expect("migrate bundled .param files into the parquet store");
    println!("wrote {} models to {}", ids.len(), out.display());
}
