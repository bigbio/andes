//! Percolator command line over andes's native engine (crate percolator-native).
//! Takes the Percolator arguments the engine accepts and writes Percolator's output files.
use std::path::PathBuf;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let limit = percolator_native::MemLimit {
        bytes: 16u64 << 30,
        spill_dir: PathBuf::from("/tmp"),
    };
    if let Err(e) = percolator_native::run(&args, limit) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}
