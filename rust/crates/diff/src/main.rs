//! msgf-diff: parity comparison tool for MS-GF+ output files.
//!
//! v0 supports a single subcommand, `compare`, which exits 0 when two files
//! are byte-identical and 1 when they differ. Schema-aware comparison and
//! per-field tolerances land in Tasks 3 and 4.

use clap::{Parser, Subcommand};
use std::fs;
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "msgf-diff", version, about = "MS-GF+ output diff tool")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Compare two .pin files at the byte level.
    Compare {
        /// First file (Java reference).
        a: PathBuf,
        /// Second file (Rust output).
        b: PathBuf,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Compare { a, b } => {
            let bytes_a = match fs::read(&a) {
                Ok(b) => b,
                Err(e) => {
                    eprintln!("cannot read {}: {e}", a.display());
                    return ExitCode::from(2);
                }
            };
            let bytes_b = match fs::read(&b) {
                Ok(b) => b,
                Err(e) => {
                    eprintln!("cannot read {}: {e}", b.display());
                    return ExitCode::from(2);
                }
            };
            if bytes_a == bytes_b {
                println!("identical");
                ExitCode::from(0)
            } else {
                println!("different ({} vs {} bytes)", bytes_a.len(), bytes_b.len());
                ExitCode::from(1)
            }
        }
    }
}
