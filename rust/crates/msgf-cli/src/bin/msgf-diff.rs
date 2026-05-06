//! msgf-diff: parity comparison tool for MS-GF+ output files.

use clap::{Parser, Subcommand};
use msgf_cli::{compare_schemas, compare_with_tolerance, PinFile, Tolerance};
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
    /// Compare two .pin files (schema + per-field tolerance).
    Compare {
        a: PathBuf,
        b: PathBuf,
        /// Comma-separated `Field:value` pairs of relative tolerance.
        /// Example: `--tolerance "SpecEValue:1e-3,EValue:1e-3"`.
        #[arg(long, default_value = "")]
        tolerance: String,
    },
}

fn parse_tolerance(spec: &str) -> Tolerance {
    let mut map = Tolerance::new();
    for entry in spec.split(',').filter(|s| !s.is_empty()) {
        if let Some((field, value)) = entry.split_once(':') {
            if let Ok(v) = value.parse::<f64>() {
                map.insert(field.to_string(), v);
            }
        }
    }
    map
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Compare { a, b, tolerance } => {
            let pin_a = match PinFile::read(&a) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!("{e}");
                    return ExitCode::from(2);
                }
            };
            let pin_b = match PinFile::read(&b) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!("{e}");
                    return ExitCode::from(2);
                }
            };
            if let Err(msg) = compare_schemas(&pin_a, &pin_b) {
                eprintln!("{msg}");
                return ExitCode::from(3);
            }
            let tol = parse_tolerance(&tolerance);
            if tol.is_empty() {
                return if pin_a.body == pin_b.body {
                    println!("identical");
                    ExitCode::from(0)
                } else {
                    println!("different (no --tolerance provided; byte-level fail)");
                    ExitCode::from(1)
                };
            }
            match compare_with_tolerance(&pin_a, &pin_b, &tol) {
                Ok(()) => {
                    println!("ok");
                    ExitCode::from(0)
                }
                Err(report) => {
                    eprintln!("{report}");
                    ExitCode::from(1)
                }
            }
        }
    }
}
