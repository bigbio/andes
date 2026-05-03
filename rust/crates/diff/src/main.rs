//! msgf-diff: parity comparison tool for MS-GF+ output files.

use clap::{Parser, Subcommand};
use msgf_diff::{compare_schemas, PinFile};
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
    /// Compare two .pin files (schema + content).
    Compare { a: PathBuf, b: PathBuf },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Compare { a, b } => {
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
            // Schema OK; bodies will be byte-compared in Task 4 with
            // tolerance support. For now, byte-equal bodies pass; otherwise
            // we still flag with exit 1.
            if pin_a.body == pin_b.body {
                println!("identical");
                ExitCode::from(0)
            } else {
                println!("different (rows differ; tolerance compare lands in Task 4)");
                ExitCode::from(1)
            }
        }
    }
}
