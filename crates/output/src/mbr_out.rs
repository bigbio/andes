//! Match-between-runs outputs.
//!
//! `<base>.mbr.tsv`: long format, one row per (precursor, run) target that was
//! identified in the run or transferred into it:
//!
//! | column | content |
//! |---|---|
//! | `precursor_key` | peptidoform and charge key |
//! | `run` | run stem |
//! | `evidence_type` | `identified`, `transfer_ms1` (transfer q-value at or below the threshold) or `below_threshold` |
//! | `quantity`, `quantity_se` | fitted MS1 quantity and its standard error |
//! | `transfer_q` | internal transfer quant-quality q-value (not an identification FDR); `NA` for identified rows |
//! | `apex_rt_s` | fitted apex, seconds |
//! | `library_q` | internal library entry gate (not an identification FDR) |
//!
//! Missing values are `NA`. A `below_threshold` row keeps its quantity so a
//! downstream tool can treat it as censored rather than missing.
//!
//! `<base>.mbr_validation.tsv`: a header and one line with the pseudo-transfer
//! validation (acceptance rate, wrong-peak rate) and the accepted
//! transfer counts of the pooled and per-family decoy estimates, always at
//! q <= 0.01 (columns ending `_q01`) and at the run's transfer q-value
//! threshold (columns ending `_at_threshold`, equal to the `_q01` ones when the
//! threshold is 0.01).
//!
//! Both files are written to `<path>.partial` and renamed into place.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// Why a quantity is reported for a (precursor, run).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceType {
    Identified,
    TransferMs1,
    BelowThreshold,
}

impl EvidenceType {
    pub fn as_str(self) -> &'static str {
        match self {
            EvidenceType::Identified => "identified",
            EvidenceType::TransferMs1 => "transfer_ms1",
            EvidenceType::BelowThreshold => "below_threshold",
        }
    }
}

/// One row of `<base>.mbr.tsv`.
#[derive(Debug, Clone, PartialEq)]
pub struct MbrRow<'a> {
    pub precursor_key: &'a str,
    pub run: &'a str,
    pub evidence_type: EvidenceType,
    pub quantity: Option<f64>,
    pub quantity_se: Option<f64>,
    pub transfer_q: Option<f64>,
    pub apex_rt_s: Option<f64>,
    pub library_q: f64,
}

pub const MBR_COLUMNS: [&str; 8] = [
    "precursor_key",
    "run",
    "evidence_type",
    "quantity",
    "quantity_se",
    "transfer_q",
    "apex_rt_s",
    "library_q",
];

/// The default transfer q-value threshold, at which the validation is
/// always reported.
pub const VALIDATION_Q: f64 = 0.01;

/// Accepted counts at one transfer q-value threshold.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Accepted {
    /// Pseudo-transfers at or below the threshold.
    pub pseudo: usize,
    /// Accepted pseudo-transfers among `MbrValidation::wrong_peak_pairs`.
    pub pairs: usize,
    /// Of those, at the wrong peak.
    pub wrong_peak: usize,
    /// Transfer targets at or below the threshold, both decoy families.
    pub pooled: usize,
    /// Mass-shift decoys only.
    pub mass_only: usize,
    /// RT-shift decoys only.
    pub rt_only: usize,
}

/// The pseudo-transfer validation and transfer counts of one finalize.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MbrValidation {
    /// The run's transfer q-value threshold (`at_threshold`).
    pub q_threshold: f64,
    /// Pseudo-transfer rows (identified pairs measured as transfers).
    pub pseudo_transfers: usize,
    /// Pseudo-transfers with a q-value.
    pub pseudo_scored: usize,
    /// Pseudo-transfers whose identified fit exists, both fitted.
    pub wrong_peak_pairs: usize,
    /// Of those, apex further than the identified width from the identified apex.
    pub wrong_peak: usize,
    /// Transfer targets.
    pub transfers: usize,
    /// Transfer targets with a q-value.
    pub transfers_scored: usize,
    /// Accepted at [`VALIDATION_Q`].
    pub at_q01: Accepted,
    /// Accepted at `q_threshold`.
    pub at_threshold: Accepted,
}

/// Columns of the accepted counts, each followed by `_q01` and by
/// `_at_threshold`.
const ACCEPTED_COLUMNS: [&str; 8] = [
    "pseudo_accepted",
    "pseudo_acceptance",
    "accepted_pairs",
    "accepted_wrong_peak",
    "accepted_wrong_peak_rate",
    "accepted_pooled",
    "accepted_mass_only",
    "accepted_rt_only",
];

/// The header of `<base>.mbr_validation.tsv`.
pub fn validation_columns() -> Vec<String> {
    let mut c: Vec<String> = [
        "q_threshold",
        "pseudo_transfers",
        "pseudo_scored",
        "wrong_peak_pairs",
        "wrong_peak",
        "wrong_peak_rate",
        "transfers",
        "transfers_scored",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for suffix in ["_q01", "_at_threshold"] {
        c.extend(ACCEPTED_COLUMNS.iter().map(|n| format!("{n}{suffix}")));
    }
    c
}

fn ratio(num: usize, den: usize) -> Option<f64> {
    (den > 0).then(|| num as f64 / den as f64)
}

fn push_opt(line: &mut String, v: Option<f64>) {
    match v {
        Some(v) if v.is_finite() => {
            let _ = write!(line, "\t{v}");
        }
        _ => line.push_str("\tNA"),
    }
}

impl MbrValidation {
    /// Pseudo-transfers accepted in `a` over all pseudo-transfers. Optimistic
    /// as a transfer sensitivity: the alignment map of a pseudo-transfer's
    /// run was fitted with that pair as an anchor.
    pub fn pseudo_acceptance(&self, a: &Accepted) -> Option<f64> {
        ratio(a.pseudo, self.pseudo_transfers)
    }

    pub fn wrong_peak_rate(&self) -> Option<f64> {
        ratio(self.wrong_peak, self.wrong_peak_pairs)
    }
}

impl Accepted {
    pub fn wrong_peak_rate(&self) -> Option<f64> {
        ratio(self.wrong_peak, self.pairs)
    }
}

fn write_atomic(path: &Path, text: &str) -> Result<(), String> {
    let mut partial = path.as_os_str().to_owned();
    partial.push(".partial");
    let partial = PathBuf::from(partial);
    let written = std::fs::write(&partial, text)
        .map_err(|e| format!("{}: {e}", partial.display()))
        .and_then(|()| {
            std::fs::rename(&partial, path).map_err(|e| format!("{}: {e}", path.display()))
        });
    if written.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}

/// Writes `<base>.mbr.tsv` rows in the given order.
pub fn write_mbr_tsv(path: &Path, rows: &[MbrRow<'_>]) -> Result<(), String> {
    let mut text = MBR_COLUMNS.join("\t");
    text.push('\n');
    for r in rows {
        text.push_str(r.precursor_key);
        text.push('\t');
        text.push_str(r.run);
        text.push('\t');
        text.push_str(r.evidence_type.as_str());
        push_opt(&mut text, r.quantity);
        push_opt(&mut text, r.quantity_se);
        push_opt(&mut text, r.transfer_q);
        push_opt(&mut text, r.apex_rt_s);
        push_opt(&mut text, Some(r.library_q));
        text.push('\n');
    }
    write_atomic(path, &text)
}

/// Writes the two-line `<base>.mbr_validation.tsv`.
pub fn write_mbr_validation(path: &Path, v: &MbrValidation) -> Result<(), String> {
    let mut text = validation_columns().join("\t");
    text.push('\n');
    let mut line = String::new();
    let _ = write!(
        line,
        "{}\t{}\t{}\t{}\t{}",
        v.q_threshold, v.pseudo_transfers, v.pseudo_scored, v.wrong_peak_pairs, v.wrong_peak
    );
    push_opt(&mut line, v.wrong_peak_rate());
    let _ = write!(line, "\t{}\t{}", v.transfers, v.transfers_scored);
    for a in [&v.at_q01, &v.at_threshold] {
        let _ = write!(line, "\t{}", a.pseudo);
        push_opt(&mut line, v.pseudo_acceptance(a));
        let _ = write!(line, "\t{}\t{}", a.pairs, a.wrong_peak);
        push_opt(&mut line, a.wrong_peak_rate());
        let _ = write!(line, "\t{}\t{}\t{}", a.pooled, a.mass_only, a.rt_only);
    }
    text.push_str(&line);
    text.push('\n');
    write_atomic(path, &text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mbr_tsv_writes_na_for_missing_values() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr.tsv");
        let rows = [
            MbrRow {
                precursor_key: "PEPTIDEK/2",
                run: "a",
                evidence_type: EvidenceType::Identified,
                quantity: Some(1.5e6),
                quantity_se: Some(2e4),
                transfer_q: None,
                apex_rt_s: Some(612.25),
                library_q: 0.001,
            },
            MbrRow {
                precursor_key: "PEPTIDEK/2",
                run: "b",
                evidence_type: EvidenceType::BelowThreshold,
                quantity: None,
                quantity_se: None,
                transfer_q: Some(0.2),
                apex_rt_s: None,
                library_q: 0.001,
            },
        ];
        write_mbr_tsv(&p, &rows).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "precursor_key\trun\tevidence_type\tquantity\tquantity_se\ttransfer_q\tapex_rt_s\tlibrary_q"
        );
        assert_eq!(
            lines[1],
            "PEPTIDEK/2\ta\tidentified\t1500000\t20000\tNA\t612.25\t0.001"
        );
        assert_eq!(
            lines[2],
            "PEPTIDEK/2\tb\tbelow_threshold\tNA\tNA\t0.2\tNA\t0.001"
        );
        assert!(!dir.path().join("x.mbr.tsv.partial").exists());
    }

    #[test]
    fn validation_is_two_lines_with_rates_at_both_thresholds() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr_validation.tsv");
        let v = MbrValidation {
            q_threshold: 0.05,
            pseudo_transfers: 10,
            pseudo_scored: 9,
            wrong_peak_pairs: 8,
            wrong_peak: 1,
            transfers: 100,
            transfers_scored: 90,
            at_q01: Accepted {
                pseudo: 8,
                pairs: 0,
                wrong_peak: 0,
                pooled: 40,
                mass_only: 41,
                rt_only: 39,
            },
            at_threshold: Accepted {
                pseudo: 9,
                pairs: 4,
                wrong_peak: 1,
                pooled: 60,
                mass_only: 61,
                rt_only: 59,
            },
        };
        write_mbr_validation(&p, &v).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        let header: Vec<&str> = lines[0].split('\t').collect();
        let values: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(header, validation_columns());
        assert_eq!(values.len(), header.len());
        let get = |name: &str| values[header.iter().position(|h| *h == name).unwrap()];
        assert_eq!(get("q_threshold"), "0.05");
        assert_eq!(get("wrong_peak_rate"), "0.125");
        assert_eq!(get("pseudo_acceptance_q01"), "0.8");
        assert_eq!(get("pseudo_acceptance_at_threshold"), "0.9");
        assert_eq!(get("accepted_wrong_peak_rate_q01"), "NA");
        assert_eq!(get("accepted_wrong_peak_rate_at_threshold"), "0.25");
        assert_eq!(get("accepted_rt_only_q01"), "39");
        assert_eq!(get("accepted_pooled_at_threshold"), "60");
    }
}
