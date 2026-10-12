//! Match-between-runs outputs.
//!
//! `<base>.mbr.tsv`: long format, one row per (precursor, run) target that was
//! identified in the run or transferred into it:
//!
//! | column | content |
//! |---|---|
//! | `precursor` | peptidoform and charge key |
//! | `run` | run stem |
//! | `evidence_type` | `identified`, `transfer_ms1` (transfer q-value at or below the threshold) or `below_threshold` |
//! | `quantity`, `quantity_se` | fitted MS1 quantity and its standard error |
//! | `transfer_q` | internal transfer quant-quality q-value (not an identification FDR); `NA` for identified rows |
//! | `apex_rt` | fitted apex, seconds |
//! | `library_q` | internal library entry gate (not an identification FDR) |
//!
//! Missing values are `NA`. A `below_threshold` row keeps its quantity so a
//! downstream tool can treat it as censored rather than missing.
//!
//! `<base>.mbr_validation.tsv`: a header and one line with the pseudo-transfer
//! validation (sensitivity lower bound, wrong-peak rate) and the transfer
//! counts at the threshold for the pooled and per-family decoy estimates.
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
    pub precursor: &'a str,
    pub run: &'a str,
    pub evidence_type: EvidenceType,
    pub quantity: Option<f64>,
    pub quantity_se: Option<f64>,
    pub transfer_q: Option<f64>,
    pub apex_rt: Option<f64>,
    pub library_q: f64,
}

pub const MBR_COLUMNS: [&str; 8] = [
    "precursor",
    "run",
    "evidence_type",
    "quantity",
    "quantity_se",
    "transfer_q",
    "apex_rt",
    "library_q",
];

/// The pseudo-transfer validation and transfer counts of one finalize.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MbrValidation {
    /// Transfer q-value threshold the counts are taken at.
    pub q_threshold: f64,
    /// Pseudo-transfer rows (identified pairs measured as transfers).
    pub pseudo_transfers: usize,
    /// Pseudo-transfers with a q-value.
    pub pseudo_scored: usize,
    /// Pseudo-transfers at or below the threshold.
    pub pseudo_accepted: usize,
    /// Pseudo-transfers whose identified fit exists, both fitted.
    pub wrong_peak_pairs: usize,
    /// Of those, apex further than the identified width from the identified apex.
    pub wrong_peak: usize,
    /// Accepted pseudo-transfers among `wrong_peak_pairs`.
    pub accepted_pairs: usize,
    /// Of those, at the wrong peak.
    pub accepted_wrong_peak: usize,
    /// Transfer targets.
    pub transfers: usize,
    /// Transfer targets with a q-value.
    pub transfers_scored: usize,
    /// Transfer targets at or below the threshold, both decoy families.
    pub accepted_pooled: usize,
    /// Mass-shift decoys only.
    pub accepted_mass_only: usize,
    /// RT-shift decoys only.
    pub accepted_rt_only: usize,
}

pub const VALIDATION_COLUMNS: [&str; 16] = [
    "q_threshold",
    "pseudo_transfers",
    "pseudo_scored",
    "pseudo_accepted",
    "sensitivity_lower_bound",
    "wrong_peak_pairs",
    "wrong_peak",
    "wrong_peak_rate",
    "accepted_pairs",
    "accepted_wrong_peak",
    "accepted_wrong_peak_rate",
    "transfers",
    "transfers_scored",
    "accepted_pooled",
    "accepted_mass_only",
    "accepted_rt_only",
];

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
    /// Accepted pseudo-transfers over all pseudo-transfers.
    pub fn sensitivity(&self) -> Option<f64> {
        ratio(self.pseudo_accepted, self.pseudo_transfers)
    }

    pub fn wrong_peak_rate(&self) -> Option<f64> {
        ratio(self.wrong_peak, self.wrong_peak_pairs)
    }

    pub fn accepted_wrong_peak_rate(&self) -> Option<f64> {
        ratio(self.accepted_wrong_peak, self.accepted_pairs)
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
        text.push_str(r.precursor);
        text.push('\t');
        text.push_str(r.run);
        text.push('\t');
        text.push_str(r.evidence_type.as_str());
        push_opt(&mut text, r.quantity);
        push_opt(&mut text, r.quantity_se);
        push_opt(&mut text, r.transfer_q);
        push_opt(&mut text, r.apex_rt);
        push_opt(&mut text, Some(r.library_q));
        text.push('\n');
    }
    write_atomic(path, &text)
}

/// Writes the two-line `<base>.mbr_validation.tsv`.
pub fn write_mbr_validation(path: &Path, v: &MbrValidation) -> Result<(), String> {
    let mut text = VALIDATION_COLUMNS.join("\t");
    text.push('\n');
    let mut line = String::new();
    let _ = write!(line, "{}", v.q_threshold);
    let _ = write!(
        line,
        "\t{}\t{}\t{}",
        v.pseudo_transfers, v.pseudo_scored, v.pseudo_accepted
    );
    push_opt(&mut line, v.sensitivity());
    let _ = write!(line, "\t{}\t{}", v.wrong_peak_pairs, v.wrong_peak);
    push_opt(&mut line, v.wrong_peak_rate());
    let _ = write!(line, "\t{}\t{}", v.accepted_pairs, v.accepted_wrong_peak);
    push_opt(&mut line, v.accepted_wrong_peak_rate());
    let _ = write!(
        line,
        "\t{}\t{}\t{}\t{}\t{}",
        v.transfers,
        v.transfers_scored,
        v.accepted_pooled,
        v.accepted_mass_only,
        v.accepted_rt_only
    );
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
                precursor: "PEPTIDEK/2",
                run: "a",
                evidence_type: EvidenceType::Identified,
                quantity: Some(1.5e6),
                quantity_se: Some(2e4),
                transfer_q: None,
                apex_rt: Some(612.25),
                library_q: 0.001,
            },
            MbrRow {
                precursor: "PEPTIDEK/2",
                run: "b",
                evidence_type: EvidenceType::BelowThreshold,
                quantity: None,
                quantity_se: None,
                transfer_q: Some(0.2),
                apex_rt: None,
                library_q: 0.001,
            },
        ];
        write_mbr_tsv(&p, &rows).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "precursor\trun\tevidence_type\tquantity\tquantity_se\ttransfer_q\tapex_rt\tlibrary_q"
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
    fn validation_is_two_lines_with_rates() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.mbr_validation.tsv");
        let v = MbrValidation {
            q_threshold: 0.01,
            pseudo_transfers: 10,
            pseudo_scored: 9,
            pseudo_accepted: 8,
            wrong_peak_pairs: 8,
            wrong_peak: 1,
            accepted_pairs: 0,
            accepted_wrong_peak: 0,
            transfers: 100,
            transfers_scored: 90,
            accepted_pooled: 40,
            accepted_mass_only: 41,
            accepted_rt_only: 39,
        };
        write_mbr_validation(&p, &v).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        let header: Vec<&str> = lines[0].split('\t').collect();
        let values: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(header, VALIDATION_COLUMNS);
        assert_eq!(values.len(), header.len());
        let get = |name: &str| values[header.iter().position(|h| *h == name).unwrap()];
        assert_eq!(get("sensitivity_lower_bound"), "0.8");
        assert_eq!(get("wrong_peak_rate"), "0.125");
        assert_eq!(get("accepted_wrong_peak_rate"), "NA");
        assert_eq!(get("accepted_rt_only"), "39");
    }
}
