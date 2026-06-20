//! End-of-run statistics: final tolerances + per-modification PSM tally.
//!
//! andes auto-resolves the scoring model, the precursor/fragment tolerances,
//! and (with `--refine`) the modification set FROM THE DATA — so the parameters
//! the search *ends* with are not necessarily the ones the user passed on the
//! CLI (e.g. precursor calibration tightens the window; a high-res model carries
//! a 20 ppm fragment tolerance even when the input named none). This module
//! reports those final, effective values together with a per-modification PSM
//! count, and writes them both to stderr and to a `statistics.log` next to the
//! PIN output so a run's true parameters are always recoverable.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use model::tolerance::Tolerance;
use search::candidate_gen::Candidate;
use search::psm::TopNQueue;
use search::search_params::SearchParams;

/// Aggregated end-of-run statistics over the rank-1 PSM of every spectrum.
///
/// Counts are taken **pre-FDR** (andes emits a PIN; Percolator applies FDR
/// downstream), over the best (rank-1) candidate of each spectrum.
pub struct RunStatistics {
    pub precursor_tolerance: String,
    pub precursor_calibration: String,
    pub fragment_tolerance: String,
    pub spectra_with_match: usize,
    pub target_psms: usize,
    pub decoy_psms: usize,
    pub unmodified_target_psms: usize,
    /// `(modification name, # rank-1 target PSMs carrying ≥1 of it)`, count-desc.
    pub ptm_counts: Vec<(String, usize)>,
}

fn fmt_tol(t: &Tolerance) -> String {
    match t {
        Tolerance::Ppm(v) => format!("{v} ppm"),
        Tolerance::Da(v) => format!("{v} Da"),
    }
}

impl RunStatistics {
    /// Tally the rank-1 PSM of each non-empty spectrum queue.
    ///
    /// `fragment_tolerance` is the resolved per-model fragment tolerance
    /// (`param.mme`); `params` carries the post-calibration precursor window.
    pub fn compute(
        queues: &[TopNQueue],
        candidates: &[Candidate],
        params: &SearchParams,
        fragment_tolerance: &Tolerance,
    ) -> Self {
        let mut spectra_with_match = 0usize;
        let mut target_psms = 0usize;
        let mut decoy_psms = 0usize;
        let mut unmodified_target_psms = 0usize;
        let mut ptm: BTreeMap<String, usize> = BTreeMap::new();

        for q in queues {
            if q.is_empty() {
                continue;
            }
            spectra_with_match += 1;
            // Rank-1 = best by rank_score (the same ordering the PIN writer uses
            // for the rank-1 row).
            let psms = q.clone().into_rank_sorted_vec();
            let Some(best) = psms.first() else {
                continue;
            };
            let cand = &candidates[best.primary_candidate_idx() as usize];
            if cand.is_decoy {
                decoy_psms += 1;
                continue;
            }
            target_psms += 1;
            // Distinct modification names on this peptide (a fixed Cam-C and a
            // variable Oxidation both count; each is tallied once per PSM).
            let mut names: BTreeSet<&str> = BTreeSet::new();
            for aa in &cand.peptide.residues {
                if let Some(m) = aa.mod_.as_ref() {
                    names.insert(m.name.as_str());
                }
            }
            if names.is_empty() {
                unmodified_target_psms += 1;
            }
            for n in names {
                *ptm.entry(n.to_string()).or_insert(0) += 1;
            }
        }

        let mut ptm_counts: Vec<(String, usize)> = ptm.into_iter().collect();
        // Most-frequent first; ties broken alphabetically for stable output.
        ptm_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        Self {
            precursor_tolerance: format!("{:?}", params.precursor_tolerance),
            precursor_calibration: format!("{:?}", params.precursor_cal_mode),
            fragment_tolerance: fmt_tol(fragment_tolerance),
            spectra_with_match,
            target_psms,
            decoy_psms,
            unmodified_target_psms,
            ptm_counts,
        }
    }

    /// Render the human-readable report (shared by stderr and `statistics.log`).
    pub fn render(&self) -> String {
        // Width for the PTM name column (incl. the "(unmodified)" sentinel).
        let label_w = self
            .ptm_counts
            .iter()
            .map(|(n, _)| n.len())
            .chain(std::iter::once("(unmodified)".len()))
            .max()
            .unwrap_or(12);

        let mut s = String::new();
        s.push_str("──────── andes run summary ────────\n");
        s.push_str(&format!(
            "  Final precursor tolerance : {} (calibration: {})\n",
            self.precursor_tolerance, self.precursor_calibration
        ));
        s.push_str(&format!(
            "  Final fragment tolerance  : {}\n",
            self.fragment_tolerance
        ));
        s.push_str(&format!(
            "  Spectra with a match      : {}\n",
            self.spectra_with_match
        ));
        s.push_str(&format!(
            "  Rank-1 PSMs (pre-FDR)     : {} target, {} decoy\n",
            self.target_psms, self.decoy_psms
        ));
        s.push_str("  PTM report (rank-1 target PSMs carrying each modification):\n");
        if self.ptm_counts.is_empty() {
            s.push_str("    (none — no modified target PSMs)\n");
        } else {
            for (name, count) in &self.ptm_counts {
                s.push_str(&format!("    {name:<label_w$} : {count}\n"));
            }
        }
        s.push_str(&format!(
            "    {:<label_w$} : {}\n",
            "(unmodified)", self.unmodified_target_psms
        ));
        s.push_str("───────────────────────────────────\n");
        s
    }
}

/// Write the rendered report to `path` (the run's `statistics.log`).
pub fn write_statistics_log(stats: &RunStatistics, path: &Path) -> io::Result<()> {
    let mut f = File::create(path)?;
    f.write_all(stats.render().as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn render_lists_ptms_desc_and_unmodified() {
        let stats = RunStatistics {
            precursor_tolerance: "Symmetric(10.0 ppm)".to_string(),
            precursor_calibration: "Auto".to_string(),
            fragment_tolerance: "0.5 Da".to_string(),
            spectra_with_match: 100,
            target_psms: 60,
            decoy_psms: 40,
            unmodified_target_psms: 25,
            ptm_counts: vec![
                ("Carbamidomethyl".to_string(), 30),
                ("Oxidation".to_string(), 10),
            ],
        };
        let out = stats.render();
        assert!(out.contains("Final fragment tolerance  : 0.5 Da"));
        assert!(out.contains("Carbamidomethyl : 30"));
        assert!(out.contains("Oxidation       : 10"));
        assert!(out.contains("(unmodified)    : 25"));
        // Most-frequent PTM listed before the rarer one.
        let cam = out.find("Carbamidomethyl").unwrap();
        let ox = out.find("Oxidation").unwrap();
        assert!(cam < ox, "PTMs should be sorted count-descending");
    }
}
