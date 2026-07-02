//! Cross-spectrum glycoform transfer (the andes-unique glyco lever).
//!
//! A single peptide is typically observed as MANY glycoforms (in serum N-glyco
//! data, ~6–7 per peptide): same peptide backbone, different glycan. Some
//! glycoforms fragment well (a strong trimannosyl-core Y-ladder → confidently
//! identified); their poorly-fragmenting siblings do not, and per-spectrum
//! candidate generation misses them because there is no core-Y ladder to anchor
//! the backbone.
//!
//! This module transfers a CONFIDENT backbone mass — learned from the well-
//! fragmented glycoforms in a first pass — to a second-pass spectrum whenever
//! `precursor − backbone` is a known glycan. The target spectrum needs NO core-Y
//! ladder of its own: the backbone is borrowed from its siblings. This is the
//! a cross-spectrum glyco engine cross-spectrum idea, and it is something per-spectrum engines
//! (the reference glyco engine) do not do.
//!
//! It is cheap: a small sorted whitelist + a binary-search glycan lookup per
//! precursor, NOT a brute force over the glycan list or a fragment-index scan.

use crate::glycan_db::GlycanComp;

/// A whitelist of confidently-identified backbone (peptide residue) masses from
/// a first pass, sorted and deduplicated.
#[derive(Debug, Clone, Default)]
pub struct GlycoformWhitelist {
    backbones: Vec<f64>,
}

impl GlycoformWhitelist {
    /// Build from confident backbone residue masses, sorting and collapsing
    /// near-duplicates within `dedup_tol` Da (many glycoforms of one peptide
    /// contribute the same backbone; keep it once).
    pub fn new(mut backbones: Vec<f64>, dedup_tol: f64) -> Self {
        backbones.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut out: Vec<f64> = Vec::with_capacity(backbones.len());
        for x in backbones {
            if out.last().map_or(true, |&l| (x - l).abs() > dedup_tol) {
                out.push(x);
            }
        }
        GlycoformWhitelist { backbones: out }
    }

    pub fn len(&self) -> usize {
        self.backbones.len()
    }
    pub fn is_empty(&self) -> bool {
        self.backbones.is_empty()
    }

    /// Propose confident backbones consistent with `precursor_neutral`: for each
    /// whitelisted backbone `bb`, if `precursor_neutral − bb ≥ min_glycan` and is
    /// within `tol` of a known glycan, emit `(bb, glycan)`. `glycan_sorted` is a
    /// `(mass, glycan_index)` view sorted ascending by mass for binary search.
    ///
    /// `tol` should be set from the PRECURSOR mass error (the measured quantity),
    /// consistent with the peptide-first path.
    pub fn transfer(
        &self,
        precursor_neutral: f64,
        glycan_sorted: &[(f64, usize)],
        glycans: &[GlycanComp],
        min_glycan: f64,
        tol: f64,
    ) -> Vec<(f64, GlycanComp)> {
        let mut out: Vec<(f64, GlycanComp)> = Vec::new();
        for &bb in &self.backbones {
            let glycan_mass = precursor_neutral - bb;
            if glycan_mass < min_glycan {
                // whitelist is sorted ascending → larger bb only shrinks the
                // implied glycan, so we can stop.
                break;
            }
            if let Some(g) = nearest_glycan(glycan_sorted, glycans, glycan_mass, tol) {
                out.push((bb, g));
            }
        }
        out
    }
}

/// Nearest known glycan to `target` within `tol`, via a binary-search start on
/// the sorted `(mass, index)` view. (Local copy so the module is self-contained;
/// mirrors the driver's `nearest_glycan_mass`.)
fn nearest_glycan(
    sorted: &[(f64, usize)],
    glycans: &[GlycanComp],
    target: f64,
    tol: f64,
) -> Option<GlycanComp> {
    let lo = target - tol;
    let hi = target + tol;
    let start = sorted.partition_point(|&(m, _)| m < lo);
    let mut best: Option<(f64, usize)> = None;
    for &(m, gi) in &sorted[start..] {
        if m > hi {
            break;
        }
        let d = (m - target).abs();
        if best.map_or(true, |(bd, _)| d < bd) {
            best = Some((d, gi));
        }
    }
    best.map(|(_, gi)| glycans[gi].clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_db::n_glycan_list;
    use crate::glycan_mass::{HEX, HEXNAC};

    fn sorted_view(glycans: &[GlycanComp]) -> Vec<(f64, usize)> {
        let mut v: Vec<(f64, usize)> = glycans.iter().enumerate().map(|(i, g)| (g.mass, i)).collect();
        v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        v
    }

    #[test]
    fn whitelist_dedups_sibling_backbones() {
        let wl = GlycoformWhitelist::new(vec![1500.0, 1500.005, 1500.01, 2000.0], 0.02);
        // The three ~1500 backbones (sibling glycoforms of one peptide) collapse to one.
        assert_eq!(wl.len(), 2);
    }

    #[test]
    fn transfer_recovers_sibling_glycoform_without_a_ladder() {
        let glycans = n_glycan_list();
        let sorted = sorted_view(&glycans);
        // A peptide backbone confidently seen on a well-fragmented glycoform.
        let backbone = 1500.0_f64;
        let wl = GlycoformWhitelist::new(vec![backbone], 0.02);

        // A DIFFERENT glycoform of the same peptide: backbone + HexNAc2Hex3.
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let precursor = backbone + glycan_mass;

        // Transfer proposes the confident backbone for this precursor — no core-Y
        // ladder from the target spectrum is needed.
        let hits = wl.transfer(precursor, &sorted, &glycans, 406.0, 0.05);
        assert!(
            hits.iter().any(|(bb, _)| (bb - backbone).abs() < 0.02),
            "transfer must propose the confident sibling backbone, got {hits:?}"
        );
    }

    #[test]
    fn transfer_rejects_precursor_with_no_valid_glycan() {
        let glycans = n_glycan_list();
        let sorted = sorted_view(&glycans);
        let wl = GlycoformWhitelist::new(vec![1500.0], 0.02);
        // Precursor implies a glycan of ~50 Da (below the N-glycan core) → nothing.
        let hits = wl.transfer(1550.0, &sorted, &glycans, 406.0, 0.05);
        assert!(hits.is_empty(), "sub-core glycan must not transfer, got {hits:?}");
    }
}
