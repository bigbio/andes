//! Glycan-Y-complementary index (glycan-Y-first candidate SELECTION).
//!
//! andes today selects glyco candidates on PEPTIDE b/y evidence and treats the
//! glycan Y-ladder as a post-truncation feature — so on stepped-HCD N-glyco
//! spectra (intensity dominated by oxonium + glycan Y-ions, weak peptide b/y)
//! the true candidate is truncated before the strong glycan signal is used. The
//! dedicated glyco engines (a glyco search engine, StrucGP) instead put the GLYCAN on the
//! search axis. This module is that primitive.
//!
//! The a glyco search engine insight: for a Y-ion (= peptide + partial glycan), the
//! **Y-complementary mass** is peptide-independent:
//! ```text
//!   Y_complement = precursor_neutral − Y_ion_neutral
//!                = glycan_mass − partial_glycan_mass
//! ```
//! So we index, per glycan, `glycan_mass − core_partial` for each core-Y rung
//! (Y0..Y5 = the trimannosyl-core ladder, present in every N-glycan). A query
//! computes `precursor_neutral − (peak)` per peak and looks it up in O(#peaks),
//! yielding the count of distinct core-Y rungs matched per glycan — WITHOUT
//! enumerating backbones or scoring peptide b/y. Glycans clearing a core-Y
//! quorum are the selected candidates (their backbone = precursor − glycan).
//!
//! Convention: pass the FULL neutral precursor (`(precursor_mz − PROTON)·z`,
//! water INCLUDED). The core partials use the peptide-neutral Y0 baseline, so
//! `Y_complement` then equals `glycan_mass − partial` exactly.

use std::collections::HashMap;

use crate::glycan_db::GlycanComp;
use crate::glycan_mass::{CORE_Y_STEPS, PROTON};

/// Cumulative core-Y partial-glycan masses: Y0 (0), Y1 (+HexNAc), Y2 (+2HexNAc),
/// Y3..Y5 (+Hex each) — the trimannosyl core, glycan-independent.
fn core_partials() -> [f64; 6] {
    [
        0.0,
        CORE_Y_STEPS[0],
        CORE_Y_STEPS[1],
        CORE_Y_STEPS[2],
        CORE_Y_STEPS[3],
        CORE_Y_STEPS[4],
    ]
}

/// Minimum composition to bear the full trimannosyl core (2×HexNAc + 3×Hex).
fn has_core(g: &GlycanComp) -> bool {
    g.hexnac >= 2 && g.hex >= 3
}

/// Inverted index: quantized `Y_complement` bin → `(glycan_idx, rung, key_mass)`.
pub struct GlycanYIndex {
    tol: f64,
    bin_width: f64,
    bins: HashMap<i64, Vec<(u32, u8, f64)>>,
}

impl GlycanYIndex {
    #[inline]
    fn bin_of(mz: f64, bin_width: f64) -> i64 {
        (mz / bin_width).round() as i64
    }

    /// Build over the glycan list, indexing each glycan's core-Y complement masses.
    /// `tol` is the match tolerance (Da); the bin is a coarse bucket, the stored
    /// key_mass is distance-checked at query time.
    pub fn build(glycans: &[GlycanComp], tol: f64) -> Self {
        let bin_width = tol.max(0.001);
        let partials = core_partials();
        let mut bins: HashMap<i64, Vec<(u32, u8, f64)>> = HashMap::new();
        for (gi, g) in glycans.iter().enumerate() {
            if !has_core(g) {
                continue;
            }
            for (rung, &p) in partials.iter().enumerate() {
                if p > g.mass {
                    continue;
                }
                let key = g.mass - p; // Y_complement for this rung
                let b = Self::bin_of(key, bin_width);
                bins.entry(b).or_default().push((gi as u32, rung as u8, key));
            }
        }
        GlycanYIndex { tol, bin_width, bins }
    }

    /// Return `(glycan_idx, distinct_core_rungs_matched)` for glycans whose core-Y
    /// ladder is supported by at least `min_core` DISTINCT rungs in the spectrum.
    ///
    /// `precursor_neutral` is the FULL neutral precursor (water included). Y-ions
    /// are tried at charges `1..=prec_z` (the core ladder is usually 1+/2+).
    pub fn query(
        &self,
        peaks: &[(f64, f32)],
        precursor_neutral: f64,
        prec_z: u8,
        min_core: u8,
    ) -> Vec<(u32, u8)> {
        // glycan_idx -> bitmask of matched rungs (6 rungs → fits u8).
        let mut matched: HashMap<u32, u8> = HashMap::new();
        let zmax = prec_z.max(1);
        for &(mz, _) in peaks {
            for z in 1..=zmax {
                let y_neutral = (mz - PROTON) * z as f64;
                if y_neutral <= 0.0 || y_neutral >= precursor_neutral {
                    continue;
                }
                let comp = precursor_neutral - y_neutral;
                let b = Self::bin_of(comp, self.bin_width);
                for nb in [b - 1, b, b + 1] {
                    if let Some(v) = self.bins.get(&nb) {
                        for &(gi, rung, key) in v {
                            if (comp - key).abs() <= self.tol {
                                *matched.entry(gi).or_insert(0) |= 1u8 << rung;
                            }
                        }
                    }
                }
            }
        }
        matched
            .into_iter()
            .map(|(gi, mask)| (gi, mask.count_ones() as u8))
            .filter(|&(_, c)| c >= min_core)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_mass::{HEX, HEXNAC};
    use crate::glycan_db::n_glycan_list;

    /// The index must select the true glycan from its core-Y ladder alone — no
    /// peptide b/y, no backbone enumeration — and reject glycans without support.
    #[test]
    fn selects_true_glycan_from_core_y_ladder() {
        let glycans = n_glycan_list();
        let idx = GlycanYIndex::build(&glycans, 0.02);

        // True glycan HexNAc2Hex3 (trimannosyl core). Backbone residue 1500.
        let true_gid = glycans
            .iter()
            .position(|g| g.hexnac == 2 && g.hex == 3 && g.fuc == 0 && g.neuac == 0 && g.neugc == 0)
            .expect("HexNAc2Hex3 in list") as u32;
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let backbone_residue = 1500.0_f64;
        let pep_neutral = backbone_residue + crate::backbone::H2O;
        let precursor_neutral = pep_neutral + glycan_mass; // FULL neutral

        // Core-Y ladder: Y0..Y5 at pep_neutral + cumulative partial + proton.
        let partials = core_partials();
        let peaks: Vec<(f64, f32)> = partials
            .iter()
            .map(|&p| (pep_neutral + p + PROTON, 100.0))
            .chain(std::iter::once((204.087, 200.0))) // oxonium noise
            .collect();

        let hits = idx.query(&peaks, precursor_neutral, 2, 2);
        let found = hits.iter().find(|&&(gid, _)| gid == true_gid);
        assert!(found.is_some(), "true glycan must be selected from its core-Y ladder");
        assert!(found.unwrap().1 >= 5, "expected most core rungs matched, got {}", found.unwrap().1);
    }

    /// A spectrum with no glycan-Y ladder selects nothing at the quorum.
    #[test]
    fn rejects_spectrum_without_core_y() {
        let glycans = n_glycan_list();
        let idx = GlycanYIndex::build(&glycans, 0.02);
        let peaks = vec![(204.087, 100.0), (138.055, 80.0), (900.1, 20.0)];
        // precursor arbitrary; no core-Y ladder present.
        let hits = idx.query(&peaks, 2400.0, 3, 2);
        assert!(hits.is_empty(), "no core-Y ladder → no selection, got {hits:?}");
    }
}
