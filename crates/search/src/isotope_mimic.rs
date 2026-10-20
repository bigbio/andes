//! Whether a precursor is better explained by whole isotope steps than by a
//! delta of nearly the same size.
//!
//! Deamidation (+0.984016) sits 19 mDa from one isotope step (+1.003355), and so
//! does an N→D or Q→E substitution. When an unmodified peptide's precursor is
//! picked on its M+1 peak, the deamidated (or substituted) candidate fits the
//! precursor within tolerance. Decoys cannot model this. The rule reads only
//! the precursor mass, so it treats targets and decoys alike.

use model::Peptide;

/// Daltons between successive precursor isotope peaks (13C minus 12C).
pub const ISOTOPE_STEP_DA: f64 = 1.003_355;
/// A delta within this many daltons of one or two isotope steps can be
/// mimicked by a precursor picked on the wrong isotope peak.
pub const ISOTOPE_MIMIC_TOL_DA: f64 = 0.05;

/// `k·step − delta` when `delta ≈ k·step` with `1 ≤ |k| ≤ 2`, else 0.0.
/// Deamidation (+0.984016) gives +0.019339.
pub fn isotope_gap_da(delta_da: f64) -> f64 {
    let k = (delta_da / ISOTOPE_STEP_DA).round();
    let gap = k * ISOTOPE_STEP_DA - delta_da;
    if k != 0.0 && k.abs() <= 2.0 && gap.abs() < ISOTOPE_MIMIC_TOL_DA {
        gap
    } else {
        0.0
    }
}

/// Total gap of the peptide's isotope-like variable mods.
pub fn mod_shift_da(peptide: &Peptide) -> f64 {
    peptide
        .residues
        .iter()
        .filter_map(|aa| aa.mod_.as_ref())
        .filter(|m| !m.fixed)
        .map(|m| isotope_gap_da(m.mass_delta))
        .sum()
}

/// Whether replacing the isotope-like deltas (total gap `shift_da`) by whole
/// isotope steps gives a smaller precursor error. `mass_error_ppm` is
/// `peptide − spectrum`.
pub fn fits_better(peptide: &Peptide, shift_da: f64, mass_error_ppm: f64) -> bool {
    if shift_da == 0.0 {
        return false;
    }
    let err_da = mass_error_ppm * 1e-6 * peptide.mass();
    (err_da + shift_da).abs() < err_da.abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::{AminoAcid, ModLocation, Modification, ResidueSpec};

    fn modded_pep(residue: u8, delta: f64) -> Peptide {
        let m = Modification {
            name: "x".into(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(residue),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
        };
        let aa = AminoAcid::standard(residue).unwrap().with_mod(m);
        Peptide::new(vec![aa, AminoAcid::standard(b'G').unwrap()], b'_', b'-')
    }

    fn plain_pep(seq: &[u8]) -> Peptide {
        Peptide::new(
            seq.iter()
                .map(|&r| AminoAcid::standard(r).unwrap())
                .collect(),
            b'_',
            b'-',
        )
    }

    #[test]
    fn shift_counts_only_isotope_like_mods() {
        let deam = mod_shift_da(&modded_pep(b'N', 0.984016));
        assert!((deam - 0.019339).abs() < 1e-5, "deamidation gap {deam}");
        // Oxidation, acetyl and pyro-Glu are not within 50 mDa of ±1–2 isotope steps.
        for d in [15.994915, 42.010565, -17.026549, -18.010565] {
            assert_eq!(mod_shift_da(&modded_pep(b'Q', d)), 0.0, "delta {d}");
        }
        assert_eq!(mod_shift_da(&plain_pep(b"PEPTIDE")), 0.0);
    }

    #[test]
    fn deamidation_is_a_mimic_at_the_mis_pick_position_only() {
        let pep = modded_pep(b'N', 0.984016);
        let shift = mod_shift_da(&pep);
        let ppm = |err_da: f64| err_da / pep.mass() * 1e6;
        // Peptide 19 mDa lighter than the precursor: an M+1 pick of the unmodified form.
        assert!(fits_better(&pep, shift, ppm(-0.0193)));
        assert!(fits_better(&pep, shift, ppm(-0.012)));
        // A real deamidation sits near zero error.
        assert!(!fits_better(&pep, shift, ppm(0.0)));
        assert!(!fits_better(&pep, shift, ppm(-0.008)));
        assert!(!fits_better(&pep, shift, ppm(0.01)));
        // No isotope-like delta: never a mimic.
        let ox = modded_pep(b'M', 15.994915);
        assert!(!fits_better(&ox, mod_shift_da(&ox), ppm(-0.0193)));
    }

    #[test]
    fn substitution_deltas_are_isotope_like_only_near_a_step() {
        assert!((isotope_gap_da(0.984016) - 0.019339).abs() < 1e-5, "N->D");
        assert!((isotope_gap_da(-0.984016) + 0.019339).abs() < 1e-5, "D->N");
        assert_eq!(isotope_gap_da(15.994915), 0.0, "A->S is not isotope-like");
        assert_eq!(isotope_gap_da(0.0), 0.0);
    }
}
