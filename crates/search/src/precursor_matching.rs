//! Precursor-mass tolerance window check.

use model::mass::{ISOTOPE, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use model::tolerance::PrecursorTolerance;

use crate::precursor_cal::adjusted_observed_neutral_mass;

#[derive(Debug, Clone, Copy)]
pub struct MassError {
    /// `peptide_mass - spectrum_neutral_mass`. Positive: peptide heavier.
    pub mass_error_da: f64,
    /// `mass_error_da / spectrum_neutral_mass * 1e6`.
    pub mass_error_ppm: f64,
    /// Isotope offset that produced this match: 0 = monoisotopic match,
    /// `+N` = spectrum's reported precursor was `N` isotope peaks above
    /// the true monoisotopic. Default range `-1..=2`.
    pub isotope_offset: i8,
}

/// Returns `Some(error)` if the peptide's neutral mass falls within
/// the tolerance window of the spectrum's neutral mass (after
/// `isotope_offset` C13 corrections) at the given charge, else `None`.
///
/// `isotope_offset = 0` is the monoisotopic match. Positive offsets
/// assume the spectrum's reported precursor m/z corresponds to the
/// `+N` isotope envelope (common when the instrument's pick missed
/// the lowest-mass peak); we subtract `N * ISOTOPE` from the spectrum's
/// neutral mass before comparing.
pub fn matches_precursor(
    spectrum: &Spectrum,
    peptide: &Peptide,
    charge: u8,
    isotope_offset: i8,
    tolerance: &PrecursorTolerance,
    shift_ppm: f64,
) -> Option<MassError> {
    if charge == 0 {
        return None;
    }
    let z = charge as f64;
    let spectrum_neutral_obs = adjusted_observed_neutral_mass(
        spectrum.precursor_mz * z - z * PROTON,
        shift_ppm,
    );
    let spectrum_neutral = spectrum_neutral_obs - (isotope_offset as f64) * ISOTOPE;
    let peptide_mass = peptide.mass();
    let mass_error_da = peptide_mass - spectrum_neutral;
    let mass_error_ppm = mass_error_da / spectrum_neutral * 1e6;

    let allowed_da = if mass_error_da < 0.0 {
        tolerance.left.as_da(spectrum_neutral)
    } else {
        tolerance.right.as_da(spectrum_neutral)
    };

    if mass_error_da.abs() <= allowed_da {
        Some(MassError { mass_error_da, mass_error_ppm, isotope_offset })
    } else {
        None
    }
}

/// Chimeric variant of [`matches_precursor`]: returns `Some(error)` when the
/// peptide's neutral mass falls within the spectrum's full **isolation window**
/// (neutral-mass space at `charge`, after `isotope_offset` C13 corrections),
/// expanded by the precursor tolerance. Used by `--chimeric` so co-isolated
/// peptides whose mass is offset from the *selected* precursor are still
/// accepted as candidates.
///
/// `lo_mz` / `hi_mz` are the isolation-window m/z bounds (selected m/z minus the
/// lower offset .. plus the upper offset). The reported `mass_error_*` is
/// measured against the nearest in-window neutral mass (clamped), so a peptide
/// matching anywhere inside the window reports a near-zero error rather than a
/// large error against the selected precursor.
#[allow(clippy::too_many_arguments, reason = "isolation-window bounds + charge + offset + tolerance are all orthogonal inputs")]
pub fn matches_isolation_window(
    peptide: &Peptide,
    charge: u8,
    isotope_offset: i8,
    lo_mz: f64,
    hi_mz: f64,
    tolerance: &PrecursorTolerance,
    shift_ppm: f64,
) -> Option<MassError> {
    if charge == 0 {
        return None;
    }
    let z = charge as f64;
    let iso = (isotope_offset as f64) * ISOTOPE;
    // neutral = mz*z - z*PROTON, monotonic in mz, so neutral_lo <= neutral_hi.
    let neutral_lo = adjusted_observed_neutral_mass(lo_mz * z - z * PROTON, shift_ppm) - iso;
    let neutral_hi = adjusted_observed_neutral_mass(hi_mz * z - z * PROTON, shift_ppm) - iso;
    let peptide_mass = peptide.mass();
    let tol_lo = tolerance.left.as_da(neutral_lo);
    let tol_hi = tolerance.right.as_da(neutral_hi);
    if peptide_mass >= neutral_lo - tol_lo && peptide_mass <= neutral_hi + tol_hi {
        let nearest = peptide_mass.clamp(neutral_lo, neutral_hi);
        let mass_error_da = peptide_mass - nearest;
        let mass_error_ppm = if peptide_mass != 0.0 {
            mass_error_da / peptide_mass * 1e6
        } else {
            0.0
        };
        Some(MassError { mass_error_da, mass_error_ppm, isotope_offset })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;
    use model::tolerance::Tolerance;

    fn peptide() -> Peptide {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        Peptide::from_str("K.PEPTIDEK.R", &aa_set).expect("valid peptide")
    }

    fn tol_20ppm() -> PrecursorTolerance {
        PrecursorTolerance { left: Tolerance::Ppm(20.0), right: Tolerance::Ppm(20.0) }
    }

    #[test]
    fn isolation_window_accepts_off_precursor_peptide_within_window() {
        // A peptide whose neutral mass corresponds to a precursor ~1.0 Da below
        // the selected m/z is REJECTED by the tight precursor check but ACCEPTED
        // by the isolation-window check (window ±1.5 Da).
        let pep = peptide();
        let mass = pep.mass();
        let charge = 2u8;
        let z = charge as f64;
        // Selected precursor placed 1.0 Da (neutral) ABOVE the peptide so the
        // peptide is the off-precursor co-isolated species.
        let selected_neutral = mass + 1.0;
        let selected_mz = (selected_neutral + z * PROTON) / z;
        let tol = tol_20ppm();

        // Tight precursor check (selected m/z) → rejected (1.0 Da >> 20 ppm).
        let mut spec = Spectrum::default();
        spec.precursor_mz = selected_mz;
        assert!(matches_precursor(&spec, &pep, charge, 0, &tol, 0.0).is_none(),
            "tight precursor check should reject the 1 Da off-precursor peptide");

        // Isolation window ±1.5 Da → accepted, near-zero error.
        let lo_mz = selected_mz - 1.5;
        let hi_mz = selected_mz + 1.5;
        let m = matches_isolation_window(&pep, charge, 0, lo_mz, hi_mz, &tol, 0.0)
            .expect("isolation-window check should accept the in-window peptide");
        assert!(m.mass_error_da.abs() <= tol.left.as_da(mass) + 1e-6,
            "in-window peptide should report near-zero mass error, got {}", m.mass_error_da);
    }

    #[test]
    fn isolation_window_rejects_peptide_outside_window() {
        let pep = peptide();
        let mass = pep.mass();
        let charge = 2u8;
        let z = charge as f64;
        // Selected precursor 5 Da above the peptide → outside a ±1.5 Da window.
        let selected_neutral = mass + 5.0;
        let selected_mz = (selected_neutral + z * PROTON) / z;
        let tol = tol_20ppm();
        let lo_mz = selected_mz - 1.5;
        let hi_mz = selected_mz + 1.5;
        assert!(matches_isolation_window(&pep, charge, 0, lo_mz, hi_mz, &tol, 0.0).is_none(),
            "peptide 5 Da outside a ±1.5 Da window must be rejected");
    }
}
