//! Precursor-mass tolerance window check.

use crate::mass::PROTON;
use crate::peptide::Peptide;
use crate::spectrum::Spectrum;
use crate::tolerance::PrecursorTolerance;

#[derive(Debug, Clone, Copy)]
pub struct MassError {
    /// `peptide_mass - spectrum_neutral_mass`. Positive: peptide heavier.
    pub mass_error_da: f64,
    /// `mass_error_da / spectrum_neutral_mass * 1e6`.
    pub mass_error_ppm: f64,
}

/// Returns `Some(error)` if the peptide's neutral mass falls within
/// the tolerance window of the spectrum's neutral mass at the given
/// charge, else `None`.
pub fn matches_precursor(
    spectrum: &Spectrum,
    peptide: &Peptide,
    charge: u8,
    tolerance: &PrecursorTolerance,
) -> Option<MassError> {
    if charge == 0 {
        return None;
    }
    let z = charge as f64;
    let spectrum_neutral = spectrum.precursor_mz * z - z * PROTON;
    let peptide_mass = peptide.mass();
    let mass_error_da = peptide_mass - spectrum_neutral;
    let mass_error_ppm = mass_error_da / spectrum_neutral * 1e6;

    // Convention: negative error → peptide lighter → check LEFT tolerance.
    //             positive error → peptide heavier → check RIGHT tolerance.
    let allowed_da = if mass_error_da < 0.0 {
        tolerance.left.as_da(spectrum_neutral)
    } else {
        tolerance.right.as_da(spectrum_neutral)
    };

    if mass_error_da.abs() <= allowed_da {
        Some(MassError { mass_error_da, mass_error_ppm })
    } else {
        None
    }
}
