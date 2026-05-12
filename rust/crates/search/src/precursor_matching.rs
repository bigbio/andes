//! Precursor-mass tolerance window check.

use model::mass::{ISOTOPE, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use model::tolerance::PrecursorTolerance;

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
) -> Option<MassError> {
    if charge == 0 {
        return None;
    }
    let z = charge as f64;
    let spectrum_neutral_obs = spectrum.precursor_mz * z - z * PROTON;
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
