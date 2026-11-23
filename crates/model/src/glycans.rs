//! Monosaccharide residue masses (a glycan's building blocks as they sit in a chain),
//! derived from their elemental compositions and the element masses in [`crate::mass`].

use crate::mass::{C, H, N, O};

/// Monoisotopic mass of a residue with `c` carbons, `h` hydrogens, `n` nitrogens and `o`
/// oxygens, in Da.
pub const fn residue_mass(c: u32, h: u32, n: u32, o: u32) -> f64 {
    c as f64 * C + h as f64 * H + n as f64 * N + o as f64 * O
}

/// N-acetylhexosamine residue, C8H13NO5, in Da.
pub const HEXNAC: f64 = residue_mass(8, 13, 1, 5);

/// Hexose residue, C6H10O5, in Da.
pub const HEX: f64 = residue_mass(6, 10, 0, 5);

/// Fucose (deoxyhexose) residue, C6H10O4, in Da.
pub const FUC: f64 = residue_mass(6, 10, 0, 4);

/// N-acetylneuraminic acid residue, C11H17NO8, in Da.
pub const NEUAC: f64 = residue_mass(11, 17, 1, 8);

/// N-glycolylneuraminic acid residue, C11H17NO9, in Da.
pub const NEUGC: f64 = residue_mass(11, 17, 1, 9);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residues_match_the_published_monoisotopic_masses() {
        let table = [
            (HEXNAC, 203.079373),
            (HEX, 162.052824),
            (FUC, 146.057909),
            (NEUAC, 291.095417),
            (NEUGC, 307.090331),
        ];
        for (got, want) in table {
            assert!((got - want).abs() < 1e-6, "{got} vs {want}");
        }
    }
}
