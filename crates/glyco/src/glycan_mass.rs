pub use model::glycan_residue::{FUC, HEX, HEXNAC, NEUAC, NEUGC};

pub const PROTON: f64 = model::mass::PROTON;

pub const WATER: f64 = model::mass::H2O;

/// Core oxonium ions, each carrying one proton: the C7H7NO2 fragment, HexNAc minus two and
/// one water, HexNAc, and HexNAc+Hex.
pub const CORE_OXONIUM_MZ: [f64; 5] = [
    model::glycan_residue::residue_mass(7, 7, 1, 2) + PROTON,
    HEXNAC + PROTON - 2.0 * WATER,
    HEXNAC + PROTON - WATER,
    HEXNAC + PROTON,
    HEXNAC + HEX + PROTON,
];

pub const CORE_Y_STEPS: [f64; 5] = [
    HEXNAC,
    2.0 * HEXNAC,
    2.0 * HEXNAC + HEX,
    2.0 * HEXNAC + 2.0 * HEX,
    2.0 * HEXNAC + 3.0 * HEX,
];

pub const O_CORE_Y_STEPS: [f64; 1] = [HEXNAC];

pub const MONO_STEPS: [f64; 5] = [HEXNAC, HEX, FUC, HEX + HEXNAC, 2.0 * HEX];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residue_masses_match_standard() {
        let table = [
            (HEXNAC, 203.079373),
            (HEX, 162.052824),
            (FUC, 146.057909),
            (NEUAC, 291.095417),
            (NEUGC, 307.090331),
            (PROTON, 1.007276),
            (WATER, 18.010565),
        ];
        for (got, want) in table {
            assert!(
                (got - want).abs() < 1e-6,
                "mass constant {got} deviates from standard {want} by {:.2e}",
                (got - want).abs()
            );
        }
    }

    #[test]
    fn core_y_steps_are_cumulative_core() {
        assert_eq!(CORE_Y_STEPS[0], HEXNAC);
        assert!((CORE_Y_STEPS[1] - (CORE_Y_STEPS[0] + HEXNAC)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[2] - (CORE_Y_STEPS[1] + HEX)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[3] - (CORE_Y_STEPS[2] + HEX)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[4] - (CORE_Y_STEPS[3] + HEX)).abs() < 1e-4);
    }

    #[test]
    fn o_core_y_steps_is_single_galnac() {
        assert_eq!(O_CORE_Y_STEPS, [HEXNAC]);
    }
}
