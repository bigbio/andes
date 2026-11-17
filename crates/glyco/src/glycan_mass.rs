pub const HEXNAC: f64 = 203.079373;

pub const HEX: f64 = 162.052824;

pub const FUC: f64 = 146.057909;

pub const NEUAC: f64 = 291.095417;

pub const NEUGC: f64 = 307.090331;

pub const PROTON: f64 = 1.007276;

pub const WATER: f64 = 18.010565;

pub const CORE_OXONIUM_MZ: [f64; 5] = [138.05496, 168.06552, 186.07608, 204.08665, 366.13947];

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
