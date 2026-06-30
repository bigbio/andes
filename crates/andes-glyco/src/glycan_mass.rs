// Monoisotopic residue masses (water already subtracted).
// Values from public/published tables (Glyco-Fragment, UniMod, etc.).

/// HexNAc (N-Acetylhexosamine) residue mass, e.g. GlcNAc/GalNAc.
pub const HEXNAC: f64 = 203.07937;

/// Hex (Hexose) residue mass, e.g. Mannose/Galactose/Glucose.
pub const HEX: f64 = 162.05282;

/// Fuc (Fucose) residue mass.
pub const FUC: f64 = 146.05791;

/// NeuAc (N-Acetylneuraminic acid / sialic acid) residue mass.
pub const NEUAC: f64 = 291.09542;

/// NeuGc (N-Glycolylneuraminic acid) residue mass.
pub const NEUGC: f64 = 307.09033;

/// Proton mass (Da).
pub const PROTON: f64 = 1.0072765;

/// Core oxonium ions: HexNAc-derived fragments + HexHexNAc; singly-charged m/z.
/// [138.05496, 168.06552, 186.07608, 204.08665, 366.13947]
pub const CORE_OXONIUM_MZ: [f64; 5] = [
    138.05496, // HexNAc fragment (loss of C2H4O2 from HexNAc+H)
    168.06552, // HexNAc fragment
    186.07608, // HexNAc fragment (ring-open)
    204.08665, // HexNAc+H (intact oxonium)
    366.13947, // HexHexNAc oxonium
];

/// Y-ion offsets above the peptide backbone for the trimannosyl core ladder (Y1..Y5).
/// Y1 = +HexNAc; Y2 = +2HexNAc; Y3 = +2HexNAc+Hex; Y4 = +2HexNAc+2Hex; Y5 = +2HexNAc+3Hex.
pub const CORE_Y_STEPS: [f64; 5] = [
    203.07937, // Y1: +HexNAc
    406.15874, // Y2: +2HexNAc
    568.21156, // Y3: +2HexNAc+Hex
    730.26438, // Y4: +2HexNAc+2Hex
    892.31720, // Y5: +2HexNAc+3Hex
];

/// Single-monosaccharide + common combo steps for the Y-ladder walk.
pub const MONO_STEPS: [f64; 5] = [
    HEXNAC,    // 203.07937
    HEX,       // 162.05282
    FUC,       // 146.05791
    365.13219, // HexHexNAc combo
    324.10565, // HexNeuAc combo
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn core_y_steps_are_cumulative_core() {
        // Y2 = Y1 + HexNAc; Y3 = Y2 + Hex; Y4 = Y3 + Hex; Y5 = Y4 + Hex
        assert!((CORE_Y_STEPS[1] - (CORE_Y_STEPS[0] + HEXNAC)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[2] - (CORE_Y_STEPS[1] + HEX)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[3] - (CORE_Y_STEPS[2] + HEX)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[4] - (CORE_Y_STEPS[3] + HEX)).abs() < 1e-4);
    }
}
