// Monoisotopic residue masses (water already subtracted): the mass a
// monosaccharide contributes when POLYMERISED (attached), = free sugar − H2O.
// Values pinned to the campaign 6-decimal standard, cross-checked vs
// Unimod/a glyco search engine — see docs/plans/glyco/30-standards/masses.md + roadmap §4, and
// the `residue_masses_match_standard` test below. A glycan attached to Asn is a
// single delta = Σ(residue masses); never add an extra +H2O (that water is only
// for a FREE/released glycan — the −18.0106 double-count is the #1 divergence).

/// HexNAc (N-Acetylhexosamine) residue mass, e.g. GlcNAc/GalNAc.
pub const HEXNAC: f64 = 203.079373;

/// Hex (Hexose) residue mass, e.g. Mannose/Galactose/Glucose.
pub const HEX: f64 = 162.052824;

/// Fuc (Fucose / dHex) residue mass.
pub const FUC: f64 = 146.057909;

/// NeuAc (N-Acetylneuraminic acid / sialic acid) residue mass.
pub const NEUAC: f64 = 291.095417;

/// NeuGc (N-Glycolylneuraminic acid) residue mass. Note NeuGc ≠ NeuAc (+15.995);
/// mixing them shifts sialylated masses by ~16 Da (mouse data, PXD005411).
pub const NEUGC: f64 = 307.090331;

/// Proton mass (Da). Use this, never neutral-H (1.007825) — ~0.5 mDa/charge drift.
pub const PROTON: f64 = 1.007276;

/// Water (monoisotopic). One H2O is added to a PEPTIDE backbone exactly once
/// (Σ residues + H2O); it is NOT added to an attached glycan.
pub const WATER: f64 = 18.010565;

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
    HEXNAC,               // Y1: +HexNAc
    2.0 * HEXNAC,         // Y2: +2HexNAc
    2.0 * HEXNAC + HEX,   // Y3: +2HexNAc+Hex
    2.0 * HEXNAC + 2.0 * HEX, // Y4: +2HexNAc+2Hex
    2.0 * HEXNAC + 3.0 * HEX, // Y5: +2HexNAc+3Hex
];

/// Single-monosaccharide + common combo steps for the Y-ladder walk. Combos are
/// computed from the residue masses so they cannot drift from the standard.
pub const MONO_STEPS: [f64; 5] = [
    HEXNAC,        // 203.079373
    HEX,           // 162.052824
    FUC,           // 146.057909
    HEX + HEXNAC,  // HexHexNAc combo (365.132197)
    2.0 * HEX,     // 2x Hex (324.105648)
];

#[cfg(test)]
mod tests {
    use super::*;

    /// G0 (mass standardization): the residue/constant table must match the
    /// campaign 6-decimal standard (docs/plans/glyco/30-standards/masses.md +
    /// roadmap §4), cross-checked vs Unimod/a glyco search engine. A >~few-µDa drift here
    /// summed over 8 monosaccharides can nudge near-isobaric compositions, so we
    /// pin to 1e-6 Da. Do not retype masses elsewhere — reference these consts.
    #[test]
    fn residue_masses_match_standard() {
        // (const, standard 6-decimal value)
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
        // Y1 must equal HEXNAC (catches literal drift in the first step)
        assert_eq!(CORE_Y_STEPS[0], HEXNAC);
        // Y2 = Y1 + HexNAc; Y3 = Y2 + Hex; Y4 = Y3 + Hex; Y5 = Y4 + Hex
        assert!((CORE_Y_STEPS[1] - (CORE_Y_STEPS[0] + HEXNAC)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[2] - (CORE_Y_STEPS[1] + HEX)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[3] - (CORE_Y_STEPS[2] + HEX)).abs() < 1e-4);
        assert!((CORE_Y_STEPS[4] - (CORE_Y_STEPS[3] + HEX)).abs() < 1e-4);
    }
}
