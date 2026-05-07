//! Pin the 20 standard AA monoisotopic residue masses to Java
//! `edu.ucsd.msjava.msutil.AminoAcid.STANDARD_AA[]`. Source-of-truth:
//! the (C, H, N, O, S) integer composition tuples copied from
//! `AminoAcid.java:163-181`. Each mass is computed in-test from those
//! tuples using the chemistry constants in `model::mass`, then
//! compared to the Rust-built `AminoAcid::standard(residue).mass`.

use model::amino_acid::AminoAcid;
use model::mass::{C, H, N, O, S};

fn java_composition_mass(c: u32, h: u32, n: u32, o: u32, s: u32) -> f64 {
    c as f64 * C + h as f64 * H + n as f64 * N + o as f64 * O + s as f64 * S
}

#[test]
fn all_20_match_java() {
    // (residue, C, H, N, O, S) — exact integer counts from
    // edu.ucsd.msjava.msutil.AminoAcid.STANDARD_AA[].
    let java: &[(u8, u32, u32, u32, u32, u32)] = &[
        (b'G',  2,  3, 1, 1, 0), (b'A',  3,  5, 1, 1, 0),
        (b'S',  3,  5, 1, 2, 0), (b'P',  5,  7, 1, 1, 0),
        (b'V',  5,  9, 1, 1, 0), (b'T',  4,  7, 1, 2, 0),
        (b'C',  3,  5, 1, 1, 1), (b'L',  6, 11, 1, 1, 0),
        (b'I',  6, 11, 1, 1, 0), (b'N',  4,  6, 2, 2, 0),
        (b'D',  4,  5, 1, 3, 0), (b'Q',  5,  8, 2, 2, 0),
        (b'K',  6, 12, 2, 1, 0), (b'E',  5,  7, 1, 3, 0),
        (b'M',  5,  9, 1, 1, 1), (b'H',  6,  7, 3, 1, 0),
        (b'F',  9,  9, 1, 1, 0), (b'R',  6, 12, 4, 1, 0),
        (b'Y',  9,  9, 1, 2, 0), (b'W', 11, 10, 2, 1, 0),
    ];

    for &(r, c, h, n, o, s) in java {
        let aa = AminoAcid::standard(r)
            .unwrap_or_else(|| panic!("residue {} missing from standard table", r as char));
        let expected = java_composition_mass(c, h, n, o, s);
        assert_eq!(
            aa.mass.to_bits(), expected.to_bits(),
            "AA {} drift: rust=0x{:016x}, java=0x{:016x}",
            r as char, aa.mass.to_bits(), expected.to_bits()
        );
    }
}

#[test]
fn exotic_residues_absent() {
    // U, O, B, Z, J, X are NOT in the standard table — Phase 1 spec
    // explicitly excludes them.
    for r in [b'U', b'O', b'B', b'Z', b'J', b'X'] {
        assert!(
            AminoAcid::standard(r).is_none(),
            "exotic residue {} unexpectedly present", r as char
        );
    }
}

#[test]
fn nominal_masses_match_java() {
    // Java AminoAcid stores nominalMass via Composition.getNominalMass()
    // = C*12 + H*1 + N*14 + O*16 + S*32. For Phase 1 we compute it via
    // `nominal_from(mass)` (the mass-based path); these happen to agree
    // for all 20 standard AAs (verified by inspection — see Composition
    // integer formulae). This test pins that agreement.
    let java: &[(u8, i32)] = &[
        (b'G', 57),  (b'A', 71),  (b'S', 87),  (b'P', 97),
        (b'V', 99),  (b'T', 101), (b'C', 103), (b'L', 113),
        (b'I', 113), (b'N', 114), (b'D', 115), (b'Q', 128),
        (b'K', 128), (b'E', 129), (b'M', 131), (b'H', 137),
        (b'F', 147), (b'R', 156), (b'Y', 163), (b'W', 186),
    ];
    for &(r, expected) in java {
        let aa = AminoAcid::standard(r).unwrap();
        assert_eq!(aa.nominal_mass(), expected, "nominal mass drift on {}", r as char);
    }
}
