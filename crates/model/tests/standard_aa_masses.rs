//! Pin the 20 standard AA monoisotopic residue masses to IUPAC atomic
//! compositions (C/H/N/O/S integer counts per residue). Each mass is
//! computed in-test from those tuples using the chemistry constants in
//! `model::mass`, then compared to the Rust-built
//! `AminoAcid::standard(residue).mass`.

use model::amino_acid::AminoAcid;
use model::mass::{C, H, N, O, S};

fn composition_mass(c: u32, h: u32, n: u32, o: u32, s: u32) -> f64 {
    c as f64 * C + h as f64 * H + n as f64 * N + o as f64 * O + s as f64 * S
}

#[test]
fn all_20_match_iupac_compositions() {
    // (residue, C, H, N, O, S) — standard 20-residue atomic formulae.
    let reference: &[(u8, u32, u32, u32, u32, u32)] = &[
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

    for &(r, c, h, n, o, s) in reference {
        let aa = AminoAcid::standard(r)
            .unwrap_or_else(|| panic!("residue {} missing from standard table", r as char));
        let expected = composition_mass(c, h, n, o, s);
        assert_eq!(
            aa.mass.to_bits(), expected.to_bits(),
            "AA {} drift: rust=0x{:016x}, expected=0x{:016x}",
            r as char, aa.mass.to_bits(), expected.to_bits()
        );
    }
}

#[test]
fn exotic_residues_absent() {
    // Note: `U` (selenocysteine) IS supported — it's a real co-translationally
    // inserted residue and selenoprotein peptides would otherwise be dropped
    // (see amino_acid::tests::selenocysteine_u_is_supported). The ambiguity
    // codes and pyrrolysine (O) remain unsupported.
    for r in [b'O', b'B', b'Z', b'J', b'X'] {
        assert!(
            AminoAcid::standard(r).is_none(),
            "exotic residue {} unexpectedly present", r as char
        );
    }
    // U is now a recognized standard residue.
    assert!(AminoAcid::standard(b'U').is_some(), "selenocysteine (U) must be supported");
}

#[test]
fn nominal_masses_match_integer_formulae() {
    // Integer formula mass: C×12 + H×1 + N×14 + O×16 + S×32.
    let reference: &[(u8, i32)] = &[
        (b'G', 57),  (b'A', 71),  (b'S', 87),  (b'P', 97),
        (b'V', 99),  (b'T', 101), (b'C', 103), (b'L', 113),
        (b'I', 113), (b'N', 114), (b'D', 115), (b'Q', 128),
        (b'K', 128), (b'E', 129), (b'M', 131), (b'H', 137),
        (b'F', 147), (b'R', 156), (b'Y', 163), (b'W', 186),
    ];
    for &(r, expected) in reference {
        let aa = AminoAcid::standard(r).unwrap();
        assert_eq!(aa.nominal_mass(), expected, "nominal mass drift on {}", r as char);
    }
}
