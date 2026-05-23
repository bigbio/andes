//! Pin per-enzyme cleavage rules to Java
//! `edu.ucsd.msjava.msutil.Enzyme` (lines 299-321). Source-of-truth
//! values copied by hand from the Java source.

use model::enzyme::Enzyme;

/// (variant, residues_cleaved_after, residues_cleaved_before)
fn java_rules() -> Vec<(Enzyme, &'static [u8], &'static [u8])> {
    vec![
        (Enzyme::Trypsin,      b"KR",      b""),
        (Enzyme::Chymotrypsin, b"FYWL",    b""),
        (Enzyme::LysC,         b"K",       b""),
        (Enzyme::AspN,         b"",        b"D"),
        (Enzyme::GluC,         b"E",       b""),
        (Enzyme::LysN,         b"",        b"K"),
        (Enzyme::ArgC,         b"R",       b""),
    ]
}

#[test]
fn cleavage_after_matches_java() {
    for (e, after, _) in java_rules() {
        for r in b'A'..=b'Z' {
            let expected = after.contains(&r);
            assert_eq!(
                e.is_cleavable_after(r), expected,
                "{:?}.is_cleavable_after({}) drift", e, r as char
            );
        }
    }
}

#[test]
fn cleavage_before_matches_java() {
    for (e, _, before) in java_rules() {
        for r in b'A'..=b'Z' {
            let expected = before.contains(&r);
            assert_eq!(
                e.is_cleavable_before(r), expected,
                "{:?}.is_cleavable_before({}) drift", e, r as char
            );
        }
    }
}

#[test]
fn no_cleavage_universal_false() {
    for r in b'A'..=b'Z' {
        assert!(!Enzyme::NoCleavage.is_cleavable_after(r));
        assert!(!Enzyme::NoCleavage.is_cleavable_before(r));
    }
}

#[test]
fn nonspecific_universal_true() {
    for r in b'A'..=b'Z' {
        assert!(Enzyme::NonSpecific.is_cleavable_after(r));
        assert!(Enzyme::NonSpecific.is_cleavable_before(r));
    }
}

#[test]
fn alphalp_universal_true() {
    for r in b'A'..=b'Z' {
        assert!(Enzyme::AlphaLP.is_cleavable_after(r));
        assert!(Enzyme::AlphaLP.is_cleavable_before(r));
    }
}
