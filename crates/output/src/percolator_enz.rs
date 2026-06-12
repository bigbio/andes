//! Percolator-style enzymatic-boundary helpers (OpenMS PercolatorInfile rules).
//! These compute the `enzN`, `enzC`, and `enzInt` PIN columns that feed
//! Percolator as enzymatic-cleavage consistency features (Kim et al., Nat
//! Commun 5:5277, 2014).
//!
//! ## Conventions
//!
//! - `n` and `c` are the two residues flanking the candidate boundary
//!   (n = the residue immediately N-terminal, c = the residue immediately
//!   C-terminal of the boundary).
//! - Protein-boundary flanking characters always count as enzymatic
//!   (OpenMS `n == '-' || c == '-'` short-circuit). Rust's `Peptide::pre`
//!   uses `_` for the protein N-terminal boundary and `-` for the protein
//!   C-terminal boundary, so both bytes are normalised to the same "boundary"
//!   semantics here.
//! - Unknown / non-builtin enzymes return `true` for any boundary —
//!   matching OpenMS's default "else" branch and Percolator's
//!   unspecific-cleavage semantics.

use model::enzyme::Enzyme;

#[inline]
fn is_protein_boundary(c: u8) -> bool {
    c == b'-' || c == b'_'
}

/// Returns `true` when the boundary between residues `n` and `c` is
/// consistent with the enzyme's cleavage rule (OpenMS PercolatorInfile).
pub(crate) fn is_enzymatic_boundary(n: u8, c: u8, enzyme: Enzyme) -> bool {
    // Protein boundaries are always enzymatic (OpenMS short-circuit),
    // generalised to Rust's `_`/`-` boundary-byte convention.
    if is_protein_boundary(n) || is_protein_boundary(c) {
        return true;
    }
    match enzyme {
        Enzyme::Trypsin => (n == b'K' || n == b'R') && c != b'P',
        Enzyme::Chymotrypsin => (n == b'F' || n == b'W' || n == b'Y' || n == b'L') && c != b'P',
        Enzyme::LysC => n == b'K' && c != b'P',
        Enzyme::LysN => c == b'K',
        Enzyme::GluC => n == b'E' && c != b'P',
        Enzyme::ArgC => n == b'R' && c != b'P',
        Enzyme::AspN => c == b'D',
        // Unknown enzymes: OpenMS default "else" branch returns true so
        // unspecific searches don't penalise every PSM as non-enzymatic.
        Enzyme::AlphaLP | Enzyme::NoCleavage | Enzyme::NonSpecific => true,
    }
}

/// Count internal boundaries consistent with the enzyme cleavage rule
/// (OpenMS PercolatorInfile).
///
/// For an empty / single-residue peptide returns `0` (no internal
/// boundaries to evaluate). For an "unknown" enzyme (universal-true
/// branch above) this returns `len - 1`.
pub(crate) fn count_internal_enzymatic(residues: &[u8], enzyme: Enzyme) -> i32 {
    if residues.len() < 2 {
        return 0;
    }
    let mut count: i32 = 0;
    for i in 1..residues.len() {
        if is_enzymatic_boundary(residues[i - 1], residues[i], enzyme) {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trypsin_cleaves_after_k_r_unless_followed_by_p() {
        // After K with non-P after: enzymatic
        assert!(is_enzymatic_boundary(b'K', b'A', Enzyme::Trypsin));
        assert!(is_enzymatic_boundary(b'R', b'A', Enzyme::Trypsin));
        // After K with P after: not enzymatic
        assert!(!is_enzymatic_boundary(b'K', b'P', Enzyme::Trypsin));
        assert!(!is_enzymatic_boundary(b'R', b'P', Enzyme::Trypsin));
        // Other letters: not enzymatic
        assert!(!is_enzymatic_boundary(b'A', b'B', Enzyme::Trypsin));
    }

    #[test]
    fn protein_boundary_short_circuits_for_all_enzymes() {
        for e in [
            Enzyme::Trypsin, Enzyme::Chymotrypsin, Enzyme::LysC, Enzyme::LysN,
            Enzyme::GluC, Enzyme::ArgC, Enzyme::AspN, Enzyme::AlphaLP,
            Enzyme::NoCleavage, Enzyme::NonSpecific,
        ] {
            // Either side `-` or `_` always cleavable.
            assert!(is_enzymatic_boundary(b'-', b'A', e), "{e:?}");
            assert!(is_enzymatic_boundary(b'A', b'-', e), "{e:?}");
            assert!(is_enzymatic_boundary(b'_', b'A', e), "{e:?}");
            assert!(is_enzymatic_boundary(b'A', b'_', e), "{e:?}");
        }
    }

    #[test]
    fn aspn_cleaves_before_d() {
        assert!(is_enzymatic_boundary(b'A', b'D', Enzyme::AspN));
        assert!(!is_enzymatic_boundary(b'D', b'A', Enzyme::AspN));
    }

    #[test]
    fn lysn_cleaves_before_k() {
        assert!(is_enzymatic_boundary(b'A', b'K', Enzyme::LysN));
        assert!(!is_enzymatic_boundary(b'K', b'A', Enzyme::LysN));
    }

    #[test]
    fn chymotrypsin_cleaves_after_fwy_l_unless_followed_by_p() {
        for n in [b'F', b'W', b'Y', b'L'] {
            assert!(is_enzymatic_boundary(n, b'A', Enzyme::Chymotrypsin));
            assert!(!is_enzymatic_boundary(n, b'P', Enzyme::Chymotrypsin));
        }
        assert!(!is_enzymatic_boundary(b'K', b'A', Enzyme::Chymotrypsin));
    }

    #[test]
    fn unspecific_enzymes_always_cleavable() {
        assert!(is_enzymatic_boundary(b'A', b'A', Enzyme::AlphaLP));
        assert!(is_enzymatic_boundary(b'A', b'A', Enzyme::NonSpecific));
        // NoCleavage follows the unknown-enzyme → true convention (OpenMS default).
        assert!(is_enzymatic_boundary(b'A', b'A', Enzyme::NoCleavage));
    }

    #[test]
    fn count_internal_handles_tryptic_peptide() {
        // PEPTIDKR has internal boundaries: PE EP PT TI ID DK KR
        // (i=1..7), only DK qualifies (after K, then R — wait, position 6 is K-R: after K with R after → enzymatic).
        // Let's verify with a concrete easy case.
        // Peptide: ABKAR → residues [A, B, K, A, R].
        // Internal boundaries at i=1,2,3,4: (A,B), (B,K), (K,A), (A,R)
        //   trypsin: only (K,A) qualifies → count = 1.
        let count = count_internal_enzymatic(b"ABKAR", Enzyme::Trypsin);
        assert_eq!(count, 1);
    }

    #[test]
    fn count_internal_zero_for_short_peptide() {
        assert_eq!(count_internal_enzymatic(b"", Enzyme::Trypsin), 0);
        assert_eq!(count_internal_enzymatic(b"A", Enzyme::Trypsin), 0);
    }

    #[test]
    fn count_internal_handles_p_block() {
        // KPKA: boundaries at i=1,2,3: (K,P), (P,K), (K,A)
        //   trypsin: (K,P) blocked, (P,K) no K/R before, (K,A) yes → count=1.
        assert_eq!(count_internal_enzymatic(b"KPKA", Enzyme::Trypsin), 1);
    }

    #[test]
    fn count_internal_universal_returns_len_minus_one() {
        assert_eq!(count_internal_enzymatic(b"ABCDE", Enzyme::NonSpecific), 4);
    }
}
