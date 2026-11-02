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
//!   (the usual `n == '-' || c == '-'` short-circuit). Rust's `Peptide::pre`
//!   uses `_` for the protein N-terminal boundary and `-` for the protein
//!   C-terminal boundary, so both bytes are normalised to the same "boundary"
//!   semantics here.
//! - A protein N-terminal peptide starts at a terminus even when the
//!   initiator Met was removed: [`enzyme_features`] reads its `M` flank as a
//!   protein boundary (OpenMS PercolatorInfile, OpenMS/OpenMS#10740).
//! - The rules keep the proline clause (no cleavage before P for trypsin,
//!   Lys-C, Arg-C, Glu-C and chymotrypsin), although the digest has none.
//!   This is deliberate: with the digest's rule, `enzN` and `enzC` are 1 on
//!   every row of a fully specific search, and the benchmarks lose
//!   identifications (DOCS.md, PIN columns).
//! - Unknown / non-builtin enzymes return `true` for any boundary —
//!   matching OpenMS's default "else" branch and Percolator's
//!   unspecific-cleavage semantics.

use model::enzyme::Enzyme;
use search::candidate_gen::Candidate;

#[inline]
fn is_protein_boundary(c: u8) -> bool {
    c == b'-' || c == b'_'
}

/// Returns `true` when the boundary between residues `n` and `c` is
/// consistent with the enzyme's cleavage rule (OpenMS PercolatorInfile).
pub(crate) fn is_enzymatic_boundary(n: u8, c: u8, enzyme: Enzyme) -> bool {
    // Protein boundaries are always enzymatic (usual short-circuit),
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
        Enzyme::AlphaLP => matches!(n, b'T' | b'A' | b'S' | b'V'),
        // Unknown enzymes: the common default branch returns true so
        // unspecific searches don't penalise every PSM as non-enzymatic.
        Enzyme::NoCleavage | Enzyme::NonSpecific => true,
    }
}

/// [`is_enzymatic_boundary`] for a digest of `primary` plus the `extra`
/// enzymes of a multi-protease `--enzyme a,b,...`: the boundary is enzymatic
/// when ANY of them cleaves it.
fn is_enzymatic_boundary_any(n: u8, c: u8, primary: Enzyme, extra: &[Enzyme]) -> bool {
    std::iter::once(&primary)
        .chain(extra)
        .any(|&e| is_enzymatic_boundary(n, c, e))
}

/// Count internal boundaries consistent with the cleavage rule of any of the
/// digest's enzymes (OpenMS PercolatorInfile, per enzyme).
///
/// For an empty / single-residue peptide returns `0` (no internal
/// boundaries to evaluate). For an "unknown" enzyme (universal-true
/// branch above) this returns `len - 1`.
fn count_internal_enzymatic_any(residues: &[u8], primary: Enzyme, extra: &[Enzyme]) -> i32 {
    residues
        .windows(2)
        .filter(|w| is_enzymatic_boundary_any(w[0], w[1], primary, extra))
        .count() as i32
}

/// The `enzN`, `enzC` and `enzInt` PIN features of `cand` under a digest of
/// `primary` plus the `extra` enzymes.
///
/// `enzN` is the boundary between the residue before the peptide and its
/// first residue. A protein N-terminal peptide starts at a terminus: at the
/// protein's own start `pre` is `_` already, and after the initiator Met was
/// removed `pre` is that `M`, which no cleavage rule accepts. Both count as
/// enzymatic, as in the search's own cleavage credit
/// (`search::match_engine::enzyme_cleavage_credit`). Decoys get the same
/// Met-cleaved walk as targets, so this changes both alike.
pub(crate) fn enzyme_features(
    cand: &Candidate,
    primary: Enzyme,
    extra: &[Enzyme],
) -> (u8, u8, i32) {
    let residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
    let first = residues.first().copied().unwrap_or(b'-');
    let last = residues.last().copied().unwrap_or(b'-');
    let pre = if cand.is_protein_n_term {
        b'_'
    } else {
        cand.peptide.pre
    };
    let enz_n = is_enzymatic_boundary_any(pre, first, primary, extra) as u8;
    let enz_c = is_enzymatic_boundary_any(last, cand.peptide.post, primary, extra) as u8;
    let enz_int = count_internal_enzymatic_any(&residues, primary, extra);
    (enz_n, enz_c, enz_int)
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
            Enzyme::Trypsin,
            Enzyme::Chymotrypsin,
            Enzyme::LysC,
            Enzyme::LysN,
            Enzyme::GluC,
            Enzyme::ArgC,
            Enzyme::AspN,
            Enzyme::AlphaLP,
            Enzyme::NoCleavage,
            Enzyme::NonSpecific,
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
        assert!(is_enzymatic_boundary(b'A', b'A', Enzyme::NonSpecific));
        // NoCleavage follows the unknown-enzyme → true convention (common default).
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
        let count = count_internal_enzymatic_any(b"ABKAR", Enzyme::Trypsin, &[]);
        assert_eq!(count, 1);
    }

    #[test]
    fn count_internal_zero_for_short_peptide() {
        assert_eq!(count_internal_enzymatic_any(b"", Enzyme::Trypsin, &[]), 0);
        assert_eq!(count_internal_enzymatic_any(b"A", Enzyme::Trypsin, &[]), 0);
    }

    #[test]
    fn count_internal_handles_p_block() {
        // KPKA: boundaries at i=1,2,3: (K,P), (P,K), (K,A)
        //   trypsin: (K,P) blocked, (P,K) no K/R before, (K,A) yes → count=1.
        assert_eq!(
            count_internal_enzymatic_any(b"KPKA", Enzyme::Trypsin, &[]),
            1
        );
    }

    #[test]
    fn count_internal_universal_returns_len_minus_one() {
        assert_eq!(
            count_internal_enzymatic_any(b"ABCDE", Enzyme::NonSpecific, &[]),
            4
        );
    }

    fn candidate(seq: &[u8], pre: u8, post: u8, n_term: bool, decoy: bool) -> Candidate {
        use model::amino_acid::AminoAcid;
        use model::peptide::Peptide;
        let residues = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Candidate {
            peptide: Peptide::new(residues, pre, post),
            protein_index: 0,
            start_offset_in_protein: usize::from(pre != b'_'),
            is_decoy: decoy,
            is_protein_n_term: n_term,
            is_protein_c_term: post == b'-',
        }
    }

    #[test]
    fn met_clipped_protein_n_term_peptide_starts_at_a_terminus() {
        // M.ASLVTEGLK.P after initiator-Met removal: M|A is no trypsin site,
        // but the peptide starts at the protein's N-terminus. Target and
        // decoy alike.
        for decoy in [false, true] {
            let met_clipped = candidate(b"ASLVTEGLK", b'M', b'P', true, decoy);
            assert_eq!(
                enzyme_features(&met_clipped, Enzyme::Trypsin, &[]),
                (1, 0, 0),
                "decoy={decoy}"
            );
        }
        // The same span after an internal M is not at a terminus.
        let internal = candidate(b"ASLVTEGLK", b'M', b'P', false, false);
        assert_eq!(enzyme_features(&internal, Enzyme::Trypsin, &[]), (0, 0, 0));
        // At the protein's own start `pre` is `_` either way.
        let at_start = candidate(b"MASLVTEGLK", b'_', b'P', true, false);
        assert_eq!(enzyme_features(&at_start, Enzyme::Trypsin, &[]).0, 1);
        // An N-terminal cleaver (Lys-N cuts before K) never cuts M|A either.
        let lysn = candidate(b"ASLVTEGL", b'M', b'K', true, false);
        assert_eq!(enzyme_features(&lysn, Enzyme::LysN, &[]), (1, 1, 0));
    }

    #[test]
    fn multi_protease_boundaries_use_the_union_of_enzymes() {
        // Trypsin + GluC: an E| boundary is enzymatic only through GluC.
        let glu = [Enzyme::GluC];
        assert!(!is_enzymatic_boundary_any(b'E', b'A', Enzyme::Trypsin, &[]));
        assert!(is_enzymatic_boundary_any(b'E', b'A', Enzyme::Trypsin, &glu));
        assert!(is_enzymatic_boundary_any(b'K', b'A', Enzyme::Trypsin, &glu));
        assert!(!is_enzymatic_boundary_any(
            b'A',
            b'A',
            Enzyme::Trypsin,
            &glu
        ));
        // AKAEAKP: (K,A) trypsin, (E,A) GluC, (K,P) neither.
        assert_eq!(
            count_internal_enzymatic_any(b"AKAEAKP", Enzyme::Trypsin, &[]),
            1
        );
        assert_eq!(
            count_internal_enzymatic_any(b"AKAEAKP", Enzyme::Trypsin, &glu),
            2
        );
    }
}
