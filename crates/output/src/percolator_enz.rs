use model::enzyme::Enzyme;
use search::candidate_gen::Candidate;

#[inline]
fn is_protein_boundary(c: u8) -> bool {
    c == b'-' || c == b'_'
}

pub(crate) fn is_enzymatic_boundary(n: u8, c: u8, enzyme: Enzyme) -> bool {
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
        Enzyme::NoCleavage | Enzyme::NonSpecific => true,
    }
}

fn is_enzymatic_boundary_any(n: u8, c: u8, primary: Enzyme, extra: &[Enzyme]) -> bool {
    std::iter::once(&primary)
        .chain(extra)
        .any(|&e| is_enzymatic_boundary(n, c, e))
}

fn count_internal_enzymatic_any(residues: &[u8], primary: Enzyme, extra: &[Enzyme]) -> i32 {
    residues
        .windows(2)
        .filter(|w| is_enzymatic_boundary_any(w[0], w[1], primary, extra))
        .count() as i32
}

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
        assert!(is_enzymatic_boundary(b'K', b'A', Enzyme::Trypsin));
        assert!(is_enzymatic_boundary(b'R', b'A', Enzyme::Trypsin));
        assert!(!is_enzymatic_boundary(b'K', b'P', Enzyme::Trypsin));
        assert!(!is_enzymatic_boundary(b'R', b'P', Enzyme::Trypsin));
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
        assert!(is_enzymatic_boundary(b'A', b'A', Enzyme::NoCleavage));
    }

    #[test]
    fn count_internal_handles_tryptic_peptide() {
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
        for decoy in [false, true] {
            let met_clipped = candidate(b"ASLVTEGLK", b'M', b'P', true, decoy);
            assert_eq!(
                enzyme_features(&met_clipped, Enzyme::Trypsin, &[]),
                (1, 0, 0),
                "decoy={decoy}"
            );
        }
        let internal = candidate(b"ASLVTEGLK", b'M', b'P', false, false);
        assert_eq!(enzyme_features(&internal, Enzyme::Trypsin, &[]), (0, 0, 0));
        let at_start = candidate(b"MASLVTEGLK", b'_', b'P', true, false);
        assert_eq!(enzyme_features(&at_start, Enzyme::Trypsin, &[]).0, 1);
        let lysn = candidate(b"ASLVTEGL", b'M', b'K', true, false);
        assert_eq!(enzyme_features(&lysn, Enzyme::LysN, &[]), (1, 1, 0));
    }

    #[test]
    fn multi_protease_boundaries_use_the_union_of_enzymes() {
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
