// Clean-room N-glycosylation sequon filter.
//
// The canonical N-glycosylation sequon is N-X-S/T where X is any amino acid
// except proline. This module exposes a single predicate operating on raw
// one-letter amino-acid bytes.

/// Returns `true` iff `residues` contains at least one N-X-S/T sequon
/// (X ≠ P).
///
/// `residues` must be a slice of uppercase one-letter amino-acid bytes (e.g.
/// `b"SVNLTK"`). The function performs a single linear scan and never
/// allocates.
pub fn has_nxst_sequon(residues: &[u8]) -> bool {
    let n = residues.len();
    if n < 3 {
        return false;
    }
    for i in 0..n.saturating_sub(2) {
        if residues[i] == b'N'
            && residues[i + 1] != b'P'
            && (residues[i + 2] == b'S' || residues[i + 2] == b'T')
        {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn svnltk_has_sequon() {
        assert!(has_nxst_sequon(b"SVNLTK"));
    }

    #[test]
    fn svnpltk_no_sequon_x_is_p() {
        assert!(!has_nxst_sequon(b"SVNPLTK"));
    }

    #[test]
    fn peptide_no_sequon() {
        assert!(!has_nxst_sequon(b"PEPTIDE"));
    }

    #[test]
    fn nst_minimal_positive() {
        assert!(has_nxst_sequon(b"NST"));
    }

    #[test]
    fn npt_no_sequon_x_is_p() {
        assert!(!has_nxst_sequon(b"NPT"));
    }

    #[test]
    fn nn_too_short_for_sequon() {
        assert!(!has_nxst_sequon(b"NN"));
    }

    #[test]
    fn empty_slice_no_sequon() {
        assert!(!has_nxst_sequon(b""));
    }
}
