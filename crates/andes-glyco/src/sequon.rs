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

/// Returns the 0-based index of the N in the FIRST N-X-S/T sequon (X ≠ P), or
/// `None` if there is none. This is the presumptive glycosylation site — used to
/// place the intact glycan on glycosite-spanning ETD c/z fragments.
pub fn first_nxst_site(residues: &[u8]) -> Option<usize> {
    let n = residues.len();
    if n < 3 {
        return None;
    }
    (0..n.saturating_sub(2)).find(|&i| {
        residues[i] == b'N'
            && residues[i + 1] != b'P'
            && (residues[i + 2] == b'S' || residues[i + 2] == b'T')
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_site_index() {
        assert_eq!(first_nxst_site(b"SVNLTK"), Some(2));
        assert_eq!(first_nxst_site(b"SVNPLTK"), None);
        assert_eq!(first_nxst_site(b"AANDSNKTQ"), Some(2)); // first of two sequons
    }

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
