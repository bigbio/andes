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

/// Boundary sequon: the N-X-S/T motif completed by the peptide's C-terminal
/// FLANKING residue `post` (the next residue in the protein). Returns `Some(n-2)`
/// — the glycosite N position inside the peptide — when `residues[n-2]=N`,
/// `residues[n-1]≠P`, and `post ∈ {S,T}`. A fully-tryptic peptide `…N-K` whose
/// completing S/T is the first residue of the next tryptic peptide carries the
/// glycan on that N but is invisible to a peptide-internal-only sequon scan
/// (round-4 candidate-gen audit: silently never generated). `None` otherwise.
pub fn boundary_nxst_site(residues: &[u8], post: u8) -> Option<usize> {
    let n = residues.len();
    if n < 2 {
        return None;
    }
    if residues[n - 2] == b'N' && residues[n - 1] != b'P' && (post == b'S' || post == b'T') {
        Some(n - 2)
    } else {
        None
    }
}

/// Returns the 0-based indices of the N in EVERY N-X-S/T sequon (X ≠ P), in
/// ascending order. Used to place the glycan at each candidate glycosite when a
/// peptide carries more than one sequon (~8% of tryptic N-glycopeptides), so the
/// c/z hyperscore can be maximized over sites instead of assuming the first —
/// [`first_nxst_site`] silently mis-places the glycan on multi-sequon peptides.
pub fn all_nxst_sites(residues: &[u8]) -> Vec<usize> {
    let n = residues.len();
    if n < 3 {
        return Vec::new();
    }
    (0..n.saturating_sub(2))
        .filter(|&i| {
            residues[i] == b'N'
                && residues[i + 1] != b'P'
                && (residues[i + 2] == b'S' || residues[i + 2] == b'T')
        })
        .collect()
}

/// Every glycosite of a peptide in its protein context: the N of each internal
/// N-X-S/T sequon ([`all_nxst_sites`]) followed by the boundary-sequon N
/// ([`boundary_nxst_site`]) when the residue after the peptide, `post`,
/// completes one. The two never coincide: an internal sequon's N sits at most
/// at `len-3`, the boundary N at `len-2`.
pub fn glycosites(residues: &[u8], post: u8) -> Vec<usize> {
    let mut sites = all_nxst_sites(residues);
    sites.extend(boundary_nxst_site(residues, post));
    sites
}

// ── O-glycosylation site helpers ─────────────────────────────────────────────
//
// O-glycans attach to Ser/Thr with NO sequon (unlike N-glycans' N-X-S/T). A
// peptide is an O-glycan candidate iff it contains any S or T; every S/T is a
// possible glycosite. There is no Pro exclusion and no boundary case (the site is
// always internal to the peptide).

/// Returns `true` iff `residues` contains at least one Ser or Thr — the O-glycan
/// analogue of [`has_nxst_sequon`] (no motif, no Pro exclusion).
pub fn has_ser_or_thr(residues: &[u8]) -> bool {
    residues.iter().any(|&r| r == b'S' || r == b'T')
}

/// Returns the 0-based index of the FIRST Ser or Thr, or `None` if there is none.
pub fn first_ser_thr_site(residues: &[u8]) -> Option<usize> {
    residues.iter().position(|&r| r == b'S' || r == b'T')
}

/// Returns the 0-based indices of EVERY Ser/Thr, ascending. Every S/T is a candidate
/// glycosite; the c/z hyperscore is maximised over them exactly as N-glycan does
/// over its sequons.
pub fn all_ser_thr_sites(residues: &[u8]) -> Vec<usize> {
    residues
        .iter()
        .enumerate()
        .filter_map(|(i, &r)| (r == b'S' || r == b'T').then_some(i))
        .collect()
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
    fn boundary_sequon_uses_c_flank() {
        // …N-K | S-…  : N at len-2, X=K≠P, post=S → glycosite at n-2.
        assert_eq!(boundary_nxst_site(b"AAANK", b'S'), Some(3));
        assert_eq!(boundary_nxst_site(b"AAANK", b'T'), Some(3));
        assert_eq!(boundary_nxst_site(b"AAANK", b'A'), None); // post not S/T
        assert_eq!(boundary_nxst_site(b"AAANP", b'S'), None); // X=P excluded
        assert_eq!(boundary_nxst_site(b"AAAKN", b'S'), None); // last is N, n-2 is K
        assert_eq!(boundary_nxst_site(b"N", b'S'), None); // too short
    }

    #[test]
    fn glycosites_include_the_boundary_sequon() {
        // Boundary only: …N-K | S-…
        assert_eq!(glycosites(b"AAANK", b'S'), vec![3]);
        // Internal plus boundary: two sites.
        assert_eq!(glycosites(b"ANLTANK", b'T'), vec![1, 5]);
        // A flank that does not complete the sequon adds nothing.
        assert_eq!(glycosites(b"ANLTANK", b'A'), vec![1]);
        // Protein C-terminus.
        assert_eq!(glycosites(b"AAANK", b'-'), Vec::<usize>::new());
    }

    #[test]
    fn all_sites_multi_and_single() {
        assert_eq!(all_nxst_sites(b"SVNLTK"), vec![2]); // single sequon
        assert_eq!(all_nxst_sites(b"AANDSNKTQ"), vec![2, 5]); // two sequons: NDS + NKT
        assert_eq!(all_nxst_sites(b"SVNPLTK"), Vec::<usize>::new()); // X=P excluded
        assert_eq!(all_nxst_sites(b"AAA"), Vec::<usize>::new()); // none
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

    #[test]
    fn ser_thr_helpers() {
        assert!(has_ser_or_thr(b"SVTK"));
        assert!(has_ser_or_thr(b"AKLT")); // T counts
        assert!(!has_ser_or_thr(b"AAAK")); // no S/T
        assert!(!has_ser_or_thr(b""));

        assert_eq!(first_ser_thr_site(b"AKST"), Some(2)); // first S
        assert_eq!(first_ser_thr_site(b"AKTA"), Some(2)); // first T
        assert_eq!(first_ser_thr_site(b"AAAA"), None);

        assert_eq!(all_ser_thr_sites(b"ASKT"), vec![1, 3]); // S and T
        assert_eq!(all_ser_thr_sites(b"AAAA"), Vec::<usize>::new());
    }
}
