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

pub fn glycosites(residues: &[u8], post: u8) -> Vec<usize> {
    let mut sites = all_nxst_sites(residues);
    sites.extend(boundary_nxst_site(residues, post));
    sites
}

pub fn has_ser_or_thr(residues: &[u8]) -> bool {
    residues.iter().any(|&r| r == b'S' || r == b'T')
}

pub fn first_ser_thr_site(residues: &[u8]) -> Option<usize> {
    residues.iter().position(|&r| r == b'S' || r == b'T')
}

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
        assert_eq!(first_nxst_site(b"AANDSNKTQ"), Some(2));
    }

    #[test]
    fn boundary_sequon_uses_c_flank() {
        assert_eq!(boundary_nxst_site(b"AAANK", b'S'), Some(3));
        assert_eq!(boundary_nxst_site(b"AAANK", b'T'), Some(3));
        assert_eq!(boundary_nxst_site(b"AAANK", b'A'), None);
        assert_eq!(boundary_nxst_site(b"AAANP", b'S'), None);
        assert_eq!(boundary_nxst_site(b"AAAKN", b'S'), None);
        assert_eq!(boundary_nxst_site(b"N", b'S'), None);
    }

    #[test]
    fn glycosites_include_the_boundary_sequon() {
        assert_eq!(glycosites(b"AAANK", b'S'), vec![3]);
        assert_eq!(glycosites(b"ANLTANK", b'T'), vec![1, 5]);
        assert_eq!(glycosites(b"ANLTANK", b'A'), vec![1]);
        assert_eq!(glycosites(b"AAANK", b'-'), Vec::<usize>::new());
    }

    #[test]
    fn all_sites_multi_and_single() {
        assert_eq!(all_nxst_sites(b"SVNLTK"), vec![2]);
        assert_eq!(all_nxst_sites(b"AANDSNKTQ"), vec![2, 5]);
        assert_eq!(all_nxst_sites(b"SVNPLTK"), Vec::<usize>::new());
        assert_eq!(all_nxst_sites(b"AAA"), Vec::<usize>::new());
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
        assert!(has_ser_or_thr(b"AKLT"));
        assert!(!has_ser_or_thr(b"AAAK"));
        assert!(!has_ser_or_thr(b""));

        assert_eq!(first_ser_thr_site(b"AKST"), Some(2));
        assert_eq!(first_ser_thr_site(b"AKTA"), Some(2));
        assert_eq!(first_ser_thr_site(b"AAAA"), None);

        assert_eq!(all_ser_thr_sites(b"ASKT"), vec![1, 3]);
        assert_eq!(all_ser_thr_sites(b"AAAA"), Vec::<usize>::new());
    }
}
