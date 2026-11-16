//! Link-site rules: where on a peptide the linker can sit.

use crate::crosslinker::LinkEnd;

/// A residue position the linker can occupy on a peptide.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LinkSite {
    /// 0-based residue index in the peptide.
    pub pos: usize,
    /// Whether the site is the N-terminal amine (on residue 0).
    pub n_term: bool,
    /// Whether the site is the protein C-terminal carboxyl (on the last
    /// residue).
    pub c_term: bool,
}

/// Every site of a peptide with residues `residues` that linker end `end` can
/// occupy.
///
/// - A reactive residue that carries no other modification (`modified[i]` is
///   false). A linked Lys blocks tryptic cleavage, so a reactive residue at the
///   peptide's C-terminus is a site only when it is also the protein's
///   C-terminus.
/// - The protein N-terminus (residue 0 of a protein N-terminal peptide), when
///   the end reacts with it and residue 0 is unmodified; with
///   [`LinkEnd::peptide_n_term`], the N-terminus of every peptide.
/// - The protein C-terminus (the last residue of a protein C-terminal
///   peptide), when the end reacts with it.
///
/// A position that qualifies both ways (a reactive residue 0 of a protein
/// N-terminal peptide, a reactive last residue of a protein C-terminal one) is
/// listed once, as the residue site, because both placements put the same mass
/// on the same residue.
pub fn link_sites(
    residues: &[u8],
    modified: &[bool],
    is_protein_n_term: bool,
    is_protein_c_term: bool,
    end: &LinkEnd,
) -> Vec<LinkSite> {
    let n = residues.len();
    let mut out = Vec::new();
    for (i, &r) in residues.iter().enumerate() {
        if !end.reacts_with(r) || modified.get(i).copied().unwrap_or(false) {
            continue;
        }
        if i + 1 == n && !is_protein_c_term {
            continue;
        }
        out.push(LinkSite {
            pos: i,
            n_term: false,
            c_term: false,
        });
    }
    if (end.peptide_n_term || (end.protein_n_term && is_protein_n_term))
        && n > 0
        && !modified.first().copied().unwrap_or(false)
        && out.first().is_none_or(|s| s.pos != 0)
    {
        out.insert(
            0,
            LinkSite {
                pos: 0,
                n_term: true,
                c_term: false,
            },
        );
    }
    if end.protein_c_term
        && is_protein_c_term
        && n > 0
        && out.last().is_none_or(|s| s.pos != n - 1 || s.n_term)
    {
        out.push(LinkSite {
            pos: n - 1,
            n_term: false,
            c_term: true,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::Crosslinker;

    fn sites(seq: &[u8], nterm: bool, cterm: bool) -> Vec<(usize, bool)> {
        let unmod = vec![false; seq.len()];
        link_sites(seq, &unmod, nterm, cterm, &Crosslinker::dsso().ends[0])
            .into_iter()
            .map(|s| (s.pos, s.n_term))
            .collect()
    }

    #[test]
    fn internal_lysines_are_sites_c_terminal_ones_are_not() {
        assert_eq!(sites(b"AKLLKR", false, false), vec![(1, false), (4, false)]);
        assert_eq!(sites(b"AKLLK", false, false), vec![(1, false)]);
        // At the protein C-terminus the last Lys can be linked.
        assert_eq!(sites(b"AKLLK", false, true), vec![(1, false), (4, false)]);
    }

    #[test]
    fn protein_n_terminus_is_a_site() {
        assert_eq!(sites(b"MAKR", true, false), vec![(0, true), (2, false)]);
        assert_eq!(sites(b"MAKR", false, false), vec![(2, false)]);
        // Lys at residue 0 of a protein N-terminal peptide: one site.
        assert_eq!(sites(b"KAR", true, false), vec![(0, false)]);
    }

    #[test]
    fn modified_residues_are_not_sites() {
        let x = Crosslinker::dsso();
        let s = link_sites(
            b"MAKAKR",
            &[true, false, true, false, false, false],
            true,
            false,
            &x.ends[0],
        );
        assert_eq!(
            s,
            vec![LinkSite {
                pos: 4,
                n_term: false,
                c_term: false,
            }]
        );
    }

    #[test]
    fn peptide_n_terminus_is_a_site_when_the_end_says_so() {
        let mut end = Crosslinker::dsso().ends[0].clone();
        end.peptide_n_term = true;
        end.protein_n_term = false;
        let s: Vec<(usize, bool)> = link_sites(b"MAKR", &[false; 4], false, false, &end)
            .into_iter()
            .map(|s| (s.pos, s.n_term))
            .collect();
        assert_eq!(s, vec![(0, true), (2, false)]);
    }

    #[test]
    fn protein_c_terminus_is_a_site_of_a_carboxyl_end() {
        let end = LinkEnd {
            residues: b"DE".to_vec(),
            protein_c_term: true,
            ..LinkEnd::default()
        };
        let s = |cterm: bool| -> Vec<(usize, bool)> {
            link_sites(b"ADAK", &[false; 4], false, cterm, &end)
                .into_iter()
                .map(|s| (s.pos, s.c_term))
                .collect()
        };
        assert_eq!(s(false), vec![(1, false)]);
        assert_eq!(s(true), vec![(1, false), (3, true)]);
        // A reactive last residue is listed once, as the residue site.
        let e = link_sites(b"AKAD", &[false; 4], false, true, &end);
        assert_eq!(e.len(), 1);
        assert!(!e[0].c_term && e[0].pos == 3);
    }
}
