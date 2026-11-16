use crate::candidate_gen::Candidate;
use crate::psm::PsmMatch;
use model::IsomerClass;

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct IsomerKey {
    precursor_mz_bits: u64,
    charge: u8,
    is_decoy: bool,
    class: IsomerClass,
}

impl IsomerKey {
    pub fn new(spectrum_precursor_mz: f64, psm: &PsmMatch, cand: &Candidate) -> Self {
        Self {
            precursor_mz_bits: psm
                .precursor_mz_override
                .unwrap_or(spectrum_precursor_mz)
                .to_bits(),
            charge: psm.charge_used,
            is_decoy: cand.is_decoy,
            class: IsomerClass::of_peptide(&cand.peptide),
        }
    }
}

pub fn collapse_isomers<'a, I>(
    spectrum_precursor_mz: f64,
    psms_best_first: I,
    candidates: &[Candidate],
) -> Vec<(usize, Vec<usize>)>
where
    I: IntoIterator<Item = &'a PsmMatch>,
{
    let mut keys: Vec<IsomerKey> = Vec::new();
    let mut kept: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, psm) in psms_best_first.into_iter().enumerate() {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let key = IsomerKey::new(spectrum_precursor_mz, psm, cand);
        match keys.iter().position(|k| *k == key) {
            Some(slot) => kept[slot].1.push(i),
            None => {
                keys.push(key);
                kept.push((i, Vec::new()));
            }
        }
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::{AminoAcid, ModLocation, Modification, Peptide, ResidueSpec};

    fn aa(c: u8) -> AminoAcid {
        AminoAcid::standard(c).unwrap()
    }

    fn modded(c: u8, delta: f64) -> AminoAcid {
        let mut a = aa(c);
        a.mod_ = Some(model::modification::leak_mod(Modification {
            name: "m".into(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(c),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }));
        a
    }

    fn pep(residues: Vec<AminoAcid>) -> Peptide {
        Peptide::new(residues, b'K', b'A')
    }

    fn cand(peptide: Peptide, protein_index: usize, is_decoy: bool) -> Candidate {
        Candidate {
            peptide,
            protein_index,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn psm(idx: u32, charge: u8) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![idx],
            charge_used: charge,
            mass_error_ppm: 0.0,
            score: 0.0,
            rank_score: 0.0,
            edge_score: 0,
            activation_method: None,
            features: Default::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        }
    }

    #[test]
    fn collapse_keeps_the_first_of_each_class() {
        let phospho = 79.96633;
        let s1 = || pep(vec![modded(b'S', phospho), aa(b'T'), aa(b'K')]);
        let s2 = || pep(vec![aa(b'S'), modded(b'T', phospho), aa(b'K')]);
        let cands = vec![
            cand(s1(), 0, false),
            cand(s2(), 0, false),
            cand(s2(), 1, true),
            cand(s2(), 2, false),
            cand(pep(vec![aa(b'S'), aa(b'T'), aa(b'K')]), 0, false),
        ];
        let psms = vec![
            psm(0, 2),
            psm(1, 2),
            psm(2, 2),
            psm(3, 2),
            psm(4, 2),
            psm(1, 3),
        ];
        let kept = collapse_isomers(500.0, &psms, &cands);
        assert_eq!(
            kept,
            vec![(0, vec![1, 3]), (2, vec![]), (4, vec![]), (5, vec![])]
        );
    }
}
