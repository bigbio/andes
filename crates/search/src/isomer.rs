//! Positional-isomer identity of a PSM.
//!
//! Two PSMs on one spectrum are positional isomers when they carry the same
//! residues and the same multiset of modification mass deltas, with the
//! modifications on different sites. Site alternatives of one peptide score
//! almost identically, so the per-spectrum rows keep one of them
//! ([`collapse_isomers`]).

use crate::candidate_gen::Candidate;
use crate::psm::PsmMatch;

/// A modification mass delta in units of 1e-5 Da, the precision of the PIN
/// `Peptide` column.
#[inline]
fn mod_delta_units(mass_delta: f64) -> i64 {
    (mass_delta * 100_000.0).round() as i64
}

/// Identity of a PSM up to the placement of its modifications: the precursor it
/// was matched at, its charge, its target/decoy side, its residues, and the
/// sorted multiset of its modification mass deltas (terminal mods included,
/// since they are carried on the terminal residue). Two PSMs with equal keys are
/// positional isomers — the same peptide with the same modifications on
/// different sites.
///
/// A target and a decoy never share a key, so a decoy that repeats a target's
/// residues keeps its own row and its `-1` label.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct IsomerKey {
    precursor_mz_bits: u64,
    charge: u8,
    is_decoy: bool,
    residues: Vec<u8>,
    mod_deltas: Vec<i64>,
}

impl IsomerKey {
    /// The key of `psm` (whose primary candidate is `cand`) on a spectrum with
    /// precursor m/z `spectrum_precursor_mz`. A chimeric secondary is keyed at
    /// its own precursor (`precursor_mz_override`).
    pub fn new(spectrum_precursor_mz: f64, psm: &PsmMatch, cand: &Candidate) -> Self {
        let residues = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
        let mut mod_deltas: Vec<i64> = cand
            .peptide
            .residues
            .iter()
            .filter_map(|aa| aa.mod_.map(|m| mod_delta_units(m.mass_delta)))
            .collect();
        mod_deltas.sort_unstable();
        Self {
            precursor_mz_bits: psm
                .precursor_mz_override
                .unwrap_or(spectrum_precursor_mz)
                .to_bits(),
            charge: psm.charge_used,
            is_decoy: cand.is_decoy,
            residues,
            mod_deltas,
        }
    }
}

/// The PSMs that represent each positional-isomer class of one spectrum, in
/// the order given: of the PSMs sharing an [`IsomerKey`], the first one is
/// kept and the rest are dropped. `psms_best_first` must be in rank order
/// (best first), so the kept PSM is the best-ranked isomer.
///
/// Returns, for each kept PSM, its index in `psms_best_first` and the indices
/// of the PSMs collapsed into it (empty when it has no isomer).
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
            cand(s2(), 0, false), // isomer of 0
            cand(s2(), 1, true),  // decoy: its own class
            cand(s2(), 2, false), // isomer of 0 in another protein
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
