use super::*;
use crate::psm::PsmMatch;
use model::amino_acid::AminoAcid;
use model::modification::{ModLocation, Modification};
use model::peptide::Peptide;
use model::ResidueSpec;

fn seq_peptide(bytes: &[u8]) -> Peptide {
    let residues: Vec<AminoAcid> = bytes
        .iter()
        .filter_map(|&b| AminoAcid::standard(b))
        .collect();
    Peptide::new(residues, b'R', b'K')
}

fn cand_with_peptide(peptide: Peptide) -> Candidate {
    Candidate {
        peptide,
        protein_index: 0,
        start_offset_in_protein: 0,
        is_decoy: false,
        is_protein_n_term: false,
        is_protein_c_term: false,
    }
}

fn credit(cand: &Candidate, enz: Enzyme) -> i32 {
    enzyme_cleavage_credit(
        cand,
        enz,
        enz.is_c_term(),
        enz.is_n_term(),
        10,
        -100,
        1,
        -1000,
    )
}

fn flanked(seq: &[u8], pre: u8, post: u8) -> Candidate {
    let residues: Vec<AminoAcid> = seq.iter().filter_map(|&b| AminoAcid::standard(b)).collect();
    cand_with_peptide(Peptide::new(residues, pre, post))
}

#[test]
fn n_term_enzyme_peptide_credit_reads_the_first_residue() {
    assert_eq!(credit(&flanked(b"DAAAG", b'G', b'D'), Enzyme::AspN), 11);
    assert_eq!(credit(&flanked(b"AAAAG", b'D', b'D'), Enzyme::AspN), -990);
    assert_eq!(credit(&flanked(b"DAAAG", b'G', b'A'), Enzyme::AspN), -99);
    assert_eq!(credit(&flanked(b"KAAAG", b'G', b'K'), Enzyme::LysN), 11);
}

#[test]
fn c_term_enzyme_credit_reads_the_last_residue() {
    assert_eq!(credit(&flanked(b"AAAAR", b'K', b'G'), Enzyme::Trypsin), 11);
    assert_eq!(
        credit(&flanked(b"AAAAG", b'G', b'G'), Enzyme::Trypsin),
        -1100
    );
}

fn psm(primary: u32, rank: f32, pin: f32) -> PsmMatch {
    PsmMatch {
        spectrum_idx: 0,
        candidate_idxs: smallvec::smallvec![primary],
        charge_used: 2,
        mass_error_ppm: 0.0,
        score: pin,
        rank_score: rank,
        edge_score: (rank - pin) as i32,
        activation_method: None,
        features: Default::default(),
        isotope_offset: 0,
        precursor_mz_override: None,
        site_evidence: None,
        tie_key: 0,
    }
}

#[test]
fn collision_decoy_relabeled_as_target_genuine_decoy_kept() {
    let cand = |seq: &[u8], is_decoy: bool| Candidate {
        peptide: seq_peptide(seq),
        protein_index: 0,
        start_offset_in_protein: 0,
        is_decoy,
        is_protein_n_term: false,
        is_protein_c_term: false,
    };
    let mut cands = vec![
        cand(b"PEPTIDEK", false),
        cand(b"PEPTIDEK", true),
        cand(b"EDITPEPK", true),
        cand(b"PEPTLDEK", true),
    ];
    let targets = folded_set(vec![b"PEPTIDEK".as_slice()]);
    relabel_collision_decoys_with(&mut cands, &targets);
    assert!(!cands[0].is_decoy, "target unchanged");
    assert!(!cands[1].is_decoy, "collision decoy relabeled to target");
    assert!(cands[2].is_decoy, "genuine decoy left as decoy");
    assert!(
        !cands[3].is_decoy,
        "a decoy equal to a target with I and L folded is relabeled to target"
    );
}

#[test]
fn dedup_uses_rank_score_not_pin_score() {
    let pep = seq_peptide(b"PEPTK");
    let cands = vec![cand_with_peptide(pep.clone())];
    let psms = vec![psm(0, 100.4, 99.6), psm(0, 120.0, 99.6)];
    let out = dedup_pepseq_score(psms, &cands);
    assert_eq!(out.len(), 2, "different rank_score keys must not merge");
}

#[test]
fn dedup_distinguishes_mod_state() {
    let mut ox = seq_peptide(b"PEPMK");
    ox.residues[3].mod_ = Some(model::modification::leak_mod(Modification {
        name: "Ox".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
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
    let cands = vec![
        cand_with_peptide(seq_peptide(b"PEPMK")),
        cand_with_peptide(ox),
    ];
    let psms = vec![psm(0, 50.0, 50.0), psm(1, 50.0, 50.0)];
    let out = dedup_pepseq_score(psms, &cands);
    assert_eq!(out.len(), 2, "mod-aware pepSeq keys must differ");
}

#[test]
fn dedup_keeps_highest_rank_score_survivor() {
    let pep = seq_peptide(b"PEPTK");
    let cands = vec![cand_with_peptide(pep)];
    let psms = vec![psm(0, 59.6, 50.0), psm(0, 60.4, 50.0)];
    let out = dedup_pepseq_score(psms, &cands);
    assert_eq!(out.len(), 1);
    assert!((out[0].rank_score - 60.4).abs() < f32::EPSILON);
}
