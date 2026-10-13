//! `RtPeptide::from_andes` on peptide strings produced by andes's own writer
//! (`model::Peptide`'s `Display`, the PIN peptide column). The strings and
//! their expected ProForma are kept in `fixtures/andes_writer_peptides.tsv`,
//! which `scripts/student/test_rt_parsing.py` reads for the Python parser.

use std::path::PathBuf;

use model::{AminoAcid, AminoAcidSet, AminoAcidSetBuilder, ModLocation, Modification, Peptide};
use onnx_infer::rt_native::RtPeptide;

fn aa_set(fixed: &[&str], variable: &[&str]) -> AminoAcidSet {
    let mut b = AminoAcidSetBuilder::new_standard();
    for line in fixed {
        b = b.add_fixed_mod(Modification::from_mods_txt_line(line).unwrap());
    }
    for line in variable {
        b = b.add_variable_mod(Modification::from_mods_txt_line(line).unwrap());
    }
    b.build().unwrap()
}

/// `seq` with residue 1 taken from the `first_loc` variants whose delta is
/// `first_delta`, and every other residue in its fixed-mod (Anywhere) form.
fn written(
    set: &AminoAcidSet,
    seq: &[u8],
    pre: u8,
    post: u8,
    first_loc: ModLocation,
    first_delta: f64,
) -> String {
    let first = set
        .variants_for(seq[0], first_loc)
        .iter()
        .find(|aa| {
            aa.mod_
                .as_ref()
                .is_some_and(|m| (m.mass_delta - first_delta).abs() < 1e-3)
        })
        .unwrap_or_else(|| panic!("no {first_delta} variant of {}", seq[0] as char))
        .clone();
    let mut residues: Vec<AminoAcid> = vec![first];
    for &r in &seq[1..] {
        let aa = set.variants_for(r, ModLocation::Anywhere)[0].clone();
        residues.push(aa);
    }
    Peptide::new(residues, pre, post).to_string()
}

fn writer_cases() -> Vec<(String, &'static str)> {
    let cam_acetyl = aa_set(
        &["57.021464,C,fix,any,Carbamidomethyl"],
        &["42.010565,*,opt,Prot-N-term,Acetyl"],
    );
    let tmtpro_acetyl = aa_set(
        &["304.207146,K,fix,any,TMTpro"],
        &["42.010565,*,opt,Prot-N-term,Acetyl"],
    );
    vec![
        (
            written(
                &cam_acetyl,
                b"CPEPTIDEK",
                b'_',
                b'R',
                ModLocation::ProtNTerm,
                99.032029,
            ),
            "[UNIMOD:1]-C[UNIMOD:4]PEPTIDEK",
        ),
        (
            written(
                &cam_acetyl,
                b"CPEPCTIDEK",
                b'M',
                b'R',
                ModLocation::ProtNTerm,
                99.032029,
            ),
            "[UNIMOD:1]-C[UNIMOD:4]PEPC[UNIMOD:4]TIDEK",
        ),
        (
            written(
                &tmtpro_acetyl,
                b"KPEPTIDEK",
                b'_',
                b'R',
                ModLocation::ProtNTerm,
                346.217711,
            ),
            "[UNIMOD:1]-K[UNIMOD:2016]PEPTIDEK[UNIMOD:2016]",
        ),
    ]
}

fn fixture() -> Vec<(String, String)> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("fixtures/andes_writer_peptides.tsv");
    std::fs::read_to_string(path)
        .expect("read fixture")
        .lines()
        .skip(1)
        .map(|l| {
            let (a, p) = l.split_once('\t').expect("two columns");
            (a.to_string(), p.to_string())
        })
        .collect()
}

#[test]
fn fixture_holds_the_strings_the_andes_writer_produces() {
    let cases: Vec<(String, String)> = writer_cases()
        .into_iter()
        .map(|(a, p)| (a, p.to_string()))
        .collect();
    assert_eq!(fixture(), cases);
}

#[test]
fn from_andes_decodes_summed_terminal_and_residue_deltas() {
    for (andes, proforma) in writer_cases() {
        let p = RtPeptide::from_andes(&andes)
            .unwrap_or_else(|| panic!("from_andes rejects writer output {andes}"));
        assert_eq!(p.to_proforma(), proforma, "{andes}");
    }
}
