//! Decoy generation on Tryp_Pig_Bov.fasta.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use input::FastaReader;
use search::decoy::{build_search_db, decoy_db, DecoyStrategy, Proteases, DEFAULT_DECOY_SEED};

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("test-fixtures/Tryp_Pig_Bov.fasta")
        .canonicalize()
        .expect("canonicalize Tryp_Pig_Bov.fasta path")
}

/// The fully cleaved tryptic segments of `seq` (cut after every K and R).
fn tryptic_segments(seq: &[u8]) -> Vec<&[u8]> {
    seq.split_inclusive(|&r| r == b'K' || r == b'R').collect()
}

#[test]
fn tryp_pig_bov_pseudo_reverses_every_tryptic_peptide() {
    let path = fixture_path();
    let target = FastaReader::load_all(BufReader::new(File::open(&path).unwrap())).unwrap();
    let decoy = decoy_db(
        &target,
        "XXX",
        DecoyStrategy::Reverse,
        DEFAULT_DECOY_SEED,
        &Proteases::default(),
    );
    assert_eq!(decoy.len(), 16);
    let mut pseudo_reversed = 0;
    let mut segments = 0;
    for (t, d) in target.iter().zip(decoy.iter()) {
        assert_eq!(d.accession, format!("XXX_{}", t.accession));
        assert_eq!(d.description, t.description);
        let (ts, ds) = (tryptic_segments(&t.sequence), tryptic_segments(&d.sequence));
        assert_eq!(ts.len(), ds.len(), "{}", t.accession);
        for (k, (a, b)) in ts.iter().zip(&ds).enumerate() {
            assert_eq!(a.len(), b.len());
            assert_eq!(a.last(), b.last(), "C-terminal residue stays");
            let (mut x, mut y) = (a.to_vec(), b.to_vec());
            x.sort_unstable();
            y.sort_unstable();
            assert_eq!(x, y, "composition stays");
            // Pseudo-reversal: the residues before the last, reversed, after
            // an initiator Met, unless a redraw replaced a palindrome.
            let lead = usize::from(k == 0 && a[0] == b'M');
            let mut want = a[lead..a.len() - 1].to_vec();
            want.reverse();
            segments += 1;
            pseudo_reversed += usize::from(b[lead..b.len() - 1] == want[..]);
        }
    }
    assert!(pseudo_reversed * 10 > segments * 9, "{pseudo_reversed}/{segments}");
}

#[test]
fn tryp_pig_bov_search_db_has_32_proteins() {
    let path = fixture_path();
    let target = FastaReader::load_all(BufReader::new(File::open(&path).unwrap())).unwrap();
    let combined = build_search_db(&target, "XXX", DecoyStrategy::Reverse, DEFAULT_DECOY_SEED);
    assert_eq!(combined.len(), 32);
    for i in 0..16 {
        assert_eq!(combined.proteins[i].accession, target.proteins[i].accession);
        assert!(combined.proteins[16 + i].accession.starts_with("XXX_"));
    }
}
