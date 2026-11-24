//! Near-identical copies of cross-linked peptides in other proteins.
//!
//! On a proteome-sized database the false inter-links a rescorer accepts are mostly a
//! correct peptide with a wrong partner that is a near-copy of the true one (a homolog one
//! or two residues away, or the same sequence up to I/L). A reversed decoy is never that
//! close, so the decoys undercount these errors. An inter-link whose peptide has such a copy
//! in a protein it does not list is estimated apart, as `LinkClass::InterAmbiguous`.

use std::collections::{HashMap, HashSet};
use std::path::Path;

/// Substitutions (same length, I = L) within which a peptide counts as a near-copy. On the
/// PXD042173 truth set searched against human + library, 44 % of the accepted false
/// inter-link partners were within 2 substitutions of the true library peptide, against
/// 18 % of the decoy partners.
pub(crate) const NEAR_COPY_MAX_SUBSTITUTIONS: usize = 2;

/// The segments of a peptide of length `len`: two peptides within
/// `NEAR_COPY_MAX_SUBSTITUTIONS` substitutions share one of them exactly.
fn segments(len: usize) -> Vec<(usize, usize)> {
    let n = NEAR_COPY_MAX_SUBSTITUTIONS + 1;
    (0..n).map(|k| (k * len / n, (k + 1) * len / n)).collect()
}

/// The bare sequence of a peptide as an andes table writes it (`K.PEPM[UNIMOD:35]K.L`),
/// with I read as L.
pub(crate) fn folded(peptide: &str) -> Vec<u8> {
    let core = match peptide.split('.').collect::<Vec<_>>().as_slice() {
        [_, mid, _] => *mid,
        _ => peptide,
    };
    let mut out = Vec::with_capacity(core.len());
    let mut depth = 0usize;
    for b in core.bytes() {
        match b {
            b'[' => depth += 1,
            b']' => depth = depth.saturating_sub(1),
            b'A'..=b'Z' if depth == 0 => out.push(if b == b'I' { b'L' } else { b }),
            _ => {}
        }
    }
    out
}

fn substitutions(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

/// For each of `peptides` (folded), the accessions of the proteins of `db` holding a
/// cleavage-bounded window within `NEAR_COPY_MAX_SUBSTITUTIONS` substitutions of it, at most
/// `cap` per peptide. A window starts at the protein N-terminus or after K/R and ends at the
/// C-terminus or on K/R (tryptic and Lys-C peptides; the cohort's cross-link searches).
pub(crate) fn near_copies(
    db: &Path,
    peptides: &[Vec<u8>],
    cap: usize,
) -> Result<Vec<Vec<String>>, String> {
    let mut index: HashMap<(usize, usize, &[u8]), Vec<u32>> = HashMap::new();
    let mut lengths: HashSet<usize> = HashSet::new();
    for (i, p) in peptides.iter().enumerate() {
        if p.len() < NEAR_COPY_MAX_SUBSTITUTIONS + 1 {
            continue;
        }
        lengths.insert(p.len());
        for (k, (a, b)) in segments(p.len()).into_iter().enumerate() {
            index
                .entry((p.len(), k, &p[a..b]))
                .or_default()
                .push(i as u32);
        }
    }
    let max_len = lengths.iter().copied().max().unwrap_or(0);
    let mut out: Vec<Vec<String>> = vec![Vec::new(); peptides.len()];
    let text =
        std::fs::read_to_string(db).map_err(|e| format!("cannot read {}: {e}", db.display()))?;
    let mut scan = |accession: &str, seq: &[u8]| {
        let seq: Vec<u8> = seq
            .iter()
            .filter(|b| b.is_ascii_alphabetic())
            .map(|b| match b.to_ascii_uppercase() {
                b'I' => b'L',
                c => c,
            })
            .collect();
        let starts = std::iter::once(0).chain(
            (0..seq.len())
                .filter(|&i| matches!(seq[i], b'K' | b'R'))
                .map(|i| i + 1),
        );
        for s in starts {
            for e in s + 1..=seq.len().min(s + max_len) {
                if !(e == seq.len() || matches!(seq[e - 1], b'K' | b'R')) {
                    continue;
                }
                let len = e - s;
                if !lengths.contains(&len) {
                    continue;
                }
                let w = &seq[s..e];
                for (k, (a, b)) in segments(len).into_iter().enumerate() {
                    let Some(hits) = index.get(&(len, k, &w[a..b])) else {
                        continue;
                    };
                    for &i in hits {
                        let o = &mut out[i as usize];
                        if o.len() < cap
                            && !o.iter().any(|x| x == accession)
                            && substitutions(&peptides[i as usize], w)
                                <= NEAR_COPY_MAX_SUBSTITUTIONS
                        {
                            o.push(accession.to_string());
                        }
                    }
                }
            }
        }
    };
    let mut accession = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !accession.is_empty() {
                scan(&accession, &seq);
            }
            accession = h.split_whitespace().next().unwrap_or("").to_string();
            seq.clear();
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !accession.is_empty() {
        scan(&accession, &seq);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_table_peptide_folds_to_its_bare_sequence_with_i_as_l() {
        assert_eq!(folded("K.PEPM[UNIMOD:35]IK.L"), b"PEPMLK");
        assert_eq!(folded("TIDEK"), b"TLDEK");
    }

    #[test]
    fn copies_within_two_substitutions_in_cleavage_bounded_windows_are_found() {
        let t = tempfile::tempdir().unwrap();
        let db = t.path().join("db.fasta");
        std::fs::write(
            &db,
            ">P1 one\nMKTLDEAGHSKR\n>P2 two\nAKTIDEAGNSKR\n>P3 three\nAKTWWWWWWSKR\n>P4 four\nATIDEAGHSKR\n",
        )
        .unwrap();
        let hits = near_copies(&db, &[b"TLDEAGHSK".to_vec()], 8).unwrap();
        assert_eq!(
            hits[0],
            ["P1", "P2"],
            "P2 is one substitution away (I = L), P3 too far, P4 not cleavage-bounded"
        );
        let capped = near_copies(&db, &[b"TLDEAGHSK".to_vec()], 1).unwrap();
        assert_eq!(capped[0].len(), 1);
    }
}
