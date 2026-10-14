//! Helpers for andes PIN peptide strings (`K.M+15.99491PEPTIDEK.R`).
//!
//! andes prints at most one mass delta per residue (`Peptide`'s `Display`).
//! When a terminal modification and a residue modification sit on the same
//! residue, the printed delta is their sum: fixed carbamidomethyl-C under a
//! variable protein N-terminal acetyl is written `C+99.03203`. N-terminal
//! deltas are printed after the first residue.

/// Split a delta printed after the first residue into an N-terminal and a
/// residue modification whose masses sum to it.
///
/// `nterm` and `residue` are the `(key, monoisotopic mass)` candidates allowed
/// at the N-terminus and on the first residue. Returns the pair with the
/// smallest mass error within `tol` Da; an exact tie keeps the first pair in
/// candidate order (N-terminal candidates outer, residue candidates inner), so
/// callers that want a deterministic tie-break pass candidates in that order.
pub fn split_nterm_residue_delta<K: Copy>(
    delta: f64,
    tol: f64,
    nterm: &[(K, f64)],
    residue: &[(K, f64)],
) -> Option<(K, K)> {
    let mut best: Option<(f64, K, K)> = None;
    for &(n, n_mass) in nterm {
        for &(r, r_mass) in residue {
            let err = (n_mass + r_mass - delta).abs();
            if err <= tol && best.is_none_or(|(e, _, _)| err < e) {
                best = Some((err, n, r));
            }
        }
    }
    best.map(|(_, n, r)| (n, r))
}

#[cfg(test)]
mod tests {
    use super::*;

    const NTERM: &[(&str, f64)] = &[("Acetyl", 42.010565), ("TMT6plex", 229.162932)];
    const RESIDUE_C: &[(&str, f64)] = &[("Carbamidomethyl", 57.021464)];

    #[test]
    fn splits_cam_plus_acetyl() {
        assert_eq!(
            split_nterm_residue_delta(99.03203, 2e-3, NTERM, RESIDUE_C),
            Some(("Acetyl", "Carbamidomethyl"))
        );
    }

    #[test]
    fn rejects_sums_outside_tolerance_or_without_candidates() {
        assert_eq!(
            split_nterm_residue_delta(99.04, 2e-3, NTERM, RESIDUE_C),
            None
        );
        assert_eq!(
            split_nterm_residue_delta(99.03203, 2e-3, &[], RESIDUE_C),
            None
        );
        assert_eq!(
            split_nterm_residue_delta::<&str>(99.03203, 2e-3, NTERM, &[]),
            None
        );
    }

    #[test]
    fn closest_sum_wins_and_exact_ties_keep_candidate_order() {
        let nterm = [(1, 10.0), (2, 10.0015)];
        let residue = [(7, 5.0)];
        assert_eq!(
            split_nterm_residue_delta(15.0014, 2e-3, &nterm, &residue),
            Some((2, 7))
        );
        assert_eq!(
            split_nterm_residue_delta(15.0001, 2e-3, &nterm, &residue),
            Some((1, 7))
        );
        let tied = [(3, 10.0), (4, 10.0)];
        assert_eq!(
            split_nterm_residue_delta(15.0, 2e-3, &tied, &residue),
            Some((3, 7))
        );
    }
}
