//! Fragment-ion prediction for a Peptide.
//!
//! Phase 5 Task 2: canonical b/y ions only, no neutral losses. Produces
//! `(PredictedIon, m/z)` pairs at every requested charge.
//!
//! Phase 6 Task 4: adds `ions_for_node` for per-nominal-mass GF DP scoring.

use std::ops::RangeInclusive;

use model::amino_acid::AminoAcid;
use model::mass::{H2O, PROTON};
use crate::param_model::{IonType, Param};
use model::peptide::Peptide;

/// For a single prefix or suffix node at `nominal_mass`, enumerate the
/// `(ion_type, theo_mz)` pairs that contribute to its node score under
/// `param`. Java reference: `NewScoredSpectrum.getNodeScore(nodeMass, isPrefix)`.
///
/// `is_prefix = true` → walk prefix ions (b-ions etc.); `false` → suffix (y-ions etc.).
/// `parent_mass` / `charge` select the segment+partition used downstream.
///
/// Returns only the `(IonType, theo_mz)` pairs whose segment, when re-derived
/// from `theo_mz`, matches the segment from which the ion was collected.
pub fn ions_for_node(
    nominal_mass: f64,
    is_prefix: bool,
    param: &Param,
    parent_mass: f64,
    charge: u8,
) -> Vec<(IonType, f64)> {
    // Java parity: per-segment iteration uses the SPECIFIC partition's ion
    // type list (mirrors `NewScoredSpectrum`'s
    // `ionTypes[seg] = scorer.getIonTypes(charge, parentMass, seg)`). The
    // earlier `param.ion_types_for_segment(seg)` returned the union across
    // all partitions in the segment, enumerating extra ion types that Java
    // doesn't score for this spectrum. Each extra ion contributed a
    // missing-ion penalty (negative score), pulling Rust's PSM scores
    // systematically below Java's.
    // (Bug fix 2026-05-09: dominant cause of the per-PSM RawScore gap on
    // PXD001819 — Java HAEHIK = 80, Rust HAEHIK was 19; expected to close
    // most of that 60-point delta.)
    let mut out = Vec::new();
    let num_segs = param.num_segments as usize;
    for seg in 0..num_segs {
        for ion in param.ion_types_for_partition(charge, parent_mass, seg) {
            let theo_mz = match (is_prefix, ion) {
                (true, IonType::Prefix { .. }) => ion.mz(nominal_mass),
                (false, IonType::Suffix { .. }) => ion.mz(nominal_mass),
                _ => continue,
            };
            // Mirror Java: verify the ion's computed mz actually falls in this segment.
            if param.segment_num(theo_mz, parent_mass) != seg {
                continue;
            }
            out.push((ion, theo_mz));
        }
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IonKind {
    /// N-terminal fragment (b-ion). Neutral mass = sum of prefix residues.
    B,
    /// C-terminal fragment (y-ion). Neutral mass = sum of suffix residues + H2O.
    Y,
}

#[derive(Debug, Clone, Copy)]
pub struct PredictedIon {
    pub kind: IonKind,
    /// 1-based: b1 = prefix length 1, y1 = suffix length 1, etc.
    pub position: u32,
    pub charge: u8,
    /// Predicted m/z value.
    pub mz: f64,
}

/// Predict every canonical b/y ion at each charge in `charge_range`.
/// For a peptide of length n, produces `2*(n-1)*|charge_range|` ions:
/// b1..b_{n-1} and y1..y_{n-1} at each charge.
pub fn predict_by_ions(peptide: &Peptide, charge_range: RangeInclusive<u8>) -> Vec<PredictedIon> {
    let residues = &peptide.residues;
    let n = residues.len();
    if n < 2 || charge_range.is_empty() {
        return Vec::new();
    }

    // Cumulative residue masses (including any mods). Index i = sum of
    // residues[0..i]. cumulative[0] = 0; cumulative[n] = total residue mass.
    let mut cumulative: Vec<f64> = Vec::with_capacity(n + 1);
    cumulative.push(0.0);
    let mut acc = 0.0;
    for aa in residues {
        acc += residue_mass_with_mod(aa);
        cumulative.push(acc);
    }
    let total_residue_mass = cumulative[n];

    let mut out = Vec::with_capacity(
        2 * (n - 1) * (charge_range.end() - charge_range.start() + 1) as usize,
    );
    for charge in charge_range.clone() {
        let z = charge as f64;
        for k in 1..n {
            // b-ion at position k: neutral mass = sum of residues 0..k
            let b_neutral = cumulative[k];
            let b_mz = (b_neutral + z * PROTON) / z;
            out.push(PredictedIon {
                kind: IonKind::B,
                position: k as u32,
                charge,
                mz: b_mz,
            });

            // y-ion at position k: neutral mass = sum of residues n-k..n + H2O
            let y_neutral = total_residue_mass - cumulative[n - k] + H2O;
            let y_mz = (y_neutral + z * PROTON) / z;
            out.push(PredictedIon {
                kind: IonKind::Y,
                position: k as u32,
                charge,
                mz: y_mz,
            });
        }
    }
    out
}

fn residue_mass_with_mod(aa: &AminoAcid) -> f64 {
    aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn empty_charge_set_produces_no_ions() {
        let peptide = pep(b"PEPTIDE");
        // Build an empty RangeInclusive without triggering the reversed_empty_ranges lint.
        let empty: RangeInclusive<u8> = RangeInclusive::new(1, 0);
        let ions = predict_by_ions(&peptide, empty);
        assert!(ions.is_empty());
    }

    #[test]
    fn short_peptide_one_charge() {
        let peptide = pep(b"AR"); // 2 residues
        let ions = predict_by_ions(&peptide, 1..=1);
        // For a 2-residue peptide, prefix lengths are 1 only (b1).
        // Suffix lengths are 1 only (y1). 2 ions total at charge 1.
        assert_eq!(ions.len(), 2);
    }

    #[test]
    fn b_ion_mz_for_alanine_at_charge_1() {
        let peptide = pep(b"AR");
        let ions = predict_by_ions(&peptide, 1..=1);
        // b1 is the A residue alone. A residue mass = 71.0371...
        // m/z = (71.0371 + 1 * PROTON) / 1 = 72.0444...
        let a_mass = AminoAcid::standard(b'A').unwrap().mass;
        let expected_b1 = (a_mass + PROTON) / 1.0;
        let b1 = ions
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 1 && p.charge == 1)
            .expect("b1+1");
        assert!(
            (b1.mz - expected_b1).abs() < 1e-9,
            "b1+1 mz drift: got {}, expected {}",
            b1.mz,
            expected_b1
        );
    }

    #[test]
    fn y_ion_mz_for_arginine_at_charge_1() {
        let peptide = pep(b"AR");
        let ions = predict_by_ions(&peptide, 1..=1);
        // y1 is the R residue + H2O. R residue mass = 156.1011...
        // y1 neutral mass = R + H2O.
        // m/z = (R + H2O + 1 * PROTON) / 1
        let r_mass = AminoAcid::standard(b'R').unwrap().mass;
        let expected_y1 = (r_mass + H2O + PROTON) / 1.0;
        let y1 = ions
            .iter()
            .find(|p| matches!(p.kind, IonKind::Y) && p.position == 1 && p.charge == 1)
            .expect("y1+1");
        assert!(
            (y1.mz - expected_y1).abs() < 1e-9,
            "y1+1 mz drift: got {}, expected {}",
            y1.mz,
            expected_y1
        );
    }

    #[test]
    fn ion_count_scales_with_peptide_length() {
        // Length-3 peptide → b1, b2 (2 b-ions) + y1, y2 (2 y-ions) = 4 ions per charge.
        let peptide = pep(b"AGR");
        let ions = predict_by_ions(&peptide, 1..=1);
        assert_eq!(ions.len(), 4);

        // Length-5 peptide → 4 b + 4 y = 8 ions per charge.
        let peptide = pep(b"PEPTR");
        let ions = predict_by_ions(&peptide, 1..=1);
        assert_eq!(ions.len(), 8);
    }

    #[test]
    fn multi_charge_doubles_ion_count() {
        let peptide = pep(b"AGR");
        let ions_1 = predict_by_ions(&peptide, 1..=1);
        let ions_12 = predict_by_ions(&peptide, 1..=2);
        assert_eq!(ions_12.len(), ions_1.len() * 2);
    }

    #[test]
    fn charge_2_mz_is_about_half_of_charge_1() {
        let peptide = pep(b"PEPTIDER");
        let ions = predict_by_ions(&peptide, 1..=2);
        // Same b/y position at charge 2 should be roughly half + small shift due to proton mass.
        let b3_z1 = ions
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 3 && p.charge == 1)
            .unwrap();
        let b3_z2 = ions
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 3 && p.charge == 2)
            .unwrap();
        // m/z2 = (neutral + 2*PROTON) / 2 vs m/z1 = (neutral + PROTON) / 1
        // m/z2 - m/z1/2 = PROTON/2 - PROTON/2 = 0... actually
        // m/z2 = neutral/2 + PROTON
        // m/z1/2 = neutral/2 + PROTON/2
        // So m/z2 = m/z1/2 + PROTON/2
        assert!((b3_z2.mz - (b3_z1.mz / 2.0 + PROTON / 2.0)).abs() < 1e-9);
    }
}
