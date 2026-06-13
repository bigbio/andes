//! Fragment-ion prediction for a Peptide.
//!
//! Canonical b/y ions, plus optional activation-gated, per-class
//! neutral-loss-shifted b/y ions for residues whose modification declares
//! `neutral_losses`. Produces `PredictedIon`s (each tagged with a
//! `loss_class`; 0 = intact) at every requested charge. Also exposes
//! `ions_for_node` for per-nominal-mass GF DP scoring.
//!
//! ## Neutral-loss emission (gated, inert by default)
//!
//! [`predict_by_ions`] never emits loss ions: it is a thin wrapper over
//! [`predict_by_ions_with_losses`] with the loss gate `false`, so every
//! existing caller is byte-identical to the canonical-only behaviour.
//!
//! [`predict_by_ions_with_losses`] emits a loss-shifted partner for each
//! intact b/y ion whose residue span includes a residue carrying
//! `neutral_losses`, ONLY when `predict_losses` is `true` (callers pass
//! `scorer.param().data_type.activation.predicts_neutral_losses()` — ETD
//! preserves labile mods, so it is `false`). Each declared loss `L` for a
//! spanned residue produces one ion at `mz_intact − L/z` with the same
//! series/charge/position and `loss_class` set to that mod's class. The
//! intact ion is always kept. If no spanned residue declares losses, or the
//! gate is off, ZERO loss ions are emitted ⇒ output identical to the
//! canonical b/y set. v1: multiple loss-bearing residues in one fragment
//! emit each residue's losses independently (no cross-products).

use std::ops::RangeInclusive;

use model::amino_acid::AminoAcid;
use model::mass::{H2O, PROTON};
use crate::param_model::{IonType, Param};
use model::peptide::Peptide;

/// For a single prefix or suffix node at `nominal_mass`, enumerate the
/// `(ion_type, theo_mz)` pairs that contribute to its node score under `param`.
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
    // Compat shim — callers in hot paths should use `for_each_ion_for_node`
    // to avoid the per-call Vec allocation.
    let mut out = Vec::new();
    for_each_ion_for_node(nominal_mass, is_prefix, param, parent_mass, charge, |ion, theo_mz, _part| {
        out.push((ion, theo_mz));
    });
    out
}

/// Callback variant of `ions_for_node`. Calls `f(ion, theo_mz, partition)`
/// once per (ion, theo_mz) pair without allocating an intermediate Vec.
/// Used by `directional_node_score` in the GF DP hot path (~5 splits ×
/// 2 directions × ~38k spectra ÷ 12 threads = millions of calls per search).
///
/// `partition` is precomputed per outer-segment iteration (constant for
/// all ions in that segment). Saves a `partition_for` binary search per
/// ion (was ~30 ns × millions of calls).
///
/// See `ions_for_node` for the per-segment / per-partition iteration
/// semantics. Produces the same set of (ion, theo_mz) pairs in the same order.
#[inline]
pub fn for_each_ion_for_node<F: FnMut(IonType, f64, crate::param_model::Partition)>(
    nominal_mass: f64,
    is_prefix: bool,
    param: &Param,
    parent_mass: f64,
    charge: u8,
    mut f: F,
) {
    let num_segs = param.num_segments as usize;
    for seg in 0..num_segs {
        // Partition is constant for all ions in this segment.
        let partition = param.partition_for(charge, parent_mass, seg);
        for &ion in param.ion_types_for_partition_slice(charge, parent_mass, seg) {
            let theo_mz = match (is_prefix, ion) {
                (true, IonType::Prefix { .. }) => ion.mz(nominal_mass),
                (false, IonType::Suffix { .. }) => ion.mz(nominal_mass),
                _ => continue,
            };
            // Verify the ion's computed mz actually falls in this segment.
            if param.segment_num(theo_mz, parent_mass) != seg {
                continue;
            }
            f(ion, theo_mz, partition);
        }
    }
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
    /// Neutral-loss class of this ion: 0 = intact (no neutral loss); 1.. = a
    /// per-mod-class loss pool (glyco=1, phospho=2, sulfo=3, generic=255),
    /// matching `Modification::loss_class`. Loss ions sit at `mz_intact − L/z`.
    pub loss_class: u8,
}

impl PredictedIon {
    /// Loss-class id: 0 = intact; 1.. = a per-mod-class neutral-loss pool.
    #[inline]
    pub fn loss_class(&self) -> u8 {
        self.loss_class
    }

    /// True if this is a neutral-loss-shifted fragment ion (any loss class).
    #[inline]
    pub fn is_loss(&self) -> bool {
        self.loss_class != 0
    }
}

/// Predict every canonical b/y ion at each charge in `charge_range`.
/// For a peptide of length n, produces `2*(n-1)*|charge_range|` ions:
/// b1..b_{n-1} and y1..y_{n-1} at each charge.
///
/// Never emits neutral-loss ions: this is [`predict_by_ions_with_losses`]
/// with the loss gate `false`, so its output is byte-identical to the
/// canonical-only behaviour. Callers on the scoring path that want
/// activation-gated loss ions call [`predict_by_ions_with_losses`] directly.
#[inline]
pub fn predict_by_ions(peptide: &Peptide, charge_range: RangeInclusive<u8>) -> Vec<PredictedIon> {
    predict_by_ions_with_losses(peptide, charge_range, false)
}

/// Predict canonical b/y ions at each charge in `charge_range`, plus optional
/// neutral-loss-shifted partners.
///
/// When `predict_losses` is `true`, for each intact b/y ion whose residue
/// span includes a residue whose modification declares `neutral_losses`, one
/// loss-shifted ion is additionally emitted per declared loss `L` at
/// `mz_intact − L/z`, with the same series/charge/position and `loss_class`
/// set to that mod's class. The intact ion is always kept.
///
/// When `predict_losses` is `false`, OR when no spanned residue declares
/// losses, ZERO loss ions are emitted and the output is byte-identical to the
/// canonical b/y set (this is the safety guarantee behind [`predict_by_ions`]).
///
/// v1: multiple loss-bearing residues within a single fragment emit each
/// residue's losses independently (no cross-products of simultaneous losses).
pub fn predict_by_ions_with_losses(
    peptide: &Peptide,
    charge_range: RangeInclusive<u8>,
    predict_losses: bool,
) -> Vec<PredictedIon> {
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

    // Per-residue declared losses, only when the gate is on. `loss_residues`
    // is empty (and never consulted) when `predict_losses` is false or no
    // residue carries `neutral_losses` ⇒ the loss-emission branches below are
    // inert and the output is byte-identical to the canonical set.
    //
    // Entry per loss-bearing residue: (residue_index, loss_class, &[loss_da]).
    let loss_residues: Vec<(usize, u8, &[f64])> = if predict_losses {
        residues
            .iter()
            .enumerate()
            .filter_map(|(i, aa)| {
                aa.mod_.as_ref().and_then(|m| {
                    (!m.neutral_losses.is_empty())
                        .then(|| (i, m.loss_class, m.neutral_losses.as_slice()))
                })
            })
            .collect()
    } else {
        Vec::new()
    };

    let mut out = Vec::with_capacity(
        2 * (n - 1) * (charge_range.end() - charge_range.start() + 1) as usize,
    );
    for charge in charge_range.clone() {
        let z = charge as f64;
        for k in 1..n {
            // b-ion at position k: neutral mass = sum of residues 0..k.
            // The prefix fragment spans residue indices [0, k).
            let b_neutral = cumulative[k];
            let b_mz = (b_neutral + z * PROTON) / z;
            out.push(PredictedIon {
                kind: IonKind::B,
                position: k as u32,
                charge,
                mz: b_mz,
                loss_class: 0,
            });
            emit_loss_ions(&mut out, &loss_residues, 0..k, IonKind::B, k as u32, charge, z, b_mz);

            // y-ion at position k: neutral mass = sum of residues n-k..n + H2O.
            // The suffix fragment spans residue indices [n-k, n).
            let y_neutral = total_residue_mass - cumulative[n - k] + H2O;
            let y_mz = (y_neutral + z * PROTON) / z;
            out.push(PredictedIon {
                kind: IonKind::Y,
                position: k as u32,
                charge,
                mz: y_mz,
                loss_class: 0,
            });
            emit_loss_ions(&mut out, &loss_residues, (n - k)..n, IonKind::Y, k as u32, charge, z, y_mz);
        }
    }
    out
}

/// Emit loss-shifted partners for an intact ion at `intact_mz`.
///
/// For each loss-bearing residue whose index lies in `span`, push one ion per
/// declared loss `L` at `intact_mz − L/z`, tagged with that residue's
/// `loss_class`. `loss_residues` is empty whenever loss prediction is
/// disabled or no residue declares losses, so this is a no-op (the loop body
/// never runs) on the standard path ⇒ output unchanged.
#[inline]
#[allow(clippy::too_many_arguments)]
fn emit_loss_ions(
    out: &mut Vec<PredictedIon>,
    loss_residues: &[(usize, u8, &[f64])],
    span: std::ops::Range<usize>,
    kind: IonKind,
    position: u32,
    charge: u8,
    z: f64,
    intact_mz: f64,
) {
    for &(idx, loss_class, losses) in loss_residues {
        if !span.contains(&idx) {
            continue;
        }
        for &loss in losses {
            out.push(PredictedIon {
                kind,
                position,
                charge,
                mz: intact_mz - loss / z,
                loss_class,
            });
        }
    }
}

fn residue_mass_with_mod(aa: &AminoAcid) -> f64 {
    aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use std::sync::Arc;

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    fn plain_peptide() -> Peptide {
        pep(b"PEPTIDE")
    }

    /// `PEPTIDE` with a loss-bearing mod (the declared losses + class) on the
    /// middle T residue (index 3), so both b and y series span it.
    fn peptide_with_loss_mod(losses: Vec<f64>, loss_class: u8) -> Peptide {
        let m = Modification {
            name: "TestLoss".to_string(),
            mass_delta: 0.0,
            residue: ResidueSpec::Specific(b'T'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: losses,
            loss_class,
        };
        let arc = Arc::new(m);
        let residues: Vec<AminoAcid> = b"PEPTIDE"
            .iter()
            .map(|&r| {
                let aa = AminoAcid::standard(r).unwrap();
                if r == b'T' {
                    aa.with_mod(arc.clone())
                } else {
                    aa
                }
            })
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn emits_loss_ions_for_loss_bearing_residue() {
        let pep = peptide_with_loss_mod(vec![162.0528, 324.1056], 1);
        let ions = predict_by_ions_with_losses(&pep, 1..=1, true);
        // intact present
        assert!(ions.iter().any(|p| !p.is_loss()));
        let losses: Vec<_> = ions.iter().filter(|p| p.is_loss()).collect();
        assert!(!losses.is_empty());
        // all loss ions tagged glyco (class 1)
        assert!(losses.iter().all(|p| p.loss_class() == 1));
        // a loss ion sits 162.0528 below some intact ion of the same series/charge:
        assert!(losses.iter().any(|l| ions.iter().any(|i| {
            !i.is_loss()
                && i.kind == l.kind
                && i.charge == l.charge
                && (i.mz - l.mz - 162.0528).abs() < 1e-4
        })));
        // the second declared loss is also represented.
        assert!(losses.iter().any(|l| ions.iter().any(|i| {
            !i.is_loss()
                && i.kind == l.kind
                && i.charge == l.charge
                && (i.mz - l.mz - 324.1056).abs() < 1e-4
        })));
    }

    #[test]
    fn no_loss_ions_when_disabled_or_no_loss_mod() {
        let plain = plain_peptide();
        // a plain peptide never produces loss ions, even with the gate on.
        assert!(predict_by_ions_with_losses(&plain, 1..=1, true)
            .iter()
            .all(|p| !p.is_loss()));
        let pep = peptide_with_loss_mod(vec![162.0528], 1);
        // ETD / gate disabled ⇒ no loss ions even with a loss mod:
        assert!(predict_by_ions_with_losses(&pep, 1..=1, false)
            .iter()
            .all(|p| !p.is_loss()));
    }

    #[test]
    fn predict_by_ions_is_byte_identical_to_gate_off() {
        // The thin wrapper must equal the explicit gate-off call, including on
        // a loss-bearing peptide (the inertness guarantee).
        let pep = peptide_with_loss_mod(vec![162.0528, 324.1056], 1);
        let via_wrapper = predict_by_ions(&pep, 1..=2);
        let via_gate_off = predict_by_ions_with_losses(&pep, 1..=2, false);
        assert_eq!(via_wrapper.len(), via_gate_off.len());
        for (a, b) in via_wrapper.iter().zip(via_gate_off.iter()) {
            assert_eq!(a.kind, b.kind);
            assert_eq!(a.position, b.position);
            assert_eq!(a.charge, b.charge);
            assert_eq!(a.mz.to_bits(), b.mz.to_bits());
            assert_eq!(a.loss_class, b.loss_class);
        }
        assert!(via_wrapper.iter().all(|p| !p.is_loss()));
    }

    #[test]
    fn activation_etd_does_not_predict_losses() {
        use model::activation::ActivationMethod;
        assert!(!ActivationMethod::ETD.predicts_neutral_losses());
        assert!(ActivationMethod::HCD.predicts_neutral_losses());
        assert!(ActivationMethod::CID.predicts_neutral_losses());
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
