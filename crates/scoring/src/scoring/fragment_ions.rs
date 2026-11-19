use std::ops::RangeInclusive;

use crate::param_model::{IonType, Param};
use model::amino_acid::AminoAcid;
use model::mass::{H2O, PROTON};
use model::peptide::Peptide;

pub fn ions_for_node(
    nominal_mass: f64,
    is_prefix: bool,
    param: &Param,
    parent_mass: f64,
    charge: u8,
) -> Vec<(IonType, f64)> {
    let mut out = Vec::new();
    for_each_ion_for_node(
        nominal_mass,
        is_prefix,
        param,
        parent_mass,
        charge,
        |ion, theo_mz, _part| {
            out.push((ion, theo_mz));
        },
    );
    out
}

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
        let partition = param.partition_for(charge, parent_mass, seg);
        for &ion in param.ion_types_for_partition_slice(charge, parent_mass, seg) {
            let theo_mz = match (is_prefix, ion) {
                (true, IonType::Prefix { .. }) => ion.mz(nominal_mass),
                (false, IonType::Suffix { .. }) => ion.mz(nominal_mass),
                _ => continue,
            };
            if param.segment_num(theo_mz, parent_mass) != seg {
                continue;
            }
            f(ion, theo_mz, partition);
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IonKind {
    B,
    Y,
    C,
    Z,
}

#[derive(Debug, Clone, Copy)]
pub struct PredictedIon {
    pub kind: IonKind,
    pub position: u32,
    pub charge: u8,
    pub mz: f64,
    pub loss_class: u8,
}

impl PredictedIon {
    #[inline]
    pub fn loss_class(&self) -> u8 {
        self.loss_class
    }

    #[inline]
    pub fn is_loss(&self) -> bool {
        self.loss_class != 0
    }
}

#[inline]
pub fn predict_by_ions(peptide: &Peptide, charge_range: RangeInclusive<u8>) -> Vec<PredictedIon> {
    predict_by_ions_with_losses(peptide, charge_range, false)
}

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

    let mut cumulative: Vec<f64> = Vec::with_capacity(n + 1);
    cumulative.push(0.0);
    let mut acc = 0.0;
    for aa in residues {
        acc += residue_mass_with_mod(aa);
        cumulative.push(acc);
    }
    let total_residue_mass = cumulative[n];

    let loss_residues: Vec<(usize, u8, &[f64])> = if predict_losses {
        residues
            .iter()
            .enumerate()
            .filter_map(|(i, aa)| {
                aa.mod_.as_ref().and_then(|m| {
                    (!m.neutral_losses.is_empty()).then_some((
                        i,
                        m.loss_class,
                        m.neutral_losses.as_slice(),
                    ))
                })
            })
            .collect()
    } else {
        Vec::new()
    };

    let mut out =
        Vec::with_capacity(2 * (n - 1) * (charge_range.end() - charge_range.start() + 1) as usize);
    for charge in charge_range.clone() {
        let z = charge as f64;
        for k in 1..n {
            let b_neutral = cumulative[k];
            let b_mz = (b_neutral + z * PROTON) / z;
            out.push(PredictedIon {
                kind: IonKind::B,
                position: k as u32,
                charge,
                mz: b_mz,
                loss_class: 0,
            });
            emit_loss_ions(
                &mut out,
                &loss_residues,
                0..k,
                IonKind::B,
                k as u32,
                charge,
                z,
                b_mz,
            );

            let y_neutral = total_residue_mass - cumulative[n - k] + H2O;
            let y_mz = (y_neutral + z * PROTON) / z;
            out.push(PredictedIon {
                kind: IonKind::Y,
                position: k as u32,
                charge,
                mz: y_mz,
                loss_class: 0,
            });
            emit_loss_ions(
                &mut out,
                &loss_residues,
                (n - k)..n,
                IonKind::Y,
                k as u32,
                charge,
                z,
                y_mz,
            );
        }
    }
    out
}

pub const NH3: f64 = model::mass::NH3;
pub const Z_DOT_OFFSET: f64 = model::mass::Z_DOT_OFFSET;

fn cz_remnant_enabled() -> bool {
    use std::sync::OnceLock;
    static F: OnceLock<bool> = OnceLock::new();
    *F.get_or_init(|| false)
}

fn cz_spanning_shifts(glycan_mass: f64, remnant: bool) -> Vec<f64> {
    if glycan_mass == 0.0 || !remnant {
        return vec![glycan_mass];
    }
    use model::glycan_residue::{HEX, HEXNAC};
    const REMNANTS: [f64; 5] = [
        0.0,
        HEXNAC,
        HEXNAC + HEX,
        2.0 * HEXNAC,
        2.0 * HEXNAC + 3.0 * HEX,
    ];

    let mut v = Vec::with_capacity(1 + REMNANTS.len());
    v.push(glycan_mass);
    for &r in &REMNANTS {
        if (r - glycan_mass).abs() > 1e-6 {
            v.push(r);
        }
    }
    v
}

pub fn predict_cz_ions(
    peptide: &Peptide,
    charge_range: RangeInclusive<u8>,
    glycan_mass: f64,
    glycosite: usize,
) -> Vec<PredictedIon> {
    let residues = &peptide.residues;
    let n = residues.len();
    if n < 2 || charge_range.is_empty() {
        return Vec::new();
    }
    let mut cumulative: Vec<f64> = Vec::with_capacity(n + 1);
    cumulative.push(0.0);
    let mut acc = 0.0;
    for aa in residues {
        acc += residue_mass_with_mod(aa);
        cumulative.push(acc);
    }
    let total_residue_mass = cumulative[n];
    let span_shifts = cz_spanning_shifts(glycan_mass, cz_remnant_enabled());
    let mut out = Vec::with_capacity(
        2 * (n - 1) * (charge_range.end() - charge_range.start() + 1) as usize * span_shifts.len(),
    );
    for charge in charge_range.clone() {
        let z = charge as f64;
        for k in 1..n {
            let c_prefix = cumulative[k] + NH3;
            if glycosite < k {
                for &shift in &span_shifts {
                    out.push(PredictedIon {
                        kind: IonKind::C,
                        position: k as u32,
                        charge,
                        mz: (c_prefix + shift + z * PROTON) / z,
                        loss_class: 0,
                    });
                }
            } else {
                out.push(PredictedIon {
                    kind: IonKind::C,
                    position: k as u32,
                    charge,
                    mz: (c_prefix + z * PROTON) / z,
                    loss_class: 0,
                });
            }
            let z_suffix = (total_residue_mass - cumulative[n - k] + H2O) + Z_DOT_OFFSET;
            if glycosite >= n - k {
                for &shift in &span_shifts {
                    out.push(PredictedIon {
                        kind: IonKind::Z,
                        position: k as u32,
                        charge,
                        mz: (z_suffix + shift + z * PROTON) / z,
                        loss_class: 0,
                    });
                }
            } else {
                out.push(PredictedIon {
                    kind: IonKind::Z,
                    position: k as u32,
                    charge,
                    mz: (z_suffix + z * PROTON) / z,
                    loss_class: 0,
                });
            }
        }
    }
    out
}

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
    for (j, &(idx, loss_class, losses)) in loss_residues.iter().enumerate() {
        if !span.contains(&idx) {
            continue;
        }
        for (k, &loss) in losses.iter().enumerate() {
            let same = |&l: &f64| l.to_bits() == loss.to_bits();
            let earlier = losses[..k].iter().any(same)
                || loss_residues[..j]
                    .iter()
                    .any(|&(i, c, ls)| span.contains(&i) && c == loss_class && ls.iter().any(same));
            if earlier {
                continue;
            }
            let mz = intact_mz - loss / z;
            if mz <= 0.0 {
                continue;
            }
            out.push(PredictedIon {
                kind,
                position,
                charge,
                mz,
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

    #[test]
    fn predict_cz_ions_offsets_and_glycan_spanning() {
        let p = pep(b"PEPTIDE");
        let by = predict_by_ions(&p, 1..=1);
        let cz = predict_cz_ions(&p, 1..=1, 0.0, 0);
        let b_at = |k: u32| {
            by.iter()
                .find(|i| i.kind == IonKind::B && i.position == k)
                .unwrap()
                .mz
        };
        let y_at = |k: u32| {
            by.iter()
                .find(|i| i.kind == IonKind::Y && i.position == k)
                .unwrap()
                .mz
        };
        let c_at = |k: u32| {
            cz.iter()
                .find(|i| i.kind == IonKind::C && i.position == k)
                .unwrap()
                .mz
        };
        let z_at = |k: u32| {
            cz.iter()
                .find(|i| i.kind == IonKind::Z && i.position == k)
                .unwrap()
                .mz
        };
        for k in 1..=6 {
            assert!((c_at(k) - (b_at(k) + NH3)).abs() < 1e-6, "c{k} != b{k}+NH3");
            assert!(
                (z_at(k) - (y_at(k) + Z_DOT_OFFSET)).abs() < 1e-6,
                "z{k} != y{k}+offset"
            );
        }
        let g = 1000.0;
        let czg = predict_cz_ions(&p, 1..=1, g, 2);
        let cg = |k: u32| {
            czg.iter()
                .find(|i| i.kind == IonKind::C && i.position == k)
                .unwrap()
                .mz
        };
        assert!(
            (cg(2) - c_at(2)).abs() < 1e-6,
            "c2 must be naked (prefix [0,2) excludes site 2)"
        );
        assert!(
            (cg(3) - (c_at(3) + g)).abs() < 1e-6,
            "c3 must carry glycan (prefix includes site 2)"
        );
    }

    #[test]
    fn cz_spanning_shifts_remnant_ladder() {
        assert_eq!(cz_spanning_shifts(1000.0, false), vec![1000.0]);
        assert_eq!(cz_spanning_shifts(0.0, true), vec![0.0]);
        let s = cz_spanning_shifts(1000.0, true);
        assert_eq!(s[0], 1000.0, "intact glycan must be first");
        assert_eq!(s.len(), 6, "intact glycan + 5 remnants, no collision");
        use model::glycan_residue::{HEX, HEXNAC};
        for r in [
            0.0,
            HEXNAC,
            HEXNAC + HEX,
            2.0 * HEXNAC,
            2.0 * HEXNAC + 3.0 * HEX,
        ] {
            assert!(
                s.iter().any(|&x| (x - r).abs() < 1e-9),
                "missing remnant {r}"
            );
        }
        let s2 = cz_spanning_shifts(HEXNAC, true);
        assert_eq!(
            s2.iter().filter(|&&x| (x - HEXNAC).abs() < 1e-6).count(),
            1,
            "remnant equal to the intact glycan must not be duplicated"
        );
    }

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
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let arc = model::modification::leak_mod(m);
        let residues: Vec<AminoAcid> = b"PEPTIDE"
            .iter()
            .map(|&r| {
                let aa = AminoAcid::standard(r).unwrap();
                if r == b'T' {
                    aa.with_mod(arc)
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
        assert!(ions.iter().any(|p| !p.is_loss()));
        let losses: Vec<_> = ions.iter().filter(|p| p.is_loss()).collect();
        assert!(!losses.is_empty());
        assert!(losses.iter().all(|p| p.loss_class() == 1));
        assert!(losses.iter().any(|l| ions.iter().any(|i| {
            !i.is_loss()
                && i.kind == l.kind
                && i.charge == l.charge
                && (i.mz - l.mz - 162.0528).abs() < 1e-4
        })));
        assert!(losses.iter().any(|l| ions.iter().any(|i| {
            !i.is_loss()
                && i.kind == l.kind
                && i.charge == l.charge
                && (i.mz - l.mz - 324.1056).abs() < 1e-4
        })));
    }

    #[test]
    fn shared_loss_of_two_sites_is_one_ion() {
        let m = model::modification::leak_mod(Modification {
            name: "Phospho".to_string(),
            mass_delta: 79.96633,
            residue: ResidueSpec::Specific(b'S'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: vec![97.97690],
            loss_class: 2,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        });
        let residues: Vec<AminoAcid> = b"SSAK"
            .iter()
            .map(|&r| {
                let aa = AminoAcid::standard(r).unwrap();
                if r == b'S' {
                    aa.with_mod(m)
                } else {
                    aa
                }
            })
            .collect();
        let pep = Peptide::new(residues, b'_', b'-');
        let ions = predict_by_ions_with_losses(&pep, 1..=1, true);
        for (kind, position) in [(IonKind::B, 2), (IonKind::Y, 3)] {
            let n = ions
                .iter()
                .filter(|p| p.is_loss() && p.kind == kind && p.position == position)
                .count();
            assert_eq!(n, 1, "{kind:?}{position}");
        }
    }

    #[test]
    fn no_loss_ions_when_disabled_or_no_loss_mod() {
        let plain = plain_peptide();
        assert!(predict_by_ions_with_losses(&plain, 1..=1, true)
            .iter()
            .all(|p| !p.is_loss()));
        let pep = peptide_with_loss_mod(vec![162.0528], 1);
        assert!(predict_by_ions_with_losses(&pep, 1..=1, false)
            .iter()
            .all(|p| !p.is_loss()));
    }

    #[test]
    fn predict_by_ions_is_byte_identical_to_gate_off() {
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
        let empty: RangeInclusive<u8> = RangeInclusive::new(1, 0);
        let ions = predict_by_ions(&peptide, empty);
        assert!(ions.is_empty());
    }

    #[test]
    fn short_peptide_one_charge() {
        let peptide = pep(b"AR");
        let ions = predict_by_ions(&peptide, 1..=1);
        assert_eq!(ions.len(), 2);
    }

    #[test]
    fn b_ion_mz_for_alanine_at_charge_1() {
        let peptide = pep(b"AR");
        let ions = predict_by_ions(&peptide, 1..=1);
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
        let peptide = pep(b"AGR");
        let ions = predict_by_ions(&peptide, 1..=1);
        assert_eq!(ions.len(), 4);

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
        let b3_z1 = ions
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 3 && p.charge == 1)
            .unwrap();
        let b3_z2 = ions
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 3 && p.charge == 2)
            .unwrap();
        assert!((b3_z2.mz - (b3_z1.mz / 2.0 + PROTON / 2.0)).abs() < 1e-9);
    }
}
