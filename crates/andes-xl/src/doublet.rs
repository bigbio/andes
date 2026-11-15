//! Reporter-doublet detection for MS-cleavable linkers.
//!
//! When an MS-cleavable linker such as DSSO breaks in MS2, each peptide is
//! released twice: once with the alkene stub and once with the thiol stub.
//! The two ions sit `delta / z` apart in m/z at fragment charge `z`
//! (`delta` = thiol − alkene = 31.9721 Da for DSSO). A doublet therefore
//! measures one peptide's neutral mass:
//!
//! ```text
//! M(peptide) = z · (mz_light − proton) − alkene
//! ```
//!
//! The detector works on the raw centroided peaks (no deconvolution): every
//! peak above an intensity floor is tried as the light member at each charge
//! `1..=max_charge` (capped at precursor charge − 1, since the other peptide
//! carries at least one charge), and the heavy member is looked up within a
//! ppm window. A doublet whose two members both have a peak one isotope
//! spacing below them is an isotope copy of a lower doublet and is dropped.
//! Doublets of the same peptide at several charges collapse into one
//! [`DoubletMass`].

use crate::crosslinker::{Cleavage, ISOTOPE_SPACING, PROTON};

/// Settings of [`find_doublets`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubletParams {
    /// m/z tolerance of the heavy member and of the isotope checks, in ppm.
    pub tol_ppm: f64,
    /// Highest fragment charge tried for a doublet.
    pub max_charge: u8,
    /// Intensity floor of both members, as a fraction of the base peak.
    pub min_rel_intensity: f64,
}

impl Default for DoubletParams {
    fn default() -> Self {
        DoubletParams {
            tol_ppm: 10.0,
            max_charge: 4,
            min_rel_intensity: 0.0,
        }
    }
}

/// One light/heavy peak pair at one charge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Doublet {
    /// Index of the light (alkene) member in the peak list.
    pub light_idx: usize,
    /// Index of the heavy (thiol) member in the peak list.
    pub heavy_idx: usize,
    /// Fragment charge.
    pub charge: u8,
    /// Neutral mass of the peptide, from the light member.
    pub peptide_mass: f64,
    /// Heavy-member m/z error against its expected position, in ppm.
    pub error_ppm: f64,
    /// Summed intensity of the two members.
    pub intensity: f64,
    /// Whether either member has its +1 isotope peak at this charge.
    pub isotope_confirmed: bool,
}

/// The doublets of one peptide mass, merged over charges.
#[derive(Debug, Clone, PartialEq)]
pub struct DoubletMass {
    /// Neutral peptide mass, from the most intense doublet of the group.
    pub mass: f64,
    /// Summed intensity of the group's doublets.
    pub intensity: f64,
    /// Number of charges at which the doublet was seen.
    pub n_charges: u8,
    /// Whether any member doublet is isotope-confirmed.
    pub isotope_confirmed: bool,
    /// Lowest charge in the group.
    pub min_charge: u8,
}

#[inline]
fn ppm_window(mz: f64, ppm: f64) -> f64 {
    mz * ppm * 1e-6
}

/// Index of the most intense peak within `tol` of `mz` in m/z-sorted `peaks`.
pub fn find_peak(peaks: &[(f64, f32)], mz: f64, tol: f64) -> Option<usize> {
    let lo = peaks.partition_point(|p| p.0 < mz - tol);
    let mut best: Option<usize> = None;
    for (i, p) in peaks.iter().enumerate().skip(lo) {
        if p.0 > mz + tol {
            break;
        }
        if best.is_none_or(|b| p.1 > peaks[b].1) {
            best = Some(i);
        }
    }
    best
}

/// Every doublet in `peaks` (sorted by m/z) for a precursor of charge
/// `precursor_charge`.
pub fn find_doublets(
    peaks: &[(f64, f32)],
    precursor_charge: u8,
    cleavage: &Cleavage,
    params: &DoubletParams,
) -> Vec<Doublet> {
    let mut out = Vec::new();
    if peaks.len() < 2 || precursor_charge < 2 {
        return out;
    }
    debug_assert!(peaks.windows(2).all(|w| w[0].0 <= w[1].0));
    let base = peaks.iter().map(|p| p.1).fold(0.0_f32, f32::max) as f64;
    let floor = base * params.min_rel_intensity;
    let max_z = params.max_charge.min(precursor_charge - 1).max(1);
    let delta = cleavage.doublet_delta();
    for (i, &(mz, inten)) in peaks.iter().enumerate() {
        if (inten as f64) < floor || inten <= 0.0 {
            continue;
        }
        for z in 1..=max_z {
            let zf = z as f64;
            let target = mz + delta / zf;
            let tol = ppm_window(target, params.tol_ppm);
            let Some(j) = find_peak(peaks, target, tol) else {
                continue;
            };
            if (peaks[j].1 as f64) < floor {
                continue;
            }
            let iso = ISOTOPE_SPACING / zf;
            // Isotope copy of a lower doublet: both members have a peak one
            // isotope spacing below.
            let below_light = find_peak(peaks, mz - iso, ppm_window(mz, params.tol_ppm));
            let below_heavy = find_peak(peaks, peaks[j].0 - iso, tol);
            if below_light.is_some() && below_heavy.is_some() {
                continue;
            }
            let above_light = find_peak(peaks, mz + iso, ppm_window(mz, params.tol_ppm));
            let above_heavy = find_peak(peaks, peaks[j].0 + iso, tol);
            out.push(Doublet {
                light_idx: i,
                heavy_idx: j,
                charge: z,
                peptide_mass: zf * (mz - PROTON) - cleavage.alkene,
                error_ppm: (peaks[j].0 - target) / target * 1e6,
                intensity: inten as f64 + peaks[j].1 as f64,
                isotope_confirmed: above_light.is_some() || above_heavy.is_some(),
            });
        }
    }
    out
}

/// Collapse doublets of the same peptide mass (within `tol_ppm`) seen at
/// several charges into one entry each, most intense first.
pub fn group_by_mass(doublets: &[Doublet], tol_ppm: f64) -> Vec<DoubletMass> {
    let mut sorted: Vec<&Doublet> = doublets.iter().collect();
    sorted.sort_by(|a, b| a.peptide_mass.total_cmp(&b.peptide_mass));
    let mut groups: Vec<Vec<&Doublet>> = Vec::new();
    for d in sorted {
        match groups.last_mut() {
            Some(g)
                if (d.peptide_mass - g[0].peptide_mass).abs()
                    <= ppm_window(d.peptide_mass, tol_ppm) =>
            {
                g.push(d)
            }
            _ => groups.push(vec![d]),
        }
    }
    let mut out: Vec<DoubletMass> = groups
        .into_iter()
        .map(|g| {
            let top = g
                .iter()
                .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
                .expect("non-empty group");
            let mut charges: Vec<u8> = g.iter().map(|d| d.charge).collect();
            charges.sort_unstable();
            charges.dedup();
            DoubletMass {
                mass: top.peptide_mass,
                intensity: g.iter().map(|d| d.intensity).sum(),
                n_charges: charges.len() as u8,
                isotope_confirmed: g.iter().any(|d| d.isotope_confirmed),
                min_charge: charges[0],
            }
        })
        .collect();
    out.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
    out
}

/// Which stub ions of a peptide of neutral mass `mass` are present at any
/// charge `1..=max_charge`: `[alkene, thiol, sulfenic]`.
pub fn stub_ions_present(
    peaks: &[(f64, f32)],
    mass: f64,
    cleavage: &Cleavage,
    max_charge: u8,
    tol_ppm: f64,
) -> [bool; 3] {
    let stubs = [cleavage.alkene, cleavage.thiol, cleavage.sulfenic];
    let mut out = [false; 3];
    for (k, stub) in stubs.iter().enumerate() {
        for z in 1..=max_charge.max(1) {
            let mz = (mass + stub) / z as f64 + PROTON;
            if find_peak(peaks, mz, ppm_window(mz, tol_ppm)).is_some() {
                out[k] = true;
                break;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crosslinker::Crosslinker;

    fn ion(mass: f64, z: u8) -> f64 {
        mass / z as f64 + PROTON
    }

    fn dsso() -> Cleavage {
        Crosslinker::dsso().cleavage.unwrap()
    }

    fn sorted(mut p: Vec<(f64, f32)>) -> Vec<(f64, f32)> {
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        p
    }

    #[test]
    fn finds_a_doublet_at_charge_two_and_recovers_the_mass() {
        let c = dsso();
        let m = 1234.5678;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 2), 100.0),
            (ion(m + c.thiol, 2), 80.0),
            (500.0, 10.0),
        ]);
        let d = find_doublets(&peaks, 4, &c, &DoubletParams::default());
        let z2: Vec<_> = d.iter().filter(|d| d.charge == 2).collect();
        assert_eq!(z2.len(), 1);
        assert!((z2[0].peptide_mass - m).abs() < 1e-6);
        assert!(!z2[0].isotope_confirmed);
        // No charge-1 doublet: the spacing at z1 is 31.97, not 15.99.
        assert!(d.iter().all(|d| d.charge != 1));
    }

    #[test]
    fn charge_is_capped_below_the_precursor_charge() {
        let c = dsso();
        let m = 2000.0;
        let peaks = sorted(vec![(ion(m + c.alkene, 3), 50.0), (ion(m + c.thiol, 3), 50.0)]);
        // A z3 doublet needs a precursor of charge >= 4.
        assert!(find_doublets(&peaks, 3, &c, &DoubletParams::default()).is_empty());
        assert_eq!(find_doublets(&peaks, 4, &c, &DoubletParams::default()).len(), 1);
    }

    #[test]
    fn isotope_copies_are_dropped_and_mark_the_monoisotope() {
        let c = dsso();
        let m = 1500.0;
        let i1 = ISOTOPE_SPACING;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 1), 100.0),
            (ion(m + c.thiol, 1), 100.0),
            (ion(m + c.alkene + i1, 1), 60.0),
            (ion(m + c.thiol + i1, 1), 60.0),
        ]);
        let d = find_doublets(&peaks, 3, &c, &DoubletParams::default());
        let z1: Vec<_> = d.iter().filter(|d| d.charge == 1).collect();
        assert_eq!(z1.len(), 1);
        assert!((z1[0].peptide_mass - m).abs() < 1e-6);
        assert!(z1[0].isotope_confirmed);
    }

    #[test]
    fn intensity_floor_removes_weak_members() {
        let c = dsso();
        let m = 1000.0;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 1), 1.0),
            (ion(m + c.thiol, 1), 1.0),
            (300.0, 1000.0),
        ]);
        let p = DoubletParams {
            min_rel_intensity: 0.01,
            ..DoubletParams::default()
        };
        assert!(find_doublets(&peaks, 3, &c, &p).is_empty());
        assert_eq!(
            find_doublets(&peaks, 3, &c, &DoubletParams::default()).len(),
            1
        );
    }

    #[test]
    fn grouping_merges_charges_of_one_peptide() {
        let c = dsso();
        let m = 1800.0;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 1), 10.0),
            (ion(m + c.thiol, 1), 10.0),
            (ion(m + c.alkene, 2), 30.0),
            (ion(m + c.thiol, 2), 30.0),
        ]);
        let d = find_doublets(&peaks, 4, &c, &DoubletParams::default());
        let g = group_by_mass(&d, 10.0);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].n_charges, 2);
        assert_eq!(g[0].min_charge, 1);
        assert!((g[0].intensity - 80.0).abs() < 1e-9);
    }

    #[test]
    fn stub_ions_are_found_at_any_charge() {
        let c = dsso();
        let m = 900.0;
        let peaks = sorted(vec![(ion(m + c.sulfenic, 2), 5.0)]);
        assert_eq!(stub_ions_present(&peaks, m, &c, 3, 10.0), [false, false, true]);
    }
}
