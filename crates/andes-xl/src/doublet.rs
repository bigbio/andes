use crate::crosslinker::{Cleavage, ISOTOPE_SPACING, PROTON};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubletParams {
    pub tol_ppm: f64,
    pub max_charge: u8,
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Doublet {
    pub light_idx: usize,
    pub heavy_idx: usize,
    pub charge: u8,
    pub peptide_mass: f64,
    pub error_ppm: f64,
    pub intensity: f64,
    pub isotope_confirmed: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DoubletMass {
    pub mass: f64,
    pub intensity: f64,
    pub n_charges: u8,
    pub isotope_confirmed: bool,
    pub min_charge: u8,
}

#[inline]
fn ppm_window(mz: f64, ppm: f64) -> f64 {
    mz * ppm * 1e-6
}

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
    let spacings: Vec<(f64, f64)> = (0..cleavage.doublets.len())
        .map(|k| (cleavage.doublet_delta(k), cleavage.doublet_light(k)))
        .collect();
    for (i, &(mz, inten)) in peaks.iter().enumerate() {
        if (inten as f64) < floor || inten <= 0.0 {
            continue;
        }
        for (z, &(delta, light)) in (1..=max_z).flat_map(|z| spacings.iter().map(move |s| (z, s))) {
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
                peptide_mass: zf * (mz - PROTON) - light,
                error_ppm: (peaks[j].0 - target) / target * 1e6,
                intensity: inten as f64 + peaks[j].1 as f64,
                isotope_confirmed: above_light.is_some() || above_heavy.is_some(),
            });
        }
    }
    out
}

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

pub fn stub_ions_present(
    peaks: &[(f64, f32)],
    mass: f64,
    cleavage: &Cleavage,
    max_charge: u8,
    tol_ppm: f64,
) -> Vec<bool> {
    let mut out = vec![false; cleavage.stubs.len()];
    for (k, stub) in cleavage.stub_masses().enumerate() {
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

    struct Dsso {
        alkene: f64,
        thiol: f64,
        sulfenic: f64,
    }

    fn masses(c: &Cleavage) -> Dsso {
        Dsso {
            alkene: c.stubs[0].mass,
            thiol: c.stubs[1].mass,
            sulfenic: c.stubs[2].mass,
        }
    }

    fn sorted(mut p: Vec<(f64, f32)>) -> Vec<(f64, f32)> {
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        p
    }

    #[test]
    fn finds_a_doublet_at_charge_two_and_recovers_the_mass() {
        let cl = dsso();
        let c = masses(&cl);
        let m = 1234.5678;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 2), 100.0),
            (ion(m + c.thiol, 2), 80.0),
            (500.0, 10.0),
        ]);
        let d = find_doublets(&peaks, 4, &cl, &DoubletParams::default());
        let z2: Vec<_> = d.iter().filter(|d| d.charge == 2).collect();
        assert_eq!(z2.len(), 1);
        assert!((z2[0].peptide_mass - m).abs() < 1e-6);
        assert!(!z2[0].isotope_confirmed);
        assert!(d.iter().all(|d| d.charge != 1));
    }

    #[test]
    fn charge_is_capped_below_the_precursor_charge() {
        let cl = dsso();
        let c = masses(&cl);
        let m = 2000.0;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 3), 50.0),
            (ion(m + c.thiol, 3), 50.0),
        ]);
        assert!(find_doublets(&peaks, 3, &cl, &DoubletParams::default()).is_empty());
        assert_eq!(
            find_doublets(&peaks, 4, &cl, &DoubletParams::default()).len(),
            1
        );
    }

    #[test]
    fn isotope_copies_are_dropped_and_mark_the_monoisotope() {
        let cl = dsso();
        let c = masses(&cl);
        let m = 1500.0;
        let i1 = ISOTOPE_SPACING;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 1), 100.0),
            (ion(m + c.thiol, 1), 100.0),
            (ion(m + c.alkene + i1, 1), 60.0),
            (ion(m + c.thiol + i1, 1), 60.0),
        ]);
        let d = find_doublets(&peaks, 3, &cl, &DoubletParams::default());
        let z1: Vec<_> = d.iter().filter(|d| d.charge == 1).collect();
        assert_eq!(z1.len(), 1);
        assert!((z1[0].peptide_mass - m).abs() < 1e-6);
        assert!(z1[0].isotope_confirmed);
    }

    #[test]
    fn intensity_floor_removes_weak_members() {
        let cl = dsso();
        let c = masses(&cl);
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
        assert!(find_doublets(&peaks, 3, &cl, &p).is_empty());
        assert_eq!(
            find_doublets(&peaks, 3, &cl, &DoubletParams::default()).len(),
            1
        );
    }

    #[test]
    fn grouping_merges_charges_of_one_peptide() {
        let cl = dsso();
        let c = masses(&cl);
        let m = 1800.0;
        let peaks = sorted(vec![
            (ion(m + c.alkene, 1), 10.0),
            (ion(m + c.thiol, 1), 10.0),
            (ion(m + c.alkene, 2), 30.0),
            (ion(m + c.thiol, 2), 30.0),
        ]);
        let d = find_doublets(&peaks, 4, &cl, &DoubletParams::default());
        let g = group_by_mass(&d, 10.0);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].n_charges, 2);
        assert_eq!(g[0].min_charge, 1);
        assert!((g[0].intensity - 80.0).abs() < 1e-9);
    }

    #[test]
    fn stub_ions_are_found_at_any_charge() {
        let cl = dsso();
        let c = masses(&cl);
        let m = 900.0;
        let peaks = sorted(vec![(ion(m + c.sulfenic, 2), 5.0)]);
        assert_eq!(
            stub_ions_present(&peaks, m, &cl, 3, 10.0),
            vec![false, false, true]
        );
    }

    fn resolve(spec: &str) -> Cleavage {
        crate::spec::XlDecl::parse(spec)
            .and_then(|d| d.resolve(model::unimod::bundled()))
            .unwrap()
            .cleavage
            .unwrap()
    }

    #[test]
    fn other_cleavable_linkers_use_their_own_spacing() {
        for name in ["DSBU", "CDI", "DSBSO"] {
            let cl = resolve(name);
            let (l, h) = cl.doublets[0];
            let (light, heavy) = (cl.stubs[l].mass, cl.stubs[h].mass);
            let m = 1400.0;
            let peaks = sorted(vec![(ion(m + light, 2), 50.0), (ion(m + heavy, 2), 40.0)]);
            let d = find_doublets(&peaks, 4, &cl, &DoubletParams::default());
            assert_eq!(d.len(), 1, "{name}");
            assert!((d[0].peptide_mass - m).abs() < 1e-6, "{name}");
            assert!(find_doublets(&peaks, 4, &dsso(), &DoubletParams::default()).is_empty());
        }
    }

    #[test]
    fn every_doublet_of_a_linker_is_searched() {
        let mut cl = dsso();
        cl.doublets.push((0, 2));
        let c = masses(&cl);
        let (m1, m2) = (1100.0, 1700.0);
        let peaks = sorted(vec![
            (ion(m1 + c.alkene, 1), 30.0),
            (ion(m1 + c.thiol, 1), 30.0),
            (ion(m2 + c.alkene, 2), 20.0),
            (ion(m2 + c.sulfenic, 2), 20.0),
        ]);
        let d = find_doublets(&peaks, 4, &cl, &DoubletParams::default());
        let g = group_by_mass(&d, 10.0);
        let mut found: Vec<f64> = g.iter().map(|g| g.mass).collect();
        found.sort_by(f64::total_cmp);
        assert_eq!(found.len(), 2);
        assert!((found[0] - m1).abs() < 1e-6 && (found[1] - m2).abs() < 1e-6);
    }
}
