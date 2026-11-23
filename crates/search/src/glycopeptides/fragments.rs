use glyco::glycan_first::PPM_FLOOR_DA;
use rustc_hash::{FxHashMap, FxHashSet};

use model::peptide::Peptide;
use scoring_crate::scoring::fragment_ions::predict_by_ions;

pub struct FragmentIndex {
    bin_width: f64,
    tol: f64,
    tol_ppm: Option<f64>,
    bins: FxHashMap<i64, Vec<(u32, f64)>>,
}

const PPM_BIN_MAX_MZ: f64 = 5000.0;

impl FragmentIndex {
    #[inline]
    fn bin_of(mz: f64, bin_width: f64) -> i64 {
        (mz / bin_width).round() as i64
    }

    #[inline]
    fn ion_id(mz: f64) -> i64 {
        (mz * 1000.0).round() as i64
    }

    pub fn build<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        tol: f64,
        max_charge: u8,
    ) -> Self {
        let bin_width = tol.max(0.001);
        let zmax = max_charge.clamp(1, 3);
        let mut bins: FxHashMap<i64, Vec<(u32, f64)>> = FxHashMap::default();
        for (idx, pep) in entries {
            for ion in predict_by_ions(pep, 1..=zmax) {
                let b = Self::bin_of(ion.mz, bin_width);
                bins.entry(b).or_default().push((idx, ion.mz));
            }
        }
        FragmentIndex {
            bin_width,
            tol,
            tol_ppm: None,
            bins,
        }
    }

    pub fn build_ppm<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        tol_ppm: f64,
        max_charge: u8,
    ) -> Self {
        Self::build_ppm_sized(entries, tol_ppm, max_charge, PPM_BIN_MAX_MZ)
    }

    pub fn build_ppm_sized<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        tol_ppm: f64,
        max_charge: u8,
        bin_max_mz: f64,
    ) -> Self {
        if tol_ppm.is_nan() || tol_ppm <= 0.0 {
            return Self::build(entries, 0.5, max_charge);
        }
        let widest = (bin_max_mz * tol_ppm * 1e-6).max(PPM_FLOOR_DA);
        let mut idx = Self::build(entries, widest, max_charge);
        idx.tol_ppm = Some(tol_ppm);
        idx
    }

    #[inline]
    fn neighbourhood(&self, mz: f64) -> i64 {
        ((self.window_da(mz) / self.bin_width).ceil() as i64).max(1)
    }

    #[inline]
    fn window_da(&self, theo: f64) -> f64 {
        match self.tol_ppm {
            Some(p) => (theo * p * 1e-6).max(PPM_FLOOR_DA),
            None => self.tol,
        }
    }

    pub fn query(&self, peaks: &[(f64, f32)], min_matches: u32) -> Vec<(u32, u32)> {
        let mut used_ion: FxHashSet<(u32, i64)> = FxHashSet::default();
        let mut matched_this_peak: FxHashSet<u32> = FxHashSet::default();
        let mut count: FxHashMap<u32, u32> = FxHashMap::default();
        for &(mz, _) in peaks {
            matched_this_peak.clear();
            let b = Self::bin_of(mz, self.bin_width);
            let k = self.neighbourhood(mz);
            for nb in (b - k)..=(b + k) {
                if let Some(v) = self.bins.get(&nb) {
                    for &(idx, theo) in v {
                        let ion_key = (idx, Self::ion_id(theo));
                        if (mz - theo).abs() <= self.window_da(theo)
                            && !matched_this_peak.contains(&idx)
                            && !used_ion.contains(&ion_key)
                        {
                            matched_this_peak.insert(idx);
                            used_ion.insert(ion_key);
                            *count.entry(idx).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
        count
            .into_iter()
            .filter(|&(_, c)| c >= min_matches)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn ppm_query_reaches_matches_beyond_the_bin_sizing_mz() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDENK.R", &aa).expect("valid");
        let theo = predict_by_ions(&pep, 1..=1)
            .iter()
            .map(|i| i.mz)
            .fold(0.0f64, f64::max);
        assert!(
            theo > 300.0,
            "test premise: heaviest ion above the sizing m/z"
        );
        let idx = FragmentIndex::build_ppm_sized([(0u32, &pep)], 20.0, 1, 300.0);
        assert!(
            idx.neighbourhood(theo) > 1,
            "must need more than one bin per side here"
        );
        let peaks = vec![(theo * (1.0 + 18e-6), 100.0f32)];
        assert_eq!(idx.query(&peaks, 1), vec![(0, 1)]);
        let far = vec![(theo * (1.0 + 30e-6), 100.0f32)];
        assert!(idx.query(&far, 1).is_empty());
    }

    #[test]
    fn ppm_window_is_per_ion_and_much_narrower_than_the_da_window() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDENK.R", &aa).expect("valid");
        let ions = predict_by_ions(&pep, 1..=1);
        let theo = ions.iter().map(|i| i.mz).fold(0.0f64, f64::max);
        let da = FragmentIndex::build([(0u32, &pep)], 0.5, 1);
        let ppm = FragmentIndex::build_ppm([(0u32, &pep)], 20.0, 1);
        let off_by_0_3 = vec![(theo + 0.3, 100.0f32)];
        assert_eq!(
            da.query(&off_by_0_3, 1),
            vec![(0, 1)],
            "0.5 Da index must accept +0.3 Da"
        );
        assert!(
            ppm.query(&off_by_0_3, 1).is_empty(),
            "20 ppm index must reject +0.3 Da"
        );
        let inside = vec![(theo * (1.0 + 10e-6), 100.0f32)];
        assert_eq!(
            ppm.query(&inside, 1),
            vec![(0, 1)],
            "20 ppm index must accept +10 ppm"
        );
        let light = ions.iter().map(|i| i.mz).fold(f64::MAX, f64::min);
        let near_light = vec![(light + 0.004, 100.0f32)];
        assert_eq!(
            ppm.query(&near_light, 1),
            vec![(0, 1)],
            "floor must admit +4 mDa on a light ion"
        );
        let zero = FragmentIndex::build_ppm([(0u32, &pep)], 0.0, 1);
        assert_eq!(zero.query(&off_by_0_3, 1), vec![(0, 1)]);
    }

    use super::*;
    use model::aa_set::AminoAcidSetBuilder;

    fn aa() -> model::aa_set::AminoAcidSet {
        AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .build()
            .unwrap()
    }

    #[test]
    fn fragment_index_selects_peptide_with_matching_by_ions() {
        let aa = aa();
        let pep_a = Peptide::from_str("K.PEPTIDESK.R", &aa).expect("valid a");
        let pep_b = Peptide::from_str("K.ELVISLIVESR.K", &aa).expect("valid b");
        let idx = FragmentIndex::build([(0u32, &pep_a), (1u32, &pep_b)], 0.02, 1);

        let peaks: Vec<(f64, f32)> = predict_by_ions(&pep_a, 1..=1)
            .iter()
            .map(|i| (i.mz, 100.0))
            .collect();

        let hits = idx.query(&peaks, 3);
        let a = hits
            .iter()
            .find(|&&(i, _)| i == 0)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        let b = hits
            .iter()
            .find(|&&(i, _)| i == 1)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        assert!(a >= 3, "pep_a must match its own ladder (got {a})");
        assert!(a > b, "pep_a ({a}) must outscore the unrelated pep_b ({b})");
    }

    #[test]
    fn fragment_index_counts_tolerance_checked_distinct_ions() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDESK.R", &aa).unwrap();
        let tol = 0.02;
        let idx = FragmentIndex::build([(0u32, &pep)], tol, 1);
        let ions = predict_by_ions(&pep, 1..=1);

        let first = ions[0].mz;
        let mut peaks: Vec<(f64, f32)> =
            vec![(first, 100.0), (first + 0.005, 100.0), (first + 0.5, 100.0)];
        peaks.push((ions[1].mz, 100.0));
        peaks.push((ions[2].mz, 100.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let hits = idx.query(&peaks, 1);
        let c = hits
            .iter()
            .find(|&&(i, _)| i == 0)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        assert_eq!(c, 3, "must count 3 distinct tolerance-valid ions, got {c}");
        assert!(
            idx.query(&peaks, 4).is_empty(),
            "must not inflate past distinct-ion count"
        );
    }

    #[test]
    fn fragment_index_is_charge_aware() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDESKELVISR.R", &aa).unwrap();
        let peaks: Vec<(f64, f32)> = predict_by_ions(&pep, 2..=2)
            .iter()
            .map(|i| (i.mz, 100.0))
            .collect();
        let cnt = |z: u8| {
            FragmentIndex::build([(0u32, &pep)], 0.02, z)
                .query(&peaks, 1)
                .iter()
                .find(|&&(i, _)| i == 0)
                .map(|&(_, c)| c)
                .unwrap_or(0)
        };
        let c1 = cnt(1);
        let c2 = cnt(2);
        assert!(
            c2 >= 3,
            "+2 index must select the peptide from its +2 ladder (got {c2})"
        );
        assert!(
            c2 > c1,
            "+2 index must match more +2 ions than +1-only ({c2} > {c1})"
        );
    }

    #[test]
    fn fragment_index_count_never_exceeds_peak_count() {
        let aa = aa();
        let pep = Peptide::from_str("K.SNSSGGVG.R", &aa).unwrap();
        let idx = FragmentIndex::build([(0u32, &pep)], 0.02, 2);

        let mut theo: Vec<f64> = predict_by_ions(&pep, 1..=2).iter().map(|i| i.mz).collect();
        theo.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for m in theo {
            if peaks
                .last()
                .map(|&(p, _)| (m - p).abs() > 0.02)
                .unwrap_or(true)
            {
                peaks.push((m, 100.0));
            }
        }

        let n_peaks = peaks.len() as u32;
        let c = idx
            .query(&peaks, 1)
            .iter()
            .find(|&&(i, _)| i == 0)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        assert!(
            c <= n_peaks,
            "matched count {c} must not exceed observed peak count {n_peaks} (cross-charge inflation)"
        );
        assert!(
            idx.query(&peaks, n_peaks + 1).is_empty(),
            "must not pass a min_matches above the number of observed peaks"
        );
    }

    #[test]
    fn fragment_index_rejects_noise() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDESK.R", &aa).unwrap();
        let idx = FragmentIndex::build([(0u32, &pep)], 0.02, 1);
        let noise: Vec<(f64, f32)> = vec![(123.456, 1.0), (777.111, 1.0), (1500.9, 1.0)];
        let hits = idx.query(&noise, 3);
        assert!(
            hits.is_empty(),
            "noise must not select the peptide, got {hits:?}"
        );
    }
}
