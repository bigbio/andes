//! A run's MS1 survey scans, indexed for extracted-ion chromatograms.
//!
//! Scans are kept whole (RT plus their m/z-sorted centroids) and ordered by
//! retention time. An XIC over an RT window is one binary search per scan,
//! which is cheap enough that no m/z binning is needed: a ±60 s window on an
//! Orbitrap run is ~100 scans, and a 10 ppm slice of each touches a handful
//! of centroids.

use std::collections::HashMap;

use model::mass::ISOTOPE;
pub use model::scan::Ms1Scan;

/// All MS1 scans of one run in RT order.
#[derive(Debug, Default)]
pub struct Ms1RunIndex {
    scans: Vec<Ms1Scan>,
    /// Native id → position, for scans that carry an id.
    by_id: HashMap<String, usize>,
}

impl Ms1RunIndex {
    /// Build from scans in any order; peaks are sorted by m/z if needed.
    pub fn new(mut scans: Vec<Ms1Scan>) -> Self {
        for s in scans.iter_mut() {
            if !s.peaks.windows(2).all(|w| w[0].0 <= w[1].0) {
                s.peaks
                    .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            }
        }
        scans.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));
        let by_id = scans
            .iter()
            .enumerate()
            .filter(|(_, s)| !s.id.is_empty())
            .map(|(i, s)| (s.id.clone(), i))
            .collect();
        Self { scans, by_id }
    }

    /// Position of the scan with native id `id`.
    pub fn position_of(&self, id: &str) -> Option<usize> {
        self.by_id.get(id).copied()
    }

    /// The scan at `idx`.
    pub fn scan(&self, idx: usize) -> &Ms1Scan {
        &self.scans[idx]
    }

    /// The distinct FAIMS compensation voltages of the run's scans, ascending
    /// (empty for a run without FAIMS).
    pub fn faims_voltages(&self) -> Vec<f32> {
        let mut cvs: Vec<f32> = self.scans.iter().filter_map(|s| s.faims_cv).collect();
        cvs.sort_by(f32::total_cmp);
        cvs.dedup();
        cvs
    }

    /// Split the run into one index per FAIMS compensation voltage. Survey scans
    /// at different voltages transmit different ions, so a chromatogram must
    /// stay within one voltage. Scans without a voltage are dropped.
    pub fn split_by_faims(self) -> Vec<(f32, Ms1RunIndex)> {
        let cvs = self.faims_voltages();
        let mut groups: Vec<Vec<Ms1Scan>> = vec![Vec::new(); cvs.len()];
        for scan in self.scans {
            if let Some(k) = scan
                .faims_cv
                .and_then(|cv| cvs.iter().position(|&c| c == cv))
            {
                groups[k].push(scan);
            }
        }
        cvs.into_iter()
            .zip(groups)
            .map(|(cv, g)| (cv, Ms1RunIndex::new(g)))
            .collect()
    }

    pub fn len(&self) -> usize {
        self.scans.len()
    }

    pub fn is_empty(&self) -> bool {
        self.scans.is_empty()
    }

    pub fn total_peaks(&self) -> usize {
        self.scans.iter().map(|s| s.peaks.len()).sum()
    }

    pub fn rt(&self, idx: usize) -> f64 {
        self.scans[idx].rt
    }

    pub fn peaks(&self, idx: usize) -> &[(f64, f32)] {
        &self.scans[idx].peaks
    }

    /// Retention-time span of the run, seconds.
    pub fn rt_range(&self) -> Option<(f64, f64)> {
        Some((self.scans.first()?.rt, self.scans.last()?.rt))
    }

    /// Indices of the scans with `rt_lo <= rt <= rt_hi`.
    pub fn scan_window(&self, rt_lo: f64, rt_hi: f64) -> std::ops::Range<usize> {
        let a = self.scans.partition_point(|s| s.rt < rt_lo);
        let b = self.scans.partition_point(|s| s.rt <= rt_hi);
        a..b.max(a)
    }

    /// Index of the scan whose RT is nearest to `rt`.
    pub fn nearest_scan(&self, rt: f64) -> Option<usize> {
        if self.scans.is_empty() {
            return None;
        }
        let i = self.scans.partition_point(|s| s.rt < rt);
        let cands = [i.checked_sub(1), (i < self.scans.len()).then_some(i)];
        cands.into_iter().flatten().min_by(|&a, &b| {
            (self.scans[a].rt - rt)
                .abs()
                .partial_cmp(&(self.scans[b].rt - rt).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Index of the last scan at or before `rt` (the survey scan an MS2 at
    /// `rt` was triggered from).
    pub fn scan_at_or_before(&self, rt: f64) -> Option<usize> {
        let i = self.scans.partition_point(|s| s.rt <= rt);
        i.checked_sub(1)
    }

    /// Most intense centroid in `[lo, hi]` of scan `idx`, 0.0 when none.
    pub fn max_in(&self, idx: usize, lo: f64, hi: f64) -> f32 {
        let peaks = &self.scans[idx].peaks;
        let start = peaks.partition_point(|&(mz, _)| mz < lo);
        let mut best = 0.0f32;
        for &(mz, inten) in &peaks[start..] {
            if mz > hi {
                break;
            }
            if inten > best {
                best = inten;
            }
        }
        best
    }

    /// Share of the most intense centroids that have an isotope partner
    /// (`+ISOTOPE/z`, z = 1..4) within `tol_ppm`, over an even sample of the
    /// run's scans. Orbitrap/TOF survey scans score 0.7–0.95 at 10 ppm; ion-trap
    /// survey scans, whose centroids are off by tenths of a Th, score far
    /// lower. `None` for an empty run.
    pub fn isotope_partner_fraction(&self, tol_ppm: f64) -> Option<f64> {
        const SAMPLE_SCANS: usize = 200;
        const TOP_PEAKS: usize = 30;
        let step = (self.scans.len() / SAMPLE_SCANS).max(1);
        let (mut hits, mut total) = (0usize, 0usize);
        for scan in self.scans.iter().step_by(step) {
            let mut order: Vec<usize> = (0..scan.peaks.len()).collect();
            order.sort_unstable_by(|&a, &b| {
                scan.peaks[b]
                    .1
                    .partial_cmp(&scan.peaks[a].1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            for &k in order.iter().take(TOP_PEAKS) {
                let mz = scan.peaks[k].0;
                total += 1;
                let partnered = (1..=4).any(|z| {
                    let target = mz + ISOTOPE / z as f64;
                    let tol = target * tol_ppm * 1e-6;
                    let i = scan.peaks.partition_point(|&(m, _)| m < target - tol);
                    scan.peaks.get(i).is_some_and(|&(m, _)| m <= target + tol)
                });
                if partnered {
                    hits += 1;
                }
            }
        }
        (total > 0).then(|| hits as f64 / total as f64)
    }

    /// Extracted-ion chromatogram of `mz ± tol_da` over `scans`: one value per
    /// scan (the most intense centroid in the slice, 0.0 when none).
    pub fn xic(&self, mz: f64, tol_da: f64, scans: std::ops::Range<usize>) -> Vec<f32> {
        scans
            .map(|i| self.max_in(i, mz - tol_da, mz + tol_da))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index() -> Ms1RunIndex {
        let scans = (0..10)
            .map(|i| Ms1Scan {
                rt: i as f64 * 2.0,
                peaks: vec![
                    (400.0, 1.0),
                    (500.0 + i as f64 * 0.001, (i * 10) as f32),
                    (600.0, 1.0),
                ],
                ..Ms1Scan::default()
            })
            .collect();
        Ms1RunIndex::new(scans)
    }

    #[test]
    fn isotope_partners_separate_high_and_low_resolution_survey_scans() {
        // High resolution: z=2 envelopes with exact isotope spacing.
        let hi = Ms1RunIndex::new(
            (0..20)
                .map(|i| {
                    let mut peaks = Vec::new();
                    for k in 0..10 {
                        let mono = 400.0 + 37.3 * k as f64 + i as f64 * 0.01;
                        for n in 0..3 {
                            peaks.push((mono + n as f64 * ISOTOPE / 2.0, 1000.0 / (n + 1) as f32));
                        }
                    }
                    Ms1Scan {
                        rt: i as f64,
                        peaks,
                        ..Ms1Scan::default()
                    }
                })
                .collect(),
        );
        assert!(hi.isotope_partner_fraction(10.0).unwrap() > 0.6);
        // Ion-trap-like: the same peaks with a tenth-of-a-Th centroid error.
        let lo = Ms1RunIndex::new(
            (0..20)
                .map(|i| Ms1Scan {
                    rt: i as f64,
                    peaks: hi
                        .peaks(i)
                        .iter()
                        .enumerate()
                        .map(|(j, &(mz, int))| (mz + 0.03 * ((j % 7) as f64 - 3.0), int))
                        .collect(),
                    ..Ms1Scan::default()
                })
                .collect(),
        );
        assert!(lo.isotope_partner_fraction(10.0).unwrap() < 0.2);
        assert_eq!(
            Ms1RunIndex::new(Vec::new()).isotope_partner_fraction(10.0),
            None
        );
    }

    #[test]
    fn windows_and_nearest() {
        let idx = index();
        assert_eq!(idx.scan_window(3.0, 9.0), 2..5);
        assert_eq!(idx.scan_window(100.0, 200.0), 10..10);
        assert_eq!(idx.nearest_scan(5.1), Some(3)); // rt 6 is nearer than 4
        assert_eq!(idx.nearest_scan(4.9), Some(2));
        assert_eq!(idx.scan_at_or_before(5.0), Some(2));
        assert_eq!(idx.scan_at_or_before(-1.0), None);
    }

    #[test]
    fn xic_follows_the_drifting_peak_within_tolerance() {
        let idx = index();
        let x = idx.xic(500.0, 0.02, 0..10);
        assert_eq!(x.len(), 10);
        assert_eq!(x[3], 30.0);
        // 500.009 at scan 9 is inside ±0.02 → found; outside a 0.005 window → 0.
        assert_eq!(x[9], 90.0);
        let tight = idx.xic(500.0, 0.005, 8..10);
        assert_eq!(tight, vec![0.0, 0.0]);
    }

    #[test]
    fn scans_are_found_by_id_and_split_by_faims_voltage() {
        let scans = (0..6)
            .map(|i| Ms1Scan {
                rt: i as f64,
                peaks: vec![(500.0, i as f32)],
                id: format!("scan={}", 10 - i),
                faims_cv: Some(if i % 2 == 0 { -40.0 } else { -60.0 }),
            })
            .collect();
        let idx = Ms1RunIndex::new(scans);
        assert_eq!(idx.position_of("scan=7"), Some(3));
        assert_eq!(idx.position_of("scan=99"), None);
        assert_eq!(idx.faims_voltages(), vec![-60.0, -40.0]);
        let split = idx.split_by_faims();
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].0, -60.0);
        let rts: Vec<f64> = (0..split[0].1.len()).map(|i| split[0].1.rt(i)).collect();
        assert_eq!(rts, vec![1.0, 3.0, 5.0]);
        assert_eq!(split[1].1.position_of("scan=8"), Some(1));
    }

    #[test]
    fn unsorted_input_is_normalised() {
        let scans = vec![
            Ms1Scan {
                rt: 5.0,
                peaks: vec![(300.0, 1.0), (200.0, 2.0)],
                ..Ms1Scan::default()
            },
            Ms1Scan {
                rt: 1.0,
                peaks: vec![],
                ..Ms1Scan::default()
            },
        ];
        let idx = Ms1RunIndex::new(scans);
        assert_eq!(idx.rt(0), 1.0);
        assert_eq!(idx.peaks(1)[0].0, 200.0);
        assert_eq!(idx.total_peaks(), 2);
    }
}
