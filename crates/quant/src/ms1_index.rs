use std::collections::HashMap;

use model::mass::ISOTOPE;
pub use model::scan::Ms1Scan;
use model::stats::{median_mut, Even};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobilityGate {
    pub lo: f32,
    pub hi: f32,
}

impl MobilityGate {
    pub fn contains(&self, k0: f32) -> bool {
        k0 >= self.lo && k0 <= self.hi
    }
}

pub fn slice_value(
    peaks: &[(f64, f32)],
    mobility: &[f32],
    start: usize,
    hi: f64,
    gate: Option<MobilityGate>,
) -> (f64, f32) {
    let (mut best_mz, mut best) = (0.0f64, 0.0f32);
    if mobility.is_empty() {
        for &(mz, inten) in &peaks[start.min(peaks.len())..] {
            if mz > hi {
                break;
            }
            if inten > best {
                best = inten;
                best_mz = mz;
            }
        }
        return (best_mz, best);
    }
    let mut sum = 0.0f32;
    for (k, &(mz, inten)) in peaks.iter().enumerate().skip(start) {
        if mz > hi {
            break;
        }
        if inten.is_nan() || inten <= 0.0 || gate.is_some_and(|g| !g.contains(mobility[k])) {
            continue;
        }
        sum += inten;
        if inten > best {
            best = inten;
            best_mz = mz;
        }
    }
    (best_mz, sum)
}

#[derive(Debug, Default)]
pub struct Ms1RunIndex {
    scans: Vec<Ms1Scan>,
    by_id: HashMap<String, usize>,
    detection: Vec<f32>,
}

const DETECTION_QUANTILE: f64 = 0.05;
const EMPTY_SCAN_DETECTION_LIMIT: f32 = 1.0;

fn scan_detection_limit(peaks: &[(f64, f32)], buf: &mut Vec<f32>) -> f32 {
    if peaks.is_empty() {
        return EMPTY_SCAN_DETECTION_LIMIT;
    }
    buf.clear();
    buf.extend(peaks.iter().map(|p| p.1));
    let rank = (DETECTION_QUANTILE * (buf.len() - 1) as f64).round() as usize;
    *buf.select_nth_unstable_by(rank, f32::total_cmp).1
}

impl Ms1RunIndex {
    pub fn new(mut scans: Vec<Ms1Scan>) -> Self {
        scans.retain(|s| s.rt.is_finite());
        for s in scans.iter_mut() {
            if !s.peaks.windows(2).all(|w| w[0].0 <= w[1].0) {
                if s.mobility.len() == s.peaks.len() {
                    let mut order: Vec<usize> = (0..s.peaks.len()).collect();
                    order.sort_by(|&a, &b| {
                        s.peaks[a]
                            .0
                            .partial_cmp(&s.peaks[b].0)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    s.peaks = order.iter().map(|&k| s.peaks[k]).collect();
                    s.mobility = order.iter().map(|&k| s.mobility[k]).collect();
                } else {
                    s.peaks
                        .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                }
            }
            if s.mobility.len() != s.peaks.len() {
                s.mobility.clear();
            }
        }
        scans.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));
        let mut buf = Vec::new();
        let detection = scans
            .iter()
            .map(|s| scan_detection_limit(&s.peaks, &mut buf))
            .collect();
        Self::assemble(scans, detection)
    }

    fn assemble(scans: Vec<Ms1Scan>, detection: Vec<f32>) -> Self {
        debug_assert_eq!(scans.len(), detection.len());
        let by_id = scans
            .iter()
            .enumerate()
            .filter(|(_, s)| !s.id.is_empty())
            .map(|(i, s)| (s.id.clone(), i))
            .collect();
        Self {
            scans,
            by_id,
            detection,
        }
    }

    pub fn position_of(&self, id: &str) -> Option<usize> {
        self.by_id.get(id).copied()
    }

    pub fn scan(&self, idx: usize) -> &Ms1Scan {
        &self.scans[idx]
    }

    pub fn faims_voltages(&self) -> Vec<f32> {
        let mut cvs: Vec<f32> = self.scans.iter().filter_map(|s| s.faims_cv).collect();
        cvs.sort_by(f32::total_cmp);
        cvs.dedup();
        cvs
    }

    pub fn split_by_faims(self) -> Vec<(f32, Ms1RunIndex)> {
        let cvs = self.faims_voltages();
        let mut groups: Vec<(Vec<Ms1Scan>, Vec<f32>)> = vec![(Vec::new(), Vec::new()); cvs.len()];
        for (scan, limit) in self.scans.into_iter().zip(self.detection) {
            if let Some(k) = scan
                .faims_cv
                .and_then(|cv| cvs.iter().position(|&c| c == cv))
            {
                groups[k].0.push(scan);
                groups[k].1.push(limit);
            }
        }
        cvs.into_iter()
            .zip(groups)
            .map(|(cv, (scans, detection))| (cv, Ms1RunIndex::assemble(scans, detection)))
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

    pub fn mobility(&self, idx: usize) -> &[f32] {
        &self.scans[idx].mobility
    }

    pub fn has_mobility(&self) -> bool {
        self.scans.iter().any(|s| !s.mobility.is_empty())
    }

    pub fn detection_limit(&self, idx: usize) -> f32 {
        self.detection[idx]
    }

    pub fn rt_range(&self) -> Option<(f64, f64)> {
        Some((self.scans.first()?.rt, self.scans.last()?.rt))
    }

    pub fn median_scan_spacing(&self) -> Option<f64> {
        let mut last: HashMap<Option<u32>, f64> = HashMap::new();
        let mut gaps: Vec<f64> = Vec::new();
        for s in &self.scans {
            if let Some(prev) = last.insert(s.faims_cv.map(f32::to_bits), s.rt) {
                if s.rt > prev {
                    gaps.push(s.rt - prev);
                }
            }
        }
        if gaps.is_empty() {
            return None;
        }
        median_mut(&mut gaps, Even::Upper)
    }

    pub fn scan_window(&self, rt_lo: f64, rt_hi: f64) -> std::ops::Range<usize> {
        let a = self.scans.partition_point(|s| s.rt < rt_lo);
        let b = self.scans.partition_point(|s| s.rt <= rt_hi);
        a..b.max(a)
    }

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

    pub fn scan_at_or_before(&self, rt: f64) -> Option<usize> {
        let i = self.scans.partition_point(|s| s.rt <= rt);
        i.checked_sub(1)
    }

    pub fn max_in(&self, idx: usize, lo: f64, hi: f64) -> f32 {
        let s = &self.scans[idx];
        let peaks = &s.peaks;
        let start = peaks.partition_point(|&(mz, _)| mz < lo);
        if !s.mobility.is_empty() {
            return slice_value(peaks, &s.mobility, start, hi, None).1;
        }
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

    pub fn slice_in(&self, idx: usize, lo: f64, hi: f64, gate: Option<MobilityGate>) -> (f64, f32) {
        let s = &self.scans[idx];
        let start = s.peaks.partition_point(|&(mz, _)| mz < lo);
        slice_value(&s.peaks, &s.mobility, start, hi, gate)
    }

    pub fn apex_centroid(&self, idx: usize, lo: f64, hi: f64) -> Option<(f64, f32)> {
        let s = &self.scans[idx];
        let peaks = &s.peaks;
        let start = peaks.partition_point(|&(mz, _)| mz < lo);
        if !s.mobility.is_empty() {
            let (mz, v) = slice_value(peaks, &s.mobility, start, hi, None);
            return (v > 0.0).then_some((mz, v));
        }
        let mut best: Option<(f64, f32)> = None;
        for &(mz, inten) in &peaks[start..] {
            if mz > hi {
                break;
            }
            if inten > best.map_or(0.0, |(_, b)| b) {
                best = Some((mz, inten));
            }
        }
        best
    }

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
    fn scans_without_a_finite_retention_time_are_dropped() {
        let scans = [1.0, f64::NAN, 3.0, f64::INFINITY, 2.0]
            .iter()
            .map(|&rt| Ms1Scan {
                rt,
                ..Ms1Scan::default()
            })
            .collect();
        let idx = Ms1RunIndex::new(scans);
        assert_eq!(idx.len(), 3);
        let rts: Vec<f64> = (0..idx.len()).map(|i| idx.scan(i).rt).collect();
        assert_eq!(rts, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn isotope_partners_separate_high_and_low_resolution_survey_scans() {
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
    fn apex_centroid_returns_position_and_intensity_of_the_most_intense() {
        let idx = Ms1RunIndex::new(vec![Ms1Scan {
            rt: 0.0,
            peaks: vec![(499.98, 5.0), (500.0, 7.0), (500.01, 7.0), (500.03, 9.0)],
            ..Ms1Scan::default()
        }]);
        assert_eq!(idx.apex_centroid(0, 499.97, 500.02), Some((500.0, 7.0)));
        assert_eq!(idx.apex_centroid(0, 499.97, 500.03), Some((500.03, 9.0)));
        assert_eq!(idx.apex_centroid(0, 600.0, 601.0), None);
        assert_eq!(idx.max_in(0, 499.97, 500.02), 7.0);
    }

    #[test]
    fn apex_centroid_skips_nan_and_non_positive_intensities() {
        let idx = Ms1RunIndex::new(vec![Ms1Scan {
            rt: 0.0,
            peaks: vec![
                (499.99, f32::NAN),
                (500.0, 0.0),
                (500.01, 3.0),
                (500.02, f32::NAN),
            ],
            ..Ms1Scan::default()
        }]);
        assert_eq!(idx.apex_centroid(0, 499.98, 500.03), Some((500.01, 3.0)));
        assert_eq!(idx.max_in(0, 499.98, 500.03), 3.0);
        assert_eq!(idx.apex_centroid(0, 499.98, 500.005), None);
    }

    #[test]
    fn detection_limit_is_per_scan_and_survives_the_faims_split() {
        let scans = vec![
            Ms1Scan {
                rt: 1.0,
                peaks: (1..=21)
                    .map(|k| (100.0 + k as f64, k as f32 * 10.0))
                    .collect(),
                faims_cv: Some(-45.0),
                ..Ms1Scan::default()
            },
            Ms1Scan {
                rt: 0.0,
                peaks: vec![],
                faims_cv: Some(-60.0),
                ..Ms1Scan::default()
            },
        ];
        let idx = Ms1RunIndex::new(scans);
        assert_eq!(idx.detection_limit(1), 20.0);
        assert_eq!(idx.detection_limit(0), 1.0);
        let split = idx.split_by_faims();
        assert_eq!(split[0].1.detection_limit(0), 1.0);
        assert_eq!(split[1].1.detection_limit(0), 20.0);
    }

    #[test]
    fn windows_and_nearest() {
        let idx = index();
        assert_eq!(idx.scan_window(3.0, 9.0), 2..5);
        assert_eq!(idx.scan_window(100.0, 200.0), 10..10);
        assert_eq!(idx.nearest_scan(5.1), Some(3));
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
                mobility: Vec::new(),
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
    fn scan_spacing_is_per_faims_voltage() {
        let scan = |rt: f64, cv: Option<f32>| Ms1Scan {
            rt,
            faims_cv: cv,
            ..Ms1Scan::default()
        };
        let alternating: Vec<Ms1Scan> = (0..8)
            .map(|i| scan(i as f64, Some(if i % 2 == 0 { -40.0 } else { -60.0 })))
            .collect();
        assert_eq!(
            Ms1RunIndex::new(alternating).median_scan_spacing(),
            Some(2.0)
        );
        let plain = vec![
            scan(0.0, None),
            scan(1.0, None),
            scan(2.5, None),
            scan(4.0, None),
        ];
        assert_eq!(Ms1RunIndex::new(plain).median_scan_spacing(), Some(1.5));
        assert_eq!(
            Ms1RunIndex::new(vec![scan(3.0, None)]).median_scan_spacing(),
            None
        );
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
