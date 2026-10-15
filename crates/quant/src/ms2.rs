//! MS2 evidence for match-between-runs transfers.
//!
//! A library entry carries an empirical fragment spectrum: the most intense
//! b/y fragments of its best donor PSM, matched on the donor scan. A probe
//! (a transfer target or one of its decoys) is checked against the MS2 scans
//! of the receiving run whose isolation window holds the probe's m/z and whose
//! retention time falls in the probe's fitted elution window. Each such scan
//! is scored against the library fragments; the best scan gives the probe's
//! [`Ms2Evidence`].
//!
//! A decoy probe is scored with the same library spectrum as its target, on
//! the scans of its own m/z and elution window. When the target peptide is
//! absent, the scans around the target's m/z hold other peptides, as the
//! decoy's do, so a decoy's score is the score of a false transfer.

use std::collections::HashMap;

use model::mass::{H2O, PROTON};
use model::{AminoAcid, Tolerance};

/// Fragments kept in a library spectrum (the most intense matched ones).
pub const LIBRARY_FRAGMENTS: usize = 12;
/// Fewest matched fragments of a library spectrum.
pub const MIN_LIBRARY_FRAGMENTS: usize = 3;
/// Shortest fragment (residues) in the ladder; shorter ions are unspecific.
pub const MIN_ION_LENGTH: usize = 3;
/// Half-width (Th) of the isolation window of a scan that does not record one.
pub const DEFAULT_ISOLATION_HALF_WIDTH: f64 = 1.0;
/// The elution window of a probe is its fitted apex ± this many fitted widths.
pub const ELUTION_WINDOW_WIDTHS: f64 = 3.0;
/// Bounds of the fitted width used for the elution window, seconds.
pub const ELUTION_WIDTH_BOUNDS_S: (f64, f64) = (1.0, 20.0);
/// MS2 support: at least this many matched library fragments ...
pub const SUPPORT_MIN_MATCHED: u8 = 4;
/// ... and at least this spectral angle.
pub const SUPPORT_MIN_ANGLE: f32 = 0.6;

/// b or y ion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IonKind {
    B,
    Y,
}

/// A fragment ion: series, length in residues and charge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Ion {
    pub kind: IonKind,
    pub length: u8,
    pub charge: u8,
}

impl Ion {
    /// `b4`, `y7`, `y9^2`.
    pub fn label(&self) -> String {
        let k = match self.kind {
            IonKind::B => 'b',
            IonKind::Y => 'y',
        };
        if self.charge == 1 {
            format!("{k}{}", self.length)
        } else {
            format!("{k}{}^{}", self.length, self.charge)
        }
    }

    /// The inverse of [`Ion::label`].
    pub fn parse(label: &str) -> Option<Ion> {
        let kind = match label.as_bytes().first()? {
            b'b' => IonKind::B,
            b'y' => IonKind::Y,
            _ => return None,
        };
        let rest = &label[1..];
        let (len, charge) = match rest.split_once('^') {
            Some((l, c)) => (l, c.parse().ok()?),
            None => (rest, 1u8),
        };
        let length: u8 = len.parse().ok()?;
        (length > 0 && charge > 0).then_some(Ion {
            kind,
            length,
            charge,
        })
    }
}

/// A theoretical fragment: its m/z and ion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fragment {
    pub mz: f64,
    pub ion: Ion,
}

/// One fragment of a library spectrum: theoretical m/z, the donor scan's
/// intensity relative to the most intense library fragment, and the ion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LibraryFragment {
    pub mz: f64,
    pub intensity: f32,
    pub ion: Ion,
}

/// An empirical fragment spectrum, fragments sorted by m/z.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LibrarySpectrum {
    pub fragments: Vec<LibraryFragment>,
}

/// Residue masses (modifications included) of a peptidoform written as
/// residues with `[...]` tags: a signed mass delta (`[+15.99491]`) or an
/// accession resolved through `mods` (`[UNIMOD:35]`). A tag before the first
/// residue goes to the first residue. Errors: an unknown residue or tag.
pub fn residue_masses(peptidoform: &str, mods: &HashMap<String, f64>) -> Result<Vec<f64>, String> {
    let mut masses: Vec<f64> = Vec::with_capacity(peptidoform.len());
    let mut pending = 0.0;
    let mut chars = peptidoform.chars();
    while let Some(c) = chars.next() {
        if c == '[' {
            let tag: String = chars.by_ref().take_while(|&c| c != ']').collect();
            let delta = if tag.starts_with('+') || tag.starts_with('-') {
                tag.parse::<f64>()
                    .map_err(|_| format!("{peptidoform}: unreadable mass tag [{tag}]"))?
            } else {
                *mods
                    .get(&tag)
                    .ok_or_else(|| format!("{peptidoform}: no mass for modification [{tag}]"))?
            };
            match masses.last_mut() {
                Some(m) => *m += delta,
                None => pending += delta,
            }
            continue;
        }
        if c == '-' && masses.is_empty() {
            continue;
        }
        let aa = u8::try_from(c)
            .ok()
            .and_then(AminoAcid::standard)
            .ok_or_else(|| format!("{peptidoform}: unknown residue {c:?}"))?;
        masses.push(aa.mass + std::mem::take(&mut pending));
    }
    if masses.is_empty() {
        return Err(format!("{peptidoform}: no residues"));
    }
    Ok(masses)
}

/// The b and y ladder of a peptide with residue masses `masses` at precursor
/// charge `precursor_charge`: fragments of [`MIN_ION_LENGTH`] residues up to
/// one residue short of the peptide, at charge 1, and also at charge 2 when
/// the precursor carries three or more charges. Sorted by m/z.
pub fn fragment_ladder(masses: &[f64], precursor_charge: u8) -> Vec<Fragment> {
    let n = masses.len();
    let max_charge = if precursor_charge >= 3 { 2 } else { 1 };
    let mut out = Vec::new();
    let mut prefix = 0.0;
    let total: f64 = masses.iter().sum();
    for (i, m) in masses.iter().enumerate().take(n.saturating_sub(1)) {
        prefix += m;
        let len_b = i + 1;
        let len_y = n - len_b;
        for charge in 1..=max_charge {
            let z = charge as f64;
            if len_b >= MIN_ION_LENGTH {
                out.push(Fragment {
                    mz: (prefix + z * PROTON) / z,
                    ion: Ion {
                        kind: IonKind::B,
                        length: len_b.min(u8::MAX as usize) as u8,
                        charge,
                    },
                });
            }
            if len_y >= MIN_ION_LENGTH {
                out.push(Fragment {
                    mz: (total - prefix + H2O + z * PROTON) / z,
                    ion: Ion {
                        kind: IonKind::Y,
                        length: len_y.min(u8::MAX as usize) as u8,
                        charge,
                    },
                });
            }
        }
    }
    out.sort_by(|a, b| a.mz.total_cmp(&b.mz));
    out
}

/// The most intense peak of `peaks` (sorted by m/z) within `tol` of `mz`.
fn best_peak(peaks: &[(f32, f32)], mz: f64, tol: Tolerance) -> Option<f32> {
    let d = tol.as_da(mz);
    let (lo, hi) = (mz - d, mz + d);
    let start = peaks.partition_point(|p| (p.0 as f64) < lo);
    peaks[start..]
        .iter()
        .take_while(|p| (p.0 as f64) <= hi)
        .map(|p| p.1)
        .filter(|i| i.is_finite() && *i > 0.0)
        .max_by(f32::total_cmp)
}

/// The library spectrum of a donor scan: the fragments of `ladder` with a
/// peak of `peaks` (sorted by m/z) within `tol`, the [`LIBRARY_FRAGMENTS`]
/// most intense of them (ties: the lower m/z), intensities relative to the
/// most intense. `None` with fewer than [`MIN_LIBRARY_FRAGMENTS`].
pub fn library_spectrum(
    ladder: &[Fragment],
    peaks: &[(f32, f32)],
    tol: Tolerance,
) -> Option<LibrarySpectrum> {
    let mut matched: Vec<(Fragment, f32)> = ladder
        .iter()
        .filter_map(|f| best_peak(peaks, f.mz, tol).map(|i| (*f, i)))
        .collect();
    if matched.len() < MIN_LIBRARY_FRAGMENTS {
        return None;
    }
    matched.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.mz.total_cmp(&b.0.mz)));
    matched.truncate(LIBRARY_FRAGMENTS);
    let top = matched[0].1;
    let mut fragments: Vec<LibraryFragment> = matched
        .into_iter()
        .map(|(f, i)| LibraryFragment {
            mz: f.mz,
            intensity: i / top,
            ion: f.ion,
        })
        .collect();
    fragments.sort_by(|a, b| a.mz.total_cmp(&b.mz));
    Some(LibrarySpectrum { fragments })
}

/// How one scan explains a library spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ScanMatch {
    /// Library fragments with a peak within the tolerance.
    pub n_matched: u8,
    /// Normalised spectral contrast angle between the square roots of the
    /// library intensities and of the matched peak intensities (0 for an
    /// unmatched fragment): `1 - 2 acos(cos) / pi`, in `[0, 1]`.
    pub angle: f32,
    /// Share of the library intensity carried by the matched fragments.
    pub explained: f32,
}

/// Scores scan `peaks` (sorted by m/z) against `lib`.
pub fn score_scan(lib: &LibrarySpectrum, peaks: &[(f32, f32)], tol: Tolerance) -> ScanMatch {
    let (mut dot, mut na, mut nb, mut total, mut hit) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut n = 0u8;
    for f in &lib.fragments {
        let a = (f.intensity.max(0.0) as f64).sqrt();
        na += a * a;
        total += f.intensity.max(0.0) as f64;
        if let Some(i) = best_peak(peaks, f.mz, tol) {
            let b = (i as f64).sqrt();
            dot += a * b;
            nb += b * b;
            hit += f.intensity.max(0.0) as f64;
            n = n.saturating_add(1);
        }
    }
    if n == 0 || na <= 0.0 || nb <= 0.0 {
        return ScanMatch::default();
    }
    let cos = (dot / (na.sqrt() * nb.sqrt())).clamp(0.0, 1.0);
    ScanMatch {
        n_matched: n,
        angle: (1.0 - 2.0 * cos.acos() / std::f64::consts::PI) as f32,
        explained: if total > 0.0 {
            (hit / total) as f32
        } else {
            0.0
        },
    }
}

/// One reduced MS2 scan of a run.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Ms2Scan {
    /// The scan's native id.
    pub native_id: String,
    /// The PSM SpecId of its rank-1 match without the `_<rank>[_<row>]`
    /// suffix: `"<title>_<scan number>"`.
    pub spec_key: String,
    pub rt: f64,
    pub precursor_mz: f64,
    /// 0 when the source records none.
    pub charge: u8,
    /// Isolation window as (target m/z, lower offset, upper offset), when the
    /// source records it.
    pub isolation: Option<(f64, f64, f64)>,
    pub faims_cv: Option<f32>,
    /// The most intense centroids, sorted by m/z.
    pub peaks: Vec<(f32, f32)>,
}

impl Ms2Scan {
    /// The isolated m/z range: the recorded window, else the precursor m/z ±
    /// [`DEFAULT_ISOLATION_HALF_WIDTH`].
    pub fn window(&self) -> (f64, f64) {
        match self.isolation {
            Some((t, lo, hi)) if t.is_finite() && lo.is_finite() && hi.is_finite() => {
                (t - lo.abs(), t + hi.abs())
            }
            _ => (
                self.precursor_mz - DEFAULT_ISOLATION_HALF_WIDTH,
                self.precursor_mz + DEFAULT_ISOLATION_HALF_WIDTH,
            ),
        }
    }
}

/// The best scan of a probe.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Ms2Evidence {
    /// Scans whose isolation window holds the probe's m/z in its elution
    /// window (the excluded ones not counted).
    pub n_scans: u16,
    /// The best of them by angle (then matched fragments); all zero without
    /// a scan.
    pub best: ScanMatch,
}

impl Ms2Evidence {
    /// The best scan has at least [`SUPPORT_MIN_MATCHED`] fragments at an
    /// angle of at least [`SUPPORT_MIN_ANGLE`].
    pub fn supports(&self) -> bool {
        self.best.n_matched >= SUPPORT_MIN_MATCHED && self.best.angle >= SUPPORT_MIN_ANGLE
    }
}

/// The elution window of a fit with apex `t0` and width `w`; `None` when
/// either is not finite.
pub fn elution_window(t0: f64, w: f64) -> Option<(f64, f64)> {
    if !t0.is_finite() || !w.is_finite() {
        return None;
    }
    let w = w.clamp(ELUTION_WIDTH_BOUNDS_S.0, ELUTION_WIDTH_BOUNDS_S.1);
    let h = ELUTION_WINDOW_WIDTHS * w;
    Some((t0 - h, t0 + h))
}

/// A run's MS2 scans in retention-time order.
#[derive(Debug, Clone, Default)]
pub struct Ms2RunIndex {
    scans: Vec<Ms2Scan>,
    windows: Vec<(f64, f64)>,
}

impl Ms2RunIndex {
    /// Sorts `scans` by retention time (stable; non-finite times dropped).
    pub fn new(mut scans: Vec<Ms2Scan>) -> Self {
        scans.retain(|s| s.rt.is_finite());
        scans.sort_by(|a, b| a.rt.total_cmp(&b.rt));
        let windows = scans.iter().map(Ms2Scan::window).collect();
        Self { scans, windows }
    }

    pub fn len(&self) -> usize {
        self.scans.len()
    }

    pub fn is_empty(&self) -> bool {
        self.scans.is_empty()
    }

    pub fn scans(&self) -> &[Ms2Scan] {
        &self.scans
    }

    /// Scan position by [`Ms2Scan::spec_key`].
    pub fn key_positions(&self) -> HashMap<&str, usize> {
        self.scans
            .iter()
            .enumerate()
            .map(|(i, s)| (s.spec_key.as_str(), i))
            .collect()
    }

    /// The best scan for `lib` among those isolating `mz` with a retention
    /// time in `rt` (inclusive), at FAIMS voltage `faims_cv` when given, and
    /// with a retention time other than those in `exclude_rt`.
    pub fn best_match(
        &self,
        lib: &LibrarySpectrum,
        mz: f64,
        rt: (f64, f64),
        faims_cv: Option<f32>,
        exclude_rt: &[f64],
        tol: Tolerance,
    ) -> Ms2Evidence {
        let mut out = Ms2Evidence::default();
        if !(mz.is_finite() && rt.0 <= rt.1) {
            return out;
        }
        let start = self.scans.partition_point(|s| s.rt < rt.0);
        for (s, w) in self.scans[start..].iter().zip(&self.windows[start..]) {
            if s.rt > rt.1 {
                break;
            }
            if !(w.0..=w.1).contains(&mz) {
                continue;
            }
            if faims_cv.is_some_and(|cv| s.faims_cv.is_some_and(|v| v != cv)) {
                continue;
            }
            if exclude_rt.iter().any(|t| t.to_bits() == s.rt.to_bits()) {
                continue;
            }
            out.n_scans = out.n_scans.saturating_add(1);
            let m = score_scan(lib, &s.peaks, tol);
            let better = m.angle > out.best.angle
                || (m.angle == out.best.angle && m.n_matched > out.best.n_matched);
            if better {
                out.best = m;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: Tolerance = Tolerance::Ppm(20.0);

    fn masses(seq: &str) -> Vec<f64> {
        residue_masses(seq, &HashMap::new()).unwrap()
    }

    fn peaks_of(frags: &[(f64, f32)]) -> Vec<(f32, f32)> {
        let mut p: Vec<(f32, f32)> = frags.iter().map(|&(m, i)| (m as f32, i)).collect();
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        p
    }

    #[test]
    fn residue_masses_resolve_deltas_and_accessions() {
        let mut mods = HashMap::new();
        mods.insert("UNIMOD:4".to_string(), 57.021464);
        let plain = masses("PEPCK");
        let m = residue_masses("PEPC[UNIMOD:4]K", &mods).unwrap();
        assert!((m[3] - plain[3] - 57.021464).abs() < 1e-9);
        let d = residue_masses("PEPCK[+8.01420]", &mods).unwrap();
        assert!((d[4] - plain[4] - 8.0142).abs() < 1e-9);
        let n = residue_masses("[+42.01057]-PEPCK", &mods).unwrap();
        assert!((n[0] - plain[0] - 42.01057).abs() < 1e-9);
        assert!(residue_masses("PEPC[UNIMOD:999]K", &mods)
            .unwrap_err()
            .contains("UNIMOD:999"));
        assert!(residue_masses("PEPZK", &mods).is_err());
    }

    #[test]
    fn ladder_matches_known_masses() {
        // PEPTIDEK: y1 of K is 147.1128; b3 of PEP is 324.1554.
        let l = fragment_ladder(&masses("PEPTIDEK"), 2);
        let b3 = l
            .iter()
            .find(|f| f.ion.label() == "b3")
            .expect("b3 in ladder");
        assert!((b3.mz - 324.155_4).abs() < 1e-3, "{}", b3.mz);
        let y3 = l.iter().find(|f| f.ion.label() == "y3").unwrap();
        // DEK + H2O + H+: 115.02694 + 129.04259 + 128.09496 + 18.01056 + 1.00728.
        assert!((y3.mz - 391.182_33).abs() < 1e-3, "{}", y3.mz);
        assert!(l.iter().all(|f| f.ion.charge == 1));
        assert!(l.iter().all(|f| f.ion.length as usize >= MIN_ION_LENGTH));
        assert!(l.windows(2).all(|w| w[0].mz <= w[1].mz));
        let l3 = fragment_ladder(&masses("PEPTIDEK"), 3);
        assert_eq!(l3.len(), 2 * l.len());
    }

    #[test]
    fn ion_labels_round_trip() {
        for label in ["b3", "y12", "y9^2", "b10^3"] {
            assert_eq!(Ion::parse(label).unwrap().label(), label);
        }
        for bad in ["", "a3", "y", "y0", "b3^0", "yx"] {
            assert!(Ion::parse(bad).is_none(), "{bad}");
        }
    }

    fn library() -> (Vec<Fragment>, LibrarySpectrum, Vec<(f32, f32)>) {
        let ladder = fragment_ladder(&masses("LGEYGFQNALIVR"), 2);
        // The donor scan: every third ladder fragment, with falling intensity,
        // plus unrelated noise.
        let mut frags: Vec<(f64, f32)> = ladder
            .iter()
            .step_by(2)
            .enumerate()
            .map(|(k, f)| (f.mz, 1000.0 / (k as f32 + 1.0)))
            .collect();
        frags.extend([(301.5, 50.0), (555.55, 20.0), (888.8, 30.0)]);
        let donor = peaks_of(&frags);
        let lib = library_spectrum(&ladder, &donor, TOL).unwrap();
        (ladder, lib, donor)
    }

    #[test]
    fn library_keeps_the_most_intense_matched_fragments() {
        let (_, lib, _) = library();
        assert!(lib.fragments.len() <= LIBRARY_FRAGMENTS);
        assert!(lib.fragments.len() >= MIN_LIBRARY_FRAGMENTS);
        assert!(lib.fragments.windows(2).all(|w| w[0].mz <= w[1].mz));
        let max = lib
            .fragments
            .iter()
            .map(|f| f.intensity)
            .fold(0.0, f32::max);
        assert_eq!(max, 1.0);
        // Too few matches: no library spectrum.
        let ladder = fragment_ladder(&masses("LGEYGFQNALIVR"), 2);
        let two = peaks_of(&[(ladder[0].mz, 10.0), (ladder[3].mz, 5.0)]);
        assert!(library_spectrum(&ladder, &two, TOL).is_none());
    }

    #[test]
    fn identical_scan_scores_one_and_unrelated_scan_zero() {
        let (_, lib, donor) = library();
        let same = score_scan(&lib, &donor, TOL);
        assert_eq!(same.n_matched as usize, lib.fragments.len());
        assert!((same.angle - 1.0).abs() < 1e-5, "{}", same.angle);
        assert!((same.explained - 1.0).abs() < 1e-6);
        let other = peaks_of(&[(400.1, 100.0), (512.3, 80.0), (777.7, 60.0)]);
        assert_eq!(score_scan(&lib, &other, TOL), ScanMatch::default());
        // Peaks just outside the tolerance do not match.
        let shifted: Vec<(f32, f32)> = donor.iter().map(|&(m, i)| (m * (1.0 + 30e-6), i)).collect();
        assert_eq!(score_scan(&lib, &shifted, TOL).n_matched, 0);
    }

    #[test]
    fn partial_and_distorted_matches_score_in_between() {
        let (_, lib, _) = library();
        let half: Vec<(f64, f32)> = lib
            .fragments
            .iter()
            .step_by(2)
            .map(|f| (f.mz, f.intensity * 100.0))
            .collect();
        let m = score_scan(&lib, &peaks_of(&half), TOL);
        assert_eq!(m.n_matched as usize, lib.fragments.len().div_ceil(2));
        assert!(m.angle > 0.0 && m.angle < 1.0, "{}", m.angle);
        assert!(m.explained > 0.0 && m.explained < 1.0);
        // All fragments, with reversed relative intensities.
        let n = lib.fragments.len();
        let mut by_int: Vec<usize> = (0..n).collect();
        by_int.sort_by(|&a, &b| {
            lib.fragments[a]
                .intensity
                .total_cmp(&lib.fragments[b].intensity)
        });
        let mut rev = vec![0.0f32; n];
        for (&lo, &hi) in by_int.iter().zip(by_int.iter().rev()) {
            rev[lo] = lib.fragments[hi].intensity;
        }
        let swapped: Vec<(f64, f32)> = lib
            .fragments
            .iter()
            .zip(&rev)
            .map(|(f, &i)| (f.mz, i))
            .collect();
        let s = score_scan(&lib, &peaks_of(&swapped), TOL);
        assert_eq!(s.n_matched as usize, n);
        assert!(s.angle < 0.95 && s.angle > 0.0, "{}", s.angle);
    }

    fn scan(rt: f64, mz: f64, peaks: Vec<(f32, f32)>) -> Ms2Scan {
        Ms2Scan {
            native_id: format!("scan={}", (rt * 10.0) as i64),
            spec_key: format!("t{rt}_1"),
            rt,
            precursor_mz: mz,
            charge: 2,
            isolation: Some((mz, 0.7, 0.7)),
            faims_cv: None,
            peaks,
        }
    }

    #[test]
    fn best_match_respects_window_isolation_faims_and_exclusion() {
        let (_, lib, donor) = library();
        let noise = vec![(410.2f32, 50.0f32), (620.3, 40.0)];
        let mut on_target = scan(100.0, 700.4, donor.clone());
        on_target.faims_cv = Some(-45.0);
        let index = Ms2RunIndex::new(vec![
            scan(130.0, 700.0, donor.clone()), // outside the RT window
            on_target,
            scan(95.0, 702.0, donor.clone()), // window does not hold 700
            scan(105.0, 700.3, noise.clone()),
        ]);
        let ev = index.best_match(&lib, 700.0, (90.0, 110.0), None, &[], TOL);
        assert_eq!(ev.n_scans, 2);
        assert!((ev.best.angle - 1.0).abs() < 1e-5);
        assert!(ev.supports());
        // The scan at another FAIMS voltage is left out.
        let ev = index.best_match(&lib, 700.0, (90.0, 110.0), Some(-60.0), &[], TOL);
        assert_eq!(ev.n_scans, 1);
        assert_eq!(ev.best, ScanMatch::default());
        assert!(!ev.supports());
        // So is an excluded scan.
        let ev = index.best_match(&lib, 700.0, (90.0, 110.0), None, &[100.0], TOL);
        assert_eq!(ev.n_scans, 1);
        assert_eq!(ev.best.n_matched, 0);
        // No scan: zeros.
        let ev = index.best_match(&lib, 900.0, (90.0, 110.0), None, &[], TOL);
        assert_eq!(ev, Ms2Evidence::default());
        let ev = index.best_match(&lib, f64::NAN, (90.0, 110.0), None, &[], TOL);
        assert_eq!(ev, Ms2Evidence::default());
    }

    #[test]
    fn decoy_and_target_scored_alike() {
        // A target probe and a mass-shift decoy probe on scans of the same
        // composition (each scan isolates its probe's m/z and holds the same
        // unrelated peaks) get the same evidence: the score depends on the
        // scans, not on which probe asks.
        let (_, lib, _) = library();
        let other = vec![(333.3f32, 70.0f32), (444.4, 60.0), (555.5, 50.0)];
        let index = Ms2RunIndex::new(vec![
            scan(100.0, 700.0, other.clone()),
            scan(100.5, 705.5, other.clone()),
        ]);
        let t = index.best_match(&lib, 700.0, (90.0, 110.0), None, &[], TOL);
        let d = index.best_match(&lib, 705.5, (90.0, 110.0), None, &[], TOL);
        assert_eq!(t, d);
        assert_eq!(t.n_scans, 1);
    }

    #[test]
    fn default_isolation_and_elution_windows() {
        let mut s = scan(10.0, 500.0, vec![]);
        s.isolation = None;
        assert_eq!(
            s.window(),
            (
                500.0 - DEFAULT_ISOLATION_HALF_WIDTH,
                500.0 + DEFAULT_ISOLATION_HALF_WIDTH
            )
        );
        s.isolation = Some((500.2, 0.6, 0.8));
        let (lo, hi) = s.window();
        assert!((lo - 499.6).abs() < 1e-9 && (hi - 501.0).abs() < 1e-9);
        assert_eq!(elution_window(100.0, 5.0), Some((85.0, 115.0)));
        assert_eq!(elution_window(100.0, 0.1), Some((97.0, 103.0)));
        assert_eq!(elution_window(100.0, 100.0), Some((40.0, 160.0)));
        assert_eq!(elution_window(f64::NAN, 5.0), None);
        // Unsorted input and non-finite RTs.
        let idx = Ms2RunIndex::new(vec![
            scan(3.0, 1.0, vec![]),
            scan(f64::NAN, 1.0, vec![]),
            scan(1.0, 1.0, vec![]),
        ]);
        assert_eq!(idx.len(), 2);
        assert!(idx.scans()[0].rt < idx.scans()[1].rt);
        assert_eq!(idx.key_positions()["t1_1"], 0);
    }
}
