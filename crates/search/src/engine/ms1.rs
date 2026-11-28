//! The DIA MS1 scans a `PreparedSearch` keeps for its MS1 gate, and the precursor checks on them.

use std::cell::Cell;

use model::mass::PROTON;

use super::{DIA_MS1_GATE_SCANS, PreparedSearch};

thread_local! {
    static GATE_MZ: Cell<Vec<f64>> = const { Cell::new(Vec::new()) };
}

/// The DIA MS1 gate of one spectrum. `Peaks` holds, sorted, the m/z of the gate-scan peaks that
/// a candidate of the spectrum can match: the peaks inside the m/z range of its isolation
/// windows widened by twice the gate tolerance, and, when the spectrum has a mobility slice, on
/// scans with per-peak mobility only the peaks inside the slice. A query is then one binary
/// search over these few peaks instead of one search per gate scan over the whole scan.
/// `Scans` keeps the per-scan search for gate scans whose peaks are not sorted by finite m/z.
pub(super) enum DiaMs1Gate {
    Peaks(Vec<f64>),
    Scans(std::ops::Range<usize>, Option<(f32, f32)>),
}

impl Drop for DiaMs1Gate {
    fn drop(&mut self) {
        if let Self::Peaks(mz) = self {
            GATE_MZ.set(std::mem::take(mz));
        }
    }
}

fn sorted_by_finite_mz(peaks: &[(f64, f32)]) -> bool {
    peaks.iter().all(|p| p.0.is_finite()) && peaks.windows(2).all(|w| w[0].0 <= w[1].0)
}

impl<'a> PreparedSearch<'a> {
    pub fn set_dia_ms1(&mut self, mut scans: Vec<model::scan::Ms1Scan>) {
        scans.retain(|s| s.rt.is_finite());
        scans.sort_by(|a, b| a.rt.total_cmp(&b.rt));
        self.dia_ms1 = scans;
        self.check_dia_ms1_order();
    }

    fn check_dia_ms1_order(&mut self) {
        self.dia_ms1_sorted = self
            .dia_ms1
            .iter()
            .map(|s| sorted_by_finite_mz(&s.peaks))
            .collect();
    }

    /// The MS1 gate of `spec`, whose candidates are matched at charge `z` within neutral masses
    /// `lo..=hi` for each `(z, lo, hi)` of `windows`, at `tol_ppm` and in mobility slice `k0`;
    /// `None` when the spectrum has no gate scans.
    pub(super) fn dia_ms1_gate(
        &self,
        rt: Option<f64>,
        windows: &[(u8, f64, f64)],
        tol_ppm: f64,
        k0: Option<(f32, f32)>,
    ) -> Option<DiaMs1Gate> {
        let scans = self.dia_ms1_scans(rt)?;
        if !self
            .dia_ms1_sorted
            .get(scans.clone())
            .is_some_and(|sorted| sorted.iter().all(|&s| s))
        {
            return Some(DiaMs1Gate::Scans(scans, k0));
        }
        let (lo, hi) = windows.iter().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(a, b), &(z, lo, hi)| {
                let z = f64::from(z);
                (a.min(lo / z + PROTON), b.max(hi / z + PROTON))
            },
        );
        let widen = 2.0 * tol_ppm * 1e-6;
        let (lo, hi) = (lo * (1.0 - widen), hi * (1.0 + widen));
        let mut mz = GATE_MZ.take();
        mz.clear();
        for s in &self.dia_ms1[scans] {
            let a = s.peaks.partition_point(|p| p.0 < lo);
            let b = a + s.peaks[a..].partition_point(|p| p.0 <= hi);
            match k0 {
                Some((klo, khi)) if s.mobility.len() == s.peaks.len() => mz.extend(
                    s.peaks[a..b]
                        .iter()
                        .zip(&s.mobility[a..b])
                        .filter(|(_, &m)| m >= klo && m <= khi)
                        .map(|(p, _)| p.0),
                ),
                _ => mz.extend(s.peaks[a..b].iter().map(|p| p.0)),
            }
        }
        mz.sort_unstable_by(f64::total_cmp);
        Some(DiaMs1Gate::Peaks(mz))
    }

    /// Whether `gate` has a peak within `tol_ppm` of m/z `mz`: [`Self::dia_ms1_has_at`] on the
    /// gate's scans and mobility slice, for an `mz` whose isolation window built the gate.
    #[inline]
    pub(super) fn dia_ms1_gate_has(&self, gate: &DiaMs1Gate, mz: f64, tol_ppm: f64) -> bool {
        match gate {
            DiaMs1Gate::Peaks(peaks) => {
                let tol = mz * tol_ppm * 1e-6;
                let i = peaks.partition_point(|&p| p < mz - tol);
                peaks.get(i).is_some_and(|&p| p <= mz + tol)
            }
            DiaMs1Gate::Scans(scans, k0) => self.dia_ms1_has_at(scans.clone(), mz, tol_ppm, *k0),
        }
    }

    /// Installs per-DIA-window 1/K0 tables for the mobility term of the MS1 gate.
    pub fn set_dia_scan_mobility(&mut self, k0: input::DiaScanMobility) {
        self.dia_scan_mobility = k0;
    }

    /// The DIA MS1 gate's ppm tolerance, for learning the gate's mobility margin (`None`
    /// outside DIA).
    pub fn dia_ms1_gate_tol_ppm(&self) -> Option<f64> {
        self.params.dia.then_some(self.params.dia_ms1_tol_ppm)
    }

    /// Installs the mobility margin learned from the calibration anchors; `None` keeps the
    /// MS1 gate on m/z alone.
    pub fn set_dia_mobility_margin(&self, margin: Option<f32>) {
        if let Ok(mut m) = self.dia_mobility_margin.write() {
            *m = margin;
        }
    }

    /// The 1/K0 range of `spec`'s mobility slice (no margin), when the run has a mobility table.
    fn slice_k0(&self, spec: &model::Spectrum) -> Option<(f32, f32)> {
        dia_scan_mobility_range(&self.dia_scan_mobility, spec)
    }

    /// For a calibration anchor at precursor m/z `mz` in `spec`: the 1/K0 distance between
    /// `spec`'s mobility slice and the nearest MS1 peak within `tol_ppm` of `mz` in the gate
    /// scans (0 inside the slice); `None` without a slice, mobility data or matching peak.
    pub fn dia_ms1_mobility_excess(
        &self,
        spec: &model::Spectrum,
        mz: f64,
        tol_ppm: f64,
    ) -> Option<f32> {
        let (lo, hi) = self.slice_k0(spec)?;
        let scans = self.dia_ms1_scans(spec.rt_seconds)?;
        let tol = mz * tol_ppm * 1e-6;
        self.dia_ms1[scans]
            .iter()
            .filter(|s| s.mobility.len() == s.peaks.len())
            .flat_map(|s| {
                let i = s.peaks.partition_point(|p| p.0 < mz - tol);
                s.peaks[i..]
                    .iter()
                    .zip(&s.mobility[i..])
                    .take_while(move |(p, _)| p.0 <= mz + tol)
                    .map(|(_, &m)| (lo - m).max(m - hi).max(0.0))
            })
            .min_by(f32::total_cmp)
    }

    /// The 1/K0 range of `spec`'s mobility slice widened on each side by the learned margin,
    /// once a margin is learned and when the run has a mobility table.
    pub(super) fn dia_slice_mobility(&self, spec: &model::Spectrum) -> Option<(f32, f32)> {
        let margin = self.dia_mobility_margin.read().ok().and_then(|g| *g)?;
        let (lo, hi) = dia_scan_mobility_range(&self.dia_scan_mobility, spec)?;
        Some((lo - margin, hi + margin))
    }

    /// [`Self::dia_ms1_has`] with the matching peak's mobility inside `k0` when given (scans
    /// without per-peak mobility pass on m/z alone).
    pub(super) fn dia_ms1_has_at(
        &self,
        scans: std::ops::Range<usize>,
        mz: f64,
        tol_ppm: f64,
        k0: Option<(f32, f32)>,
    ) -> bool {
        let Some((lo, hi)) = k0 else {
            return self.dia_ms1_has(scans, mz, tol_ppm);
        };
        let tol = mz * tol_ppm * 1e-6;
        self.dia_ms1[scans].iter().any(|s| {
            let i = s.peaks.partition_point(|p| p.0 < mz - tol);
            if s.mobility.len() != s.peaks.len() {
                return s.peaks.get(i).is_some_and(|p| p.0 <= mz + tol);
            }
            s.peaks[i..]
                .iter()
                .zip(&s.mobility[i..])
                .take_while(|(p, _)| p.0 <= mz + tol)
                .any(|(_, &m)| m >= lo && m <= hi)
        })
    }

    pub fn take_dia_ms1(&mut self) -> Vec<model::scan::Ms1Scan> {
        self.dia_ms1_sorted = Vec::new();
        std::mem::take(&mut self.dia_ms1)
    }

    pub fn shift_dia_ms1_ppm(&mut self, ppm: f64) {
        let f = 1.0 - ppm * 1e-6;
        for s in &mut self.dia_ms1 {
            for p in &mut s.peaks {
                p.0 *= f;
            }
        }
        self.check_dia_ms1_order();
    }

    pub fn dia_ms1_error_ppm(&self, rt: Option<f64>, mz: f64, window_ppm: f64) -> Option<f64> {
        let tol = mz * window_ppm * 1e-6;
        self.dia_ms1_scans(rt)?
            .filter_map(|k| {
                let s = &self.dia_ms1[k];
                let i = s.peaks.partition_point(|p| p.0 < mz - tol);
                s.peaks[i..]
                    .iter()
                    .take_while(|p| p.0 <= mz + tol)
                    .max_by(|a, b| a.1.total_cmp(&b.1))
                    .copied()
            })
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|p| (p.0 - mz) / mz * 1e6)
    }

    pub(super) fn dia_ms1_has(&self, scans: std::ops::Range<usize>, mz: f64, tol_ppm: f64) -> bool {
        let tol = mz * tol_ppm * 1e-6;
        self.dia_ms1[scans].iter().any(|s| {
            let i = s.peaks.partition_point(|p| p.0 < mz - tol);
            s.peaks.get(i).is_some_and(|p| p.0 <= mz + tol)
        })
    }

    pub(super) fn dia_ms1_scans(&self, rt: Option<f64>) -> Option<std::ops::Range<usize>> {
        let rt = rt?;
        if self.dia_ms1.is_empty() {
            return None;
        }
        let i = self.dia_ms1.partition_point(|s| s.rt < rt);
        let near = match (i.checked_sub(1), (i < self.dia_ms1.len()).then_some(i)) {
            (Some(a), Some(b)) => {
                if rt - self.dia_ms1[a].rt <= self.dia_ms1[b].rt - rt {
                    a
                } else {
                    b
                }
            }
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => return None,
        };
        Some(
            near.saturating_sub(DIA_MS1_GATE_SCANS)
                ..(near + DIA_MS1_GATE_SCANS + 1).min(self.dia_ms1.len()),
        )
    }
}

fn dia_scan_mobility_range(
    mobility: &input::DiaScanMobility,
    spec: &model::Spectrum,
) -> Option<(f32, f32)> {
    let (a, b) = spec.mobility_window?;
    let scan = usize::try_from(spec.scan?).ok()?;
    let table_index = mobility.table_for_scan.get(scan)?.as_ref()?;
    let table = mobility.tables.get(*table_index)?;
    let ka = *table.get(a as usize)?;
    let kb = *table.get(b as usize)?;
    Some((ka.min(kb), ka.max(kb)))
}

#[cfg(test)]
mod tests {
    use super::dia_scan_mobility_range;
    use model::Spectrum;

    fn spectrum(scan: i32) -> Spectrum {
        Spectrum {
            scan: Some(scan),
            mobility_window: Some((1, 3)),
            ..Spectrum::default()
        }
    }

    #[test]
    fn dia_mobility_uses_the_spectrum_scan_table() {
        let mobility = input::DiaScanMobility {
            tables: vec![vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 10.0, 20.0, 30.0]],
            table_for_scan: vec![None, Some(0), Some(1)],
        };
        assert_eq!(
            dia_scan_mobility_range(&mobility, &spectrum(1)),
            Some((1.0, 3.0))
        );
        assert_eq!(
            dia_scan_mobility_range(&mobility, &spectrum(2)),
            Some((10.0, 30.0))
        );
        assert_eq!(dia_scan_mobility_range(&mobility, &spectrum(0)), None);
        assert_eq!(dia_scan_mobility_range(&mobility, &spectrum(3)), None);
    }
}
