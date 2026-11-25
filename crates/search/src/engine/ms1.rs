//! The DIA MS1 scans a `PreparedSearch` keeps for its MS1 gate, and the precursor checks on them.

use super::{PreparedSearch, DIA_MS1_GATE_SCANS};

impl<'a> PreparedSearch<'a> {
    pub fn set_dia_ms1(&mut self, mut scans: Vec<model::scan::Ms1Scan>) {
        scans.retain(|s| s.rt.is_finite());
        scans.sort_by(|a, b| a.rt.total_cmp(&b.rt));
        self.dia_ms1 = scans;
    }

    /// Installs the 1/K0 of every TIMS scan number, for the mobility term of the MS1 gate.
    pub fn set_dia_scan_mobility(&mut self, k0: Vec<f32>) {
        self.dia_scan_mobility = k0;
    }

    /// Whether the search asks for the MS1 gate's mobility margin learned from the run, and
    /// the gate's ppm tolerance.
    pub fn wants_learned_mobility_margin(&self) -> Option<f64> {
        self.params
            .dia_ms1_mobility_margin
            .is_some_and(|m| m < 0.0)
            .then_some(self.params.dia_ms1_tol_ppm)
    }

    /// Installs the mobility margin learned from the calibration anchors (used when the
    /// search parameters ask for the learned margin).
    pub fn set_dia_mobility_margin(&self, margin: Option<f32>) {
        *self.dia_mobility_margin.write().expect("margin lock") = margin;
    }

    /// The 1/K0 range of `spec`'s mobility slice (no margin), when the run has a scan
    /// mobility table.
    fn slice_k0(&self, spec: &model::Spectrum) -> Option<(f32, f32)> {
        let (a, b) = spec.mobility_window?;
        let ka = *self.dia_scan_mobility.get(a as usize)?;
        let kb = *self.dia_scan_mobility.get(b as usize)?;
        Some((ka.min(kb), ka.max(kb)))
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

    /// The 1/K0 range of `spec`'s mobility slice widened by `margin` on each side (a negative
    /// `margin`: the learned one, none until it is learned), when the run has a scan mobility
    /// table and the MS1 scans carry per-peak mobility.
    pub(super) fn dia_slice_mobility(
        &self,
        spec: &model::Spectrum,
        margin: Option<f64>,
    ) -> Option<(f32, f32)> {
        let margin = match margin? {
            m if m < 0.0 => (*self.dia_mobility_margin.read().expect("margin lock"))?,
            m => m as f32,
        };
        let (a, b) = spec.mobility_window?;
        let ka = *self.dia_scan_mobility.get(a as usize)?;
        let kb = *self.dia_scan_mobility.get(b as usize)?;
        Some((ka.min(kb) - margin, ka.max(kb) + margin))
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
        std::mem::take(&mut self.dia_ms1)
    }

    pub fn shift_dia_ms1_ppm(&mut self, ppm: f64) {
        let f = 1.0 - ppm * 1e-6;
        for s in &mut self.dia_ms1 {
            for p in &mut s.peaks {
                p.0 *= f;
            }
        }
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
