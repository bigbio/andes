//! The DIA MS1 scans a `PreparedSearch` keeps for its MS1 gate, and the precursor checks on them.

use super::{DIA_MS1_GATE_SCANS, PreparedSearch};

impl<'a> PreparedSearch<'a> {
    pub fn set_dia_ms1(&mut self, mut scans: Vec<model::scan::Ms1Scan>) {
        scans.retain(|s| s.rt.is_finite());
        scans.sort_by(|a, b| a.rt.total_cmp(&b.rt));
        self.dia_ms1 = scans;
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
