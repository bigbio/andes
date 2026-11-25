//! The DIA MS1 scans a `PreparedSearch` keeps for its MS1 gate, and the precursor checks on them.

use super::{PreparedSearch, DIA_MS1_GATE_SCANS};

impl<'a> PreparedSearch<'a> {
    pub fn set_dia_ms1(&mut self, mut scans: Vec<model::scan::Ms1Scan>) {
        scans.retain(|s| s.rt.is_finite());
        scans.sort_by(|a, b| a.rt.total_cmp(&b.rt));
        self.dia_ms1 = scans;
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
