#[derive(Debug, Clone, Default, PartialEq)]
pub struct Ms1Scan {
    pub rt: f64,
    pub peaks: Vec<(f64, f32)>,
    pub id: String,
    pub faims_cv: Option<f32>,
    pub mobility: Vec<f32>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ProductScan {
    pub id: String,
    pub scan: Option<i32>,
    pub rt_seconds: Option<f64>,
    pub ms_level: u8,
    pub parent_id: Option<String>,
    pub precursor_mz: f64,
    pub scan_window: Option<(f64, f64)>,
    pub peaks: Vec<(f64, f32)>,
    pub noise: Option<NoiseProfile>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ScanContext {
    pub id: String,
    pub parent_id: Option<String>,
    pub scan_window: Option<(f64, f64)>,
    pub faims_cv: Option<f32>,
    pub noise: Option<NoiseProfile>,
    pub ion_mobility: Option<f32>,
}

pub const REPORTER_NOISE_MAX_MZ: f64 = 140.0;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct NoiseProfile {
    points: Vec<(f64, f32)>,
}

impl NoiseProfile {
    pub fn reporter_region(points: impl IntoIterator<Item = (f64, f32)>) -> Option<Self> {
        let mut points: Vec<(f64, f32)> = points
            .into_iter()
            .filter(|&(mz, n)| mz.is_finite() && mz > 0.0 && n.is_finite() && n > 0.0)
            .collect();
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        let above = points.partition_point(|&(mz, _)| mz <= REPORTER_NOISE_MAX_MZ);
        points.truncate((above + 1).min(points.len()));
        (!points.is_empty()).then_some(Self { points })
    }

    pub fn at(&self, mz: f64) -> f32 {
        let p = &self.points;
        let i = p.partition_point(|&(x, _)| x < mz);
        if i == 0 {
            return p[0].1;
        }
        if i == p.len() {
            return p[p.len() - 1].1;
        }
        let ((x0, y0), (x1, y1)) = (p[i - 1], p[i]);
        if x1 <= x0 {
            return y1;
        }
        let t = ((mz - x0) / (x1 - x0)) as f32;
        y0 + t * (y1 - y0)
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noise_profile_keeps_the_reporter_region_and_interpolates() {
        let p = NoiseProfile::reporter_region([
            (130.0, 200.0),
            (120.0, 100.0),
            (125.0, f32::NAN),
            (150.0, 400.0),
            (300.0, 900.0),
            (110.0, 0.0),
        ])
        .unwrap();
        assert_eq!(p.len(), 3);
        assert_eq!(p.at(100.0), 100.0);
        assert_eq!(p.at(120.0), 100.0);
        assert!((p.at(125.0) - 150.0).abs() < 1e-4);
        assert!((p.at(140.0) - 300.0).abs() < 1e-4);
        assert_eq!(p.at(500.0), 400.0);
        assert!(NoiseProfile::reporter_region([(126.0, 0.0)]).is_none());
        assert!(NoiseProfile::reporter_region([]).is_none());
    }
}
