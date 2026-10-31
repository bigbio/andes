//! Scans that are not searched but quantified: MS1 survey scans (label-free
//! precursor features) and MS3 product scans (SPS-MS3 reporter ions), plus the
//! acquisition context of each searched MS2 that the [`crate::spectrum::Spectrum`]
//! itself does not carry.
//!
//! The readers capture these on request next to the MS2 stream; the
//! quantification layer consumes them. They carry only what quantification
//! needs, so a whole run's MS1 scans fit in memory (an Orbitrap Astral run is
//! ~30 M centroids, ~480 MB at 16 bytes each: `(f64, f32)` pads to 16).

/// One MS1 scan: retention time in seconds and its centroids, m/z ascending.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Ms1Scan {
    pub rt: f64,
    pub peaks: Vec<(f64, f32)>,
    /// Native id (`<spectrum id>` in mzML, the controller id in `.raw`); empty
    /// when the source has none.
    pub id: String,
    /// FAIMS compensation voltage of the scan, when the run used FAIMS. Survey
    /// scans at different voltages sample different ion populations and are
    /// not one chromatographic series.
    pub faims_cv: Option<f32>,
    /// Inverse reduced ion mobility (1/K0, V·s/cm²) of each centroid, parallel
    /// to `peaks`, when the source resolves the survey scan by ion mobility
    /// (Bruker timsTOF `.d` read with mobility-resolved MS1); empty otherwise.
    /// A mobility-resolved scan may hold several centroids at one m/z: the
    /// ion-mobility components of one time-of-flight peak.
    pub mobility: Vec<f32>,
}

/// An MS3 (or higher) product scan kept for reporter-ion quantification.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ProductScan {
    /// Native id (`<spectrum id>` in mzML, the controller id in `.raw`).
    pub id: String,
    pub scan: Option<i32>,
    pub rt_seconds: Option<f64>,
    pub ms_level: u8,
    /// Native id of the scan this one was produced from (mzML
    /// `<precursor spectrumRef>`; the master scan in `.raw`). The MS2 whose
    /// PSM this product scan quantifies.
    pub parent_id: Option<String>,
    /// Selected precursor m/z of this scan (the last SPS notch listed for an
    /// MS3; informational).
    pub precursor_mz: f64,
    /// Acquired m/z range (`scan window lower/upper limit`), when recorded.
    pub scan_window: Option<(f64, f64)>,
    /// Centroids, m/z ascending.
    pub peaks: Vec<(f64, f32)>,
    /// Noise level in the reporter-ion region, when the source records it
    /// and noise capture is on (reporter signal-to-noise).
    pub noise: Option<NoiseProfile>,
}

/// Acquisition context of one emitted MS2 scan, keyed by its native id.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ScanContext {
    /// Native id of the MS2 (the emitted spectrum's title).
    pub id: String,
    /// Native id of the survey scan the MS2 was triggered from (mzML
    /// `<precursor spectrumRef>`; the `Master Scan Number` trailer in `.raw`).
    pub parent_id: Option<String>,
    /// Acquired m/z range (`scan window lower/upper limit`), when recorded.
    pub scan_window: Option<(f64, f64)>,
    /// FAIMS compensation voltage, when the run used FAIMS.
    pub faims_cv: Option<f32>,
    /// Noise level in the reporter-ion region of the MS2, when the source
    /// records it and noise capture is on (reporter signal-to-noise at MS2).
    pub noise: Option<NoiseProfile>,
    /// Inverse reduced ion mobility (1/K0, V·s/cm²) of the precursor, when
    /// the source records it (Bruker timsTOF `.d`).
    pub ion_mobility: Option<f32>,
}

/// Upper end of the m/z range a [`NoiseProfile`] keeps: every isobaric
/// reporter ion (iTRAQ 113-121, TMT 126-135) lies below it.
pub const REPORTER_NOISE_MAX_MZ: f64 = 140.0;

/// Noise level along m/z of one scan, as the instrument estimated it: Orbitrap
/// noise per centroid (Thermo `.raw`) or the sampled noise packets an mzML
/// converter writes (`sampled noise m/z array` / `sampled noise intensity
/// array`). Reporter signal-to-noise is a reporter's intensity over the noise
/// at its m/z.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NoiseProfile {
    /// `(m/z, noise)`, m/z ascending, noise > 0.
    points: Vec<(f64, f32)>,
}

impl NoiseProfile {
    /// Keep the finite, positive points up to [`REPORTER_NOISE_MAX_MZ`] (and
    /// the first point above it, so the last reporter still interpolates).
    /// `None` when no usable point is left.
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

    /// The noise at `mz`: linear between the two nearest points, the end value
    /// outside the sampled range.
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
        // 110 (zero) and 125 (NaN) dropped; 300 is past the first point above 140.
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
