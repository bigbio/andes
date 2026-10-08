//! Scans that are not searched but quantified: MS1 survey scans (label-free
//! precursor features) and MS3 product scans (SPS-MS3 reporter ions).
//!
//! The readers capture these on request next to the MS2 stream; the
//! quantification layer consumes them. They carry only what quantification
//! needs, so a whole run's MS1 scans fit in memory (an Orbitrap Astral run is
//! ~30 M centroids, ~400 MB at 12 bytes each).

/// One MS1 scan: retention time in seconds and its centroids, m/z ascending.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Ms1Scan {
    pub rt: f64,
    pub peaks: Vec<(f64, f32)>,
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
    /// Centroids, m/z ascending.
    pub peaks: Vec<(f64, f32)>,
}
