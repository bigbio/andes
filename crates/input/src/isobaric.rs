//! Auto-detection of isobaric labeling (TMT / iTRAQ) from MS2 reporter ions.
//!
//! Reporter ions sit in the low-m/z region of (almost) every MS2 scan of an
//! isobaric-labeled run: TMT in the 126–131 cluster, iTRAQ in the 114–117
//! cluster. We sample MS2 spectra and, per spectrum, count how many of a label's
//! reporter channels carry a peak above an intensity floor; if a large fraction
//! of the sampled spectra show the cluster, the run is that label.
//!
//! This drives **zero-config** protocol selection: when the user leaves
//! `--protocol auto` (the default), a detected label engages the isobaric
//! windowed dense-peak filter automatically (the same path `--protocol TMT`
//! triggers today). It is **never** used to override an explicit `--protocol`,
//! and it returns `None` for label-free data, so behavior is unchanged unless a
//! genuine reporter cluster is present.
//!
//! Detection is intentionally presence-based and tolerant in m/z so it works on
//! both high-res MS2 (Orbitrap, reporters resolved to mDa) and low-res CID-MS2
//! (ion trap, reporters at ~unit resolution).

use model::spectrum::Spectrum;

/// An isobaric label detected from reporter ions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsobaricLabel {
    /// TMT (6/10/11/16/18-plex share the 126–131 reporter cluster).
    Tmt,
    /// iTRAQ (4/8-plex; the 114–117 reporter cluster is the 4-plex core).
    Itraq,
}

// Reporter-ion m/z, one representative per nominal channel.
const TMT_CHANNELS: &[f64] = &[126.1277, 127.1248, 128.1344, 129.1378, 130.1411, 131.1382];
const ITRAQ_CHANNELS: &[f64] = &[114.1112, 115.1083, 116.1116, 117.1150];

// Match window, resolution-dependent. High-res analyzers (Orbitrap/TOF) resolve
// reporters to a few mDa, so a TIGHT window is essential: the low-m/z region of
// real spectra is crowded with immonium ions and (for glyco) oxonium ions —
// e.g. HexNAc oxonium 126.0550 sits only 73 mDa from TMT 126.1277, and Arg
// immonium 129.1135 only 24 mDa from TMT 129.1378. A loose window lets those
// near-misses masquerade as reporters and false-trigger on label-free data.
// Low-res ion-trap reporters are only ~unit-accurate, so they need a loose
// window — but there the near-misses are unavoidable and the user can pass an
// explicit `--protocol`.
const HIGH_RES_TOL: f64 = 0.01;
const LOW_RES_TOL: f64 = 0.3;

/// Minimum number of MS2 (with peaks) needed before we trust a detection.
pub const MIN_SAMPLE: usize = 50;
/// A channel counts as present if its peak is ≥ this fraction of the base peak.
const REL_FLOOR: f32 = 0.01;
/// A spectrum "shows TMT" if ≥ this many of the 6 TMT channels are present.
const TMT_MIN_CHANNELS: usize = 4;
/// A spectrum "shows iTRAQ" if ≥ this many of the 4 iTRAQ channels are present.
const ITRAQ_MIN_CHANNELS: usize = 3;
/// The run is labeled if ≥ this fraction of sampled MS2 show the cluster.
const FRACTION_THRESHOLD: f64 = 0.30;

/// Count how many of `channels` have a peak within `tol`, at or above
/// `rel_floor × base_peak`, in `spec`.
fn channels_present(spec: &Spectrum, channels: &[f64], rel_floor: f32, tol: f64) -> usize {
    if spec.peaks.is_empty() {
        return 0;
    }
    let base = spec.peaks.iter().map(|&(_, i)| i).fold(0.0_f32, f32::max);
    if base <= 0.0 {
        return 0;
    }
    let floor = base * rel_floor;
    channels
        .iter()
        .filter(|&&c| {
            spec.peaks
                .iter()
                .any(|&(mz, i)| (mz - c).abs() <= tol && i >= floor)
        })
        .count()
}

/// Detect the isobaric label from a sample of MS2 spectra, or `None` if the run
/// is label-free or the sample is too small / ambiguous. `high_res` selects the
/// reporter match window (tight for Orbitrap/TOF, loose for ion-trap) — pass the
/// detected analyzer's resolution class.
///
/// TMT and iTRAQ reporter regions do not overlap (126–131 vs 114–117), so the
/// stronger cluster wins when (improbably) both clear the threshold.
pub fn detect_isobaric(spectra: &[Spectrum], high_res: bool) -> Option<IsobaricLabel> {
    let tol = if high_res { HIGH_RES_TOL } else { LOW_RES_TOL };
    let ms2: Vec<&Spectrum> = spectra.iter().filter(|s| !s.peaks.is_empty()).collect();
    if ms2.len() < MIN_SAMPLE {
        return None;
    }
    let n = ms2.len() as f64;

    let tmt_hits = ms2
        .iter()
        .filter(|s| channels_present(s, TMT_CHANNELS, REL_FLOOR, tol) >= TMT_MIN_CHANNELS)
        .count();
    let itraq_hits = ms2
        .iter()
        .filter(|s| channels_present(s, ITRAQ_CHANNELS, REL_FLOOR, tol) >= ITRAQ_MIN_CHANNELS)
        .count();

    let tmt_frac = tmt_hits as f64 / n;
    let itraq_frac = itraq_hits as f64 / n;

    if tmt_frac >= FRACTION_THRESHOLD && tmt_frac >= itraq_frac {
        Some(IsobaricLabel::Tmt)
    } else if itraq_frac >= FRACTION_THRESHOLD {
        Some(IsobaricLabel::Itraq)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec_with(peaks: Vec<(f64, f32)>) -> Spectrum {
        Spectrum {
            title: "t".into(),
            precursor_mz: 600.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }

    /// Backbone peaks shared by all synthetic spectra (away from reporter region).
    fn backbone() -> Vec<(f64, f32)> {
        vec![(300.0, 500.0), (450.0, 800.0), (620.0, 1000.0), (780.0, 400.0)]
    }

    fn tmt_spectrum() -> Spectrum {
        let mut p = backbone();
        // all 6 TMT reporters present, modest intensity
        for &c in TMT_CHANNELS {
            p.push((c, 200.0));
        }
        p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        spec_with(p)
    }

    fn itraq_spectrum() -> Spectrum {
        let mut p = backbone();
        for &c in ITRAQ_CHANNELS {
            p.push((c, 200.0));
        }
        p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        spec_with(p)
    }

    fn labelfree_spectrum() -> Spectrum {
        spec_with(backbone())
    }

    #[test]
    fn detects_tmt() {
        let run: Vec<Spectrum> = (0..200).map(|_| tmt_spectrum()).collect();
        assert_eq!(detect_isobaric(&run, true), Some(IsobaricLabel::Tmt));
    }

    #[test]
    fn detects_itraq() {
        let run: Vec<Spectrum> = (0..200).map(|_| itraq_spectrum()).collect();
        assert_eq!(detect_isobaric(&run, true), Some(IsobaricLabel::Itraq));
    }

    #[test]
    fn label_free_is_none() {
        let run: Vec<Spectrum> = (0..200).map(|_| labelfree_spectrum()).collect();
        assert_eq!(detect_isobaric(&run, true), None);
    }

    #[test]
    fn low_res_reporters_still_detected() {
        // Ion-trap reporters land at ~nominal masses (126.1, 127.1, ...); the
        // loose low-res window catches them.
        let run: Vec<Spectrum> = (0..200)
            .map(|_| {
                let mut p = backbone();
                for c in [126.1, 127.1, 128.1, 129.1, 130.1, 131.1] {
                    p.push((c, 150.0));
                }
                p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                spec_with(p)
            })
            .collect();
        assert_eq!(detect_isobaric(&run, false), Some(IsobaricLabel::Tmt));
    }

    #[test]
    fn sample_too_small_is_none() {
        let run: Vec<Spectrum> = (0..10).map(|_| tmt_spectrum()).collect();
        assert_eq!(detect_isobaric(&run, true), None);
    }

    #[test]
    fn sparse_coincidental_reporters_do_not_trigger() {
        // Only ~10% of spectra carry a couple of reporter-region peaks — below
        // both the per-spectrum channel count and the run fraction threshold.
        let mut run: Vec<Spectrum> = (0..180).map(|_| labelfree_spectrum()).collect();
        for _ in 0..20 {
            let mut p = backbone();
            p.push((126.1277, 100.0));
            p.push((127.1248, 100.0));
            p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            run.push(spec_with(p));
        }
        assert_eq!(detect_isobaric(&run, true), None);
    }

    #[test]
    fn reporters_below_floor_not_counted() {
        // Reporters present but at 0.1% of base peak (< 1% floor) → not counted.
        let run: Vec<Spectrum> = (0..200)
            .map(|_| {
                let mut p = backbone(); // base peak 1000
                for &c in TMT_CHANNELS {
                    p.push((c, 0.5)); // 0.05% of base → below floor
                }
                p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                spec_with(p)
            })
            .collect();
        assert_eq!(detect_isobaric(&run, true), None);
    }

    /// Regression for the real-data false positive: a high-res GLYCO run is
    /// label-free but peak-dense in the low-m/z region (HexNAc oxonium 126.055,
    /// Arg immonium 129.114, etc.). With a loose window those near-misses faked
    /// a TMT cluster; the tight high-res window must reject them.
    #[test]
    fn dense_highres_glyco_lowmz_is_not_tmt() {
        let run: Vec<Spectrum> = (0..300)
            .map(|_| {
                let mut p = backbone();
                // A peak in EACH nominal reporter slot (126–131) but offset ~50 mDa
                // off the exact TMT centroids — mimicking a dense glyco low-m/z
                // region. Loose (low-res) window catches all 6 → false TMT; tight
                // (high-res) window rejects all 6.
                for n in 126..=131 {
                    p.push((n as f64 + 0.05, 600.0));
                }
                // plus real glyco oxonium + Arg immonium (also off-centroid).
                for m in [126.0550, 138.0550, 144.0655, 204.0867, 129.1135] {
                    p.push((m, 600.0));
                }
                p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                spec_with(p)
            })
            .collect();
        assert_eq!(detect_isobaric(&run, true), None, "high-res glyco must not be read as TMT");
        // And with a loose (low-res) window this WOULD have false-triggered —
        // documents why high-res must use the tight window.
        assert_eq!(detect_isobaric(&run, false), Some(IsobaricLabel::Tmt));
    }
}
