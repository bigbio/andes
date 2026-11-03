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
//!
//! Presence is judged against the run's own background. The low-m/z region of
//! a label-free ion-trap HCD run holds a peak at most nominal masses, so at the
//! loose low-res window a 4-of-6 "TMT cluster" turns up by chance in over half
//! of its spectra. Each rule is therefore also evaluated on the channel set
//! shifted by whole Daltons ([`CONTROL_SHIFTS_DA`], positions where no reporter
//! sits): their mean is the chance rate `c`, and a fraction `f` counts as
//! `(f - c) / (1 - c)`, the share of the spectra not explained by chance. On
//! high-res data, where `c` is about zero, this is the plain fraction. The
//! same correction is applied per channel, and a cluster also needs as many
//! channels above chance as a spectrum needs present, so label-free ions that
//! sit at one or two reporter masses do not make a cluster.

use model::spectrum::Spectrum;

/// An isobaric label detected from reporter ions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsobaricLabel {
    /// TMT 6/10/11-plex: the 126–131 reporter cluster only.
    Tmt,
    /// TMTpro 16/18/32/35-plex: the 126–131 cluster plus the 132–135 reporters.
    TmtPro,
    /// iTRAQ 4-plex: the 114–117 reporter cluster.
    Itraq,
    /// iTRAQ 8-plex: the 114–117 cluster plus 113, 118, 119 and 121.
    Itraq8,
}

impl IsobaricLabel {
    /// The Unimod name and monoisotopic mass of the label added to K and the
    /// peptide N-terminus.
    pub fn tag(self) -> (&'static str, f64) {
        match self {
            IsobaricLabel::Tmt => ("TMT6plex", 229.162932),
            IsobaricLabel::TmtPro => ("TMTpro", 304.207146),
            IsobaricLabel::Itraq => ("iTRAQ4plex", 144.102063),
            IsobaricLabel::Itraq8 => ("iTRAQ8plex", 304.205360),
        }
    }
}

// Reporter-ion m/z, one representative per nominal channel.
const TMT_CHANNELS: &[f64] = &[126.1277, 127.1248, 128.1344, 129.1378, 130.1411, 131.1382];
const ITRAQ_CHANNELS: &[f64] = &[114.1112, 115.1083, 116.1116, 117.1150];
// Reporters only the larger kits carry. TMTpro 132N/132C/133N/133C/134N
// (TMT 6–11 stop at 131); iTRAQ 8-plex 113/118/119/121 (4-plex is 114–117).
// At the loose low-res window the N/C pairs merge, so the count still reaches
// the minimum below from the nominal 132–134 peaks.
const TMTPRO_EXTRA: &[f64] = &[132.1415, 132.1479, 133.1449, 133.1512, 134.1482];
const ITRAQ8_EXTRA: &[f64] = &[113.1079, 118.1120, 119.1153, 121.1220];
/// A labeled spectrum shows the larger kit if ≥ this many extra reporters are present.
const EXTRA_MIN_CHANNELS: usize = 2;
/// The run is the larger kit if ≥ this fraction of its labeled spectra show it.
const EXTRA_FRACTION: f64 = 0.5;

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

/// Whole-Dalton offsets of the control channel sets that measure the run's
/// chance rate. At least 9 Da on either side, so a shifted set never overlaps
/// the 8-Da span of a kit's own reporters (TMTpro 126–134, iTRAQ 8-plex
/// 113–121); on both sides, so the local peak density is sampled above and
/// below the reporters.
const CONTROL_SHIFTS_DA: [i32; 24] = [
    -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    18, 19, 20,
];

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

/// Count how many of `channels`, moved by `shift` Da, have a peak within `tol`,
/// at or above `rel_floor × base_peak`, in `spec`.
fn channels_present_at(
    spec: &Spectrum,
    channels: &[f64],
    shift: f64,
    rel_floor: f32,
    tol: f64,
) -> usize {
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
                .any(|&(mz, i)| (mz - (c + shift)).abs() <= tol && i >= floor)
        })
        .count()
}

fn channels_present(spec: &Spectrum, channels: &[f64], rel_floor: f32, tol: f64) -> usize {
    channels_present_at(spec, channels, 0.0, rel_floor, tol)
}

/// The share of `obs` not explained by the chance rate `chance`:
/// `(obs - chance) / (1 - chance)`, floored at 0.
fn above_chance(obs: f64, chance: f64) -> f64 {
    if chance >= 1.0 {
        return 0.0;
    }
    ((obs - chance) / (1.0 - chance)).max(0.0)
}

/// Fraction of `spectra` in which at least `min` of `channels` (moved by
/// `shift` Da) are present.
fn fraction_showing(
    spectra: &[&Spectrum],
    channels: &[f64],
    min: usize,
    shift: f64,
    tol: f64,
) -> f64 {
    if spectra.is_empty() {
        return 0.0;
    }
    let n = spectra
        .iter()
        .filter(|s| channels_present_at(s, channels, shift, REL_FLOOR, tol) >= min)
        .count();
    n as f64 / spectra.len() as f64
}

/// The `min`-th highest presence rate above chance among `channels`: each
/// channel's rate in `spectra`, corrected for the chance rate at its own
/// [`CONTROL_SHIFTS_DA`] positions. A reporter cluster raises every channel;
/// label-free ions that happen to sit at one or two reporter masses (an
/// immonium or internal ion near 127 or 129) raise only those.
fn channel_rate_above_chance(spectra: &[&Spectrum], channels: &[f64], min: usize, tol: f64) -> f64 {
    let mut rates: Vec<f64> = channels
        .iter()
        .map(|c| {
            let one = std::slice::from_ref(c);
            let obs = fraction_showing(spectra, one, 1, 0.0, tol);
            let chance = CONTROL_SHIFTS_DA
                .iter()
                .map(|&d| fraction_showing(spectra, one, 1, f64::from(d), tol))
                .sum::<f64>()
                / CONTROL_SHIFTS_DA.len() as f64;
            above_chance(obs, chance)
        })
        .collect();
    rates.sort_by(|a, b| b.total_cmp(a));
    rates.get(min.saturating_sub(1)).copied().unwrap_or(0.0)
}

/// [`fraction_showing`] at the reporter positions, corrected for the chance
/// rate measured at the [`CONTROL_SHIFTS_DA`] positions of the same spectra.
fn fraction_above_chance(spectra: &[&Spectrum], channels: &[f64], min: usize, tol: f64) -> f64 {
    let obs = fraction_showing(spectra, channels, min, 0.0, tol);
    let chance = CONTROL_SHIFTS_DA
        .iter()
        .map(|&d| fraction_showing(spectra, channels, min, f64::from(d), tol))
        .sum::<f64>()
        / CONTROL_SHIFTS_DA.len() as f64;
    above_chance(obs, chance)
}

/// Detect the isobaric label from a sample of MS2 spectra, or `None` if the run
/// is label-free or the sample is too small / ambiguous. `high_res` selects the
/// reporter match window (tight for Orbitrap/TOF, loose for ion-trap) — pass the
/// detected analyzer's resolution class.
///
/// Every fraction is taken above the run's chance rate (see the module docs).
/// TMT and iTRAQ reporter regions do not overlap (126–131 vs 114–117), so the
/// stronger cluster wins when (improbably) both clear the threshold.
pub fn detect_isobaric(spectra: &[Spectrum], high_res: bool) -> Option<IsobaricLabel> {
    let tol = if high_res { HIGH_RES_TOL } else { LOW_RES_TOL };
    let ms2: Vec<&Spectrum> = spectra.iter().filter(|s| !s.peaks.is_empty()).collect();
    if ms2.len() < MIN_SAMPLE {
        return None;
    }
    let n = ms2.len() as f64;

    let tmt_hits: Vec<&&Spectrum> = ms2
        .iter()
        .filter(|s| channels_present(s, TMT_CHANNELS, REL_FLOOR, tol) >= TMT_MIN_CHANNELS)
        .collect();
    let itraq_hits: Vec<&&Spectrum> = ms2
        .iter()
        .filter(|s| channels_present(s, ITRAQ_CHANNELS, REL_FLOOR, tol) >= ITRAQ_MIN_CHANNELS)
        .collect();

    // A cluster counts as far as both the spectra that show it and its
    // weakest required channel stand above chance.
    let cluster = |hits: usize, channels: &[f64], min: usize| {
        if (hits as f64) < FRACTION_THRESHOLD * n {
            return 0.0;
        }
        fraction_above_chance(&ms2, channels, min, tol)
            .min(channel_rate_above_chance(&ms2, channels, min, tol))
    };
    let tmt_frac = cluster(tmt_hits.len(), TMT_CHANNELS, TMT_MIN_CHANNELS);
    let itraq_frac = cluster(itraq_hits.len(), ITRAQ_CHANNELS, ITRAQ_MIN_CHANNELS);

    // Among the spectra that carry the core cluster, how many also carry the
    // larger kit's extra reporters, above the chance rate in those spectra.
    let shows_extra = |hits: &[&&Spectrum], extra: &[f64]| {
        let hits: Vec<&Spectrum> = hits.iter().map(|s| **s).collect();
        fraction_above_chance(&hits, extra, EXTRA_MIN_CHANNELS, tol) >= EXTRA_FRACTION
    };

    if tmt_frac >= FRACTION_THRESHOLD && tmt_frac >= itraq_frac {
        Some(if shows_extra(&tmt_hits, TMTPRO_EXTRA) {
            IsobaricLabel::TmtPro
        } else {
            IsobaricLabel::Tmt
        })
    } else if itraq_frac >= FRACTION_THRESHOLD {
        Some(if shows_extra(&itraq_hits, ITRAQ8_EXTRA) {
            IsobaricLabel::Itraq8
        } else {
            IsobaricLabel::Itraq
        })
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
            isolation_target_mz: None,
        }
    }

    /// Backbone peaks shared by all synthetic spectra (away from reporter region).
    fn backbone() -> Vec<(f64, f32)> {
        vec![
            (300.0, 500.0),
            (450.0, 800.0),
            (620.0, 1000.0),
            (780.0, 400.0),
        ]
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

    fn run_with(channels: &[f64]) -> Vec<Spectrum> {
        (0..200)
            .map(|_| {
                let mut p = backbone();
                for &c in channels {
                    p.push((c, 200.0));
                }
                p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                spec_with(p)
            })
            .collect()
    }

    #[test]
    fn detects_tmtpro_and_itraq8_from_their_extra_reporters() {
        let tmtpro: Vec<f64> = TMT_CHANNELS.iter().chain(TMTPRO_EXTRA).copied().collect();
        assert_eq!(
            detect_isobaric(&run_with(&tmtpro), true),
            Some(IsobaricLabel::TmtPro)
        );
        let itraq8: Vec<f64> = ITRAQ_CHANNELS.iter().chain(ITRAQ8_EXTRA).copied().collect();
        assert_eq!(
            detect_isobaric(&run_with(&itraq8), true),
            Some(IsobaricLabel::Itraq8)
        );
        // Low-res: nominal 132–134 peaks still read as TMTpro.
        let lowres: Vec<f64> = (126..=134).map(|c| c as f64 + 0.1).collect();
        assert_eq!(
            detect_isobaric(&run_with(&lowres), false),
            Some(IsobaricLabel::TmtPro)
        );
        // The 6–11-plex cluster alone stays TMT.
        assert_eq!(
            detect_isobaric(&run_with(TMT_CHANNELS), true),
            Some(IsobaricLabel::Tmt)
        );
    }

    /// A dense low-m/z region: a peak at each nominal mass 100–160 with
    /// probability `p` (deterministic pseudo-random), at the +0.1 Da offset of
    /// fragment ions, at 5% of the base peak.
    fn dense_lowres_background(seed: u64, p: f64) -> Vec<(f64, f32)> {
        let mut x = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut peaks = backbone();
        for m in 100..=160 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            if ((x >> 33) as f64 / (1u64 << 31) as f64) < p {
                peaks.push((m as f64 + 0.1, 50.0));
            }
        }
        peaks
    }

    /// Regression for a label-free ion-trap HCD run read as TMTpro: its
    /// low-m/z region holds a peak at most nominal masses, so at the low-res
    /// window the 4-of-6 cluster and the 132–134 extras appear by chance. The
    /// chance rate at the control positions explains them away. A TMTpro run
    /// over the same background is still detected, and a TMT 6–11 run over it
    /// is not promoted to TMTpro.
    #[test]
    fn dense_lowres_label_free_is_not_tmt() {
        let label_free: Vec<Spectrum> = (0..300)
            .map(|i| {
                let mut p = dense_lowres_background(i, 0.65);
                p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                spec_with(p)
            })
            .collect();
        let raw_hits = label_free
            .iter()
            .filter(|s| channels_present(s, TMT_CHANNELS, REL_FLOOR, LOW_RES_TOL) >= 4)
            .count();
        assert!(
            raw_hits as f64 >= FRACTION_THRESHOLD * label_free.len() as f64,
            "the background alone clears the raw fraction ({raw_hits})"
        );
        assert_eq!(detect_isobaric(&label_free, false), None);

        let with = |channels: &[f64]| -> Vec<Spectrum> {
            (0..300)
                .map(|i| {
                    let mut p = dense_lowres_background(i, 0.65);
                    for &c in channels {
                        p.push((c, 300.0));
                    }
                    p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                    spec_with(p)
                })
                .collect()
        };
        let tmtpro: Vec<f64> = TMT_CHANNELS.iter().chain(TMTPRO_EXTRA).copied().collect();
        assert_eq!(
            detect_isobaric(&with(&tmtpro), false),
            Some(IsobaricLabel::TmtPro)
        );
        assert_eq!(
            detect_isobaric(&with(TMT_CHANNELS), false),
            Some(IsobaricLabel::Tmt)
        );
    }

    /// The shape of the measured false positive: over a dense background,
    /// label-free ions at two reporter masses (nominal 127 and 129) lift the
    /// 4-of-6 cluster above the chance rate of the control positions. Two
    /// channels above chance are not a cluster.
    #[test]
    fn two_label_free_ions_at_reporter_masses_are_not_tmt() {
        let run: Vec<Spectrum> = (0..400)
            .map(|i| {
                let mut p = dense_lowres_background(i, 0.5);
                if i % 10 != 0 {
                    p.push((127.09, 80.0));
                    p.push((129.10, 80.0));
                }
                p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                spec_with(p)
            })
            .collect();
        let refs: Vec<&Spectrum> = run.iter().collect();
        assert!(
            fraction_above_chance(&refs, TMT_CHANNELS, TMT_MIN_CHANNELS, LOW_RES_TOL)
                >= FRACTION_THRESHOLD,
            "the spectrum-level rule alone clears the threshold"
        );
        assert_eq!(detect_isobaric(&run, false), None);
    }

    #[test]
    fn above_chance_is_the_unexplained_share() {
        assert_eq!(above_chance(0.8, 0.0), 0.8);
        assert!((above_chance(0.65, 0.6) - 0.125).abs() < 1e-12);
        assert_eq!(above_chance(0.5, 0.6), 0.0);
        assert_eq!(above_chance(1.0, 1.0), 0.0);
    }

    #[test]
    fn label_masses() {
        assert_eq!(IsobaricLabel::Tmt.tag().1, 229.162932);
        assert_eq!(IsobaricLabel::TmtPro.tag().1, 304.207146);
        assert_eq!(IsobaricLabel::Itraq.tag().1, 144.102063);
        assert_eq!(IsobaricLabel::Itraq8.tag().1, 304.205360);
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
        assert_eq!(
            detect_isobaric(&run, true),
            None,
            "high-res glyco must not be read as TMT"
        );
        // And with a loose (low-res) window this WOULD have false-triggered —
        // documents why high-res must use the tight window.
        assert_eq!(detect_isobaric(&run, false), Some(IsobaricLabel::Tmt));
    }
}
