use model::reporter;
use model::spectrum::Spectrum;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsobaricLabel {
    Tmt,
    TmtPro,
    Itraq,
    Itraq8,
}

impl IsobaricLabel {
    pub fn tag(self) -> (&'static str, f64) {
        match self {
            IsobaricLabel::Tmt => ("TMT6plex", 229.162932),
            IsobaricLabel::TmtPro => ("TMTpro", 304.207146),
            IsobaricLabel::Itraq => ("iTRAQ4plex", 144.102063),
            IsobaricLabel::Itraq8 => ("iTRAQ8plex", 304.205360),
        }
    }
}

const TMT_CHANNELS: &[f64] = &[
    reporter::TMT6[0].mz,
    reporter::TMT6[1].mz,
    reporter::TMT6[2].mz,
    reporter::TMT6[3].mz,
    reporter::TMT6[4].mz,
    reporter::TMT6[5].mz,
];
const ITRAQ_CHANNELS: &[f64] = &[
    reporter::ITRAQ4[0].mz,
    reporter::ITRAQ4[1].mz,
    reporter::ITRAQ4[2].mz,
    reporter::ITRAQ4[3].mz,
];
const TMTPRO_EXTRA: &[f64] = &[
    reporter::TMT16[11].mz,
    reporter::TMT16[12].mz,
    reporter::TMT16[13].mz,
    reporter::TMT16[14].mz,
    reporter::TMT16[15].mz,
];
const ITRAQ8_EXTRA: &[f64] = &[
    reporter::ITRAQ8[0].mz,
    reporter::ITRAQ8[5].mz,
    reporter::ITRAQ8[6].mz,
    reporter::ITRAQ8[7].mz,
];
const EXTRA_MIN_CHANNELS: usize = 2;
const EXTRA_FRACTION: f64 = 0.5;

const HIGH_RES_TOL: f64 = 0.01;
const LOW_RES_TOL: f64 = 0.3;

const CONTROL_SHIFTS_DA: [i32; 24] = [
    -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    18, 19, 20,
];

pub const MIN_SAMPLE: usize = 50;
const REL_FLOOR: f32 = 0.01;
const TMT_MIN_CHANNELS: usize = 4;
const ITRAQ_MIN_CHANNELS: usize = 3;
const FRACTION_THRESHOLD: f64 = 0.30;

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

fn above_chance(obs: f64, chance: f64) -> f64 {
    if chance >= 1.0 {
        return 0.0;
    }
    ((obs - chance) / (1.0 - chance)).max(0.0)
}

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

fn fraction_above_chance(spectra: &[&Spectrum], channels: &[f64], min: usize, tol: f64) -> f64 {
    let obs = fraction_showing(spectra, channels, min, 0.0, tol);
    let chance = CONTROL_SHIFTS_DA
        .iter()
        .map(|&d| fraction_showing(spectra, channels, min, f64::from(d), tol))
        .sum::<f64>()
        / CONTROL_SHIFTS_DA.len() as f64;
    above_chance(obs, chance)
}

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

    let cluster = |hits: usize, channels: &[f64], min: usize| {
        if (hits as f64) < FRACTION_THRESHOLD * n {
            return 0.0;
        }
        fraction_above_chance(&ms2, channels, min, tol)
            .min(channel_rate_above_chance(&ms2, channels, min, tol))
    };
    let tmt_frac = cluster(tmt_hits.len(), TMT_CHANNELS, TMT_MIN_CHANNELS);
    let itraq_frac = cluster(itraq_hits.len(), ITRAQ_CHANNELS, ITRAQ_MIN_CHANNELS);

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
            mobility_window: None,
        }
    }

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
        let lowres: Vec<f64> = (126..=134).map(|c| c as f64 + 0.1).collect();
        assert_eq!(
            detect_isobaric(&run_with(&lowres), false),
            Some(IsobaricLabel::TmtPro)
        );
        assert_eq!(
            detect_isobaric(&run_with(TMT_CHANNELS), true),
            Some(IsobaricLabel::Tmt)
        );
    }

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
        let run: Vec<Spectrum> = (0..200)
            .map(|_| {
                let mut p = backbone();
                for &c in TMT_CHANNELS {
                    p.push((c, 0.5));
                }
                p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                spec_with(p)
            })
            .collect();
        assert_eq!(detect_isobaric(&run, true), None);
    }

    #[test]
    fn dense_highres_glyco_lowmz_is_not_tmt() {
        let run: Vec<Spectrum> = (0..300)
            .map(|_| {
                let mut p = backbone();
                for n in 126..=131 {
                    p.push((n as f64 + 0.05, 600.0));
                }
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
        assert_eq!(detect_isobaric(&run, false), Some(IsobaricLabel::Tmt));
    }
}
