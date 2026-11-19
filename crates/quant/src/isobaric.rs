pub use model::reporter::{Channel, DEUTERATED_TAG_SHIFT_DA};
use model::reporter::{ITRAQ4, ITRAQ8, TMT10, TMT11, TMT16, TMT18, TMT32, TMT35, TMT6};
use model::scan::NoiseProfile;
use model::Tolerance;

use crate::nnls::nnls;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Plex {
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
    Tmt18,
    Tmt32,
    Tmt35,
    Itraq4,
    Itraq8,
}

impl Plex {
    pub fn parse(s: &str) -> Result<Self, String> {
        let lower = s.trim().to_ascii_lowercase();
        let key = lower.strip_suffix("plex").unwrap_or(&lower);
        Ok(match key {
            "tmt6" => Plex::Tmt6,
            "tmt10" => Plex::Tmt10,
            "tmt11" => Plex::Tmt11,
            "tmt16" | "tmtpro16" | "tmtpro" => Plex::Tmt16,
            "tmt18" | "tmtpro18" => Plex::Tmt18,
            "tmt32" | "tmtpro32" => Plex::Tmt32,
            "tmt35" | "tmtpro35" => Plex::Tmt35,
            "itraq4" => Plex::Itraq4,
            "itraq8" => Plex::Itraq8,
            _ => {
                return Err(format!(
                    "unknown isobaric plex '{s}' (expected tmt6, tmt10, tmt11, tmt16, tmt18, \
                     tmt32, tmt35, itraq4 or itraq8)"
                ))
            }
        })
    }

    pub fn name(&self) -> &'static str {
        match self {
            Plex::Tmt6 => "tmt6plex",
            Plex::Tmt10 => "tmt10plex",
            Plex::Tmt11 => "tmt11plex",
            Plex::Tmt16 => "tmt16plex",
            Plex::Tmt18 => "tmt18plex",
            Plex::Tmt32 => "tmt32plex",
            Plex::Tmt35 => "tmt35plex",
            Plex::Itraq4 => "itraq4plex",
            Plex::Itraq8 => "itraq8plex",
        }
    }

    pub fn label_prefix(&self) -> &'static str {
        match self {
            Plex::Itraq4 | Plex::Itraq8 => "ITRAQ",
            _ => "TMT",
        }
    }

    pub fn channels(&self) -> &'static [Channel] {
        match self {
            Plex::Tmt6 => TMT6,
            Plex::Tmt10 => TMT10,
            Plex::Tmt11 => TMT11,
            Plex::Tmt16 => TMT16,
            Plex::Tmt18 => TMT18,
            Plex::Tmt32 => TMT32,
            Plex::Tmt35 => TMT35,
            Plex::Itraq4 => ITRAQ4,
            Plex::Itraq8 => ITRAQ8,
        }
    }

    pub fn n_channels(&self) -> usize {
        self.channels().len()
    }

    pub fn min_channel_spacing(&self) -> f64 {
        let mut mz: Vec<f64> = self.channels().iter().map(|c| c.mz).collect();
        mz.sort_by(f64::total_cmp);
        mz.windows(2)
            .map(|w| w[1] - w[0])
            .fold(f64::INFINITY, f64::min)
    }

    pub fn default_reporter_ppm(&self) -> f64 {
        if self.has_deuterated_channels() {
            10.0
        } else {
            20.0
        }
    }

    pub fn has_deuterated_channels(&self) -> bool {
        matches!(self, Plex::Tmt32 | Plex::Tmt35)
    }

    pub fn check_tolerance(&self, tol: Tolerance) -> Result<(), String> {
        let half_spacing = self.min_channel_spacing() / 2.0;
        let widest = self
            .channels()
            .iter()
            .map(|c| tol.as_da(c.mz))
            .fold(0.0f64, f64::max);
        if widest >= half_spacing {
            return Err(format!(
                "reporter tolerance {:.4} Da is not below half the {:.4} Da spacing between \
                 the closest {} channels, so neighbouring channels would read the same peak",
                widest,
                2.0 * half_spacing,
                self.name()
            ));
        }
        Ok(())
    }

    pub fn label_mass(&self) -> f64 {
        match self {
            Plex::Tmt6 | Plex::Tmt10 | Plex::Tmt11 => 229.162932,
            Plex::Tmt16 | Plex::Tmt18 | Plex::Tmt32 | Plex::Tmt35 => 304.207146,
            Plex::Itraq4 => 144.102063,
            Plex::Itraq8 => 304.205360,
        }
    }

    pub fn correction_columns(&self) -> Option<usize> {
        match self {
            Plex::Tmt6 | Plex::Tmt10 | Plex::Tmt11 | Plex::Itraq4 | Plex::Itraq8 => Some(4),
            Plex::Tmt16 | Plex::Tmt18 => Some(8),
            Plex::Tmt32 | Plex::Tmt35 => None,
        }
    }

    fn affected(&self, from: usize, col: usize) -> Option<usize> {
        let chans = self.channels();
        let (nominal, series) = parse_channel(chans[from].name, *self);
        let (delta, series_change): (i32, SeriesChange) = match self.correction_columns()? {
            4 => {
                let d = [-2, -1, 1, 2][col];
                (d, SeriesChange::Same)
            }
            _ => match col {
                0 => (-2, SeriesChange::Same),
                1 => (-2, SeriesChange::NToC),
                2 => (-1, SeriesChange::Same),
                3 => (-1, SeriesChange::NToC),
                4 => (1, SeriesChange::CToN),
                5 => (1, SeriesChange::Same),
                6 => (2, SeriesChange::CToN),
                _ => (2, SeriesChange::Same),
            },
        };
        let target_series = match (series, series_change) {
            (s, SeriesChange::Same) => s,
            (Series::N, SeriesChange::NToC) => Series::C,
            (Series::C, SeriesChange::CToN) => Series::N,
            (Series::Plain, SeriesChange::CToN) => Series::N,
            _ => return None,
        };
        let target_nominal = nominal as i32 + delta;
        if target_nominal <= 0 {
            return None;
        }
        chans.iter().position(|c| {
            let (n, s) = parse_channel(c.name, *self);
            n as i32 == target_nominal && series_compatible(s, target_series, *self)
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Series {
    Plain,
    N,
    C,
}

#[derive(Debug, Clone, Copy)]
enum SeriesChange {
    Same,
    NToC,
    CToN,
}

fn parse_channel(name: &str, plex: Plex) -> (u32, Series) {
    let (digits, suffix) = name.split_at(name.trim_end_matches(['N', 'C']).len());
    let nominal: u32 = digits.parse().unwrap_or(0);
    let series = match (suffix, plex) {
        ("N", _) => Series::N,
        ("C", _) => Series::C,
        ("", Plex::Tmt10 | Plex::Tmt11 | Plex::Tmt16 | Plex::Tmt18) => {
            if nominal == 126 {
                Series::C
            } else if nominal == 131 {
                Series::N
            } else {
                Series::Plain
            }
        }
        _ => Series::Plain,
    };
    (nominal, series)
}

fn series_compatible(have: Series, want: Series, plex: Plex) -> bool {
    match plex {
        Plex::Tmt6 | Plex::Itraq4 | Plex::Itraq8 => true,
        _ => have == want,
    }
}

pub fn extract_reporters(peaks: &[(f64, f32)], plex: Plex, tol: Tolerance) -> Vec<f32> {
    plex.channels()
        .iter()
        .map(|c| {
            let half = tol.as_da(c.mz);
            let lo = c.mz - half;
            let hi = c.mz + half;
            let start = peaks.partition_point(|&(mz, _)| mz < lo);
            let mut best = 0.0f32;
            for &(mz, inten) in &peaks[start..] {
                if mz > hi {
                    break;
                }
                if inten > best {
                    best = inten;
                }
            }
            best
        })
        .collect()
}

pub fn reporter_noise(noise: &NoiseProfile, plex: Plex) -> Vec<f32> {
    plex.channels().iter().map(|c| noise.at(c.mz)).collect()
}

pub fn signal_to_noise(values: &[f32], noise: &[f32]) -> Vec<f32> {
    values
        .iter()
        .zip(noise)
        .map(|(&v, &n)| if n > 0.0 { v.max(0.0) / n } else { 0.0 })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
pub struct CorrectionMatrix {
    plex: Plex,
    n: usize,
    m: Vec<f64>,
}

impl CorrectionMatrix {
    pub fn parse(plex: Plex, text: &str) -> Result<Self, String> {
        let n = plex.n_channels();
        let cols = Self::check_supported(plex)?;
        let rows: Vec<&str> = text
            .lines()
            .map(|l| l.split('#').next().unwrap_or("").trim())
            .filter(|l| !l.is_empty())
            .collect();
        if rows.len() != n {
            return Err(format!(
                "impurity matrix for {} needs {n} channel rows, found {}",
                plex.name(),
                rows.len()
            ));
        }
        let mut m = vec![0.0f64; n * n];
        for (from, row) in rows.iter().enumerate() {
            let (label, body) = match row.find([':', ',']) {
                Some(p) if !row[..p].contains('/') => (Some(row[..p].trim()), row[p + 1..].trim()),
                _ => match row.split_once(char::is_whitespace) {
                    Some((head, rest))
                        if !head.contains('/') && !rest.trim_start().starts_with('/') =>
                    {
                        (Some(head), rest.trim())
                    }
                    _ => (None, *row),
                },
            };
            if let Some(label) = label {
                let expected = plex.channels()[from].name;
                if label != expected {
                    return Err(format!(
                        "impurity row {} is labelled '{label}', expected '{expected}' in kit order",
                        from + 1
                    ));
                }
            }
            let values: Vec<&str> = body.split('/').map(str::trim).collect();
            if values.len() != cols {
                return Err(format!(
                    "impurity row {} ({}) must have {cols} '/'-separated entries, found {}",
                    from + 1,
                    plex.channels()[from].name,
                    values.len()
                ));
            }
            let mut self_share = 100.0f64;
            for (col, v) in values.iter().enumerate() {
                let upper = v.to_ascii_uppercase();
                if upper.is_empty() || upper == "NA" || upper == "-1" {
                    continue;
                }
                let pct: f64 = upper.parse().map_err(|_| {
                    format!(
                        "impurity row {} entry {}: '{v}' is not a percentage",
                        from + 1,
                        col + 1
                    )
                })?;
                if !(0.0..=100.0).contains(&pct) {
                    return Err(format!(
                        "impurity row {} entry {}: {pct} is outside 0..100",
                        from + 1,
                        col + 1
                    ));
                }
                if let Some(target) = plex.affected(from, col) {
                    m[target * n + from] = pct / 100.0;
                }
                self_share -= pct;
            }
            if self_share < 0.0 {
                return Err(format!(
                    "impurity row {} ({}) percentages sum to more than 100",
                    from + 1,
                    plex.channels()[from].name
                ));
            }
            m[from * n + from] = self_share / 100.0;
        }
        Ok(Self { plex, n, m })
    }

    pub fn check_supported(plex: Plex) -> Result<usize, String> {
        plex.correction_columns().ok_or_else(|| {
            format!(
                "isotope impurity correction is not supported for {} (the reference \
                 isobaric implementation applies none for this kit either)",
                plex.name()
            )
        })
    }

    pub fn plex(&self) -> Plex {
        self.plex
    }

    pub fn is_identity(&self) -> bool {
        (0..self.n).all(|i| {
            (0..self.n).all(|j| {
                let v = self.m[i * self.n + j];
                if i == j {
                    (v - 1.0).abs() < 1e-12
                } else {
                    v.abs() < 1e-12
                }
            })
        })
    }

    pub fn get(&self, target: usize, source: usize) -> f64 {
        self.m[target * self.n + source]
    }

    pub fn correct(&self, observed: &[f32]) -> Vec<f32> {
        if observed.len() != self.n || observed.iter().all(|&v| v <= 0.0) {
            return observed.to_vec();
        }
        let b: Vec<f64> = observed.iter().map(|&v| v.max(0.0) as f64).collect();
        let x = nnls(&self.m, self.n, self.n, &b, 6 * self.n + 10);
        x.into_iter().map(|v| v as f32).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL: [Plex; 9] = [
        Plex::Tmt6,
        Plex::Tmt10,
        Plex::Tmt11,
        Plex::Tmt16,
        Plex::Tmt18,
        Plex::Tmt32,
        Plex::Tmt35,
        Plex::Itraq4,
        Plex::Itraq8,
    ];

    #[test]
    fn plex_parse_accepts_common_spellings() {
        assert_eq!(Plex::parse("tmt10").unwrap(), Plex::Tmt10);
        assert_eq!(Plex::parse("TMT11plex").unwrap(), Plex::Tmt11);
        assert_eq!(Plex::parse("tmtpro").unwrap(), Plex::Tmt16);
        assert_eq!(Plex::parse("iTRAQ8").unwrap(), Plex::Itraq8);
        assert_eq!(Plex::parse("tmt32").unwrap(), Plex::Tmt32);
        assert_eq!(Plex::parse("TMT35plex").unwrap(), Plex::Tmt35);
        assert_eq!(Plex::parse("tmtpro35").unwrap(), Plex::Tmt35);
        assert!(Plex::parse("tmt12").is_err());
        assert!(Plex::parse("tmt36").is_err());
        for plex in ALL {
            assert_eq!(Plex::parse(plex.name()).unwrap(), plex);
        }
    }

    #[test]
    fn channel_tables_are_sorted_and_sized() {
        for plex in ALL {
            let c = plex.channels();
            if matches!(plex, Plex::Tmt32 | Plex::Tmt35) {
                assert!(
                    c.windows(2).all(|w| w[0].mz.floor() <= w[1].mz.floor()),
                    "{plex:?}"
                );
            } else {
                assert!(c.windows(2).all(|w| w[0].mz < w[1].mz), "{plex:?}");
            }
            assert!(plex.min_channel_spacing() > 0.0029, "{plex:?}");
            let names: std::collections::HashSet<_> = c.iter().map(|c| c.name).collect();
            assert_eq!(names.len(), c.len(), "{plex:?}");
        }
        assert_eq!(Plex::Tmt10.n_channels(), 10);
        assert_eq!(Plex::Tmt18.n_channels(), 18);
        assert_eq!(Plex::Tmt32.n_channels(), 32);
        assert_eq!(Plex::Tmt35.n_channels(), 35);
        assert_eq!(Plex::Itraq8.n_channels(), 8);
    }

    #[test]
    fn deuterated_kits_extend_tmtpro() {
        let in_35 = |c: &Channel| {
            Plex::Tmt35
                .channels()
                .iter()
                .any(|d| d.name == c.name && (d.mz - c.mz).abs() < 2e-6)
        };
        assert!(Plex::Tmt18.channels().iter().all(in_35));
        assert!(Plex::Tmt32.channels().iter().all(in_35));
        let mz = |name: &str| {
            Plex::Tmt35
                .channels()
                .iter()
                .find(|c| c.name == name)
                .unwrap()
                .mz
        };
        for (d, twin) in [
            ("127D", "127C"),
            ("128ND", "128N"),
            ("128CD", "128C"),
            ("133CD", "133C"),
            ("134ND", "134N"),
            ("135ND", "135N"),
        ] {
            assert!(
                (mz(d) - mz(twin) - DEUTERATED_TAG_SHIFT_DA).abs() < 5e-6,
                "{d}"
            );
        }
        assert_eq!(Plex::Tmt35.label_mass(), Plex::Tmt18.label_mass());
    }

    #[test]
    fn extraction_takes_most_intense_peak_in_window() {
        let peaks = vec![
            (126.1270, 50.0f32),
            (126.1278, 900.0),
            (127.1311, 300.0),
            (131.1382, 10.0),
        ];
        let r = extract_reporters(&peaks, Plex::Tmt10, Tolerance::Ppm(20.0));
        assert_eq!(r.len(), 10);
        assert_eq!(r[0], 900.0);
        assert_eq!(r[1], 0.0);
        assert_eq!(r[2], 300.0);
        assert_eq!(r[9], 10.0);
        let r2 = extract_reporters(&peaks, Plex::Tmt10, Tolerance::Da(0.3));
        assert_eq!(r2[1], 300.0);
        assert_eq!(r2[0], 900.0);
    }

    #[test]
    fn signal_to_noise_divides_by_the_noise_at_each_reporter() {
        let profile =
            NoiseProfile::reporter_region([(125.0, 100.0), (127.0, 300.0), (132.0, 300.0)])
                .unwrap();
        let noise = reporter_noise(&profile, Plex::Tmt6);
        assert_eq!(noise.len(), 6);
        assert!(
            (noise[0] - (100.0 + 200.0 * 1.127725 / 2.0)).abs() < 0.1,
            "{noise:?}"
        );
        assert_eq!(noise[5], 300.0);
        let sn = signal_to_noise(&[3000.0, 0.0, 600.0, -1.0, 300.0, 900.0], &noise);
        assert!((sn[0] - 3000.0 / noise[0]).abs() < 1e-3);
        assert_eq!(sn[1], 0.0);
        assert_eq!(sn[3], 0.0);
        assert_eq!(sn[5], 3.0);
    }

    #[test]
    fn tolerance_must_separate_the_closest_channels() {
        assert!((Plex::Tmt10.min_channel_spacing() - 0.00632).abs() < 1e-5);
        assert!(Plex::Tmt10.check_tolerance(Tolerance::Ppm(20.0)).is_ok());
        assert!(Plex::Tmt18.check_tolerance(Tolerance::Ppm(20.0)).is_ok());
        assert!(Plex::Tmt10.check_tolerance(Tolerance::Da(0.3)).is_err());
        assert!(Plex::Tmt16.check_tolerance(Tolerance::Da(0.005)).is_err());
        assert!(Plex::Tmt6.check_tolerance(Tolerance::Da(0.3)).is_ok());
        assert!(Plex::Itraq8.check_tolerance(Tolerance::Da(0.3)).is_ok());
        for plex in [Plex::Tmt32, Plex::Tmt35] {
            assert!(plex.has_deuterated_channels());
            assert!((plex.min_channel_spacing() - DEUTERATED_TAG_SHIFT_DA).abs() < 5e-6);
            assert!(plex.check_tolerance(Tolerance::Ppm(20.0)).is_err());
            assert!(plex.check_tolerance(Tolerance::Ppm(10.0)).is_ok());
        }
        for plex in ALL {
            let tol = Tolerance::Ppm(plex.default_reporter_ppm());
            assert!(plex.check_tolerance(tol).is_ok(), "{plex:?}");
        }
    }

    #[test]
    fn deuterated_channels_are_read_separately() {
        let peaks = vec![
            (128.128116, 100.0f32),
            (128.131038, 200.0),
            (128.134436, 300.0),
            (128.137358, 400.0),
        ];
        let r = extract_reporters(&peaks, Plex::Tmt35, Tolerance::Ppm(10.0));
        let idx = |name: &str| {
            Plex::Tmt35
                .channels()
                .iter()
                .position(|c| c.name == name)
                .unwrap()
        };
        assert_eq!(r[idx("128N")], 100.0);
        assert_eq!(r[idx("128ND")], 200.0);
        assert_eq!(r[idx("128C")], 300.0);
        assert_eq!(r[idx("128CD")], 400.0);
        assert_eq!(r[idx("127D")], 0.0);
    }

    #[test]
    fn tmt10_affected_channels_follow_the_c13_series() {
        let p = Plex::Tmt10;
        let idx = |name: &str| p.channels().iter().position(|c| c.name == name).unwrap();
        assert_eq!(p.affected(idx("126"), 2), Some(idx("127C")));
        assert_eq!(p.affected(idx("126"), 3), Some(idx("128C")));
        assert_eq!(p.affected(idx("126"), 1), None);
        assert_eq!(p.affected(idx("127N"), 2), Some(idx("128N")));
        assert_eq!(p.affected(idx("128N"), 1), Some(idx("127N")));
        assert_eq!(p.affected(idx("130N"), 2), Some(idx("131")));
        assert_eq!(p.affected(idx("130C"), 2), None);
        assert_eq!(p.affected(idx("127C"), 1), Some(idx("126")));
    }

    #[test]
    fn tmt6_and_itraq_use_plain_nominal_neighbours() {
        let p = Plex::Tmt6;
        assert_eq!(p.affected(0, 2), Some(1));
        assert_eq!(p.affected(5, 2), None);
        let q = Plex::Itraq8;
        let idx = |name: &str| q.channels().iter().position(|c| c.name == name).unwrap();
        assert_eq!(q.affected(idx("119"), 2), None);
        assert_eq!(q.affected(idx("119"), 3), Some(idx("121")));
    }

    #[test]
    fn tmtpro_eight_column_mapping() {
        let p = Plex::Tmt16;
        let idx = |name: &str| p.channels().iter().position(|c| c.name == name).unwrap();
        assert_eq!(p.affected(idx("127C"), 4), Some(idx("128N")));
        assert_eq!(p.affected(idx("127C"), 5), Some(idx("128C")));
        assert_eq!(p.affected(idx("127C"), 6), Some(idx("129N")));
        assert_eq!(p.affected(idx("127C"), 2), Some(idx("126")));
        assert_eq!(p.affected(idx("127C"), 3), None);
        assert_eq!(p.affected(idx("128N"), 3), Some(idx("127C")));
        assert_eq!(p.affected(idx("128N"), 1), Some(idx("126")));
    }

    #[test]
    fn parses_text_matrix_and_corrects() {
        let text = "126:0.0/0.0/5.0/0.0\n127:0.0/0.0/5.0/0.0\n128:0.0/0.0/5.0/0.0\n\
                    129:0.0/0.0/5.0/0.0\n130:0.0/0.0/5.0/0.0\n131:0.0/0.0/0.0/0.0\n";
        let m = CorrectionMatrix::parse(Plex::Tmt6, text).unwrap();
        assert!(!m.is_identity());
        assert!((m.get(0, 0) - 0.95).abs() < 1e-12);
        assert!((m.get(1, 0) - 0.05).abs() < 1e-12);
        assert!((m.get(5, 5) - 1.0).abs() < 1e-12);
        let truth = [1000.0f32, 2000.0, 500.0, 0.0, 100.0, 300.0];
        let observed: Vec<f32> = (0..6)
            .map(|t| (0..6).map(|s| m.get(t, s) as f32 * truth[s]).sum())
            .collect();
        let corrected = m.correct(&observed);
        for (c, t) in corrected.iter().zip(truth.iter()) {
            assert!((c - t).abs() < 1e-2, "{corrected:?}");
        }
    }

    #[test]
    fn deuterated_kits_refuse_an_impurity_sheet() {
        for plex in [Plex::Tmt32, Plex::Tmt35] {
            let na = "NA/NA/NA/NA/NA/NA/NA/NA\n".repeat(plex.n_channels());
            let error = CorrectionMatrix::parse(plex, &na).unwrap_err();
            assert!(error.contains("not supported"), "{error}");
        }
    }

    #[test]
    fn matrix_row_count_and_width_are_checked() {
        assert!(CorrectionMatrix::parse(Plex::Tmt6, "0/0/1/0\n").is_err());
        let bad_width = "0/0/1\n".repeat(6);
        assert!(CorrectionMatrix::parse(Plex::Tmt6, &bad_width).is_err());
        let na = "NA/NA/NA/NA\n".repeat(4);
        let m = CorrectionMatrix::parse(Plex::Itraq4, &na).unwrap();
        assert!(m.is_identity());
    }

    #[test]
    fn correction_rejects_swapped_labels_and_impossible_spill() {
        let swapped = "115:0/0/5/0\n114:0/0/1/0\n116:0/0/0/0\n117:0/0/0/0\n";
        let error = CorrectionMatrix::parse(Plex::Itraq4, swapped).unwrap_err();
        assert!(error.contains("expected '114'"), "{error}");
        let impossible = "0/0/60/50\n0/0/0/0\n0/0/0/0\n0/0/0/0\n";
        let error = CorrectionMatrix::parse(Plex::Itraq4, impossible).unwrap_err();
        assert!(error.contains("more than 100"), "{error}");
    }

    #[test]
    fn unlabelled_percentages_can_have_spaces_around_separators() {
        let spaced = CorrectionMatrix::parse(Plex::Itraq4, &"0 / 0 / 5 / 0\n".repeat(4)).unwrap();
        let compact = CorrectionMatrix::parse(Plex::Itraq4, &"0/0/5/0\n".repeat(4)).unwrap();
        assert_eq!(spaced, compact);
        let labelled = "114 0 / 0 / 5 / 0\n115,0/0/5/0\n116:0/0/5/0\n117:0/0/5/0\n";
        assert_eq!(
            CorrectionMatrix::parse(Plex::Itraq4, labelled).unwrap(),
            compact
        );
    }
}
