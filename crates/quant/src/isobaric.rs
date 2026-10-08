//! Isobaric reporter-ion quantification: plex tables, extraction and
//! isotope-impurity correction.
//!
//! Reporter m/z values are the ones OpenMS's `IsobaricQuantitationMethod`
//! subclasses carry (TMT 6/10/11/16/18, iTRAQ 4/8). Extraction takes the most
//! intense centroid within a tolerance of each reporter (Sage's rule; OpenMS
//! takes the nearest within `reporter_mass_shift`). Impurity correction reads
//! the vendor lot sheet in the OpenMS text form (one line per channel,
//! `-2/-1/+1/+2` percentages, eight columns for TMTpro) and solves the
//! observed = M · true system by non-negative least squares, as
//! `IsobaricIsotopeCorrector` does.

use model::Tolerance;

use crate::nnls::nnls;

/// One reporter channel: its name (as OpenMS labels it) and exact m/z.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Channel {
    pub name: &'static str,
    pub mz: f64,
}

const fn ch(name: &'static str, mz: f64) -> Channel {
    Channel { name, mz }
}

const TMT6: &[Channel] = &[
    ch("126", 126.127725),
    ch("127", 127.124760),
    ch("128", 128.134433),
    ch("129", 129.131468),
    ch("130", 130.141141),
    ch("131", 131.138176),
];

const TMT10: &[Channel] = &[
    ch("126", 126.127726),
    ch("127N", 127.124761),
    ch("127C", 127.131081),
    ch("128N", 128.128116),
    ch("128C", 128.134436),
    ch("129N", 129.131471),
    ch("129C", 129.137790),
    ch("130N", 130.134825),
    ch("130C", 130.141145),
    ch("131", 131.138180),
];

const TMT11: &[Channel] = &[
    ch("126", 126.127726),
    ch("127N", 127.124761),
    ch("127C", 127.131081),
    ch("128N", 128.128116),
    ch("128C", 128.134436),
    ch("129N", 129.131471),
    ch("129C", 129.137790),
    ch("130N", 130.134825),
    ch("130C", 130.141145),
    ch("131N", 131.138180),
    ch("131C", 131.144500),
];

const TMT16: &[Channel] = &[
    ch("126", 126.127726),
    ch("127N", 127.124761),
    ch("127C", 127.131081),
    ch("128N", 128.128116),
    ch("128C", 128.134436),
    ch("129N", 129.131471),
    ch("129C", 129.137790),
    ch("130N", 130.134825),
    ch("130C", 130.141145),
    ch("131N", 131.138180),
    ch("131C", 131.144500),
    ch("132N", 132.141535),
    ch("132C", 132.147855),
    ch("133N", 133.144890),
    ch("133C", 133.151210),
    ch("134N", 134.148245),
];

const TMT18: &[Channel] = &[
    ch("126", 126.127726),
    ch("127N", 127.124761),
    ch("127C", 127.131081),
    ch("128N", 128.128116),
    ch("128C", 128.134436),
    ch("129N", 129.131471),
    ch("129C", 129.137790),
    ch("130N", 130.134825),
    ch("130C", 130.141145),
    ch("131N", 131.138180),
    ch("131C", 131.144500),
    ch("132N", 132.141535),
    ch("132C", 132.147855),
    ch("133N", 133.144890),
    ch("133C", 133.151210),
    ch("134N", 134.148245),
    ch("134C", 134.154565),
    ch("135N", 135.151600),
];

const ITRAQ4: &[Channel] = &[
    ch("114", 114.1112),
    ch("115", 115.1082),
    ch("116", 116.1116),
    ch("117", 117.1149),
];

const ITRAQ8: &[Channel] = &[
    ch("113", 113.1078),
    ch("114", 114.1112),
    ch("115", 115.1082),
    ch("116", 116.1116),
    ch("117", 117.1149),
    ch("118", 118.1120),
    ch("119", 119.1153),
    ch("121", 121.1220),
];

/// Supported isobaric labelling kits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Plex {
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
    Tmt18,
    Itraq4,
    Itraq8,
}

impl Plex {
    /// Parse a CLI name: `tmt6`, `tmt10`, `tmt11`, `tmt16`, `tmt18`, `itraq4`,
    /// `itraq8` (case-insensitive, `plex` suffix optional, e.g. `tmt10plex`).
    pub fn parse(s: &str) -> Result<Self, String> {
        let lower = s.trim().to_ascii_lowercase();
        let key = lower.strip_suffix("plex").unwrap_or(&lower);
        Ok(match key {
            "tmt6" => Plex::Tmt6,
            "tmt10" => Plex::Tmt10,
            "tmt11" => Plex::Tmt11,
            "tmt16" | "tmtpro16" | "tmtpro" => Plex::Tmt16,
            "tmt18" | "tmtpro18" => Plex::Tmt18,
            "itraq4" => Plex::Itraq4,
            "itraq8" => Plex::Itraq8,
            _ => {
                return Err(format!(
                    "unknown isobaric plex '{s}' (expected tmt6, tmt10, tmt11, tmt16, tmt18, \
                     itraq4 or itraq8)"
                ))
            }
        })
    }

    /// The kit's name in OpenMS form (`tmt10plex`, `itraq4plex`, ...).
    pub fn name(&self) -> &'static str {
        match self {
            Plex::Tmt6 => "tmt6plex",
            Plex::Tmt10 => "tmt10plex",
            Plex::Tmt11 => "tmt11plex",
            Plex::Tmt16 => "tmt16plex",
            Plex::Tmt18 => "tmt18plex",
            Plex::Itraq4 => "itraq4plex",
            Plex::Itraq8 => "itraq8plex",
        }
    }

    /// The QPX / mzTab label prefix for this kit's channels: `TMT126`,
    /// `iTRAQ114`, ...
    pub fn label_prefix(&self) -> &'static str {
        match self {
            Plex::Itraq4 | Plex::Itraq8 => "iTRAQ",
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
            Plex::Itraq4 => ITRAQ4,
            Plex::Itraq8 => ITRAQ8,
        }
    }

    pub fn n_channels(&self) -> usize {
        self.channels().len()
    }

    /// The label mass added to K and the peptide N-terminus by this kit
    /// (Unimod monoisotopic delta), for the run log and the default fixed mod.
    pub fn label_mass(&self) -> f64 {
        match self {
            Plex::Tmt6 | Plex::Tmt10 | Plex::Tmt11 => 229.162932,
            Plex::Tmt16 | Plex::Tmt18 => 304.207146,
            Plex::Itraq4 => 144.102063,
            Plex::Itraq8 => 304.205360,
        }
    }

    /// Number of impurity columns per channel in the vendor sheet: four
    /// (`-2/-1/+1/+2`) for TMT ≤ 11-plex and iTRAQ, eight for TMTpro
    /// (`-2C13/-N15-C13/-C13/-N15/+N15/+C13/+N15+C13/+2C13`).
    pub fn correction_columns(&self) -> usize {
        match self {
            Plex::Tmt16 | Plex::Tmt18 => 8,
            _ => 4,
        }
    }

    /// Index of the channel that column `col` of channel `from`'s impurity
    /// row spills into, or `None` when that isotopologue is not a channel of
    /// the kit (the spill is still subtracted from the channel's own share).
    fn affected(&self, from: usize, col: usize) -> Option<usize> {
        let chans = self.channels();
        let (nominal, series) = parse_channel(chans[from].name, *self);
        let (delta, series_change): (i32, SeriesChange) = match self.correction_columns() {
            4 => {
                let d = [-2, -1, 1, 2][col];
                (d, SeriesChange::Same)
            }
            _ => match col {
                0 => (-2, SeriesChange::Same), // -2C13
                1 => (-2, SeriesChange::NToC), // -N15-C13
                2 => (-1, SeriesChange::Same), // -C13
                3 => (-1, SeriesChange::NToC), // -N15
                4 => (1, SeriesChange::CToN),  // +N15
                5 => (1, SeriesChange::Same),  // +C13
                6 => (2, SeriesChange::CToN),  // +N15+C13
                _ => (2, SeriesChange::Same),  // +2C13
            },
        };
        let target_series = match (series, series_change) {
            (s, SeriesChange::Same) => s,
            (Series::N, SeriesChange::NToC) => Series::C,
            (Series::C, SeriesChange::CToN) => Series::N,
            (Series::Plain, SeriesChange::CToN) => Series::N, // 126 behaves as C
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

/// Nominal mass and isotope series of a channel name. TMT 10/11-plex's `126`
/// carries no 15N and sits in the 13C series; its `131` (TMT10) is `131N`.
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

/// Extract the reporter intensities of `plex` from an m/z-sorted centroid
/// list: the most intense peak within `tol` of each channel m/z, 0.0 when no
/// peak lies in the window. The result has one entry per channel, in kit order.
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

/// Isotope impurity matrix `M` with `observed = M · true`: `m[target][source]`
/// is the fraction of channel `source`'s signal that is read in channel
/// `target`; the diagonal is what remains after every listed spill.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrectionMatrix {
    plex: Plex,
    n: usize,
    /// Row-major `n × n`.
    m: Vec<f64>,
}

impl CorrectionMatrix {
    /// Parse the OpenMS / vendor sheet text: one line per channel in kit
    /// order, `-2/-1/+1/+2` percentages (eight entries for TMTpro) separated
    /// by `/`. An optional `<channel>:` or `<channel>,` prefix names the row;
    /// `#` starts a comment; `NA`, `-1` and blank entries mean no spill.
    ///
    /// ```text
    /// # TMT10plex lot XY
    /// 126:0.0/0.0/5.09/0.0
    /// 127N:0.0/0.17/5.57/0.0
    /// ...
    /// ```
    pub fn parse(plex: Plex, text: &str) -> Result<Self, String> {
        let n = plex.n_channels();
        let cols = plex.correction_columns();
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
            // Strip a leading channel label ("126:" / "127N," / "127N ").
            let body = match row.find([':', ',']) {
                Some(p) if !row[..p].contains('/') => row[p + 1..].trim(),
                _ => match row.split_once(char::is_whitespace) {
                    Some((head, rest)) if !head.contains('/') => rest.trim(),
                    _ => row,
                },
            };
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
            m[from * n + from] = self_share / 100.0;
        }
        Ok(Self { plex, n, m })
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

    /// Entry `(target, source)`.
    pub fn get(&self, target: usize, source: usize) -> f64 {
        self.m[target * self.n + source]
    }

    /// Corrected (true) channel intensities for one spectrum's observed
    /// reporter vector, by NNLS. All-zero input stays all-zero.
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

    #[test]
    fn plex_parse_accepts_common_spellings() {
        assert_eq!(Plex::parse("tmt10").unwrap(), Plex::Tmt10);
        assert_eq!(Plex::parse("TMT11plex").unwrap(), Plex::Tmt11);
        assert_eq!(Plex::parse("tmtpro").unwrap(), Plex::Tmt16);
        assert_eq!(Plex::parse("iTRAQ8").unwrap(), Plex::Itraq8);
        assert!(Plex::parse("tmt12").is_err());
    }

    #[test]
    fn channel_tables_are_sorted_and_sized() {
        for plex in [
            Plex::Tmt6,
            Plex::Tmt10,
            Plex::Tmt11,
            Plex::Tmt16,
            Plex::Tmt18,
            Plex::Itraq4,
            Plex::Itraq8,
        ] {
            let c = plex.channels();
            assert!(c.windows(2).all(|w| w[0].mz < w[1].mz), "{plex:?}");
        }
        assert_eq!(Plex::Tmt10.n_channels(), 10);
        assert_eq!(Plex::Tmt18.n_channels(), 18);
        assert_eq!(Plex::Itraq8.n_channels(), 8);
    }

    #[test]
    fn extraction_takes_most_intense_peak_in_window() {
        // 126 has two peaks inside 20 ppm; 127N none; 127C exactly on target.
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
        // A 0.3 Da ion-trap window around 127N (126.82..127.42) picks up the 127C peak.
        let r2 = extract_reporters(&peaks, Plex::Tmt10, Tolerance::Da(0.3));
        assert_eq!(r2[1], 300.0);
        assert_eq!(r2[0], 900.0);
    }

    #[test]
    fn tmt10_affected_channels_follow_the_c13_series() {
        let p = Plex::Tmt10;
        let idx = |name: &str| p.channels().iter().position(|c| c.name == name).unwrap();
        // 126 +1 → 127C, +2 → 128C; 126 −1 → none.
        assert_eq!(p.affected(idx("126"), 2), Some(idx("127C")));
        assert_eq!(p.affected(idx("126"), 3), Some(idx("128C")));
        assert_eq!(p.affected(idx("126"), 1), None);
        // 127N +1 → 128N; 128N −1 → 127N.
        assert_eq!(p.affected(idx("127N"), 2), Some(idx("128N")));
        assert_eq!(p.affected(idx("128N"), 1), Some(idx("127N")));
        // 130N +1 → 131 (which is the N series in TMT10); 130C +1 → none.
        assert_eq!(p.affected(idx("130N"), 2), Some(idx("131")));
        assert_eq!(p.affected(idx("130C"), 2), None);
        // 127C −1 → 126.
        assert_eq!(p.affected(idx("127C"), 1), Some(idx("126")));
    }

    #[test]
    fn tmt6_and_itraq_use_plain_nominal_neighbours() {
        let p = Plex::Tmt6;
        assert_eq!(p.affected(0, 2), Some(1)); // 126 +1 → 127
        assert_eq!(p.affected(5, 2), None); // 131 +1 → nothing
        let q = Plex::Itraq8;
        let idx = |name: &str| q.channels().iter().position(|c| c.name == name).unwrap();
        assert_eq!(q.affected(idx("119"), 2), None); // 120 is not a channel
        assert_eq!(q.affected(idx("119"), 3), Some(idx("121")));
    }

    #[test]
    fn tmtpro_eight_column_mapping() {
        let p = Plex::Tmt16;
        let idx = |name: &str| p.channels().iter().position(|c| c.name == name).unwrap();
        // 127C: +N15 → 128N, +C13 → 128C, +N15+C13 → 129N, −C13 → 126, −N15 → none.
        assert_eq!(p.affected(idx("127C"), 4), Some(idx("128N")));
        assert_eq!(p.affected(idx("127C"), 5), Some(idx("128C")));
        assert_eq!(p.affected(idx("127C"), 6), Some(idx("129N")));
        assert_eq!(p.affected(idx("127C"), 2), Some(idx("126")));
        assert_eq!(p.affected(idx("127C"), 3), None);
        // 128N: −N15 → 127C, −N15−C13 → 126 … wait: 128N = 13C+15N; −N15 → 127C,
        // −N15−C13 → 126.
        assert_eq!(p.affected(idx("128N"), 3), Some(idx("127C")));
        assert_eq!(p.affected(idx("128N"), 1), Some(idx("126")));
    }

    #[test]
    fn parses_openms_style_matrix_and_corrects() {
        // TMT6 with a 5 % +1 spill from every channel except the last.
        let text = "126:0.0/0.0/5.0/0.0\n127:0.0/0.0/5.0/0.0\n128:0.0/0.0/5.0/0.0\n\
                    129:0.0/0.0/5.0/0.0\n130:0.0/0.0/5.0/0.0\n131:0.0/0.0/0.0/0.0\n";
        let m = CorrectionMatrix::parse(Plex::Tmt6, text).unwrap();
        assert!(!m.is_identity());
        assert!((m.get(0, 0) - 0.95).abs() < 1e-12);
        assert!((m.get(1, 0) - 0.05).abs() < 1e-12);
        assert!((m.get(5, 5) - 1.0).abs() < 1e-12);
        let truth = [1000.0f32, 2000.0, 500.0, 0.0, 100.0, 300.0];
        // observed = M · truth
        let observed: Vec<f32> = (0..6)
            .map(|t| (0..6).map(|s| m.get(t, s) as f32 * truth[s]).sum())
            .collect();
        let corrected = m.correct(&observed);
        for (c, t) in corrected.iter().zip(truth.iter()) {
            assert!((c - t).abs() < 1e-2, "{corrected:?}");
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
}
