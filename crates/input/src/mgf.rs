//! Streaming, regex-based MGF reader yielding andes's Spectrum.
//! Sync I/O — MGF is line-oriented.

use std::io::BufRead;

use model::Spectrum;

pub struct MgfReader<R: BufRead> {
    reader: R,
    line_no: usize,
    /// Reusable line buffer to avoid per-line allocations.
    buf: String,
    /// Number of spectra parsed so far.
    spectrum_index: usize,
}

impl<R: BufRead> MgfReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_no: 0,
            buf: String::new(),
            spectrum_index: 0,
        }
    }

    /// Read the next non-blank, non-comment line. Returns `Ok(None)`
    /// at EOF. Advances `line_no`.
    fn next_significant_line(&mut self) -> Result<Option<String>, MgfParseError> {
        loop {
            self.buf.clear();
            let n = self
                .reader
                .read_line(&mut self.buf)
                .map_err(|source| MgfParseError::Io {
                    line: self.line_no + 1,
                    source,
                })?;
            if n == 0 {
                return Ok(None);
            }
            self.line_no += 1;
            let trimmed = self.buf.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            return Ok(Some(trimmed.to_string()));
        }
    }
}

impl<R: BufRead> Iterator for MgfReader<R> {
    type Item = Result<Spectrum, MgfParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let begin_line = match self.next_significant_line() {
            Ok(None) => return None,
            Ok(Some(line)) => line,
            Err(e) => return Some(Err(e)),
        };

        if begin_line != "BEGIN IONS" {
            return Some(Err(MgfParseError::ExpectedBeginIons {
                line: self.line_no,
                got: begin_line,
            }));
        }

        let begin_line_no = self.line_no;
        self.spectrum_index += 1;

        let mut title = String::new();
        let mut precursor_mz: Option<f64> = None;
        let mut precursor_intensity: Option<f32> = None;
        let mut precursor_charge: Option<i32> = None;
        let mut rt_seconds: Option<f64> = None;
        let mut scan: Option<i32> = None;
        let mut peaks: Vec<(f64, f32)> = Vec::new();

        loop {
            let line = match self.next_significant_line() {
                Ok(None) => {
                    return Some(Err(MgfParseError::UnterminatedSpectrum {
                        line: begin_line_no,
                    }));
                }
                Ok(Some(l)) => l,
                Err(e) => return Some(Err(e)),
            };

            if line == "END IONS" {
                break;
            }

            if let Some(eq) = line.find('=') {
                let key = line[..eq].to_ascii_uppercase();
                let value = line[eq + 1..].trim().to_string();
                match key.as_str() {
                    "TITLE" => title = value,
                    "PEPMASS" => match parse_pepmass(&value) {
                        Ok((mz, intensity)) => {
                            precursor_mz = Some(mz);
                            precursor_intensity = intensity;
                        }
                        Err(()) => {
                            return Some(Err(MgfParseError::BadPepmass {
                                line: self.line_no,
                                got: value,
                            }))
                        }
                    },
                    "CHARGE" => match parse_charge(&value) {
                        // `CHARGE=0` means "unknown", as in mzML: leave it absent
                        // so the search enumerates the charge range; so does an
                        // out-of-range value or a list of several charge states.
                        Ok(z) => {
                            precursor_charge = z
                                .and_then(|z| model::spectrum::known_precursor_charge(i64::from(z)))
                        }
                        Err(ChargeError::Malformed) => {
                            return Some(Err(MgfParseError::BadCharge {
                                line: self.line_no,
                                got: value,
                            }))
                        }
                        Err(ChargeError::Negative) => {
                            return Some(Err(MgfParseError::NegativeCharge {
                                line: self.line_no,
                                got: value,
                            }))
                        }
                    },
                    "RTINSECONDS" => {
                        rt_seconds = value.parse().ok();
                    }
                    "SCANS" => {
                        scan = value.parse().ok();
                    }
                    _ => { /* ignore unknown keys */ }
                }
                continue;
            }

            match parse_peak(&line) {
                // Drop non-finite / non-positive peaks (finding 4.4), mirroring
                // the mzML / timsTOF readers; a syntactically-malformed line is
                // still a hard error.
                Ok(Some((mz, intensity))) => peaks.push((mz, intensity)),
                Ok(None) => { /* peak out of domain: skip */ }
                Err(()) => {
                    return Some(Err(MgfParseError::BadPeak {
                        line: self.line_no,
                        got: line,
                    }))
                }
            }
        }

        let precursor_mz = match precursor_mz {
            Some(v) => v,
            None => {
                return Some(Err(MgfParseError::MissingPepmass {
                    line: begin_line_no,
                }))
            }
        };

        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // `SCANS=` wins; msconvert omits it and puts the scan in the TITLE
        // instead. Without either, the 1-based position in the file keeps
        // the PIN `ScanNr` unique.
        let scan = scan
            .or_else(|| scan_from_title(&title))
            .or_else(|| i32::try_from(self.spectrum_index).ok());

        Some(Ok(Spectrum {
            title,
            precursor_mz,
            precursor_intensity,
            precursor_charge,
            rt_seconds,
            scan,
            peaks,
            // MGF doesn't carry an activation method in the standard
            // header set; could be extended via a custom `ACTIVATION=`
            // field if needed. For now: leave it absent.
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        }))
    }
}

fn parse_pepmass(value: &str) -> Result<(f64, Option<f32>), ()> {
    let mut iter = value.split_ascii_whitespace();
    let mz: f64 = iter.next().ok_or(())?.parse().map_err(|_| ())?;
    // Reject a non-finite or non-positive precursor m/z (e.g. `PEPMASS=0` or
    // `inf`) — it would produce nonsense search windows. The Thermo and timsTOF
    // readers gate on `precursor_mz > 0`; do the same here for consistency.
    if !mz.is_finite() || mz <= 0.0 {
        return Err(());
    }
    let intensity = iter
        .next()
        .map(|s| s.parse::<f32>())
        .transpose()
        .map_err(|_| ())?;
    Ok((mz, intensity))
}

/// Scan number carried in an MGF `TITLE`: a native-id `scan=N` token
/// (`... scan=1234`, `NativeID:"controllerType=0 controllerNumber=1
/// scan=1234"`), or the TPP/msconvert `name.START.END.CHARGE` form.
fn scan_from_title(title: &str) -> Option<i32> {
    if let Some(pos) = title.find("scan=") {
        let digits: String = title[pos + "scan=".len()..]
            .chars()
            .take_while(char::is_ascii_digit)
            .collect();
        if let Ok(scan) = digits.parse::<i32>() {
            return Some(scan);
        }
    }
    // `name.START.END.CHARGE`, possibly followed by a space and free text.
    let head = title.split_whitespace().next()?;
    let mut fields = head.rsplit('.');
    let charge = fields.next()?;
    let end = fields.next()?;
    let start = fields.next()?;
    fields.next()?;
    if charge.is_empty() || !charge.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let (start, end) = (start.parse::<i32>().ok()?, end.parse::<i32>().ok()?);
    (start > 0 && end >= start).then_some(start)
}

/// Why a `CHARGE=` value cannot be used.
#[derive(Debug, PartialEq, Eq)]
enum ChargeError {
    Malformed,
    /// A negative charge state: the search computes the precursor mass of a
    /// protonated ion, which is wrong for a negative-ion spectrum.
    Negative,
}

/// Parse a `CHARGE=` value: one charge state (`2+`, `2`, `+2`) or several
/// (`2+ and 3+`, `2+,3+`), as Mascot and msconvert write when the charge was
/// not determined.
///
/// `Ok(Some(z))` when the value names one charge state. `Ok(None)` when it
/// names several different ones: the charge is then unknown, and the search
/// enumerates its configured charge range for the spectrum.
fn parse_charge(value: &str) -> Result<Option<i32>, ChargeError> {
    let mut charge: Option<i32> = None;
    let mut several = false;
    let tokens = value
        .split(|c: char| c.is_whitespace() || c == ',')
        .filter(|t| !t.is_empty() && !t.eq_ignore_ascii_case("and"));
    for token in tokens {
        let (digits, negative) = if let Some(d) = token.strip_suffix('-') {
            (d, true)
        } else if let Some(d) = token.strip_prefix('-') {
            (d, true)
        } else {
            let d = token.strip_suffix('+').unwrap_or(token);
            (d.strip_prefix('+').unwrap_or(d), false)
        };
        if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
            return Err(ChargeError::Malformed);
        }
        let z: i32 = digits.parse().map_err(|_| ChargeError::Malformed)?;
        if negative && z != 0 {
            return Err(ChargeError::Negative);
        }
        match charge {
            None => charge = Some(z),
            Some(first) if first != z => several = true,
            Some(_) => {}
        }
    }
    match charge {
        None => Err(ChargeError::Malformed),
        Some(_) if several => Ok(None),
        Some(z) => Ok(Some(z)),
    }
}

/// Parse one `mz intensity` peak line.
///
/// `Err(())` is a syntactically malformed line (hard error upstream). `Ok(None)`
/// is a syntactically-valid but out-of-domain peak (non-finite / non-positive
/// m/z, or non-finite / negative intensity) that must be dropped rather than
/// scored (finding 4.4). `Ok(Some(..))` is a usable peak.
fn parse_peak(line: &str) -> Result<Option<(f64, f32)>, ()> {
    let mut iter = line.split_ascii_whitespace();
    let mz: f64 = iter.next().ok_or(())?.parse().map_err(|_| ())?;
    let intensity: f32 = iter.next().ok_or(())?.parse().map_err(|_| ())?;
    if !mz.is_finite() || mz <= 0.0 || !intensity.is_finite() || intensity < 0.0 {
        return Ok(None);
    }
    Ok(Some((mz, intensity)))
}

#[derive(thiserror::Error, Debug)]
pub enum MgfParseError {
    #[error("I/O error at line {line}: {source}")]
    Io {
        line: usize,
        #[source]
        source: std::io::Error,
    },

    #[error("expected `BEGIN IONS` at line {line}, got {got:?}")]
    ExpectedBeginIons { line: usize, got: String },

    #[error("unterminated spectrum starting at line {line} (no `END IONS` before EOF)")]
    UnterminatedSpectrum { line: usize },

    #[error("malformed PEPMASS at line {line}: {got:?}")]
    BadPepmass { line: usize, got: String },

    #[error("malformed CHARGE at line {line}: {got:?}")]
    BadCharge { line: usize, got: String },

    #[error(
        "negative CHARGE at line {line}: {got:?}; negative-ion spectra are not supported \
         (the precursor mass is computed for a protonated ion)"
    )]
    NegativeCharge { line: usize, got: String },

    #[error("malformed peak line at line {line}: expected `mz intensity`, got {got:?}")]
    BadPeak { line: usize, got: String },

    #[error("missing PEPMASS in spectrum starting at line {line}")]
    MissingPepmass { line: usize },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_pepmass_with_intensity() {
        assert_eq!(
            parse_pepmass("500.5 1000.0").unwrap(),
            (500.5, Some(1000.0))
        );
    }

    #[test]
    fn parse_pepmass_without_intensity() {
        assert_eq!(parse_pepmass("500.5").unwrap(), (500.5, None));
    }

    #[test]
    fn parse_pepmass_rejects_nonpositive_and_nonfinite() {
        // A non-positive or non-finite precursor m/z would produce nonsense
        // search windows — reject it, matching the Thermo/timsTOF readers.
        assert!(parse_pepmass("0").is_err());
        assert!(parse_pepmass("0.0 1000.0").is_err());
        assert!(parse_pepmass("-5.0").is_err());
        assert!(parse_pepmass("inf").is_err());
        assert!(parse_pepmass("nan").is_err());
    }

    #[test]
    fn parse_pepmass_garbage_errors() {
        assert!(parse_pepmass("garbage").is_err());
    }

    #[test]
    fn parse_charge_single_value_forms() {
        assert_eq!(parse_charge("2+"), Ok(Some(2)));
        assert_eq!(parse_charge("3+"), Ok(Some(3)));
        assert_eq!(parse_charge("4"), Ok(Some(4)));
        assert_eq!(parse_charge("+2"), Ok(Some(2)));
        assert_eq!(parse_charge(" 2+ "), Ok(Some(2)));
        assert_eq!(parse_charge("0"), Ok(Some(0)));
    }

    #[test]
    fn parse_charge_several_states_is_unknown() {
        assert_eq!(parse_charge("2+ and 3+"), Ok(None));
        assert_eq!(parse_charge("2+,3+"), Ok(None));
        assert_eq!(parse_charge("3 and 4"), Ok(None));
        assert_eq!(parse_charge("2+, 3+ and 4+"), Ok(None));
    }

    #[test]
    fn parse_charge_repeated_state_is_that_state() {
        assert_eq!(parse_charge("2+ and 2+"), Ok(Some(2)));
        assert_eq!(parse_charge("3,3"), Ok(Some(3)));
    }

    #[test]
    fn parse_charge_negative_is_rejected() {
        assert_eq!(parse_charge("1-"), Err(ChargeError::Negative));
        assert_eq!(parse_charge("-2"), Err(ChargeError::Negative));
        assert_eq!(parse_charge("2+ and 2-"), Err(ChargeError::Negative));
    }

    #[test]
    fn parse_charge_malformed_is_rejected() {
        for value in [
            "", "and", "banana", "2+ or 3+", "+", "2++", "+-2", "2.0", "2+and3+",
        ] {
            assert_eq!(
                parse_charge(value),
                Err(ChargeError::Malformed),
                "{value:?}"
            );
        }
    }

    fn read_all(mgf: &str) -> Vec<Spectrum> {
        MgfReader::new(std::io::Cursor::new(mgf.as_bytes()))
            .map(|r| r.unwrap())
            .collect()
    }

    #[test]
    fn charge_zero_is_unknown() {
        let s = read_all("BEGIN IONS\nPEPMASS=500.5\nCHARGE=0\n100 1\nEND IONS\n");
        assert_eq!(s[0].precursor_charge, None);
        let s = read_all("BEGIN IONS\nPEPMASS=500.5\nCHARGE=2+\n100 1\nEND IONS\n");
        assert_eq!(s[0].precursor_charge, Some(2));
    }

    #[test]
    fn several_charge_states_leave_the_charge_unknown() {
        let s = read_all("BEGIN IONS\nPEPMASS=500.5\nCHARGE=2+ and 3+\n100 1\nEND IONS\n");
        assert_eq!(s[0].precursor_charge, None);
    }

    #[test]
    fn negative_charge_is_an_error_naming_the_line() {
        let mgf = "BEGIN IONS\nPEPMASS=500.5\nCHARGE=2-\n100 1\nEND IONS\n";
        let err = MgfReader::new(std::io::Cursor::new(mgf.as_bytes()))
            .next()
            .unwrap()
            .unwrap_err();
        assert!(
            matches!(err, MgfParseError::NegativeCharge { line: 3, .. }),
            "{err}"
        );
        assert!(
            err.to_string()
                .contains("negative-ion spectra are not supported"),
            "{err}"
        );
    }

    #[test]
    fn out_of_range_charge_is_unknown() {
        let s = read_all("BEGIN IONS\nPEPMASS=500.5\nCHARGE=300+\n100 1\nEND IONS\n");
        assert_eq!(s[0].precursor_charge, None);
    }

    #[test]
    fn scan_from_title_forms() {
        assert_eq!(scan_from_title("run01.1234.1234.2"), Some(1234));
        assert_eq!(
            scan_from_title("run01.1234.1234.2 File:\"run01.raw\", NativeID:\"x\""),
            Some(1234)
        );
        assert_eq!(
            scan_from_title(
                "run01.1234.1234.2 File:\"run01.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=1234\""
            ),
            Some(1234)
        );
        assert_eq!(
            scan_from_title("controllerType=0 controllerNumber=1 scan=77"),
            Some(77)
        );
        assert_eq!(scan_from_title("my.sample.v2.3"), None);
        assert_eq!(scan_from_title("spectrum 5"), None);
        assert_eq!(scan_from_title(""), None);
    }

    #[test]
    fn scan_number_precedence() {
        // SCANS= wins over the title; the title wins over the file position;
        // the 1-based file position keeps ScanNr unique otherwise.
        let mgf = "BEGIN IONS\nTITLE=a.10.10.2\nSCANS=99\nPEPMASS=500\n100 1\nEND IONS\n\
                   BEGIN IONS\nTITLE=a.20.20.2\nPEPMASS=500\n100 1\nEND IONS\n\
                   BEGIN IONS\nTITLE=no scan here\nPEPMASS=500\n100 1\nEND IONS\n";
        let scans: Vec<_> = read_all(mgf).iter().map(|s| s.scan).collect();
        assert_eq!(scans, vec![Some(99), Some(20), Some(3)]);
    }

    #[test]
    fn parse_peak_space_separator() {
        assert_eq!(parse_peak("100.0 1.5").unwrap(), Some((100.0, 1.5)));
    }

    #[test]
    fn parse_peak_tab_separator() {
        assert_eq!(parse_peak("100.0\t1.5").unwrap(), Some((100.0, 1.5)));
    }

    #[test]
    fn parse_peak_garbage_errors() {
        assert!(parse_peak("not a peak").is_err());
    }

    #[test]
    fn parse_peak_out_of_domain_is_dropped() {
        // Finding 4.4: syntactically valid but out-of-domain peaks are dropped
        // (Ok(None)), not errored, mirroring the mzML / timsTOF readers.
        assert_eq!(parse_peak("0 100").unwrap(), None, "zero m/z");
        assert_eq!(parse_peak("-5 100").unwrap(), None, "negative m/z");
        assert_eq!(parse_peak("inf 100").unwrap(), None, "inf m/z");
        assert_eq!(parse_peak("nan 100").unwrap(), None, "nan m/z");
        assert_eq!(parse_peak("100 -1").unwrap(), None, "negative intensity");
        assert_eq!(parse_peak("100 inf").unwrap(), None, "inf intensity");
        assert_eq!(parse_peak("100 nan").unwrap(), None, "nan intensity");
        // Zero intensity is allowed (a real, if uninformative, peak).
        assert_eq!(parse_peak("100 0").unwrap(), Some((100.0, 0.0)));
    }
}
