//! Streaming MGF reader. Sage's regex-based pattern adapted to msgf-rust's
//! Spectrum shape. Sync I/O — MGF is line-oriented, no async benefit.

use std::io::BufRead;

use model::Spectrum;

pub struct MgfReader<R: BufRead> {
    reader: R,
    line_no: usize,
    /// Reusable line buffer to avoid per-line allocations.
    buf: String,
}

impl<R: BufRead> MgfReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader, line_no: 0, buf: String::new() }
    }

    /// Read the next non-blank, non-comment line. Returns `Ok(None)`
    /// at EOF. Advances `line_no`.
    fn next_significant_line(&mut self) -> Result<Option<String>, MgfParseError> {
        loop {
            self.buf.clear();
            let n = self.reader.read_line(&mut self.buf)
                .map_err(|source| MgfParseError::Io { line: self.line_no + 1, source })?;
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
                line: self.line_no, got: begin_line,
            }));
        }

        let begin_line_no = self.line_no;

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
                    return Some(Err(MgfParseError::UnterminatedSpectrum { line: begin_line_no }));
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
                    "TITLE"       => title = value,
                    "PEPMASS"     => {
                        match parse_pepmass(&value) {
                            Ok((mz, intensity)) => {
                                precursor_mz = Some(mz);
                                precursor_intensity = intensity;
                            }
                            Err(()) => return Some(Err(MgfParseError::BadPepmass {
                                line: self.line_no, got: value,
                            })),
                        }
                    }
                    "CHARGE"      => {
                        match parse_charge(&value) {
                            Ok(z) => precursor_charge = Some(z),
                            Err(()) => return Some(Err(MgfParseError::BadCharge {
                                line: self.line_no, got: value,
                            })),
                        }
                    }
                    "RTINSECONDS" => {
                        rt_seconds = value.parse().ok();
                    }
                    "SCANS"       => {
                        scan = value.parse().ok();
                    }
                    _ => { /* ignore unknown keys */ }
                }
                continue;
            }

            match parse_peak(&line) {
                Ok((mz, intensity)) => peaks.push((mz, intensity)),
                Err(()) => return Some(Err(MgfParseError::BadPeak {
                    line: self.line_no, got: line,
                })),
            }
        }

        let precursor_mz = match precursor_mz {
            Some(v) => v,
            None => return Some(Err(MgfParseError::MissingPepmass { line: begin_line_no })),
        };

        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        Some(Ok(Spectrum {
            title,
            precursor_mz,
            precursor_intensity,
            precursor_charge,
            rt_seconds,
            scan,
            peaks,
        }))
    }
}

fn parse_pepmass(value: &str) -> Result<(f64, Option<f32>), ()> {
    let mut iter = value.split_ascii_whitespace();
    let mz: f64 = iter.next().ok_or(())?.parse().map_err(|_| ())?;
    let intensity = iter.next().map(|s| s.parse::<f32>()).transpose().map_err(|_| ())?;
    Ok((mz, intensity))
}

fn parse_charge(value: &str) -> Result<i32, ()> {
    let trimmed = value.trim();
    let stripped = trimmed
        .strip_suffix('+')
        .or_else(|| trimmed.strip_suffix('-'))
        .unwrap_or(trimmed);
    stripped.parse().map_err(|_| ())
}

fn parse_peak(line: &str) -> Result<(f64, f32), ()> {
    let mut iter = line.split_ascii_whitespace();
    let mz: f64 = iter.next().ok_or(())?.parse().map_err(|_| ())?;
    let intensity: f32 = iter.next().ok_or(())?.parse().map_err(|_| ())?;
    Ok((mz, intensity))
}

#[derive(thiserror::Error, Debug)]
pub enum MgfParseError {
    #[error("I/O error at line {line}: {source}")]
    Io { line: usize, #[source] source: std::io::Error },

    #[error("expected `BEGIN IONS` at line {line}, got {got:?}")]
    ExpectedBeginIons { line: usize, got: String },

    #[error("unterminated spectrum starting at line {line} (no `END IONS` before EOF)")]
    UnterminatedSpectrum { line: usize },

    #[error("malformed PEPMASS at line {line}: {got:?}")]
    BadPepmass { line: usize, got: String },

    #[error("malformed CHARGE at line {line}: {got:?}")]
    BadCharge { line: usize, got: String },

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
        assert_eq!(parse_pepmass("500.5 1000.0").unwrap(), (500.5, Some(1000.0)));
    }

    #[test]
    fn parse_pepmass_without_intensity() {
        assert_eq!(parse_pepmass("500.5").unwrap(), (500.5, None));
    }

    #[test]
    fn parse_pepmass_garbage_errors() {
        assert!(parse_pepmass("garbage").is_err());
    }

    #[test]
    fn parse_charge_strips_plus() {
        assert_eq!(parse_charge("2+").unwrap(), 2);
        assert_eq!(parse_charge("3+").unwrap(), 3);
    }

    #[test]
    fn parse_charge_strips_minus() {
        assert_eq!(parse_charge("1-").unwrap(), 1);
    }

    #[test]
    fn parse_charge_no_sign_ok() {
        assert_eq!(parse_charge("4").unwrap(), 4);
    }

    #[test]
    fn parse_peak_space_separator() {
        assert_eq!(parse_peak("100.0 1.5").unwrap(), (100.0, 1.5));
    }

    #[test]
    fn parse_peak_tab_separator() {
        assert_eq!(parse_peak("100.0\t1.5").unwrap(), (100.0, 1.5));
    }

    #[test]
    fn parse_peak_garbage_errors() {
        assert!(parse_peak("not a peak").is_err());
    }
}
