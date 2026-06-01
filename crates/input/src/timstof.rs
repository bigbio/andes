//! Native Bruker timsTOF `.d` reader (`feature = "timstof"`).
//!
//! Wraps the pure-Rust [`timsrust`] crate (the same reader Sage uses) and
//! yields the same [`Spectrum`] model as the mzML/MGF/Thermo readers, so the
//! search path is format-agnostic. A Bruker `.d` is a DIRECTORY holding a TDF
//! SQLite database (`analysis.tdf`) plus its binary blob (`analysis.tdf_bin`);
//! `timsrust` reads both natively with NO vendor runtime and NO bundling
//! (unlike the Thermo `.raw` reader, which needs the hosted .NET runtime).
//!
//! Built only under `--features timstof`. mzML/MGF reading never pulls in
//! `timsrust`.
//!
//! Scope: DDA-PASEF, MS2 only. `timsrust`'s [`SpectrumReader`] already groups
//! the TIMS frames into centroided MS2 fragment spectra, each carrying the DDA
//! precursor it was selected from — exactly the unit the search needs. The ion
//! mobility dimension is extra metadata MS-GF+ scoring does not use; RT and
//! precursor are carried as usual and mobility is ignored for the base search
//! (a future Percolator-feature idea, not implemented here). There is no MS1
//! stream here, so `--chimeric` on `.d` degrades gracefully to a normal search
//! (handled by the binary, like MGF).

use std::path::Path;

use timsrust::readers::{SpectrumReader, SpectrumReaderError};

use crate::Spectrum;

/// Error opening or reading a Bruker timsTOF `.d` directory.
#[derive(Debug, thiserror::Error)]
pub enum TimsTofParseError {
    #[error("failed to open Bruker .d '{path}': {source}")]
    Open {
        path: String,
        source: SpectrumReaderError,
    },
    #[error("failed to read spectrum {index} from Bruker .d: {source}")]
    Read {
        index: usize,
        source: SpectrumReaderError,
    },
}

/// Reader over a Bruker timsTOF `.d` directory, yielding [`Spectrum`] values.
///
/// `timsrust`'s [`SpectrumReader`] exposes a flat, index-addressed list of
/// centroided MS2 spectra (`len` / `get`), so this reader iterates by index and
/// owns the underlying handle without a self-referential borrow — matching the
/// `ThermoRawReader` structure. Every emitted spectrum is MS2 by construction
/// (the reader only produces fragment spectra), so there is no MS-level filter.
pub struct TimsTofReader {
    reader: SpectrumReader,
    next: usize,
    len: usize,
}

impl TimsTofReader {
    /// Open a `.d` directory. The path is the `.d` folder itself (which holds
    /// `analysis.tdf` + `analysis.tdf_bin`), not a file inside it.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, TimsTofParseError> {
        let p = path.as_ref();
        let reader = SpectrumReader::new(p).map_err(|source| TimsTofParseError::Open {
            path: p.display().to_string(),
            source,
        })?;
        let len = reader.len();
        Ok(Self {
            reader,
            next: 0,
            len,
        })
    }

    /// Number of MS2 spectra in the `.d`.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the `.d` has no spectra.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Convert one `timsrust` spectrum to a [`Spectrum`], or `None` if the scan
    /// has no usable DDA precursor m/z — matching the mzML/Thermo readers, which
    /// skip precursor-less MS2 rather than forwarding a `precursor_mz = 0.0`
    /// spectrum that would search against a nonsense mass window.
    // `!(mz > 0.0)` deliberately rejects NaN as well as 0 / negative m/z — the
    // same precursor-sanity idiom the Thermo reader uses.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn convert(raw: &timsrust::Spectrum) -> Option<Spectrum> {
        // The 0-based `timsrust` index becomes a 1-based scan number, matching
        // the convention the mzML/Thermo readers use for `SCANS=`/scan.
        let scan_number = raw.index as i32 + 1;

        // A DDA precursor with a positive m/z is required to search an MS2.
        let precursor = raw.precursor.as_ref()?;
        let precursor_mz = precursor.mz;
        if !(precursor_mz > 0.0) {
            return None;
        }

        let precursor_charge = precursor.charge.and_then(|z| {
            // 0 (or an absurd value) means "unknown"; let the engine sweep the
            // configured charge range instead of trusting a bogus state.
            if z == 0 {
                None
            } else {
                Some(z as i32)
            }
        });

        let precursor_intensity = precursor
            .intensity
            .and_then(|i| if i > 0.0 { Some(i as f32) } else { None });

        // The quadrupole isolation window is reported as a center m/z + total
        // width. Convert to the symmetric lower/upper offsets the model stores
        // (mirroring how the mzML reader records `isolation window lower/upper
        // offset`). Only `--chimeric` consumes these, and chimeric on `.d` runs
        // as a normal search, but they are filled in for completeness when the
        // width is positive.
        let (isolation_lower_offset, isolation_upper_offset) = if raw.isolation_width > 0.0 {
            let half = raw.isolation_width / 2.0;
            (Some(half), Some(half))
        } else {
            (None, None)
        };

        Some(Spectrum {
            // The 1-based scan index uniquely identifies the spectrum within the
            // `.d` and becomes the PIN `SpecID` column.
            title: format!("scan={scan_number}"),
            precursor_mz,
            precursor_intensity,
            precursor_charge,
            // `timsrust` retention time is in SECONDS (Sage divides it by 60 to
            // get minutes); the model stores seconds, so it is carried as-is.
            rt_seconds: Some(precursor.rt),
            scan: Some(scan_number),
            peaks: extract_peaks(raw),
            // timsTOF DDA-PASEF is CID/HCD-style beam-type fragmentation; the
            // `.d` does not record a discrete activation cvParam the resolver
            // keys on, so leave this `None` and let the param resolver use its
            // default (or an explicit `--fragmentation`/`--instrument`).
            activation_method: None,
            isolation_lower_offset,
            isolation_upper_offset,
        })
    }
}

/// Extract centroided peaks `(m/z, intensity)` from a `timsrust` spectrum,
/// ascending by m/z. `timsrust` stores parallel `mz_values` / `intensities`
/// vectors of `f64`; intensity is narrowed to `f32` to match the model.
///
/// Downstream consumers require m/z-ascending peaks (the scorer and chimeric
/// co-isolation detection binary-search the list). `timsrust` returns peaks
/// sorted by m/z, so the list is only re-sorted if an inversion is present.
fn extract_peaks(raw: &timsrust::Spectrum) -> Vec<(f64, f32)> {
    let mut peaks: Vec<(f64, f32)> = raw
        .mz_values
        .iter()
        .zip(raw.intensities.iter())
        .map(|(&mz, &intensity)| (mz, intensity as f32))
        .collect();
    if peaks.windows(2).any(|w| w[0].0 > w[1].0) {
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    peaks
}

impl Iterator for TimsTofReader {
    type Item = Result<Spectrum, TimsTofParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next < self.len {
            let i = self.next;
            self.next += 1;
            match self.reader.get(i) {
                Ok(raw) => match Self::convert(&raw) {
                    // Skip MS2 scans with no usable precursor (mirrors mzML).
                    Some(spec) => return Some(Ok(spec)),
                    None => continue,
                },
                Err(source) => return Some(Err(TimsTofParseError::Read { index: i, source })),
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsrust::Precursor;

    /// Build a minimal `timsrust::Spectrum` for the conversion tests.
    fn spectrum(
        index: usize,
        precursor: Option<Precursor>,
        isolation_width: f64,
        mz_values: Vec<f64>,
        intensities: Vec<f64>,
    ) -> timsrust::Spectrum {
        timsrust::Spectrum {
            mz_values,
            intensities,
            precursor,
            index,
            isolation_width,
            ..Default::default()
        }
    }

    #[test]
    fn convert_maps_precursor_and_peaks() {
        let precursor = Precursor {
            mz: 654.32,
            rt: 1800.0, // seconds
            charge: Some(2),
            intensity: Some(5000.0),
            ..Default::default()
        };
        let raw = spectrum(
            41,
            Some(precursor),
            2.0,
            vec![100.0, 200.5, 300.0],
            vec![10.0, 20.0, 30.0],
        );

        let spec = TimsTofReader::convert(&raw).expect("convertible spectrum");
        // 0-based index becomes a 1-based scan number.
        assert_eq!(spec.scan, Some(42));
        assert_eq!(spec.title, "scan=42");
        assert_eq!(spec.precursor_mz, 654.32);
        assert_eq!(spec.precursor_charge, Some(2));
        assert_eq!(spec.precursor_intensity, Some(5000.0));
        // RT is carried through in seconds (NOT divided by 60).
        assert_eq!(spec.rt_seconds, Some(1800.0));
        // Symmetric offsets are half the total isolation width.
        assert_eq!(spec.isolation_lower_offset, Some(1.0));
        assert_eq!(spec.isolation_upper_offset, Some(1.0));
        assert_eq!(
            spec.peaks,
            vec![(100.0, 10.0_f32), (200.5, 20.0_f32), (300.0, 30.0_f32)]
        );
    }

    #[test]
    fn convert_skips_precursorless_spectrum() {
        let raw = spectrum(0, None, 2.0, vec![100.0], vec![1.0]);
        assert!(TimsTofReader::convert(&raw).is_none());
    }

    #[test]
    fn convert_skips_nonpositive_precursor_mz() {
        let precursor = Precursor { mz: 0.0, ..Default::default() };
        let raw = spectrum(0, Some(precursor), 2.0, vec![100.0], vec![1.0]);
        assert!(TimsTofReader::convert(&raw).is_none());
    }

    #[test]
    fn convert_treats_zero_charge_as_unknown() {
        let precursor = Precursor { mz: 500.0, charge: Some(0), ..Default::default() };
        let raw = spectrum(0, Some(precursor), 0.0, vec![100.0], vec![1.0]);
        let spec = TimsTofReader::convert(&raw).expect("positive m/z is convertible");
        assert_eq!(spec.precursor_charge, None);
        // A non-positive isolation width yields no offsets.
        assert_eq!(spec.isolation_lower_offset, None);
        assert_eq!(spec.isolation_upper_offset, None);
    }

    #[test]
    fn extract_peaks_sorts_unordered_input() {
        let raw = spectrum(
            0,
            None,
            0.0,
            vec![300.0, 100.0, 200.0],
            vec![3.0, 1.0, 2.0],
        );
        let peaks = extract_peaks(&raw);
        assert_eq!(
            peaks,
            vec![(100.0, 1.0_f32), (200.0, 2.0_f32), (300.0, 3.0_f32)]
        );
    }
}
