//! Native Bruker timsTOF `.d` reader (`feature = "timstof"`).
//!
//! Wraps the pure-Rust [`timsrust`] crate and yields the same [`Spectrum`]
//! model as the mzML/MGF/Thermo readers, so the search path is
//! format-agnostic. A Bruker `.d` is a DIRECTORY holding a TDF SQLite database
//! (`analysis.tdf`) plus its binary blob (`analysis.tdf_bin`); `timsrust` reads
//! both natively with NO vendor runtime and NO bundling.
//!
//! Built only under `--features timstof`. mzML/MGF reading never pulls in
//! `timsrust`.
//!
//! Scope: DDA-PASEF, MS2 only. `timsrust`'s [`SpectrumReader`] already groups
//! the TIMS frames into centroided MS2 fragment spectra, each carrying the DDA
//! precursor it was selected from — exactly the unit the search needs. The ion
//! mobility dimension is extra metadata the rank-based scorer (Kim et al.,
//! Nat Commun 5:5277, 2014) does not use; RT and
//! precursor are carried as usual and mobility is ignored for the base search
//! (a future Percolator-feature idea, not implemented here). There is no MS1
//! stream here, so `--chimeric` on `.d` degrades gracefully to a normal search
//! (handled by the binary, like MGF).

use std::path::{Path, PathBuf};

use timsrust::readers::{SpectrumReader, SpectrumReaderError};

use crate::Spectrum;

/// Error opening or reading a Bruker timsTOF `.d` directory.
#[derive(Debug, thiserror::Error)]
pub enum TimsTofParseError {
    #[error("failed to open Bruker .d '{path}': {source}{hint}")]
    Open {
        path: String,
        source: SpectrumReaderError,
        /// Extra guidance for a known unsupported layout (empty otherwise).
        hint: &'static str,
    },
    #[error(
        "'{path}' is not a Bruker .d run folder: it has no analysis.tdf, but the folder \
         '{nested}' inside it has one; pass that folder"
    )]
    Nested { path: String, nested: String },
    #[error("failed to read spectrum {index} from Bruker .d: {source}")]
    Read {
        index: usize,
        source: SpectrumReaderError,
    },
}

/// Guidance appended to an open error whose cause is a known unsupported `.d`
/// layout. Frame compression type 1 is the scheme written by early timsTOF
/// acquisition software (data from about 2018 and before); the reader only
/// decodes the later type 2.
fn open_error_hint(message: &str) -> &'static str {
    if message.contains("Compression type 1") {
        " (this .d uses the frame compression of early timsTOF acquisition software, \
         which the native reader does not decode; convert the run to mzML and search that)"
    } else {
        ""
    }
}

/// The run folder inside `dir` when `dir` has no `analysis.tdf` of its own but
/// a `.d` folder inside it has one: archives that wrap the run in a folder of
/// the same name unpack to `X.d/X.d/`. The reader panics on the outer folder,
/// so it is caught here first.
fn nested_run_folder(dir: &Path) -> Option<PathBuf> {
    if !dir.is_dir() || dir.join("analysis.tdf").exists() {
        return None;
    }
    let mut nested: Vec<PathBuf> = std::fs::read_dir(dir)
        .ok()?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| {
            p.is_dir()
                && p.extension().is_some_and(|x| x.eq_ignore_ascii_case("d"))
                && p.join("analysis.tdf").exists()
        })
        .collect();
    nested.sort();
    nested.into_iter().next()
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
        if let Some(nested) = nested_run_folder(p) {
            return Err(TimsTofParseError::Nested {
                path: p.display().to_string(),
                nested: nested.display().to_string(),
            });
        }
        let reader = SpectrumReader::new(p).map_err(|source| {
            let hint = open_error_hint(&source.to_string());
            TimsTofParseError::Open {
                path: p.display().to_string(),
                source,
                hint,
            }
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

        // 0 (or an absurd value) means "unknown"; let the engine sweep the
        // configured charge range instead of trusting a bogus state.
        let precursor_charge = precursor
            .charge
            .and_then(|z| i64::try_from(z).ok())
            .and_then(model::spectrum::known_precursor_charge);

        let precursor_intensity =
            precursor
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
            // `timsrust` retention time is in SECONDS; the model stores
            // seconds, so it is carried as-is (RT seconds → minutes downstream).
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
            isolation_target_mz: None,
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
///
/// The two arrays are paired only over their overlapping prefix, so a malformed
/// length mismatch can't read past either vector, and non-finite or
/// non-positive-m/z / negative-intensity points are dropped before they can
/// reach the scorer.
fn extract_peaks(raw: &timsrust::Spectrum) -> Vec<(f64, f32)> {
    let n = raw.mz_values.len().min(raw.intensities.len());
    let mut peaks: Vec<(f64, f32)> = (0..n)
        .map(|i| (raw.mz_values[i], raw.intensities[i] as f32))
        .filter(|&(mz, intensity)| {
            mz.is_finite() && mz > 0.0 && intensity.is_finite() && intensity >= 0.0
        })
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

    #[test]
    fn compression_type_1_gets_a_conversion_hint() {
        assert!(open_error_hint("Compression type 1 not understood").contains("mzML"));
        assert_eq!(open_error_hint("No such file or directory"), "");
    }

    #[test]
    fn a_run_folder_nested_in_its_own_name_is_an_error_naming_the_inner_folder() {
        let tmp = tempfile::tempdir().unwrap();
        let outer = tmp.path().join("run.d");
        let inner = outer.join("run.d");
        std::fs::create_dir_all(&inner).unwrap();
        std::fs::write(inner.join("analysis.tdf"), b"").unwrap();
        assert_eq!(nested_run_folder(&outer), Some(inner.clone()));
        assert_eq!(nested_run_folder(&inner), None);
        let err = TimsTofReader::open(&outer)
            .err()
            .expect("the outer folder is not a run");
        assert!(matches!(err, TimsTofParseError::Nested { .. }), "{err}");
        assert!(
            err.to_string().contains(&inner.display().to_string()),
            "{err}"
        );
    }

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
        let precursor = Precursor {
            mz: 0.0,
            ..Default::default()
        };
        let raw = spectrum(0, Some(precursor), 2.0, vec![100.0], vec![1.0]);
        assert!(TimsTofReader::convert(&raw).is_none());
    }

    #[test]
    fn convert_treats_zero_charge_as_unknown() {
        let precursor = Precursor {
            mz: 500.0,
            charge: Some(0),
            ..Default::default()
        };
        let raw = spectrum(0, Some(precursor), 0.0, vec![100.0], vec![1.0]);
        let spec = TimsTofReader::convert(&raw).expect("positive m/z is convertible");
        assert_eq!(spec.precursor_charge, None);
        // A non-positive isolation width yields no offsets.
        assert_eq!(spec.isolation_lower_offset, None);
        assert_eq!(spec.isolation_upper_offset, None);
    }

    #[test]
    fn extract_peaks_sorts_unordered_input() {
        let raw = spectrum(0, None, 0.0, vec![300.0, 100.0, 200.0], vec![3.0, 1.0, 2.0]);
        let peaks = extract_peaks(&raw);
        assert_eq!(
            peaks,
            vec![(100.0, 1.0_f32), (200.0, 2.0_f32), (300.0, 3.0_f32)]
        );
    }

    #[test]
    fn extract_peaks_drops_nonfinite_and_nonpositive() {
        let raw = spectrum(
            0,
            None,
            0.0,
            vec![100.0, f64::NAN, 0.0, 200.0, f64::INFINITY, 300.0],
            vec![1.0, 5.0, 5.0, f64::NAN, 5.0, -1.0],
        );
        // Drops: NaN m/z, 0.0 m/z, +inf m/z, NaN intensity, negative intensity.
        let peaks = extract_peaks(&raw);
        assert_eq!(peaks, vec![(100.0, 1.0_f32)]);
    }

    #[test]
    fn extract_peaks_pairs_only_overlapping_prefix() {
        // A malformed length mismatch must not read past the shorter vector.
        let raw = spectrum(0, None, 0.0, vec![100.0, 200.0, 300.0], vec![1.0, 2.0]);
        let peaks = extract_peaks(&raw);
        assert_eq!(peaks, vec![(100.0, 1.0_f32), (200.0, 2.0_f32)]);
    }
}
