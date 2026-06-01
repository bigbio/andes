//! Native Thermo `.raw` reader (`feature = "thermo"`).
//!
//! Wraps Thermo's official RawFileReader (hosted .NET runtime) via the
//! [`thermorawfilereader`] crate and yields the same [`Spectrum`] model as the
//! mzML/MGF readers, so the search path is format-agnostic.
//!
//! Built only under `--features thermo`. The build needs no .NET SDK (the
//! RawFileReader assemblies are vendored in the dependency); opening a `.raw`
//! at runtime requires the .NET 8 runtime, auto-discovered via hostfxr
//! (`DOTNET_ROOT` or a system install). mzML/MGF reading never loads .NET.

use std::path::Path;

use thermorawfilereader::{RawFileReader, RawSpectrum};

use crate::Spectrum;

/// Error opening a Thermo `.raw` file.
#[derive(Debug, thiserror::Error)]
pub enum ThermoParseError {
    #[error("failed to open Thermo .raw '{path}': {source}")]
    Open {
        path: String,
        source: std::io::Error,
    },
}

/// Reader over a Thermo `.raw`, yielding [`Spectrum`] values.
///
/// Iterates by index through the underlying `RawFileReader::get`, so it owns
/// the file handle without a self-referential borrow. By default only MS2
/// scans are emitted (the search path); call [`with_all_ms_levels`] to include
/// MS1 as well (required by the chimeric cascade's MS1 gating).
///
/// [`with_all_ms_levels`]: ThermoRawReader::with_all_ms_levels
pub struct ThermoRawReader {
    handle: RawFileReader,
    next: usize,
    len: usize,
    /// `Some(level)` emits only that MS level; `None` emits all levels.
    ms_level_filter: Option<u8>,
}

impl ThermoRawReader {
    /// Open a `.raw`. Requests centroided peaks, matching the centroid data the
    /// mzML path feeds the scorer.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, ThermoParseError> {
        let p = path.as_ref();
        let mut handle = RawFileReader::open(p).map_err(|source| ThermoParseError::Open {
            path: p.display().to_string(),
            source,
        })?;
        handle.set_centroid_spectra(true);
        let len = handle.len();
        Ok(Self {
            handle,
            next: 0,
            len,
            ms_level_filter: Some(2),
        })
    }

    /// Emit all MS levels (MS1 + MS2) instead of MS2-only. Needed for `--chimeric`.
    pub fn with_all_ms_levels(mut self) -> Self {
        self.ms_level_filter = None;
        self
    }

    /// Restrict to a single MS level, or `None` for all levels.
    pub fn with_ms_level(mut self, level: Option<u8>) -> Self {
        self.ms_level_filter = level;
        self
    }

    /// Number of spectra in the file (all MS levels).
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the file has no spectra.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn convert(raw: &RawSpectrum) -> Spectrum {
        let scan_number = raw.index() as i32 + 1;

        let (precursor_mz, precursor_charge, precursor_intensity) = match raw.precursor() {
            Some(p) => {
                let charge = match p.charge() {
                    0 => None,
                    z => Some(z),
                };
                let inten = p.intensity();
                let inten = if inten > 0.0 { Some(inten) } else { None };
                (p.mz(), charge, inten)
            }
            None => (0.0, None, None),
        };

        // `data_raw()` exposes the FlatBuffers vectors directly; `.mz()` /
        // `.intensity()` are `Option<Vector<_>>`, and a `Vector`'s `iter()`
        // already yields owned scalars (no `.copied()`).
        let peaks: Vec<(f64, f32)> = match raw.data_raw() {
            Some(data) => match (data.mz(), data.intensity()) {
                (Some(mz), Some(intensity)) => mz.iter().zip(intensity.iter()).collect(),
                _ => Vec::new(),
            },
            None => Vec::new(),
        };

        Spectrum {
            // The Thermo controller native id (e.g.
            // `controllerType=0 controllerNumber=1 scan=N`), matching the
            // `<spectrum id=...>` an mzML produced from the same `.raw` carries,
            // so the PIN `SpecID` column lines up across the two input formats.
            title: raw.native_id(),
            precursor_mz,
            precursor_intensity,
            precursor_charge,
            // RawFileReader retention time is in minutes; the model uses seconds.
            rt_seconds: Some(raw.time() * 60.0),
            scan: Some(scan_number),
            peaks,
            // Activation + isolation window are mapped in a later milestone
            // (the CLI accepts `--fragmentation`/`--instrument`; the chimeric
            // cascade is the only consumer of the isolation window).
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }
}

impl Iterator for ThermoRawReader {
    type Item = Result<Spectrum, ThermoParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next < self.len {
            let i = self.next;
            self.next += 1;
            let raw = match self.handle.get(i) {
                Some(raw) => raw,
                None => continue,
            };
            if let Some(want) = self.ms_level_filter {
                if raw.ms_level() != want {
                    continue;
                }
            }
            return Some(Ok(Self::convert(&raw)));
        }
        None
    }
}
