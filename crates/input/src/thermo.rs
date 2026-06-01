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

use crate::{Ms1Link, Spectrum};

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

        let mut precursor_mz = 0.0;
        let mut precursor_charge = None;
        let mut precursor_intensity = None;
        let mut isolation_lower_offset = None;
        let mut isolation_upper_offset = None;
        if let Some(p) = raw.precursor() {
            precursor_mz = p.mz();
            precursor_charge = match p.charge() {
                0 => None,
                z => Some(z),
            };
            let inten = p.intensity();
            precursor_intensity = if inten > 0.0 { Some(inten) } else { None };
            // Isolation window: RawFileReader gives absolute m/z bounds
            // (lower, target, upper); the model stores offsets from the target,
            // matching mzML's `isolation window lower/upper offset`. Skip a
            // degenerate/absent window (the cascade is its only consumer).
            let w = p.isolation_window();
            let (lo, tg, up) = (w.lower(), w.target(), w.upper());
            if tg > 0.0 && up > lo {
                isolation_lower_offset = Some((tg - lo).abs());
                isolation_upper_offset = Some((up - tg).abs());
            }
        }

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
            peaks: extract_peaks(raw),
            // Activation is left to the CLI (`--fragmentation`/`--instrument`).
            activation_method: None,
            isolation_lower_offset,
            isolation_upper_offset,
        }
    }

    /// Stream MS2 in `chunk_size` batches, each paired with a bounded
    /// [`Ms1Link`] (the MS1 scans seen in that chunk plus a per-MS2 link to the
    /// most-recent preceding MS1). Mirrors `MzMLReader::read_with_ms1_chunked`
    /// so the chimeric cascade consumes `.raw` exactly like mzML. The last MS1
    /// of a chunk is carried into the next chunk (index 0) so the first MS2 of
    /// a chunk still links to its true preceding MS1 across the boundary.
    ///
    /// Returns `(error_count, sample_errors)` for API symmetry; `.raw` reads do
    /// not surface per-scan parse errors, so these are always `(0, vec![])`.
    pub fn read_with_ms1_chunked<F>(
        self,
        chunk_size: usize,
        cap: usize,
        mut on_chunk: F,
    ) -> (usize, Vec<String>)
    where
        F: FnMut(Vec<Spectrum>, Ms1Link),
    {
        let mut chunk: Vec<Spectrum> = Vec::new();
        let mut ms1_peaks: Vec<Vec<(f64, f32)>> = Vec::new();
        let mut ms2_to_ms1: Vec<Option<usize>> = Vec::new();
        let mut current_ms1: Option<usize> = None;
        let mut emitted_ms2 = 0usize;

        for i in 0..self.len {
            if cap > 0 && emitted_ms2 >= cap {
                break;
            }
            let raw = match self.handle.get(i) {
                Some(raw) => raw,
                None => continue,
            };
            match raw.ms_level() {
                1 => {
                    let idx = ms1_peaks.len();
                    ms1_peaks.push(extract_peaks(&raw));
                    current_ms1 = Some(idx);
                }
                2 => {
                    chunk.push(Self::convert(&raw));
                    ms2_to_ms1.push(current_ms1);
                    emitted_ms2 += 1;
                    if chunk.len() >= chunk_size {
                        // Carry the current MS1's peaks into the next chunk so
                        // its first MS2 still links across the boundary.
                        let carry = current_ms1.map(|idx| ms1_peaks[idx].clone());
                        let link = Ms1Link {
                            ms1_peaks: std::mem::take(&mut ms1_peaks),
                            ms2_to_ms1: std::mem::take(&mut ms2_to_ms1),
                        };
                        on_chunk(std::mem::take(&mut chunk), link);
                        current_ms1 = carry.map(|p| {
                            ms1_peaks.push(p);
                            0
                        });
                    }
                }
                _ => {}
            }
        }

        if !chunk.is_empty() {
            on_chunk(chunk, Ms1Link { ms1_peaks, ms2_to_ms1 });
        }
        (0, Vec::new())
    }
}

/// Extract centroided peaks `(m/z, intensity)` from a raw spectrum. The
/// FlatBuffers vectors yield owned scalars, so no `.copied()` is needed.
fn extract_peaks(raw: &RawSpectrum) -> Vec<(f64, f32)> {
    match raw.data_raw() {
        Some(data) => match (data.mz(), data.intensity()) {
            (Some(mz), Some(intensity)) => mz.iter().zip(intensity.iter()).collect(),
            _ => Vec::new(),
        },
        None => Vec::new(),
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
