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

use std::collections::HashMap;

use thermorawfilereader::schema::{DissociationMethod, MassAnalyzer};
use thermorawfilereader::{RawFileReader, RawSpectrum};

use model::{ActivationMethod, InstrumentType};

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

    /// Convert one MS2 `RawSpectrum` to a [`Spectrum`], or `None` if the scan
    /// has no usable precursor m/z — matching the mzML reader, which skips
    /// precursor-less MS2 rather than forwarding a `precursor_mz = 0.0` spectrum
    /// that would search against a nonsense mass window.
    fn convert(raw: &RawSpectrum) -> Option<Spectrum> {
        let scan_number = raw.index() as i32 + 1;

        // A precursor record with a positive m/z is required to search an MS2.
        let p = raw.precursor()?;
        let precursor_mz = p.mz();
        if !(precursor_mz > 0.0) {
            return None;
        }
        let precursor_charge = match p.charge() {
            0 => None,
            z => Some(z),
        };
        let inten = p.intensity();
        let precursor_intensity = if inten > 0.0 { Some(inten) } else { None };

        // Isolation window: RawFileReader gives absolute m/z bounds
        // (lower, target, upper); the model stores offsets from the target,
        // matching mzML's `isolation window lower/upper offset`. Skip a
        // degenerate/absent window (the cascade is its only consumer).
        let w = p.isolation_window();
        let (lo, tg, up) = (w.lower(), w.target(), w.upper());
        let (isolation_lower_offset, isolation_upper_offset) = if tg > 0.0 && up > lo {
            (Some((tg - lo).abs()), Some((up - tg).abs()))
        } else {
            (None, None)
        };

        Some(Spectrum {
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
            // Real activation from the .raw, used by the param auto-router so a
            // CID/ETD .raw is not silently scored with the HCD model.
            activation_method: map_dissociation(p.activation().dissociation_method()),
            isolation_lower_offset,
            isolation_upper_offset,
        })
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
        let mut skipped_higher_ms = 0usize;

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
                // MS3+ scans (e.g. TMT SPS-MS3 reporter-quant scans) are NOT
                // identification spectra — only MS2 is searched and only MS1
                // links it, in chimeric OR normal mode. Skipping them here keeps
                // the cascade from ever scoring an MS3 in TMT-MS3 acquisitions.
                level if level >= 3 => {
                    skipped_higher_ms += 1;
                }
                2 => {
                    // Skip MS2 with no usable precursor (mirrors mzML); such a
                    // scan must not consume a chunk slot or an Ms1Link entry.
                    let spec = match Self::convert(&raw) {
                        Some(spec) => spec,
                        None => continue,
                    };
                    chunk.push(spec);
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
        if skipped_higher_ms > 0 {
            eprintln!(
                "Thermo .raw: skipped {skipped_higher_ms} MS3+ scans (not searched; \
                 e.g. TMT SPS-MS3 reporter-quant scans)"
            );
        }
        (0, Vec::new())
    }
}

/// Extract centroided peaks `(m/z, intensity)` from a raw spectrum, ascending
/// by m/z. The FlatBuffers vectors yield owned scalars, so no `.copied()`.
///
/// Downstream consumers require m/z-ascending peaks: the chimeric co-isolation
/// detection binary-searches the list (`partition_point`), and the scorer/`Spectrum`
/// contract documents the invariant. RawFileReader centroid streams are normally
/// already sorted, so the list is only sorted when an inversion is present.
fn extract_peaks(raw: &RawSpectrum) -> Vec<(f64, f32)> {
    let mut peaks: Vec<(f64, f32)> = match raw.data_raw() {
        Some(data) => match (data.mz(), data.intensity()) {
            (Some(mz), Some(intensity)) => mz.iter().zip(intensity.iter()).collect(),
            _ => Vec::new(),
        },
        None => Vec::new(),
    };
    if peaks.windows(2).any(|w| w[0].0 > w[1].0) {
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    peaks
}

/// Map a Thermo `DissociationMethod` to the model's [`ActivationMethod`].
/// Electron-based variants (ETD/ETCID/ETHCD/ECD/…) collapse to `ETD`; unknown
/// or unsupported methods return `None` (the resolver then keeps its default).
fn map_dissociation(d: DissociationMethod) -> Option<ActivationMethod> {
    match d.0 {
        1 => Some(ActivationMethod::CID),
        2 => Some(ActivationMethod::HCD),
        // ETD, ETCID, ETHCD, ECD, ECCID, ECHCD, NETD
        4 | 5 | 6 | 8 | 9 | 10 | 16 => Some(ActivationMethod::ETD),
        _ => None,
    }
}

/// Map a Thermo `MassAnalyzer` to the model's [`InstrumentType`]. FTMS
/// (Orbitrap) maps to `QExactive`; ASTMS (Astral) maps to the explicit
/// `OrbitrapAstral` class; TOFMS is `TOF`; ITMS (ion trap) is `LowRes`.
/// Others return `None`.
fn map_analyzer(a: MassAnalyzer) -> Option<InstrumentType> {
    match a.0 {
        1 => Some(InstrumentType::LowRes),         // ITMS (ion trap)
        4 => Some(InstrumentType::TOF),            // TOFMS
        5 => Some(InstrumentType::QExactive),      // FTMS (Orbitrap)
        7 => Some(InstrumentType::OrbitrapAstral), // ASTMS (Astral)
        _ => None,
    }
}

/// Peek the first `peek_n` MS2 scans of a `.raw` and return the dominant
/// `(activation, instrument)` for param auto-routing, mirroring the mzML
/// `detect_dominant_activation` / `detect_instrument_type_for_path` helpers so a
/// CID/ETD/ion-trap `.raw` is not silently scored with the HCD/QExactive model.
/// Returns `None` if the file can't be opened or no MS2 activation is found.
pub fn detect_activation_instrument<P: AsRef<Path>>(
    path: P,
    peek_n: usize,
) -> Option<(ActivationMethod, Option<InstrumentType>)> {
    let reader = RawFileReader::open(path.as_ref()).ok()?;
    let mut act: HashMap<u8, usize> = HashMap::new();
    let mut inst: HashMap<u8, usize> = HashMap::new();
    let mut seen = 0usize;
    for i in 0..reader.len() {
        if seen >= peek_n {
            break;
        }
        let raw = match reader.get(i) {
            Some(r) => r,
            None => continue,
        };
        if raw.ms_level() != 2 {
            continue;
        }
        seen += 1;
        if let Some(p) = raw.precursor() {
            if let Some(m) = map_dissociation(p.activation().dissociation_method()) {
                *act.entry(m as u8).or_insert(0) += 1;
            }
        }
        if let Some(a) = raw.acquisition() {
            if let Some(t) = map_analyzer(a.mass_analyzer()) {
                *inst.entry(t as u8).or_insert(0) += 1;
            }
        }
    }
    let dominant_act = act.into_iter().max_by_key(|&(_, n)| n).map(|(k, _)| k)?;
    let activation = match dominant_act {
        x if x == ActivationMethod::CID as u8 => ActivationMethod::CID,
        x if x == ActivationMethod::HCD as u8 => ActivationMethod::HCD,
        x if x == ActivationMethod::ETD as u8 => ActivationMethod::ETD,
        x if x == ActivationMethod::PQD as u8 => ActivationMethod::PQD,
        _ => ActivationMethod::UVPD,
    };
    let instrument = inst
        .into_iter()
        .max_by_key(|&(_, n)| n)
        .map(|(k, _)| match k {
            x if x == InstrumentType::LowRes as u8         => InstrumentType::LowRes,
            x if x == InstrumentType::TOF as u8            => InstrumentType::TOF,
            x if x == InstrumentType::HighRes as u8        => InstrumentType::HighRes,
            x if x == InstrumentType::OrbitrapAstral as u8 => InstrumentType::OrbitrapAstral,
            x if x == InstrumentType::TimsTOF as u8        => InstrumentType::TimsTOF,
            _ => InstrumentType::QExactive,
        });
    Some((activation, instrument))
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
            match Self::convert(&raw) {
                // Skip MS2 scans with no usable precursor (mirrors mzML).
                Some(spec) => return Some(Ok(spec)),
                None => continue,
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dissociation_maps_to_activation() {
        assert_eq!(map_dissociation(DissociationMethod(1)), Some(ActivationMethod::CID));
        assert_eq!(map_dissociation(DissociationMethod(2)), Some(ActivationMethod::HCD));
        // ETD / ETCID / ETHCD / ECD-family / NETD all collapse to ETD.
        for v in [4u8, 5, 6, 8, 9, 10, 16] {
            assert_eq!(map_dissociation(DissociationMethod(v)), Some(ActivationMethod::ETD));
        }
        assert_eq!(map_dissociation(DissociationMethod(0)), None); // Unknown
    }

    #[test]
    fn analyzer_maps_to_instrument() {
        assert_eq!(map_analyzer(MassAnalyzer(1)), Some(InstrumentType::LowRes));         // ITMS
        assert_eq!(map_analyzer(MassAnalyzer(4)), Some(InstrumentType::TOF));            // TOFMS
        assert_eq!(map_analyzer(MassAnalyzer(5)), Some(InstrumentType::QExactive));      // FTMS (Orbitrap)
        assert_eq!(map_analyzer(MassAnalyzer(7)), Some(InstrumentType::OrbitrapAstral)); // ASTMS (Astral)
        assert_eq!(map_analyzer(MassAnalyzer(0)), None); // Unknown
    }
}
