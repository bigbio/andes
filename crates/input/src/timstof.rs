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
//! Scope: DDA-PASEF. Each precursor of the run becomes one centroided MS2
//! fragment spectrum, summed over all its PASEF selections and converted to
//! m/z with the run's calibration (see [`crate::timstof_ms2`]); it carries the
//! precursor's m/z, charge, retention time, parent MS1 frame and ion mobility.
//! MS2 spectra are built in parallel batches. On request the reader also yields the run's MS1 survey frames,
//! summed over ion mobility and centroided (see [`crate::timstof_ms1`]), with
//! the same capture interface as the mzML and Thermo readers: every MS1 scan
//! for quantification ([`TimsTofReader::with_run_ms1_capture`]), the
//! acquisition context of each MS2 with its parent frame and precursor
//! mobility ([`TimsTofReader::with_scan_context_capture`]), and MS2 chunks
//! linked to their parent MS1 for the MS1-linked search paths
//! ([`TimsTofReader::read_with_ms1_chunked`]). MS1 frames are decoded in
//! parallel batches as the MS2 stream reaches them, so only one batch of
//! frames is held besides the captured centroids. A TDF holds MS1 and PASEF
//! MS2 frames only; there are no MS3 scans.

use std::collections::{HashMap, VecDeque};
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use timsrust::converters::{ConvertableDomain, Scan2ImConverter};
use timsrust::readers::{FrameReader, FrameReaderError, MetadataReader};
use timsrust::MSLevel;

use model::scan::{Ms1Scan, ScanContext};

use crate::timstof_ms1::{
    centroid_profile, fit_calibration, match_precursors, mode_offset_ppm, sum_over_mobility,
    CalibrationFit, FrameConditions, MzCalibration, MIN_REFINE_WINDOW_PPM, NARROW_WINDOW_PPM,
    REFINE_ROUNDS, REFINE_SIGMAS, WIDE_WINDOW_PPM,
};
use crate::timstof_ms2::{
    ms2_peaks, push_scans, read_frame_rows, read_pasef_precursors, FrameRow, PasefPrecursor,
};
use crate::{Ms1Link, Spectrum};

/// Error opening or reading a Bruker timsTOF `.d` directory.
#[derive(Debug, thiserror::Error)]
pub enum TimsTofParseError {
    #[error("failed to open Bruker .d '{path}': {source}{hint}")]
    Open {
        path: String,
        source: Box<FrameReaderError>,
        /// Extra guidance for a known unsupported layout (empty otherwise).
        hint: &'static str,
    },
    #[error(
        "'{path}' is not a Bruker .d run folder: it has no analysis.tdf, but the folder \
         '{nested}' inside it has one; pass that folder"
    )]
    Nested { path: String, nested: String },
    #[error("failed to read the MS2 frames of precursor {precursor} from Bruker .d: {message}")]
    Read { precursor: usize, message: String },
    #[error("failed to open the MS1 frames of Bruker .d '{path}': {source}")]
    Frames {
        path: String,
        source: Box<FrameReaderError>,
    },
    #[error("failed to read the precursor and frame tables of Bruker .d '{path}': {source}")]
    Tables {
        path: String,
        source: rusqlite::Error,
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

/// MS1 frames decoded per parallel batch.
const MS1_BATCH: usize = 64;
/// MS1 frames sampled across the run for the m/z calibration.
const CALIBRATION_FRAMES: usize = 400;

/// The MS1 survey frames of a `.d`, decoded lazily in file order.
struct Ms1Frames {
    frames: FrameReader,
    /// Frame-reader positions and frame ids of the MS1 frames, in file order.
    positions: Vec<usize>,
    ids: Vec<usize>,
    /// Frame id → `(T1, T2)` temperatures.
    temperatures: HashMap<usize, [f64; 2]>,
    calibration: MzCalibration,
    /// Next entry of `positions` to decode.
    cursor: usize,
    /// Decoded frames not yet handed out: frame id and scan.
    ready: VecDeque<(usize, Ms1Scan)>,
}

/// The `Frames` temperatures and the vendor precursor peaks of a `.d`.
struct TdfTables {
    temperatures: HashMap<usize, [f64; 2]>,
    /// Parent frame id → monoisotopic m/z of its precursors
    /// (`MonoisotopicMz`, else `LargestPeakMz` when the monoisotopic peak was
    /// not determined).
    precursor_peaks: HashMap<usize, Vec<f64>>,
}

fn read_tables(dir: &Path) -> Result<TdfTables, rusqlite::Error> {
    let conn = rusqlite::Connection::open_with_flags(
        dir.join("analysis.tdf"),
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )?;
    let mut temperatures = HashMap::new();
    // Older schemas may lack the temperature columns; the calibration then
    // has no temperature terms.
    if let Ok(mut stmt) = conn.prepare("SELECT Id, T1, T2 FROM Frames") {
        let rows = stmt.query_map([], |r| {
            Ok((
                r.get::<_, i64>(0)?,
                r.get::<_, Option<f64>>(1)?,
                r.get::<_, Option<f64>>(2)?,
            ))
        })?;
        for row in rows {
            let (id, t1, t2) = row?;
            temperatures.insert(id as usize, [t1.unwrap_or(0.0), t2.unwrap_or(0.0)]);
        }
    }
    let mut precursor_peaks: HashMap<usize, Vec<f64>> = HashMap::new();
    let mut stmt = conn.prepare("SELECT Parent, MonoisotopicMz, LargestPeakMz FROM Precursors")?;
    let rows = stmt.query_map([], |r| {
        Ok((
            r.get::<_, Option<i64>>(0)?,
            r.get::<_, Option<f64>>(1)?,
            r.get::<_, Option<f64>>(2)?,
        ))
    })?;
    for row in rows {
        let (parent, mono, largest) = row?;
        let mz = mono.filter(|m| m.is_finite() && *m > 0.0).or(largest);
        if let (Some(parent), Some(mz)) = (parent, mz) {
            precursor_peaks.entry(parent as usize).or_default().push(mz);
        }
    }
    Ok(TdfTables {
        temperatures,
        precursor_peaks,
    })
}

/// The TOF-domain centroids of one decoded frame.
fn frame_centroids(frame: &timsrust::Frame) -> Vec<crate::timstof_ms1::TofCentroid> {
    centroid_profile(&sum_over_mobility(&frame.tof_indices, &frame.intensities))
}

/// Centroids converted to `(m/z, intensity)`, m/z ascending.
fn to_peaks(
    centroids: &[crate::timstof_ms1::TofCentroid],
    calibration: &MzCalibration,
    frame: FrameConditions,
) -> Vec<(f64, f32)> {
    let mut peaks: Vec<(f64, f32)> = centroids
        .iter()
        .map(|c| {
            let cond = crate::timstof_ms1::conditions(frame, c.intensity);
            (calibration.mz(c.tof, cond), c.intensity as f32)
        })
        .filter(|&(mz, i)| mz.is_finite() && mz > 0.0 && i.is_finite() && i > 0.0)
        .collect();
    if peaks.windows(2).any(|w| w[0].0 > w[1].0) {
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    peaks
}

impl Ms1Frames {
    fn open(dir: &Path) -> Result<Self, TimsTofParseError> {
        let frames_err = |source| TimsTofParseError::Frames {
            path: dir.display().to_string(),
            source: Box::new(source),
        };
        let frames = FrameReader::new(dir).map_err(frames_err)?;
        let metadata = MetadataReader::new(dir).map_err(|e| frames_err(e.into()))?;
        let tables = read_tables(dir).map_err(|source| TimsTofParseError::Tables {
            path: dir.display().to_string(),
            source,
        })?;
        let mut positions = Vec::new();
        let mut ids = Vec::new();
        for pos in 0..frames.len() {
            let f = frames
                .get_frame_without_coordinates(pos)
                .map_err(frames_err)?;
            if f.ms_level == MSLevel::MS1 {
                positions.push(pos);
                ids.push(f.index);
            }
        }
        // The acquisition-range model of timsrust as `sqrt(m/z) = a + b * tof`.
        let conv = metadata.mz_converter;
        let intercept = conv.convert(0.0).sqrt();
        let slope = conv.convert(1.0).sqrt() - intercept;
        let scale = conv.invert(metadata.upper_mz).max(1.0);
        let base = MzCalibration::linear(intercept, slope, scale);
        let mut ms1 = Self {
            frames,
            positions,
            ids,
            temperatures: tables.temperatures,
            calibration: base.clone(),
            cursor: 0,
            ready: VecDeque::new(),
        };
        match ms1.calibrate(&base, &tables.precursor_peaks) {
            Some(fit) => {
                eprintln!(
                    "Bruker .d: MS1 m/z calibration from {} precursor peaks: residual to the \
                     vendor m/z median {:+.2} ppm, robust SD {:.2} ppm (acquisition-range \
                     model: {:+.2} ppm, SD {:.2} ppm)",
                    fit.n_used,
                    fit.median_ppm,
                    fit.sigma_ppm,
                    fit.base_median_ppm,
                    fit.base_sigma_ppm
                );
                ms1.calibration = fit.calibration;
            }
            None => eprintln!(
                "WARN: Bruker .d: too few precursor peaks matched their MS1 frame for an m/z \
                 calibration; MS1 m/z uses the acquisition-range model"
            ),
        }
        Ok(ms1)
    }

    fn temperature(&self, id: usize) -> [f64; 2] {
        self.temperatures.get(&id).copied().unwrap_or([0.0; 2])
    }

    /// Fit the m/z calibration on up to [`CALIBRATION_FRAMES`] MS1 frames
    /// spread over the run that have precursors.
    fn calibrate(
        &self,
        base: &MzCalibration,
        peaks: &HashMap<usize, Vec<f64>>,
    ) -> Option<CalibrationFit> {
        let with_precursors: Vec<usize> = (0..self.positions.len())
            .filter(|&k| peaks.contains_key(&self.ids[k]))
            .collect();
        let stride = (with_precursors.len() / CALIBRATION_FRAMES).max(1);
        let sample: Vec<usize> = with_precursors.into_iter().step_by(stride).collect();
        type Decoded = (usize, FrameConditions, Vec<crate::timstof_ms1::TofCentroid>);
        let decoded: Vec<Decoded> = sample
            .par_iter()
            .filter_map(|&k| {
                let id = self.ids[k];
                let frame = self.frames.get(self.positions[k]).ok()?;
                let cond = FrameConditions {
                    t: self.temperature(id),
                    rt_s: frame.rt_in_seconds,
                };
                Some((id, cond, frame_centroids(&frame)))
            })
            .collect();
        let pass = |calibration: &MzCalibration, offset: f64, window: f64| {
            let mut points = Vec::new();
            for (id, cond, centroids) in &decoded {
                match_precursors(
                    centroids,
                    &peaks[id],
                    *cond,
                    calibration,
                    offset,
                    window,
                    &mut points,
                );
            }
            points
        };
        let wide = pass(base, 0.0, WIDE_WINDOW_PPM);
        let offset = mode_offset_ppm(&wide, base)?;
        let mut fit = fit_calibration(&pass(base, offset, NARROW_WINDOW_PPM), base)?;
        for _ in 0..REFINE_ROUNDS {
            let window = (REFINE_SIGMAS * fit.sigma_ppm).max(MIN_REFINE_WINDOW_PPM);
            match fit_calibration(&pass(&fit.calibration, 0.0, window), base) {
                Some(next) => fit = next,
                None => break,
            }
        }
        Some(fit)
    }

    /// Decode the next batch of MS1 frames in parallel.
    fn decode_batch(&mut self) {
        let end = (self.cursor + MS1_BATCH).min(self.positions.len());
        let batch: Vec<(usize, usize)> = (self.cursor..end)
            .map(|k| (self.positions[k], self.ids[k]))
            .collect();
        self.cursor = end;
        let decoded: Vec<(usize, Ms1Scan)> = batch
            .par_iter()
            .map(|&(pos, id)| {
                let (rt, peaks) = match self.frames.get(pos) {
                    Ok(frame) => (
                        frame.rt_in_seconds,
                        to_peaks(
                            &frame_centroids(&frame),
                            &self.calibration,
                            FrameConditions {
                                t: self.temperature(id),
                                rt_s: frame.rt_in_seconds,
                            },
                        ),
                    ),
                    // A frame that cannot be decoded is kept as an empty scan
                    // so the survey series keeps its time points.
                    Err(_) => (
                        self.frames
                            .get_frame_without_coordinates(pos)
                            .map(|f| f.rt_in_seconds)
                            .unwrap_or(f64::NAN),
                        Vec::new(),
                    ),
                };
                (
                    id,
                    Ms1Scan {
                        rt,
                        peaks,
                        id: ms1_native_id(id),
                        faims_cv: None,
                    },
                )
            })
            .collect();
        self.ready.extend(decoded);
    }

    /// Hand the decoded MS1 frames with an id up to `frame_id` (all remaining
    /// frames for `None`) to `sink` with their frame id, in file order.
    fn take_up_to(&mut self, frame_id: Option<usize>, mut sink: impl FnMut(usize, Ms1Scan)) {
        loop {
            while let Some((id, _)) = self.ready.front() {
                if frame_id.is_some_and(|limit| *id > limit) {
                    return;
                }
                let (id, scan) = self.ready.pop_front().expect("front exists");
                sink(id, scan);
            }
            if self.cursor >= self.positions.len() {
                return;
            }
            self.decode_batch();
        }
    }
}

/// An emitted MS2, the id of its parent MS1 frame and its precursor 1/K0.
type Ms2Item = (Spectrum, Option<usize>, f64);

/// Native id of an MS1 frame: `frame=<TDF frame id>`.
fn ms1_native_id(frame_id: usize) -> String {
    format!("frame={frame_id}")
}

/// Precursors whose MS2 spectra are built per parallel batch.
const MS2_BATCH: usize = 1024;

/// The PASEF MS2 side of a `.d`: the precursors with their selections and the
/// frames they were fragmented in.
struct Ms2Source {
    frames: FrameReader,
    /// Frame id → position in `frames`.
    positions: HashMap<usize, usize>,
    rows: HashMap<usize, FrameRow>,
    precursors: Vec<PasefPrecursor>,
    im: Scan2ImConverter,
}

impl Ms2Source {
    fn open(dir: &Path) -> Result<Self, TimsTofParseError> {
        let open_err = |source: FrameReaderError| {
            let hint = open_error_hint(&source.to_string());
            TimsTofParseError::Open {
                path: dir.display().to_string(),
                source: Box::new(source),
                hint,
            }
        };
        let frames = FrameReader::new(dir).map_err(open_err)?;
        let metadata = MetadataReader::new(dir).map_err(|e| open_err(e.into()))?;
        let tables_err = |source| TimsTofParseError::Tables {
            path: dir.display().to_string(),
            source,
        };
        let conn = rusqlite::Connection::open_with_flags(
            dir.join("analysis.tdf"),
            rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
        )
        .map_err(tables_err)?;
        let precursors = read_pasef_precursors(&conn).map_err(tables_err)?;
        let rows = read_frame_rows(&conn).map_err(tables_err)?;
        let mut positions = HashMap::with_capacity(frames.len());
        for pos in 0..frames.len() {
            let f = frames
                .get_frame_without_coordinates(pos)
                .map_err(open_err)?;
            positions.insert(f.index, pos);
        }
        Ok(Self {
            frames,
            positions,
            rows,
            precursors,
            im: metadata.im_converter,
        })
    }

    /// Build the MS2 spectra of `precursors[range]`: the frames they need are
    /// decoded once each, in parallel, then the spectra are assembled in
    /// parallel. Precursors without a usable m/z yield nothing.
    fn build(
        &self,
        range: std::ops::Range<usize>,
        calibration: &MzCalibration,
    ) -> Vec<Result<Ms2Item, TimsTofParseError>> {
        let batch = &self.precursors[range];
        let mut needed: Vec<usize> = batch
            .iter()
            .flat_map(|p| p.selections.iter().map(|s| s.frame_id))
            .filter(|id| self.rows.get(id).is_none_or(|r| r.num_peaks > 0))
            .collect();
        needed.sort_unstable();
        needed.dedup();
        let decoded: HashMap<usize, Result<timsrust::Frame, String>> = needed
            .par_iter()
            .map(|&id| {
                let frame = match self.positions.get(&id) {
                    Some(&pos) => self.frames.get(pos).map_err(|e| e.to_string()),
                    None => Err(format!("frame {id} is not in the frame table")),
                };
                (id, frame)
            })
            .collect();
        batch
            .par_iter()
            .filter_map(|p| self.assemble(p, &decoded, calibration))
            .collect()
    }

    /// One precursor's spectrum from its decoded frames (see
    /// [`crate::timstof_ms2`]).
    fn assemble(
        &self,
        p: &PasefPrecursor,
        decoded: &HashMap<usize, Result<timsrust::Frame, String>>,
        calibration: &MzCalibration,
    ) -> Option<Result<Ms2Item, TimsTofParseError>> {
        let mut tof = Vec::new();
        let mut intensity = Vec::new();
        let mut n_frames = 0usize;
        let (mut t, mut rt) = ([0.0f64; 2], 0.0f64);
        let mut error: Option<String> = None;
        for s in &p.selections {
            let Some(result) = decoded.get(&s.frame_id) else {
                // A frame recorded with no peaks contributes nothing.
                continue;
            };
            match result {
                Ok(frame) => {
                    push_scans(
                        &frame.scan_offsets,
                        &frame.tof_indices,
                        &frame.intensities,
                        s.scan_begin,
                        s.scan_end,
                        &mut tof,
                        &mut intensity,
                    );
                    let row = self.rows.get(&s.frame_id);
                    let ft = row.map_or([0.0; 2], |r| r.t);
                    t[0] += ft[0];
                    t[1] += ft[1];
                    rt += frame.rt_in_seconds;
                    n_frames += 1;
                }
                Err(message) => {
                    error.get_or_insert_with(|| message.clone());
                }
            }
        }
        if n_frames == 0 {
            if let Some(message) = error {
                return Some(Err(TimsTofParseError::Read {
                    precursor: p.id,
                    message,
                }));
            }
        }
        let parent_rt = p
            .parent
            .and_then(|id| self.rows.get(&id))
            .map(|r| r.time)
            .filter(|t| t.is_finite());
        let peaks = if n_frames > 0 {
            let n = n_frames as f64;
            let cond = FrameConditions {
                t: [t[0] / n, t[1] / n],
                rt_s: rt / n,
            };
            ms2_peaks(&tof, &intensity, calibration, cond, n_frames)
        } else {
            Vec::new()
        };
        let rt_seconds = parent_rt.or_else(|| {
            p.selections
                .first()
                .and_then(|s| self.rows.get(&s.frame_id))
                .map(|r| r.time)
        });
        let spec = to_spectrum(p, rt_seconds, peaks)?;
        let im = p.scan_number.map_or(f64::NAN, |s| self.im.convert(s));
        Some(Ok((spec, p.parent, im)))
    }
}

/// The [`Spectrum`] of a precursor, or `None` without a positive precursor m/z
/// (the mzML and Thermo readers skip such scans too).
// `!(mz > 0.0)` deliberately rejects NaN as well as 0 / negative m/z — the
// same precursor-sanity idiom the Thermo reader uses.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn to_spectrum(
    p: &PasefPrecursor,
    rt_seconds: Option<f64>,
    peaks: Vec<(f64, f32)>,
) -> Option<Spectrum> {
    if !(p.mz > 0.0) {
        return None;
    }
    let scan_number = i32::try_from(p.id).ok()?;
    // 0 (or an absurd value) means "unknown"; the engine then sweeps the
    // configured charge range.
    let precursor_charge = p.charge.and_then(model::spectrum::known_precursor_charge);
    let precursor_intensity = p
        .intensity
        .filter(|i| i.is_finite() && *i > 0.0)
        .map(|i| i as f32);
    // The quadrupole isolation window is recorded as a centre m/z and a total
    // width; the model stores symmetric lower / upper offsets.
    let (isolation_lower_offset, isolation_upper_offset) = if p.isolation_width > 0.0 {
        let half = p.isolation_width / 2.0;
        (Some(half), Some(half))
    } else {
        (None, None)
    };
    Some(Spectrum {
        // The precursor id identifies the spectrum within the `.d` and becomes
        // the PIN `SpecID` / scan number.
        title: format!("scan={scan_number}"),
        precursor_mz: p.mz,
        precursor_intensity,
        precursor_charge,
        // Retention time of the parent MS1 frame, in seconds.
        rt_seconds,
        scan: Some(scan_number),
        peaks,
        // timsTOF DDA-PASEF is CID/HCD-style beam-type fragmentation; the
        // `.d` does not record a discrete activation cvParam the resolver
        // keys on, so this stays `None` and the param resolver uses its
        // default (or an explicit `--fragmentation`/`--instrument`).
        activation_method: None,
        isolation_lower_offset,
        isolation_upper_offset,
        isolation_target_mz: None,
    })
}

/// Reader over a Bruker timsTOF `.d` directory, yielding [`Spectrum`] values.
///
/// Every emitted spectrum is a PASEF MS2 fragment spectrum, one per precursor
/// in precursor-id order, so there is no MS-level filter. MS1 survey frames
/// are read alongside on request.
pub struct TimsTofReader {
    ms2: Ms2Source,
    dir: PathBuf,
    /// Next precursor (index into the precursor list) to build.
    next: usize,
    /// Built MS2 not yet handed out.
    ready: VecDeque<Result<Ms2Item, TimsTofParseError>>,
    /// The run's TOF → m/z calibration, fitted on the first MS2 request.
    calibration: Option<MzCalibration>,
    /// Quantification capture (see the mzML reader): every MS1 frame of the
    /// run, and the acquisition context of every emitted MS2.
    capture_run_ms1: bool,
    run_ms1: Vec<Ms1Scan>,
    capture_scan_context: bool,
    scan_contexts: Vec<ScanContext>,
    /// Opened on the first MS1 or MS2 request.
    ms1: Option<Ms1Frames>,
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
        let ms2 = Ms2Source::open(p)?;
        Ok(Self {
            ms2,
            dir: p.to_path_buf(),
            next: 0,
            ready: VecDeque::new(),
            calibration: None,
            capture_run_ms1: false,
            run_ms1: Vec::new(),
            capture_scan_context: false,
            scan_contexts: Vec::new(),
            ms1: None,
        })
    }

    /// Keep every MS1 frame (retention time + centroids summed over ion
    /// mobility) for label-free quantification; take them with
    /// [`Self::take_run_ms1`] once the stream is exhausted.
    pub fn with_run_ms1_capture(mut self, capture: bool) -> Self {
        self.capture_run_ms1 = capture;
        self
    }

    /// Record the acquisition context of every emitted MS2 (parent MS1 frame
    /// and precursor ion mobility); take them with [`Self::take_scan_contexts`].
    pub fn with_scan_context_capture(mut self, capture: bool) -> Self {
        self.capture_scan_context = capture;
        self
    }

    pub fn take_run_ms1(&mut self) -> Vec<Ms1Scan> {
        std::mem::take(&mut self.run_ms1)
    }

    pub fn take_scan_contexts(&mut self) -> Vec<ScanContext> {
        std::mem::take(&mut self.scan_contexts)
    }

    /// Number of MS2 spectra in the `.d` (precursors with a PASEF selection).
    pub fn len(&self) -> usize {
        self.ms2.precursors.len()
    }

    /// Whether the `.d` has no spectra.
    pub fn is_empty(&self) -> bool {
        self.ms2.precursors.is_empty()
    }

    fn ms1_frames(&mut self) -> Result<&mut Ms1Frames, TimsTofParseError> {
        if self.ms1.is_none() {
            self.ms1 = Some(Ms1Frames::open(&self.dir)?);
        }
        Ok(self.ms1.as_mut().expect("just opened"))
    }

    /// The next MS2 with its parent MS1 frame id and precursor ion mobility,
    /// or `None` at the end of the run. The m/z calibration is fitted on the
    /// first call (from the MS1 frames, see [`Ms1Frames::open`]).
    fn next_ms2(&mut self) -> Option<Result<Ms2Item, TimsTofParseError>> {
        loop {
            if let Some(item) = self.ready.pop_front() {
                return Some(item);
            }
            let total = self.ms2.precursors.len();
            if self.next >= total {
                return None;
            }
            if self.calibration.is_none() {
                match self.ms1_frames() {
                    Ok(ms1) => self.calibration = Some(ms1.calibration.clone()),
                    Err(e) => {
                        self.next = total;
                        return Some(Err(e));
                    }
                }
            }
            let end = (self.next + MS2_BATCH).min(total);
            let calibration = self.calibration.as_ref().expect("set above");
            let built = self.ms2.build(self.next..end, calibration);
            self.ready.extend(built);
            self.next = end;
        }
    }

    /// Record the context of an emitted MS2 when context capture is on.
    fn record_context(&mut self, spec: &Spectrum, parent: Option<usize>, im: f64) {
        if self.capture_scan_context {
            self.scan_contexts.push(ScanContext {
                id: spec.title.clone(),
                parent_id: parent.map(ms1_native_id),
                scan_window: None,
                faims_cv: None,
                noise: None,
                ion_mobility: (im.is_finite() && im > 0.0).then_some(im as f32),
            });
        }
    }

    /// Stream MS2 in `chunk_size` batches, each paired with a bounded
    /// [`Ms1Link`] (the MS1 frames reached in that chunk plus, per MS2, the
    /// index of its parent MS1 frame). Mirrors `MzMLReader::read_with_ms1_chunked`
    /// so the MS1-linked search paths consume a `.d` like mzML. The last MS1
    /// of a chunk is carried into the next chunk (index 0) so the first MS2 of
    /// a chunk still links to its parent across the boundary. Stops after
    /// `cap` MS2.
    ///
    /// Returns `(error_count, sample_errors)`: unreadable MS2 are counted and
    /// skipped. The quantification captures work on this path too.
    pub fn read_with_ms1_chunked<F>(
        &mut self,
        chunk_size: usize,
        cap: usize,
        mut on_chunk: F,
    ) -> (usize, Vec<String>)
    where
        F: FnMut(Vec<Spectrum>, Ms1Link),
    {
        let mut err_count = 0usize;
        let mut first_errors: Vec<String> = Vec::new();
        if let Err(e) = self.ms1_frames() {
            return (1, vec![e.to_string()]);
        }
        let mut chunk: Vec<Spectrum> = Vec::new();
        let mut ms1_peaks: Vec<Vec<(f64, f32)>> = Vec::new();
        let mut ms2_to_ms1: Vec<Option<usize>> = Vec::new();
        // (frame id, index into `ms1_peaks`) of the latest MS1 reached.
        let mut current: Option<(usize, usize)> = None;
        let mut emitted = 0usize;
        while emitted < cap {
            let (spec, parent, im) = match self.next_ms2() {
                None => break,
                Some(Ok(x)) => x,
                Some(Err(e)) => {
                    err_count += 1;
                    if first_errors.len() < 3 {
                        first_errors.push(e.to_string());
                    }
                    continue;
                }
            };
            if let Some(p) = parent {
                let keep = self.capture_run_ms1;
                let run_ms1 = &mut self.run_ms1;
                let ms1 = self.ms1.as_mut().expect("opened above");
                ms1.take_up_to(Some(p), |id, scan| {
                    let peaks = if keep {
                        let peaks = scan.peaks.clone();
                        run_ms1.push(scan);
                        peaks
                    } else {
                        scan.peaks
                    };
                    ms1_peaks.push(peaks);
                    current = Some((id, ms1_peaks.len() - 1));
                });
            }
            self.record_context(&spec, parent, im);
            let link = match (parent, current) {
                (Some(p), Some((id, idx))) if id == p => Some(idx),
                _ => None,
            };
            chunk.push(spec);
            ms2_to_ms1.push(link);
            emitted += 1;
            if chunk.len() >= chunk_size {
                let carry = current.map(|(id, idx)| (id, ms1_peaks[idx].clone()));
                on_chunk(
                    std::mem::take(&mut chunk),
                    Ms1Link {
                        ms1_peaks: std::mem::take(&mut ms1_peaks),
                        ms2_to_ms1: std::mem::take(&mut ms2_to_ms1),
                    },
                );
                current = carry.map(|(id, peaks)| {
                    ms1_peaks.push(peaks);
                    (id, 0)
                });
            }
        }
        if !chunk.is_empty() {
            on_chunk(
                chunk,
                Ms1Link {
                    ms1_peaks,
                    ms2_to_ms1,
                },
            );
        }
        if self.capture_run_ms1 {
            self.flush_ms1();
        }
        (err_count, first_errors)
    }

    /// Capture the MS1 frames after the last MS2.
    fn flush_ms1(&mut self) {
        let run_ms1 = &mut self.run_ms1;
        if let Some(ms1) = self.ms1.as_mut() {
            ms1.take_up_to(None, |_, scan| run_ms1.push(scan));
        }
    }

    /// Collision energy (eV) and isolation width of the first `n` MS2 spectra.
    /// The `.d` records no FAIMS voltage.
    pub fn sample_ms2_meta(&self, n: usize) -> crate::Ms2MetaSample {
        let mut out = crate::Ms2MetaSample::default();
        for p in self.ms2.precursors.iter().take(n) {
            out.n_ms2 += 1;
            if p.collision_energy.is_finite() && p.collision_energy > 0.0 {
                out.collision_energies.push(p.collision_energy);
            }
            if p.isolation_width.is_finite() && p.isolation_width > 0.0 {
                out.isolation_widths.push(p.isolation_width);
            }
        }
        out
    }
}

impl Iterator for TimsTofReader {
    type Item = Result<Spectrum, TimsTofParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_ms2() {
            Some(Ok((spec, parent, im))) => {
                if self.capture_run_ms1 {
                    if let Err(e) = self.ms1_frames() {
                        self.capture_run_ms1 = false;
                        return Some(Err(e));
                    }
                    if let Some(p) = parent {
                        let run_ms1 = &mut self.run_ms1;
                        let ms1 = self.ms1.as_mut().expect("opened above");
                        ms1.take_up_to(Some(p), |_, scan| run_ms1.push(scan));
                    }
                }
                self.record_context(&spec, parent, im);
                Some(Ok(spec))
            }
            Some(Err(e)) => Some(Err(e)),
            None => {
                if self.capture_run_ms1 {
                    self.flush_ms1();
                }
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::timstof_ms2::Selection;

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

    fn precursor(id: usize, mz: f64, charge: Option<i64>, width: f64) -> PasefPrecursor {
        PasefPrecursor {
            id,
            mz,
            charge,
            intensity: Some(5000.0),
            parent: Some(7),
            scan_number: Some(700.0),
            isolation_width: width,
            collision_energy: 28.0,
            selections: vec![Selection {
                frame_id: 8,
                scan_begin: 690,
                scan_end: 715,
            }],
        }
    }

    #[test]
    fn to_spectrum_maps_precursor_and_peaks() {
        let p = precursor(42, 654.32, Some(2), 2.0);
        let peaks = vec![(100.0, 10.0_f32), (200.5, 20.0_f32)];
        let spec = to_spectrum(&p, Some(1800.0), peaks.clone()).expect("convertible");
        // The precursor id is the scan number.
        assert_eq!(spec.scan, Some(42));
        assert_eq!(spec.title, "scan=42");
        assert_eq!(spec.precursor_mz, 654.32);
        assert_eq!(spec.precursor_charge, Some(2));
        assert_eq!(spec.precursor_intensity, Some(5000.0));
        // RT is carried in seconds.
        assert_eq!(spec.rt_seconds, Some(1800.0));
        // Symmetric offsets are half the total isolation width.
        assert_eq!(spec.isolation_lower_offset, Some(1.0));
        assert_eq!(spec.isolation_upper_offset, Some(1.0));
        assert_eq!(spec.peaks, peaks);
    }

    #[test]
    fn to_spectrum_skips_a_nonpositive_precursor_mz() {
        assert!(to_spectrum(&precursor(1, 0.0, Some(2), 2.0), None, Vec::new()).is_none());
        assert!(to_spectrum(&precursor(1, f64::NAN, Some(2), 2.0), None, Vec::new()).is_none());
    }

    #[test]
    fn to_spectrum_treats_zero_or_missing_charge_as_unknown() {
        let spec = to_spectrum(&precursor(1, 500.0, Some(0), 0.0), None, Vec::new()).unwrap();
        assert_eq!(spec.precursor_charge, None);
        // A non-positive isolation width yields no offsets.
        assert_eq!(spec.isolation_lower_offset, None);
        assert_eq!(spec.isolation_upper_offset, None);
        let spec = to_spectrum(&precursor(1, 500.0, None, 2.0), None, Vec::new()).unwrap();
        assert_eq!(spec.precursor_charge, None);
    }
}
