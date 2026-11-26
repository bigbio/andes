pub mod archive;
pub(crate) mod ms1;
pub(crate) mod ms2;

use std::collections::{HashMap, VecDeque};
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use timsrust::converters::{ConvertableDomain, Scan2ImConverter};
use timsrust::readers::{FrameReader, FrameReaderError, MetadataReader};
use timsrust::MSLevel;

use model::scan::{Ms1Scan, ScanContext};

use crate::timstof::ms1::{
    centroid_frame, centroid_profile, fit_calibration, match_precursors, offset_curve,
    sum_over_mobility, CalibrationFit, FrameConditions, MobilityCentroid, MzCalibration,
    OffsetCurve, MIN_REFINE_WINDOW_PPM, NARROW_WINDOW_PPM, REFINE_ROUNDS, REFINE_SIGMAS,
    WIDE_WINDOW_PPM,
};
use crate::timstof::ms2::{
    has_table, ms2_peaks, push_scans, read_dia_windows, read_frame_rows, read_pasef_precursors,
    slice_dia_windows, FrameRow, PasefPrecursor,
};
use crate::{Ms1Link, Spectrum};

#[derive(Debug, thiserror::Error)]
pub enum TimsTofParseError {
    #[error("failed to open Bruker .d '{path}': {source}{hint}")]
    Open {
        path: String,
        source: Box<FrameReaderError>,
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

fn open_error_hint(message: &str) -> &'static str {
    if message.contains("Compression type 1") {
        " (this .d uses the frame compression of early timsTOF acquisition software, \
         which the native reader does not decode; convert the run to mzML and search that)"
    } else {
        ""
    }
}

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

const MS1_BATCH: usize = 64;
const CALIBRATION_FRAMES: usize = 400;

struct Ms1Frames {
    frames: FrameReader,
    positions: Vec<usize>,
    ids: Vec<usize>,
    temperatures: HashMap<usize, [f64; 2]>,
    calibration: MzCalibration,
    im: timsrust::converters::Scan2ImConverter,
    resolve: bool,
    cursor: usize,
    ready: VecDeque<DecodedMs1>,
}

struct DecodedMs1 {
    id: usize,
    scan: Ms1Scan,
    summed: Option<Vec<(f64, f32)>>,
}

struct TdfTables {
    temperatures: HashMap<usize, [f64; 2]>,
    precursor_peaks: HashMap<usize, Vec<f64>>,
}

fn read_tables(dir: &Path) -> Result<TdfTables, rusqlite::Error> {
    let conn = rusqlite::Connection::open_with_flags(
        dir.join("analysis.tdf"),
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )?;
    let mut temperatures = HashMap::new();
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
    if !has_table(&conn, "Precursors")? {
        return Ok(TdfTables {
            temperatures,
            precursor_peaks,
        });
    }
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

fn frame_centroids(frame: &timsrust::Frame) -> Vec<crate::timstof::ms1::TofCentroid> {
    centroid_profile(&sum_over_mobility(&frame.tof_indices, &frame.intensities))
}

fn to_peaks(
    centroids: &[crate::timstof::ms1::TofCentroid],
    calibration: &MzCalibration,
    frame: FrameConditions,
) -> Vec<(f64, f32)> {
    finite_sorted_peaks(centroids.iter().map(|c| {
        let cond = crate::timstof::ms1::conditions(frame, c.intensity);
        (calibration.mz(c.tof, cond), c.intensity as f32)
    }))
}

/// The peaks with a finite positive m/z and intensity, in m/z order.
pub(crate) fn finite_sorted_peaks(peaks: impl Iterator<Item = (f64, f32)>) -> Vec<(f64, f32)> {
    let mut peaks: Vec<(f64, f32)> = peaks
        .filter(|&(mz, i)| mz.is_finite() && mz > 0.0 && i.is_finite() && i > 0.0)
        .collect();
    if peaks.windows(2).any(|w| w[0].0 > w[1].0) {
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    peaks.shrink_to_fit();
    peaks
}

fn to_resolved_peaks(
    centroids: &[MobilityCentroid],
    calibration: &MzCalibration,
    im: &timsrust::converters::Scan2ImConverter,
    frame: FrameConditions,
) -> (Vec<(f64, f32)>, Vec<f32>) {
    let mut rows: Vec<(f64, f32, f32)> = centroids
        .iter()
        .map(|c| {
            let cond = crate::timstof::ms1::conditions(frame, c.intensity);
            (
                calibration.mz(c.tof, cond),
                c.intensity as f32,
                im.convert(c.scan) as f32,
            )
        })
        .filter(|&(mz, i, k0)| {
            mz.is_finite() && mz > 0.0 && i.is_finite() && i > 0.0 && k0.is_finite()
        })
        .collect();
    if rows.windows(2).any(|w| w[0].0 > w[1].0) {
        rows.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    (
        rows.iter().map(|r| (r.0, r.1)).collect(),
        rows.iter().map(|r| r.2).collect(),
    )
}

impl Ms1Frames {
    fn open(dir: &Path, resolve: bool) -> Result<Self, TimsTofParseError> {
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
            im: metadata.im_converter,
            resolve,
            cursor: 0,
            ready: VecDeque::new(),
        };
        match ms1.calibrate(&base, &tables.precursor_peaks) {
            Some(fit) => {
                log::info!(
                    "MS1 m/z calibration from {} precursor peaks: residual to the \
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
            None => log::warn!(
                "too few precursor peaks matched their MS1 frame for an m/z \
                 calibration; MS1 m/z uses the acquisition-range model"
            ),
        }
        Ok(ms1)
    }

    fn temperature(&self, id: usize) -> [f64; 2] {
        self.temperatures.get(&id).copied().unwrap_or([0.0; 2])
    }

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
        type Decoded = (usize, FrameConditions, Vec<crate::timstof::ms1::TofCentroid>);
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
        let pass = |calibration: &MzCalibration, offset: &OffsetCurve, window: f64| {
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
        let wide = pass(base, &OffsetCurve::constant(0.0), WIDE_WINDOW_PPM);
        let offset = offset_curve(&wide, base)?;
        let mut fit = fit_calibration(&pass(base, &offset, NARROW_WINDOW_PPM), base)?;
        for _ in 0..REFINE_ROUNDS {
            let window = (REFINE_SIGMAS * fit.sigma_ppm).max(MIN_REFINE_WINDOW_PPM);
            match fit_calibration(
                &pass(&fit.calibration, &OffsetCurve::constant(0.0), window),
                base,
            ) {
                Some(next) => fit = next,
                None => break,
            }
        }
        Some(fit)
    }

    fn decode_batch(&mut self) {
        let end = (self.cursor + MS1_BATCH).min(self.positions.len());
        let batch: Vec<(usize, usize)> = (self.cursor..end)
            .map(|k| (self.positions[k], self.ids[k]))
            .collect();
        self.cursor = end;
        let decoded: Vec<DecodedMs1> = batch
            .par_iter()
            .map(|&(pos, id)| {
                let (rt, peaks, mobility, summed) = match self.frames.get(pos) {
                    Ok(frame) => {
                        let cond = FrameConditions {
                            t: self.temperature(id),
                            rt_s: frame.rt_in_seconds,
                        };
                        let (summed, resolved) = centroid_frame(
                            &frame.scan_offsets,
                            &frame.tof_indices,
                            &frame.intensities,
                            self.resolve,
                        );
                        let summed = to_peaks(&summed, &self.calibration, cond);
                        if self.resolve {
                            let (peaks, mobility) =
                                to_resolved_peaks(&resolved, &self.calibration, &self.im, cond);
                            (frame.rt_in_seconds, peaks, mobility, Some(summed))
                        } else {
                            (frame.rt_in_seconds, summed, Vec::new(), None)
                        }
                    }
                    Err(_) => (
                        self.frames
                            .get_frame_without_coordinates(pos)
                            .map(|f| f.rt_in_seconds)
                            .unwrap_or(f64::NAN),
                        Vec::new(),
                        Vec::new(),
                        self.resolve.then(Vec::new),
                    ),
                };
                DecodedMs1 {
                    id,
                    scan: Ms1Scan {
                        rt,
                        peaks,
                        id: ms1_native_id(id),
                        faims_cv: None,
                        mobility,
                    },
                    summed,
                }
            })
            .collect();
        self.ready.extend(decoded);
    }

    fn take_up_to(&mut self, frame_id: Option<usize>, mut sink: impl FnMut(DecodedMs1)) {
        loop {
            while let Some(d) = self.ready.front() {
                if frame_id.is_some_and(|limit| d.id > limit) {
                    return;
                }
                sink(self.ready.pop_front().expect("front exists"));
            }
            if self.cursor >= self.positions.len() {
                return;
            }
            self.decode_batch();
        }
    }
}

type Ms2Item = (Spectrum, Option<usize>, f64);

fn ms1_native_id(frame_id: usize) -> String {
    format!("frame={frame_id}")
}

const MS2_BATCH: usize = 1024;

struct Ms2Source {
    frames: FrameReader,
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
        let mut precursors = read_pasef_precursors(&conn).map_err(tables_err)?;
        if precursors.is_empty() {
            precursors = read_dia_windows(&conn).map_err(tables_err)?;
        }
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
            .filter(|_| !p.dia)
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
    let precursor_charge = p.charge.and_then(model::spectrum::known_precursor_charge);
    let precursor_intensity = p
        .intensity
        .filter(|i| i.is_finite() && *i > 0.0)
        .map(|i| i as f32);
    let (isolation_lower_offset, isolation_upper_offset) = if p.isolation_width > 0.0 {
        let half = p.isolation_width / 2.0;
        (Some(half), Some(half))
    } else {
        (None, None)
    };
    Some(Spectrum {
        title: format!("scan={scan_number}"),
        precursor_mz: p.mz,
        precursor_intensity,
        precursor_charge,
        rt_seconds,
        scan: Some(scan_number),
        peaks,
        activation_method: None,
        isolation_lower_offset,
        isolation_upper_offset,
        isolation_target_mz: p.dia.then_some(p.mz),
        mobility_window: p
            .selections
            .first()
            .filter(|_| p.dia)
            .map(|s| (s.scan_begin as u32, s.scan_end as u32)),
    })
}

pub struct TimsTofReader {
    ms2: Ms2Source,
    dir: PathBuf,
    next: usize,
    ready: VecDeque<Result<Ms2Item, TimsTofParseError>>,
    calibration: Option<MzCalibration>,
    capture_run_ms1: bool,
    run_ms1: Vec<Ms1Scan>,
    capture_scan_context: bool,
    scan_contexts: Vec<ScanContext>,
    ms1: Option<Ms1Frames>,
    ms1_mobility: bool,
}

impl TimsTofReader {
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
            ms1_mobility: false,
        })
    }

    pub fn with_run_ms1_capture(mut self, capture: bool) -> Self {
        self.capture_run_ms1 = capture;
        self
    }

    pub fn with_ms1_mobility(mut self, resolve: bool) -> Self {
        self.ms1_mobility = resolve;
        self
    }

    pub fn with_dia_mobility_slices(mut self, n: usize) -> Self {
        let p = std::mem::take(&mut self.ms2.precursors);
        self.ms2.precursors = slice_dia_windows(p, n);
        self
    }

    pub fn with_scan_context_capture(mut self, capture: bool) -> Self {
        self.capture_scan_context = capture;
        self
    }

    /// 1/K0 of every TIMS scan number up to the last scan a DIA window selects; empty without
    /// DIA windows. Index `s` is scan `s` of a frame.
    pub fn dia_scan_mobility(&self) -> Vec<f32> {
        let last = self
            .ms2
            .precursors
            .iter()
            .filter(|p| p.dia)
            .flat_map(|p| p.selections.iter().map(|s| s.scan_end))
            .max();
        last.map_or_else(Vec::new, |n| {
            (0..=n).map(|s| self.ms2.im.convert(s as f64) as f32).collect()
        })
    }

    pub fn take_run_ms1(&mut self) -> Vec<Ms1Scan> {
        std::mem::take(&mut self.run_ms1)
    }

    pub fn take_scan_contexts(&mut self) -> Vec<ScanContext> {
        std::mem::take(&mut self.scan_contexts)
    }

    pub fn len(&self) -> usize {
        self.ms2.precursors.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ms2.precursors.is_empty()
    }

    fn ms1_frames(&mut self) -> Result<&mut Ms1Frames, TimsTofParseError> {
        if self.ms1.is_none() {
            self.ms1 = Some(Ms1Frames::open(
                &self.dir,
                self.ms1_mobility && self.capture_run_ms1,
            )?);
        }
        Ok(self.ms1.as_mut().expect("just opened"))
    }

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
        let mut current: Option<(usize, usize)> = None;
        let mut current_linked = false;
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
                ms1.take_up_to(Some(p), |d| {
                    let peaks = match (keep, d.summed) {
                        (true, Some(summed)) => {
                            run_ms1.push(d.scan);
                            summed
                        }
                        (true, None) => {
                            let peaks = d.scan.peaks.clone();
                            run_ms1.push(d.scan);
                            peaks
                        }
                        (false, summed) => summed.unwrap_or(d.scan.peaks),
                    };
                    match current {
                        Some((_, idx)) if !current_linked => {
                            ms1_peaks[idx] = peaks;
                            current = Some((d.id, idx));
                        }
                        _ => {
                            ms1_peaks.push(peaks);
                            current = Some((d.id, ms1_peaks.len() - 1));
                        }
                    }
                    current_linked = false;
                });
            }
            self.record_context(&spec, parent, im);
            let link = match (parent, current) {
                (Some(p), Some((id, idx))) if id == p => Some(idx),
                _ => None,
            };
            current_linked |= link.is_some();
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
                current_linked = false;
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

    fn flush_ms1(&mut self) {
        let run_ms1 = &mut self.run_ms1;
        if let Some(ms1) = self.ms1.as_mut() {
            ms1.take_up_to(None, |d| run_ms1.push(d.scan));
        }
    }

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
                        ms1.take_up_to(Some(p), |d| run_ms1.push(d.scan));
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
    use crate::timstof::ms2::Selection;

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
            dia: false,
        }
    }

    #[test]
    fn to_spectrum_maps_precursor_and_peaks() {
        let p = precursor(42, 654.32, Some(2), 2.0);
        let peaks = vec![(100.0, 10.0_f32), (200.5, 20.0_f32)];
        let spec = to_spectrum(&p, Some(1800.0), peaks.clone()).expect("convertible");
        assert_eq!(spec.scan, Some(42));
        assert_eq!(spec.title, "scan=42");
        assert_eq!(spec.precursor_mz, 654.32);
        assert_eq!(spec.precursor_charge, Some(2));
        assert_eq!(spec.precursor_intensity, Some(5000.0));
        assert_eq!(spec.rt_seconds, Some(1800.0));
        assert_eq!(spec.isolation_lower_offset, Some(1.0));
        assert_eq!(spec.isolation_upper_offset, Some(1.0));
        assert_eq!(spec.peaks, peaks);
        assert_eq!(spec.isolation_target_mz, None);
        assert_eq!(spec.mobility_window, None);
        let dia = PasefPrecursor { dia: true, ..p };
        let spec = to_spectrum(&dia, Some(1800.0), Vec::new()).unwrap();
        assert_eq!(spec.isolation_target_mz, Some(654.32));
        assert_eq!(spec.mobility_window, Some((690, 715)));
    }

    #[test]
    fn read_tables_needs_no_precursors_table_but_reports_a_broken_one() {
        let dir = tempfile::tempdir().unwrap();
        let conn = rusqlite::Connection::open(dir.path().join("analysis.tdf")).unwrap();
        conn.execute_batch(
            "CREATE TABLE Frames (Id INTEGER, T1 REAL, T2 REAL);
             INSERT INTO Frames VALUES (1, 25.0, 26.0);",
        )
        .unwrap();
        let t = read_tables(dir.path()).unwrap();
        assert_eq!(t.temperatures[&1], [25.0, 26.0]);
        assert!(t.precursor_peaks.is_empty());
        conn.execute_batch("CREATE TABLE Precursors (Parent INTEGER);")
            .unwrap();
        assert!(read_tables(dir.path()).is_err());
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
        assert_eq!(spec.isolation_lower_offset, None);
        assert_eq!(spec.isolation_upper_offset, None);
        let spec = to_spectrum(&precursor(1, 500.0, None, 2.0), None, Vec::new()).unwrap();
        assert_eq!(spec.precursor_charge, None);
    }
}
