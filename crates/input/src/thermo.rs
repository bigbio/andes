use std::path::Path;

use std::collections::HashMap;

use thermorawfilereader::schema::{DissociationMethod, MassAnalyzer};
use thermorawfilereader::{RawFileReader, RawSpectrum};

use model::scan::{Ms1Scan, NoiseProfile, ProductScan, ScanContext};
use model::{ActivationMethod, InstrumentType};

use crate::{Ms1Link, Spectrum};

#[derive(Debug, thiserror::Error)]
pub enum ThermoParseError {
    #[error("failed to open Thermo .raw '{path}': {source}{}", runtime_hint(.source))]
    Open {
        path: String,
        source: std::io::Error,
    },
}

fn runtime_hint(source: &std::io::Error) -> &'static str {
    match source.kind() {
        std::io::ErrorKind::NotFound | std::io::ErrorKind::InvalidData => "",
        _ => {
            " (opening a .raw needs the .NET 8 runtime: install it, or set DOTNET_ROOT \
              to its directory)"
        }
    }
}

pub struct ThermoRawReader {
    handle: RawFileReader,
    path: std::path::PathBuf,
    /// Further readers of the same file: a batch of scans is read and converted by `handle` and
    /// these, each over its own contiguous part of the batch.
    helpers: Vec<RawFileReader>,
    /// Scans read and converted ahead, in scan order, not yet taken.
    fetched: std::collections::VecDeque<FetchedScan>,
    next: usize,
    len: usize,
    ms_level_filter: Option<u8>,
    capture_run_ms1: bool,
    run_ms1: Vec<Ms1Scan>,
    capture_product_scans: bool,
    product_scans: Vec<ProductScan>,
    capture_scan_context: bool,
    scan_contexts: Vec<ScanContext>,
    capture_noise: bool,
}

/// What a read converts from each scan besides the MS2 spectra.
#[derive(Clone, Copy)]
struct Capture {
    /// The MS1 scans, for the run MS1 or (`linked`) the MS1 link of the MS2 scans.
    ms1: bool,
    product: bool,
    context: bool,
    noise: bool,
    level: Option<u8>,
}

/// One scan, read and converted: what it gives the reader's outputs.
struct FetchedScan {
    level: u8,
    faims_cv: Option<f32>,
    ms1: Option<Ms1Scan>,
    product: Option<ProductScan>,
    ms2: Option<(Spectrum, Option<ScanContext>)>,
}

/// Scans each reader of a parallel read converts per batch.
const FETCH_SCANS_PER_READER: usize = 2048;

fn fetch_scan(handle: &RawFileReader, index: usize, cap: Capture) -> Option<FetchedScan> {
    let raw = handle.get(index)?;
    let level = raw.ms_level();
    let faims = faims_cv(&raw);
    let ms1 = (level == 1 && cap.ms1).then(|| Ms1Scan {
        rt: raw.time() * 60.0,
        peaks: extract_peaks(&raw),
        id: raw.native_id(),
        faims_cv: faims,
        mobility: Vec::new(),
    });
    if level >= 3 && cap.product {
        return Some(FetchedScan {
            level,
            faims_cv: faims,
            ms1,
            product: Some(product_scan(handle, index, &raw, cap.noise)),
            ms2: None,
        });
    }
    let ms2 = cap
        .level
        .is_none_or(|want| want == level)
        .then(|| ThermoRawReader::convert(&raw, handle))
        .flatten()
        .map(|spec| {
            let context = cap
                .context
                .then(|| scan_context(handle, index, &raw, &spec, cap.noise));
            (spec, context)
        });
    Some(FetchedScan {
        level,
        faims_cv: faims,
        ms1,
        product: None,
        ms2,
    })
}

fn scan_context(
    handle: &RawFileReader,
    index: usize,
    raw: &RawSpectrum,
    spec: &Spectrum,
    noise: bool,
) -> ScanContext {
    ScanContext {
        id: spec.title.clone(),
        parent_id: master_scan_number(handle, index)
            .map(|scan| with_scan_number(&spec.title, scan)),
        scan_window: scan_window(raw),
        faims_cv: faims_cv(raw),
        noise: noise.then(|| reporter_noise(handle, index)).flatten(),
        ion_mobility: None,
    }
}

/// The MS3+ scan at `index`, linked to its parent scan by the master scan number.
fn product_scan(
    handle: &RawFileReader,
    index: usize,
    raw: &RawSpectrum,
    noise: bool,
) -> ProductScan {
    let parent_id =
        master_scan_number(handle, index).map(|scan| with_scan_number(&raw.native_id(), scan));
    ProductScan {
        id: raw.native_id(),
        scan: Some(raw.index() as i32 + 1),
        rt_seconds: Some(raw.time() * 60.0),
        ms_level: raw.ms_level(),
        parent_id,
        precursor_mz: raw
            .filter_string()
            .and_then(last_reaction_mz)
            .or_else(|| raw.precursor().map(|p| p.mz()))
            .unwrap_or(0.0),
        scan_window: scan_window(raw),
        peaks: extract_peaks(raw),
        noise: noise.then(|| reporter_noise(handle, index)).flatten(),
    }
}

impl ThermoRawReader {
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
            path: p.to_path_buf(),
            helpers: Vec::new(),
            fetched: std::collections::VecDeque::new(),
            next: 0,
            len,
            ms_level_filter: Some(2),
            capture_run_ms1: false,
            run_ms1: Vec::new(),
            capture_product_scans: false,
            product_scans: Vec::new(),
            capture_scan_context: false,
            scan_contexts: Vec::new(),
            capture_noise: false,
        })
    }

    /// Reads and converts the scans with `readers` readers of the file at once (one per
    /// thread), each over its own part of a batch; the scans come out in scan order as with
    /// one reader. A further reader that fails to open leaves the read with fewer readers.
    pub fn with_parallel_reads(mut self, readers: usize) -> Self {
        for _ in 1..readers {
            match RawFileReader::open(&self.path) {
                Ok(mut h) => {
                    h.set_centroid_spectra(true);
                    self.helpers.push(h);
                }
                Err(e) => {
                    log::debug!(
                        "{}: a further reader did not open ({e}); reading with {}",
                        self.path.display(),
                        1 + self.helpers.len()
                    );
                    break;
                }
            }
        }
        self
    }

    pub fn with_run_ms1_capture(mut self, capture: bool) -> Self {
        self.capture_run_ms1 = capture;
        self
    }

    pub fn with_product_scan_capture(mut self, capture: bool) -> Self {
        self.capture_product_scans = capture;
        self
    }

    pub fn with_scan_context_capture(mut self, capture: bool) -> Self {
        self.capture_scan_context = capture;
        self
    }

    pub fn with_noise_capture(mut self, capture: bool) -> Self {
        self.capture_noise = capture;
        self
    }

    pub fn take_scan_contexts(&mut self) -> Vec<ScanContext> {
        std::mem::take(&mut self.scan_contexts)
    }

    pub fn take_run_ms1(&mut self) -> Vec<Ms1Scan> {
        std::mem::take(&mut self.run_ms1)
    }

    pub fn take_product_scans(&mut self) -> Vec<ProductScan> {
        std::mem::take(&mut self.product_scans)
    }

    pub fn with_all_ms_levels(mut self) -> Self {
        self.ms_level_filter = None;
        self
    }

    pub fn with_ms_level(mut self, level: Option<u8>) -> Self {
        self.ms_level_filter = level;
        self
    }

    fn capture(&self, linked: bool) -> Capture {
        Capture {
            ms1: linked || self.capture_run_ms1,
            product: self.capture_product_scans,
            context: self.capture_scan_context,
            noise: self.capture_noise,
            level: if linked { Some(2) } else { self.ms_level_filter },
        }
    }

    /// Reads and converts the next batch of scans into `fetched`: one scan with one reader,
    /// else a batch split over the readers, each converting a contiguous part on its own thread.
    fn fetch_batch(&mut self, cap: Capture) {
        let lo = self.next;
        if self.helpers.is_empty() {
            self.next += 1;
            self.fetched.extend(fetch_scan(&self.handle, lo, cap));
            return;
        }
        let readers: Vec<&RawFileReader> =
            std::iter::once(&self.handle).chain(&self.helpers).collect();
        let hi = (lo + FETCH_SCANS_PER_READER * readers.len()).min(self.len);
        self.next = hi;
        let per = (hi - lo).div_ceil(readers.len());
        let parts: Vec<Vec<FetchedScan>> = std::thread::scope(|s| {
            let workers: Vec<_> = readers
                .iter()
                .enumerate()
                .map(|(j, &reader)| {
                    let (a, b) = ((lo + j * per).min(hi), (lo + (j + 1) * per).min(hi));
                    s.spawn(move || {
                        (a..b)
                            .filter_map(|i| fetch_scan(reader, i, cap))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            workers
                .into_iter()
                .map(|w| w.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect()
        });
        for part in parts {
            self.fetched.extend(part);
        }
    }

    /// The next scan read and converted, in scan order; `None` after the last scan.
    fn next_fetched(&mut self, cap: Capture) -> Option<FetchedScan> {
        loop {
            if let Some(f) = self.fetched.pop_front() {
                return Some(f);
            }
            if self.next >= self.len {
                return None;
            }
            self.fetch_batch(cap);
        }
    }

    pub fn sample_level(&mut self, level: u8, n: usize) -> Vec<Spectrum> {
        let stride = (self.len / n.max(1) / 4).max(1);
        self.handle.set_signal_loading(false);
        let kept: Vec<usize> = (0..self.len)
            .step_by(stride)
            .filter(|&i| {
                self.handle
                    .get(i)
                    .is_some_and(|raw| raw.ms_level() == level)
            })
            .take(n)
            .collect();
        self.handle.set_signal_loading(true);
        kept.into_iter()
            .filter_map(|i| self.handle.get(i))
            .map(|raw| Spectrum {
                title: raw.native_id(),
                precursor_mz: 0.0,
                precursor_intensity: None,
                precursor_charge: None,
                rt_seconds: Some(raw.time() * 60.0),
                scan: Some(raw.index() as i32 + 1),
                peaks: extract_peaks(&raw),
                activation_method: None,
                isolation_lower_offset: None,
                isolation_upper_offset: None,
                isolation_target_mz: None,
                mobility_window: None,
            })
            .collect()
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn convert(raw: &RawSpectrum, handle: &RawFileReader) -> Option<Spectrum> {
        let scan_number = raw.index() as i32 + 1;

        let p = raw.precursor()?;
        let precursor_mz = p.mz();
        if precursor_mz.is_nan() || precursor_mz <= 0.0 {
            return None;
        }
        let precursor_charge = model::spectrum::known_precursor_charge(i64::from(p.charge()));
        let inten = p.intensity();
        let precursor_intensity = if inten > 0.0 { Some(inten) } else { None };

        let w = p.isolation_window();
        let (lo, tg) = (w.lower(), w.target());
        let (lo, up) = full_width_isolation_bounds(
            lo,
            w.upper(),
            isolation_width_trailer(handle, raw.index(), raw.ms_level()),
        );
        let (isolation_lower_offset, isolation_upper_offset, isolation_target_mz) =
            if tg > 0.0 && up > lo {
                (Some((tg - lo).abs()), Some((up - tg).abs()), Some(tg))
            } else {
                (None, None, None)
            };

        Some(Spectrum {
            title: raw.native_id(),
            precursor_mz,
            precursor_intensity,
            precursor_charge,
            rt_seconds: Some(raw.time() * 60.0),
            scan: Some(scan_number),
            peaks: extract_peaks(raw),
            activation_method: scan_activation(raw),
            isolation_lower_offset,
            isolation_upper_offset,
            isolation_target_mz,
            mobility_window: None,
        })
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
        let mut chunk: Vec<Spectrum> = Vec::new();
        let mut ms1_peaks: Vec<Vec<(f64, f32)>> = Vec::new();
        let mut ms2_to_ms1: Vec<Option<usize>> = Vec::new();
        let mut tracker = crate::ms1::Ms1Tracker::default();
        let mut emitted_ms2 = 0usize;
        let mut skipped_higher_ms = 0usize;
        let capture = self.capture(true);

        loop {
            let capped = cap > 0 && emitted_ms2 >= cap;
            if capped && !self.capture_product_scans {
                break;
            }
            let Some(scan) = self.next_fetched(capture) else {
                break;
            };
            if capped && scan.level < 3 {
                break;
            }
            match scan.level {
                1 => {
                    let Some(ms1) = scan.ms1 else { continue };
                    let idx = ms1_peaks.len();
                    if self.capture_run_ms1 {
                        ms1_peaks.push(ms1.peaks.clone());
                        self.run_ms1.push(ms1);
                    } else {
                        ms1_peaks.push(ms1.peaks);
                    }
                    tracker.push(idx, scan.faims_cv);
                }
                level if level >= 3 => match scan.product {
                    Some(p) => self.product_scans.push(p),
                    None => skipped_higher_ms += 1,
                },
                2 => {
                    let Some((spec, context)) = scan.ms2 else {
                        continue;
                    };
                    self.scan_contexts.extend(context);
                    chunk.push(spec);
                    ms2_to_ms1.push(tracker.link(scan.faims_cv));
                    emitted_ms2 += 1;
                    if chunk.len() >= chunk_size {
                        // Carry the MS1 scans later MS2 scans can still link to.
                        let carried: Vec<Vec<(f64, f32)>> = tracker
                            .live()
                            .iter()
                            .map(|&j| ms1_peaks[j].clone())
                            .collect();
                        let link = Ms1Link {
                            ms1_peaks: std::mem::take(&mut ms1_peaks),
                            ms2_to_ms1: std::mem::take(&mut ms2_to_ms1),
                        };
                        on_chunk(std::mem::take(&mut chunk), link);
                        ms1_peaks = carried;
                        tracker.renumber_to_live();
                    }
                }
                _ => {}
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
        if skipped_higher_ms > 0 {
            log::info!(
                "skipped {skipped_higher_ms} MS3+ scans (not searched; \
                 e.g. TMT SPS-MS3 reporter-quant scans)"
            );
        }
        (0, Vec::new())
    }
}

fn master_scan_number(handle: &RawFileReader, index: usize) -> Option<usize> {
    let trailers = handle.get_raw_trailers_for(index)?;
    let value = trailers
        .iter()
        .find_map(|t| is_master_scan_label(t.label).then(|| t.value.trim().to_string()))?;
    value.parse::<usize>().ok().filter(|&scan| scan > 0)
}

fn isolation_width_trailer(handle: &RawFileReader, index: usize, ms_level: u8) -> Option<f64> {
    if ms_level < 2 {
        return None;
    }
    let label = format!("MS{ms_level} Isolation Width");
    let trailers = handle.get_raw_trailers_for(index)?;
    let value = trailers.iter().find_map(|t| {
        (t.label.trim().trim_end_matches(':').trim() == label).then(|| t.value.trim().to_string())
    })?;
    value
        .parse::<f64>()
        .ok()
        .filter(|w| w.is_finite() && *w > 0.0)
}

fn full_width_isolation_bounds(lower: f64, upper: f64, full_width: Option<f64>) -> (f64, f64) {
    match full_width {
        Some(w) if ((upper - lower) - 2.0 * w).abs() <= 1e-6 => {
            let center = 0.5 * (lower + upper);
            (center - 0.5 * w, center + 0.5 * w)
        }
        _ => (lower, upper),
    }
}

fn scan_window(raw: &RawSpectrum) -> Option<(f64, f64)> {
    let a = raw.acquisition()?;
    let (lo, hi) = (a.low_mz(), a.high_mz());
    (lo.is_finite() && hi > lo).then_some((lo, hi))
}

fn faims_cv(raw: &RawSpectrum) -> Option<f32> {
    raw.acquisition()?
        .compensation_voltages()
        .and_then(|v| v.iter().next())
}

fn is_master_scan_label(label: &str) -> bool {
    label
        .trim()
        .trim_end_matches(':')
        .trim()
        .eq_ignore_ascii_case("Master Scan Number")
}

/// The precursor m/z of the last reaction in a scan filter: `375.5128` in
/// `ITMS + c NSI r d Full ms3 405.5227@cid25.00 375.5128@cid35.00 [99.0000-1137.0000]`, the
/// ion an MS3 scan fragments (the first reaction is its MS2 precursor).
fn last_reaction_mz(filter: &str) -> Option<f64> {
    filter
        .split_whitespace()
        .filter_map(|t| t.split_once('@'))
        .filter_map(|(mz, _)| mz.parse::<f64>().ok())
        .rfind(|mz| mz.is_finite() && *mz > 0.0)
}

fn with_scan_number(native_id: &str, scan: usize) -> String {
    match native_id.rfind("scan=") {
        Some(pos) => {
            let rest = &native_id[pos + "scan=".len()..];
            let digits = rest.chars().take_while(char::is_ascii_digit).count();
            format!("{}scan={scan}{}", &native_id[..pos], &rest[digits..])
        }
        None => format!("controllerType=0 controllerNumber=1 scan={scan}"),
    }
}

fn reporter_noise(handle: &RawFileReader, index: usize) -> Option<NoiseProfile> {
    let data = handle.get_extended_spectrum_data(index, false)?;
    let view = data.view();
    let (mass, noise) = (view.mass()?, view.noise()?);
    let n = mass.len().min(noise.len());
    NoiseProfile::reporter_region((0..n).map(|j| (mass.get(j), noise.get(j))))
}

fn extract_peaks(raw: &RawSpectrum) -> Vec<(f64, f32)> {
    let (mz, intensity) = match raw.data_raw() {
        Some(data) => match (data.mz(), data.intensity()) {
            (Some(mz), Some(intensity)) => (mz, intensity),
            _ => return Vec::new(),
        },
        None => return Vec::new(),
    };

    if mz.len() != intensity.len() {
        static WARNED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
        if !WARNED.swap(true, std::sync::atomic::Ordering::Relaxed) {
            log::warn!(
                "a scan has mismatched peak arrays (m/z {} vs intensity {}); pairing only the \
                 common prefix (reported once; --verbose lists every such scan)",
                mz.len(),
                intensity.len()
            );
        } else {
            log::debug!(
                "scan has mismatched peak arrays (m/z {} vs intensity {}); \
                 pairing only the common prefix",
                mz.len(),
                intensity.len()
            );
        }
    }
    let n = mz.len().min(intensity.len());

    let raw_pairs: Vec<(f64, f32)> = (0..n).map(|i| (mz.get(i), intensity.get(i))).collect();
    sanitize_thermo_peaks(raw_pairs)
}

fn sanitize_thermo_peaks(pairs: Vec<(f64, f32)>) -> Vec<(f64, f32)> {
    let mut peaks: Vec<(f64, f32)> = pairs
        .into_iter()
        .filter(|&(m, it)| m.is_finite() && m > 0.0 && it.is_finite() && it >= 0.0)
        .collect();
    if peaks.windows(2).any(|w| w[0].0 > w[1].0) {
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    peaks
}

fn scan_activation(raw: &RawSpectrum) -> Option<ActivationMethod> {
    let d = raw.precursor()?.activation().dissociation_method();
    let electron = matches!(d.0, 4 | 8 | 16);
    if electron
        && raw
            .filter_string()
            .is_some_and(crate::mzml::filter_has_supplemental_collision)
    {
        return map_dissociation(DissociationMethod::ETHCD);
    }
    map_dissociation(d)
}

fn map_dissociation(d: DissociationMethod) -> Option<ActivationMethod> {
    match d.0 {
        1 => Some(ActivationMethod::CID),
        2 => Some(ActivationMethod::HCD),
        4 | 8 | 16 => Some(ActivationMethod::ETD),
        5 | 6 | 9 | 10 => {
            crate::mzml::log_ethcd_once();
            if crate::mzml::ethcd_as_etd_policy() {
                Some(ActivationMethod::ETD)
            } else {
                Some(ActivationMethod::EThcD)
            }
        }
        _ => None,
    }
}

fn map_analyzer(a: MassAnalyzer) -> Option<InstrumentType> {
    match a.0 {
        1 => Some(InstrumentType::LowRes),
        4 => Some(InstrumentType::TOF),
        5 => Some(InstrumentType::QExactive),
        7 => Some(InstrumentType::OrbitrapAstral),
        _ => None,
    }
}

pub fn detect_activation_instrument<P: AsRef<Path>>(
    path: P,
    peek_n: usize,
) -> (Option<ActivationMethod>, Option<InstrumentType>) {
    let Ok(mut reader) = RawFileReader::open(path.as_ref()) else {
        return (None, None);
    };
    reader.set_signal_loading(false);
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
        if let Some(m) = scan_activation(&raw) {
            *act.entry(m as u8).or_insert(0) += 1;
        }
        if let Some(a) = raw.acquisition() {
            if let Some(t) = map_analyzer(a.mass_analyzer()) {
                *inst.entry(t as u8).or_insert(0) += 1;
            }
        }
    }
    warn_mixed_electron_activation(&act, seen);
    dominant_acquisition(act, inst)
}

fn dominant_acquisition(
    act: HashMap<u8, usize>,
    inst: HashMap<u8, usize>,
) -> (Option<ActivationMethod>, Option<InstrumentType>) {
    let activation = act
        .into_iter()
        .max_by_key(|&(k, n)| (n, std::cmp::Reverse(k)))
        .map(|(k, _)| match k {
            x if x == ActivationMethod::CID as u8 => ActivationMethod::CID,
            x if x == ActivationMethod::HCD as u8 => ActivationMethod::HCD,
            x if x == ActivationMethod::ETD as u8 => ActivationMethod::ETD,
            x if x == ActivationMethod::PQD as u8 => ActivationMethod::PQD,
            x if x == ActivationMethod::UVPD as u8 => ActivationMethod::UVPD,
            x if x == ActivationMethod::EThcD as u8 => ActivationMethod::EThcD,
            _ => ActivationMethod::HCD,
        });
    let instrument = inst
        .into_iter()
        .max_by_key(|&(k, n)| (n, std::cmp::Reverse(k)))
        .map(|(k, _)| match k {
            x if x == InstrumentType::LowRes as u8 => InstrumentType::LowRes,
            x if x == InstrumentType::TOF as u8 => InstrumentType::TOF,
            x if x == InstrumentType::HighRes as u8 => InstrumentType::HighRes,
            x if x == InstrumentType::OrbitrapAstral as u8 => InstrumentType::OrbitrapAstral,
            x if x == InstrumentType::TimsTOF as u8 => InstrumentType::TimsTOF,
            _ => InstrumentType::QExactive,
        });
    (activation, instrument)
}

fn warn_mixed_electron_activation(act: &HashMap<u8, usize>, seen: usize) {
    let n = |m: ActivationMethod| act.get(&(m as u8)).copied().unwrap_or(0);
    let (etd, ethcd) = (n(ActivationMethod::ETD), n(ActivationMethod::EThcD));
    if etd > 0 && ethcd > 0 {
        log::warn!(
            "mixed electron activation in the first {seen} MS2 scans: EThcD/ETciD {ethcd}, \
             ETD {etd}. The run is routed by the dominant activation; ETD scans are scored \
             with the run's model (the fragment-ion index cannot serve c/z-only spectra). \
             Pass --model or --ethcd-activation etd to choose explicitly."
        );
    }
}

pub fn sample_ms2_meta<P: AsRef<Path>>(path: P, n: usize) -> Option<crate::Ms2MetaSample> {
    let mut reader = RawFileReader::open(path.as_ref()).ok()?;
    reader.set_signal_loading(false);
    let mut out = crate::Ms2MetaSample::default();
    for i in 0..reader.len() {
        if out.n_ms2 >= n {
            break;
        }
        let Some(raw) = reader.get(i) else {
            continue;
        };
        if raw.ms_level() != 2 {
            continue;
        }
        out.n_ms2 += 1;
        if let Some(p) = raw.precursor() {
            let ce = p.activation().collision_energy();
            if ce.is_finite() && ce > 0.0 {
                out.collision_energies.push(ce);
            }
            let w = p.isolation_window();
            let (lo, up) = full_width_isolation_bounds(
                w.lower(),
                w.upper(),
                isolation_width_trailer(&reader, raw.index(), raw.ms_level()),
            );
            if up > lo && lo > 0.0 {
                out.isolation_widths.push(up - lo);
            }
        }
        if let Some(cv) = faims_cv(&raw) {
            out.faims_cvs.push(cv);
        }
    }
    Some(out)
}

impl Iterator for ThermoRawReader {
    type Item = Result<Spectrum, ThermoParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let capture = self.capture(false);
        while let Some(scan) = self.next_fetched(capture) {
            if let Some(ms1) = scan.ms1 {
                self.run_ms1.push(ms1);
            }
            if let Some(p) = scan.product {
                self.product_scans.push(p);
                continue;
            }
            if let Some((spec, context)) = scan.ms2 {
                self.scan_contexts.extend(context);
                return Some(Ok(spec));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_ms3_precursor_is_the_last_reaction_of_its_filter() {
        let f = "ITMS + c NSI cv=-75.00 r d Full ms3 405.5227@cid25.00 375.5128@cid35.00 \
                 [99.0000-1137.0000]";
        assert_eq!(last_reaction_mz(f), Some(375.5128));
        assert_eq!(
            last_reaction_mz("FTMS + p NSI d Full ms2 784.37@hcd28.00 [120.00-2000.00]"),
            Some(784.37)
        );
        assert_eq!(
            last_reaction_mz("FTMS + p NSI Full ms [375.0000-1600.0000]"),
            None
        );
    }

    fn tally(pairs: &[(u8, usize)]) -> HashMap<u8, usize> {
        pairs.iter().copied().collect()
    }

    #[test]
    fn analyzer_is_kept_when_no_scan_names_an_activation() {
        let (activation, instrument) = dominant_acquisition(
            HashMap::new(),
            tally(&[(InstrumentType::OrbitrapAstral as u8, 64)]),
        );
        assert_eq!(activation, None);
        assert_eq!(instrument, Some(InstrumentType::OrbitrapAstral));
    }

    #[test]
    fn activation_is_kept_when_no_scan_names_an_analyzer() {
        let (activation, instrument) =
            dominant_acquisition(tally(&[(ActivationMethod::ETD as u8, 3)]), HashMap::new());
        assert_eq!(activation, Some(ActivationMethod::ETD));
        assert_eq!(instrument, None);
    }

    #[test]
    fn dominant_acquisition_is_the_most_frequent_of_each_tally() {
        let (activation, instrument) = dominant_acquisition(
            tally(&[
                (ActivationMethod::CID as u8, 10),
                (ActivationMethod::HCD as u8, 54),
            ]),
            tally(&[
                (InstrumentType::LowRes as u8, 10),
                (InstrumentType::QExactive as u8, 54),
            ]),
        );
        assert_eq!(activation, Some(ActivationMethod::HCD));
        assert_eq!(instrument, Some(InstrumentType::QExactive));
    }

    #[test]
    fn open_failure_names_the_file_and_the_runtime_only_when_it_may_be_missing() {
        let host = ThermoParseError::Open {
            path: "run.raw".into(),
            source: std::io::Error::other("One of the dependent libraries is missing."),
        }
        .to_string();
        assert!(host.contains("'run.raw'"), "{host}");
        assert!(host.contains("DOTNET_ROOT"), "{host}");

        for kind in [
            std::io::ErrorKind::NotFound,
            std::io::ErrorKind::InvalidData,
        ] {
            let file = ThermoParseError::Open {
                path: "run.raw".into(),
                source: std::io::Error::new(kind, "bad file"),
            }
            .to_string();
            assert!(
                file.contains("'run.raw'") && !file.contains("DOTNET_ROOT"),
                "{file}"
            );
        }
    }

    #[test]
    fn master_scan_links_to_the_native_id_of_that_scan() {
        assert_eq!(
            with_scan_number("controllerType=0 controllerNumber=1 scan=4780", 4773),
            "controllerType=0 controllerNumber=1 scan=4773"
        );
        assert_eq!(with_scan_number("scan=12", 3), "scan=3");
        assert!(is_master_scan_label("Master Scan Number"));
        assert!(is_master_scan_label("Master Scan Number:"));
        assert!(!is_master_scan_label("Master Index"));
    }

    #[test]
    fn sanitize_thermo_peaks_filters_and_sorts() {
        let pairs = vec![
            (300.0, 1000.0f32),
            (f64::NAN, 500.0),
            (0.0, 100.0),
            (-5.0, 100.0),
            (150.0, 2000.0),
            (400.0, f32::NAN),
            (450.0, f32::INFINITY),
            (500.0, -1.0),
            (250.0, 0.0),
        ];
        let out = super::sanitize_thermo_peaks(pairs);
        assert_eq!(out, vec![(150.0, 2000.0), (250.0, 0.0), (300.0, 1000.0)]);
    }

    #[test]
    fn dissociation_maps_to_activation() {
        assert_eq!(
            map_dissociation(DissociationMethod(1)),
            Some(ActivationMethod::CID)
        );
        assert_eq!(
            map_dissociation(DissociationMethod(2)),
            Some(ActivationMethod::HCD)
        );
        for v in [4u8, 8, 16] {
            assert_eq!(
                map_dissociation(DissociationMethod(v)),
                Some(ActivationMethod::ETD)
            );
        }
        for v in [5u8, 6, 9, 10] {
            assert_eq!(
                map_dissociation(DissociationMethod(v)),
                Some(ActivationMethod::EThcD)
            );
        }
        assert_eq!(map_dissociation(DissociationMethod(0)), None);
    }

    #[test]
    fn trailer_isolation_width_is_the_full_width() {
        let (lo, up) = full_width_isolation_bounds(498.4, 501.6, Some(1.6));
        assert!(
            (lo - 499.2).abs() < 1e-9 && (up - 500.8).abs() < 1e-9,
            "{lo} {up}"
        );
        let (lo, up) = full_width_isolation_bounds(498.6, 501.8, Some(1.6));
        assert!(
            (lo - 499.4).abs() < 1e-9 && (up - 501.0).abs() < 1e-9,
            "{lo} {up}"
        );
        assert_eq!(
            full_width_isolation_bounds(499.2, 500.8, Some(1.6)),
            (499.2, 500.8)
        );
        assert_eq!(
            full_width_isolation_bounds(498.4, 501.6, None),
            (498.4, 501.6)
        );
    }

    #[test]
    fn analyzer_maps_to_instrument() {
        assert_eq!(map_analyzer(MassAnalyzer(1)), Some(InstrumentType::LowRes));
        assert_eq!(map_analyzer(MassAnalyzer(4)), Some(InstrumentType::TOF));
        assert_eq!(
            map_analyzer(MassAnalyzer(5)),
            Some(InstrumentType::QExactive)
        );
        assert_eq!(
            map_analyzer(MassAnalyzer(7)),
            Some(InstrumentType::OrbitrapAstral)
        );
        assert_eq!(map_analyzer(MassAnalyzer(0)), None);
    }
}
