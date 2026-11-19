use std::collections::HashMap;
use std::io::BufRead;

use base64::{engine::general_purpose::STANDARD, Engine as _};
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::ZlibDecoder;
use quick_xml::{events::Event, Reader};

use model::scan::{Ms1Scan, NoiseProfile, ProductScan, ScanContext};
use model::{ActivationMethod, InstrumentType, Spectrum};

const CV_ANALYZER_ION_TRAP: &str = "MS:1000264";
const CV_ANALYZER_QUAD_ION_TRAP: &str = "MS:1000082";
const CV_ANALYZER_RADIAL_LIT: &str = "MS:1000083";
const CV_ANALYZER_LINEAR_ION_TRAP: &str = "MS:1000291";
const CV_ANALYZER_ORBITRAP: &str = "MS:1000484";
const CV_ANALYZER_FTICR: &str = "MS:1000079";
const CV_ANALYZER_TOF: &str = "MS:1000084";

const CV_MODEL_Q_EXACTIVE: &str = "MS:1001911";
const CV_MODEL_Q_EXACTIVE_HF: &str = "MS:1002523";
const CV_MODEL_Q_EXACTIVE_HF_X: &str = "MS:1002634";
const CV_MODEL_Q_EXACTIVE_PLUS: &str = "MS:1002877";
const CV_MODEL_ORBITRAP_FUSION: &str = "MS:1002416";
const CV_MODEL_ORBITRAP_ASTRAL: &str = "MS:1003378";

const CV_MS_LEVEL: &str = "MS:1000511";
const CV_SCAN_TIME: &str = "MS:1000016";
const CV_PROFILE: &str = "MS:1000128";
const CV_SELECTED_ION_MZ: &str = "MS:1000744";
const CV_ISOLATION_LOWER_OFFSET: &str = "MS:1000828";
const CV_ISOLATION_TARGET_MZ: &str = "MS:1000827";
const CV_ISOLATION_UPPER_OFFSET: &str = "MS:1000829";
const CV_SCAN_WINDOW_LOWER: &str = "MS:1000501";
const CV_SCAN_WINDOW_UPPER: &str = "MS:1000500";
const CV_FAIMS_CV: &str = "MS:1001581";
const CV_FILTER_STRING: &str = "MS:1000512";
const CV_MZ_PLAIN: &str = "MS:1000040";
const CV_CHARGE_STATE: &str = "MS:1000041";
const CV_PEAK_INTENSITY: &str = "MS:1000042";
const CV_MZ_ARRAY: &str = "MS:1000514";
const CV_INTENSITY_ARRAY: &str = "MS:1000515";
const CV_NOISE_MZ_ARRAY: &str = "MS:1002743";
const CV_NOISE_ARRAY: &str = "MS:1002744";
const CV_64BIT: &str = "MS:1000523";
const CV_32BIT: &str = "MS:1000521";
const CV_ZLIB: &str = "MS:1000574";
const CV_TRUNCATION_ZLIB: &str = "MS:1003088";

const UNSUPPORTED_COMPRESSION: &[(&str, &str)] = &[
    ("MS:1002312", "MS-Numpress linear prediction"),
    ("MS:1002313", "MS-Numpress positive integer"),
    ("MS:1002314", "MS-Numpress short logged float"),
    ("MS:1002746", "MS-Numpress linear prediction + zlib"),
    ("MS:1002747", "MS-Numpress positive integer + zlib"),
    ("MS:1002748", "MS-Numpress short logged float + zlib"),
    ("MS:1003089", "truncation, delta prediction + zlib"),
    ("MS:1003090", "truncation, linear prediction + zlib"),
    ("MS:1003780", "zstd"),
    ("MS:1003781", "byte-shuffled zstd"),
    ("MS:1003782", "dictionary-encoded zstd"),
    ("MS:1003783", "MS-Numpress linear prediction + zstd"),
    ("MS:1003784", "MS-Numpress positive integer + zstd"),
    ("MS:1003785", "MS-Numpress short logged float + zstd"),
    ("MS:1003826", "coordinate grid encoding"),
];

const CV_CID: &str = "MS:1000133";
const CV_HCD: &str = "MS:1000422";
const CV_ETD: &str = "MS:1000598";
const CV_PQD: &str = "MS:1000599";
const CV_UVPD: &str = "MS:1000435";
const CV_ECD: &str = "MS:1000250";

const CV_TRAP_CID: &str = "MS:1002472";
const CV_HIGH_ENERGY_BEAM_CID: &str = "MS:1002481";
const CV_UVPD_CHILD: &str = "MS:1003246";

const CV_SUPPLEMENTAL_HCD: &str = "MS:1002678";
const CV_SUPPLEMENTAL_CID: &str = "MS:1002679";

const CV_ETHCD: &str = "MS:1002631";
const CV_ETCID: &str = "MS:1003182";

const CV_UNIT_MINUTE: &str = "UO:0000031";
const CV_UNIT_MINUTE_LEGACY: &str = "MS:1000038";

fn is_minute_unit(unit_accession: &str, unit_name: &str) -> bool {
    match unit_accession {
        CV_UNIT_MINUTE | CV_UNIT_MINUTE_LEGACY => true,
        "" => matches!(
            unit_name.trim().to_ascii_lowercase().as_str(),
            "minute" | "min" | "minutes"
        ),
        _ => false,
    }
}

#[derive(Debug, thiserror::Error)]
pub enum MzMLParseError {
    #[error("XML parse error: {0}")]
    Xml(#[from] quick_xml::Error),

    #[error("base64 decode error: {0}")]
    Base64(#[from] base64::DecodeError),

    #[error("zlib decode error: {0}")]
    Zlib(std::io::Error),

    #[error("mzML structure: {0}")]
    Structure(String),

    #[error("mismatched binary array lengths: m/z {mz_len} vs intensity {int_len}")]
    LengthMismatch { mz_len: usize, int_len: usize },

    #[error(
        "unsupported binary encoding {accession} ({name}); re-convert without it \
         (e.g. msconvert without --numpress*/--zstd)"
    )]
    UnsupportedCompression {
        accession: &'static str,
        name: &'static str,
    },

    #[error("truncated binary array: decoded {decoded_len} bytes is not a multiple of the {elem_bytes}-byte element width")]
    TruncatedBinary {
        decoded_len: usize,
        elem_bytes: usize,
    },
}

impl From<std::io::Error> for MzMLParseError {
    fn from(e: std::io::Error) -> Self {
        MzMLParseError::Zlib(e)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
enum State {
    #[default]
    Outside,
    Spectrum,
    Scan,
    SelectedIon,
    IsolationWindow,
    Activation,
    ProductList,
    BinaryDataArray,
    Binary,
}

#[derive(Debug)]
struct BinaryArrayCtx {
    is_mz: bool,
    is_intensity: bool,
    is_noise_mz: bool,
    is_noise: bool,
    precision_bits: u8,
    zlib: bool,
    unsupported: Option<(&'static str, &'static str)>,
    b64_text: String,
}

impl BinaryArrayCtx {
    fn new() -> Self {
        BinaryArrayCtx {
            is_mz: false,
            is_intensity: false,
            is_noise_mz: false,
            is_noise: false,
            precision_bits: 64,
            zlib: false,
            unsupported: None,
            b64_text: String::new(),
        }
    }
}

pub(crate) fn log_ethcd_once() {
    use std::sync::atomic::{AtomicBool, Ordering};
    static WARNED: AtomicBool = AtomicBool::new(false);
    if WARNED
        .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
        .is_ok()
    {
        if ethcd_as_etd() {
            log::info!(
                "EThcD/ETciD detected; --ethcd-activation etd is set, so these spectra \
                 stay classified as ETD (c/z stack active)."
            );
        } else {
            log::info!(
                "EThcD/ETciD detected (electron transfer + supplemental collisional \
                 activation); these spectra are labelled EThcD and scored with b/y and c/z \
                 ions. IN --glyco the c/z stack (c/z generation and scoring, the c/z \
                 truncation gate, the ETD DB fallback, --glyco-hcd-pair) runs only on ETD \
                 spectra; pass --ethcd-activation etd to label these spectra ETD."
            );
        }
    }
}

fn log_profile_centroid_once() {
    use std::sync::atomic::{AtomicBool, Ordering};
    static WARNED: AtomicBool = AtomicBool::new(false);
    if WARNED
        .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
        .is_ok()
    {
        log::info!(
            "profile-mode spectra detected (MS:1000128) — centroiding peaks on the fly \
             (local-max + intensity-weighted m/z). For best results pre-centroid during \
             conversion (e.g. the converter's peak-picking option)."
        );
    }
}

fn centroid_profile(profile: &[(f64, f32)]) -> Vec<(f64, f32)> {
    let n = profile.len();
    if n < 3 {
        return profile.to_vec();
    }
    let mut out: Vec<(f64, f32)> = Vec::new();
    let mut i = 1usize;
    while i < n - 1 {
        let (prev, cur, next) = (profile[i - 1].1, profile[i].1, profile[i + 1].1);
        if cur > 0.0 && cur >= prev && cur > next {
            let mut l = i;
            while l > 0 && profile[l - 1].1 <= profile[l].1 && profile[l - 1].1 > 0.0 {
                l -= 1;
            }
            let mut r = i;
            while r + 1 < n && profile[r + 1].1 <= profile[r].1 && profile[r + 1].1 > 0.0 {
                r += 1;
            }
            let mut wsum = 0f64;
            let mut msum = 0f64;
            let mut apex = 0f32;
            for &(m, it) in &profile[l..=r] {
                wsum += it as f64;
                msum += m * it as f64;
                if it > apex {
                    apex = it;
                }
            }
            if wsum > 0.0 {
                out.push((msum / wsum, apex));
            }
            i = r + 1;
        } else {
            i += 1;
        }
    }
    out
}

#[derive(Debug, Default)]
struct SpectrumBuilder {
    id: String,
    index: Option<i32>,
    ms_level: Option<u32>,
    rt_seconds: Option<f64>,
    precursor_mz: Option<f64>,
    monoisotopic_mz_override: Option<f64>,
    precursor_charge: Option<i32>,
    precursor_intensity: Option<f32>,
    selected_ions: u32,
    activation_method: Option<ActivationMethod>,
    act_saw_electron: bool,
    act_saw_collisional: bool,
    filter_supplemental: bool,
    is_profile: bool,
    isolation_lower_offset: Option<f64>,
    isolation_upper_offset: Option<f64>,
    isolation_target_mz: Option<f64>,
    parent_ref: Option<String>,
    scan_window_lower: Option<f64>,
    scan_window_upper: Option<f64>,
    faims_cv: Option<f32>,
    mz_array: Option<Vec<f64>>,
    intensity_array: Option<Vec<f64>>,
    noise_mz_array: Option<Vec<f64>>,
    noise_array: Option<Vec<f64>>,
}

impl SpectrumBuilder {
    fn take_noise(&mut self) -> Option<NoiseProfile> {
        let mz = self.noise_mz_array.take()?;
        let noise = self.noise_array.take()?;
        if mz.len() != noise.len() {
            return None;
        }
        NoiseProfile::reporter_region(mz.into_iter().zip(noise.into_iter().map(|n| n as f32)))
    }
}

fn activation_user_param_accession(name: &str) -> Option<&'static str> {
    match name.trim().to_ascii_lowercase().as_str() {
        "cid" | "collision-induced dissociation" => Some(CV_CID),
        "hcd" | "beam-type collision-induced dissociation" => Some(CV_HCD),
        "etd" | "electron transfer dissociation" => Some(CV_ETD),
        _ => None,
    }
}

struct CvParamInfo {
    accession: String,
    value: String,
    unit_accession: String,
    unit_name: String,
}

impl CvParamInfo {
    fn from_bytes_start(e: &quick_xml::events::BytesStart<'_>) -> Option<Self> {
        let accession = attr_str(e, b"accession")?;
        let value = attr_str(e, b"value").unwrap_or_default();
        let unit_accession = attr_str(e, b"unitAccession").unwrap_or_default();
        let unit_name = attr_str(e, b"unitName").unwrap_or_default();
        Some(CvParamInfo {
            accession,
            value,
            unit_accession,
            unit_name,
        })
    }
}

#[derive(Debug, Default, Clone)]
pub struct Ms1Link {
    pub ms1_peaks: Vec<Vec<(f64, f32)>>,
    pub ms2_to_ms1: Vec<Option<usize>>,
}

pub struct MzMLReader<R: BufRead> {
    xml: Reader<R>,
    buf: Vec<u8>,
    ms_level_min: u32,
    ms_level_max: u32,
    state: State,
    current: Option<SpectrumBuilder>,
    binary_ctx: Option<BinaryArrayCtx>,
    done: bool,
    capture_ms1: bool,
    captured_ms1: Vec<Vec<(f64, f32)>>,
    ms1_tracker: crate::ms1_link::Ms1Tracker,
    last_faims_cv: Option<f32>,
    capture_run_ms1: bool,
    run_ms1: Vec<Ms1Scan>,
    capture_product_scans: bool,
    product_scans: Vec<ProductScan>,
    capture_scan_context: bool,
    scan_contexts: Vec<ScanContext>,
    capture_noise: bool,
    strict: bool,
    skipped: usize,
    skipped_errors: Vec<String>,
}

impl<R: BufRead> MzMLReader<R> {
    pub fn new(reader: R) -> Self {
        let mut xml = Reader::from_reader(reader);
        xml.trim_text(true);
        Self {
            xml,
            buf: Vec::with_capacity(4096),
            ms_level_min: 2,
            ms_level_max: 2,
            state: State::Outside,
            current: None,
            binary_ctx: None,
            done: false,
            capture_ms1: false,
            captured_ms1: Vec::new(),
            ms1_tracker: Default::default(),
            last_faims_cv: None,
            capture_run_ms1: false,
            run_ms1: Vec::new(),
            capture_product_scans: false,
            product_scans: Vec::new(),
            capture_scan_context: false,
            scan_contexts: Vec::new(),
            capture_noise: false,
            strict: false,
            skipped: 0,
            skipped_errors: Vec::new(),
        }
    }

    pub fn strict(mut self) -> Self {
        self.strict = true;
        self
    }

    pub fn skipped_count(&self) -> usize {
        self.skipped
    }

    pub fn skipped_errors(&self) -> &[String] {
        &self.skipped_errors
    }

    pub fn with_ms_level_range(mut self, min: u32, max: u32) -> Self {
        self.ms_level_min = min;
        self.ms_level_max = max;
        self
    }

    pub fn with_ms1_capture(mut self, capture: bool) -> Self {
        self.capture_ms1 = capture;
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

    fn finish_spectrum(&self, sb: SpectrumBuilder) -> Result<Option<Spectrum>, MzMLParseError> {
        let level = sb.ms_level.unwrap_or(0);
        if level < self.ms_level_min || level > self.ms_level_max {
            return Ok(None);
        }
        Self::build_spectrum(sb)
    }

    fn build_spectrum(sb: SpectrumBuilder) -> Result<Option<Spectrum>, MzMLParseError> {
        let precursor_mz = match (sb.monoisotopic_mz_override, sb.precursor_mz) {
            (Some(m), _) => m,
            (None, Some(v)) => v,
            (None, None) => return Ok(None),
        };

        if !precursor_mz.is_finite() || precursor_mz <= 0.0 {
            return Ok(None);
        }

        let peaks = Self::build_peaks(sb.mz_array, sb.intensity_array, sb.is_profile)?;
        let scan = scan_number(&sb.id, sb.index);

        Ok(Some(Spectrum {
            title: sb.id,
            precursor_mz,
            precursor_charge: sb.precursor_charge,
            precursor_intensity: sb.precursor_intensity,
            rt_seconds: sb.rt_seconds,
            scan,
            peaks,
            activation_method: sb.activation_method,
            isolation_lower_offset: sb.isolation_lower_offset,
            isolation_upper_offset: sb.isolation_upper_offset,
            isolation_target_mz: sb.isolation_target_mz,
            mobility_window: None,
        }))
    }

    fn build_peaks(
        mz_array: Option<Vec<f64>>,
        intensity_array: Option<Vec<f64>>,
        is_profile: bool,
    ) -> Result<Vec<(f64, f32)>, MzMLParseError> {
        let mz_vals = mz_array.unwrap_or_default();
        let int_vals = intensity_array.unwrap_or_default();

        if mz_vals.len() != int_vals.len() {
            return Err(MzMLParseError::LengthMismatch {
                mz_len: mz_vals.len(),
                int_len: int_vals.len(),
            });
        }

        let mut peaks: Vec<(f64, f32)> = mz_vals
            .into_iter()
            .zip(int_vals)
            .map(|(mz, inten)| (mz, inten as f32))
            .filter(|&(mz, inten)| mz.is_finite() && mz > 0.0 && inten.is_finite() && inten >= 0.0)
            .collect();

        if !peaks.windows(2).all(|w| w[0].0 <= w[1].0) {
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        }

        if is_profile {
            log_profile_centroid_once();
            peaks = centroid_profile(&peaks);
        }

        Ok(peaks)
    }

    fn in_first_selected_ion(&self) -> bool {
        self.state == State::SelectedIon
            && self
                .current
                .as_ref()
                .is_some_and(|sb| sb.selected_ions == 1)
    }

    fn apply_cv_param(&mut self, cv: CvParamInfo) {
        let (state, first_ion) = (self.state, self.in_first_selected_ion());
        let (in_scan, activation) = (
            matches!(state, State::Scan | State::Spectrum),
            state == State::Activation,
        );
        let (sb, ctx) = (self.current.as_mut(), self.binary_ctx.as_mut());
        let value = cv.value.as_str();
        match cv.accession.as_str() {
            CV_MS_LEVEL => {
                if let (Ok(lvl), Some(sb)) = (value.parse::<u32>(), sb) {
                    sb.ms_level = Some(lvl);
                }
            }

            CV_PROFILE => {
                if let Some(sb) = sb {
                    sb.is_profile = true;
                }
            }

            CV_SCAN_TIME if in_scan => {
                if let (Ok(t), Some(sb)) = (value.parse::<f64>(), sb) {
                    let minutes = is_minute_unit(&cv.unit_accession, &cv.unit_name);
                    sb.rt_seconds = Some(if minutes { t * 60.0 } else { t });
                }
            }

            CV_SELECTED_ION_MZ if first_ion => {
                if let (Ok(mz), Some(sb)) = (value.parse::<f64>(), sb) {
                    sb.precursor_mz = Some(mz);
                }
            }

            CV_MZ_PLAIN if first_ion => {
                if let (Ok(mz), Some(sb)) = (value.parse::<f64>(), sb) {
                    sb.precursor_mz.get_or_insert(mz);
                }
            }

            CV_CHARGE_STATE if first_ion => {
                if let (Ok(z), Some(sb)) = (value.parse::<i32>(), sb) {
                    sb.precursor_charge = model::spectrum::known_precursor_charge(i64::from(z));
                }
            }

            CV_PEAK_INTENSITY if first_ion => {
                if let (Ok(inten), Some(sb)) = (value.parse::<f32>(), sb) {
                    sb.precursor_intensity = Some(inten);
                }
            }

            CV_ISOLATION_LOWER_OFFSET if state == State::IsolationWindow => {
                if let (Ok(off), Some(sb)) = (value.parse::<f64>(), sb) {
                    sb.isolation_lower_offset = Some(off);
                }
            }
            CV_ISOLATION_UPPER_OFFSET if state == State::IsolationWindow => {
                if let (Ok(off), Some(sb)) = (value.parse::<f64>(), sb) {
                    sb.isolation_upper_offset = Some(off);
                }
            }
            CV_ISOLATION_TARGET_MZ if state == State::IsolationWindow => {
                if let (Ok(mz), Some(sb)) = (value.parse::<f64>(), sb) {
                    sb.isolation_target_mz = (mz.is_finite() && mz > 0.0).then_some(mz);
                }
            }

            CV_SCAN_WINDOW_LOWER if state == State::Scan => {
                if let (Ok(v), Some(sb)) = (value.parse::<f64>(), sb) {
                    sb.scan_window_lower.get_or_insert(v);
                }
            }
            CV_SCAN_WINDOW_UPPER if state == State::Scan => {
                if let (Ok(v), Some(sb)) = (value.parse::<f64>(), sb) {
                    sb.scan_window_upper.get_or_insert(v);
                }
            }
            CV_FILTER_STRING if in_scan => {
                if let Some(sb) = sb {
                    sb.filter_supplemental = filter_has_supplemental_collision(value);
                }
            }
            CV_FAIMS_CV if in_scan => {
                if let (Ok(v), Some(sb)) = (value.parse::<f32>(), sb) {
                    sb.faims_cv.get_or_insert(v);
                }
            }

            CV_CID | CV_TRAP_CID if activation => {
                if let Some(sb) = sb {
                    sb.act_saw_collisional = true;
                    sb.activation_method.get_or_insert(ActivationMethod::CID);
                }
            }
            CV_HCD | CV_HIGH_ENERGY_BEAM_CID if activation => {
                if let Some(sb) = sb {
                    sb.act_saw_collisional = true;
                    sb.activation_method.get_or_insert(ActivationMethod::HCD);
                }
            }
            CV_ETD if activation => {
                if let Some(sb) = sb {
                    sb.act_saw_electron = true;
                    sb.activation_method = Some(ActivationMethod::ETD);
                }
            }
            CV_SUPPLEMENTAL_HCD | CV_SUPPLEMENTAL_CID if activation => {
                if let Some(sb) = sb {
                    sb.act_saw_collisional = true;
                }
            }
            CV_ETHCD | CV_ETCID if activation => {
                if let Some(sb) = sb {
                    sb.act_saw_electron = true;
                    sb.act_saw_collisional = true;
                    sb.activation_method = Some(ActivationMethod::ETD);
                }
            }
            CV_ECD if activation => {
                if let Some(sb) = sb {
                    sb.act_saw_electron = true;
                    sb.activation_method.get_or_insert(ActivationMethod::ETD);
                }
            }
            CV_PQD if activation => {
                if let Some(sb) = sb {
                    sb.activation_method.get_or_insert(ActivationMethod::PQD);
                }
            }
            CV_UVPD | CV_UVPD_CHILD if activation => {
                if let Some(sb) = sb {
                    sb.activation_method.get_or_insert(ActivationMethod::UVPD);
                }
            }

            acc if state == State::BinaryDataArray => {
                let Some(ctx) = ctx else { return };
                match acc {
                    CV_MZ_ARRAY => ctx.is_mz = true,
                    CV_INTENSITY_ARRAY => ctx.is_intensity = true,
                    CV_NOISE_MZ_ARRAY => ctx.is_noise_mz = true,
                    CV_NOISE_ARRAY => ctx.is_noise = true,
                    CV_64BIT => ctx.precision_bits = 64,
                    CV_32BIT => ctx.precision_bits = 32,
                    CV_ZLIB | CV_TRUNCATION_ZLIB => ctx.zlib = true,
                    _ => {
                        if let Some(&known) =
                            UNSUPPORTED_COMPRESSION.iter().find(|(a, _)| *a == acc)
                        {
                            ctx.unsupported = Some(known);
                        }
                    }
                }
            }

            _ => {}
        }
    }

    fn end_spectrum(&mut self) -> Result<Option<Spectrum>, MzMLParseError> {
        self.state = State::Outside;
        let Some(mut sb) = self.current.take() else {
            return Ok(None);
        };
        if self.capture_run_ms1 && sb.ms_level == Some(1) {
            let emits_ms1 = self.ms_level_min <= 1 && !self.capture_ms1;
            let (mz, inten) = if emits_ms1 {
                (sb.mz_array.clone(), sb.intensity_array.clone())
            } else {
                (sb.mz_array.take(), sb.intensity_array.take())
            };
            let peaks = Self::build_peaks(mz, inten, sb.is_profile)?;
            self.run_ms1.push(Ms1Scan {
                rt: sb.rt_seconds.unwrap_or(0.0),
                peaks,
                id: sb.id.clone(),
                faims_cv: sb.faims_cv,
                mobility: Vec::new(),
            });
            if !emits_ms1 && !self.capture_ms1 {
                return Ok(None);
            }
            if self.capture_ms1 {
                let last = self.run_ms1.last().expect("just pushed");
                self.captured_ms1.push(last.peaks.clone());
                self.ms1_tracker
                    .push(self.captured_ms1.len() - 1, sb.faims_cv);
                return Ok(None);
            }
        }
        if self.capture_product_scans && sb.ms_level.unwrap_or(0) >= 3 {
            let level = sb.ms_level.unwrap_or(3) as u8;
            let noise = sb.take_noise();
            let peaks =
                Self::build_peaks(sb.mz_array.take(), sb.intensity_array.take(), sb.is_profile)?;
            self.product_scans.push(ProductScan {
                scan: scan_number(&sb.id, sb.index),
                id: sb.id,
                rt_seconds: sb.rt_seconds,
                ms_level: level,
                parent_id: sb.parent_ref,
                precursor_mz: sb
                    .monoisotopic_mz_override
                    .or(sb.precursor_mz)
                    .unwrap_or(0.0),
                scan_window: sb.scan_window_lower.zip(sb.scan_window_upper),
                peaks,
                noise,
            });
            return Ok(None);
        }
        if self.capture_ms1 && sb.ms_level == Some(1) {
            let peaks = Self::build_peaks(sb.mz_array, sb.intensity_array, sb.is_profile)?;
            self.captured_ms1.push(peaks);
            self.ms1_tracker
                .push(self.captured_ms1.len() - 1, sb.faims_cv);
            return Ok(None);
        }
        let noise = sb.take_noise();
        let context = self.capture_scan_context.then(|| ScanContext {
            id: sb.id.clone(),
            parent_id: sb.parent_ref.clone(),
            scan_window: sb.scan_window_lower.zip(sb.scan_window_upper),
            faims_cv: sb.faims_cv,
            noise,
            ion_mobility: None,
        });
        let faims_cv = sb.faims_cv;
        if let Some(s) = self.finish_spectrum(sb)? {
            self.last_faims_cv = faims_cv;
            if let Some(c) = context {
                self.scan_contexts.push(c);
            }
            return Ok(Some(s));
        }
        Ok(None)
    }

    fn pump(&mut self) -> Result<Option<Spectrum>, MzMLParseError> {
        loop {
            self.buf.clear();
            let event = self.xml.read_event_into(&mut self.buf)?;

            match event {
                Event::Eof => {
                    self.done = true;
                    return Ok(None);
                }

                Event::Start(ref e) => {
                    let tag = e.local_name().as_ref().to_owned();
                    match tag.as_slice() {
                        b"spectrum" => {
                            let id = attr_str(e, b"id").unwrap_or_default();
                            let index = attr_str(e, b"index").and_then(|v| v.parse().ok());
                            self.current = Some(SpectrumBuilder {
                                id,
                                index,
                                ..Default::default()
                            });
                            self.state = State::Spectrum;
                        }
                        b"scan" if self.state == State::Spectrum => {
                            self.state = State::Scan;
                        }
                        b"selectedIon" if self.state == State::Spectrum => {
                            self.state = State::SelectedIon;
                            if let Some(sb) = self.current.as_mut() {
                                sb.selected_ions += 1;
                            }
                        }
                        b"precursor" if self.state == State::Spectrum => {
                            if let Some(sb) = self.current.as_mut() {
                                if sb.parent_ref.is_none() {
                                    sb.parent_ref = attr_str(e, b"spectrumRef");
                                }
                            }
                        }
                        b"productList" if self.state == State::Spectrum => {
                            self.state = State::ProductList;
                        }
                        b"isolationWindow" if self.state == State::Spectrum => {
                            self.state = State::IsolationWindow;
                        }
                        b"activation" if self.state == State::Spectrum => {
                            self.state = State::Activation;
                            if let Some(sb) = self.current.as_mut() {
                                sb.act_saw_electron = false;
                                sb.act_saw_collisional = false;
                            }
                        }
                        b"binaryDataArray" if self.state == State::Spectrum => {
                            let dropped = self
                                .current
                                .as_ref()
                                .and_then(|sb| sb.ms_level)
                                .is_some_and(|level| {
                                    (level < self.ms_level_min || level > self.ms_level_max)
                                        && !((self.capture_ms1 || self.capture_run_ms1)
                                            && level == 1)
                                        && !(self.capture_product_scans && level >= 3)
                                });
                            self.binary_ctx = (!dropped).then(BinaryArrayCtx::new);
                            self.state = State::BinaryDataArray;
                        }
                        b"binary" if self.state == State::BinaryDataArray => {
                            self.state = State::Binary;
                        }
                        _ => {}
                    }
                }

                Event::Empty(ref e) => {
                    let tag = e.local_name().as_ref().to_owned();
                    if tag == b"cvParam" {
                        if let Some(cv) = CvParamInfo::from_bytes_start(e) {
                            self.apply_cv_param(cv);
                        }
                    } else if tag == b"userParam" && self.state == State::Activation {
                        if let Some(accession) =
                            attr_str(e, b"name").and_then(|n| activation_user_param_accession(&n))
                        {
                            self.apply_cv_param(CvParamInfo {
                                accession: accession.to_string(),
                                value: String::new(),
                                unit_accession: String::new(),
                                unit_name: String::new(),
                            });
                        }
                    } else if tag == b"userParam"
                        && matches!(self.state, State::Scan | State::Spectrum)
                    {
                        if let (Some(name), Some(val)) =
                            (attr_str(e, b"name"), attr_str(e, b"value"))
                        {
                            let normalized = name.to_lowercase();
                            if normalized.contains("monoisotopic m/z")
                                || normalized.contains("monoisotopic mz")
                            {
                                if let Ok(mz) = val.parse::<f64>() {
                                    if mz > 0.0 {
                                        if let Some(sb) = self.current.as_mut() {
                                            sb.monoisotopic_mz_override = Some(mz);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                Event::Text(ref e) if self.state == State::Binary && self.binary_ctx.is_some() => {
                    let chunk = e.unescape()?;
                    if let Some(ctx) = self.binary_ctx.as_mut() {
                        ctx.b64_text.push_str(chunk.as_ref());
                    }
                }

                Event::End(ref e) => {
                    let tag = e.local_name().as_ref().to_owned();
                    match tag.as_slice() {
                        b"spectrum" => {
                            if let Some(s) = self.end_spectrum()? {
                                return Ok(Some(s));
                            }
                        }
                        b"scan" if self.state == State::Scan => {
                            self.state = State::Spectrum;
                        }
                        b"selectedIon" if self.state == State::SelectedIon => {
                            self.state = State::Spectrum;
                        }
                        b"isolationWindow" if self.state == State::IsolationWindow => {
                            self.state = State::Spectrum;
                        }
                        b"productList" if self.state == State::ProductList => {
                            self.state = State::Spectrum;
                        }
                        b"activation" if self.state == State::Activation => {
                            if let Some(sb) = self.current.as_mut() {
                                if sb.act_saw_electron
                                    && (sb.act_saw_collisional || sb.filter_supplemental)
                                {
                                    log_ethcd_once();
                                    if !ethcd_as_etd() {
                                        sb.activation_method = Some(ActivationMethod::EThcD);
                                    }
                                }
                            }
                            self.state = State::Spectrum;
                        }
                        b"binary" if self.state == State::Binary => {
                            self.state = State::BinaryDataArray;
                        }
                        b"binaryDataArray" if self.state == State::BinaryDataArray => {
                            if let Some(ctx) = self.binary_ctx.take() {
                                let wanted = ctx.is_mz
                                    || ctx.is_intensity
                                    || (self.capture_noise && (ctx.is_noise_mz || ctx.is_noise));
                                if wanted {
                                    let vals = decode_binary_array(&ctx)?;
                                    if let Some(sb) = self.current.as_mut() {
                                        if ctx.is_mz {
                                            sb.mz_array = Some(vals);
                                        } else if ctx.is_intensity {
                                            sb.intensity_array = Some(vals);
                                        } else if ctx.is_noise_mz {
                                            sb.noise_mz_array = Some(vals);
                                        } else {
                                            sb.noise_array = Some(vals);
                                        }
                                    }
                                }
                            }
                            self.state = State::Spectrum;
                        }
                        _ => {}
                    }
                }

                _ => {}
            }
        }
    }

    pub fn read_with_ms1(mut self) -> std::io::Result<(Vec<Spectrum>, Ms1Link)> {
        if self.capture_ms1 {
            self.ms_level_min = 1;
        }

        let mut spectra: Vec<Spectrum> = Vec::new();
        let mut ms2_to_ms1: Vec<Option<usize>> = Vec::new();

        loop {
            match self.pump() {
                Ok(Some(s)) => {
                    ms2_to_ms1.push(self.ms1_tracker.link(self.last_faims_cv));
                    spectra.push(s);
                }
                Ok(None) => break,
                Err(_e) => match self.resync_to_next_spectrum() {
                    Ok(true) => continue,
                    Ok(false) | Err(_) => break,
                },
            }
        }

        let link = Ms1Link {
            ms1_peaks: std::mem::take(&mut self.captured_ms1),
            ms2_to_ms1,
        };
        Ok((spectra, link))
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
        self.capture_ms1 = true;
        self.ms_level_min = 1;

        let mut err_count = 0usize;
        let mut first_errors: Vec<String> = Vec::new();
        let mut total = 0usize;

        let mut chunk: Vec<Spectrum> = Vec::with_capacity(chunk_size);
        let mut chunk_ms1: Vec<Vec<(f64, f32)>> = Vec::new();
        let mut links: Vec<Option<usize>> = Vec::with_capacity(chunk_size);
        // Captured MS1 index -> its copy in this chunk.
        let mut in_chunk: Vec<(usize, usize)> = Vec::new();
        let mut kept = 0usize;

        loop {
            if total >= cap {
                break;
            }
            match self.pump() {
                Ok(Some(ms2)) => {
                    if self.captured_ms1.len() > kept {
                        // New MS1 scans arrived: keep only those a later MS2 can link to.
                        let live = self.ms1_tracker.live();
                        let mut all = std::mem::take(&mut self.captured_ms1);
                        self.captured_ms1 =
                            live.iter().map(|&i| std::mem::take(&mut all[i])).collect();
                        in_chunk = in_chunk
                            .into_iter()
                            .filter_map(|(g, c)| live.binary_search(&g).ok().map(|n| (n, c)))
                            .collect();
                        self.ms1_tracker.renumber_to_live();
                        kept = self.captured_ms1.len();
                    }
                    let link = self.ms1_tracker.link(self.last_faims_cv).map(|g| {
                        match in_chunk.iter().find(|(k, _)| *k == g) {
                            Some(&(_, c)) => c,
                            None => {
                                chunk_ms1.push(self.captured_ms1[g].clone());
                                in_chunk.push((g, chunk_ms1.len() - 1));
                                chunk_ms1.len() - 1
                            }
                        }
                    });
                    chunk.push(ms2);
                    links.push(link);
                    total += 1;

                    if chunk.len() >= chunk_size {
                        on_chunk(
                            std::mem::take(&mut chunk),
                            Ms1Link {
                                ms1_peaks: std::mem::take(&mut chunk_ms1),
                                ms2_to_ms1: std::mem::take(&mut links),
                            },
                        );
                        chunk = Vec::with_capacity(chunk_size);
                        chunk_ms1 = Vec::new();
                        links = Vec::with_capacity(chunk_size);
                        in_chunk.clear();
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    err_count += 1;
                    if first_errors.len() < 3 {
                        first_errors.push(format!("{e}"));
                    }
                    match self.resync_to_next_spectrum() {
                        Ok(true) => continue,
                        Ok(false) | Err(_) => break,
                    }
                }
            }
        }

        if !chunk.is_empty() {
            on_chunk(
                chunk,
                Ms1Link {
                    ms1_peaks: chunk_ms1,
                    ms2_to_ms1: links,
                },
            );
        }
        (err_count, first_errors)
    }

    fn resync_to_next_spectrum(&mut self) -> Result<bool, MzMLParseError> {
        self.current = None;
        self.binary_ctx = None;
        self.state = State::Outside;
        loop {
            self.buf.clear();
            match self.xml.read_event_into(&mut self.buf)? {
                Event::Eof => {
                    self.done = true;
                    return Ok(false);
                }
                Event::Start(ref e) if e.local_name().as_ref() == b"spectrum" => {
                    let id = attr_str(e, b"id").unwrap_or_default();
                    let index = attr_str(e, b"index").and_then(|v| v.parse().ok());
                    self.current = Some(SpectrumBuilder {
                        id,
                        index,
                        ..Default::default()
                    });
                    self.state = State::Spectrum;
                    return Ok(true);
                }
                _ => {}
            }
        }
    }
}

pub const MAX_SKIPPED_ERRORS: usize = 3;

impl<R: BufRead> Iterator for MzMLReader<R> {
    type Item = Result<Spectrum, MzMLParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            match self.pump() {
                Ok(Some(s)) => return Some(Ok(s)),
                Ok(None) => {
                    self.done = true;
                    return None;
                }
                Err(_e) if self.strict => {
                    self.done = true;
                    return Some(Err(_e));
                }
                Err(e) => {
                    self.skipped += 1;
                    if self.skipped_errors.len() < MAX_SKIPPED_ERRORS {
                        self.skipped_errors.push(e.to_string());
                    }
                    match self.resync_to_next_spectrum() {
                        Ok(true) => continue,
                        Ok(false) => {
                            self.done = true;
                            return None;
                        }
                        Err(e) => {
                            self.done = true;
                            return Some(Err(e));
                        }
                    }
                }
            }
        }
    }
}

const PARALLEL_PIECE_BYTES: u64 = 16 << 20;

fn next_spectrum_start(file: &std::fs::File, from: u64) -> std::io::Result<Option<u64>> {
    use std::io::{Read, Seek, SeekFrom};
    const PAT: &[u8] = b"<spectrum";
    let mut f = file;
    f.seek(SeekFrom::Start(from))?;
    let mut block = vec![0u8; 1 << 20];
    let keep = PAT.len();
    let mut carried = 0usize;
    let mut block_start = from;
    loop {
        let n = f.read(&mut block[carried..])?;
        let len = carried + n;
        if n == 0 {
            return Ok(None);
        }
        let hay = &block[..len];
        let mut i = 0;
        while i + PAT.len() < len {
            if &hay[i..i + PAT.len()] == PAT && hay[i + PAT.len()].is_ascii_whitespace() {
                return Ok(Some(block_start + i as u64));
            }
            i += 1;
        }
        let tail = keep.min(len);
        block.copy_within(len - tail..len, 0);
        block_start += (len - tail) as u64;
        carried = tail;
    }
}

pub fn read_spectra_parallel(
    path: &std::path::Path,
    ms_level_min: u32,
    ms_level_max: u32,
    threads: usize,
) -> std::io::Result<(Vec<Spectrum>, Vec<String>)> {
    read_spectra_in_pieces(
        path,
        ms_level_min,
        ms_level_max,
        threads,
        PARALLEL_PIECE_BYTES,
    )
}

fn spectrum_pieces(path: &std::path::Path, piece_bytes: u64) -> std::io::Result<Vec<(u64, u64)>> {
    let file = std::fs::File::open(path)?;
    let len = file.metadata()?.len();
    let n_cuts = len / piece_bytes;
    let mut starts: Vec<u64> = vec![0];
    for i in 1..=n_cuts {
        let from = (i * piece_bytes).max(*starts.last().unwrap() + 1);
        match next_spectrum_start(&file, from)? {
            Some(off) if off > *starts.last().unwrap() => starts.push(off),
            Some(_) => {}
            None => break,
        }
    }
    let n_pieces = starts.len();
    Ok((0..n_pieces)
        .map(|i| (starts[i], starts.get(i + 1).copied().unwrap_or(len)))
        .collect())
}

fn piece_bytes(
    path: &std::path::Path,
    pieces: &[(u64, u64)],
    i: usize,
) -> std::io::Result<Vec<u8>> {
    use std::io::{Read, Seek, SeekFrom};
    let (a, b) = pieces[i];
    let mut f = std::fs::File::open(path)?;
    f.seek(SeekFrom::Start(a))?;
    let mut bytes = Vec::with_capacity((b - a) as usize);
    f.take(b - a).read_to_end(&mut bytes)?;
    if i + 1 == pieces.len() {
        const END: &[u8] = b"</spectrumList";
        if let Some(p) = bytes.windows(END.len()).position(|w| w == END) {
            bytes.truncate(p);
        }
    }
    Ok(bytes)
}

fn map_pieces<T: Send>(
    n_pieces: usize,
    threads: usize,
    parse: impl Fn(usize) -> T + Sync,
) -> Vec<T> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let next = AtomicUsize::new(0);
    let mut results: Vec<(usize, T)> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..threads.clamp(1, n_pieces.max(1)))
            .map(|_| {
                scope.spawn(|| {
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n_pieces {
                            break;
                        }
                        out.push((i, parse(i)));
                    }
                    out
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| w.join().expect("mzML parse worker panicked"))
            .collect()
    });
    results.sort_by_key(|(i, _)| *i);
    results.into_iter().map(|(_, r)| r).collect()
}

fn read_spectra_in_pieces(
    path: &std::path::Path,
    ms_level_min: u32,
    ms_level_max: u32,
    threads: usize,
    piece_bytes_target: u64,
) -> std::io::Result<(Vec<Spectrum>, Vec<String>)> {
    let pieces = spectrum_pieces(path, piece_bytes_target)?;
    type PieceResult = std::io::Result<(Vec<Spectrum>, Vec<String>)>;
    let parse_piece = |i: usize| -> PieceResult {
        let bytes = piece_bytes(path, &pieces, i)?;
        let mut reader =
            MzMLReader::new(bytes.as_slice()).with_ms_level_range(ms_level_min, ms_level_max);
        let mut spectra = Vec::new();
        let mut errors = Vec::new();
        for r in reader.by_ref() {
            match r {
                Ok(s) => spectra.push(s),
                Err(e) => errors.push(e.to_string()),
            }
        }
        let kept = reader.skipped_errors().to_vec();
        let unlisted = reader.skipped_count() - kept.len();
        errors.extend(kept);
        errors.extend(std::iter::repeat_n(
            "malformed spectrum skipped".to_string(),
            unlisted,
        ));
        Ok((spectra, errors))
    };
    let mut spectra = Vec::new();
    let mut errors = Vec::new();
    for r in map_pieces(pieces.len(), threads, parse_piece) {
        let (s, e) = r?;
        spectra.extend(s);
        errors.extend(e);
    }
    Ok((spectra, errors))
}

#[derive(Debug, Default)]
pub struct Ms1LinkedRead {
    pub spectra: Vec<Spectrum>,
    pub ms2_to_ms1: Vec<Option<usize>>,
    pub run_ms1: Vec<Ms1Scan>,
    pub contexts: Vec<ScanContext>,
}

pub fn read_ms1_linked_parallel(
    path: &std::path::Path,
    threads: usize,
    capture_context: bool,
    capture_noise: bool,
) -> std::io::Result<Option<Ms1LinkedRead>> {
    read_ms1_linked_in_pieces(
        path,
        threads,
        capture_context,
        capture_noise,
        PARALLEL_PIECE_BYTES,
    )
}

/// Per MS2 scan of a piece: the MS1 scans before it and its FAIMS compensation voltage.
type Ms2Seen = (usize, Option<f32>);

fn read_ms1_linked_in_pieces(
    path: &std::path::Path,
    threads: usize,
    capture_context: bool,
    capture_noise: bool,
    piece_bytes_target: u64,
) -> std::io::Result<Option<Ms1LinkedRead>> {
    let pieces = spectrum_pieces(path, piece_bytes_target)?;
    let parse_piece = |i: usize| -> std::io::Result<Option<(Ms1LinkedRead, Vec<Ms2Seen>)>> {
        let bytes = piece_bytes(path, &pieces, i)?;
        let mut reader = MzMLReader::new(bytes.as_slice())
            .with_ms_level_range(2, 2)
            .with_run_ms1_capture(true)
            .with_scan_context_capture(capture_context)
            .with_noise_capture(capture_noise);
        let mut read = Ms1LinkedRead::default();
        // Per MS2 scan: the MS1 scans of the piece before it, and its compensation voltage.
        let mut seen: Vec<Ms2Seen> = Vec::new();
        loop {
            match reader.pump() {
                Ok(Some(spectrum)) => {
                    seen.push((reader.run_ms1.len(), reader.last_faims_cv));
                    read.spectra.push(spectrum);
                }
                Ok(None) => break,
                Err(_) => return Ok(None),
            }
        }
        read.run_ms1 = reader.take_run_ms1();
        read.contexts = reader.take_scan_contexts();
        Ok(Some((read, seen)))
    };
    let mut out = Ms1LinkedRead::default();
    let mut seen: Vec<Ms2Seen> = Vec::new();
    for r in map_pieces(pieces.len(), threads, parse_piece) {
        let Some((piece, piece_seen)) = r? else {
            return Ok(None);
        };
        let base = out.run_ms1.len();
        seen.extend(piece_seen.into_iter().map(|(n, cv)| (base + n, cv)));
        out.spectra.extend(piece.spectra);
        out.run_ms1.extend(piece.run_ms1);
        out.contexts.extend(piece.contexts);
    }
    // Link in acquisition order across the pieces.
    let mut tracker = crate::ms1_link::Ms1Tracker::default();
    let mut next = 0usize;
    out.ms2_to_ms1 = seen
        .into_iter()
        .map(|(n, cv)| {
            while next < n {
                tracker.push(next, out.run_ms1[next].faims_cv);
                next += 1;
            }
            tracker.link(cv)
        })
        .collect();
    Ok(Some(out))
}

pub fn filter_has_supplemental_collision(filter: &str) -> bool {
    let head = filter.split('[').next().unwrap_or("");
    let mut tokens = head.split_whitespace();
    if !tokens.any(|t| {
        let t = t.to_ascii_lowercase();
        t.len() > 2 && t.starts_with("ms") && t[2..].bytes().all(|b| b.is_ascii_digit())
    }) {
        return false;
    }
    let stages: Vec<&str> = tokens.filter(|t| t.contains('@')).collect();
    let Some(last) = stages.last() else {
        return false;
    };
    let precursor = |t: &str| t.split('@').next().unwrap_or("").to_string();
    let last_precursor = precursor(last);
    let methods: Vec<String> = stages
        .iter()
        .rev()
        .take_while(|t| precursor(t) == last_precursor)
        .flat_map(|t| t.split('@').skip(1))
        .map(|a| {
            a.chars()
                .take_while(char::is_ascii_alphabetic)
                .collect::<String>()
                .to_ascii_lowercase()
        })
        .collect();
    let electron = methods
        .iter()
        .any(|m| matches!(m.as_str(), "etd" | "ecd" | "netd"));
    let collisional = methods.iter().any(|m| matches!(m.as_str(), "hcd" | "cid"));
    electron && collisional
}

static ETHCD_AS_ETD: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

static ETHCD_POLICY_READ: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

pub fn init_ethcd_as_etd(as_etd: bool) {
    debug_assert!(
        !ETHCD_POLICY_READ.load(std::sync::atomic::Ordering::Relaxed) || !as_etd,
        "init_ethcd_as_etd(true) called after the policy was already read: every \
         spectrum read before this point was labelled under the default (EThcD)"
    );
    let _ = ETHCD_AS_ETD.set(as_etd);
}

#[inline]
pub fn ethcd_as_etd_policy() -> bool {
    ethcd_as_etd()
}

#[inline]
fn ethcd_as_etd() -> bool {
    ETHCD_POLICY_READ.store(true, std::sync::atomic::Ordering::Relaxed);
    *ETHCD_AS_ETD.get().unwrap_or(&false)
}

pub fn detect_instrument_type<R: BufRead>(reader: R) -> Option<InstrumentType> {
    let mut xml = Reader::from_reader(reader);
    xml.trim_text(true);

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum S {
        Outside,
        ReferenceableParamGroupList,
        ReferenceableParamGroup,
        InstrumentConfigurationList,
        InstrumentConfiguration,
        ComponentListAnalyzer,
        Run,
        Spectrum,
        Scan,
    }

    let mut state = S::Outside;
    let mut buf: Vec<u8> = Vec::with_capacity(4096);

    let mut ic_map: HashMap<String, InstrumentType> = HashMap::new();
    let mut current_ic_id: Option<String> = None;
    let mut current_ic_type: Option<InstrumentType> = None;
    let mut current_ic_analyzer: Option<InstrumentType> = None;

    let mut rpg_map: HashMap<String, InstrumentType> = HashMap::new();
    let mut current_rpg_id: Option<String> = None;
    let mut current_rpg_type: Option<InstrumentType> = None;

    let mut default_ic_ref: Option<String> = None;

    let mut ms2_counts: HashMap<InstrumentType, usize> = HashMap::new();
    let mut current_spec_is_ms2: Option<bool> = None;
    let mut current_spec_ic_ref: Option<String> = None;
    let mut ms2_seen: usize = 0;

    const MAX_PEEK: usize = 64;

    loop {
        buf.clear();
        let event = match xml.read_event_into(&mut buf) {
            Ok(e) => e,
            Err(_) => break,
        };
        match event {
            Event::Eof => break,

            Event::Start(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                match tag.as_slice() {
                    b"referenceableParamGroupList" if state == S::Outside => {
                        state = S::ReferenceableParamGroupList;
                    }
                    b"referenceableParamGroup" if state == S::ReferenceableParamGroupList => {
                        current_rpg_id = attr_str(e, b"id");
                        current_rpg_type = None;
                        state = S::ReferenceableParamGroup;
                    }
                    b"instrumentConfigurationList" if state == S::Outside => {
                        state = S::InstrumentConfigurationList;
                    }
                    b"instrumentConfiguration" if state == S::InstrumentConfigurationList => {
                        current_ic_id = attr_str(e, b"id");
                        current_ic_type = None;
                        current_ic_analyzer = None;
                        state = S::InstrumentConfiguration;
                    }
                    b"analyzer" if state == S::InstrumentConfiguration => {
                        state = S::ComponentListAnalyzer;
                    }
                    b"run" if state == S::Outside => {
                        default_ic_ref = attr_str(e, b"defaultInstrumentConfigurationRef");
                        state = S::Run;
                    }
                    b"spectrum" if state == S::Run => {
                        current_spec_is_ms2 = None;
                        current_spec_ic_ref = None;
                        state = S::Spectrum;
                    }
                    b"scan" if state == S::Spectrum => {
                        if let Some(r) = attr_str(e, b"instrumentConfigurationRef") {
                            current_spec_ic_ref = Some(r);
                        }
                        state = S::Scan;
                    }
                    _ => {}
                }
            }

            Event::Empty(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                if tag == b"scan" && state == S::Spectrum {
                    if let Some(r) = attr_str(e, b"instrumentConfigurationRef") {
                        current_spec_ic_ref = Some(r);
                    }
                }
                if tag == b"referenceableParamGroupRef" && state == S::InstrumentConfiguration {
                    if let Some(r) = attr_str(e, b"ref") {
                        if let Some(&t) = rpg_map.get(&r) {
                            current_ic_type = Some(t);
                        }
                    }
                }
                if tag == b"cvParam" {
                    let acc = attr_str(e, b"accession").unwrap_or_default();
                    match state {
                        S::ComponentListAnalyzer => {
                            let typ = match acc.as_str() {
                                CV_ANALYZER_ORBITRAP => Some(InstrumentType::QExactive),
                                CV_ANALYZER_FTICR => Some(InstrumentType::HighRes),
                                CV_ANALYZER_TOF => Some(InstrumentType::TOF),
                                CV_ANALYZER_ION_TRAP
                                | CV_ANALYZER_QUAD_ION_TRAP
                                | CV_ANALYZER_RADIAL_LIT
                                | CV_ANALYZER_LINEAR_ION_TRAP => Some(InstrumentType::LowRes),
                                _ => None,
                            };
                            if let Some(t) = typ {
                                if current_ic_type.is_none() {
                                    current_ic_type = Some(t);
                                }
                                current_ic_analyzer.get_or_insert(t);
                            }
                        }
                        S::InstrumentConfiguration => {
                            let cv_name = attr_str(e, b"name").unwrap_or_default();
                            if let Some(t) = instrument_model_type(acc.as_str(), &cv_name) {
                                current_ic_type = Some(t);
                            }
                        }
                        S::ReferenceableParamGroup => {
                            let cv_name = attr_str(e, b"name").unwrap_or_default();
                            if current_rpg_type.is_none() {
                                if let Some(t) = instrument_model_type(acc.as_str(), &cv_name) {
                                    current_rpg_type = Some(t);
                                }
                            }
                        }
                        S::Spectrum => {
                            if acc == CV_MS_LEVEL {
                                let val = attr_str(e, b"value").unwrap_or_default();
                                current_spec_is_ms2 = Some(val == "2");
                            }
                        }
                        _ => {}
                    }
                }
            }

            Event::End(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                match tag.as_slice() {
                    b"referenceableParamGroup" if state == S::ReferenceableParamGroup => {
                        if let (Some(id), Some(t)) =
                            (current_rpg_id.take(), current_rpg_type.take())
                        {
                            rpg_map.insert(id, t);
                        }
                        state = S::ReferenceableParamGroupList;
                    }
                    b"referenceableParamGroupList" if state == S::ReferenceableParamGroupList => {
                        state = S::Outside;
                    }
                    b"analyzer" if state == S::ComponentListAnalyzer => {
                        state = S::InstrumentConfiguration;
                    }
                    b"instrumentConfiguration" if state == S::InstrumentConfiguration => {
                        let ic_type = match current_ic_analyzer.take() {
                            Some(InstrumentType::LowRes) => Some(InstrumentType::LowRes),
                            _ => current_ic_type.take(),
                        };
                        if let (Some(id), Some(t)) = (current_ic_id.take(), ic_type) {
                            ic_map.insert(id, t);
                        }
                        state = S::InstrumentConfigurationList;
                    }
                    b"instrumentConfigurationList" if state == S::InstrumentConfigurationList => {
                        state = S::Outside;
                    }
                    b"scan" if state == S::Scan => {
                        state = S::Spectrum;
                    }
                    b"spectrum" if state == S::Spectrum => {
                        let is_ms2 = current_spec_is_ms2.unwrap_or(false);
                        if is_ms2 {
                            let ic_ref = current_spec_ic_ref
                                .clone()
                                .or_else(|| default_ic_ref.clone());
                            if let Some(r) = ic_ref {
                                if let Some(&t) = ic_map.get(&r) {
                                    *ms2_counts.entry(t).or_insert(0) += 1;
                                }
                            }
                            ms2_seen += 1;
                            if ms2_seen >= MAX_PEEK {
                                break;
                            }
                        }
                        current_spec_is_ms2 = None;
                        current_spec_ic_ref = None;
                        state = S::Run;
                    }
                    b"run" if state == S::Run => {
                        state = S::Outside;
                    }
                    _ => {}
                }
            }

            _ => {}
        }
    }

    if !ms2_counts.is_empty() {
        return ms2_counts
            .iter()
            .max_by_key(|(&k, &n)| (n, std::cmp::Reverse(k as u8)))
            .map(|(&t, _)| t);
    }

    if let Some(r) = default_ic_ref.as_ref() {
        if let Some(&t) = ic_map.get(r) {
            return Some(t);
        }
    }

    if ic_map.len() == 1 {
        return ic_map.into_values().next();
    }

    None
}

fn instrument_model_type(acc: &str, name: &str) -> Option<InstrumentType> {
    match acc {
        CV_MODEL_ORBITRAP_ASTRAL => Some(InstrumentType::OrbitrapAstral),
        CV_MODEL_Q_EXACTIVE
        | CV_MODEL_Q_EXACTIVE_HF
        | CV_MODEL_Q_EXACTIVE_HF_X
        | CV_MODEL_Q_EXACTIVE_PLUS
        | CV_MODEL_ORBITRAP_FUSION => Some(InstrumentType::QExactive),
        _ => {
            if name.to_ascii_lowercase().contains("astral") {
                Some(InstrumentType::OrbitrapAstral)
            } else {
                None
            }
        }
    }
}

fn attr_str(e: &quick_xml::events::BytesStart<'_>, name: &[u8]) -> Option<String> {
    e.attributes()
        .filter_map(|a| a.ok())
        .find(|a| a.key.local_name().as_ref() == name)
        .and_then(|a| {
            std::str::from_utf8(a.value.as_ref())
                .ok()
                .map(str::to_owned)
        })
}

fn extract_scan_from_id(id: &str) -> Option<i32> {
    ["scan=", "index=", "scanId="].iter().find_map(|key| {
        id.split_whitespace()
            .find_map(|tok| tok.strip_prefix(key)?.parse::<i32>().ok())
    })
}

fn scan_number(id: &str, index: Option<i32>) -> Option<i32> {
    extract_scan_from_id(id).or_else(|| index.map(|i| i.saturating_add(1)))
}

fn decode_binary_array(ctx: &BinaryArrayCtx) -> Result<Vec<f64>, MzMLParseError> {
    if let Some((accession, name)) = ctx.unsupported {
        return Err(MzMLParseError::UnsupportedCompression { accession, name });
    }
    let trimmed = ctx.b64_text.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }

    let raw = STANDARD.decode(trimmed)?;

    let bytes: Vec<u8> = if ctx.zlib {
        let mut decoder = ZlibDecoder::new(&raw[..]);
        let mut out = Vec::with_capacity(raw.len() * 2);
        std::io::Read::read_to_end(&mut decoder, &mut out).map_err(MzMLParseError::Zlib)?;
        out
    } else {
        raw
    };

    let elem_bytes = if ctx.precision_bits == 64 {
        8usize
    } else {
        4usize
    };
    if bytes.len() % elem_bytes != 0 {
        return Err(MzMLParseError::TruncatedBinary {
            decoded_len: bytes.len(),
            elem_bytes,
        });
    }

    let mut cur = std::io::Cursor::new(&bytes);
    let mut out: Vec<f64> = Vec::with_capacity(bytes.len() / elem_bytes);

    if ctx.precision_bits == 64 {
        while let Ok(v) = cur.read_f64::<LittleEndian>() {
            out.push(v);
        }
    } else {
        while let Ok(v) = cur.read_f32::<LittleEndian>() {
            out.push(v as f64);
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn collect_ok(xml: &str) -> Vec<Spectrum> {
        MzMLReader::new(Cursor::new(xml))
            .map(|r| r.expect("parse error"))
            .collect()
    }

    fn wrap_spectra(spectra: &str) -> String {
        format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectra}
    </spectrumList>
  </run>
</mzML>"#
        )
    }

    fn encode_f64_b64(vals: &[f64]) -> String {
        use byteorder::WriteBytesExt;
        let mut buf: Vec<u8> = Vec::with_capacity(vals.len() * 8);
        for &v in vals {
            buf.write_f64::<LittleEndian>(v).unwrap();
        }
        STANDARD.encode(&buf)
    }

    fn encode_f64_zlib_b64(vals: &[f64]) -> String {
        use byteorder::WriteBytesExt;
        use flate2::{write::ZlibEncoder, Compression};
        use std::io::Write;

        let mut raw: Vec<u8> = Vec::with_capacity(vals.len() * 8);
        for &v in vals {
            raw.write_f64::<LittleEndian>(v).unwrap();
        }
        let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
        enc.write_all(&raw).unwrap();
        STANDARD.encode(enc.finish().unwrap())
    }

    fn bda_block(cv_array: &str, compression_cv: &str, b64: &str) -> String {
        format!(
            r#"<binaryDataArray>
              <cvParam accession="MS:1000523" name="64-bit float" value=""/>
              <cvParam accession="{compression_cv}" name="compression" value=""/>
              <cvParam accession="{cv_array}" name="" value=""/>
              <binary>{b64}</binary>
            </binaryDataArray>"#
        )
    }

    fn bda_plain(cv_array: &str, b64: &str) -> String {
        bda_block(cv_array, "MS:1000576", b64)
    }

    fn bda_zlib(cv_array: &str, b64: &str) -> String {
        bda_block(cv_array, "MS:1000574", b64)
    }

    fn ms2_spectrum_xml(
        id: &str,
        mz_bda: &str,
        int_bda: &str,
        precursor_mz: f64,
        charge: Option<i32>,
    ) -> String {
        let charge_param = match charge {
            Some(z) => {
                format!(r#"<cvParam accession="MS:1000041" name="charge state" value="{z}"/>"#)
            }
            None => String::new(),
        };
        format!(
            r#"<spectrum index="0" id="{id}" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <cvParam accession="MS:1000016" name="scan start time" value="1.5"
                           unitAccession="UO:0000031" unitName="minute"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="{precursor_mz}"/>
                      {charge_param}
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz_bda}
                {int_bda}
              </binaryDataArrayList>
            </spectrum>"#
        )
    }

    #[test]
    fn parses_minimal_mzml_with_one_ms2_spectrum() {
        let mz_b64 = encode_f64_b64(&[100.0, 200.0]);
        let int_b64 = encode_f64_b64(&[1000.0, 500.0]);

        let spec = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(2),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));

        assert_eq!(spectra.len(), 1, "expected exactly one MS2 spectrum");
        assert_eq!(spectra[0].peaks.len(), 2, "expected two peaks");
    }

    #[test]
    fn decodes_zlib_compressed_peaks() {
        let mz_vals = [150.0_f64, 300.0, 450.0];
        let int_vals = [2000.0_f64, 1000.0, 500.0];

        let spec = format!(
            r#"<spectrum index="0" id="scan=7" defaultArrayLength="3">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.0"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_zlib("MS:1000514", &encode_f64_zlib_b64(&mz_vals)),
            int = bda_zlib("MS:1000515", &encode_f64_zlib_b64(&int_vals)),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));

        assert_eq!(spectra.len(), 1);
        let peaks = &spectra[0].peaks;
        assert_eq!(peaks.len(), 3);
        assert!((peaks[0].0 - 150.0).abs() < 1e-6, "first m/z");
        assert!((peaks[1].0 - 300.0).abs() < 1e-6, "second m/z");
        assert!((peaks[2].0 - 450.0).abs() < 1e-6, "third m/z");
    }

    #[test]
    fn decodes_uncompressed_64bit_peaks() {
        let mz_b64 = encode_f64_b64(&[200.0, 400.0]);
        let int_b64 = encode_f64_b64(&[5000.0, 2500.0]);

        let spec = ms2_spectrum_xml(
            "scan=3",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            600.0,
            None,
        );
        let spectra = collect_ok(&wrap_spectra(&spec));

        assert_eq!(spectra.len(), 1);
        let peaks = &spectra[0].peaks;
        assert_eq!(peaks.len(), 2);
        assert!((peaks[0].0 - 200.0).abs() < 1e-6);
        assert!((peaks[1].0 - 400.0).abs() < 1e-6);
        assert!((peaks[0].1 - 5000.0_f32).abs() < 1.0);
        assert!((peaks[1].1 - 2500.0_f32).abs() < 1.0);
    }

    #[test]
    fn filters_out_ms1_spectra() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[100.0]);

        let ms1 = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );

        let ms2_mz_b64 = encode_f64_b64(&[200.0, 300.0]);
        let ms2_int_b64 = encode_f64_b64(&[800.0, 400.0]);
        let ms2 = ms2_spectrum_xml(
            "scan=2",
            &bda_plain("MS:1000514", &ms2_mz_b64),
            &bda_plain("MS:1000515", &ms2_int_b64),
            500.0,
            Some(2),
        );

        let xml = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="2" defaultDataProcessingRef="dp">
      {ms1}
      {ms2}
    </spectrumList>
  </run>
</mzML>"#
        );

        let spectra = collect_ok(&xml);
        assert_eq!(spectra.len(), 1, "only the MS2 should be emitted");
        assert_eq!(spectra[0].scan, Some(2));
    }

    #[test]
    fn extracts_scan_number_from_id_attr() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);

        let spec = ms2_spectrum_xml(
            "controllerType=0 controllerNumber=1 scan=1234",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.0,
            None,
        );
        let spectra = collect_ok(&wrap_spectra(&spec));

        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].scan, Some(1234));
    }

    #[test]
    fn extracts_precursor_mz_and_charge() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);

        let spec = ms2_spectrum_xml(
            "scan=10",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(2),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));

        assert_eq!(spectra.len(), 1);
        assert!((spectra[0].precursor_mz - 500.5).abs() < 1e-6);
        assert_eq!(spectra[0].precursor_charge, Some(2));
    }

    #[test]
    fn first_selected_ion_gives_mz_and_charge_together() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let ions = |second: &str| {
            format!(
                r#"<spectrum index="0" id="scan=7" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <cvParam accession="MS:1000016" name="scan start time" value="1.5"
                           unitAccession="UO:0000031" unitName="minute"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="2">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.25"/>
                      <cvParam accession="MS:1000041" name="charge state" value="2"/>
                      <cvParam accession="MS:1000042" name="peak intensity" value="1000"/>
                    </selectedIon>
                    <selectedIon>{second}</selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {}
                {}
              </binaryDataArrayList>
            </spectrum>"#,
                bda_plain("MS:1000514", &mz_b64),
                bda_plain("MS:1000515", &int_b64),
            )
        };
        for second in [
            r#"<cvParam accession="MS:1000744" name="selected ion m/z" value="733.70"/>"#,
            r#"<cvParam accession="MS:1000744" name="selected ion m/z" value="733.70"/>
               <cvParam accession="MS:1000041" name="charge state" value="3"/>
               <cvParam accession="MS:1000042" name="peak intensity" value="50"/>"#,
            r#"<cvParam accession="MS:1000041" name="charge state" value="4"/>"#,
        ] {
            let spectra = collect_ok(&wrap_spectra(&ions(second)));
            assert_eq!(spectra.len(), 1);
            assert!((spectra[0].precursor_mz - 500.25).abs() < 1e-9, "{second}");
            assert_eq!(spectra[0].precursor_charge, Some(2), "{second}");
            assert_eq!(spectra[0].precursor_intensity, Some(1000.0), "{second}");
        }
    }

    #[test]
    fn out_of_range_charge_state_is_treated_as_absent() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        for z in [300, -2] {
            let spec = ms2_spectrum_xml(
                "scan=10",
                &bda_plain("MS:1000514", &mz_b64),
                &bda_plain("MS:1000515", &int_b64),
                500.5,
                Some(z),
            );
            let spectra = collect_ok(&wrap_spectra(&spec));
            assert_eq!(spectra[0].precursor_charge, None, "charge {z}");
        }
    }

    #[test]
    fn charge_state_zero_is_treated_as_absent() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let spec = ms2_spectrum_xml(
            "scan=10",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(0),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].precursor_charge, None);
    }

    #[test]
    fn non_finite_and_nonpositive_peaks_are_filtered() {
        let mz_b64 = encode_f64_b64(&[f64::NAN, 100.0, -50.0, f64::INFINITY, 200.0]);
        let int_b64 = encode_f64_b64(&[10.0, 1000.0, 10.0, 10.0, 500.0]);
        let spec = ms2_spectrum_xml(
            "scan=11",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(2),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));
        assert_eq!(spectra.len(), 1);
        let mzs: Vec<f64> = spectra[0].peaks.iter().map(|&(mz, _)| mz).collect();
        assert_eq!(
            mzs,
            vec![100.0, 200.0],
            "only finite positive-m/z peaks survive"
        );
    }

    #[test]
    fn nonpositive_precursor_mz_is_skipped() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let spec = ms2_spectrum_xml(
            "scan=12",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            0.0,
            Some(2),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));
        assert!(
            spectra.is_empty(),
            "spectrum with precursor m/z 0 must be skipped"
        );
    }

    #[test]
    fn peaks_sorted_ascending_by_mz() {
        let mz_b64 = encode_f64_b64(&[300.0, 100.0, 200.0]);
        let int_b64 = encode_f64_b64(&[3.0, 1.0, 2.0]);

        let spec = format!(
            r#"<spectrum index="0" id="scan=5" defaultArrayLength="3">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="600.0"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));

        assert_eq!(spectra.len(), 1);
        let mzs: Vec<f64> = spectra[0].peaks.iter().map(|p| p.0).collect();
        assert_eq!(mzs, vec![100.0, 200.0, 300.0]);
    }

    #[test]
    fn parses_real_test_fixture() {
        let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../test-fixtures/tiny.pwiz.mzML");

        if !fixture.exists() {
            eprintln!("SKIP: fixture not found at {}", fixture.display());
            return;
        }

        let file = std::fs::File::open(&fixture).expect("failed to open tiny.pwiz.mzML");
        let spectra: Vec<Spectrum> = MzMLReader::new(std::io::BufReader::new(file))
            .map(|r| r.expect("parse error"))
            .collect();

        assert!(!spectra.is_empty(), "expected at least one MS2 spectrum");
        let s = &spectra[0];
        assert!(!s.peaks.is_empty(), "MS2 spectrum should have peaks");
        assert!(s.precursor_mz > 0.0, "precursor m/z should be positive");
    }

    #[test]
    fn minute_unit_by_accession_legacy_term_or_name() {
        assert!(is_minute_unit("UO:0000031", "minute"));
        assert!(is_minute_unit("MS:1000038", ""));
        assert!(is_minute_unit("", "minute"));
        assert!(is_minute_unit("", "Min"));
        assert!(!is_minute_unit("UO:0000010", "second"));
        assert!(!is_minute_unit("", "second"));
        assert!(!is_minute_unit("", ""));
    }

    #[test]
    fn numpress_array_is_an_error_not_garbage_peaks() {
        let b64 = encode_f64_b64(&[100.0]);
        let xml = ms2_spectrum_xml(
            "scan=1",
            &bda_block("MS:1000514", "MS:1002312", &b64),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1000.0])),
            500.5,
            Some(2),
        );
        let mut reader = MzMLReader::new(Cursor::new(wrap_spectra(&xml)));
        assert_eq!(reader.by_ref().filter_map(Result::ok).count(), 0);
        assert_eq!(reader.skipped_count(), 1);
        let err = &reader.skipped_errors()[0];
        assert!(
            err.contains("MS:1002312") && err.contains("Numpress"),
            "{err}"
        );
    }

    #[test]
    fn truncation_zlib_decodes_like_zlib() {
        let mz: Vec<u8> = [100.0f64, 200.0]
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        let mut enc = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut enc, &mz).unwrap();
        let b64 = STANDARD.encode(enc.finish().unwrap());
        let xml = ms2_spectrum_xml(
            "scan=1",
            &bda_block("MS:1000514", "MS:1003088", &b64),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0, 2.0])),
            500.5,
            Some(2),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra[0].peaks, vec![(100.0, 1.0), (200.0, 2.0)]);
    }

    #[test]
    fn extract_scan_plain() {
        assert_eq!(extract_scan_from_id("scan=1234"), Some(1234));
    }

    #[test]
    fn extract_scan_pwiz_format() {
        assert_eq!(
            extract_scan_from_id("controllerType=0 controllerNumber=1 scan=42"),
            Some(42)
        );
    }

    #[test]
    fn thermo_trailer_monoisotopic_overrides_selected_ion_mz() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <cvParam accession="MS:1000016" name="scan start time"
                           value="1.5" unitAccession="UO:0000031"
                           unitName="minute"/>
                  <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:"
                             type="xsd:float" value="625.0037"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="625.338134765625"/>
                      <cvParam accession="MS:1000041" name="charge state" value="3"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {}
                {}
              </binaryDataArrayList>
            </spectrum>"#,
            bda_plain("MS:1000514", &mz_b64),
            bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra.len(), 1);
        assert!(
            (spectra[0].precursor_mz - 625.0037).abs() < 1e-6,
            "expected Thermo trailer monoisotopic m/z (625.0037), got {}",
            spectra[0].precursor_mz
        );
    }

    #[test]
    fn precursor_mz_falls_back_to_selected_ion_without_trailer() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let spec = ms2_spectrum_xml(
            "scan=42",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(2),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));
        assert_eq!(spectra.len(), 1);
        assert!((spectra[0].precursor_mz - 500.5).abs() < 1e-6);
    }

    #[test]
    fn zero_thermo_trailer_does_not_override() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:"
                             type="xsd:float" value="0"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="700.25"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {}
                {}
              </binaryDataArrayList>
            </spectrum>"#,
            bda_plain("MS:1000514", &mz_b64),
            bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra.len(), 1);
        assert!(
            (spectra[0].precursor_mz - 700.25).abs() < 1e-6,
            "zero-trailer must fall back to selectedIon m/z; got {}",
            spectra[0].precursor_mz
        );
    }

    #[test]
    fn extract_scan_missing() {
        assert_eq!(extract_scan_from_id("spectrum=1"), None);
        assert_eq!(extract_scan_from_id("index=101"), Some(101));
        assert_eq!(extract_scan_from_id("scanId=7"), Some(7));
        assert_eq!(
            extract_scan_from_id("controllerType=0 controllerNumber=1 scan=42 index=9"),
            Some(42)
        );
        assert_eq!(
            scan_number("sample=1 period=1 cycle=4 experiment=2", Some(16)),
            Some(17)
        );
        assert_eq!(scan_number("index=101", Some(0)), Some(101));
    }

    fn spectrum_xml_with_activation(activation_cv: Option<&str>) -> String {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let act_block = match activation_cv {
            Some(cv) => format!(
                r#"<activation>
                     <cvParam accession="{cv}" name="" value=""/>
                   </activation>"#
            ),
            None => String::new(),
        };
        format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="500.5"/>
                      <cvParam accession="MS:1000041" name="charge state" value="2"/>
                    </selectedIon>
                  </selectedIonList>
                  {act_block}
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        )
    }

    #[test]
    fn parses_cid_activation() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
            "MS:1000133",
        ))));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::CID));
    }

    #[test]
    fn parses_hcd_activation() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
            "MS:1000422",
        ))));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::HCD));
    }

    #[test]
    fn parses_etd_activation() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
            "MS:1000598",
        ))));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
    }

    #[test]
    fn parses_ecd_as_etd() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
            "MS:1000250",
        ))));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
    }

    #[test]
    fn missing_activation_block_yields_none() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(None)));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, None);
    }

    fn spectrum_xml_with_two_activations(cv1: &str, cv2: &str) -> String {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                      <cvParam accession="MS:1000041" name="charge state" value="2"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="{cv1}" name="" value=""/>
                    <cvParam accession="{cv2}" name="" value=""/>
                  </activation>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        )
    }

    #[test]
    fn ethcd_etd_plus_hcd_is_labelled_ethcd() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(
            "MS:1000598",
            "MS:1000422",
        )));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
    }

    #[test]
    fn ethcd_detection_is_order_independent() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(
            "MS:1000422",
            "MS:1000598",
        )));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
    }

    #[test]
    fn supplemental_activation_terms_mark_ethcd_and_etcid() {
        for supp in ["MS:1002678", "MS:1002679"] {
            for (a, b) in [("MS:1000598", supp), (supp, "MS:1000598")] {
                let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(a, b)));
                assert_eq!(
                    spectra[0].activation_method,
                    Some(ActivationMethod::EThcD),
                    "{a} + {b}"
                );
            }
        }
    }

    #[test]
    fn supplemental_collision_read_from_the_filter_chain() {
        for f in [
            "FTMS + p NSI d Full ms2 784.37@etd50.00@hcd25.00 [120.00-2000.00]",
            "ITMS + c NSI d Full ms2 652.10@etd100.00@cid35.00 [50.00-1315.00]",
            "FTMS + c NSI Full ms3 900.00@cid35.00 450.20@etd50.00@hcd20.00 [100.00-1000.00]",
            "FTMS + p NSI d Full ms2 391.1936@etd55.49 391.1936@hcd20.00 [100.0000-1184.0000]",
            "FTMS + c NSI d sa Full ms2 361.9050@etd143.26@hcd28.00 [100.0000-734.0000]",
        ] {
            assert!(filter_has_supplemental_collision(f), "{f}");
        }
        for f in [
            "FTMS + p NSI d Full ms2 784.37@etd50.00 [120.00-2000.00]",
            "FTMS + p NSI d Full ms2 784.37@hcd28.00 [120.00-2000.00]",
            "FTMS + c NSI Full ms3 900.00@cid35.00 450.20@etd50.00 [100.00-1000.00]",
            "FTMS + c NSI Full ms3 900.00@hcd35.00 450.20@etd50.00 [100.00-1000.00]",
            "FTMS + p NSI Full ms [350.00-1800.00]",
            "",
        ] {
            assert!(!filter_has_supplemental_collision(f), "{f}");
        }
    }

    #[test]
    fn etd_term_with_supplemental_filter_string_is_labelled_ethcd() {
        let xml = spectrum_xml_with_activation(Some("MS:1000598")).replace(
            "<scanList count=\"1\"><scan/></scanList>",
            "<scanList count=\"1\"><scan><cvParam accession=\"MS:1000512\" name=\"filter string\" \
             value=\"FTMS + p NSI d Full ms2 500.50@etd55.49 500.50@hcd20.00 [100.00-1000.00]\"/>\
             </scan></scanList>",
        );
        assert!(
            xml.contains("MS:1000512"),
            "fixture must carry the filter string"
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
        let xml = spectrum_xml_with_activation(Some("MS:1000598")).replace(
            "<scanList count=\"1\"><scan/></scanList>",
            "<scanList count=\"1\"><scan><cvParam accession=\"MS:1000512\" name=\"filter string\" \
             value=\"FTMS + p NSI d Full ms2 500.50@etd55.49 [100.00-1000.00]\"/></scan></scanList>",
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
    }

    #[test]
    fn combined_and_child_activation_terms() {
        for (cv, want) in [
            ("MS:1002631", ActivationMethod::EThcD),
            ("MS:1003182", ActivationMethod::EThcD),
            ("MS:1002472", ActivationMethod::CID),
            ("MS:1002481", ActivationMethod::HCD),
            ("MS:1003246", ActivationMethod::UVPD),
        ] {
            let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(cv))));
            assert_eq!(spectra[0].activation_method, Some(want), "{cv}");
        }
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
            "MS:1002678",
        ))));
        assert_eq!(spectra[0].activation_method, None);
    }

    #[test]
    fn etcid_etd_plus_cid_is_labelled_ethcd() {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(
            "MS:1000598",
            "MS:1000133",
        )));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
    }

    #[test]
    fn centroid_profile_collapses_one_peak_to_intensity_weighted_centroid() {
        let prof = [
            (100.0, 1.0f32),
            (100.1, 4.0),
            (100.2, 9.0),
            (100.3, 4.0),
            (100.4, 1.0),
        ];
        let c = centroid_profile(&prof);
        assert_eq!(c.len(), 1, "one peak ⇒ one centroid");
        assert!(
            (c[0].0 - 100.2).abs() < 1e-9,
            "centroid m/z {} != 100.2",
            c[0].0
        );
        assert_eq!(c[0].1, 9.0, "centroid intensity is the apex");
    }

    #[test]
    fn centroid_profile_separates_two_peaks() {
        let prof = [
            (200.0, 1.0f32),
            (200.1, 5.0),
            (200.2, 1.0),
            (200.3, 0.0),
            (200.4, 2.0),
            (200.5, 8.0),
            (200.6, 2.0),
        ];
        let c = centroid_profile(&prof);
        assert_eq!(c.len(), 2, "two peaks ⇒ two centroids");
        assert!((c[0].0 - 200.1).abs() < 1e-9);
        assert!((c[1].0 - 200.5).abs() < 1e-9);
    }

    fn profile_spectrum_xml(profile: bool) -> String {
        let mz_b64 = encode_f64_b64(&[150.0, 150.1, 150.2]);
        let int_b64 = encode_f64_b64(&[2.0, 9.0, 2.0]);
        let rep = if profile {
            r#"<cvParam accession="MS:1000128" name="profile spectrum" value=""/>"#
        } else {
            r#"<cvParam accession="MS:1000127" name="centroid spectrum" value=""/>"#
        };
        format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="3">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              {rep}
              <scanList count="1"><scan/></scanList>
              <precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>
                <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                <cvParam accession="MS:1000041" name="charge state" value="2"/>
              </selectedIon></selectedIonList></precursor></precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        )
    }

    #[test]
    fn profile_mode_spectrum_is_centroided_on_parse() {
        let prof = collect_ok(&wrap_spectra(&profile_spectrum_xml(true)));
        assert_eq!(prof.len(), 1);
        assert_eq!(prof[0].peaks.len(), 1, "profile 3-point peak ⇒ 1 centroid");
        assert!((prof[0].peaks[0].0 - 150.1).abs() < 1e-9);

        let cent = collect_ok(&wrap_spectra(&profile_spectrum_xml(false)));
        assert_eq!(
            cent[0].peaks.len(),
            3,
            "centroided spectrum is not re-centroided"
        );
    }

    #[test]
    fn multiple_activations_first_wins() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="2">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="MS:1000133" name="CID" value=""/>
                  </activation>
                </precursor>
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="350.0"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="MS:1000422" name="HCD" value=""/>
                  </activation>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );

        let wrapped = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {xml}
    </spectrumList>
  </run>
</mzML>"#
        );
        let spectra: Vec<Spectrum> = MzMLReader::new(Cursor::new(wrapped))
            .with_ms_level_range(2, 3)
            .map(|r| r.expect("parse error"))
            .collect();
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::CID));
    }

    fn activation_of(activation_xml: &str) -> Option<ActivationMethod> {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  </selectedIon></selectedIonList>
                  <activation>{activation_xml}</activation>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let wrapped = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run><spectrumList count="1" defaultDataProcessingRef="dp">
      {xml}
  </spectrumList></run>
</mzML>"#
        );
        let spectra: Vec<Spectrum> = MzMLReader::new(Cursor::new(wrapped))
            .map(|r| r.expect("parse error"))
            .collect();
        assert_eq!(spectra.len(), 1);
        spectra[0].activation_method
    }

    #[test]
    fn activation_user_param_sets_activation() {
        assert_eq!(
            activation_of(r#"<userParam name="CID" value=""/>"#),
            Some(ActivationMethod::CID)
        );
        assert_eq!(
            activation_of(r#"<userParam name="beam-type collision-induced dissociation"/>"#),
            Some(ActivationMethod::HCD)
        );
        assert_eq!(
            activation_of(r#"<userParam name="ETD"/>"#),
            Some(ActivationMethod::ETD)
        );
        assert_eq!(
            activation_of(r#"<userParam name="collision energy" value="35"/>"#),
            None
        );
        assert_eq!(
            activation_of(
                r#"<cvParam accession="MS:1000422" name="HCD" value=""/><userParam name="CID"/>"#
            ),
            Some(ActivationMethod::HCD)
        );
    }

    #[test]
    fn separate_activation_blocks_do_not_falsely_trigger_ethcd() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="2">
                <precursor>
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  </selectedIon></selectedIonList>
                  <activation><cvParam accession="MS:1000133" name="CID" value=""/></activation>
                </precursor>
                <precursor>
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="350.0"/>
                  </selectedIon></selectedIonList>
                  <activation><cvParam accession="MS:1000598" name="ETD" value=""/></activation>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let wrapped = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run><spectrumList count="1" defaultDataProcessingRef="dp">
      {xml}
  </spectrumList></run>
</mzML>"#
        );
        let spectra: Vec<Spectrum> = MzMLReader::new(Cursor::new(wrapped))
            .with_ms_level_range(2, 3)
            .map(|r| r.expect("parse error"))
            .collect();
        assert_eq!(spectra.len(), 1);
        assert_ne!(
            spectra[0].activation_method,
            Some(ActivationMethod::EThcD),
            "separate CID + ETD blocks must not be mis-detected as EThcD"
        );
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
    }

    #[test]
    fn etd_plus_supplemental_collisional_is_labelled_ethcd() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="MS:1000133" name="CID" value=""/>
                    <cvParam accession="MS:1000598" name="ETD" value=""/>
                  </activation>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
    }

    #[test]
    fn parses_isolation_window_offsets() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <isolationWindow>
                    <cvParam accession="MS:1000827" name="isolation window target m/z"
                             value="500.5"/>
                    <cvParam accession="MS:1000828" name="isolation window lower offset"
                             value="1.5"/>
                    <cvParam accession="MS:1000829" name="isolation window upper offset"
                             value="1.5"/>
                  </isolationWindow>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].isolation_lower_offset, Some(1.5));
        assert_eq!(spectra[0].isolation_upper_offset, Some(1.5));
        assert_eq!(spectra[0].isolation_target_mz, Some(500.5));
        assert!((spectra[0].precursor_mz - 500.5).abs() < 1e-6);
    }

    #[test]
    fn product_isolation_window_does_not_overwrite_precursor_window() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let xml = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <isolationWindow>
                    <cvParam accession="MS:1000827" name="isolation window target m/z" value="500.5"/>
                    <cvParam accession="MS:1000828" name="isolation window lower offset" value="0.7"/>
                    <cvParam accession="MS:1000829" name="isolation window upper offset" value="0.7"/>
                  </isolationWindow>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <productList count="1">
                <product>
                  <isolationWindow>
                    <cvParam accession="MS:1000827" name="isolation window target m/z" value="800.0"/>
                    <cvParam accession="MS:1000828" name="isolation window lower offset" value="400.0"/>
                    <cvParam accession="MS:1000829" name="isolation window upper offset" value="400.0"/>
                  </isolationWindow>
                </product>
              </productList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].isolation_lower_offset, Some(0.7));
        assert_eq!(spectra[0].isolation_upper_offset, Some(0.7));
        assert_eq!(spectra[0].isolation_target_mz, Some(500.5));
        assert_eq!(spectra[0].peaks.len(), 1);
    }

    #[test]
    fn missing_isolation_window_yields_none() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let spec = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(2),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].isolation_lower_offset, None);
        assert_eq!(spectra[0].isolation_upper_offset, None);
        assert_eq!(spectra[0].isolation_target_mz, None);
    }

    fn wrap_with_instrument_configs(
        instrument_configs: &str,
        default_ic_ref: &str,
        spectra_xml: &str,
    ) -> String {
        format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <instrumentConfigurationList count="1">
    {instrument_configs}
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="{default_ic_ref}">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectra_xml}
    </spectrumList>
  </run>
</mzML>"#
        )
    }

    fn ic_block(id: &str, analyzer_cv: &str) -> String {
        format!(
            r#"<instrumentConfiguration id="{id}">
              <componentList count="3">
                <source order="1">
                  <cvParam accession="MS:1000398" name="nanoelectrospray" value=""/>
                </source>
                <analyzer order="2">
                  <cvParam accession="{analyzer_cv}" name="" value=""/>
                </analyzer>
                <detector order="3">
                  <cvParam accession="MS:1000624" name="inductive detector" value=""/>
                </detector>
              </componentList>
            </instrumentConfiguration>"#
        )
    }

    fn ms2_spectrum_with_ic_ref(ic_ref: &str) -> String {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan instrumentConfigurationRef="{ic_ref}"/>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        )
    }

    #[test]
    fn detect_instrument_orbitrap_analyzer_to_qexactive() {
        let xml = wrap_with_instrument_configs(
            &ic_block("IC1", "MS:1000484"),
            "IC1",
            &ms2_spectrum_with_ic_ref("IC1"),
        );
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::QExactive));
    }

    #[test]
    fn detect_instrument_ion_trap_analyzer_to_lowres() {
        let xml = wrap_with_instrument_configs(
            &ic_block("IC1", "MS:1000291"),
            "IC1",
            &ms2_spectrum_with_ic_ref("IC1"),
        );
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::LowRes));
    }

    #[test]
    fn detect_instrument_quad_ion_trap_to_lowres() {
        let xml = wrap_with_instrument_configs(
            &ic_block("IC1", "MS:1000082"),
            "IC1",
            &ms2_spectrum_with_ic_ref("IC1"),
        );
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::LowRes));
    }

    #[test]
    fn detect_instrument_fticr_to_highres() {
        let xml = wrap_with_instrument_configs(
            &ic_block("IC1", "MS:1000079"),
            "IC1",
            &ms2_spectrum_with_ic_ref("IC1"),
        );
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::HighRes));
    }

    #[test]
    fn detect_instrument_tof_analyzer() {
        let xml = wrap_with_instrument_configs(
            &ic_block("IC1", "MS:1000084"),
            "IC1",
            &ms2_spectrum_with_ic_ref("IC1"),
        );
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::TOF));
    }

    #[test]
    fn detect_instrument_ms2_referenced_ic_wins_pxd001819_pattern() {
        let ics = format!(
            "{}\n{}",
            ic_block("IC1", "MS:1000484"),
            ic_block("IC2", "MS:1000264"),
        );
        let xml = wrap_with_instrument_configs(&ics, "IC1", &ms2_spectrum_with_ic_ref("IC2"));
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::LowRes));
    }

    #[test]
    fn detect_instrument_falls_back_to_default_ic_when_scan_lacks_ref() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let spec = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let xml = wrap_with_instrument_configs(&ic_block("IC1", "MS:1000484"), "IC1", &spec);
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::QExactive));
    }

    #[test]
    fn detect_instrument_returns_none_when_no_ic_info() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        let spec = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            None,
        );
        let xml = wrap_spectra(&spec);
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, None);
    }

    #[test]
    fn ms1_capture_links_two_ms2_to_preceding_ms1() {
        let ms1_mz_b64 = encode_f64_b64(&[111.0, 222.0]);
        let ms1_int_b64 = encode_f64_b64(&[10.0, 20.0]);
        let ms1 = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &ms1_mz_b64),
            int = bda_plain("MS:1000515", &ms1_int_b64),
        );

        let ms2a = ms2_spectrum_xml(
            "scan=2",
            &bda_plain("MS:1000514", &encode_f64_b64(&[300.0, 400.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0, 2.0])),
            500.0,
            Some(2),
        );
        let ms2b = ms2_spectrum_xml(
            "scan=3",
            &bda_plain("MS:1000514", &encode_f64_b64(&[500.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[3.0])),
            600.0,
            Some(2),
        );

        let xml = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="3" defaultDataProcessingRef="dp">
      {ms1}
      {ms2a}
      {ms2b}
    </spectrumList>
  </run>
</mzML>"#
        );

        let (spectra, link) = MzMLReader::new(Cursor::new(xml))
            .with_ms1_capture(true)
            .read_with_ms1()
            .expect("read_with_ms1 failed");

        assert_eq!(spectra.len(), 2, "expected exactly two MS2 spectra");
        assert_eq!(spectra[0].scan, Some(2));
        assert_eq!(spectra[1].scan, Some(3));

        assert_eq!(link.ms1_peaks.len(), 1, "expected one captured MS1");
        let ms1_peaks = &link.ms1_peaks[0];
        assert_eq!(ms1_peaks.len(), 2);
        assert!((ms1_peaks[0].0 - 111.0).abs() < 1e-6);
        assert!((ms1_peaks[1].0 - 222.0).abs() < 1e-6);
        assert!((ms1_peaks[0].1 - 10.0_f32).abs() < 1e-3);
        assert!((ms1_peaks[1].1 - 20.0_f32).abs() < 1e-3);

        assert_eq!(link.ms2_to_ms1, vec![Some(0), Some(0)]);
    }

    #[test]
    fn run_ms1_capture_keeps_rt_and_peaks_without_emitting_ms1() {
        let ms1 = |id: &str, rt_min: f64, mz: f64| {
            format!(
                r#"<spectrum index="0" id="{id}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="{rt_min}"
                         unitAccession="UO:0000031" unitName="minute"/>
              </scan></scanList>
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
                mzb = bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
                intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
            )
        };
        let ms2 = ms2_spectrum_xml(
            "scan=2",
            &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            350.0,
            Some(2),
        );
        let xml = wrap_spectra(&format!(
            "{}{}{}",
            ms1("scan=1", 1.0, 400.0),
            ms2,
            ms1("scan=3", 1.5, 401.0)
        ));
        let mut reader = MzMLReader::new(Cursor::new(xml)).with_run_ms1_capture(true);
        let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
        assert_eq!(spectra.len(), 1, "only the MS2 is emitted");
        let ms1s = reader.take_run_ms1();
        assert_eq!(ms1s.len(), 2);
        assert_eq!(ms1s[0].rt, 60.0);
        assert_eq!(ms1s[1].rt, 90.0);
        assert_eq!(ms1s[1].peaks, vec![(401.0, 7.0)]);
        assert!(reader.take_run_ms1().is_empty(), "take drains");
    }

    #[test]
    fn scan_context_records_survey_reference_scan_window_and_faims_voltage() {
        let spectrum = |id: &str, level: u32, extra: &str| {
            format!(
                r#"<spectrum index="0" id="{id}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="{level}"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="1.0"
                         unitAccession="UO:0000031" unitName="minute"/>
                <cvParam accession="MS:1001581" name="FAIMS compensation voltage" value="-45"/>
                <scanWindowList count="1"><scanWindow>
                  <cvParam accession="MS:1000501" name="scan window lower limit" value="148.0"/>
                  <cvParam accession="MS:1000500" name="scan window upper limit" value="1866.0"/>
                </scanWindow></scanWindowList>
              </scan></scanList>
              {extra}
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
                mzb = bda_plain("MS:1000514", &encode_f64_b64(&[400.0])),
                intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
            )
        };
        let precursor = r#"<precursorList count="1"><precursor spectrumRef="scan=1">
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    <cvParam accession="MS:1000041" name="charge state" value="2"/>
                  </selectedIon></selectedIonList>
                </precursor></precursorList>"#;
        let xml = wrap_spectra(&format!(
            "{}{}",
            spectrum("scan=1", 1, ""),
            spectrum("scan=2", 2, precursor)
        ));
        let mut reader = MzMLReader::new(Cursor::new(xml.clone()))
            .with_run_ms1_capture(true)
            .with_scan_context_capture(true);
        let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
        assert_eq!(spectra.len(), 1);
        let ms1 = reader.take_run_ms1();
        assert_eq!(ms1[0].id, "scan=1");
        assert_eq!(ms1[0].faims_cv, Some(-45.0));
        let contexts = reader.take_scan_contexts();
        assert_eq!(
            contexts,
            vec![ScanContext {
                id: "scan=2".into(),
                parent_id: Some("scan=1".into()),
                scan_window: Some((148.0, 1866.0)),
                faims_cv: Some(-45.0),
                noise: None,
                ion_mobility: None,
            }]
        );
        let mut plain = MzMLReader::new(Cursor::new(xml));
        assert_eq!(plain.by_ref().count(), 1);
        assert!(plain.take_scan_contexts().is_empty());
    }

    #[test]
    fn product_scan_capture_links_ms3_to_its_ms2_and_hides_it() {
        let ms2 = ms2_spectrum_xml(
            "controllerType=0 controllerNumber=1 scan=10",
            &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            350.0,
            Some(2),
        );
        let ms3 = format!(
            r#"<spectrum index="1" id="controllerType=0 controllerNumber=1 scan=11" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="2.0"
                         unitAccession="UO:0000031" unitName="minute"/>
              </scan></scanList>
              <precursorList count="2">
                <precursor spectrumRef="controllerType=0 controllerNumber=1 scan=10">
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  </selectedIon></selectedIonList>
                  <activation><cvParam accession="MS:1000422" name="HCD" value=""/></activation>
                </precursor>
                <precursor spectrumRef="controllerType=0 controllerNumber=1 scan=10">
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="612.2"/>
                  </selectedIon></selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
            mzb = bda_plain("MS:1000514", &encode_f64_b64(&[126.1277, 127.1248])),
            intb = bda_plain("MS:1000515", &encode_f64_b64(&[1000.0, 2000.0])),
        );
        let xml = wrap_spectra(&format!("{ms2}{ms3}"));
        let mut plain = MzMLReader::new(Cursor::new(xml.clone()));
        let n = plain
            .by_ref()
            .collect::<Result<Vec<_>, _>>()
            .expect("parse")
            .len();
        assert_eq!(n, 1);
        assert!(plain.take_product_scans().is_empty());
        let mut reader = MzMLReader::new(Cursor::new(xml)).with_product_scan_capture(true);
        let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].scan, Some(10));
        let ms3s = reader.take_product_scans();
        assert_eq!(ms3s.len(), 1);
        let p = &ms3s[0];
        assert_eq!(p.scan, Some(11));
        assert_eq!(p.ms_level, 3);
        assert_eq!(
            p.parent_id.as_deref(),
            Some("controllerType=0 controllerNumber=1 scan=10")
        );
        assert_eq!(p.rt_seconds, Some(120.0));
        assert_eq!(
            p.precursor_mz, 500.5,
            "the first selected ion gives the precursor m/z, as the first \
             precursor gives the parent reference"
        );
        assert_eq!(p.peaks.len(), 2);
        assert_eq!(p.peaks[1], (127.1248, 2000.0));
    }

    #[test]
    fn noise_capture_reads_the_sampled_noise_arrays() {
        let noise_arrays = format!(
            "{}{}",
            bda_plain("MS:1002743", &encode_f64_b64(&[120.0, 130.0, 200.0, 900.0])),
            bda_plain("MS:1002744", &encode_f64_b64(&[100.0, 200.0, 300.0, 50.0])),
        );
        let ms2 = format!(
            r#"<spectrum index="0" id="scan=10" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1"><precursor>
                <selectedIonList count="1"><selectedIon>
                  <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  <cvParam accession="MS:1000041" name="charge state" value="2"/>
                </selectedIon></selectedIonList>
              </precursor></precursorList>
              <binaryDataArrayList count="4">{mz}{int}{noise_arrays}</binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &encode_f64_b64(&[126.1277])),
            int = bda_plain("MS:1000515", &encode_f64_b64(&[1000.0])),
        );
        let ms3 = format!(
            r#"<spectrum index="1" id="scan=11" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1"><precursor spectrumRef="scan=10">
                <selectedIonList count="1"><selectedIon>
                  <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                </selectedIon></selectedIonList>
              </precursor></precursorList>
              <binaryDataArrayList count="4">{mz}{int}{noise_arrays}</binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &encode_f64_b64(&[126.1277])),
            int = bda_plain("MS:1000515", &encode_f64_b64(&[1000.0])),
        );
        let xml = wrap_spectra(&format!("{ms2}{ms3}"));
        let mut reader = MzMLReader::new(Cursor::new(xml.clone()))
            .with_product_scan_capture(true)
            .with_scan_context_capture(true)
            .with_noise_capture(true);
        let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].peaks, vec![(126.1277, 1000.0)]);
        let contexts = reader.take_scan_contexts();
        let ms3s = reader.take_product_scans();
        for noise in [&contexts[0].noise, &ms3s[0].noise] {
            let noise = noise.as_ref().expect("noise captured");
            assert_eq!(noise.len(), 3);
            assert!((noise.at(125.0) - 150.0).abs() < 1e-4);
        }
        let mut plain = MzMLReader::new(Cursor::new(xml))
            .with_product_scan_capture(true)
            .with_scan_context_capture(true);
        let _: Vec<Spectrum> = plain.by_ref().map(|r| r.expect("parse")).collect();
        assert!(plain.take_scan_contexts()[0].noise.is_none());
        assert!(plain.take_product_scans()[0].noise.is_none());
    }

    #[test]
    fn chunked_ms1_linkage_matches_batch_across_boundaries() {
        let ms1a = format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">{mz}{int}</binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &encode_f64_b64(&[111.0, 222.0])),
            int = bda_plain("MS:1000515", &encode_f64_b64(&[10.0, 20.0])),
        );
        let ms1b = format!(
            r#"<spectrum index="4" id="scan=5" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">{mz}{int}</binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &encode_f64_b64(&[333.0])),
            int = bda_plain("MS:1000515", &encode_f64_b64(&[30.0])),
        );
        let mk_ms2 = |id: &str, mz: f64| {
            ms2_spectrum_xml(
                id,
                &bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
                &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
                mz + 50.0,
                Some(2),
            )
        };
        let xml = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="6" defaultDataProcessingRef="dp">
{a}{m1}{m2}{m3}{b}{m4}
</spectrumList></run></mzML>"#,
            a = ms1a,
            m1 = mk_ms2("scan=2", 300.0),
            m2 = mk_ms2("scan=3", 400.0),
            m3 = mk_ms2("scan=4", 410.0),
            b = ms1b,
            m4 = mk_ms2("scan=6", 600.0),
        );

        let (b_spectra, b_link) = MzMLReader::new(Cursor::new(xml.clone()))
            .with_ms1_capture(true)
            .read_with_ms1()
            .expect("batch read");
        let batch_resolved: Vec<Option<&Vec<(f64, f32)>>> = b_link
            .ms2_to_ms1
            .iter()
            .map(|o| o.map(|i| &b_link.ms1_peaks[i]))
            .collect();

        let mut chunked_spectra: Vec<Spectrum> = Vec::new();
        let mut chunked_resolved: Vec<Option<Vec<(f64, f32)>>> = Vec::new();
        let (errc, _errs) = MzMLReader::new(Cursor::new(xml))
            .with_ms1_capture(true)
            .read_with_ms1_chunked(2, usize::MAX, |chunk, link| {
                for (i, s) in chunk.into_iter().enumerate() {
                    chunked_resolved.push(link.ms2_to_ms1[i].map(|j| link.ms1_peaks[j].clone()));
                    chunked_spectra.push(s);
                }
            });

        assert_eq!(errc, 0, "clean input → no parse errors");
        assert_eq!(
            chunked_spectra.len(),
            b_spectra.len(),
            "same MS2 count as batch"
        );
        assert_eq!(chunked_resolved.len(), batch_resolved.len());
        for (i, (cr, br)) in chunked_resolved
            .iter()
            .zip(batch_resolved.iter())
            .enumerate()
        {
            match (cr, br) {
                (Some(c), Some(b)) => {
                    assert_eq!(c, *b, "MS2 #{i} linked to a different MS1 than batch")
                }
                (None, None) => {}
                _ => panic!("MS2 #{i} linkage presence differs from batch: {cr:?} vs {br:?}"),
            }
        }
        assert_eq!(chunked_resolved[2].as_ref().map(|v| v.len()), Some(2));
        assert_eq!(chunked_resolved[3].as_ref().map(|v| v.len()), Some(1));
    }

    #[test]
    fn chunked_read_keeps_run_ms1_and_scan_contexts_for_quantification() {
        let ms1 = |id: &str, rt_min: f64, mz: f64| {
            format!(
                r#"<spectrum index="0" id="{id}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="{rt_min}"
                         unitAccession="UO:0000031" unitName="minute"/>
              </scan></scanList>
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
                mzb = bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
                intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
            )
        };
        let ms2 = |id: &str, mz: f64| {
            ms2_spectrum_xml(
                id,
                &bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
                &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
                mz + 50.0,
                Some(2),
            )
        };
        let xml = wrap_spectra(&format!(
            "{}{}{}{}{}{}",
            ms1("scan=1", 1.0, 400.0),
            ms2("scan=2", 300.0),
            ms2("scan=3", 310.0),
            ms1("scan=4", 1.5, 401.0),
            ms2("scan=5", 320.0),
            ms1("scan=6", 2.0, 402.0),
        ));
        let mut reader = MzMLReader::new(Cursor::new(xml.clone()))
            .with_run_ms1_capture(true)
            .with_scan_context_capture(true);
        let mut linked: Vec<Option<Vec<(f64, f32)>>> = Vec::new();
        let mut n_ms2 = 0usize;
        let (errc, _) = reader.read_with_ms1_chunked(2, usize::MAX, |chunk, link| {
            n_ms2 += chunk.len();
            for i in 0..chunk.len() {
                linked.push(link.ms2_to_ms1[i].map(|j| link.ms1_peaks[j].clone()));
            }
        });
        assert_eq!(errc, 0);
        assert_eq!(n_ms2, 3, "only MS2 is emitted");
        assert_eq!(
            linked,
            vec![
                Some(vec![(400.0, 7.0)]),
                Some(vec![(400.0, 7.0)]),
                Some(vec![(401.0, 7.0)])
            ]
        );
        let run = reader.take_run_ms1();
        assert_eq!(
            run.iter()
                .map(|s| (s.id.as_str(), s.rt))
                .collect::<Vec<_>>(),
            vec![("scan=1", 60.0), ("scan=4", 90.0), ("scan=6", 120.0)],
            "every MS1 of the run, the trailing one included"
        );
        assert_eq!(run[2].peaks, vec![(402.0, 7.0)]);
        let contexts = reader.take_scan_contexts();
        assert_eq!(
            contexts.iter().map(|c| c.id.as_str()).collect::<Vec<_>>(),
            vec!["scan=2", "scan=3", "scan=5"]
        );
        let mut plain = MzMLReader::new(Cursor::new(xml))
            .with_run_ms1_capture(true)
            .with_scan_context_capture(true);
        assert_eq!(plain.by_ref().count(), 3);
        assert_eq!(plain.take_run_ms1(), run);
        assert_eq!(plain.take_scan_contexts(), contexts);
    }

    #[test]
    fn iterator_resyncs_past_a_malformed_spectrum() {
        let good1 = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            350.0,
            Some(2),
        );
        let bad = ms2_spectrum_xml(
            "scan=2",
            &bda_plain("MS:1000514", "@@@not-valid-base64@@@"),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            450.0,
            Some(2),
        );
        let good2 = ms2_spectrum_xml(
            "scan=3",
            &bda_plain("MS:1000514", &encode_f64_b64(&[500.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            550.0,
            Some(2),
        );
        let xml = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="3" defaultDataProcessingRef="dp">
{good1}{bad}{good2}
</spectrumList></run></mzML>"#
        );

        let got = collect_ok(&xml);
        assert_eq!(
            got.len(),
            2,
            "both good spectra must survive the resync (no truncation)"
        );
        assert_eq!(got[0].scan, Some(1));
        assert_eq!(
            got[1].scan,
            Some(3),
            "the post-error spectrum must still be parsed"
        );
    }

    #[test]
    fn chunked_resyncs_past_a_malformed_spectrum() {
        let good1 = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            350.0,
            Some(2),
        );
        let bad = ms2_spectrum_xml(
            "scan=2",
            &bda_plain("MS:1000514", "@@@not-valid-base64@@@"),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            450.0,
            Some(2),
        );
        let good2 = ms2_spectrum_xml(
            "scan=3",
            &bda_plain("MS:1000514", &encode_f64_b64(&[500.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            550.0,
            Some(2),
        );
        let xml = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="3" defaultDataProcessingRef="dp">
{good1}{bad}{good2}
</spectrumList></run></mzML>"#
        );

        let mut got: Vec<Spectrum> = Vec::new();
        let (errc, errs) = MzMLReader::new(Cursor::new(xml))
            .with_ms1_capture(true)
            .read_with_ms1_chunked(5000, usize::MAX, |chunk, _link| got.extend(chunk));

        assert_eq!(errc, 1, "exactly one malformed spectrum should be counted");
        assert_eq!(errs.len(), 1, "the error message should be captured");
        assert_eq!(
            got.len(),
            2,
            "both good spectra must survive the resync (no truncation)"
        );
        assert_eq!(got[0].scan, Some(1));
        assert_eq!(
            got[1].scan,
            Some(3),
            "the post-error spectrum must still be parsed"
        );
    }

    #[test]
    fn ms1_capture_ms2_before_any_ms1_links_to_none() {
        let ms2a = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            500.0,
            Some(2),
        );
        let ms1 = format!(
            r#"<spectrum index="1" id="scan=2" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
            mz = bda_plain("MS:1000514", &encode_f64_b64(&[999.0])),
            int = bda_plain("MS:1000515", &encode_f64_b64(&[5.0])),
        );
        let ms2b = ms2_spectrum_xml(
            "scan=3",
            &bda_plain("MS:1000514", &encode_f64_b64(&[400.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[2.0])),
            600.0,
            Some(2),
        );

        let xml = format!(
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="3" defaultDataProcessingRef="dp">
      {ms2a}
      {ms1}
      {ms2b}
    </spectrumList>
  </run>
</mzML>"#
        );

        let (spectra, link) = MzMLReader::new(Cursor::new(xml))
            .with_ms1_capture(true)
            .read_with_ms1()
            .expect("read_with_ms1 failed");

        assert_eq!(spectra.len(), 2, "expected two MS2 spectra");
        assert_eq!(link.ms1_peaks.len(), 1, "one captured MS1");
        assert_eq!(link.ms2_to_ms1, vec![None, Some(0)]);
    }

    #[test]
    fn ms1_capture_off_yields_empty_link() {
        let mz_b64 = encode_f64_b64(&[100.0, 200.0]);
        let int_b64 = encode_f64_b64(&[1000.0, 500.0]);
        let spec = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(2),
        );
        let (spectra, link) = MzMLReader::new(Cursor::new(wrap_spectra(&spec)))
            .read_with_ms1()
            .expect("read_with_ms1 failed");
        assert_eq!(spectra.len(), 1);
        assert!(
            link.ms1_peaks.is_empty(),
            "no MS1 captured when capture off"
        );
        assert_eq!(link.ms2_to_ms1, vec![None]);
    }

    #[test]
    fn detect_instrument_qexactive_model_cv_param() {
        let ic = r#"<instrumentConfiguration id="IC1">
            <cvParam accession="MS:1001911" name="Q Exactive" value=""/>
            <componentList count="3">
              <source order="1">
                <cvParam accession="MS:1000398" name="nanoelectrospray" value=""/>
              </source>
              <analyzer order="2"/>
              <detector order="3"/>
            </componentList>
          </instrumentConfiguration>"#;
        let xml = wrap_with_instrument_configs(ic, "IC1", &ms2_spectrum_with_ic_ref("IC1"));
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::QExactive));
    }

    #[test]
    fn detect_instrument_astral_model_in_referenceable_param_group() {
        let xml = r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <referenceableParamGroupList count="1">
    <referenceableParamGroup id="commonInstrumentParams">
      <cvParam accession="MS:1003378" name="Orbitrap Astral" value=""/>
      <cvParam accession="MS:1000529" name="instrument serial number" value="FSN123"/>
    </referenceableParamGroup>
  </referenceableParamGroupList>
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
      <referenceableParamGroupRef ref="commonInstrumentParams"/>
      <componentList count="3">
        <source order="1">
          <cvParam accession="MS:1000398" name="nanoelectrospray" value=""/>
        </source>
        <analyzer order="2">
          <cvParam accession="MS:1000484" name="orbitrap" value=""/>
        </analyzer>
        <detector order="3"/>
      </componentList>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="IC1">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectrum}
    </spectrumList>
  </run>
</mzML>"#
            .replace("{spectrum}", &ms2_spectrum_with_ic_ref("IC1"));
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::OrbitrapAstral));
    }

    #[test]
    fn detect_instrument_fusion_ion_trap_ms2_is_low_res() {
        let doc = |ms2_ic: &str| {
            r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <referenceableParamGroupList count="1">
    <referenceableParamGroup id="common">
      <cvParam accession="MS:1002416" name="Orbitrap Fusion" value=""/>
    </referenceableParamGroup>
  </referenceableParamGroupList>
  <instrumentConfigurationList count="2">
    <instrumentConfiguration id="IC1">
      <referenceableParamGroupRef ref="common"/>
      <componentList count="1">
        <analyzer order="2"><cvParam accession="MS:1000484" name="orbitrap" value=""/></analyzer>
      </componentList>
    </instrumentConfiguration>
    <instrumentConfiguration id="IC2">
      <referenceableParamGroupRef ref="common"/>
      <componentList count="1">
        <analyzer order="2"><cvParam accession="MS:1000083" name="radial ejection linear ion trap" value=""/></analyzer>
      </componentList>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="IC1">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectrum}
    </spectrumList>
  </run>
</mzML>"#
                .replace("{spectrum}", &ms2_spectrum_with_ic_ref(ms2_ic))
        };
        assert_eq!(
            detect_instrument_type(Cursor::new(doc("IC2"))),
            Some(InstrumentType::LowRes)
        );
        assert_eq!(
            detect_instrument_type(Cursor::new(doc("IC1"))),
            Some(InstrumentType::QExactive)
        );
    }

    #[test]
    fn truncated_binary_array_is_rejected() {
        let truncated_b64 = STANDARD.encode([0u8; 12]);
        let spec = ms2_spectrum_xml(
            "scan=1",
            &bda_plain("MS:1000514", &truncated_b64),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1000.0, 2000.0])),
            500.5,
            Some(2),
        );
        let mut reader = MzMLReader::new(Cursor::new(wrap_spectra(&spec))).strict();
        let first = reader.next().expect("should yield a result");
        assert!(
            matches!(
                first,
                Err(MzMLParseError::TruncatedBinary { elem_bytes: 8, .. })
            ),
            "truncated f64 array must be rejected, got {first:?}"
        );
    }

    fn malformed_spectrum() -> String {
        ms2_spectrum_xml(
            "scan=bad",
            &bda_plain("MS:1000514", "!!!not-base64!!!"),
            &bda_plain("MS:1000515", "!!!not-base64!!!"),
            500.5,
            Some(2),
        )
    }

    #[test]
    fn strict_iterator_surfaces_malformed_spectrum_as_err() {
        let xml = wrap_spectra(&malformed_spectrum());
        let mut reader = MzMLReader::new(Cursor::new(xml)).strict();
        let first = reader.next().expect("should yield a result");
        assert!(first.is_err(), "strict mode must surface the parse error");
    }

    #[test]
    fn tolerant_iterator_skips_and_counts_malformed_spectrum() {
        let good = ms2_spectrum_xml(
            "scan=2",
            &bda_plain("MS:1000514", &encode_f64_b64(&[150.0, 300.0])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[2000.0, 1000.0])),
            500.5,
            Some(2),
        );
        let xml = wrap_spectra(&format!("{}{}", malformed_spectrum(), good));
        let mut reader = MzMLReader::new(Cursor::new(xml));
        let collected: Vec<_> = reader.by_ref().filter_map(|r| r.ok()).collect();
        assert_eq!(collected.len(), 1, "tolerant mode keeps the good spectrum");
        assert!(
            reader.skipped_count() >= 1,
            "tolerant mode records the dropped scan"
        );
    }

    fn spectrum_fingerprint(s: &Spectrum) -> String {
        format!(
            "{}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}",
            s.title,
            s.precursor_mz.to_bits(),
            s.precursor_charge,
            s.rt_seconds.map(f64::to_bits),
            s.scan,
            s.activation_method,
            s.isolation_lower_offset.map(f64::to_bits),
            s.isolation_upper_offset.map(f64::to_bits),
            s.isolation_target_mz.map(f64::to_bits),
            s.precursor_intensity.map(f32::to_bits),
            s.peaks
                .iter()
                .map(|p| (p.0.to_bits(), p.1.to_bits()))
                .collect::<Vec<_>>()
        )
    }

    #[test]
    fn parallel_read_matches_the_sequential_reader() {
        use std::io::Read;
        let gz = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../test-fixtures/orbitrap_lumos_120.mzML.gz");
        let mut xml = Vec::new();
        flate2::read::GzDecoder::new(std::fs::File::open(&gz).expect("open fixture"))
            .read_to_end(&mut xml)
            .expect("decompress fixture");
        let path =
            std::env::temp_dir().join(format!("andes-parallel-mzml-{}.mzML", std::process::id()));
        std::fs::write(&path, &xml).expect("write fixture copy");
        for (min, max) in [(2, 2), (1, 2), (1, 1)] {
            let sequential: Vec<String> = MzMLReader::new(xml.as_slice())
                .with_ms_level_range(min, max)
                .map(|r| spectrum_fingerprint(&r.expect("parse")))
                .collect();
            if (min, max) == (2, 2) {
                assert!(!sequential.is_empty(), "the fixture carries MS2 spectra");
            }
            for piece in [1_000u64, 7_919, 65_536, 1 << 30] {
                for threads in [1, 3] {
                    let (spectra, errors) =
                        read_spectra_in_pieces(&path, min, max, threads, piece).expect("read");
                    assert!(errors.is_empty(), "{errors:?}");
                    let parallel: Vec<String> = spectra.iter().map(spectrum_fingerprint).collect();
                    assert_eq!(
                        parallel, sequential,
                        "levels {min}-{max}, piece {piece}, {threads} threads"
                    );
                }
            }
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn an_ms2_scan_links_to_the_latest_ms1_at_its_faims_voltage() {
        let spectrum = |n: u32, level: u32, cv: Option<f64>, mz: f64| {
            let cv = cv.map_or(String::new(), |v| {
                format!(r#"<cvParam accession="MS:1001581" name="FAIMS compensation voltage" value="{v}"/>"#)
            });
            let precursor = if level == 2 {
                r#"<precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    <cvParam accession="MS:1000041" name="charge state" value="2"/>
                  </selectedIon></selectedIonList></precursor></precursorList>"#
            } else {
                ""
            };
            format!(
                r#"<spectrum index="{n}" id="scan={n}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="{level}"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="1.0"
                         unitAccession="UO:0000031" unitName="minute"/>
                {cv}
              </scan></scanList>
              {precursor}
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
                mzb = bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
                intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
            )
        };
        // MS1 at -50 (m/z 100), MS1 at -60 (200), MS2 -50, MS2 -60, MS1 at -50 (300),
        // MS2 -60, MS2 without a voltage.
        let xml = wrap_spectra(
            &[
                spectrum(1, 1, Some(-50.0), 100.0),
                spectrum(2, 1, Some(-60.0), 200.0),
                spectrum(3, 2, Some(-50.0), 400.0),
                spectrum(4, 2, Some(-60.0), 400.0),
                spectrum(5, 1, Some(-50.0), 300.0),
                spectrum(6, 2, Some(-60.0), 400.0),
                spectrum(7, 2, None, 400.0),
            ]
            .concat(),
        );
        let want = vec![Some(100.0), Some(200.0), Some(200.0), Some(300.0)];
        let first_mz = |link: &Ms1Link| -> Vec<Option<f64>> {
            link.ms2_to_ms1
                .iter()
                .map(|l| l.map(|m| link.ms1_peaks[m][0].0))
                .collect()
        };
        let (_, link) = MzMLReader::new(Cursor::new(xml.clone()))
            .with_ms1_capture(true)
            .read_with_ms1()
            .expect("read");
        assert_eq!(first_mz(&link), want);
        for chunk in [1, 2, 3] {
            let mut got = Vec::new();
            MzMLReader::new(Cursor::new(xml.clone())).read_with_ms1_chunked(
                chunk,
                usize::MAX,
                |_, link| got.extend(first_mz(&link)),
            );
            assert_eq!(got, want, "chunk {chunk}");
        }
        let path =
            std::env::temp_dir().join(format!("andes-faims-link-{}.mzML", std::process::id()));
        std::fs::write(&path, &xml).expect("write");
        for piece in [300u64, 1 << 30] {
            let read = read_ms1_linked_in_pieces(&path, 2, false, false, piece)
                .expect("read")
                .expect("no parse error");
            let got: Vec<Option<f64>> = read
                .ms2_to_ms1
                .iter()
                .map(|l| l.map(|m| read.run_ms1[m].peaks[0].0))
                .collect();
            assert_eq!(got, want, "piece {piece}");
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn ms1_linked_parallel_read_matches_the_chunked_reader() {
        use std::io::Read;
        let gz = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../test-fixtures/mbr_run1.mzML.gz");
        let mut xml = Vec::new();
        flate2::read::GzDecoder::new(std::fs::File::open(&gz).expect("open fixture"))
            .read_to_end(&mut xml)
            .expect("decompress fixture");
        let path = std::env::temp_dir().join(format!(
            "andes-parallel-ms1-mzml-{}.mzML",
            std::process::id()
        ));
        std::fs::write(&path, &xml).expect("write fixture copy");
        for (context, noise) in [(true, false), (false, true)] {
            let mut reader = MzMLReader::new(xml.as_slice())
                .with_ms_level_range(2, 2)
                .with_ms1_capture(true)
                .with_run_ms1_capture(true)
                .with_scan_context_capture(context)
                .with_noise_capture(noise);
            let mut spectra = Vec::new();
            let mut linked = Vec::new();
            let (errors, _) = reader.read_with_ms1_chunked(7, usize::MAX, |chunk, link| {
                for (s, l) in chunk.into_iter().zip(&link.ms2_to_ms1) {
                    spectra.push(format!("{s:?}"));
                    linked.push(l.map(|m| link.ms1_peaks[m].clone()));
                }
            });
            assert_eq!(errors, 0);
            let run_ms1 = reader.take_run_ms1();
            let contexts = reader.take_scan_contexts();
            assert!(
                run_ms1.len() > 1 && !spectra.is_empty(),
                "MS1 and MS2 scans"
            );
            for piece in [1_000u64, 7_919, 65_536, 1 << 30] {
                for threads in [1, 3] {
                    let read = read_ms1_linked_in_pieces(&path, threads, context, noise, piece)
                        .expect("read")
                        .expect("no parse error");
                    let what = format!("piece {piece}, {threads} threads");
                    let got: Vec<String> = read.spectra.iter().map(|s| format!("{s:?}")).collect();
                    assert_eq!(got, spectra, "{what}");
                    assert_eq!(read.run_ms1, run_ms1, "{what}");
                    assert_eq!(read.contexts, contexts, "{what}");
                    let got_linked: Vec<_> = read
                        .ms2_to_ms1
                        .iter()
                        .map(|l| l.map(|m| read.run_ms1[m].peaks.clone()))
                        .collect();
                    assert_eq!(got_linked, linked, "{what}");
                }
            }
        }
        let _ = std::fs::remove_file(&path);
    }
}
