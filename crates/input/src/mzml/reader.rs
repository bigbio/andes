//! The streaming mzML reader: one spectrum at a time from any `BufRead`, with optional
//! capture of MS1 scans, MS3 product scans, scan contexts and noise arrays.

use std::io::BufRead;

use quick_xml::{events::Event, Reader};

use model::scan::{Ms1Scan, NoiseProfile, ProductScan, ScanContext};
use model::{ActivationMethod, Spectrum};

use super::activation::{
    activation_user_param_accession, ethcd_as_etd, filter_has_supplemental_collision,
    log_ethcd_once,
};
use super::binary::{decode_binary_array, BinaryArrayCtx};
use super::centroid::{centroid_profile, log_profile_centroid_once};
use super::cv::{
    is_minute_unit, CV_32BIT, CV_64BIT, CV_CHARGE_STATE, CV_CID, CV_ECD, CV_ETCID, CV_ETD,
    CV_ETHCD, CV_FAIMS_CV, CV_FILTER_STRING, CV_HCD, CV_HIGH_ENERGY_BEAM_CID, CV_INTENSITY_ARRAY,
    CV_ISOLATION_LOWER_OFFSET, CV_ISOLATION_TARGET_MZ, CV_ISOLATION_UPPER_OFFSET, CV_MS_LEVEL,
    CV_MZ_ARRAY, CV_MZ_PLAIN, CV_NOISE_ARRAY, CV_NOISE_MZ_ARRAY, CV_PEAK_INTENSITY, CV_PQD,
    CV_PROFILE, CV_SCAN_TIME, CV_SCAN_WINDOW_LOWER, CV_SCAN_WINDOW_UPPER, CV_SELECTED_ION_MZ,
    CV_SUPPLEMENTAL_CID, CV_SUPPLEMENTAL_HCD, CV_TRAP_CID, CV_TRUNCATION_ZLIB, CV_UVPD,
    CV_UVPD_CHILD, CV_ZLIB, UNSUPPORTED_COMPRESSION,
};
use super::{attr_str, MzMLParseError};

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
    ms1_tracker: crate::ms1::Ms1Tracker,
    pub(super) last_faims_cv: Option<f32>,
    capture_run_ms1: bool,
    pub(super) run_ms1: Vec<Ms1Scan>,
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

    pub(super) fn pump(&mut self) -> Result<Option<Spectrum>, MzMLParseError> {
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

pub(super) fn extract_scan_from_id(id: &str) -> Option<i32> {
    ["scan=", "index=", "scanId="].iter().find_map(|key| {
        id.split_whitespace()
            .find_map(|tok| tok.strip_prefix(key)?.parse::<i32>().ok())
    })
}

pub(super) fn scan_number(id: &str, index: Option<i32>) -> Option<i32> {
    extract_scan_from_id(id).or_else(|| index.map(|i| i.saturating_add(1)))
}
