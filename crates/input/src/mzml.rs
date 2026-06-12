//! Streaming mzML reader. Event-driven via quick-xml; no serde tree.
//!
//! By default only MS2 spectra are emitted (ms level == 2). The parser
//! decodes base64 peak arrays (32-bit or 64-bit float, little-endian)
//! with optional zlib compression and zips (m/z, intensity) pairs into
//! `Vec<(f64, f32)>` sorted ascending by m/z.

use std::collections::HashMap;
use std::io::BufRead;

use base64::{engine::general_purpose::STANDARD, Engine as _};
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::ZlibDecoder;
use quick_xml::{events::Event, Reader};

use model::{ActivationMethod, InstrumentType, Spectrum};

// ── CV accessions we care about ─────────────────────────────────────────────

// Mass-analyzer cvParams used by `detect_instrument_type`. Sourced from the
// PSI-MS controlled vocabulary (HUPO-PSI MS ontology). When no `--instrument`
// flag is given, the CLI defaults to low-resolution ion-trap routing; per-file
// auto-detection reads these terms to pick a sensible bundled `.param` file
// (LTQ Velos data → CID_LowRes; Orbitrap CID → CID_HighRes).
//
// Ion-trap family → InstrumentType::LowRes.
const CV_ANALYZER_ION_TRAP:           &str = "MS:1000264"; // ion trap (generic)
const CV_ANALYZER_QUAD_ION_TRAP:      &str = "MS:1000082"; // quadrupole ion trap
const CV_ANALYZER_RADIAL_LIT:         &str = "MS:1000083"; // radial ejection linear ion trap
const CV_ANALYZER_LINEAR_ION_TRAP:    &str = "MS:1000291"; // linear ion trap
// Orbitrap / FT family → InstrumentType::QExactive / HighRes.
const CV_ANALYZER_ORBITRAP:           &str = "MS:1000484"; // orbitrap
const CV_ANALYZER_FTICR:              &str = "MS:1000079"; // Fourier transform ion cyclotron resonance
// TOF.
const CV_ANALYZER_TOF:                &str = "MS:1000084"; // time-of-flight

// Instrument-model cvParams in `<instrument>` / `<referenceableParamGroup>`
// that explicitly identify a QExactive-family box. We don't enumerate every
// Orbitrap model — falling back to "MS:1000484 orbitrap analyzer ⇒ QExactive"
// covers the typical case. These exist for cases where the analyzer cvParam
// is absent but the instrument model is recorded.
const CV_MODEL_Q_EXACTIVE:            &str = "MS:1001911";
const CV_MODEL_Q_EXACTIVE_HF:         &str = "MS:1002523";
const CV_MODEL_Q_EXACTIVE_HF_X:       &str = "MS:1002634";
const CV_MODEL_Q_EXACTIVE_PLUS:       &str = "MS:1002877";
const CV_MODEL_ORBITRAP_FUSION:       &str = "MS:1002416";
// Orbitrap Astral instrument model (PSI-MS MS:1003378). Checked BEFORE the
// generic Orbitrap-analyzer → QExactive mapping so Astral wins.
const CV_MODEL_ORBITRAP_ASTRAL:       &str = "MS:1003378";

const CV_MS_LEVEL: &str = "MS:1000511";
const CV_SCAN_TIME: &str = "MS:1000016";
const CV_SELECTED_ION_MZ: &str = "MS:1000744";
/// Isolation-window lower offset in Da (selected m/z − lower = window start).
const CV_ISOLATION_LOWER_OFFSET: &str = "MS:1000828";
/// Isolation-window upper offset in Da (selected m/z + upper = window end).
const CV_ISOLATION_UPPER_OFFSET: &str = "MS:1000829";
/// Older mzML files sometimes use plain m/z accession in selectedIon.
const CV_MZ_PLAIN: &str = "MS:1000040";
const CV_CHARGE_STATE: &str = "MS:1000041";
const CV_PEAK_INTENSITY: &str = "MS:1000042";
const CV_MZ_ARRAY: &str = "MS:1000514";
const CV_INTENSITY_ARRAY: &str = "MS:1000515";
const CV_64BIT: &str = "MS:1000523";
const CV_32BIT: &str = "MS:1000521";
const CV_ZLIB: &str = "MS:1000574";

// Activation-method CV accessions (inside <precursor><activation>).
// Mapped to the five canonical ActivationMethod variants per PSI-MS
// (HUPO-PSI MS ontology). Unknown / unhandled child terms fall through
// and the spectrum's activation_method stays None.
const CV_CID: &str  = "MS:1000133"; // collision-induced dissociation
const CV_HCD: &str  = "MS:1000422"; // beam-type CID = HCD
const CV_ETD: &str  = "MS:1000598"; // electron transfer dissociation
const CV_PQD: &str  = "MS:1000599"; // pulsed Q dissociation
const CV_UVPD: &str = "MS:1000435"; // photodissociation (PSI-MS UVPD term)
// ECD is MS:1000250; we don't have a dedicated variant for it — callers
// that need ECD usually look up either ETD or treat as electron-based.
// We map ECD → ETD for electron-based activation grouping when ECD is
// the only signal (canonical table covers ETD/CID/HCD/PQD/UVPD).
const CV_ECD: &str  = "MS:1000250"; // electron capture dissociation

/// Unit: minutes → multiply by 60 to get seconds.
const CV_UNIT_MINUTE: &str = "UO:0000031";

// ── Error type ───────────────────────────────────────────────────────────────

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
}

// io::Error → MzMLParseError via the Zlib variant.
// Cannot use #[from] because quick_xml::Error already wraps io::Error and that
// would introduce an overlapping From impl.
impl From<std::io::Error> for MzMLParseError {
    fn from(e: std::io::Error) -> Self {
        MzMLParseError::Zlib(e)
    }
}

// ── State machine ────────────────────────────────────────────────────────────

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
enum State {
    #[default]
    Outside,
    Spectrum,
    Scan,
    SelectedIon,
    /// Inside `<precursor><isolationWindow>` — we read the lower/upper
    /// offset cvParams here and set the `SpectrumBuilder` isolation fields.
    IsolationWindow,
    /// Inside `<precursor><activation>` — we read activation-method
    /// cvParams here and set `SpectrumBuilder::activation_method`.
    Activation,
    BinaryDataArray,
    Binary,
}

#[derive(Debug)]
struct BinaryArrayCtx {
    is_mz: bool,
    is_intensity: bool,
    /// 32 or 64 bits.
    precision_bits: u8,
    zlib: bool,
    b64_text: String,
}

impl BinaryArrayCtx {
    fn new() -> Self {
        BinaryArrayCtx {
            is_mz: false,
            is_intensity: false,
            precision_bits: 64,
            zlib: false,
            b64_text: String::new(),
        }
    }
}

#[derive(Debug, Default)]
struct SpectrumBuilder {
    id: String,
    ms_level: Option<u32>,
    rt_seconds: Option<f64>,
    precursor_mz: Option<f64>,
    /// Thermo-specific monoisotopic-corrected precursor m/z, when the mzML
    /// file is produced from a Thermo .raw and the instrument firmware ran
    /// its on-board deisotoping. Lives under `<scan>` as a userParam:
    ///   `<userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" value="..."/>`
    /// When present, this is preferred over `selectedIon.MS:1000744` because
    /// the raw isolation m/z may be off-by-one-or-more C13 isotopes for
    /// Orbitrap-style data — Thermo firmware deisotoping is preferred over
    /// the raw isolation m/z when the trailer userParam is present.
    monoisotopic_mz_override: Option<f64>,
    precursor_charge: Option<i32>,
    precursor_intensity: Option<f32>,
    /// Activation method recorded under `<precursor><activation>` — set
    /// when we see a known cvParam (CID/HCD/ETD/PQD/UVPD/ECD). Stays
    /// `None` when no `<activation>` block is present or the term is
    /// unknown.
    activation_method: Option<ActivationMethod>,
    /// Isolation-window lower offset in Da, from `<isolationWindow>`
    /// `MS:1000828`. `None` when the mzML omits the isolation window.
    isolation_lower_offset: Option<f64>,
    /// Isolation-window upper offset in Da, from `<isolationWindow>`
    /// `MS:1000829`. `None` when the mzML omits the isolation window.
    isolation_upper_offset: Option<f64>,
    mz_array: Option<Vec<f64>>,
    intensity_array: Option<Vec<f64>>,
}

// ── Extracted cv-param info (avoids borrow-checker conflicts) ────────────────

/// What we extract from a `<cvParam>` element without holding a reference
/// into the event buffer.
struct CvParamInfo {
    accession: String,
    value: String,
    unit_accession: String,
}

impl CvParamInfo {
    fn from_bytes_start(e: &quick_xml::events::BytesStart<'_>) -> Option<Self> {
        let accession = attr_str(e, b"accession")?;
        let value = attr_str(e, b"value").unwrap_or_default();
        let unit_accession = attr_str(e, b"unitAccession").unwrap_or_default();
        Some(CvParamInfo { accession, value, unit_accession })
    }
}

// ── MS1 linkage side structure ───────────────────────────────────────────────

/// Side structure produced by the opt-in MS1-capture reader path
/// ([`MzMLReader::read_with_ms1`]). It records the peak lists of the captured
/// MS1 scans (stored ONCE) and, for each emitted MS2 spectrum, the index of
/// its most-recent preceding MS1.
///
/// The chimeric MS1 isotope filter consumes this to subtract co-isolated
/// precursors from the MS2 candidate space. Producing it is opt-in; the
/// default reader path never builds one and stays MS2-only.
#[derive(Debug, Default, Clone)]
pub struct Ms1Link {
    /// Peak lists of the captured MS1 scans, in file order. Stored ONCE
    /// (not duplicated per MS2).
    pub ms1_peaks: Vec<Vec<(f64, f32)>>,
    /// For each emitted MS2 (by its index in the returned `Vec<Spectrum>`),
    /// the index into `ms1_peaks` of the most-recent preceding MS1, or
    /// `None` if no MS1 preceded it. Length equals the number of emitted MS2.
    pub ms2_to_ms1: Vec<Option<usize>>,
}

// ── Public reader ────────────────────────────────────────────────────────────

/// Streaming mzML reader. Emits MS2 spectra by default.
pub struct MzMLReader<R: BufRead> {
    xml: Reader<R>,
    buf: Vec<u8>,
    ms_level_min: u32,
    ms_level_max: u32,
    state: State,
    current: Option<SpectrumBuilder>,
    binary_ctx: Option<BinaryArrayCtx>,
    done: bool,
    /// When `true`, [`Self::read_with_ms1`] additionally captures MS1 scans
    /// (level 1) into an [`Ms1Link`] instead of discarding them. Has NO effect
    /// on the default [`Iterator`] path or on `read_with_ms1` emitting MS2-only
    /// — it only controls whether MS1 peaks are retained for linkage.
    capture_ms1: bool,
    /// MS1 peak lists captured during a `read_with_ms1` run (capture on).
    /// Stored once; the default path never touches this.
    captured_ms1: Vec<Vec<(f64, f32)>>,
    /// Index (into `captured_ms1`) of the most-recent preceding MS1, or
    /// `None` if no MS1 has been seen yet. Read when emitting each MS2 to
    /// build `Ms1Link::ms2_to_ms1`.
    latest_ms1_idx: Option<usize>,
}

impl<R: BufRead> MzMLReader<R> {
    /// Create a reader that emits MS2 spectra (ms level == 2).
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
            latest_ms1_idx: None,
        }
    }

    /// Widen or narrow the ms-level filter (e.g. `with_ms_level_range(1, 2)`
    /// emits both MS1 and MS2).
    pub fn with_ms_level_range(mut self, min: u32, max: u32) -> Self {
        self.ms_level_min = min;
        self.ms_level_max = max;
        self
    }

    /// Opt in to MS1 capture for the [`Self::read_with_ms1`] path. When
    /// enabled, MS1 scans (level 1) are retained in the returned [`Ms1Link`]
    /// and each emitted MS2 is linked to its most-recent preceding MS1.
    ///
    /// This flag has NO effect on the default [`Iterator`] path, which always
    /// emits MS2-only. Default is `false` (no capture).
    pub fn with_ms1_capture(mut self, capture: bool) -> Self {
        self.capture_ms1 = capture;
        self
    }

    // ── Build a Spectrum from a completed SpectrumBuilder ────────────────────

    fn finish_spectrum(&self, sb: SpectrumBuilder) -> Result<Option<Spectrum>, MzMLParseError> {
        let level = sb.ms_level.unwrap_or(0);
        if level < self.ms_level_min || level > self.ms_level_max {
            return Ok(None);
        }
        Self::build_spectrum(sb)
    }

    /// Build a [`Spectrum`] from a completed builder, applying the
    /// Thermo-trailer precursor preference and the ascending-m/z invariant.
    /// Returns `Ok(None)` for spectra without any precursor m/z (e.g. an MS2
    /// that never recorded a precursor). Caller is responsible for the
    /// ms-level filter (see [`Self::finish_spectrum`]).
    fn build_spectrum(sb: SpectrumBuilder) -> Result<Option<Spectrum>, MzMLParseError> {
        // Prefer the Thermo Trailer Extra monoisotopic m/z when available —
        // the instrument firmware's deisotoping is more accurate than the
        // raw isolation m/z (selectedIon.MS:1000744) for Orbitrap-class
        // data. Falls back to the selected ion when the trailer is absent.
        // Kim et al. (Nat Commun 5:5277, 2014) use the deisotoped precursor
        // mass for Orbitrap-class data when available.
        let precursor_mz = match (sb.monoisotopic_mz_override, sb.precursor_mz) {
            (Some(m), _) => m,
            (None, Some(v)) => v,
            // MS2 without any precursor m/z: skip rather than error.
            (None, None) => return Ok(None),
        };

        // Reject a non-finite or non-positive precursor m/z (garbage from a
        // malformed file) — it would produce nonsense search windows. The
        // Thermo and timsTOF readers already gate on `precursor_mz > 0`; do the
        // same here so all formats behave consistently.
        if !precursor_mz.is_finite() || precursor_mz <= 0.0 {
            return Ok(None);
        }

        let peaks = Self::build_peaks(sb.mz_array, sb.intensity_array)?;
        let scan = extract_scan_from_id(&sb.id);

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
        }))
    }

    /// Zip the decoded m/z and intensity arrays into ascending-by-m/z peaks.
    /// Shared by [`Self::build_spectrum`] and the MS1-capture path.
    fn build_peaks(
        mz_array: Option<Vec<f64>>,
        intensity_array: Option<Vec<f64>>,
    ) -> Result<Vec<(f64, f32)>, MzMLParseError> {
        let mz_vals = mz_array.unwrap_or_default();
        let int_vals = intensity_array.unwrap_or_default();

        if mz_vals.len() != int_vals.len() {
            return Err(MzMLParseError::LengthMismatch {
                mz_len: mz_vals.len(),
                int_len: int_vals.len(),
            });
        }

        // Drop non-finite or non-positive-m/z / negative-intensity points
        // before they can reach the scorer (matches the timsTOF reader, which
        // already filters; mzML previously passed NaN/Inf/garbage through).
        let mut peaks: Vec<(f64, f32)> = mz_vals
            .into_iter()
            .zip(int_vals)
            .map(|(mz, inten)| (mz, inten as f32))
            .filter(|&(mz, inten)| mz.is_finite() && mz > 0.0 && inten.is_finite() && inten >= 0.0)
            .collect();

        // Enforce ascending-by-m/z invariant required by downstream consumers.
        if !peaks.windows(2).all(|w| w[0].0 <= w[1].0) {
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        }

        Ok(peaks)
    }

    // ── Apply a CvParamInfo to current parse state ───────────────────────────

    fn apply_cv_param(&mut self, cv: CvParamInfo) {
        match cv.accession.as_str() {
            CV_MS_LEVEL => {
                if let Ok(lvl) = cv.value.parse::<u32>() {
                    if let Some(sb) = self.current.as_mut() {
                        sb.ms_level = Some(lvl);
                    }
                }
            }

            CV_SCAN_TIME if matches!(self.state, State::Scan | State::Spectrum) => {
                if let Ok(t) = cv.value.parse::<f64>() {
                    let secs = if cv.unit_accession == CV_UNIT_MINUTE {
                        t * 60.0
                    } else {
                        t
                    };
                    if let Some(sb) = self.current.as_mut() {
                        sb.rt_seconds = Some(secs);
                    }
                }
            }

            CV_SELECTED_ION_MZ if self.state == State::SelectedIon => {
                if let Ok(mz) = cv.value.parse::<f64>() {
                    if let Some(sb) = self.current.as_mut() {
                        sb.precursor_mz = Some(mz);
                    }
                }
            }

            CV_MZ_PLAIN if self.state == State::SelectedIon => {
                if let Ok(mz) = cv.value.parse::<f64>() {
                    if let Some(sb) = self.current.as_mut() {
                        if sb.precursor_mz.is_none() {
                            sb.precursor_mz = Some(mz);
                        }
                    }
                }
            }

            CV_CHARGE_STATE if self.state == State::SelectedIon => {
                if let Ok(z) = cv.value.parse::<i32>() {
                    if let Some(sb) = self.current.as_mut() {
                        // A charge state of 0 means "unknown" — store it as
                        // absent (None) rather than Some(0), which would later
                        // be treated as a real charge.
                        sb.precursor_charge = (z != 0).then_some(z);
                    }
                }
            }

            CV_PEAK_INTENSITY if self.state == State::SelectedIon => {
                if let Ok(inten) = cv.value.parse::<f32>() {
                    if let Some(sb) = self.current.as_mut() {
                        sb.precursor_intensity = Some(inten);
                    }
                }
            }

            // Isolation-window offsets under `<precursor><isolationWindow>`.
            // These are only consumed by the `--chimeric` path; absent in
            // most fixtures, so the builder fields stay `None`.
            CV_ISOLATION_LOWER_OFFSET if self.state == State::IsolationWindow => {
                if let Ok(off) = cv.value.parse::<f64>() {
                    if let Some(sb) = self.current.as_mut() {
                        sb.isolation_lower_offset = Some(off);
                    }
                }
            }
            CV_ISOLATION_UPPER_OFFSET if self.state == State::IsolationWindow => {
                if let Ok(off) = cv.value.parse::<f64>() {
                    if let Some(sb) = self.current.as_mut() {
                        sb.isolation_upper_offset = Some(off);
                    }
                }
            }

            // Activation-method cvParams under <precursor><activation>.
            // PSI-MS accessions map to the five canonical variants. ECD
            // (MS:1000250) is grouped with ETD for electron-based param routing.
            //
            // Selection rule:
            //   - ETD always wins (unconditional override).
            //   - Other methods: first-wins. A spectrum with multiple
            //     `<precursor><activation>` blocks (MS3 SPS, supplementary
            //     activation) records the first activation we see.
            //
            // Why first-wins matters: TMT SPS-MS3 mzMLs chain CID (MS2
            // isolation) → HCD (MS3 fragmentation). First-wins routes those
            // to a CID-trained model.
            CV_CID  if self.state == State::Activation => {
                if let Some(sb) = self.current.as_mut() {
                    if sb.activation_method.is_none() {
                        sb.activation_method = Some(ActivationMethod::CID);
                    }
                }
            }
            CV_HCD  if self.state == State::Activation => {
                if let Some(sb) = self.current.as_mut() {
                    if sb.activation_method.is_none() {
                        sb.activation_method = Some(ActivationMethod::HCD);
                    }
                }
            }
            CV_ETD  if self.state == State::Activation => {
                // ETD wins unconditionally over other activation methods.
                if let Some(sb) = self.current.as_mut() {
                    sb.activation_method = Some(ActivationMethod::ETD);
                }
            }
            CV_ECD  if self.state == State::Activation => {
                // ECD is electron-based — group with ETD for param routing.
                if let Some(sb) = self.current.as_mut() {
                    if sb.activation_method.is_none() {
                        sb.activation_method = Some(ActivationMethod::ETD);
                    }
                }
            }
            CV_PQD  if self.state == State::Activation => {
                if let Some(sb) = self.current.as_mut() {
                    if sb.activation_method.is_none() {
                        sb.activation_method = Some(ActivationMethod::PQD);
                    }
                }
            }
            CV_UVPD if self.state == State::Activation => {
                if let Some(sb) = self.current.as_mut() {
                    if sb.activation_method.is_none() {
                        sb.activation_method = Some(ActivationMethod::UVPD);
                    }
                }
            }

            CV_MZ_ARRAY if self.state == State::BinaryDataArray => {
                if let Some(ctx) = self.binary_ctx.as_mut() {
                    ctx.is_mz = true;
                }
            }
            CV_INTENSITY_ARRAY if self.state == State::BinaryDataArray => {
                if let Some(ctx) = self.binary_ctx.as_mut() {
                    ctx.is_intensity = true;
                }
            }
            CV_64BIT if self.state == State::BinaryDataArray => {
                if let Some(ctx) = self.binary_ctx.as_mut() {
                    ctx.precision_bits = 64;
                }
            }
            CV_32BIT if self.state == State::BinaryDataArray => {
                if let Some(ctx) = self.binary_ctx.as_mut() {
                    ctx.precision_bits = 32;
                }
            }
            CV_ZLIB if self.state == State::BinaryDataArray => {
                if let Some(ctx) = self.binary_ctx.as_mut() {
                    ctx.zlib = true;
                }
            }

            _ => {}
        }
    }

    // ── Event pump ───────────────────────────────────────────────────────────

    fn pump(&mut self) -> Result<Option<Spectrum>, MzMLParseError> {
        loop {
            self.buf.clear();
            // Read the next event. The lifetime of `event` is tied to `self.buf`,
            // so we must *not* hold onto it across a `&mut self` method call.
            // We extract what we need (as owned Strings) before calling helpers.
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
                            self.current =
                                Some(SpectrumBuilder { id, ..Default::default() });
                            self.state = State::Spectrum;
                        }
                        b"scan" if self.state == State::Spectrum => {
                            self.state = State::Scan;
                        }
                        b"selectedIon" if self.state == State::Spectrum => {
                            self.state = State::SelectedIon;
                        }
                        b"isolationWindow" if self.state == State::Spectrum => {
                            // `<isolationWindow>` is a sibling of
                            // `<selectedIon>` / `<activation>` under
                            // `<precursor>`. We don't track the intermediate
                            // `<precursor>` / `<precursorList>` elements, so
                            // we transition from Spectrum here. The closing
                            // tag pops us back to Spectrum.
                            self.state = State::IsolationWindow;
                        }
                        b"activation" if self.state == State::Spectrum => {
                            // `<activation>` lives under
                            // `<precursorList><precursor>…</precursor></precursorList>`.
                            // We don't track the intermediate `<precursor>` /
                            // `<precursorList>` elements, so we transition
                            // from Spectrum here. The closing tag pops us
                            // back to Spectrum.
                            self.state = State::Activation;
                        }
                        b"binaryDataArray" if self.state == State::Spectrum => {
                            self.binary_ctx = Some(BinaryArrayCtx::new());
                            self.state = State::BinaryDataArray;
                        }
                        b"binary" if self.state == State::BinaryDataArray => {
                            self.state = State::Binary;
                        }
                        _ => {}
                    }
                }

                // Self-closing elements — mostly cvParam and userParam.
                Event::Empty(ref e) => {
                    let tag = e.local_name().as_ref().to_owned();
                    if tag == b"cvParam" {
                        // Extract info before any &mut self call.
                        if let Some(cv) = CvParamInfo::from_bytes_start(e) {
                            self.apply_cv_param(cv);
                        }
                    } else if tag == b"userParam"
                        && matches!(self.state, State::Scan | State::Spectrum)
                    {
                        // The only userParam we care about is the Thermo
                        // monoisotopic-correction recorded by the instrument
                        // firmware. It lives under <scan> and is preferred
                        // over selectedIon.MS:1000744 (raw isolation m/z) when
                        // present. Orbitrap precursors are routinely
                        // mis-isotoped by the isolation logic, and the
                        // Trailer Extra value carries the deisotoped C0 peak.
                        if let (Some(name), Some(val)) =
                            (attr_str(e, b"name"), attr_str(e, b"value"))
                        {
                            // Match either the canonical Thermo string or a
                            // few near-equivalent forms seen in older
                            // proteomics workflows (case-insensitive on the
                            // "Monoisotopic" word). Strict accept-list — no
                            // unrelated userParams sneak through.
                            let normalized = name.to_lowercase();
                            if normalized.contains("monoisotopic m/z")
                                || normalized.contains("monoisotopic mz")
                            {
                                if let Ok(mz) = val.parse::<f64>() {
                                    // mzML files sometimes emit "0" or
                                    // negative sentinels when the firmware
                                    // couldn't decide. Treat as absent.
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

                Event::Text(ref e) if self.state == State::Binary => {
                    let chunk = e.unescape()?;
                    if let Some(ctx) = self.binary_ctx.as_mut() {
                        ctx.b64_text.push_str(chunk.as_ref());
                    }
                }

                Event::End(ref e) => {
                    let tag = e.local_name().as_ref().to_owned();
                    match tag.as_slice() {
                        b"spectrum" => {
                            let sb = self.current.take();
                            self.state = State::Outside;
                            if let Some(sb) = sb {
                                // Capture path: intercept MS1 scans before they
                                // reach the level filter. MS1 peaks are stored
                                // once in `captured_ms1` and never emitted as a
                                // scorable Spectrum; we just remember the latest
                                // MS1 index so the next MS2 can link to it.
                                // This branch is inert when `capture_ms1` is
                                // false (default), keeping that path byte-exact.
                                if self.capture_ms1 && sb.ms_level == Some(1) {
                                    let peaks =
                                        Self::build_peaks(sb.mz_array, sb.intensity_array)?;
                                    self.captured_ms1.push(peaks);
                                    self.latest_ms1_idx = Some(self.captured_ms1.len() - 1);
                                    continue;
                                }
                                if let Some(s) = self.finish_spectrum(sb)? {
                                    return Ok(Some(s));
                                }
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
                        b"activation" if self.state == State::Activation => {
                            self.state = State::Spectrum;
                        }
                        b"binary" if self.state == State::Binary => {
                            self.state = State::BinaryDataArray;
                        }
                        b"binaryDataArray" if self.state == State::BinaryDataArray => {
                            if let Some(ctx) = self.binary_ctx.take() {
                                let vals = decode_binary_array(&ctx)?;
                                if let Some(sb) = self.current.as_mut() {
                                    if ctx.is_mz {
                                        sb.mz_array = Some(vals);
                                    } else if ctx.is_intensity {
                                        sb.intensity_array = Some(vals);
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

    /// Drain the reader, returning the emitted MS2 spectra plus an
    /// [`Ms1Link`] linking each MS2 to its most-recent preceding MS1.
    ///
    /// The returned `Vec<Spectrum>` is exactly what the default [`Iterator`]
    /// path produces (MS2-only, in file order). When [`Self::with_ms1_capture`]
    /// is `false` (default), `Ms1Link::ms1_peaks` is empty and every entry in
    /// `ms2_to_ms1` is `None`; only the linkage bookkeeping is added on top of
    /// the unchanged MS2 stream. When capture is `true`, MS1 scans are
    /// retained once in `ms1_peaks` and each MS2 links to its preceding MS1.
    ///
    /// `ms2_to_ms1.len()` always equals the number of returned MS2 spectra.
    pub fn read_with_ms1(mut self) -> std::io::Result<(Vec<Spectrum>, Ms1Link)> {
        // To capture MS1 (level 1) the parser must let level-1 spectra reach
        // the spectrum-End handler rather than being dropped by the level
        // filter. We widen the internal min level to 1 ONLY when capturing;
        // MS1 is intercepted before `finish_spectrum`, so the effective output
        // is still MS2-only. With capture off, the filter is untouched.
        if self.capture_ms1 {
            self.ms_level_min = 1;
        }

        let mut spectra: Vec<Spectrum> = Vec::new();
        let mut ms2_to_ms1: Vec<Option<usize>> = Vec::new();

        loop {
            match self.pump() {
                Ok(Some(s)) => {
                    // Each spectrum returned by `pump` here is an emitted MS2
                    // (MS1 is intercepted inside `pump` and never returned).
                    // Link it to whatever MS1 most recently preceded it.
                    ms2_to_ms1.push(self.latest_ms1_idx);
                    spectra.push(s);
                }
                Ok(None) => break,
                Err(_e) => {
                    // Resync past the malformed spectrum and keep parsing (skip the
                    // bad scan, not the rest of the file). Only an unreadable XML
                    // stream stops us.
                    match self.resync_to_next_spectrum() {
                        Ok(true) => continue,
                        Ok(false) | Err(_) => break,
                    }
                }
            }
        }

        let link = Ms1Link {
            ms1_peaks: std::mem::take(&mut self.captured_ms1),
            ms2_to_ms1,
        };
        Ok((spectra, link))
    }

    /// Streaming, bounded-memory, tolerant variant of [`Self::read_with_ms1`] for
    /// the chimeric cascade.
    ///
    /// Calls `on_chunk(ms2_spectra, ms1_link)` for each batch of up to
    /// `chunk_size` MS2 spectra, where `ms1_link` covers ONLY that chunk. RSS
    /// stays bounded by the chunk size: at most the MS1 scans referenced by the
    /// in-flight chunk are retained, never the whole file (each MS2 links to its
    /// most-recent preceding MS1, so only that carry-over scan crosses a chunk
    /// boundary). Stops after `cap` total MS2 (`usize::MAX` = unbounded).
    ///
    /// Tolerant: a malformed spectrum does NOT abort the run. The first parse
    /// error stops streaming and the successfully-parsed spectra so far are still
    /// delivered (mirroring the MS2-only streaming path); the error count and the
    /// first few messages are returned for reporting.
    pub fn read_with_ms1_chunked<F>(
        mut self,
        chunk_size: usize,
        cap: usize,
        mut on_chunk: F,
    ) -> (usize, Vec<String>)
    where
        F: FnMut(Vec<Spectrum>, Ms1Link),
    {
        self.capture_ms1 = true;
        self.ms_level_min = 1; // let MS1 reach the capture hook; output stays MS2-only

        let mut err_count = 0usize;
        let mut first_errors: Vec<String> = Vec::new();
        let mut total = 0usize;

        let mut chunk: Vec<Spectrum> = Vec::with_capacity(chunk_size);
        let mut chunk_ms1: Vec<Vec<(f64, f32)>> = Vec::new();
        let mut links: Vec<Option<usize>> = Vec::with_capacity(chunk_size);
        // Most-recent MS1 peaks, carried across chunk boundaries. `carry_in_chunk`
        // is its index within the CURRENT `chunk_ms1`, or `None` until copied in.
        let mut carry: Option<Vec<(f64, f32)>> = None;
        let mut carry_in_chunk: Option<usize> = None;

        loop {
            if total >= cap {
                break;
            }
            match self.pump() {
                Ok(Some(ms2)) => {
                    // A newer MS1 was captured since the last MS2: adopt the most
                    // recent as the carry and clear the reader's buffer to bound RSS.
                    if self.latest_ms1_idx.is_some() {
                        carry = self.captured_ms1.pop(); // newest; older ones unreferenced
                        self.captured_ms1.clear();
                        self.latest_ms1_idx = None;
                        carry_in_chunk = None;
                    }
                    let link = if carry.is_some() {
                        if carry_in_chunk.is_none() {
                            chunk_ms1.push(carry.clone().expect("carry is Some"));
                            carry_in_chunk = Some(chunk_ms1.len() - 1);
                        }
                        carry_in_chunk
                    } else {
                        None
                    };
                    chunk.push(ms2);
                    links.push(link);
                    total += 1;

                    if chunk.len() >= chunk_size {
                        on_chunk(
                            std::mem::take(&mut chunk),
                            Ms1Link { ms1_peaks: std::mem::take(&mut chunk_ms1), ms2_to_ms1: std::mem::take(&mut links) },
                        );
                        chunk = Vec::with_capacity(chunk_size);
                        chunk_ms1 = Vec::new();
                        links = Vec::with_capacity(chunk_size);
                        carry_in_chunk = None; // re-copy carry into the next chunk on demand
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    err_count += 1;
                    if first_errors.len() < 3 {
                        first_errors.push(format!("{e}"));
                    }
                    // Resync past the bad scan and keep parsing the rest of the file
                    // (skip-one-bad-spectrum, not truncate-the-tail). Only a broken
                    // XML stream — where resync itself can't read — stops us.
                    match self.resync_to_next_spectrum() {
                        Ok(true) => continue,
                        Ok(false) | Err(_) => break,
                    }
                }
            }
        }

        if !chunk.is_empty() {
            on_chunk(chunk, Ms1Link { ms1_peaks: chunk_ms1, ms2_to_ms1: links });
        }
        (err_count, first_errors)
    }

    /// After a recoverable per-spectrum parse error, discard the partial spectrum
    /// and skip events until the next `<spectrum>` start, so streaming RESYNCS and
    /// continues past one bad scan instead of truncating the rest of the file.
    /// Returns `Ok(true)` positioned at a fresh spectrum (caller continues),
    /// `Ok(false)` at EOF. Propagates `Err` only when the XML stream itself is
    /// unreadable (broken markup) — that is genuinely unrecoverable.
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
                    self.current = Some(SpectrumBuilder { id, ..Default::default() });
                    self.state = State::Spectrum;
                    return Ok(true);
                }
                _ => {}
            }
        }
    }
}

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
                Err(_e) => {
                    // Resync past the malformed spectrum and keep parsing (skip the
                    // bad scan, not the rest of the file). Only an unreadable XML
                    // stream stops us — mirrors `read_with_ms1` / chunked paths.
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

// ── Instrument-type detection (separate, lightweight pass) ──────────────────

/// Quick mzML scan that returns the dominant
/// [`InstrumentType`] of MS2 spectra in the file.
///
/// Strategy:
/// 1. Parse `<instrumentConfigurationList>` and build a map from
///    `id` → analyzer [`InstrumentType`] using the analyzer / instrument-model
///    cvParams listed at the top of this module.
/// 2. As `<spectrum>` elements stream by, inspect their `<scan>`'s
///    `instrumentConfigurationRef=` attribute. Tally analyzer types for MS2
///    spectra only, stop after `MAX_PEEK` MS2 scans (early exit).
/// 3. Return the most-common analyzer mapped through `InstrumentType`. If no
///    MS2 scan referenced a known IC, fall back to the run-level
///    `defaultInstrumentConfigurationRef`. If nothing resolves, return `None`.
///
/// This intentionally does *not* mutate `MzMLReader`. We keep the
/// instrument-detection path as a separate, one-shot pre-pass so the main
/// streaming reader stays focused on per-spectrum data and remains
/// peak-memory-friendly.
pub fn detect_instrument_type<R: BufRead>(reader: R) -> Option<InstrumentType> {
    let mut xml = Reader::from_reader(reader);
    xml.trim_text(true);

    /// Internal scan state. Mirrors the structure of the streaming reader
    /// without sharing it, since the instrument-type detection cares about
    /// a different subset of the mzML schema.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum S {
        Outside,
        InstrumentConfigurationList,
        InstrumentConfiguration, // inside <instrumentConfiguration id="X">
        ComponentListAnalyzer,   // inside <componentList><analyzer>
        Run,
        Spectrum,
        Scan,
    }

    let mut state = S::Outside;
    let mut buf: Vec<u8> = Vec::with_capacity(4096);

    // IC id → detected InstrumentType.
    let mut ic_map: HashMap<String, InstrumentType> = HashMap::new();
    // Stored under the IC currently being parsed.
    let mut current_ic_id: Option<String> = None;
    let mut current_ic_type: Option<InstrumentType> = None;

    // run-level defaultInstrumentConfigurationRef.
    let mut default_ic_ref: Option<String> = None;

    // Tally of InstrumentType for MS2 spectra (via per-scan ref).
    let mut ms2_counts: HashMap<InstrumentType, usize> = HashMap::new();
    let mut current_spec_is_ms2: Option<bool> = None;
    let mut current_spec_ic_ref: Option<String> = None;
    let mut ms2_seen: usize = 0;

    const MAX_PEEK: usize = 64;

    loop {
        buf.clear();
        let event = match xml.read_event_into(&mut buf) {
            Ok(e) => e,
            // On parse error we just return whatever we've found so far —
            // detection is best-effort, never load-bearing for correctness.
            Err(_) => break,
        };
        match event {
            Event::Eof => break,

            Event::Start(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                match tag.as_slice() {
                    b"instrumentConfigurationList" if state == S::Outside => {
                        state = S::InstrumentConfigurationList;
                    }
                    b"instrumentConfiguration" if state == S::InstrumentConfigurationList => {
                        current_ic_id = attr_str(e, b"id");
                        current_ic_type = None;
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
                // A self-closing `<scan instrumentConfigurationRef="..."/>`
                // doesn't fire a Start event. Capture the IC ref attribute
                // here so files that emit empty `<scan/>` elements still
                // route correctly. Common in trimmed test fixtures.
                if tag == b"scan" && state == S::Spectrum {
                    if let Some(r) = attr_str(e, b"instrumentConfigurationRef") {
                        current_spec_ic_ref = Some(r);
                    }
                    // Don't transition state — the spectrum tag is still
                    // open; the End handler for `<spectrum>` consumes it.
                }
                if tag == b"cvParam" {
                    let acc = attr_str(e, b"accession").unwrap_or_default();
                    match state {
                        // Within <analyzer>: pick up the mass-analyzer cvParam.
                        S::ComponentListAnalyzer => {
                            let typ = match acc.as_str() {
                                CV_ANALYZER_ORBITRAP        => Some(InstrumentType::QExactive),
                                CV_ANALYZER_FTICR           => Some(InstrumentType::HighRes),
                                CV_ANALYZER_TOF             => Some(InstrumentType::TOF),
                                CV_ANALYZER_ION_TRAP
                                | CV_ANALYZER_QUAD_ION_TRAP
                                | CV_ANALYZER_RADIAL_LIT
                                | CV_ANALYZER_LINEAR_ION_TRAP => Some(InstrumentType::LowRes),
                                _ => None,
                            };
                            if let Some(t) = typ {
                                // First analyzer wins for a given IC when
                                // mzMLs declare more than one (PSI-MS practice).
                                if current_ic_type.is_none() {
                                    current_ic_type = Some(t);
                                }
                            }
                        }
                        // Within <instrumentConfiguration> at the top level
                        // (not inside <analyzer>): an instrument-model cvParam
                        // may be present and gives us a stronger signal for
                        // Orbitrap-class boxes than analyzer alone.
                        // OrbitrapAstral is checked FIRST so it wins over the
                        // generic QExactive group.
                        S::InstrumentConfiguration => {
                            let model = match acc.as_str() {
                                CV_MODEL_ORBITRAP_ASTRAL => Some(InstrumentType::OrbitrapAstral),
                                CV_MODEL_Q_EXACTIVE
                                | CV_MODEL_Q_EXACTIVE_HF
                                | CV_MODEL_Q_EXACTIVE_HF_X
                                | CV_MODEL_Q_EXACTIVE_PLUS
                                | CV_MODEL_ORBITRAP_FUSION => Some(InstrumentType::QExactive),
                                _ => {
                                    // Name-substring fallback: some mzML files record
                                    // the instrument model name even when using an
                                    // unregistered or future accession. If the cvParam
                                    // `name` attribute contains "astral" (case-insensitive)
                                    // treat it as OrbitrapAstral so detection is robust.
                                    let cv_name = attr_str(e, b"name").unwrap_or_default();
                                    if cv_name.to_ascii_lowercase().contains("astral") {
                                        Some(InstrumentType::OrbitrapAstral)
                                    } else {
                                        None
                                    }
                                }
                            };
                            if let Some(t) = model {
                                // Model wins outright if seen.
                                current_ic_type = Some(t);
                            }
                        }
                        // Within <spectrum>: pick up ms-level.
                        S::Spectrum => {
                            if acc == CV_MS_LEVEL {
                                let val = attr_str(e, b"value").unwrap_or_default();
                                if val == "2" {
                                    current_spec_is_ms2 = Some(true);
                                } else {
                                    current_spec_is_ms2 = Some(false);
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }

            Event::End(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                match tag.as_slice() {
                    b"analyzer" if state == S::ComponentListAnalyzer => {
                        state = S::InstrumentConfiguration;
                    }
                    b"instrumentConfiguration" if state == S::InstrumentConfiguration => {
                        if let (Some(id), Some(t)) = (current_ic_id.take(), current_ic_type.take()) {
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
                        // Tally if this was MS2 and we know its IC ref (or the
                        // file-wide default IC).
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

    // Prefer the dominant analyzer across MS2 scans.
    if !ms2_counts.is_empty() {
        return ms2_counts
            .iter()
            .max_by_key(|(_, &n)| n)
            .map(|(&t, _)| t);
    }

    // No MS2-referenced IC info — fall back to default IC if it's known.
    if let Some(r) = default_ic_ref.as_ref() {
        if let Some(&t) = ic_map.get(r) {
            return Some(t);
        }
    }

    // No default-IC info either — use the first IC we found (some mzMLs only
    // declare one IC and don't reference it from each scan).
    if ic_map.len() == 1 {
        return ic_map.into_values().next();
    }

    None
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Extract a named attribute value as an owned String.
fn attr_str(e: &quick_xml::events::BytesStart<'_>, name: &[u8]) -> Option<String> {
    e.attributes()
        .filter_map(|a| a.ok())
        .find(|a| a.key.local_name().as_ref() == name)
        .and_then(|a| std::str::from_utf8(a.value.as_ref()).ok().map(str::to_owned))
}

/// Parse the scan number from a spectrum id attribute.
///
/// Handles ProteoWizard format: `"controllerType=0 controllerNumber=1 scan=1234"`
/// and plain `"scan=1234"`.
fn extract_scan_from_id(id: &str) -> Option<i32> {
    id.split_whitespace()
        .find_map(|tok| tok.strip_prefix("scan=")?.parse::<i32>().ok())
}

/// Decode a `<binaryDataArray>` payload: base64 → optional zlib → f64 values.
fn decode_binary_array(ctx: &BinaryArrayCtx) -> Result<Vec<f64>, MzMLParseError> {
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

    let mut cur = std::io::Cursor::new(&bytes);
    let mut out: Vec<f64> = Vec::new();

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

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn collect_ok(xml: &str) -> Vec<Spectrum> {
        MzMLReader::new(Cursor::new(xml))
            .map(|r| r.expect("parse error"))
            .collect()
    }

    /// Minimal valid mzML wrapper around raw `<spectrum>` XML.
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

    // ── Encoding helpers ──────────────────────────────────────────────────────

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
            Some(z) => format!(
                r#"<cvParam accession="MS:1000041" name="charge state" value="{z}"/>"#
            ),
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

    // ── Test 1 ────────────────────────────────────────────────────────────────

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

    // ── Test 2 ────────────────────────────────────────────────────────────────

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

    // ── Test 3 ────────────────────────────────────────────────────────────────

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

    // ── Test 4 ────────────────────────────────────────────────────────────────

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

    // ── Test 5 ────────────────────────────────────────────────────────────────

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

    // ── Test 6 ────────────────────────────────────────────────────────────────

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
    fn charge_state_zero_is_treated_as_absent() {
        // A `charge state` of 0 means "unknown" and must be stored as None, not
        // Some(0), which downstream code would treat as a real charge.
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
        // NaN / Inf / negative-m/z points must be dropped before reaching the
        // scorer (the timsTOF reader already does this; mzML used to pass them
        // through). Only the two finite, positive-m/z peaks should survive.
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
        assert_eq!(mzs, vec![100.0, 200.0], "only finite positive-m/z peaks survive");
    }

    #[test]
    fn nonpositive_precursor_mz_is_skipped() {
        // A precursor m/z <= 0 is garbage and would yield nonsense search
        // windows; the spectrum must be dropped (Ok(None)).
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
        assert!(spectra.is_empty(), "spectrum with precursor m/z 0 must be skipped");
    }

    // ── Test 7 ────────────────────────────────────────────────────────────────

    #[test]
    fn peaks_sorted_ascending_by_mz() {
        // Provide peaks deliberately out of order.
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

    // ── Test 8: integration — real tiny.pwiz.mzML fixture ────────────────────

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

        // tiny.pwiz.mzML: scan=19 MS1, scan=20 MS2, scan=21 MS1, scan=22 MS1.
        // Only scan=20 should pass the default MS2 filter.
        assert!(!spectra.is_empty(), "expected at least one MS2 spectrum");
        let s = &spectra[0];
        assert!(!s.peaks.is_empty(), "MS2 spectrum should have peaks");
        assert!(s.precursor_mz > 0.0, "precursor m/z should be positive");
    }

    // ── Unit helpers for extract_scan_from_id ────────────────────────────────

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

    /// Thermo Trailer Extra `Monoisotopic M/Z` userParam under `<scan>`
    /// overrides the raw isolation m/z (`selectedIon.MS:1000744`). Kim et al.
    /// (Nat Commun 5:5277, 2014) prefer the firmware-deisotoped precursor for
    /// Thermo Orbitrap data; load-bearing for TMT / Orbitrap recall.
    #[test]
    fn thermo_trailer_monoisotopic_overrides_selected_ion_mz() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        // Same shape as ms2_spectrum_xml but with a Thermo trailer under
        // <scan>. selectedIon m/z = 625.338 (raw isolation), trailer
        // monoisotopic m/z = 625.004 (firmware deisotoping, off by 1 C13/3).
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

    /// When the Thermo trailer is absent, the reader still falls back to
    /// `selectedIon.MS:1000744`. Regression test for the existing path.
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

    /// A zero or negative trailer value (firmware "no decision" sentinel)
    /// must not override a real selectedIon m/z — otherwise we'd plant a
    /// nonsense precursor mass.
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
    }

    // ── Activation-method parsing ────────────────────────────────────────────

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
            mz  = bda_plain("MS:1000514", &mz_b64),
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
        // ECD is electron-based; we collapse to ETD for param routing.
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

    /// SPS-MS3 mzMLs chain `<precursor><activation>` blocks (CID then HCD).
    /// First-wins selection (modulo ETD precedence) routes TMT SPS data to a
    /// CID-trained model.
    #[test]
    fn multiple_activations_first_wins() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        // Two `<precursor>` blocks: first CID (MS:1000133), second HCD
        // (MS:1000422). First-wins → CID.
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
            mz  = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );

        // Wrap and widen to MS3 so the spectrum isn't filtered out.
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

    /// ETD has unconditional precedence over CID/HCD within a single
    /// `<activation>` block (PSI-MS activation precedence rule).
    #[test]
    fn etd_precedence_over_other_methods() {
        let mz_b64 = encode_f64_b64(&[100.0]);
        let int_b64 = encode_f64_b64(&[1000.0]);
        // Activation has CID first, then ETD. ETD must win.
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
            mz  = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
    }

    // ── Isolation-window parsing ─────────────────────────────────────────────

    /// `<precursor><isolationWindow>` carries the lower/upper offsets
    /// (`MS:1000828` / `MS:1000829`). The parser must capture them on the
    /// emitted `Spectrum` without disturbing the sibling `<selectedIon>`
    /// precursor m/z. Load-bearing for the `--chimeric` co-isolation path.
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
            mz  = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let spectra = collect_ok(&wrap_spectra(&xml));
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].isolation_lower_offset, Some(1.5));
        assert_eq!(spectra[0].isolation_upper_offset, Some(1.5));
        // The isolation-window target (MS:1000827) must NOT clobber the
        // selectedIon precursor m/z.
        assert!((spectra[0].precursor_mz - 500.5).abs() < 1e-6);
    }

    /// When the mzML omits `<isolationWindow>`, both offsets stay `None`.
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
    }

    // ── Instrument-type detection ────────────────────────────────────────────

    /// Build an mzML wrapper with one or more `<instrumentConfiguration>`
    /// blocks and `<run>`-level `defaultInstrumentConfigurationRef`.
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
            mz  = bda_plain("MS:1000514", &mz_b64),
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
        // Linear ion trap (MS:1000291) — LTQ Velos and similar.
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
        // Mimics PXD001819: MS1 uses IC1 (orbitrap) but MS2 uses IC2 (ion trap).
        // The MS2-referenced IC must win → LowRes.
        let ics = format!(
            "{}\n{}",
            ic_block("IC1", "MS:1000484"), // orbitrap
            ic_block("IC2", "MS:1000264"), // ion trap
        );
        // MS2 references IC2.
        let xml = wrap_with_instrument_configs(&ics, "IC1", &ms2_spectrum_with_ic_ref("IC2"));
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::LowRes));
    }

    #[test]
    fn detect_instrument_falls_back_to_default_ic_when_scan_lacks_ref() {
        // Spectrum's <scan> has no instrumentConfigurationRef — falls back to
        // run-level defaultInstrumentConfigurationRef.
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
            mz  = bda_plain("MS:1000514", &mz_b64),
            int = bda_plain("MS:1000515", &int_b64),
        );
        let xml = wrap_with_instrument_configs(&ic_block("IC1", "MS:1000484"), "IC1", &spec);
        let result = detect_instrument_type(Cursor::new(xml));
        assert_eq!(result, Some(InstrumentType::QExactive));
    }

    #[test]
    fn detect_instrument_returns_none_when_no_ic_info() {
        // No instrumentConfigurationList block at all.
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

    // ── MS1 capture + MS2→MS1 linkage (Ms1Link) ──────────────────────────────

    /// One MS1 followed by two MS2 spectra: the capture path stores the MS1
    /// peaks once, emits only the two MS2 spectra, and links both MS2 to the
    /// single preceding MS1 (index 0).
    #[test]
    fn ms1_capture_links_two_ms2_to_preceding_ms1() {
        // MS1 with two peaks.
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

        // Only the two MS2 spectra are emitted (the MS1 is NOT in the Vec).
        assert_eq!(spectra.len(), 2, "expected exactly two MS2 spectra");
        assert_eq!(spectra[0].scan, Some(2));
        assert_eq!(spectra[1].scan, Some(3));

        // MS1 peaks stored once.
        assert_eq!(link.ms1_peaks.len(), 1, "expected one captured MS1");
        let ms1_peaks = &link.ms1_peaks[0];
        assert_eq!(ms1_peaks.len(), 2);
        assert!((ms1_peaks[0].0 - 111.0).abs() < 1e-6);
        assert!((ms1_peaks[1].0 - 222.0).abs() < 1e-6);
        assert!((ms1_peaks[0].1 - 10.0_f32).abs() < 1e-3);
        assert!((ms1_peaks[1].1 - 20.0_f32).abs() < 1e-3);

        // Both MS2 link to MS1 index 0.
        assert_eq!(link.ms2_to_ms1, vec![Some(0), Some(0)]);
    }

    #[test]
    fn chunked_ms1_linkage_matches_batch_across_boundaries() {
        // Layout exercising both an MS1 carry-over across a chunk boundary and a
        // new MS1 mid-chunk:  MS1#a, MS2_1, MS2_2, MS2_3, MS1#b, MS2_4.
        // With chunk_size=2: chunk1=[1,2]→a; chunk2=[3,4] where 3 carries over a
        // and 4 links to the new b.
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

        // Batch reference: resolve each MS2 to its MS1 peak list.
        let (b_spectra, b_link) = MzMLReader::new(Cursor::new(xml.clone()))
            .with_ms1_capture(true)
            .read_with_ms1()
            .expect("batch read");
        let batch_resolved: Vec<Option<&Vec<(f64, f32)>>> = b_link
            .ms2_to_ms1
            .iter()
            .map(|o| o.map(|i| &b_link.ms1_peaks[i]))
            .collect();

        // Chunked read with a small chunk size → multiple chunks + carry-over.
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
        assert_eq!(chunked_spectra.len(), b_spectra.len(), "same MS2 count as batch");
        assert_eq!(chunked_resolved.len(), batch_resolved.len());
        for (i, (cr, br)) in chunked_resolved.iter().zip(batch_resolved.iter()).enumerate() {
            match (cr, br) {
                (Some(c), Some(b)) => assert_eq!(c, *b, "MS2 #{i} linked to a different MS1 than batch"),
                (None, None) => {}
                _ => panic!("MS2 #{i} linkage presence differs from batch: {cr:?} vs {br:?}"),
            }
        }
        // Sanity: MS2 #2 (chunk-2 first) carried over MS1#a; MS2 #3 linked to MS1#b.
        assert_eq!(chunked_resolved[2].as_ref().map(|v| v.len()), Some(2)); // a has 2 peaks
        assert_eq!(chunked_resolved[3].as_ref().map(|v| v.len()), Some(1)); // b has 1 peak
    }

    #[test]
    fn iterator_resyncs_past_a_malformed_spectrum() {
        // The default Iterator path must skip a bad scan and still yield the next
        // good spectrum — same resync semantics as `read_with_ms1` / chunked.
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
        assert_eq!(got.len(), 2, "both good spectra must survive the resync (no truncation)");
        assert_eq!(got[0].scan, Some(1));
        assert_eq!(got[1].scan, Some(3), "the post-error spectrum must still be parsed");
    }

    #[test]
    fn chunked_resyncs_past_a_malformed_spectrum() {
        // A bad MS2 (invalid base64 in its m/z array) sits between two good MS2s.
        // The reader must resync past the bad scan and still deliver BOTH good ones
        // — NOT truncate the file at the first error.
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
        assert_eq!(got.len(), 2, "both good spectra must survive the resync (no truncation)");
        assert_eq!(got[0].scan, Some(1));
        assert_eq!(got[1].scan, Some(3), "the post-error spectrum must still be parsed");
    }

    /// An MS2 that appears BEFORE any MS1 links to `None`.
    #[test]
    fn ms1_capture_ms2_before_any_ms1_links_to_none() {
        // MS2 first, then an MS1, then a second MS2.
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
        // First MS2 had no preceding MS1 → None; second is after the MS1 → Some(0).
        assert_eq!(link.ms2_to_ms1, vec![None, Some(0)]);
    }

    /// With capture OFF (default), the new path is not engaged and the plain
    /// iterator still emits MS2-only with no MS1 side effects.
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
        assert!(link.ms1_peaks.is_empty(), "no MS1 captured when capture off");
        // ms2_to_ms1 still has one entry per emitted MS2 (all None).
        assert_eq!(link.ms2_to_ms1, vec![None]);
    }

    #[test]
    fn detect_instrument_qexactive_model_cv_param() {
        // No analyzer cvParam, but a Q Exactive instrument-model cvParam
        // appears at the top of the IC block.
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
}
