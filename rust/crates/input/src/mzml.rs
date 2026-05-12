//! Streaming mzML reader. Event-driven via quick-xml; no serde tree.
//!
//! By default only MS2 spectra are emitted (ms level == 2). The parser
//! decodes base64 peak arrays (32-bit or 64-bit float, little-endian)
//! with optional zlib compression and zips (m/z, intensity) pairs into
//! `Vec<(f64, f32)>` sorted ascending by m/z.

use std::io::BufRead;

use base64::{engine::general_purpose::STANDARD, Engine as _};
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::ZlibDecoder;
use quick_xml::{events::Event, Reader};

use model::Spectrum;

// ── CV accessions we care about ─────────────────────────────────────────────

const CV_MS_LEVEL: &str = "MS:1000511";
const CV_SCAN_TIME: &str = "MS:1000016";
const CV_SELECTED_ION_MZ: &str = "MS:1000744";
/// Older mzML files sometimes use plain m/z accession in selectedIon.
const CV_MZ_PLAIN: &str = "MS:1000040";
const CV_CHARGE_STATE: &str = "MS:1000041";
const CV_PEAK_INTENSITY: &str = "MS:1000042";
const CV_MZ_ARRAY: &str = "MS:1000514";
const CV_INTENSITY_ARRAY: &str = "MS:1000515";
const CV_64BIT: &str = "MS:1000523";
const CV_32BIT: &str = "MS:1000521";
const CV_ZLIB: &str = "MS:1000574";

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
    precursor_charge: Option<i32>,
    precursor_intensity: Option<f32>,
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
        }
    }

    /// Widen or narrow the ms-level filter (e.g. `with_ms_level_range(1, 2)`
    /// emits both MS1 and MS2).
    pub fn with_ms_level_range(mut self, min: u32, max: u32) -> Self {
        self.ms_level_min = min;
        self.ms_level_max = max;
        self
    }

    // ── Build a Spectrum from a completed SpectrumBuilder ────────────────────

    fn finish_spectrum(&self, sb: SpectrumBuilder) -> Result<Option<Spectrum>, MzMLParseError> {
        let level = sb.ms_level.unwrap_or(0);
        if level < self.ms_level_min || level > self.ms_level_max {
            return Ok(None);
        }

        let precursor_mz = match sb.precursor_mz {
            Some(v) => v,
            // MS2 without a precursor m/z: skip rather than error.
            None => return Ok(None),
        };

        let mz_vals = sb.mz_array.unwrap_or_default();
        let int_vals = sb.intensity_array.unwrap_or_default();

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
            .collect();

        // Enforce ascending-by-m/z invariant required by downstream consumers.
        if !peaks.windows(2).all(|w| w[0].0 <= w[1].0) {
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        }

        let scan = extract_scan_from_id(&sb.id);

        Ok(Some(Spectrum {
            title: sb.id,
            precursor_mz,
            precursor_charge: sb.precursor_charge,
            precursor_intensity: sb.precursor_intensity,
            rt_seconds: sb.rt_seconds,
            scan,
            peaks,
        }))
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
                        sb.precursor_charge = Some(z);
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

                // Self-closing elements — mostly cvParam.
                Event::Empty(ref e) => {
                    let tag = e.local_name().as_ref().to_owned();
                    if tag == b"cvParam" {
                        // Extract info before any &mut self call.
                        if let Some(cv) = CvParamInfo::from_bytes_start(e) {
                            self.apply_cv_param(cv);
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
}

impl<R: BufRead> Iterator for MzMLReader<R> {
    type Item = Result<Spectrum, MzMLParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.pump() {
            Ok(Some(s)) => Some(Ok(s)),
            Ok(None) => None,
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
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
            .join("../../../../src/test/resources/tiny.pwiz.mzML");

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

    #[test]
    fn extract_scan_missing() {
        assert_eq!(extract_scan_from_id("spectrum=1"), None);
    }
}
