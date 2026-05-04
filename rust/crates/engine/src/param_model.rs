//! Loader for Java MS-GF+'s `.param` binary format. See
//! `docs/superpowers/2026-05-03-phase2-param-loader-design.md` and the
//! Java reference at `src/main/java/edu/ucsd/msjava/msscorer/NewRankScorer.java`
//! lines 197-425.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::io::Cursor;
use std::path::Path;

use byteorder::{BigEndian, ReadBytesExt};

use crate::activation::ActivationMethod;
use crate::enzyme::Enzyme;
use crate::instrument::InstrumentType;
use crate::protocol::Protocol;
use crate::tolerance::Tolerance;

#[derive(Debug, Clone)]
pub struct Param {
    pub version: i32,
    pub data_type: SpecDataType,
    pub mme: Tolerance,
    pub apply_deconvolution: bool,
    pub deconvolution_error_tolerance: f32,
    pub charge_hist: Vec<(i32, i32)>,
    pub min_charge: i32,
    pub max_charge: i32,
    pub num_segments: i32,
    pub partitions: Vec<Partition>,
    pub num_precursor_off: i32,
    pub precursor_off_map: HashMap<i32, Vec<PrecursorOffsetFrequency>>,
    pub frag_off_table: HashMap<Partition, Vec<FragmentOffsetFrequency>>,
    pub max_rank: i32,
    pub rank_dist_table: HashMap<Partition, HashMap<IonType, Vec<f32>>>,
    pub error_scaling_factor: i32,
    pub ion_err_dist_table: HashMap<Partition, Vec<f32>>,
    pub noise_err_dist_table: HashMap<Partition, Vec<f32>>,
    pub ion_existence_table: HashMap<Partition, Vec<f32>>,
}

impl Param {
    /// Parse a complete `.param` byte stream produced by Java's
    /// `DataOutputStream`. Errors on buffer underruns, unknown enum
    /// names, missing validation marker, or trailing bytes.
    pub fn load_from_bytes(bytes: &[u8]) -> Result<Self, ParamParseError> {
        let mut cursor = Cursor::new(bytes);
        let param = read_param(&mut cursor)?;

        let validation = cursor.read_i32::<BigEndian>()
            .map_err(|_| ParamParseError::UnexpectedEof {
                offset: cursor.position() as usize, needed: 4,
            })?;
        if validation != i32::MAX {
            return Err(ParamParseError::ValidationMarker { got: validation });
        }
        let unread = (bytes.len() as u64).saturating_sub(cursor.position()) as usize;
        if unread != 0 {
            return Err(ParamParseError::TrailingBytes { unread });
        }
        Ok(param)
    }

    pub fn load_from_file(path: &Path) -> Result<Self, ParamParseError> {
        let bytes = std::fs::read(path)?;
        Self::load_from_bytes(&bytes)
    }
}

fn read_param(cursor: &mut Cursor<&[u8]>) -> Result<Param, ParamParseError> {
    // -- Section 1: header --
    let version = read_i32(cursor)?;

    let len_act = read_i8_as_u8(cursor)?;
    let act_str = read_utf16be_string(cursor, len_act)?;
    let activation = ActivationMethod::from_name(&act_str)
        .ok_or_else(|| ParamParseError::BadEnum { kind: "ActivationMethod", value: act_str })?;

    let len_inst = read_i8_as_u8(cursor)?;
    let inst_str = read_utf16be_string(cursor, len_inst)?;
    let instrument = InstrumentType::from_name(&inst_str)
        .ok_or_else(|| ParamParseError::BadEnum { kind: "InstrumentType", value: inst_str })?;

    let len_enz = read_i8_as_u8(cursor)?;
    let enzyme = if len_enz == 0 {
        None
    } else {
        let enz_str = read_utf16be_string(cursor, len_enz)?;
        Some(Enzyme::from_name(&enz_str)
            .ok_or_else(|| ParamParseError::BadEnum { kind: "Enzyme", value: enz_str })?)
    };

    let len_prot = read_i8_as_u8(cursor)?;
    let protocol = if len_prot == 0 {
        Protocol::Automatic
    } else {
        let prot_str = read_utf16be_string(cursor, len_prot)?;
        Protocol::from_name(&prot_str)
            .ok_or_else(|| ParamParseError::BadEnum { kind: "Protocol", value: prot_str })?
    };

    let data_type = SpecDataType { activation, instrument, enzyme, protocol };

    // -- Section 2: tolerance --
    let is_tol_ppm = read_bool(cursor)?;
    let mme_val = read_f32(cursor)?;
    let mme = if is_tol_ppm { Tolerance::Ppm(mme_val as f64) } else { Tolerance::Da(mme_val as f64) };

    // -- Section 3: deconvolution --
    let apply_deconvolution = read_bool(cursor)?;
    let deconvolution_error_tolerance = read_f32(cursor)?;

    // -- Section 4: charge histogram --
    let size = read_i32(cursor)?;
    let mut charge_hist = Vec::with_capacity(size as usize);
    let mut min_charge = i32::MAX;
    let mut max_charge = i32::MIN;
    for _ in 0..size {
        let charge = read_i32(cursor)?;
        let num_specs = read_i32(cursor)?;
        if charge < min_charge { min_charge = charge; }
        if charge > max_charge { max_charge = charge; }
        charge_hist.push((charge, num_specs));
    }
    let (min_charge, max_charge) = if size == 0 { (0, 0) } else { (min_charge, max_charge) };

    // -- Section 5: partition info --
    let part_size = read_i32(cursor)?;
    let num_segments = read_i32(cursor)?;
    let mut partitions = Vec::with_capacity(part_size as usize);
    for _ in 0..part_size {
        let charge = read_i32(cursor)?;
        let parent_mass = read_f32(cursor)?;
        let seg_num = read_i32(cursor)?;
        partitions.push(Partition { charge, parent_mass, seg_num });
    }
    // Java uses TreeSet for partition ordering — sort to match.
    partitions.sort();

    // -- Section 6: precursor offset frequency --
    let num_precursor_off = read_i32(cursor)?;
    let mut precursor_off_map: HashMap<i32, Vec<PrecursorOffsetFrequency>> = HashMap::new();
    for _ in 0..num_precursor_off {
        let charge = read_i32(cursor)?;
        let reduced_charge = read_i32(cursor)?;
        let offset = read_f32(cursor)?;
        let is_tol_ppm = read_bool(cursor)?;
        let tol_val = read_f32(cursor)?;
        let frequency = read_f32(cursor)?;
        let tolerance = if is_tol_ppm {
            Tolerance::Ppm(tol_val as f64)
        } else {
            Tolerance::Da(tol_val as f64)
        };
        precursor_off_map.entry(charge).or_default().push(PrecursorOffsetFrequency {
            reduced_charge, offset, tolerance, frequency,
        });
    }

    // -- Section 7: fragment offset frequency (per partition, in TreeSet order) --
    let mut frag_off_table: HashMap<Partition, Vec<FragmentOffsetFrequency>> = HashMap::new();
    for &partition in &partitions {
        let size = read_i32(cursor)?;
        let mut frags = Vec::with_capacity(size as usize);
        for _ in 0..size {
            let is_prefix = read_bool(cursor)?;
            let charge = read_i32(cursor)?;
            let offset = read_f32(cursor)?;
            let frequency = read_f32(cursor)?;
            let ion_type = if is_prefix {
                IonType::Prefix { charge, offset_bits: offset.to_bits() }
            } else {
                IonType::Suffix { charge, offset_bits: offset.to_bits() }
            };
            frags.push(FragmentOffsetFrequency { ion_type, frequency });
        }
        frag_off_table.insert(partition, frags);
    }

    // -- Section 8: rank distributions (per partition × per ion type incl. NOISE) --
    let max_rank = read_i32(cursor)?;
    let mut rank_dist_table: HashMap<Partition, HashMap<IonType, Vec<f32>>> = HashMap::new();
    for &partition in &partitions {
        let frag_list = frag_off_table.get(&partition);
        // Java skips partitions with no ion types; mirror that.
        if frag_list.map_or(true, |v| v.is_empty()) {
            continue;
        }
        let mut table: HashMap<IonType, Vec<f32>> = HashMap::new();
        let mut ion_types: Vec<IonType> = frag_list.unwrap().iter().map(|f| f.ion_type).collect();
        ion_types.push(IonType::Noise);
        for ion in ion_types {
            let mut frequencies = Vec::with_capacity((max_rank + 1) as usize);
            for _ in 0..(max_rank + 1) {
                frequencies.push(read_f32(cursor)?);
            }
            table.insert(ion, frequencies);
        }
        rank_dist_table.insert(partition, table);
    }

    // -- Section 9: error distributions (conditional) --
    let error_scaling_factor = read_i32(cursor)?;
    let mut ion_err_dist_table: HashMap<Partition, Vec<f32>> = HashMap::new();
    let mut noise_err_dist_table: HashMap<Partition, Vec<f32>> = HashMap::new();
    let mut ion_existence_table: HashMap<Partition, Vec<f32>> = HashMap::new();
    if error_scaling_factor > 0 {
        let dist_len = (error_scaling_factor as usize) * 2 + 1;
        for &partition in &partitions {
            let mut ion_err = Vec::with_capacity(dist_len);
            for _ in 0..dist_len { ion_err.push(read_f32(cursor)?); }
            ion_err_dist_table.insert(partition, ion_err);

            let mut noise_err = Vec::with_capacity(dist_len);
            for _ in 0..dist_len { noise_err.push(read_f32(cursor)?); }
            noise_err_dist_table.insert(partition, noise_err);

            let mut ion_ex = Vec::with_capacity(4);
            for _ in 0..4 { ion_ex.push(read_f32(cursor)?); }
            ion_existence_table.insert(partition, ion_ex);
        }
    }

    Ok(Param {
        version,
        data_type,
        mme,
        apply_deconvolution,
        deconvolution_error_tolerance,
        charge_hist,
        min_charge,
        max_charge,
        num_segments,
        partitions,
        num_precursor_off,
        precursor_off_map,
        frag_off_table,
        max_rank,
        rank_dist_table,
        error_scaling_factor,
        ion_err_dist_table,
        noise_err_dist_table,
        ion_existence_table,
    })
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SpecDataType {
    pub activation: ActivationMethod,
    pub instrument: InstrumentType,
    pub enzyme: Option<Enzyme>,
    pub protocol: Protocol,
}

#[derive(Debug, Clone, Copy)]
pub struct Partition {
    pub charge: i32,
    pub parent_mass: f32,
    pub seg_num: i32,
}

impl PartialEq for Partition {
    fn eq(&self, other: &Self) -> bool {
        self.charge == other.charge
            && self.parent_mass.to_bits() == other.parent_mass.to_bits()
            && self.seg_num == other.seg_num
    }
}

impl Eq for Partition {}

impl Hash for Partition {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.charge.hash(state);
        self.parent_mass.to_bits().hash(state);
        self.seg_num.hash(state);
    }
}

impl Ord for Partition {
    fn cmp(&self, other: &Self) -> Ordering {
        self.charge.cmp(&other.charge)
            .then_with(|| self.parent_mass.to_bits().cmp(&other.parent_mass.to_bits()))
            .then_with(|| self.seg_num.cmp(&other.seg_num))
    }
}

impl PartialOrd for Partition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IonType {
    /// Java: `IonType.PrefixIon`. `offset_bits` is `f32::to_bits` so the
    /// type can derive Eq/Hash; recover the float via `offset()`.
    Prefix { charge: i32, offset_bits: u32 },
    Suffix { charge: i32, offset_bits: u32 },
    Noise,
}

impl IonType {
    pub fn offset(&self) -> Option<f32> {
        match self {
            IonType::Prefix { offset_bits, .. } | IonType::Suffix { offset_bits, .. } => {
                Some(f32::from_bits(*offset_bits))
            }
            IonType::Noise => None,
        }
    }

    pub fn charge(&self) -> Option<i32> {
        match self {
            IonType::Prefix { charge, .. } | IonType::Suffix { charge, .. } => Some(*charge),
            IonType::Noise => None,
        }
    }

    pub fn is_prefix(&self) -> bool { matches!(self, IonType::Prefix { .. }) }
    pub fn is_suffix(&self) -> bool { matches!(self, IonType::Suffix { .. }) }
    pub fn is_noise(&self) -> bool { matches!(self, IonType::Noise) }
}

#[derive(Debug, Clone, Copy)]
pub struct PrecursorOffsetFrequency {
    pub reduced_charge: i32,
    pub offset: f32,
    pub tolerance: Tolerance,
    pub frequency: f32,
}

#[derive(Debug, Clone, Copy)]
pub struct FragmentOffsetFrequency {
    pub ion_type: IonType,
    pub frequency: f32,
}

#[derive(thiserror::Error, Debug)]
pub enum ParamParseError {
    #[error("I/O error reading param file: {source}")]
    Io { #[from] source: std::io::Error },
    #[error("buffer underrun at offset {offset}: needed {needed} more bytes")]
    UnexpectedEof { offset: usize, needed: usize },
    #[error("unknown {kind} {value:?} (Java enum failed to resolve)")]
    BadEnum { kind: &'static str, value: String },
    #[error("validation marker mismatch: got {got}, expected i32::MAX")]
    ValidationMarker { got: i32 },
    #[error("trailing bytes after validation marker: {unread} bytes left")]
    TrailingBytes { unread: usize },
    #[error("bad string length {got} (negative)")]
    BadStringLength { got: i8 },
    #[error("param reader path not yet implemented (Phase 2 stub)")]
    Unimplemented,
}

/// Read a UTF-16BE string of the given length (in 2-byte code units).
/// Length 0 → empty string. Non-ASCII code units are rejected.
fn read_utf16be_string(cursor: &mut Cursor<&[u8]>, len: u8) -> Result<String, ParamParseError> {
    let mut buf = String::with_capacity(len as usize);
    for _ in 0..len {
        let pos = cursor.position() as usize;
        let hi = cursor.read_u8()
            .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 1 })?;
        let lo = cursor.read_u8()
            .map_err(|_| ParamParseError::UnexpectedEof { offset: pos + 1, needed: 1 })?;
        let code_unit = ((hi as u16) << 8) | (lo as u16);
        if code_unit > 0x7F {
            return Err(ParamParseError::BadEnum {
                kind: "string",
                value: format!("non-ASCII u+{:04X}", code_unit),
            });
        }
        buf.push(code_unit as u8 as char);
    }
    Ok(buf)
}

// --- low-level read helpers ---

fn read_i32(cursor: &mut Cursor<&[u8]>) -> Result<i32, ParamParseError> {
    let pos = cursor.position() as usize;
    cursor.read_i32::<BigEndian>()
        .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 4 })
}

fn read_f32(cursor: &mut Cursor<&[u8]>) -> Result<f32, ParamParseError> {
    let pos = cursor.position() as usize;
    cursor.read_f32::<BigEndian>()
        .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 4 })
}

fn read_bool(cursor: &mut Cursor<&[u8]>) -> Result<bool, ParamParseError> {
    let pos = cursor.position() as usize;
    let b = cursor.read_u8()
        .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 1 })?;
    Ok(b != 0)
}

/// Read a single signed byte as the length prefix for a UTF-16BE string.
/// Java's `readByte` returns `i8`; values < 0 are illegal here.
fn read_i8_as_u8(cursor: &mut Cursor<&[u8]>) -> Result<u8, ParamParseError> {
    let pos = cursor.position() as usize;
    let b = cursor.read_i8()
        .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 1 })?;
    if b < 0 {
        return Err(ParamParseError::BadStringLength { got: b });
    }
    Ok(b as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_eq_via_to_bits() {
        let a = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let b = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        assert_eq!(a, b);
        let c = Partition { charge: 2, parent_mass: 1000.0001, seg_num: 0 };
        assert_ne!(a, c);
    }

    #[test]
    fn partition_ord_lex_order() {
        let a = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let b = Partition { charge: 2, parent_mass: 1000.0, seg_num: 1 };
        let c = Partition { charge: 3, parent_mass: 500.0,  seg_num: 0 };
        assert!(a < b);
        assert!(b < c);
    }

    #[test]
    fn partition_hash_consistent_with_eq() {
        use std::collections::HashSet;
        let a = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let b = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let set: HashSet<_> = [a, b].into_iter().collect();
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn ion_type_helpers() {
        let p = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let s = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let n = IonType::Noise;
        assert!(p.is_prefix());  assert!(!p.is_suffix()); assert!(!p.is_noise());
        assert!(!s.is_prefix()); assert!(s.is_suffix());  assert!(!s.is_noise());
        assert!(!n.is_prefix()); assert!(!n.is_suffix()); assert!(n.is_noise());
        assert_eq!(p.charge(), Some(1));
        assert_eq!(n.charge(), None);
    }

    #[test]
    fn ion_type_offset_round_trip() {
        let i = IonType::Prefix { charge: 2, offset_bits: 1.5_f32.to_bits() };
        assert_eq!(i.offset(), Some(1.5));
    }

    /// Build a minimal `.param`-style byte buffer that exercises sections
    /// 1-4 (header + tolerance + deconvolution + charge histogram).
    /// Tasks 7-9 extend this fixture as their tests are added.
    fn buf_sections_1_to_4() -> Vec<u8> {
        let mut b = Vec::new();
        // version
        b.extend(&10001_i32.to_be_bytes());
        // activation method "CID" — len 3, then 3 UTF-16BE chars
        b.push(3);
        for c in b"CID" { b.push(0); b.push(*c); }
        // instrument type "LowRes" — len 6
        b.push(6);
        for c in b"LowRes" { b.push(0); b.push(*c); }
        // enzyme "Tryp" — len 4 (Java's short name for Trypsin)
        b.push(4);
        for c in b"Tryp" { b.push(0); b.push(*c); }
        // protocol "Standard" — len 8
        b.push(8);
        for c in b"Standard" { b.push(0); b.push(*c); }
        // tolerance: is_ppm=true, mmeVal=20.0
        b.push(1);
        b.extend(&20.0_f32.to_be_bytes());
        // deconvolution: apply=false, errTol=0.5
        b.push(0);
        b.extend(&0.5_f32.to_be_bytes());
        // charge histogram: size=2, then 2 × (charge, num_specs)
        b.extend(&2_i32.to_be_bytes());
        b.extend(&2_i32.to_be_bytes()); b.extend(&100_i32.to_be_bytes());
        b.extend(&3_i32.to_be_bytes()); b.extend(&50_i32.to_be_bytes());
        b
    }

    #[test]
    fn reader_header_through_charge_hist() {
        // Append zero-content stubs for sections 5-9 + validation marker
        let mut b = buf_sections_1_to_4();
        b.extend(&0_i32.to_be_bytes()); b.extend(&1_i32.to_be_bytes());  // partition: size=0, num_segments=1
        b.extend(&0_i32.to_be_bytes());  // precursor OFF: size=0
        // fragment OFF: zero partitions => zero iterations (no bytes)
        b.extend(&0_i32.to_be_bytes());  // max_rank
        b.extend(&0_i32.to_be_bytes());  // error_scaling_factor=0
        b.extend(&i32::MAX.to_be_bytes());  // validation

        let param = Param::load_from_bytes(&b).unwrap();
        assert_eq!(param.version, 10001);
        assert_eq!(param.data_type.activation, ActivationMethod::CID);
        assert_eq!(param.data_type.instrument, InstrumentType::LowRes);
        assert_eq!(param.data_type.enzyme, Some(Enzyme::Trypsin));
        assert_eq!(param.data_type.protocol, Protocol::Standard);
        match param.mme {
            Tolerance::Ppm(v) => assert_eq!(v, 20.0),
            _ => panic!("expected Ppm"),
        }
        assert!(!param.apply_deconvolution);
        assert_eq!(param.deconvolution_error_tolerance, 0.5);
        assert_eq!(param.charge_hist.len(), 2);
        assert_eq!(param.min_charge, 2);
        assert_eq!(param.max_charge, 3);
    }
}
