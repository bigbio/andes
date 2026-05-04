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

/// Tasks 6-9 fill this in section by section. Phase 2/Task 5 leaves a
/// stub that returns `Unimplemented` so the type-level wiring can be
/// tested without the full reader.
fn read_param(_cursor: &mut Cursor<&[u8]>) -> Result<Param, ParamParseError> {
    Err(ParamParseError::Unimplemented)
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
#[allow(dead_code)]
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
}
