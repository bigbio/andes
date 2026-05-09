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

use model::activation::ActivationMethod;
use model::enzyme::Enzyme;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;

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
    /// Pre-filtered ion-type list per partition (Noise excluded), populated
    /// at load time. Used by `ion_types_for_partition_slice` to avoid
    /// per-call Vec allocation in the GF DP hot path.
    pub(crate) partition_ion_types_cache: HashMap<Partition, Vec<IonType>>,
}

impl Param {
    /// Find the partition matching `(charge, parent_mass, seg_num)` via
    /// the same lookup Java uses (`partitionSet.floor(target)`):
    /// returns the largest partition ≤ target by lex order on
    /// (charge, parent_mass.to_bits(), seg_num).
    ///
    /// Falls back gracefully:
    /// - If no partition matches the requested charge: use the smallest
    ///   charge available with the requested mass + segment.
    /// - If charge > all available: use the largest available charge.
    pub fn find_partition(&self, charge: i32, parent_mass: f32, seg_num: i32) -> Option<Partition> {
        if self.partitions.is_empty() {
            return None;
        }

        // Build the target partition for the floor lookup.
        let target = Partition { charge, parent_mass, seg_num };

        // partitions is already sorted (Phase 2 guarantee). Find the largest
        // partition <= target.
        // Use binary search via partition_point.
        let pos = self.partitions.partition_point(|p| p <= &target);
        if pos > 0 {
            // partitions[pos - 1] is the largest <= target.
            let candidate = self.partitions[pos - 1];
            if candidate.charge == charge {
                return Some(candidate);
            }
            // Floor returned a partition with smaller charge; Java's logic
            // says: if no exact-charge match, find smallest available charge,
            // then floor on (smallest_charge, parent_mass, seg_num).
        }

        // Fall back: find smallest charge in partitions, retry.
        let min_charge = self.partitions.iter().map(|p| p.charge).min()?;
        let max_charge = self.partitions.iter().map(|p| p.charge).max()?;
        let fallback_charge = if charge < min_charge {
            min_charge
        } else if charge > max_charge {
            max_charge
        } else {
            // charge is in range but had no exact match — already handled above.
            return self.partitions.last().copied();
        };
        let fallback_target = Partition { charge: fallback_charge, parent_mass, seg_num };
        let fallback_pos = self.partitions.partition_point(|p| p <= &fallback_target);
        if fallback_pos > 0 {
            let candidate = self.partitions[fallback_pos - 1];
            if candidate.charge == fallback_charge {
                return Some(candidate);
            }
        }
        // Last resort: just return any partition with the fallback charge.
        self.partitions.iter().find(|p| p.charge == fallback_charge).copied()
    }

    /// Compute the segment number for a peak m/z relative to the peptide's
    /// parent mass. Mirrors Java's `getSegmentNum`.
    pub fn segment_num_for(&self, peak_mz: f64, parent_mass: f64) -> i32 {
        if parent_mass <= 0.0 || self.num_segments <= 0 {
            return 0;
        }
        let seg = (peak_mz / parent_mass * self.num_segments as f64) as i32;
        seg.min(self.num_segments - 1).max(0)
    }

    /// Alias for `segment_num_for` matching the name used by the GF DP code
    /// (`param.segment_num(theo_mz, parent_mass)`).
    #[inline]
    pub fn segment_num(&self, peak_mz: f64, parent_mass: f64) -> usize {
        self.segment_num_for(peak_mz, parent_mass) as usize
    }

    /// Mirrors Java NewRankScorer.getIonTypes(charge, parentMass, seg) — ion-type membership comes from frag_off_table, not rank_dist_table.
    ///
    /// Collect the unique ion types (Prefix and Suffix, not Noise) whose
    /// partition has `seg_num == seg`. Derived from `frag_off_table` keys.
    ///
    /// Returned in stable insertion order; duplicates suppressed.
    pub fn ion_types_for_segment(&self, seg: usize) -> Vec<IonType> {
        let mut seen: std::collections::HashSet<IonType> = std::collections::HashSet::new();
        let mut out: Vec<IonType> = Vec::new();
        for (partition, frag_list) in &self.frag_off_table {
            if partition.seg_num as usize != seg {
                continue;
            }
            for fof in frag_list {
                let ion = fof.ion_type;
                if matches!(ion, IonType::Noise) {
                    continue;
                }
                if seen.insert(ion) {
                    out.push(ion);
                }
            }
        }
        out
    }

    /// Find the partition for `(charge, parent_mass, seg_num)` using the
    /// floor-lookup semantics of `find_partition`. Returns a synthetic
    /// partition if none is found (so callers don't need to unwrap).
    pub fn partition_for(&self, charge: u8, parent_mass: f64, seg_num: usize) -> Partition {
        self.find_partition(charge as i32, parent_mass as f32, seg_num as i32)
            .unwrap_or(Partition {
                charge: charge as i32,
                parent_mass: parent_mass as f32,
                seg_num: seg_num as i32,
            })
    }

    /// Ion types for the SPECIFIC partition `(charge, parent_mass, seg)`.
    /// Mirrors Java `NewRankScorer.getIonTypes(charge, parentMass, segNum)`,
    /// which selects the partition's ion list from `fragOFFTable` rather
    /// than the segment-wide union returned by `ion_types_for_segment`.
    /// Used in the per-node scoring path so that Rust enumerates the
    /// same ion set as Java for a given spectrum.
    pub fn ion_types_for_partition(&self, charge: u8, parent_mass: f64, seg: usize) -> Vec<IonType> {
        // Compat shim — callers in hot paths should use
        // `ion_types_for_partition_slice` to avoid the allocation.
        self.ion_types_for_partition_slice(charge, parent_mass, seg).to_vec()
    }

    /// Slice-borrowing version of `ion_types_for_partition`. Reads from the
    /// pre-filtered `partition_ion_types_cache` populated at param-load time.
    /// Zero allocations per call. Used by the GF DP hot path.
    pub fn ion_types_for_partition_slice(&self, charge: u8, parent_mass: f64, seg: usize) -> &[IonType] {
        let part = self.partition_for(charge, parent_mass, seg);
        self.partition_ion_types_cache
            .get(&part)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Parse a complete `.param` byte stream produced by Java's
    /// `DataOutputStream`. Errors on buffer underruns, unknown enum
    /// names, missing validation marker, or trailing bytes.
    pub fn load_from_bytes(bytes: &[u8]) -> Result<Self> {
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

    pub fn load_from_file(path: &Path) -> Result<Self> {
        let bytes = std::fs::read(path)?;
        Self::load_from_bytes(&bytes)
    }
}

fn read_param(cursor: &mut Cursor<&[u8]>) -> Result<Param> {
    // -- Section 1: header --
    let version = read_i32(cursor)?;

    let len_act = read_i8_as_u8(cursor)?;
    let act_str = read_utf16be_string(cursor, len_act)?;
    let activation = ActivationMethod::from_name(&act_str)
        .ok_or(ParamParseError::BadEnum { kind: "ActivationMethod", value: act_str })?;

    let len_inst = read_i8_as_u8(cursor)?;
    let inst_str = read_utf16be_string(cursor, len_inst)?;
    let instrument = InstrumentType::from_name(&inst_str)
        .ok_or(ParamParseError::BadEnum { kind: "InstrumentType", value: inst_str })?;

    let len_enz = read_i8_as_u8(cursor)?;
    let enzyme = if len_enz == 0 {
        None
    } else {
        let enz_str = read_utf16be_string(cursor, len_enz)?;
        Some(Enzyme::from_name(&enz_str)
            .ok_or(ParamParseError::BadEnum { kind: "Enzyme", value: enz_str })?)
    };

    let len_prot = read_i8_as_u8(cursor)?;
    let protocol = if len_prot == 0 {
        Protocol::Automatic
    } else {
        let prot_str = read_utf16be_string(cursor, len_prot)?;
        Protocol::from_name(&prot_str)
            .ok_or(ParamParseError::BadEnum { kind: "Protocol", value: prot_str })?
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

    // Pre-build per-partition ion-type cache (Noise excluded), so the GF
    // DP hot path can borrow a slice instead of allocating a Vec per call.
    let mut partition_ion_types_cache: HashMap<Partition, Vec<IonType>> = HashMap::new();
    for (&part, frag_list) in &frag_off_table {
        let mut ions: Vec<IonType> = Vec::with_capacity(frag_list.len());
        for fof in frag_list {
            if !matches!(fof.ion_type, IonType::Noise) {
                ions.push(fof.ion_type);
            }
        }
        partition_ion_types_cache.insert(part, ions);
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
        partition_ion_types_cache,
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
        // Mirror Java `Partition.compareTo`: charge → segIndex → parentMass.
        // (Bug fix 2026-05-09: previously charge → parent_mass → seg_num,
        // which produced different floor-lookup results — `find_partition`
        // for seg=0 queries returned a seg=1 partition with the same
        // parent_mass tier, looking up the WRONG rank distribution table.)
        self.charge.cmp(&other.charge)
            .then_with(|| self.seg_num.cmp(&other.seg_num))
            .then_with(|| self.parent_mass.to_bits().cmp(&other.parent_mass.to_bits()))
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

    /// Compute the predicted m/z for this ion type given a node mass (in Da).
    ///
    /// For both `Prefix` and `Suffix`, the formula mirrors Java's
    /// `IonType.PrefixIon.getMz(mass)` and `IonType.SuffixIon.getMz(mass)`:
    ///   `mz = (node_mass + offset + charge * PROTON) / charge`
    ///
    /// For `Noise`, returns 0.0.
    ///
    /// Note: for suffix ions, `node_mass` is the suffix mass already in
    /// MS-GF+ convention; callers are responsible for supplying the correct
    /// directional mass.
    pub fn mz(&self, node_mass: f64) -> f64 {
        match self {
            IonType::Prefix { charge, offset_bits } | IonType::Suffix { charge, offset_bits } => {
                let offset = f32::from_bits(*offset_bits) as f64;
                let c = *charge as f64;
                (node_mass + offset + c * model::mass::PROTON) / c
            }
            IonType::Noise => 0.0,
        }
    }

    /// Inverse of `mz`: given an observed peak m/z, recover the node mass.
    ///
    /// For `Prefix { charge, offset }`: `mass = mz * charge - charge * PROTON - offset`.
    /// For `Suffix`: same formula.
    /// For `Noise`: returns 0.0.
    pub fn mass_from_mz(&self, mz: f64) -> f64 {
        match self {
            IonType::Prefix { charge, offset_bits } | IonType::Suffix { charge, offset_bits } => {
                let offset = f32::from_bits(*offset_bits) as f64;
                let c = *charge as f64;
                mz * c - c * model::mass::PROTON - offset
            }
            IonType::Noise => 0.0,
        }
    }
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

/// Module-local Result alias to reduce signature noise.
pub type Result<T> = std::result::Result<T, ParamParseError>;

/// Read a UTF-16BE string of the given length (in 2-byte code units).
/// Length 0 → empty string. Non-ASCII code units are rejected.
fn read_utf16be_string(cursor: &mut Cursor<&[u8]>, len: u8) -> Result<String> {
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

fn read_i32(cursor: &mut Cursor<&[u8]>) -> Result<i32> {
    let pos = cursor.position() as usize;
    cursor.read_i32::<BigEndian>()
        .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 4 })
}

fn read_f32(cursor: &mut Cursor<&[u8]>) -> Result<f32> {
    let pos = cursor.position() as usize;
    cursor.read_f32::<BigEndian>()
        .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 4 })
}

fn read_bool(cursor: &mut Cursor<&[u8]>) -> Result<bool> {
    let pos = cursor.position() as usize;
    let b = cursor.read_u8()
        .map_err(|_| ParamParseError::UnexpectedEof { offset: pos, needed: 1 })?;
    Ok(b != 0)
}

/// Read a single signed byte as the length prefix for a UTF-16BE string.
/// Java's `readByte` returns `i8`; values < 0 are illegal here.
fn read_i8_as_u8(cursor: &mut Cursor<&[u8]>) -> Result<u8> {
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

    #[test]
    fn reader_partitions_and_precursor_off() {
        let mut b = buf_sections_1_to_4();
        // Partition info: size=2, num_segments=4
        b.extend(&2_i32.to_be_bytes()); b.extend(&4_i32.to_be_bytes());
        // Partition 1: charge=2, parentMass=500.0, segNum=0
        b.extend(&2_i32.to_be_bytes());
        b.extend(&500.0_f32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        // Partition 2: charge=2, parentMass=1500.0, segNum=1
        b.extend(&2_i32.to_be_bytes());
        b.extend(&1500.0_f32.to_be_bytes());
        b.extend(&1_i32.to_be_bytes());
        // Precursor OFF: size=1
        b.extend(&1_i32.to_be_bytes());
        // entry: charge=2, reducedCharge=1, offset=0.0, isTolPpm=false, tolVal=0.5, freq=0.8
        b.extend(&2_i32.to_be_bytes());
        b.extend(&1_i32.to_be_bytes());
        b.extend(&0.0_f32.to_be_bytes());
        b.push(0);
        b.extend(&0.5_f32.to_be_bytes());
        b.extend(&0.8_f32.to_be_bytes());
        // Fragment OFF for both partitions: each empty (size=0)
        b.extend(&0_i32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        // Rank distributions: max_rank=0; partitions skip because frag_list empty
        b.extend(&0_i32.to_be_bytes());
        // Error distributions: error_scaling_factor=0
        b.extend(&0_i32.to_be_bytes());
        // Validation
        b.extend(&i32::MAX.to_be_bytes());

        let p = Param::load_from_bytes(&b).unwrap();
        assert_eq!(p.partitions.len(), 2);
        // Sorted by (charge, parent_mass.to_bits(), seg_num)
        assert_eq!(p.partitions[0].seg_num, 0);
        assert_eq!(p.partitions[1].seg_num, 1);
        assert_eq!(p.num_segments, 4);
        assert_eq!(p.num_precursor_off, 1);

        let off_list = p.precursor_off_map.get(&2).unwrap();
        assert_eq!(off_list.len(), 1);
        assert_eq!(off_list[0].reduced_charge, 1);
        match off_list[0].tolerance {
            Tolerance::Da(v) => assert_eq!(v, 0.5),
            _ => panic!("expected Da"),
        }
    }

    #[test]
    fn reader_fragment_off_and_rank_dist() {
        let mut b = buf_sections_1_to_4();
        // Partition info: 1 partition, num_segments=1
        b.extend(&1_i32.to_be_bytes()); b.extend(&1_i32.to_be_bytes());
        b.extend(&2_i32.to_be_bytes());
        b.extend(&1000.0_f32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        // Precursor OFF: 0 entries
        b.extend(&0_i32.to_be_bytes());
        // Fragment OFF for partition 1: size=2 (1 prefix + 1 suffix)
        b.extend(&2_i32.to_be_bytes());
        // Frag entry 1: prefix, charge=1, offset=1.00782, freq=0.7
        b.push(1);
        b.extend(&1_i32.to_be_bytes());
        b.extend(&1.00782_f32.to_be_bytes());
        b.extend(&0.7_f32.to_be_bytes());
        // Frag entry 2: suffix, charge=1, offset=18.01057, freq=0.6
        b.push(0);
        b.extend(&1_i32.to_be_bytes());
        b.extend(&18.01057_f32.to_be_bytes());
        b.extend(&0.6_f32.to_be_bytes());
        // Rank distributions: max_rank=2, so 3 floats per ion type.
        b.extend(&2_i32.to_be_bytes());
        // 3 ion types: prefix, suffix, NOISE; 3 floats each
        for &v in &[0.5_f32, 0.4, 0.3] { b.extend(&v.to_be_bytes()); }
        for &v in &[0.45_f32, 0.35, 0.25] { b.extend(&v.to_be_bytes()); }
        for &v in &[0.05_f32, 0.05, 0.05] { b.extend(&v.to_be_bytes()); }
        // Error distributions: error_scaling_factor=0
        b.extend(&0_i32.to_be_bytes());
        // Validation
        b.extend(&i32::MAX.to_be_bytes());

        let p = Param::load_from_bytes(&b).unwrap();
        assert_eq!(p.partitions.len(), 1);
        let part = p.partitions[0];
        let frags = p.frag_off_table.get(&part).unwrap();
        assert_eq!(frags.len(), 2);
        assert!(frags[0].ion_type.is_prefix());
        assert!(frags[1].ion_type.is_suffix());
        assert_eq!(p.max_rank, 2);
        let rank_table = p.rank_dist_table.get(&part).unwrap();
        // 2 ion types + NOISE = 3 entries
        assert_eq!(rank_table.len(), 3);
        for freqs in rank_table.values() {
            assert_eq!(freqs.len(), 3);
        }
    }

    #[test]
    fn reader_error_distributions() {
        let mut b = buf_sections_1_to_4();
        // 1 partition
        b.extend(&1_i32.to_be_bytes()); b.extend(&1_i32.to_be_bytes());
        b.extend(&2_i32.to_be_bytes());
        b.extend(&1000.0_f32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        // 0 precursor OFF
        b.extend(&0_i32.to_be_bytes());
        // Fragment OFF: 1 prefix entry
        b.extend(&1_i32.to_be_bytes());
        b.push(1);
        b.extend(&1_i32.to_be_bytes());
        b.extend(&1.0_f32.to_be_bytes());
        b.extend(&0.5_f32.to_be_bytes());
        // Rank dist max_rank=0; 2 ion types (prefix + NOISE) × 1 float each
        b.extend(&0_i32.to_be_bytes());
        b.extend(&0.5_f32.to_be_bytes());
        b.extend(&0.1_f32.to_be_bytes());
        // Error distributions: error_scaling_factor=2 → 2*2+1 = 5 floats per dist
        b.extend(&2_i32.to_be_bytes());
        // ionErr: 5 floats
        for v in [0.1_f32, 0.2, 0.4, 0.2, 0.1] { b.extend(&v.to_be_bytes()); }
        // noiseErr: 5 floats
        for v in [0.05_f32, 0.10, 0.70, 0.10, 0.05] { b.extend(&v.to_be_bytes()); }
        // ionExistence: 4 floats
        for v in [0.9_f32, 0.8, 0.7, 0.6] { b.extend(&v.to_be_bytes()); }
        // Validation
        b.extend(&i32::MAX.to_be_bytes());

        let p = Param::load_from_bytes(&b).unwrap();
        assert_eq!(p.error_scaling_factor, 2);
        let part = p.partitions[0];
        assert_eq!(p.ion_err_dist_table.get(&part).unwrap().len(), 5);
        assert_eq!(p.noise_err_dist_table.get(&part).unwrap().len(), 5);
        assert_eq!(p.ion_existence_table.get(&part).unwrap().len(), 4);
    }

    #[test]
    fn reader_rejects_bad_validation_marker() {
        let mut b = buf_sections_1_to_4();
        b.extend(&0_i32.to_be_bytes()); b.extend(&1_i32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        // BAD validation marker
        b.extend(&0_i32.to_be_bytes());

        let err = Param::load_from_bytes(&b).unwrap_err();
        match err {
            ParamParseError::ValidationMarker { got } => assert_eq!(got, 0),
            other => panic!("expected ValidationMarker, got {:?}", other),
        }
    }

    #[test]
    fn reader_rejects_trailing_bytes() {
        let mut b = buf_sections_1_to_4();
        b.extend(&0_i32.to_be_bytes()); b.extend(&1_i32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        b.extend(&0_i32.to_be_bytes());
        b.extend(&i32::MAX.to_be_bytes());
        // Trailing junk
        b.extend(&[1u8, 2, 3, 4]);

        let err = Param::load_from_bytes(&b).unwrap_err();
        match err {
            ParamParseError::TrailingBytes { unread } => assert_eq!(unread, 4),
            other => panic!("expected TrailingBytes, got {:?}", other),
        }
    }

    #[test]
    fn reader_rejects_unknown_activation() {
        let mut b = Vec::new();
        b.extend(&10001_i32.to_be_bytes());
        // activation: "GARBAGE"
        b.push(7);
        for c in b"GARBAGE" { b.push(0); b.push(*c); }
        let err = Param::load_from_bytes(&b).unwrap_err();
        match err {
            ParamParseError::BadEnum { kind, value } => {
                assert_eq!(kind, "ActivationMethod");
                assert_eq!(value, "GARBAGE");
            }
            other => panic!("expected BadEnum, got {:?}", other),
        }
    }

    fn make_param() -> Param {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use std::collections::HashMap;

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![],
            min_charge: 2,
            max_charge: 3,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 0,
            precursor_off_map: HashMap::new(),
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        }
    }

    #[test]
    fn find_partition_exact_charge_match() {
        let mut param = make_param();
        param.partitions = vec![
            Partition { charge: 2, parent_mass: 500.0, seg_num: 0 },
            Partition { charge: 2, parent_mass: 500.0, seg_num: 1 },
            Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 500.0, seg_num: 0 },
        ];
        // Sort matches the Phase 2 invariant.
        param.partitions.sort();

        // Target (2, 800.0, 0) → lex floor:
        // Sorted charge-2 partitions: (2,500.0,0), (2,500.0,1), (2,1000.0,0).
        // Largest ≤ (2,800.0,0): (2,500.0,1) because 500.0 < 800.0 and seg_num
        // doesn't matter once parent_mass is strictly less.
        let p = param.find_partition(2, 800.0, 0).expect("find");
        assert_eq!(p.charge, 2);
        assert_eq!(p.parent_mass, 500.0);
        assert_eq!(p.seg_num, 1);
    }

    #[test]
    fn find_partition_low_charge_fallback() {
        let mut param = make_param();
        param.partitions = vec![
            Partition { charge: 2, parent_mass: 500.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 500.0, seg_num: 0 },
        ];
        param.partitions.sort();

        // Target charge 1 (below all): falls back to smallest charge = 2.
        let p = param.find_partition(1, 500.0, 0).expect("find with fallback");
        assert_eq!(p.charge, 2);
    }

    #[test]
    fn find_partition_high_charge_fallback() {
        let mut param = make_param();
        param.partitions = vec![
            Partition { charge: 2, parent_mass: 500.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 500.0, seg_num: 0 },
        ];
        param.partitions.sort();

        // Target charge 5 (above all): falls back to largest = 3.
        let p = param.find_partition(5, 500.0, 0).expect("find with fallback");
        assert_eq!(p.charge, 3);
    }

    #[test]
    fn segment_num_clamps_to_max() {
        let mut param = make_param();
        param.num_segments = 3;
        // peak_mz / parent_mass × num_segments = floor calculation
        assert_eq!(param.segment_num_for(50.0, 100.0), 1);
        assert_eq!(param.segment_num_for(99.0, 100.0), 2);
        assert_eq!(param.segment_num_for(100.0, 100.0), 2);  // clamped
        assert_eq!(param.segment_num_for(120.0, 100.0), 2);  // clamped
    }

    #[test]
    fn ion_type_mz_prefix_charge1_offset0() {
        use model::mass::PROTON;
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let node_mass = 100.0_f64;
        let expected = (node_mass + 0.0 + PROTON) / 1.0;
        assert!((ion.mz(node_mass) - expected).abs() < 1e-9);
    }

    #[test]
    fn ion_type_mz_prefix_charge2() {
        use model::mass::PROTON;
        let ion = IonType::Prefix { charge: 2, offset_bits: 0.0_f32.to_bits() };
        let node_mass = 200.0_f64;
        let expected = (node_mass + 0.0 + 2.0 * PROTON) / 2.0;
        assert!((ion.mz(node_mass) - expected).abs() < 1e-9);
    }

    #[test]
    fn ion_type_mz_suffix_same_formula_as_prefix() {
        // Java's SuffixIon.getMz uses the same formula as PrefixIon.getMz.
        let offset = 18.01_f32;
        let prefix = IonType::Prefix { charge: 1, offset_bits: offset.to_bits() };
        let suffix = IonType::Suffix { charge: 1, offset_bits: offset.to_bits() };
        let node_mass = 150.0_f64;
        assert!((prefix.mz(node_mass) - suffix.mz(node_mass)).abs() < 1e-9);
    }

    #[test]
    fn ion_type_mz_noise_returns_zero() {
        assert_eq!(IonType::Noise.mz(100.0), 0.0);
    }

    #[test]
    fn ion_type_mass_from_mz_round_trips() {
        use model::mass::PROTON;
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let node_mass = 100.0_f64;
        let mz = ion.mz(node_mass);
        let recovered = ion.mass_from_mz(mz);
        assert!((recovered - node_mass).abs() < 1e-9,
            "round-trip failed: original {node_mass}, recovered {recovered}, mz={mz}, PROTON={PROTON}");
    }

    #[test]
    fn ion_types_for_segment_returns_unique() {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use std::collections::HashMap;

        let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let prefix = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let suffix = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() };

        // Populate frag_off_table (the source of truth for ion_types_for_segment,
        // mirroring Java NewRankScorer.getIonTypes which reads from fragOFFTable).
        let mut frag_off_table: HashMap<Partition, Vec<FragmentOffsetFrequency>> = HashMap::new();
        frag_off_table.insert(part, vec![
            FragmentOffsetFrequency { ion_type: prefix, frequency: 0.7 },
            FragmentOffsetFrequency { ion_type: suffix, frequency: 0.6 },
        ]);

        let param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Da(0.5),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: HashMap::new(),
            frag_off_table,
            max_rank: 2,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        };

        let seg0 = param.ion_types_for_segment(0);
        // Should return prefix and suffix (not noise), no duplicates.
        assert_eq!(seg0.len(), 2);
        assert!(seg0.iter().all(|i| !i.is_noise()));
        assert!(seg0.iter().any(|i| i.is_prefix()));
        assert!(seg0.iter().any(|i| i.is_suffix()));

        // Segment 1 has no partitions → empty.
        let seg1 = param.ion_types_for_segment(1);
        assert!(seg1.is_empty());
    }
}
