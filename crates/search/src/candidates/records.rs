use std::io::{self, Write};
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

use rayon::prelude::*;

use super::{candidate_has_nxst_sequon, enumerate_protein, Candidate};
use crate::decoy::is_decoy_accession_affix;
use crate::index::SearchIndex;
use crate::params::SearchParams;

pub const INDEX_RECORD_SIZE: usize = 20;

const INDEX_MAGIC: &[u8; 8] = b"ANDESCI\x00";
const INDEX_FORMAT_VERSION: u32 = 1;
const CANDIDATE_SET_REVISION: u32 = 2;
const INDEX_HEADER_SIZE: usize = INDEX_MAGIC.len() + 4;

pub mod flags {
    pub const IS_DECOY: u16 = 1 << 0;
    pub const IS_PROTEIN_N_TERM: u16 = 1 << 1;
    pub const IS_PROTEIN_C_TERM: u16 = 1 << 2;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexRecord {
    pub mass_milli: u64,
    pub protein_index: u32,
    pub start_offset: u32,
    pub length: u16,
    pub flags: u16,
}

impl IndexRecord {
    pub fn is_decoy(&self) -> bool {
        self.flags & flags::IS_DECOY != 0
    }

    pub fn to_bytes(&self) -> [u8; INDEX_RECORD_SIZE] {
        let mut buf = [0u8; INDEX_RECORD_SIZE];
        buf[0..8].copy_from_slice(&self.mass_milli.to_le_bytes());
        buf[8..12].copy_from_slice(&self.protein_index.to_le_bytes());
        buf[12..16].copy_from_slice(&self.start_offset.to_le_bytes());
        buf[16..18].copy_from_slice(&self.length.to_le_bytes());
        buf[18..20].copy_from_slice(&self.flags.to_le_bytes());
        buf
    }

    pub fn from_bytes(buf: &[u8; INDEX_RECORD_SIZE]) -> Self {
        let mass_milli = u64::from_le_bytes(buf[0..8].try_into().unwrap());
        let protein_index = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        let start_offset = u32::from_le_bytes(buf[12..16].try_into().unwrap());
        let length = u16::from_le_bytes(buf[16..18].try_into().unwrap());
        let flags = u16::from_le_bytes(buf[18..20].try_into().unwrap());
        Self {
            mass_milli,
            protein_index,
            start_offset,
            length,
            flags,
        }
    }
}

/// Proteins digested per parallel block of `base_peptide_records`: each block's records join the
/// result before the next block starts, so the records are held about once, not twice.
const DIGEST_BLOCK_PROTEINS: usize = 1024;

/// The unmodified candidates of `idx` as index records, in protein order: the order
/// `enumerate_candidates` yields them in. The proteins of each block are digested in parallel
/// and their records are appended in protein order.
pub fn base_peptide_records(
    idx: &SearchIndex,
    params: &SearchParams,
    decoy_prefix: &str,
) -> Vec<IndexRecord> {
    let base_params = SearchParams {
        max_variable_mods_per_peptide: 0,
        max_label_sites: 0,
        ..params.clone()
    };
    let suffix = idx.decoy_suffix.as_deref();
    let require_sequon = base_params.require_nxst_sequon;
    let mut records = Vec::new();
    for (block_idx, block) in idx.db.proteins.chunks(DIGEST_BLOCK_PROTEINS).enumerate() {
        let first = block_idx * DIGEST_BLOCK_PROTEINS;
        let block_records: Vec<IndexRecord> = block
            .par_iter()
            .enumerate()
            .flat_map_iter(|(i, protein)| {
                let is_decoy = is_decoy_accession_affix(&protein.accession, decoy_prefix, suffix);
                enumerate_protein(protein, first + i, is_decoy, &base_params)
                    .into_iter()
                    .filter(|c| !require_sequon || candidate_has_nxst_sequon(c))
                    .map(|c| index_record(&c))
            })
            .collect();
        records.extend(block_records);
    }
    records
}

fn index_record(c: &Candidate) -> IndexRecord {
    let mass_milli = (c.peptide.mass() * 1000.0).round() as u64;
    let mut f: u16 = 0;
    if c.is_decoy {
        f |= flags::IS_DECOY;
    }
    if c.is_protein_n_term {
        f |= flags::IS_PROTEIN_N_TERM;
    }
    if c.is_protein_c_term {
        f |= flags::IS_PROTEIN_C_TERM;
    }
    debug_assert!(
        c.protein_index <= u32::MAX as usize,
        "protein_index overflows u32"
    );
    debug_assert!(
        c.start_offset_in_protein <= u32::MAX as usize,
        "start_offset_in_protein overflows u32"
    );
    IndexRecord {
        mass_milli,
        protein_index: c.protein_index as u32,
        start_offset: c.start_offset_in_protein as u32,
        length: c.peptide.residues.len() as u16,
        flags: f,
    }
}

/// Records encoded per write: bounds the encode buffer at 20 MB.
const ENCODE_BLOCK_RECORDS: usize = 1 << 20;

pub fn build_base_peptide_index<W: Write>(
    idx: &SearchIndex,
    params: &SearchParams,
    decoy_prefix: &str,
    out: &mut W,
) -> io::Result<usize> {
    let mut records = base_peptide_records(idx, params, decoy_prefix);
    records.par_sort_by_key(|r| r.mass_milli);

    out.write_all(INDEX_MAGIC)?;
    out.write_all(&INDEX_FORMAT_VERSION.to_le_bytes())?;

    let n = records.len();
    let mut buf = vec![0u8; ENCODE_BLOCK_RECORDS.min(n) * INDEX_RECORD_SIZE];
    for block in records.chunks(ENCODE_BLOCK_RECORDS) {
        let bytes = &mut buf[..block.len() * INDEX_RECORD_SIZE];
        bytes
            .par_chunks_mut(INDEX_RECORD_SIZE)
            .zip(block.par_iter())
            .for_each(|(dst, rec)| dst.copy_from_slice(&rec.to_bytes()));
        out.write_all(bytes)?;
    }
    Ok(n)
}

struct KeyHasher(Sha256);

impl KeyHasher {
    fn u64(&mut self, v: u64) {
        self.0.update(v.to_le_bytes());
    }

    fn bytes(&mut self, b: &[u8]) {
        self.u64(b.len() as u64);
        self.0.update(b);
    }

    fn str(&mut self, s: &str) {
        self.bytes(s.as_bytes());
    }

    fn finish_hex(self) -> String {
        use std::fmt::Write;
        let digest = self.0.finalize();
        let mut out = String::with_capacity(32);
        for b in &digest[..16] {
            let _ = write!(out, "{b:02x}");
        }
        out
    }
}

pub fn index_cache_path(idx: &SearchIndex, params: &SearchParams, cache_dir: &Path) -> PathBuf {
    let mut h = KeyHasher(Sha256::new());

    h.u64(u64::from(INDEX_FORMAT_VERSION));
    h.u64(u64::from(CANDIDATE_SET_REVISION));
    h.str(env!("CARGO_PKG_VERSION"));

    for protein in &idx.db.proteins {
        h.str(&protein.accession);
        h.bytes(protein.sequence.as_ref());
        let variants = protein.variant_parts();
        h.u64(variants.len() as u64);
        for v in variants {
            h.u64(u64::from(v.start));
            h.u64(u64::from(v.end));
        }
        match protein.chain_ends() {
            Some(c) => {
                h.u64(1 + u64::from(c.n_term) + 2 * u64::from(c.c_term));
            }
            None => h.u64(0),
        }
    }
    h.str(&idx.decoy_prefix);
    h.str(
        idx.decoy_suffix
            .as_deref()
            .filter(|s| !s.is_empty())
            .unwrap_or(""),
    );

    let fixed = params.aa_set.fixed_mod_fingerprint();
    h.u64(fixed.len() as u64);
    for (res, loc_rank, delta_bits) in fixed {
        h.u64(u64::from(res));
        h.u64(u64::from(loc_rank));
        h.u64(delta_bits);
    }

    h.str(&format!("{:?}", params.enzyme));
    h.u64(params.extra_enzymes.len() as u64);
    for extra in &params.extra_enzymes {
        h.str(&format!("{extra:?}"));
    }
    h.u64(u64::from(params.max_missed_cleavages));
    h.u64(u64::from(params.min_length));
    h.u64(u64::from(params.max_length));
    h.u64(u64::from(params.num_tolerable_termini));
    h.u64(u64::from(params.require_nxst_sequon));

    cache_dir.join(format!("andes-candidx-{}.bin", h.finish_hex()))
}

fn cache_write_error(dir: &Path, e: io::Error) -> io::Error {
    io::Error::new(
        e.kind(),
        format!(
            "cannot write the candidate-index cache in {} ({e}); pass --index-cache-dir with a \
             writable directory",
            dir.display()
        ),
    )
}

fn lock_path_for(path: &Path) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(".lock");
    PathBuf::from(name)
}

const LOCK_HEARTBEAT: std::time::Duration = std::time::Duration::from_secs(5);
const LOCK_STALE: std::time::Duration =
    std::time::Duration::from_secs(if cfg!(test) { 2 } else { 60 });
const LOCK_POLL: std::time::Duration = std::time::Duration::from_millis(250);

struct BuildLock {
    path: PathBuf,
    stop: Option<std::sync::mpsc::Sender<()>>,
    beat: Option<std::thread::JoinHandle<()>>,
}

impl BuildLock {
    fn try_acquire(path: &Path) -> io::Result<Option<Self>> {
        let file = match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
        {
            Ok(f) => f,
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => return Ok(None),
            Err(e) => return Err(e),
        };
        let (stop, rx) = std::sync::mpsc::channel::<()>();
        let beat = std::thread::spawn(move || {
            while let Err(std::sync::mpsc::RecvTimeoutError::Timeout) =
                rx.recv_timeout(LOCK_HEARTBEAT)
            {
                let _ = file.set_modified(std::time::SystemTime::now());
            }
        });
        Ok(Some(Self {
            path: path.to_path_buf(),
            stop: Some(stop),
            beat: Some(beat),
        }))
    }
}

impl Drop for BuildLock {
    fn drop(&mut self) {
        drop(self.stop.take());
        if let Some(beat) = self.beat.take() {
            let _ = beat.join();
        }
        let _ = std::fs::remove_file(&self.path);
    }
}

fn wait_for_builder(path: &Path, lock_path: &Path) -> Option<MmapCandidateIndex> {
    log::info!(
        target: "candidate-index",
        "another process is building {}; waiting for it",
        path.display()
    );
    let mut seen_mtime = None;
    let mut seen_at = std::time::Instant::now();
    loop {
        if let Ok(mmap_index) = MmapCandidateIndex::open(path) {
            return Some(mmap_index);
        }
        let mtime = match std::fs::metadata(lock_path) {
            Ok(m) => m.modified().ok(),
            Err(_) => return None,
        };
        if mtime != seen_mtime {
            seen_mtime = mtime;
            seen_at = std::time::Instant::now();
        } else if seen_at.elapsed() > LOCK_STALE {
            log::warn!(
                target: "candidate-index",
                "the lock {} has not been refreshed for {} s; taking over",
                lock_path.display(),
                LOCK_STALE.as_secs()
            );
            let _ = std::fs::remove_file(lock_path);
            return None;
        }
        std::thread::sleep(LOCK_POLL);
    }
}

pub struct MmapCandidateIndex {
    mmap: memmap2::Mmap,
    len: usize,
    in_memory: bool,
}

impl MmapCandidateIndex {
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let byte_len = file.metadata()?.len() as usize;
        if byte_len < INDEX_HEADER_SIZE {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index file size {byte_len} is smaller than the {INDEX_HEADER_SIZE}-byte header"),
            ));
        }
        // SAFETY: The file is opened read-only and its contents are valid for
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        Self::from_mmap(mmap, false)
    }

    fn from_mmap(mmap: memmap2::Mmap, in_memory: bool) -> io::Result<Self> {
        let byte_len = mmap.len();
        if byte_len < INDEX_HEADER_SIZE {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "index size {byte_len} is smaller than the {INDEX_HEADER_SIZE}-byte header"
                ),
            ));
        }
        if &mmap[..INDEX_MAGIC.len()] != INDEX_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "index file is not an andes candidate index (bad magic)",
            ));
        }
        let version = u32::from_le_bytes(
            mmap[INDEX_MAGIC.len()..INDEX_HEADER_SIZE]
                .try_into()
                .expect("4 version bytes"),
        );
        if version != INDEX_FORMAT_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index file format version {version} != expected {INDEX_FORMAT_VERSION}"),
            ));
        }
        let data_len = byte_len - INDEX_HEADER_SIZE;
        if data_len % INDEX_RECORD_SIZE != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index payload {data_len} is not a multiple of {INDEX_RECORD_SIZE} bytes"),
            ));
        }
        let len = data_len / INDEX_RECORD_SIZE;
        Ok(Self {
            mmap,
            len,
            in_memory,
        })
    }

    pub fn build_in_memory(
        idx: &SearchIndex,
        params: &SearchParams,
        decoy_prefix: &str,
    ) -> io::Result<Self> {
        let mut image = Vec::new();
        build_base_peptide_index(idx, params, decoy_prefix, &mut image)?;
        let mut anon = memmap2::MmapMut::map_anon(image.len())?;
        anon.copy_from_slice(&image);
        drop(image);
        Self::from_mmap(anon.make_read_only()?, true)
    }

    pub fn is_in_memory(&self) -> bool {
        self.in_memory
    }

    pub fn open_or_build(
        path: &Path,
        idx: &SearchIndex,
        params: &SearchParams,
        decoy_prefix: &str,
    ) -> io::Result<(Self, bool)> {
        if let Ok(mmap_index) = Self::open(path) {
            return Ok((mmap_index, false));
        }

        let parent = path.parent().unwrap_or_else(|| Path::new("."));
        std::fs::create_dir_all(parent).map_err(|e| cache_write_error(parent, e))?;
        let lock_path = lock_path_for(path);
        let _lock = loop {
            match BuildLock::try_acquire(&lock_path) {
                Ok(Some(lock)) => break lock,
                Ok(None) => {
                    if let Some(mmap_index) = wait_for_builder(path, &lock_path) {
                        return Ok((mmap_index, false));
                    }
                }
                Err(e) => return Err(cache_write_error(parent, e)),
            }
        };
        if let Ok(mmap_index) = Self::open(path) {
            return Ok((mmap_index, false));
        }
        static BUILD_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let seq = BUILD_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let tmp_path = parent.join(format!(
            ".andes-candidx-build-{}-{}.tmp",
            std::process::id(),
            seq,
        ));
        let file = std::fs::File::create(&tmp_path).map_err(|e| cache_write_error(parent, e))?;
        let written = {
            use std::io::{BufWriter, Write};
            let mut bw = BufWriter::new(file);
            build_base_peptide_index(idx, params, decoy_prefix, &mut bw).and_then(|_| bw.flush())
        };
        if let Err(e) = written.and_then(|()| std::fs::rename(&tmp_path, path)) {
            let _ = std::fs::remove_file(&tmp_path);
            return Err(cache_write_error(parent, e));
        }

        let mmap_index = Self::open(path)?;
        Ok((mmap_index, true))
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn record_at(&self, i: usize) -> IndexRecord {
        assert!(
            i < self.len,
            "record_at index {i} out of bounds (len={})",
            self.len
        );
        let start = INDEX_HEADER_SIZE + i * INDEX_RECORD_SIZE;
        let chunk: &[u8; INDEX_RECORD_SIZE] = self.mmap[start..start + INDEX_RECORD_SIZE]
            .try_into()
            .expect("slice length is always INDEX_RECORD_SIZE");
        IndexRecord::from_bytes(chunk)
    }

    pub fn records(&self) -> Vec<IndexRecord> {
        (0..self.len).map(|i| self.record_at(i)).collect()
    }

    pub fn mass_window(&self, lo_milli: u64, hi_milli: u64) -> Vec<IndexRecord> {
        if self.len == 0 || lo_milli > hi_milli {
            return Vec::new();
        }
        let lo_idx = self.partition_point(|mass| mass < lo_milli);
        let hi_idx = self.partition_point(|mass| mass <= hi_milli);
        (lo_idx..hi_idx).map(|i| self.record_at(i)).collect()
    }

    fn partition_point<F>(&self, predicate: F) -> usize
    where
        F: Fn(u64) -> bool,
    {
        let mut lo = 0usize;
        let mut hi = self.len;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let mass = self.mass_milli_at(mid);
            if predicate(mass) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    fn mass_milli_at(&self, i: usize) -> u64 {
        let start = INDEX_HEADER_SIZE + i * INDEX_RECORD_SIZE;
        u64::from_le_bytes(
            self.mmap[start..start + 8]
                .try_into()
                .expect("always 8 bytes"),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidates::enumerate_candidates;
    use crate::index::SearchIndex;
    use crate::params::SearchParams;
    use model::amino_acid::set::AminoAcidSetBuilder;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use model::protein::{Protein, ProteinDb};

    fn build_fixture_index() -> tempfile::NamedTempFile {
        use std::io::BufWriter;
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        let tmp = tempfile::NamedTempFile::new().expect("tempfile");
        {
            let mut bw = BufWriter::new(tmp.as_file());
            build_base_peptide_index(&idx, &params, "XXX", &mut bw)
                .expect("build_base_peptide_index failed");
        }
        tmp
    }

    #[test]
    fn mmap_window_returns_in_range_records() {
        let tmp = build_fixture_index();
        let mi = MmapCandidateIndex::open(tmp.path()).unwrap();
        let all = mi.records();
        let target = all[all.len() / 2].mass_milli;
        let win = mi.mass_window(target - 5, target + 5);
        assert!(
            win.iter()
                .all(|r| r.mass_milli >= target - 5 && r.mass_milli <= target + 5),
            "mass_window returned out-of-range records"
        );
        assert!(
            win.iter().any(|r| r.mass_milli == target),
            "target mass not found in window"
        );
        assert!(
            mi.mass_window(0, all[0].mass_milli - 1).is_empty(),
            "expected empty window below minimum mass"
        );
    }

    fn read_records(data: &[u8]) -> Vec<IndexRecord> {
        assert!(
            data.len() >= INDEX_HEADER_SIZE && &data[..INDEX_MAGIC.len()] == INDEX_MAGIC,
            "binary blob missing the index header"
        );
        let payload = &data[INDEX_HEADER_SIZE..];
        assert_eq!(
            payload.len() % INDEX_RECORD_SIZE,
            0,
            "payload is not a multiple of {INDEX_RECORD_SIZE} bytes"
        );
        payload
            .chunks_exact(INDEX_RECORD_SIZE)
            .map(|chunk| {
                let arr: &[u8; INDEX_RECORD_SIZE] = chunk.try_into().unwrap();
                IndexRecord::from_bytes(arr)
            })
            .collect()
    }

    fn make_toy_index() -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "toy protein".into(),
                sequence: b"MKWVRPK".to_vec(),
                annotation: None,
            }],
        };
        SearchIndex::from_target_db(&target, "XXX")
    }

    #[test]
    fn builder_writes_mass_sorted_base_peptides() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        let mut buf: Vec<u8> = Vec::new();
        let n = build_base_peptide_index(&idx, &params, "XXX", &mut buf).expect("write failed");

        assert!(n > 0, "expected at least one record, got 0");
        assert_eq!(
            buf.len(),
            INDEX_HEADER_SIZE + n * INDEX_RECORD_SIZE,
            "buffer size mismatch"
        );

        let records = read_records(&buf);
        assert_eq!(records.len(), n);

        for window in records.windows(2) {
            assert!(
                window[0].mass_milli <= window[1].mass_milli,
                "records not sorted: {} > {}",
                window[0].mass_milli,
                window[1].mass_milli
            );
        }

        assert!(
            records.iter().any(|r| !r.is_decoy()),
            "no target records found"
        );

        assert!(
            records.iter().any(|r| r.is_decoy()),
            "no decoy records found"
        );

        assert!(
            records.iter().any(|r| r.length == 3),
            "expected at least one length-3 record (WVR from MKWVRPK)"
        );
    }

    #[test]
    fn index_record_round_trips_through_bytes() {
        let rec = IndexRecord {
            mass_milli: 1_234_567_890,
            protein_index: 42,
            start_offset: 7,
            length: 9,
            flags: flags::IS_DECOY | flags::IS_PROTEIN_C_TERM,
        };
        let bytes = rec.to_bytes();
        assert_eq!(bytes.len(), INDEX_RECORD_SIZE);
        let back = IndexRecord::from_bytes(&bytes);
        assert_eq!(back, rec);
        assert!(back.is_decoy());
    }

    #[test]
    fn base_params_forces_zero_variable_mods() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 3;

        let records = base_peptide_records(&idx, &params, "XXX");
        assert!(!records.is_empty(), "base records should not be empty");
    }

    #[test]
    fn base_peptide_records_matches_unmodified_enumerate_candidates() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P2".into(),
                description: "test protein".into(),
                sequence: b"WMKDALERQPK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");

        let ox_m = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(ox_m)
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 1;

        let base_records = base_peptide_records(&idx, &params, "XXX");
        let base_params = SearchParams {
            max_variable_mods_per_peptide: 0,
            max_label_sites: 0,
            ..params.clone()
        };
        let base_candidates: Vec<_> = enumerate_candidates(&idx, &base_params, "XXX").collect();
        assert_eq!(
            base_records.len(),
            base_candidates.len(),
            "base_peptide_records and enumerate_candidates(NumMods=0) must yield same count"
        );

        let mut lhs: Vec<(Vec<u8>, u64)> = base_candidates
            .iter()
            .map(|c| {
                let bytes: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
                let mass_milli = (c.peptide.mass() * 1000.0).round() as u64;
                (bytes, mass_milli)
            })
            .collect();
        lhs.sort_unstable();

        let all_candidates: Vec<_> = enumerate_candidates(&idx, &params, "XXX").collect();
        let mut rhs: Vec<(Vec<u8>, u64)> = all_candidates
            .iter()
            .filter(|c| {
                c.peptide
                    .residues
                    .iter()
                    .all(|aa| aa.mod_.as_ref().map(|m| m.fixed).unwrap_or(true))
            })
            .map(|c| {
                let bytes: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
                let mass_milli = (c.peptide.mass() * 1000.0).round() as u64;
                (bytes, mass_milli)
            })
            .collect();
        rhs.sort_unstable();

        assert_eq!(
            lhs, rhs,
            "base_peptide_records must produce exactly the unmodified subset of enumerate_candidates"
        );

        assert!(!lhs.is_empty(), "expected at least one base peptide");

        let full_count = all_candidates.len();
        let base_count = lhs.len();
        assert!(
            full_count > base_count,
            "full enumeration ({full_count}) should exceed base ({base_count}) when variable mods present"
        );
    }

    #[test]
    fn unwritable_cache_dir_is_an_error_naming_the_flag() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        let dir = tempfile::tempdir().expect("tempdir");
        let blocker = dir.path().join("blocker");
        std::fs::write(&blocker, b"x").unwrap();
        let path = blocker.join("sub").join("cache.idx");
        let err = match MmapCandidateIndex::open_or_build(&path, &idx, &params, "XXX") {
            Ok(_) => panic!("an unwritable cache directory must fail"),
            Err(e) => e,
        };
        assert!(
            err.to_string().contains("--index-cache-dir"),
            "error must name the flag: {err}"
        );
        assert!(!path.exists());
    }

    #[test]
    fn missing_cache_dir_is_created() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        let dir = tempfile::tempdir().expect("tempdir");
        let cache_dir = dir.path().join("a").join("b");
        let path = index_cache_path(&idx, &params, &cache_dir);
        let (mi, built) = MmapCandidateIndex::open_or_build(&path, &idx, &params, "XXX").unwrap();
        assert!(built && !mi.is_in_memory());
        assert!(path.is_file());
    }

    #[test]
    fn index_cache_path_is_in_the_given_dir() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = SearchParams::default_tryptic(aa_set);
        let a = index_cache_path(&idx, &params, Path::new("/one"));
        let b = index_cache_path(&idx, &params, Path::new("/two"));
        assert_eq!(a.parent(), Some(Path::new("/one")));
        assert_eq!(b.parent(), Some(Path::new("/two")));
        assert_eq!(a.file_name(), b.file_name());
    }

    #[test]
    fn open_or_build_reuses_existing_cache_file() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        let tmp = tempfile::NamedTempFile::new().expect("NamedTempFile");
        let cache_path = tmp.path();

        let (mi1, was_built1) = MmapCandidateIndex::open_or_build(cache_path, &idx, &params, "XXX")
            .expect("first open_or_build failed");
        assert!(was_built1, "first call must build the index (cache miss)");
        assert!(!mi1.is_empty(), "built index must be non-empty");

        let mtime_after_build = std::fs::metadata(cache_path)
            .expect("metadata after build")
            .modified()
            .expect("mtime after build");

        let (mi2, was_built2) = MmapCandidateIndex::open_or_build(cache_path, &idx, &params, "XXX")
            .expect("second open_or_build failed");
        assert!(!was_built2, "second call must reuse the cache (cache hit)");

        let mtime_after_reuse = std::fs::metadata(cache_path)
            .expect("metadata after reuse")
            .modified()
            .expect("mtime after reuse");
        assert_eq!(
            mtime_after_build, mtime_after_reuse,
            "cache file mtime changed on reuse — index was rebuilt unnecessarily"
        );

        let records1 = mi1.records();
        let records2 = mi2.records();
        assert_eq!(
            records1, records2,
            "cached records must be identical to freshly-built records"
        );
    }

    #[test]
    fn open_rejects_header_less_file_and_rebuilds() {
        use std::io::Write;
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        let tmp = tempfile::NamedTempFile::new().expect("NamedTempFile");
        tmp.as_file()
            .write_all(&[0u8; 2 * INDEX_RECORD_SIZE])
            .unwrap();
        tmp.as_file().sync_all().unwrap();

        assert!(
            MmapCandidateIndex::open(tmp.path()).is_err(),
            "open must reject a file with no magic header"
        );
        let (mi, was_built) = MmapCandidateIndex::open_or_build(tmp.path(), &idx, &params, "XXX")
            .expect("open_or_build");
        assert!(was_built, "header-less file must force a rebuild");
        assert!(!mi.is_empty(), "rebuilt index must be non-empty");
    }

    #[test]
    fn index_cache_path_distinguishes_different_params() {
        let idx = make_toy_index();
        let aa_set_a = AminoAcidSetBuilder::new_standard().build().unwrap();
        let aa_set_b = AminoAcidSetBuilder::new_standard().build().unwrap();

        let mut params_a = SearchParams::default_tryptic(aa_set_a);
        params_a.min_length = 3;

        let mut params_b = SearchParams::default_tryptic(aa_set_b);
        params_b.min_length = 6;

        let path_a = index_cache_path(&idx, &params_a, Path::new("cache"));
        let path_b = index_cache_path(&idx, &params_b, Path::new("cache"));

        assert_ne!(
            path_a,
            path_b,
            "cache paths must differ when min_length differs ({} vs {})",
            path_a.display(),
            path_b.display()
        );
    }

    #[test]
    fn index_cache_path_distinguishes_fixed_mods() {
        let idx = make_toy_index();
        let mut p_nomod =
            SearchParams::default_tryptic(AminoAcidSetBuilder::new_standard().build().unwrap());
        p_nomod.min_length = 3;
        let mut p_cam = SearchParams::default_tryptic(
            AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
                .build()
                .unwrap(),
        );
        p_cam.min_length = 3;
        assert_ne!(
            index_cache_path(&idx, &p_nomod, Path::new("cache")),
            index_cache_path(&idx, &p_cam, Path::new("cache")),
            "cache key must differ when the fixed-mod set differs"
        );
    }

    #[test]
    fn index_cache_path_keys_on_the_decoy_suffix() {
        let db = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: String::new(),
                    sequence: b"MNGSK".to_vec(),
                    annotation: None,
                },
                Protein {
                    accession: "P1_rev".into(),
                    description: String::new(),
                    sequence: b"ASNVTR".to_vec(),
                    annotation: None,
                },
            ],
        };
        let strategy = crate::decoy::DecoyStrategy::None;
        let seed = crate::decoy::DEFAULT_DECOY_SEED;
        let plain = SearchIndex::from_target_db_with_strategy(&db, "XXX", None, strategy, seed);
        let suffixed =
            SearchIndex::from_target_db_with_strategy(&db, "XXX", Some("_rev"), strategy, seed);
        let mut params =
            SearchParams::default_tryptic(AminoAcidSetBuilder::new_standard().build().unwrap());
        params.min_length = 3;
        params.require_nxst_sequon = true;

        let dir = tempfile::tempdir().unwrap();
        let plain_path = index_cache_path(&plain, &params, dir.path());
        let suffixed_path = index_cache_path(&suffixed, &params, dir.path());
        let (first, _) = MmapCandidateIndex::open_or_build(&plain_path, &plain, &params, "XXX")
            .expect("build the suffix-less index");
        assert!(first.records().iter().all(|r| !r.is_decoy()));
        let (second, _) =
            MmapCandidateIndex::open_or_build(&suffixed_path, &suffixed, &params, "XXX")
                .expect("open the suffixed index");
        let decoys = second.records().iter().filter(|r| r.is_decoy()).count();
        assert_eq!(decoys, 1, "ASNVTR from P1_rev is the one sequon decoy");
        assert_ne!(
            plain_path, suffixed_path,
            "the decoy suffix must key the cache"
        );
    }

    #[test]
    fn index_cache_path_is_stable_for_identical_params() {
        let idx = make_toy_index();
        let aa_set_a = AminoAcidSetBuilder::new_standard().build().unwrap();
        let aa_set_b = AminoAcidSetBuilder::new_standard().build().unwrap();

        let mut params_a = SearchParams::default_tryptic(aa_set_a);
        params_a.min_length = 4;
        let mut params_b = SearchParams::default_tryptic(aa_set_b);
        params_b.min_length = 4;

        let path_a = index_cache_path(&idx, &params_a, Path::new("cache"));
        let path_b = index_cache_path(&idx, &params_b, Path::new("cache"));

        assert_eq!(
            path_a, path_b,
            "cache paths must be identical for identical params"
        );
    }

    #[test]
    fn index_cache_key_is_pinned() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 4;
        let path = index_cache_path(&idx, &params, Path::new("cache"));
        assert_eq!(
            path.file_name().and_then(|n| n.to_str()),
            Some("andes-candidx-db3e57f5c212f6468961501684c9872b.bin")
        );
    }

    fn toy_params() -> SearchParams {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params
    }

    #[test]
    fn concurrent_open_or_build_builds_once() {
        let idx = make_toy_index();
        let params = toy_params();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("cache.bin");
        let built: Vec<bool> = std::thread::scope(|s| {
            let handles: Vec<_> = (0..6)
                .map(|_| {
                    s.spawn(|| {
                        MmapCandidateIndex::open_or_build(&path, &idx, &params, "XXX")
                            .unwrap()
                            .1
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        assert_eq!(built.iter().filter(|b| **b).count(), 1, "{built:?}");
        assert!(!lock_path_for(&path).exists(), "lock must be released");
    }

    #[test]
    fn waiter_opens_the_file_the_lock_holder_builds() {
        let idx = make_toy_index();
        let params = toy_params();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("cache.bin");
        let lock = BuildLock::try_acquire(&lock_path_for(&path))
            .unwrap()
            .expect("lock is free");
        std::thread::scope(|s| {
            let waiter =
                s.spawn(|| MmapCandidateIndex::open_or_build(&path, &idx, &params, "XXX").unwrap());
            std::thread::sleep(std::time::Duration::from_millis(600));
            assert!(!waiter.is_finished(), "waiter must not build while locked");
            let mut f = std::fs::File::create(&path).unwrap();
            build_base_peptide_index(&idx, &params, "XXX", &mut f).unwrap();
            drop(f);
            drop(lock);
            let (mi, was_built) = waiter.join().unwrap();
            assert!(!was_built, "waiter must reuse the holder's file");
            assert!(!mi.is_empty());
        });
    }

    #[test]
    fn stale_lock_is_taken_over() {
        let idx = make_toy_index();
        let params = toy_params();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("cache.bin");
        std::fs::File::create(lock_path_for(&path)).unwrap();
        let (mi, was_built) =
            MmapCandidateIndex::open_or_build(&path, &idx, &params, "XXX").unwrap();
        assert!(was_built);
        assert!(!mi.is_empty());
        assert!(!lock_path_for(&path).exists());
    }
}
