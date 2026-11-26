use std::sync::Arc;

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use model::enzyme::Enzyme;
use model::hash::{fnv1a64, splitmix64};
use model::protein::{ModSite, Protein, ProteinAnnotation, ProteinDb, SiteAnchor, VariantRegion};

use crate::candidates::{compute_cleavage_positions, has_met_cleaved_walk};
use crate::params::SearchParams;

pub const DEFAULT_DECOY_PREFIX: &str = "XXX";

pub const DEFAULT_DECOY_SEED: u64 = 42;

const MAX_REDRAWS: u64 = 20;

mod collisions;
pub mod entrapment;

pub use collisions::LengthOverlap;

/// The digestion the decoys are built for: the enzymes and specificity that place the cleavage
/// sites, and the peptide length and missed-cleavage limits of the search, so every decoy peptide
/// the search enumerates is checked against the target peptides it enumerates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Proteases {
    pub primary: Enzyme,
    pub extra: Vec<Enzyme>,
    pub num_tolerable_termini: u8,
    pub min_length: u32,
    pub max_length: u32,
    pub max_missed_cleavages: u32,
}

impl Default for Proteases {
    /// Trypsin, fully specific, with the search command's default lengths (6 to 50) and
    /// missed cleavages (2).
    fn default() -> Self {
        Self {
            primary: Enzyme::Trypsin,
            extra: Vec::new(),
            num_tolerable_termini: 2,
            min_length: 6,
            max_length: 50,
            max_missed_cleavages: 2,
        }
    }
}

impl Proteases {
    pub fn of(params: &SearchParams) -> Self {
        Self {
            primary: params.enzyme,
            extra: params.extra_enzymes.clone(),
            num_tolerable_termini: params.num_tolerable_termini,
            min_length: params.min_length,
            max_length: params.max_length,
            max_missed_cleavages: params.max_missed_cleavages,
        }
    }

    fn all(&self) -> impl Iterator<Item = Enzyme> + '_ {
        std::iter::once(self.primary).chain(self.extra.iter().copied())
    }

    pub fn is_non_specific(&self) -> bool {
        self.num_tolerable_termini == 0 || self.all().any(|e| e == Enzyme::NonSpecific)
    }

    fn makes_cut(&self, r: u8) -> bool {
        self.all()
            .any(|e| e.is_cleavable_after(r) || e.is_cleavable_before(r))
    }

    /// The enzymes, the primary one first: `Trypsin` or `Trypsin,LysC`.
    pub fn enzymes_name(&self) -> String {
        self.all().map(Enzyme::name).collect::<Vec<_>>().join(",")
    }

    /// The specificity as the search's `--enzyme-specificity` names it.
    pub fn specificity_name(&self) -> &'static str {
        match self.num_tolerable_termini {
            2 => "fully",
            1 => "semi",
            _ => "non-specific",
        }
    }
}

/// SHA-256 (hex) of the target proteins a decoy database is built from: the accession and the
/// residues of every entry that carries neither decoy affix, in order. Descriptions are left out;
/// the decoys do not depend on them.
pub fn targets_sha256(db: &ProteinDb, prefix: &str, suffix: Option<&str>) -> String {
    use sha2::{Digest, Sha256};
    let mut h = Sha256::new();
    for p in db
        .proteins
        .iter()
        .filter(|p| !is_decoy_accession_affix(&p.accession, prefix, suffix))
    {
        h.update((p.accession.len() as u64).to_le_bytes());
        h.update(p.accession.as_bytes());
        h.update((p.sequence.len() as u64).to_le_bytes());
        h.update(&p.sequence);
    }
    hex(&h.finalize())
}

/// SHA-256 (hex) of a file's bytes.
pub fn file_sha256(path: &std::path::Path) -> std::io::Result<String> {
    use sha2::{Digest, Sha256};
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(hex(&h.finalize()))
}

fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write;
    bytes.iter().fold(String::with_capacity(2 * bytes.len()), |mut s, b| {
        let _ = write!(s, "{b:02x}");
        s
    })
}

/// The decoy key: what the decoys of a database depend on and nothing else, the target proteins
/// (`targets_sha256`), the strategy, the seed, the decoy prefix and suffix, and the digests the
/// peptide-level check covered (the build digest first). Mods, tolerances, charges and every
/// other search setting are left out. 32 hex characters.
pub fn decoy_key(
    targets_sha256: &str,
    strategy: DecoyStrategy,
    seed: u64,
    prefix: &str,
    suffix: Option<&str>,
    digests: &[Proteases],
) -> String {
    use sha2::{Digest, Sha256};
    let mut h = Sha256::new();
    let mut put = |s: &str| {
        h.update((s.len() as u64).to_le_bytes());
        h.update(s.as_bytes());
    };
    put("andes.decoy_key.1");
    put(targets_sha256);
    put(strategy.name());
    put(&seed.to_string());
    put(&normalize_decoy_prefix(prefix));
    put(suffix.unwrap_or(""));
    for d in digests {
        put(&d.enzymes_name());
        put(d.specificity_name());
        put(&format!("{}-{}/{}", d.min_length, d.max_length, d.max_missed_cleavages));
    }
    hex(&h.finalize()[..16])
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecoyStrategy {
    #[default]
    Shuffle,
    Reverse,
    Sequon,
    Mutate,
    None,
}

impl DecoyStrategy {
    /// The `--decoy-strategy` name.
    pub fn name(self) -> &'static str {
        match self {
            Self::Shuffle => "shuffle",
            Self::Reverse => "reverse",
            Self::Sequon => "sequon",
            Self::Mutate => "mutate",
            Self::None => "none",
        }
    }

    pub fn from_name(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "reverse" | "rev" => Some(Self::Reverse),
            "shuffle" | "shuf" => Some(Self::Shuffle),
            "sequon" => Some(Self::Sequon),
            "mutate" => Some(Self::Mutate),
            "none" | "off" => Some(Self::None),
            _ => None,
        }
    }
}

pub fn decoy_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
) -> ProteinDb {
    if strategy == DecoyStrategy::None {
        return ProteinDb {
            proteins: Vec::new(),
        };
    }
    checked_decoy_db(target, prefix, strategy, seed, proteases, &[]).decoys
}

/// Decoys built for one digest and checked against the target peptides of further digests.
#[derive(Debug, Clone)]
pub struct CheckedDecoys {
    pub decoys: ProteinDb,
    /// Decoy peptides that equalled a target peptide and were redrawn.
    pub redrawn: usize,
    /// Decoy peptides that still equal a target peptide after the redraws.
    pub left: usize,
    /// Distinct target peptides (I and L folded) of all digests together.
    pub target_peptides: usize,
}

/// The decoys of `strategy` built for the `build` digest, as [`decoy_db`] builds them, with the
/// peptide-level check run against the target peptides of `build` and of every digest in
/// `also`: a cohort whose runs search the same decoys with different enzymes, specificities,
/// lengths or missed cleavages lists its other digests there. With `also` empty the decoys are
/// those of [`decoy_db`].
pub fn checked_decoy_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    build: &Proteases,
    also: &[Proteases],
) -> CheckedDecoys {
    if strategy == DecoyStrategy::None {
        return CheckedDecoys {
            decoys: ProteinDb {
                proteins: Vec::new(),
            },
            redrawn: 0,
            left: 0,
            target_peptides: 0,
        };
    }
    let t0 = std::time::Instant::now();
    let digests: Vec<Proteases> = std::iter::once(build.clone())
        .chain(also.iter().filter(|d| *d != build).cloned())
        .collect();
    let (proteins, repairs, target_peptides) =
        build_decoys(target, prefix, strategy, seed, &digests);
    collisions::log_repairs(&repairs, target_peptides, t0.elapsed().as_secs_f64());
    CheckedDecoys {
        decoys: ProteinDb { proteins },
        redrawn: repairs.iter().map(|r| r.found).sum(),
        left: repairs.iter().map(|r| r.left).sum(),
        target_peptides,
    }
}

/// Per peptide length, the distinct target and decoy peptides that `digest` enumerates and the
/// decoy peptides equal to a target peptide (I and L folded): the overlap report of a decoy
/// database.
pub fn peptide_overlap(
    targets: &ProteinDb,
    decoys: &ProteinDb,
    digest: &Proteases,
) -> Vec<LengthOverlap> {
    collisions::peptide_overlap(targets, decoys, digest)
}

/// The decoy proteins of `strategy` (not `None`) built for `digests[0]`, the outcome of the
/// peptide-level check against the target peptides of every digest for each, and the number of
/// distinct target peptides checked against.
fn build_decoys(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    digests: &[Proteases],
) -> (Vec<Protein>, Vec<collisions::Repair>, usize) {
    let proteases = &digests[0];
    let prefix = normalize_decoy_prefix(prefix);
    let peptides = collisions::TargetPeptides::union(target, digests);
    let built = match strategy {
        DecoyStrategy::Mutate => mutated_decoys(target, &prefix, seed, proteases, &peptides),
        _ if proteases.is_non_specific() => {
            whole_sequence_decoys(target, &prefix, strategy, seed, proteases, Some(&peptides))
        }
        _ => peptide_level_decoys(target, &prefix, strategy, seed, proteases, &peptides),
    };
    let (proteins, repairs) = built.into_iter().unzip();
    (proteins, repairs, peptides.len())
}

pub(crate) fn whole_sequence_decoy_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
) -> ProteinDb {
    let prefix = normalize_decoy_prefix(prefix);
    let proteins = match strategy {
        DecoyStrategy::None => Vec::new(),
        DecoyStrategy::Mutate => whole_sequence_decoys(
            target,
            &prefix,
            DecoyStrategy::Reverse,
            seed,
            &Proteases::default(),
            None,
        ),
        _ => whole_sequence_decoys(target, &prefix, strategy, seed, &Proteases::default(), None),
    };
    ProteinDb {
        proteins: proteins.into_iter().map(|(p, _)| p).collect(),
    }
}

const MUTATE_MIN_SEGMENT: usize = 7;
const MUTATE_POOL: &[u8] = b"GASVTLNDQEHFYW";

fn mutated_decoys(
    target: &ProteinDb,
    prefix: &str,
    seed: u64,
    proteases: &Proteases,
    peptides: &collisions::TargetPeptides,
) -> Vec<(Protein, collisions::Repair)> {
    let segments = TargetSegments::new(target, proteases);
    let mut protected: FxHashMap<(u64, bool), Vec<usize>> = FxHashMap::default();
    for p in &target.proteins {
        let sites = p.mod_sites();
        if sites.is_empty() {
            continue;
        }
        let met = has_met_cleaved_walk(p);
        let cuts = compute_cleavage_positions(&p.sequence, proteases.primary, &proteases.extra);
        for w in cuts.windows(2) {
            let (a, b) = (w[0] as usize, w[1] as usize);
            let local = sites
                .iter()
                .map(|s| s.position as usize)
                .filter(|&i| (a..b).contains(&i))
                .map(|i| i - a);
            let key = (fnv1a64(&p.sequence[a..b]), met && a == 0);
            protected.entry(key).or_default().extend(local);
        }
    }
    target
        .proteins
        .par_iter()
        .map(|p| {
            let mut sequence = p.sequence.clone();
            let met = has_met_cleaved_walk(p);
            let cuts = compute_cleavage_positions(&sequence, proteases.primary, &proteases.extra);
            for w in cuts.windows(2) {
                let (a, b) = (w[0] as usize, w[1] as usize);
                let seg = &p.sequence[a..b];
                let fixed = protected
                    .get(&(fnv1a64(seg), met && a == 0))
                    .map_or(&[][..], Vec::as_slice);
                if let Some((at, r)) =
                    segment_mutation(seg, met && a == 0, fixed, seed, proteases, &segments)
                {
                    sequence[a + at] = r;
                }
            }
            let decoy = Protein {
                accession: format!("{}_{}", prefix, p.accession),
                description: p.description.clone(),
                sequence,
                annotation: p.annotation.clone(),
            };
            collisions::repair_mutated(p, decoy, seed, peptides)
        })
        .collect()
}

fn segment_mutation(
    seg: &[u8],
    lead_met: bool,
    fixed: &[usize],
    seed: u64,
    proteases: &Proteases,
    targets: &TargetSegments,
) -> Option<(usize, u8)> {
    let m = seg.len();
    if m < MUTATE_MIN_SEGMENT {
        return None;
    }
    let mass = |r: u8| model::AminoAcid::standard(r).map(|a| a.mass);
    let open: Vec<usize> = (1..m - 1)
        .filter(|&i| {
            let fixed = (lead_met && i == 1)
                || fixed.contains(&i)
                || proteases.makes_cut(seg[i])
                || mass(seg[i]).is_none();
            !fixed
        })
        .collect();
    if open.is_empty() {
        return None;
    }
    let mut fallback = None;
    for k in 0..=MAX_REDRAWS {
        let h = segment_seed(seed, seg, lead_met, k);
        let at = open[(h % open.len() as u64) as usize];
        let old = mass(seg[at])?;
        let pool: Vec<u8> = MUTATE_POOL
            .iter()
            .copied()
            .filter(|&r| mass(r).is_some_and(|m| (m - old).abs() > 1.1) && !proteases.makes_cut(r))
            .collect();
        if pool.is_empty() {
            continue;
        }
        let r = pool[(splitmix64(h) % pool.len() as u64) as usize];
        let mut s = seg.to_vec();
        s[at] = r;
        let identity: Vec<u32> = (0..m as u32).collect();
        if !targets.collides(&s, &identity, lead_met) {
            return Some((at, r));
        }
        fallback.get_or_insert((at, r));
    }
    fallback
}

fn decoy_protein(p: &Protein, prefix: &str, source: Vec<u32>, proteases: &Proteases) -> Protein {
    let accession = format!("{}_{}", prefix, p.accession);
    if p.annotation.is_some() {
        return annotated_decoy(p, accession, source, proteases);
    }
    Protein {
        accession,
        description: p.description.clone(),
        sequence: source.iter().map(|&i| p.sequence[i as usize]).collect(),
        annotation: None,
    }
}

fn peptide_level_decoys(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
    peptides: &collisions::TargetPeptides,
) -> Vec<(Protein, collisions::Repair)> {
    let segments = TargetSegments::new(target, proteases);
    target
        .proteins
        .par_iter()
        .map(|p| {
            let source = peptide_level_permutation(p, strategy, seed, proteases, &segments);
            collisions::repair_peptide_level(p, prefix, source, strategy, seed, proteases, peptides)
        })
        .collect()
}

/// Whole-sequence decoys of a non-specific search. With `peptides`, every decoy peptide the
/// search enumerates that equals a target peptide is redrawn; without, the decoys are returned
/// as built (the peptide-anchored refinement index, whose entries are single peptides).
fn whole_sequence_decoys(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
    peptides: Option<&collisions::TargetPeptides>,
) -> Vec<(Protein, collisions::Repair)> {
    target
        .proteins
        .par_iter()
        .enumerate()
        .map(|(i, p)| {
            let n = p.sequence.len() as u32;
            let source: Vec<u32> = match strategy {
                DecoyStrategy::Shuffle | DecoyStrategy::Sequon => {
                    let mut source: Vec<u32> = (0..n).collect();
                    shuffle_in_place(
                        &mut source,
                        seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15),
                    );
                    if strategy == DecoyStrategy::Sequon {
                        restore_sequons(&p.sequence, &mut source, &vec![true; n as usize]);
                    }
                    source
                }
                _ => (0..n).rev().collect(),
            };
            match peptides {
                Some(peptides) => {
                    collisions::repair_whole_sequence(p, prefix, source, seed, proteases, peptides)
                }
                None => (
                    decoy_protein(p, prefix, source, proteases),
                    collisions::Repair::default(),
                ),
            }
        })
        .collect()
}

fn peptide_level_permutation(
    p: &Protein,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
    targets: &TargetSegments,
) -> Vec<u32> {
    let seq = &p.sequence;
    let met = has_met_cleaved_walk(p);
    let mut source: Vec<u32> = (0..seq.len() as u32).collect();
    let cuts = compute_cleavage_positions(seq, proteases.primary, &proteases.extra);
    for w in cuts.windows(2) {
        let (a, b) = (w[0] as usize, w[1] as usize);
        let local = segment_permutation(
            &seq[a..b],
            met && a == 0,
            strategy,
            seed,
            proteases,
            targets,
        );
        for (k, from) in local.into_iter().enumerate() {
            source[a + k] = a as u32 + from;
        }
    }
    source
}

fn segment_permutation(
    seg: &[u8],
    lead_met: bool,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
    targets: &TargetSegments,
) -> Vec<u32> {
    let Some(first) = draw_segment(seg, lead_met, strategy, seed, proteases, 0) else {
        return (0..seg.len() as u32).collect();
    };
    if !targets.collides(seg, &first, lead_met) {
        return first;
    }
    (1..=MAX_REDRAWS)
        .filter_map(|k| draw_segment(seg, lead_met, strategy, seed, proteases, k))
        .find(|local| !targets.collides(seg, local, lead_met))
        .unwrap_or(first)
}

/// Draw `k` of a segment's permutation (segment-local source indices): the cleavage residues,
/// the cleavage-side terminal residue and a leading initiator Met stay in place; draw 0 of
/// `Reverse` reverses the movable residues, every other draw shuffles them with a seed of the
/// segment's residues and `k`, so a segment gets the same draw wherever it occurs. `None` when
/// fewer than two residues can move.
fn draw_segment(
    seg: &[u8],
    lead_met: bool,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
    k: u64,
) -> Option<Vec<u32>> {
    let m = seg.len();
    if m == 0 {
        return None;
    }
    let anchor = if proteases.primary.is_n_term() {
        0
    } else {
        m - 1
    };
    let movable: Vec<bool> = (0..m)
        .map(|i| i != anchor && !(lead_met && i == 0) && !proteases.makes_cut(seg[i]))
        .collect();
    let free: Vec<u32> = (0..m as u32).filter(|&i| movable[i as usize]).collect();
    if free.len() < 2 {
        return None;
    }
    let mut order = free.clone();
    if k == 0 && strategy == DecoyStrategy::Reverse {
        order.reverse();
    } else {
        shuffle_in_place(&mut order, segment_seed(seed, seg, lead_met, k));
    }
    let mut local: Vec<u32> = (0..m as u32).collect();
    for (&at, &from) in free.iter().zip(&order) {
        local[at as usize] = from;
    }
    if strategy == DecoyStrategy::Sequon {
        restore_sequons(seg, &mut local, &movable);
    }
    Some(local)
}

struct TargetSegments(Vec<u64>);

impl TargetSegments {
    fn new(db: &ProteinDb, proteases: &Proteases) -> Self {
        let mut hashes: Vec<u64> = db
            .proteins
            .iter()
            .flat_map(|p| {
                let seq = &p.sequence;
                let met = has_met_cleaved_walk(p);
                let cuts = compute_cleavage_positions(seq, proteases.primary, &proteases.extra);
                let mut out = Vec::with_capacity(cuts.len() + 1);
                for w in cuts.windows(2) {
                    let s = &seq[w[0] as usize..w[1] as usize];
                    out.push(folded_hash(s.iter().copied()));
                    if met && w[0] == 0 {
                        out.push(folded_hash(s[1..].iter().copied()));
                    }
                }
                out
            })
            .collect();
        hashes.sort_unstable();
        hashes.dedup();
        Self(hashes)
    }

    fn contains(&self, residues: impl Iterator<Item = u8>) -> bool {
        self.0.binary_search(&folded_hash(residues)).is_ok()
    }

    fn collides(&self, seg: &[u8], local: &[u32], lead_met: bool) -> bool {
        let residues = || local.iter().map(|&i| seg[i as usize]);
        self.contains(residues()) || (lead_met && self.contains(residues().skip(1)))
    }
}

fn folded_hash(residues: impl Iterator<Item = u8>) -> u64 {
    fnv1a64(residues.map(|r| if r == b'I' { b'L' } else { r }))
}

fn segment_seed(seed: u64, seg: &[u8], lead_met: bool, k: u64) -> u64 {
    let h = fnv1a64(seg);
    splitmix64(h ^ splitmix64(seed ^ splitmix64(k.wrapping_mul(2) | u64::from(lead_met))))
}

fn shuffle_in_place<T>(v: &mut [T], mut state: u64) {
    if state == 0 {
        state = 0xD1B5_4A32_D192_ED03;
    }
    for j in (1..v.len()).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let k = (state % (j as u64 + 1)) as usize;
        v.swap(j, k);
    }
}

fn annotated_decoy(
    p: &Protein,
    accession: String,
    mut source: Vec<u32>,
    proteases: &Proteases,
) -> Protein {
    fit_permutation(p, &mut source, proteases);
    Protein {
        accession,
        description: p.description.clone(),
        sequence: source.iter().map(|&i| p.sequence[i as usize]).collect(),
        annotation: decoy_annotation(p, &source),
    }
}

fn fit_permutation(target: &Protein, source: &mut [u32], proteases: &Proteases) {
    let seq = &target.sequence;
    let n = seq.len();
    if n == 0 {
        return;
    }
    let chain = target.chain_ends().is_some();
    if chain && !proteases.all().any(|e| e == Enzyme::NonSpecific) {
        let mut image = vec![0u32; n];
        for (i, &from) in source.iter().enumerate() {
            image[from as usize] = i as u32;
        }
        for (i, from) in source.iter_mut().enumerate() {
            *from = i as u32;
        }
        let cuts = compute_cleavage_positions(seq, proteases.primary, &proteases.extra);
        for w in cuts.windows(2) {
            let free: Vec<u32> = (w[0]..w[1])
                .filter(|&i| !proteases.makes_cut(seq[i as usize]))
                .collect();
            let mut order = free.clone();
            order.sort_by_key(|&i| image[i as usize]);
            if order == free {
                order.reverse();
            }
            for (&at, &from) in free.iter().zip(&order) {
                source[at as usize] = from;
            }
        }
    }
    let pin = |source: &mut [u32], q: usize| {
        if let Some(i) = source.iter().position(|&from| from as usize == q) {
            source.swap(i, q);
        }
    };
    let after_met = !chain && seq[0] == b'M';
    let (mut first, mut second, mut last) = (false, false, false);
    for s in target.mod_sites() {
        match s.anchor {
            SiteAnchor::NTerminus if s.position == 0 => first = true,
            SiteAnchor::NTerminus if s.position == 1 && after_met => second = true,
            SiteAnchor::CTerminus if s.position as usize + 1 == n => last = true,
            _ => {}
        }
    }
    if first || second {
        pin(source, 0);
    }
    if second {
        pin(source, 1);
    }
    if last {
        pin(source, n - 1);
    }
}

pub fn decoy_annotation(target: &Protein, source: &[u32]) -> Option<Arc<ProteinAnnotation>> {
    let ann = target.annotation.as_deref()?;
    let n = source.len();
    debug_assert_eq!(n, target.sequence.len());
    let mut image = vec![0u32; n];
    for (i, &src) in source.iter().enumerate() {
        image[src as usize] = i as u32;
    }
    let mod_sites: Vec<ModSite> = ann
        .mod_sites
        .iter()
        .filter(|s| (s.position as usize) < n)
        .map(|s| ModSite {
            position: image[s.position as usize],
            ..*s
        })
        .collect();
    let variant = ann.variant.iter().map(|&v| map_region(v, &image)).collect();
    Some(Arc::new(
        ProteinAnnotation::new(mod_sites, None)
            .with_variant(variant)
            .with_chain(ann.chain)
            .with_origin(ann.origin.clone())
            .with_substitutions(
                ann.substitutions
                    .iter()
                    .filter(|s| (s.position as usize) < n)
                    .map(|s| model::Substitution {
                        position: image[s.position as usize],
                        ..*s
                    })
                    .collect(),
            ),
    ))
}

fn map_region(v: VariantRegion, image: &[u32]) -> VariantRegion {
    let n = image.len();
    let (a, b) = (v.start as usize, v.end as usize);
    if a < b && b <= n {
        let lo = image[a..b].iter().copied().min().unwrap_or(0);
        let hi = image[a..b].iter().copied().max().unwrap_or(0);
        let width = (b - a) as u32;
        if hi - lo + 1 == width {
            VariantRegion {
                start: lo,
                end: hi + 1,
            }
        } else {
            let start = image[a].min(n as u32 - width);
            VariantRegion {
                start,
                end: start + width,
            }
        }
    } else if a >= b && b >= 1 && a < n {
        let (lo, hi) = (b - 1, a);
        let width = (hi - lo) as u32;
        let min = image[lo..=hi].iter().copied().min().unwrap_or(0);
        let max = image[lo..=hi].iter().copied().max().unwrap_or(0);
        let first = if max - min == width {
            min
        } else {
            image[lo].min(n as u32 - 1 - width)
        };
        VariantRegion {
            start: first + width,
            end: first + 1,
        }
    } else {
        v
    }
}

fn count_sequon_starts(s: &[u8]) -> usize {
    (0..s.len().saturating_sub(2))
        .filter(|&i| s[i] == b'N' && s[i + 1] != b'P' && (s[i + 2] == b'S' || s[i + 2] == b'T'))
        .count()
}

fn restore_sequons(orig: &[u8], source: &mut [u32], movable: &[bool]) {
    let n = orig.len();
    if n < 3 {
        return;
    }
    let target = count_sequon_starts(orig);
    let mut r: Vec<u8> = source.iter().map(|&i| orig[i as usize]).collect();
    let mut cur = count_sequon_starts(&r);
    let mut i = 0usize;
    while cur < target && i + 2 < n {
        if movable[i + 2]
            && r[i] == b'N'
            && r[i + 1] != b'P'
            && r[i + 2] != b'S'
            && r[i + 2] != b'T'
        {
            if let Some(j) = (i + 3..n).find(|&j| {
                movable[j]
                    && (r[j] == b'S' || r[j] == b'T')
                    && !(r[j - 2] == b'N' && r[j - 1] != b'P')
            }) {
                r.swap(i + 2, j);
                let new_count = count_sequon_starts(&r);
                if new_count > cur && new_count <= target {
                    cur = new_count;
                    source.swap(i + 2, j);
                } else {
                    r.swap(i + 2, j);
                }
            }
        }
        i += 1;
    }
}

pub fn build_search_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
) -> ProteinDb {
    build_search_db_for(target, prefix, strategy, seed, &Proteases::default())
}

pub fn build_search_db_for(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
) -> ProteinDb {
    let decoy = decoy_db(target, prefix, strategy, seed, proteases);
    let mut proteins = target.proteins.clone();
    proteins.extend(decoy.proteins);
    ProteinDb { proteins }
}

pub fn normalize_decoy_prefix(prefix: &str) -> String {
    let trimmed = prefix.trim().trim_end_matches('_');
    if trimmed.is_empty() {
        DEFAULT_DECOY_PREFIX.to_string()
    } else {
        trimmed.to_string()
    }
}

pub fn decoy_accession_needle(prefix: &str) -> String {
    format!("{}_", normalize_decoy_prefix(prefix))
}

pub fn is_decoy_accession(accession: &str, prefix: &str) -> bool {
    accession.starts_with(&decoy_accession_needle(prefix))
}

pub fn is_decoy_accession_affix(accession: &str, prefix: &str, suffix: Option<&str>) -> bool {
    accession.starts_with(&decoy_accession_needle(prefix))
        || suffix.is_some_and(|s| !s.is_empty() && accession.ends_with(s))
}

const FOREIGN_DECOY_TAGS: [&str; 8] = [
    "rev", "reverse", "reversed", "decoy", "xxx", "shuffle", "shuffled", "random",
];

fn strip_foreign_decoy_tag(accession: &str) -> Option<&str> {
    let is_tag = |t: &str| FOREIGN_DECOY_TAGS.iter().any(|k| t.eq_ignore_ascii_case(k));
    if let Some(i) = accession.find('_') {
        let rest = accession[i + 1..].trim_start_matches('_');
        if is_tag(&accession[..i]) && !rest.is_empty() {
            return Some(rest);
        }
    }
    if let Some(i) = accession.rfind('_') {
        let rest = accession[..i].trim_end_matches('_');
        if is_tag(&accession[i + 1..]) && !rest.is_empty() {
            return Some(rest);
        }
    }
    None
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ForeignDecoys {
    pub n_records: usize,
    pub n_reversed: usize,
    pub n_named_paired: usize,
    pub n_named_unpaired: usize,
    pub example_reversed: Option<String>,
    pub example_named: Option<String>,
}

fn min_reversible_len(n_records: usize) -> usize {
    let pairs = (n_records as f64) * (n_records as f64);
    let mut len = 1;
    let mut sequences = 20.0_f64;
    while sequences <= pairs {
        sequences *= 20.0;
        len += 1;
    }
    len
}

impl ForeignDecoys {
    pub fn scan(db: &ProteinDb) -> Self {
        let min_len = min_reversible_len(db.proteins.len());
        let sequences: FxHashSet<&[u8]> = db.proteins.iter().map(|p| &p.sequence[..]).collect();
        let accessions: FxHashSet<&str> = db.proteins.iter().map(|p| &p.accession[..]).collect();
        let mut found = ForeignDecoys {
            n_records: db.proteins.len(),
            ..Default::default()
        };
        let mut last_reversed: Option<&Protein> = None;
        let mut reversed: Vec<u8> = Vec::new();
        for p in &db.proteins {
            reversed.clear();
            reversed.extend(p.sequence.iter().rev());
            if p.sequence.len() >= min_len
                && reversed != p.sequence
                && sequences.contains(&reversed[..])
            {
                found.n_reversed += 1;
                last_reversed = Some(p);
            }
            if let Some(rest) = strip_foreign_decoy_tag(&p.accession) {
                if accessions.contains(rest) {
                    found.n_named_paired += 1;
                    found.example_named = Some(p.accession.clone());
                } else {
                    found.n_named_unpaired += 1;
                    found
                        .example_named
                        .get_or_insert_with(|| p.accession.clone());
                }
            }
        }
        found.example_reversed = last_reversed.map(|p| p.accession.clone());
        found
    }

    pub fn verdict(&self) -> Result<Option<String>, String> {
        const FIX: &str = "andes would search them as targets and add its own decoys, which \
                           underestimates the FDR. Pass a target-only FASTA, or declare them: \
                           `--decoy-prefix <tag> --decoy-strategy none` (`--decoy-suffix` for a \
                           tag written after the accession).";
        if self.n_named_paired > 0 {
            return Err(format!(
                "the database already contains decoys: {} of {} records are named as the decoy \
                 of another record (e.g. {}). {FIX}",
                self.n_named_paired,
                self.n_records,
                self.example_named.as_deref().unwrap_or("?"),
            ));
        }
        if self.n_reversed > 0 && 2 * self.n_reversed >= self.n_records {
            return Err(format!(
                "the database already contains decoys: {} of {} records are the exact reverse \
                 of another record (e.g. {}). {FIX}",
                self.n_reversed,
                self.n_records,
                self.example_reversed.as_deref().unwrap_or("?"),
            ));
        }
        if self.n_reversed > 0 || self.n_named_unpaired > 0 {
            return Ok(Some(format!(
                "WARN: the database may contain decoys: {} of {} records are the exact reverse \
                 of another record{}, {} carry a decoy tag in the accession{}. If they are \
                 decoys, {FIX}",
                self.n_reversed,
                self.n_records,
                self.example_reversed
                    .as_deref()
                    .map_or(String::new(), |a| format!(" (e.g. {a})")),
                self.n_named_unpaired,
                self.example_named
                    .as_deref()
                    .map_or(String::new(), |a| format!(" (e.g. {a})")),
            )));
        }
        Ok(None)
    }
}

#[cfg(test)]
#[path = "decoy/tests/decoy.rs"]
mod tests;
