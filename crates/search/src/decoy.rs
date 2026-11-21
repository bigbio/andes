use std::sync::Arc;

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use model::enzyme::Enzyme;
use model::hash::{fnv1a64, splitmix64};
use model::protein::{ModSite, Protein, ProteinAnnotation, ProteinDb, SiteAnchor, VariantRegion};

use crate::candidate_gen::{compute_cleavage_positions, has_met_cleaved_walk};
use crate::search_params::SearchParams;

pub const DEFAULT_DECOY_PREFIX: &str = "XXX";

pub const DEFAULT_DECOY_SEED: u64 = 42;

const MAX_REDRAWS: u64 = 20;

mod collisions;

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
    let t0 = std::time::Instant::now();
    let (proteins, repairs, target_peptides) =
        build_decoys(target, prefix, strategy, seed, proteases);
    collisions::log_repairs(&repairs, target_peptides, t0.elapsed().as_secs_f64());
    ProteinDb { proteins }
}

/// The decoy proteins of `strategy` (not `None`), the outcome of the peptide-level check for
/// each, and the number of distinct target peptides checked against.
fn build_decoys(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
) -> (Vec<Protein>, Vec<collisions::Repair>, usize) {
    let prefix = normalize_decoy_prefix(prefix);
    let peptides = collisions::TargetPeptides::new(target, proteases);
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
            collisions::repair_mutated(p, decoy, seed, proteases, peptides)
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
mod tests {
    use super::*;

    fn make_db(proteins: &[(&str, &[u8])]) -> ProteinDb {
        ProteinDb {
            proteins: proteins
                .iter()
                .map(|(acc, seq)| Protein {
                    accession: acc.to_string(),
                    description: String::new(),
                    sequence: seq.to_vec(),
                    annotation: None,
                })
                .collect(),
        }
    }

    #[test]
    fn foreign_decoy_tag_is_the_whole_first_or_last_field() {
        assert_eq!(
            strip_foreign_decoy_tag("rev_sp|P1|A_HUMAN"),
            Some("sp|P1|A_HUMAN")
        );
        assert_eq!(strip_foreign_decoy_tag("REV__P1"), Some("P1"));
        assert_eq!(strip_foreign_decoy_tag("DECOY_P1"), Some("P1"));
        assert_eq!(strip_foreign_decoy_tag("P1_rev"), Some("P1"));
        assert_eq!(strip_foreign_decoy_tag("sp|P04618|REV_HV1H2"), None);
        assert_eq!(strip_foreign_decoy_tag("sp|Q9UBZ9|REV1_HUMAN"), None);
        assert_eq!(strip_foreign_decoy_tag("rev_"), None);
        assert_eq!(strip_foreign_decoy_tag("P1"), None);
    }

    #[test]
    fn a_target_only_database_has_no_foreign_decoys() {
        let db = make_db(&[
            ("sp|P1|A_HUMAN", b"MKTAYIAKQR"),
            ("sp|P2|B_HUMAN", b"MSEQNNTEMTFQIQR"),
            ("sp|P3|C_HUMAN", b"MKAKM"),
        ]);
        let found = ForeignDecoys::scan(&db);
        assert_eq!(found.n_reversed, 0);
        assert_eq!(found.n_named_paired + found.n_named_unpaired, 0);
        assert_eq!(found.verdict(), Ok(None));
    }

    #[test]
    fn reversed_records_under_an_unknown_name_are_an_error() {
        let db = make_db(&[
            ("sp|P1|A_HUMAN", b"MKTAYIAKQR"),
            ("sp|P2|B_HUMAN", b"MSEQNNTEMTFQIQR"),
            ("zz1", b"RQKAIYATKM"),
            ("zz2", b"RQIQFTMETNNQESM"),
        ]);
        let found = ForeignDecoys::scan(&db);
        assert_eq!(found.n_reversed, 4);
        assert_eq!(found.n_named_paired, 0);
        let err = found.verdict().unwrap_err();
        assert!(
            err.contains("4 of 4 records are the exact reverse"),
            "{err}"
        );
        assert!(err.contains("e.g. zz2"), "{err}");
    }

    #[test]
    fn short_reversed_sequences_are_not_evidence() {
        assert_eq!(min_reversible_len(4), 1);
        assert_eq!(min_reversible_len(20), 3);
        assert_eq!(min_reversible_len(63_778), 8);
        let filler: Vec<(String, Vec<u8>)> = (0..18u8)
            .map(|i| (format!("p{i}"), vec![b'M', b'A' + i, b'G', b'K', b'W']))
            .collect();
        let mut records: Vec<(&str, &[u8])> =
            filler.iter().map(|(a, s)| (&a[..], &s[..])).collect();
        records.push(("d1", b"EI"));
        records.push(("d2", b"IE"));
        records.push(("t1", b"MKTAYIAKQR"));
        records.push(("t2", b"MSEQNNTEMTFQIQR"));
        assert_eq!(ForeignDecoys::scan(&make_db(&records)).n_reversed, 0);
        *records.last_mut().unwrap() = ("t2", &b"RQKAIYATKM"[..]);
        assert_eq!(ForeignDecoys::scan(&make_db(&records)).n_reversed, 2);
    }

    #[test]
    fn a_record_named_as_the_decoy_of_another_is_an_error() {
        let db = make_db(&[
            ("sp|P1|A_HUMAN", b"MKTAYIAKQR"),
            ("sp|P2|B_HUMAN", b"MSEQNNTEMTFQIQR"),
            ("DECOY_sp|P1|A_HUMAN", b"MAKYTIKAQR"),
        ]);
        let found = ForeignDecoys::scan(&db);
        assert_eq!(found.n_reversed, 0);
        assert_eq!(found.n_named_paired, 1);
        let err = found.verdict().unwrap_err();
        assert!(
            err.contains("1 of 3 records are named as the decoy"),
            "{err}"
        );
        assert!(err.contains("DECOY_sp|P1|A_HUMAN"), "{err}");
    }

    #[test]
    fn weaker_evidence_is_a_warning() {
        let db = make_db(&[
            ("p1", b"MKTAYIAKQR"),
            ("p2", b"AKM"),
            ("p3", b"MKA"),
            ("p4", b"MSEQNNTEMTFQIQR"),
            ("p5", b"MLLPEPTIDEK"),
            ("p6", b"MGGSSTTWWYK"),
            ("REV_HV1H2", b"MAGRSGDSDEELLK"),
        ]);
        let found = ForeignDecoys::scan(&db);
        assert_eq!(found.n_reversed, 2);
        assert_eq!(found.n_named_unpaired, 1);
        let warning = found.verdict().unwrap().unwrap();
        assert!(
            warning.starts_with("WARN: the database may contain decoys"),
            "{warning}"
        );
        assert!(warning.contains("REV_HV1H2"), "{warning}");
    }

    fn decoys(db: &ProteinDb, strategy: DecoyStrategy, proteases: &Proteases) -> Vec<Vec<u8>> {
        decoy_db(db, "XXX", strategy, DEFAULT_DECOY_SEED, proteases)
            .proteins
            .into_iter()
            .map(|p| p.sequence)
            .collect()
    }

    fn reversed(db: &ProteinDb) -> Vec<Vec<u8>> {
        decoys(db, DecoyStrategy::Reverse, &Proteases::default())
    }

    fn enzyme(primary: Enzyme) -> Proteases {
        Proteases {
            primary,
            ..Proteases::default()
        }
    }

    #[test]
    fn mutate_changes_one_interior_residue_per_long_segment() {
        let db = make_db(&[("P1", b"PEPTIDELSAGRVLSK"), ("P2", b"MAGCTFYHRQGILW")]);
        let d = decoys(&db, DecoyStrategy::Mutate, &Proteases::default());
        let mass = |r: u8| model::AminoAcid::standard(r).unwrap().mass;
        for (t, s) in db.proteins.iter().zip(&d) {
            assert_eq!(t.sequence.len(), s.len());
            let cuts = compute_cleavage_positions(&t.sequence, Enzyme::Trypsin, &[]);
            for w in cuts.windows(2) {
                let (a, b) = (w[0] as usize, w[1] as usize);
                let diff: Vec<usize> = (a..b).filter(|&i| t.sequence[i] != s[i]).collect();
                if b - a < MUTATE_MIN_SEGMENT {
                    assert!(diff.is_empty());
                    continue;
                }
                assert_eq!(diff.len(), 1, "{:?}", std::str::from_utf8(s));
                let i = diff[0];
                assert!(i > a && i + 1 < b);
                assert!(!(t.sequence[0] == b'M' && a == 0 && i == 1));
                assert!(!b"KRP".contains(&s[i]));
                assert!((mass(s[i]) - mass(t.sequence[i])).abs() > 1.1);
            }
        }
        assert_eq!(d, decoys(&db, DecoyStrategy::Mutate, &Proteases::default()));
        assert_eq!(
            DecoyStrategy::from_name("mutate"),
            Some(DecoyStrategy::Mutate)
        );
    }

    #[test]
    fn mutate_gives_a_shared_segment_one_decoy() {
        let seg = b"AGLSTVEFK";
        let site = |position: u32| ModSite {
            position,
            mass_delta: 79.96633,
            anchor: SiteAnchor::Residue,
        };
        let entry = |acc: &str, open: usize| Protein {
            accession: acc.to_string(),
            description: String::new(),
            sequence: [b"WWK".as_slice(), seg].concat(),
            annotation: Some(Arc::new(ProteinAnnotation::new(
                (1..seg.len() - 1)
                    .filter(|&i| i != open && i != 4)
                    .map(|i| site(3 + i as u32))
                    .collect(),
                None,
            ))),
        };
        let db = ProteinDb {
            proteins: vec![entry("P1", 2), entry("P2", 6)],
        };
        let d = decoys(&db, DecoyStrategy::Mutate, &Proteases::default());
        let tail = |s: &[u8]| s[3..].to_vec();
        assert_eq!(tail(&d[0]), tail(&d[1]));
        let changed: Vec<usize> = (0..seg.len())
            .filter(|&i| tail(&d[0])[i] != seg[i])
            .collect();
        assert_eq!(changed, vec![4]);
    }

    #[test]
    fn refinement_mutate_decoys_are_reversed() {
        let db = make_db(&[("P1", b"PEPKTIDE"), ("P2", b"ACDEFGHIK")]);
        let d = whole_sequence_decoy_db(&db, "XXX", DecoyStrategy::Mutate, 7);
        let r = whole_sequence_decoy_db(&db, "XXX", DecoyStrategy::Reverse, 7);
        let seqs = |db: &ProteinDb| -> Vec<Vec<u8>> {
            db.proteins.iter().map(|p| p.sequence.clone()).collect()
        };
        assert_eq!(seqs(&d), seqs(&r));
        for (t, s) in db.proteins.iter().zip(&d.proteins) {
            assert_ne!(t.sequence, s.sequence);
        }
    }

    #[test]
    fn pseudo_reverses_each_tryptic_segment() {
        let db = make_db(&[("P1", b"PEPTIDEKVLSAGR"), ("P2", b"AGCTFYHRQGIL")]);
        assert_eq!(
            reversed(&db),
            [b"EDITPEPKGASLVR".to_vec(), b"HYFTCGARIGQL".to_vec()]
        );
    }

    #[test]
    fn initiator_met_stays() {
        let db = make_db(&[("P1", b"MAGPEKR")]);
        assert_eq!(reversed(&db), [b"MEPGAKR".to_vec()]);
    }

    #[test]
    fn n_terminal_enzyme_keeps_the_first_residue_of_each_segment() {
        let db = make_db(&[("P1", b"DAAGLKDPEW")]);
        let d = decoys(&db, DecoyStrategy::Reverse, &enzyme(Enzyme::AspN));
        assert_eq!(d, [b"DKLGAADWEP".to_vec()]);
    }

    #[test]
    fn non_specific_search_reverses_the_whole_sequence() {
        let db = make_db(&[("P1", b"MKWVLR")]);
        let want = [b"RLVWKM".to_vec()];
        assert_eq!(
            decoys(&db, DecoyStrategy::Reverse, &enzyme(Enzyme::NonSpecific)),
            want
        );
        let ntt0 = Proteases {
            num_tolerable_termini: 0,
            ..Proteases::default()
        };
        assert_eq!(decoys(&db, DecoyStrategy::Reverse, &ntt0), want);
        let semi = Proteases {
            num_tolerable_termini: 1,
            ..Proteases::default()
        };
        assert_eq!(
            decoys(&db, DecoyStrategy::Reverse, &semi),
            [b"MKLVWR".to_vec()]
        );
    }

    #[test]
    fn a_segment_has_one_decoy_wherever_it_occurs() {
        let db = make_db(&[("P1", b"GGHWKLLSAGNDKPPR"), ("P2", b"LLSAGNDKAEFRLLSAGNDK")]);
        for strategy in [
            DecoyStrategy::Reverse,
            DecoyStrategy::Shuffle,
            DecoyStrategy::Sequon,
        ] {
            let d = decoys(&db, strategy, &Proteases::default());
            let seg = d[0][5..13].to_vec();
            assert_eq!(seg[7], b'K', "{strategy:?}");
            assert_eq!(d[1][0..8], seg[..], "{strategy:?}");
            assert_eq!(d[1][12..20], seg[..], "{strategy:?}");
        }
    }

    #[test]
    fn palindromic_segment_is_redrawn() {
        let db = make_db(&[("P1", b"ADAKGLVLGR")]);
        let d = &reversed(&db)[0];
        assert_ne!(&d[0..4], b"ADAK");
        assert_eq!(d[3], b'K');
        let mut seg = d[0..4].to_vec();
        seg.sort_unstable();
        assert_eq!(seg, b"AADK");
        assert_ne!(&d[4..10], b"GLVLGR");
        assert_eq!(d[9], b'R');
    }

    #[test]
    fn decoy_segment_never_spells_another_target_segment() {
        let db = make_db(&[
            ("P1", b"GASVK"),
            ("P2", b"VSAGK"),
            ("P3", b"LIAGK"),
            ("P4", b"GALLK"),
        ]);
        let d = reversed(&db);
        let targets: Vec<Vec<u8>> = db
            .proteins
            .iter()
            .map(|p| {
                p.sequence
                    .iter()
                    .map(|&r| if r == b'I' { b'L' } else { r })
                    .collect()
            })
            .collect();
        for (t, s) in db.proteins.iter().zip(&d) {
            let folded: Vec<u8> = s
                .iter()
                .map(|&r| if r == b'I' { b'L' } else { r })
                .collect();
            assert!(!targets.contains(&folded), "{}", String::from_utf8_lossy(s));
            let (mut a, mut b) = (t.sequence.clone(), s.clone());
            a.sort_unstable();
            b.sort_unstable();
            assert_eq!(a, b);
            assert_eq!(s[4], b'K');
        }
    }

    fn random_db(count: usize) -> ProteinDb {
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let mut state = 0x1234_5678_9ABC_DEF1u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        ProteinDb {
            proteins: (0..count)
                .map(|i| {
                    let len = 30 + (next() % 300) as usize;
                    Protein {
                        accession: format!("P{i}"),
                        description: String::new(),
                        sequence: (0..len).map(|_| aa[(next() % 20) as usize]).collect(),
                        annotation: None,
                    }
                })
                .collect(),
        }
    }

    #[test]
    fn every_segment_keeps_composition_cuts_and_cleavage_side_residue() {
        let db = random_db(200);
        let digests = [
            enzyme(Enzyme::Trypsin),
            enzyme(Enzyme::LysC),
            enzyme(Enzyme::GluC),
            enzyme(Enzyme::Chymotrypsin),
            enzyme(Enzyme::AspN),
            enzyme(Enzyme::LysN),
            enzyme(Enzyme::NoCleavage),
            Proteases {
                primary: Enzyme::Trypsin,
                extra: vec![Enzyme::AspN],
                num_tolerable_termini: 2,
                ..Proteases::default()
            },
        ];
        for proteases in &digests {
            for strategy in [
                DecoyStrategy::Reverse,
                DecoyStrategy::Shuffle,
                DecoyStrategy::Sequon,
            ] {
                let d = decoys(&db, strategy, proteases);
                let mut changed = 0usize;
                let mut segments = 0usize;
                for (t, s) in db.proteins.iter().zip(&d) {
                    let cuts = compute_cleavage_positions(
                        &t.sequence,
                        proteases.primary,
                        &proteases.extra,
                    );
                    assert_eq!(
                        cuts,
                        compute_cleavage_positions(s, proteases.primary, &proteases.extra),
                        "{proteases:?} {strategy:?}"
                    );
                    for w in cuts.windows(2) {
                        let (a, b) = (w[0] as usize, w[1] as usize);
                        let (mut x, mut y) = (t.sequence[a..b].to_vec(), s[a..b].to_vec());
                        let anchor = if proteases.primary.is_n_term() {
                            a
                        } else {
                            b - 1
                        };
                        assert_eq!(t.sequence[anchor], s[anchor], "{proteases:?} {strategy:?}");
                        segments += 1;
                        changed += usize::from(x != y);
                        x.sort_unstable();
                        y.sort_unstable();
                        assert_eq!(x, y, "{proteases:?} {strategy:?}");
                    }
                    if t.sequence[0] == b'M' {
                        assert_eq!(s[0], b'M');
                    }
                }
                assert!(
                    changed * 2 > segments,
                    "{proteases:?} {strategy:?} {changed}/{segments}"
                );
            }
        }
    }

    #[test]
    fn shuffle_is_deterministic_for_a_fixed_seed() {
        let db = random_db(20);
        let p = Proteases::default();
        let a = decoy_db(&db, "XXX", DecoyStrategy::Shuffle, DEFAULT_DECOY_SEED, &p);
        let b = decoy_db(&db, "XXX", DecoyStrategy::Shuffle, DEFAULT_DECOY_SEED, &p);
        let c = decoy_db(
            &db,
            "XXX",
            DecoyStrategy::Shuffle,
            DEFAULT_DECOY_SEED + 1,
            &p,
        );
        let seqs = |d: &ProteinDb| -> Vec<Vec<u8>> {
            d.proteins.iter().map(|p| p.sequence.clone()).collect()
        };
        assert_eq!(seqs(&a), seqs(&b));
        assert_ne!(seqs(&a), seqs(&c));
        assert_ne!(
            seqs(&a),
            seqs(&decoy_db(&db, "XXX", DecoyStrategy::Reverse, 0, &p))
        );
    }

    #[test]
    fn decoy_db_prefixes_accessions() {
        let mut db = make_db(&[("P1", b"AB")]);
        db.proteins[0].description = "Some description".into();
        for prefix in ["XXX", "XXX_", ""] {
            let d = decoy_db(
                &db,
                prefix,
                DecoyStrategy::Reverse,
                0,
                &Proteases::default(),
            );
            assert_eq!(d.proteins[0].accession, "XXX_P1");
            assert_eq!(d.proteins[0].description, "Some description");
        }
    }

    #[test]
    fn build_search_db_appends_decoys_after_targets() {
        let target = make_db(&[("P1", b"PEPTIDEK"), ("P2", b"AGCTR")]);
        let built = build_search_db(&target, "XXX", DecoyStrategy::Reverse, DEFAULT_DECOY_SEED);
        assert_eq!(built.len(), 4);
        assert_eq!(built.proteins[0].accession, "P1");
        assert_eq!(built.proteins[2].accession, "XXX_P1");
        assert_eq!(built.proteins[2].sequence, b"EDITPEPK");
        assert_eq!(built.proteins[3].accession, "XXX_P2");
    }

    #[test]
    fn whole_sequence_sequon_shuffle_is_a_permutation_and_recovers_sequons() {
        let db = make_db(&[("P1", b"NISCDENLTKRSTNVTAAKNGSLLNPTVRQ")]);
        let orig = &db.proteins[0].sequence;
        let non_specific = enzyme(Enzyme::NonSpecific);
        let d = &decoys(&db, DecoyStrategy::Sequon, &non_specific)[0];
        let plain = &decoys(&db, DecoyStrategy::Shuffle, &non_specific)[0];
        let (mut a, mut b) = (orig.clone(), d.clone());
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "composition preserved");
        let (t, s, p) = (
            count_sequon_starts(orig),
            count_sequon_starts(d),
            count_sequon_starts(plain),
        );
        assert!(s <= t.max(p), "no overshoot ({s} <= {t})");
        assert!(s >= p, "at least the plain shuffle ({s} >= {p})");
    }

    #[test]
    fn peptide_level_sequon_shuffle_raises_the_sequon_count() {
        let db = random_db(200);
        let p = Proteases::default();
        let plain = decoys(&db, DecoyStrategy::Shuffle, &p);
        let sequon = decoys(&db, DecoyStrategy::Sequon, &p);
        let count = |seqs: &[Vec<u8>]| -> usize {
            db.proteins
                .iter()
                .zip(seqs)
                .map(|(target, s)| {
                    let cuts = compute_cleavage_positions(&target.sequence, p.primary, &p.extra);
                    cuts.windows(2)
                        .map(|w| count_sequon_starts(&s[w[0] as usize..w[1] as usize]))
                        .sum::<usize>()
                })
                .sum()
        };
        let targets: Vec<Vec<u8>> = db.proteins.iter().map(|p| p.sequence.clone()).collect();
        let (t, r, s) = (count(&targets), count(&plain), count(&sequon));
        assert!(s > r, "target {t}, shuffle {r}, sequon {s}");
    }

    #[test]
    fn build_search_db_none_is_target_only() {
        let target = make_db(&[("P1", b"ABCD"), ("P2", b"EFGH")]);
        let built = build_search_db(&target, "XXX", DecoyStrategy::None, DEFAULT_DECOY_SEED);
        assert_eq!(built.len(), 2, "None appends no decoys");
        assert_eq!(built.proteins[0].accession, "P1");
        assert_eq!(built.proteins[1].accession, "P2");
    }

    #[test]
    fn decoy_accession_needle_requires_underscore_delimiter() {
        assert_eq!(decoy_accession_needle("XXX"), "XXX_");
        assert_eq!(decoy_accession_needle("XXX_"), "XXX_");
        assert_eq!(decoy_accession_needle(""), "XXX_");
        assert!(is_decoy_accession("XXX_P12345", "XXX"));
        assert!(
            !is_decoy_accession("XXXP12345", "XXX"),
            "no '_' delimiter ⇒ target, not decoy"
        );
        assert!(!is_decoy_accession("P12345", "XXX"));
        assert!(is_decoy_accession("rev_sp|P1", "rev"));
        assert!(!is_decoy_accession("reverse_kinase", "rev"));
    }

    #[test]
    fn affix_matches_prefix_or_suffix() {
        assert!(is_decoy_accession_affix("XXX_P12345", "XXX", None));
        assert!(!is_decoy_accession_affix("P12345", "XXX", None));
        assert!(is_decoy_accession_affix(
            "sp|P02769|ALBU_BOVIN_rev",
            "XXX",
            Some("rev")
        ));
        assert!(!is_decoy_accession_affix(
            "sp|P02769|ALBU_BOVIN",
            "XXX",
            Some("rev")
        ));
        assert!(is_decoy_accession_affix("XXX_P1", "XXX", Some("rev")));
        assert!(!is_decoy_accession_affix("P1", "XXX", Some("")));
    }

    #[test]
    fn decoy_strategy_parses_names() {
        assert_eq!(
            DecoyStrategy::from_name("reverse"),
            Some(DecoyStrategy::Reverse)
        );
        assert_eq!(
            DecoyStrategy::from_name("Shuffle"),
            Some(DecoyStrategy::Shuffle)
        );
        assert_eq!(DecoyStrategy::from_name("NONE"), Some(DecoyStrategy::None));
        assert_eq!(DecoyStrategy::from_name("bogus"), None);
        assert_eq!(
            DecoyStrategy::from_name("sequon"),
            Some(DecoyStrategy::Sequon)
        );
        assert_eq!(DecoyStrategy::default(), DecoyStrategy::Shuffle);
    }

    fn digest(min_length: u32, max_length: u32) -> Proteases {
        Proteases {
            min_length,
            max_length,
            ..Proteases::default()
        }
    }

    fn folded(s: &[u8]) -> Vec<u8> {
        s.iter().map(|&r| if r == b'I' { b'L' } else { r }).collect()
    }

    fn peptides_of(seqs: &[Vec<u8>], proteases: &Proteases) -> FxHashSet<Vec<u8>> {
        let mut out = FxHashSet::default();
        for seq in seqs {
            let walk = seq.len() > 1 && seq[0] == b'M';
            collisions::for_each_span(seq, walk, proteases, |s, e| {
                out.insert(folded(&seq[s..e]));
            });
        }
        out
    }

    fn sequences(db: &ProteinDb) -> Vec<Vec<u8>> {
        db.proteins.iter().map(|p| p.sequence.clone()).collect()
    }

    #[test]
    fn a_missed_cleavage_join_of_short_segments_is_redrawn() {
        let db = make_db(&[("P1", b"ASKGVR"), ("P2", b"SAKVGR")]);
        let proteases = digest(6, 30);
        let segments = TargetSegments::new(&db, &proteases);
        let unchecked: Vec<Vec<u8>> = db
            .proteins
            .iter()
            .map(|p| {
                let source = peptide_level_permutation(
                    p,
                    DecoyStrategy::Reverse,
                    DEFAULT_DECOY_SEED,
                    &proteases,
                    &segments,
                );
                decoy_protein(p, "XXX", source, &proteases).sequence
            })
            .collect();
        assert_eq!(unchecked, [b"SAKVGR".to_vec(), b"ASKGVR".to_vec()]);
        let targets = peptides_of(&sequences(&db), &proteases);
        let d = decoys(&db, DecoyStrategy::Reverse, &proteases);
        for (t, s) in db.proteins.iter().zip(&d) {
            assert!(!targets.contains(&folded(s)), "{}", String::from_utf8_lossy(s));
            for (a, b) in [(0, 3), (3, 6)] {
                let (mut x, mut y) = (t.sequence[a..b].to_vec(), s[a..b].to_vec());
                assert_eq!(x[2], y[2]);
                x.sort_unstable();
                y.sort_unstable();
                assert_eq!(x, y);
            }
        }
    }

    #[test]
    fn a_non_specific_decoy_holds_no_target_peptide() {
        let fwd = b"ACDEFGHKLMNPQRSTVWYAGLK".to_vec();
        let rev: Vec<u8> = fwd.iter().rev().copied().collect();
        let db = make_db(&[("P1", &fwd), ("P2", &rev)]);
        let proteases = Proteases {
            primary: Enzyme::NonSpecific,
            num_tolerable_termini: 0,
            ..digest(8, 12)
        };
        let targets = peptides_of(&sequences(&db), &proteases);
        for strategy in [
            DecoyStrategy::Reverse,
            DecoyStrategy::Shuffle,
            DecoyStrategy::Sequon,
        ] {
            let d = decoys(&db, strategy, &proteases);
            let hits = peptides_of(&d, &proteases)
                .into_iter()
                .filter(|p| targets.contains(p))
                .count();
            assert_eq!(hits, 0, "{strategy:?}");
            for (t, s) in db.proteins.iter().zip(&d) {
                let (mut x, mut y) = (t.sequence.clone(), s.clone());
                x.sort_unstable();
                y.sort_unstable();
                assert_eq!(x, y, "{strategy:?}");
            }
        }
    }

    #[test]
    fn mutate_changes_joins_of_segments_shorter_than_seven() {
        let db = make_db(&[("P1", b"GASVKLLTGRDEFGHKW")]);
        let proteases = digest(6, 30);
        let targets = peptides_of(&sequences(&db), &proteases);
        let d = decoys(&db, DecoyStrategy::Mutate, &proteases);
        let hits: Vec<Vec<u8>> = peptides_of(&d, &proteases)
            .into_iter()
            .filter(|p| targets.contains(p))
            .collect();
        assert!(hits.is_empty(), "{hits:?}");
        let (t, s) = (&db.proteins[0].sequence, &d[0]);
        assert_eq!(
            compute_cleavage_positions(t, Enzyme::Trypsin, &[]),
            compute_cleavage_positions(s, Enzyme::Trypsin, &[])
        );
        let mass = |r: u8| model::AminoAcid::standard(r).unwrap().mass;
        for i in (0..t.len()).filter(|&i| t[i] != s[i]) {
            assert!((mass(s[i]) - mass(t[i])).abs() > 1.1);
            assert!(!b"KRP".contains(&s[i]) && !b"KR".contains(&t[i]));
        }
    }

    fn ecoli(count: usize) -> ProteinDb {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test-fixtures/ecoli.fasta");
        let file = std::fs::File::open(path).expect("open ecoli.fasta");
        let mut db = input::fasta::FastaReader::load_all(std::io::BufReader::new(file))
            .expect("read ecoli.fasta");
        db.proteins.truncate(count);
        db
    }

    fn search_peptides(
        db: &ProteinDb,
        proteases: &Proteases,
    ) -> Vec<(usize, usize, Vec<u8>)> {
        let aa = model::aa_set::AminoAcidSetBuilder::new_standard()
            .build()
            .expect("standard amino acids");
        let mut params = SearchParams::default_tryptic(aa);
        params.enzyme = proteases.primary;
        params.extra_enzymes = proteases.extra.clone();
        params.num_tolerable_termini = proteases.num_tolerable_termini;
        params.min_length = proteases.min_length;
        params.max_length = proteases.max_length;
        params.max_missed_cleavages = proteases.max_missed_cleavages;
        params.max_variable_mods_per_peptide = 0;
        let idx = crate::search_index::SearchIndex::from_target_db_for(
            db,
            "NONE",
            None,
            DecoyStrategy::None,
            DEFAULT_DECOY_SEED,
            proteases,
        );
        let mut out: Vec<(usize, usize, Vec<u8>)> =
            crate::candidate_gen::enumerate_candidates(&idx, &params, "NONE")
                .map(|c| {
                    let r: Vec<u8> = c.peptide.residues.iter().map(|a| a.residue).collect();
                    (c.protein_index, c.start_offset_in_protein, folded(&r))
                })
                .collect();
        out.sort();
        out.dedup();
        out
    }

    #[test]
    fn decoy_peptides_equal_to_targets_are_exactly_the_logged_remainder() {
        let db = ecoli(300);
        let tryptic = digest(7, 30);
        let semi = Proteases {
            num_tolerable_termini: 1,
            ..digest(7, 30)
        };
        let non_specific = Proteases {
            primary: Enzyme::NonSpecific,
            num_tolerable_termini: 0,
            ..digest(8, 12)
        };
        let permutations = [
            DecoyStrategy::Reverse,
            DecoyStrategy::Shuffle,
            DecoyStrategy::Sequon,
        ];
        let mut with_mutate = permutations.to_vec();
        with_mutate.push(DecoyStrategy::Mutate);
        for (label, proteases, strategies) in [
            ("tryptic", &tryptic, with_mutate),
            ("semi", &semi, permutations.to_vec()),
            ("non-specific", &non_specific, permutations.to_vec()),
        ] {
            let targets: FxHashSet<Vec<u8>> = search_peptides(&db, proteases)
                .into_iter()
                .map(|(_, _, r)| r)
                .collect();
            for strategy in strategies {
                let (proteins, repairs, _) =
                    build_decoys(&db, "XXX", strategy, DEFAULT_DECOY_SEED, proteases);
                let found: usize = repairs.iter().map(|r| r.found).sum();
                let left: usize = repairs.iter().map(|r| r.left).sum();
                let still = search_peptides(&ProteinDb { proteins }, proteases)
                    .into_iter()
                    .filter(|(_, _, r)| targets.contains(r))
                    .count();
                assert_eq!(still, left, "[{label} {strategy:?}] found {found}");
                assert!(left <= found, "[{label} {strategy:?}] {left} of {found} left");
                if label == "tryptic" && strategy == DecoyStrategy::Shuffle {
                    assert!(found > 0, "the fixture exercises no collision");
                }
            }
        }
    }
}
