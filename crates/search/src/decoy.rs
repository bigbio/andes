//! Decoy database generation.
//!
//! Decoys are peptide-level. Each target protein is cut into the segments
//! between the cleavage sites of the search's proteases (a missed-cleavage
//! peptide is a run of consecutive segments), and the residues of each
//! segment are permuted among themselves while three kinds of residue stay in
//! place: every residue that makes a cut (K and R for trypsin), the residue on
//! the segment's cleavage side (its last residue for an enzyme that cuts after
//! its residue, its first for AspN and LysN), and an initiator Met. Every
//! enzymatic peptide of a target therefore has a decoy peptide at the same
//! coordinates with the same composition, mass, length, cleavage-side residue
//! and missed cleavages, the way a peptide from a random wrong sequence of
//! the same digest would look. A segment's decoy depends only on its residues
//! (and on whether it starts with an initiator Met), so a peptide shared by
//! several entries has one decoy, and the decoys of overlapping
//! missed-cleavage spans agree with each other. Modification sites keep their
//! counts per peptide because residues only move within a peptide.
//!
//! A segment whose permutation reads as a target segment after folding I and
//! L (a palindrome, or a permutation that spells another protein's segment)
//! is redrawn with a shuffle seeded by its residues; a segment too short to
//! change keeps the target's residues, and candidate generation relabels such
//! a decoy peptide as a target.
//!
//! A non-specific search (protease `nonspecific`, or enzyme specificity
//! `non-specific`) treats every span as a peptide, so there are no segments
//! to permute within; its decoys are whole-sequence permutations.
//!
//! A target's database annotation (UniProt XML: annotated modification sites,
//! sequence-variant parts, mature-chain termini) is carried onto its decoy
//! through the same permutation that builds the decoy sequence
//! ([`decoy_annotation`]), so the decoy offers site-restricted modifications
//! on the same residues and, for a variant entry, variant parts of the same
//! shape. For an annotated target that permutation is first fitted
//! ([`fit_permutation`]): a terminal site keeps its residue at the decoy's
//! terminus, and a mature chain keeps its cleavage sites. The decoy side thus
//! mirrors the search space the annotation adds on the target side.

use std::sync::Arc;

use model::enzyme::Enzyme;
use model::protein::{ModSite, Protein, ProteinAnnotation, ProteinDb, SiteAnchor, VariantRegion};

use crate::candidate_gen::{compute_cleavage_positions, has_met_cleaved_walk};
use crate::search_params::SearchParams;

/// Default decoy accession prefix.
pub const DEFAULT_DECOY_PREFIX: &str = "XXX";

/// Default seed of the shuffles (the shuffle strategy, and the redraw of a
/// segment whose decoy reads as a target segment).
pub const DEFAULT_DECOY_SEED: u64 = 42;

/// Seeded redraws tried for a segment whose decoy reads as a target segment.
const MAX_REDRAWS: u64 = 20;

/// The digest of a search, which decides the segments a decoy permutes
/// within.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Proteases {
    pub primary: Enzyme,
    pub extra: Vec<Enzyme>,
    /// Termini of a span that must be cleavage sites: 2 fully specific, 1
    /// semi-specific, 0 non-specific (as
    /// [`SearchParams::num_tolerable_termini`]).
    pub num_tolerable_termini: u8,
}

impl Default for Proteases {
    /// Fully specific trypsin, the default search digest.
    fn default() -> Self {
        Self {
            primary: Enzyme::Trypsin,
            extra: Vec::new(),
            num_tolerable_termini: 2,
        }
    }
}

impl Proteases {
    /// The digest of `params`.
    pub fn of(params: &SearchParams) -> Self {
        Self {
            primary: params.enzyme,
            extra: params.extra_enzymes.clone(),
            num_tolerable_termini: params.num_tolerable_termini,
        }
    }

    fn all(&self) -> impl Iterator<Item = Enzyme> + '_ {
        std::iter::once(self.primary).chain(self.extra.iter().copied())
    }

    /// Whether every span is a peptide: a non-specific protease, or a search
    /// that requires no cleavage terminus.
    pub fn is_non_specific(&self) -> bool {
        self.num_tolerable_termini == 0 || self.all().any(|e| e == Enzyme::NonSpecific)
    }

    /// Whether residue `r` makes a cut under any of the proteases.
    fn makes_cut(&self, r: u8) -> bool {
        self.all()
            .any(|e| e.is_cleavable_after(r) || e.is_cleavable_before(r))
    }
}

/// How decoy proteins are generated for the target+decoy search database.
/// Each strategy permutes within cleavage segments (see the module
/// documentation), or the whole sequence for a non-specific search.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecoyStrategy {
    /// Reverse each segment's movable residues (the pseudo-reversed peptide:
    /// `PEPTIDEK` → `EDITPEPK`). Default.
    #[default]
    Reverse,
    /// Shuffle each segment's movable residues with a shuffle seeded by the
    /// segment's residues and the decoy seed.
    Shuffle,
    /// [`Self::Reverse`], then swaps among the movable residues that raise
    /// the segment's N-X-S/T sequon count toward (never above) the target
    /// segment's, so the glyco sequon gate admits decoys at the target rate.
    /// Reversal alone turns N-X-S/T into S/T-X-N and leaves fewer sequons.
    SequonReverse,
    /// Generate no decoys: the search database is the input FASTA verbatim. For
    /// inputs that already contain decoys (accessions carrying the decoy prefix)
    /// or when FDR is computed externally downstream.
    None,
}

impl DecoyStrategy {
    /// Parse the `--decoy-strategy` CLI value. Case-insensitive.
    pub fn from_name(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "reverse" | "rev" => Some(Self::Reverse),
            "shuffle" | "shuf" => Some(Self::Shuffle),
            "sequon-reverse" | "sequon" | "sequon_reverse" => Some(Self::SequonReverse),
            "none" | "off" => Some(Self::None),
            _ => None,
        }
    }
}

/// The decoy proteins of `target` under `strategy` for a search with
/// `proteases`, in target order, each accession `<prefix>_<target
/// accession>` (`prefix` normalized by [`normalize_decoy_prefix`]). Empty for
/// [`DecoyStrategy::None`]. `seed` drives the shuffles.
pub fn decoy_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
) -> ProteinDb {
    let prefix = normalize_decoy_prefix(prefix);
    let proteins = match strategy {
        DecoyStrategy::None => Vec::new(),
        _ if proteases.is_non_specific() => {
            whole_sequence_decoys(target, &prefix, strategy, seed, proteases)
        }
        _ => peptide_level_decoys(target, &prefix, strategy, seed, proteases),
    };
    ProteinDb { proteins }
}

/// Decoys that permute each target sequence as a whole: reversal, a shuffle
/// seeded per protein, or a sequon-preserving reversal. The decoys of a
/// non-specific search, and of sequences that are a single peptide each.
pub(crate) fn whole_sequence_decoy_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
) -> ProteinDb {
    let prefix = normalize_decoy_prefix(prefix);
    let proteins = match strategy {
        DecoyStrategy::None => Vec::new(),
        _ => whole_sequence_decoys(target, &prefix, strategy, seed, &Proteases::default()),
    };
    ProteinDb { proteins }
}

/// The decoy of `p` built by the permutation `source` (`source[i]`: the
/// target position decoy residue `i` comes from).
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

/// Runs on the calling thread: the search configures the global thread pool
/// after it has built the database.
fn peptide_level_decoys(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
) -> Vec<Protein> {
    let segments = TargetSegments::new(target, proteases);
    target
        .proteins
        .iter()
        .map(|p| {
            let source = peptide_level_permutation(p, strategy, seed, proteases, &segments);
            decoy_protein(p, prefix, source, proteases)
        })
        .collect()
}

fn whole_sequence_decoys(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
) -> Vec<Protein> {
    target
        .proteins
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let n = p.sequence.len() as u32;
            let source: Vec<u32> = match strategy {
                DecoyStrategy::Shuffle => {
                    // Per-protein stream: the seed mixed with the index, so
                    // identical sequences do not shuffle alike.
                    let mut source: Vec<u32> = (0..n).collect();
                    shuffle_in_place(
                        &mut source,
                        seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15),
                    );
                    source
                }
                DecoyStrategy::SequonReverse => {
                    let mut source: Vec<u32> = (0..n).rev().collect();
                    restore_sequons(&p.sequence, &mut source, &vec![true; n as usize]);
                    source
                }
                _ => (0..n).rev().collect(),
            };
            decoy_protein(p, prefix, source, proteases)
        })
        .collect()
}

/// The peptide-level permutation of target `p`: each cleavage segment
/// permuted by [`segment_permutation`].
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

/// The decoy of cleavage segment `seg` as a permutation of its positions
/// (`local[i]`: the position decoy residue `i` comes from). `lead_met`: the
/// segment starts with an initiator Met, which stays. The residues that make
/// a cut and the cleavage-side residue stay too; the others are reversed
/// (or shuffled), redrawn by seeded shuffles while the result reads as a
/// target segment.
fn segment_permutation(
    seg: &[u8],
    lead_met: bool,
    strategy: DecoyStrategy,
    seed: u64,
    proteases: &Proteases,
    targets: &TargetSegments,
) -> Vec<u32> {
    let m = seg.len();
    let identity: Vec<u32> = (0..m as u32).collect();
    if m == 0 {
        return identity;
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
        return identity;
    }
    let draw = |k: u64| -> Vec<u32> {
        let mut order = free.clone();
        if k == 0 && strategy != DecoyStrategy::Shuffle {
            order.reverse();
        } else {
            shuffle_in_place(&mut order, segment_seed(seed, seg, lead_met, k));
        }
        let mut local = identity.clone();
        for (&at, &from) in free.iter().zip(&order) {
            local[at as usize] = from;
        }
        if strategy == DecoyStrategy::SequonReverse {
            restore_sequons(seg, &mut local, &movable);
        }
        local
    };
    let first = draw(0);
    if !targets.collides(seg, &first, lead_met) {
        return first;
    }
    (1..=MAX_REDRAWS)
        .map(draw)
        .find(|local| !targets.collides(seg, local, lead_met))
        .unwrap_or(first)
}

/// Sorted hashes of every cleavage segment of the targets with I and L
/// folded, and of each initiator-Met segment without its Met: the residue
/// strings a decoy segment must not spell.
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

    /// Whether the decoy `local` makes of `seg` reads as a target segment,
    /// also after removing a leading initiator Met.
    fn collides(&self, seg: &[u8], local: &[u32], lead_met: bool) -> bool {
        let residues = || local.iter().map(|&i| seg[i as usize]);
        self.contains(residues()) || (lead_met && self.contains(residues().skip(1)))
    }
}

/// FNV-1a hash of `residues` with I read as L.
fn folded_hash(residues: impl Iterator<Item = u8>) -> u64 {
    let mut h: u64 = 0xCBF2_9CE4_8422_2325;
    for r in residues {
        h ^= u64::from(if r == b'I' { b'L' } else { r });
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// Shuffle seed of draw `k` of segment `seg`: depends only on the residues,
/// so identical segments shuffle alike wherever they occur.
fn segment_seed(seed: u64, seg: &[u8], lead_met: bool, k: u64) -> u64 {
    let mut h: u64 = 0xCBF2_9CE4_8422_2325;
    for &r in seg {
        h ^= u64::from(r);
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    splitmix64(h ^ splitmix64(seed ^ splitmix64(k.wrapping_mul(2) | u64::from(lead_met))))
}

/// Fisher-Yates from the back, drawing each index from a xorshift64 stream
/// started at `state` (never zero, where xorshift would stick).
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

/// The decoy of the annotated target `p` under the strategy's permutation
/// `source` (`source[i]`: the target position decoy residue `i` comes from),
/// fitted to the annotation first ([`fit_permutation`]).
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

/// Fit the strategy's permutation `source` of an annotated target so that its
/// decoy offers the search space the annotation adds:
///
/// - A mature-chain entry yields only its terminal spans, so its decoy keeps
///   their cleavage structure. Residues are permuted only within the
///   cleavage segments of `proteases`, in the order the strategy gave them
///   (reversed where that is no change), and every residue that makes a cut
///   stays in place. Cleavage rules look only at the residues around a cut,
///   so the decoy has the same cut sites, and each terminal span of the
///   target has a decoy span with the same length, mass and missed
///   cleavages. A peptide-level permutation already has that form and is
///   kept. A non-specific digest cuts everywhere, so any permutation keeps
///   that.
/// - A residue carrying a site tied to the protein terminus where the target
///   has it (the first residue; the second after an initiator Met, which then
///   stays too; the last residue) stays at the decoy's terminus. A
///   residue-specific terminal mod restricted to annotated sites (pyro-Glu on
///   Q) can then be placed on the decoy as on the target.
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
            // A shuffle can leave a short segment as it was, which would make
            // the decoy peptide the target's; reverse it then.
            if order == free {
                order.reverse();
            }
            for (&at, &from) in free.iter().zip(&order) {
                source[at as usize] = from;
            }
        }
    }
    // Keep target residue `q` at decoy position `q`. Within a chain both
    // positions lie in the first or last segment, so the cuts stay too.
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

/// The annotation of the decoy built from `target` by a permutation of its
/// residues: `source[i]` is the target position decoy residue `i` came from.
///
/// Modification sites follow their residue; a mature chain keeps its
/// new-terminus flags. Builders fit the permutation first
/// ([`fit_permutation`]), so a terminal site's residue is already at the
/// decoy's terminus and a chain's terminal spans keep their cleavage sites.
///
/// Each variant part maps to its image when that is contiguous (always,
/// under reversal); otherwise (shuffling scatters it) to a part of the same
/// length at the image of its first residue, moved back from the end if
/// needed. The same holds for the stretch a deletion junction or an indel in
/// a repeat requires a span to hold (see [`VariantRegion`]). Either way the
/// decoy yields variant spans of the same shape as the target, which is what
/// keeps the two sides' search spaces matched. `None` when the target has no
/// annotation.
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

/// The decoy image of variant part `v` under the permutation whose inverse
/// is `image` (see [`decoy_annotation`]).
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
        // A covering span holds all of [b - 1, a]: the two residues of a
        // deletion junction (a == b) or a repeat stretch (a > b). Keep that
        // requirement on the image, or on a stretch of the same length at
        // the image of its first residue when shuffling scattered it.
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

/// Count N-X-S/T sequon start positions (X ≠ P) in a residue slice.
fn count_sequon_starts(s: &[u8]) -> usize {
    (0..s.len().saturating_sub(2))
        .filter(|&i| s[i] == b'N' && s[i + 1] != b'P' && (s[i + 2] == b'S' || s[i + 2] == b'T'))
        .count()
}

/// Raise the N-X-S/T sequon count of the decoy that `source` makes of `orig`
/// (`source[i]`: the position of `orig` decoy residue `i` comes from) toward,
/// never above, the count in `orig`. Only positions with `movable` set take
/// part in a swap.
///
/// Reversal does not drive the sequon count to zero: on a human proteome it
/// leaves a different set of sequons at ~87% of the target's density, a
/// deficit that makes the glyco sequon gate admit decoys slightly less often
/// than targets. Restoring every original sequon on top of that set
/// over-corrects (~160%). This closes only the deficit: scanning windows in
/// order, it pulls a later S/T into the third slot of an `N-X-·` window. A
/// swap can also break or create a sequon at the donor's window, so each swap
/// is kept only if the recounted total rises and stays within the target's.
/// The result is still a permutation of `orig`. Deterministic.
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
        // A window that is one S/T away from being a sequon: N, X≠P, non-S/T.
        if movable[i + 2]
            && r[i] == b'N'
            && r[i + 1] != b'P'
            && r[i + 2] != b'S'
            && r[i + 2] != b'T'
        {
            // Donor: a later movable S/T that is not the third residue of a
            // sequon (moving it would break that sequon, a net-zero swap).
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

/// The searched database: `target` followed by its decoys ([`decoy_db`]) for
/// a fully specific trypsin search; the target alone for
/// [`DecoyStrategy::None`]. Decoy membership downstream is decided by the
/// accession prefix, not the position.
pub fn build_search_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
) -> ProteinDb {
    build_search_db_for(target, prefix, strategy, seed, &Proteases::default())
}

/// [`build_search_db`] for a search with `proteases`, whose cleavage
/// segments the decoys permute within.
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

/// Normalize a user-supplied decoy accession prefix: trim whitespace,
/// strip trailing `_`, and fall back to [`DEFAULT_DECOY_PREFIX`] when empty.
pub fn normalize_decoy_prefix(prefix: &str) -> String {
    let trimmed = prefix.trim().trim_end_matches('_');
    if trimmed.is_empty() {
        DEFAULT_DECOY_PREFIX.to_string()
    } else {
        trimmed.to_string()
    }
}

/// The decoy-membership needle for `prefix`: `"<normalized_prefix>_"`.
///
/// SINGLE SOURCE OF TRUTH for target/decoy membership. Decoy accessions are
/// built as `"<prefix>_<orig>"` (see [`decoy_db`]), so a protein
/// is a decoy iff its accession starts with this needle. Both candidate
/// generation and [`crate::SearchIndex`] MUST use this helper so the two paths
/// can never diverge: a target accession that merely starts with the BARE prefix
/// (e.g. a real protein named `"XXXfoo"`, or any accession under a short custom
/// prefix like `rev`) is NOT a decoy — the `_` delimiter is required.
pub fn decoy_accession_needle(prefix: &str) -> String {
    format!("{}_", normalize_decoy_prefix(prefix))
}

/// Returns `true` iff `accession` is a decoy accession under `prefix`
/// (i.e. starts with [`decoy_accession_needle`]).
pub fn is_decoy_accession(accession: &str, prefix: &str) -> bool {
    accession.starts_with(&decoy_accession_needle(prefix))
}

/// Returns `true` iff `accession` is a decoy under the AFFIX convention: it
/// starts with the prefix needle `"<prefix>_"` (andes' native decoy naming) OR
/// — when `suffix` is set and non-empty — ends with `suffix`.
///
/// The suffix form lets andes recognize decoys produced by an external pipeline
/// whose pre-built target+decoy DB tags decoys with an accession SUFFIX (e.g.
/// quantms: `<orig>_rev`) rather than andes' prefix. Used
/// with `--decoy-strategy none` so andes consumes the shared decoy DB verbatim
/// instead of generating its own (which would double-decoy and bias FDR).
pub fn is_decoy_accession_affix(accession: &str, prefix: &str, suffix: Option<&str>) -> bool {
    accession.starts_with(&decoy_accession_needle(prefix))
        || suffix.is_some_and(|s| !s.is_empty() && accession.ends_with(s))
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
        // AspN cuts before D: segments DAAGLK and DPEW keep their D.
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
        // Semi-specific searches keep the segments.
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
        // LLSAGNDK in two contexts, once inside a missed-cleavage stretch.
        let db = make_db(&[("P1", b"GGHWKLLSAGNDKPPR"), ("P2", b"LLSAGNDKAEFRLLSAGNDK")]);
        for strategy in [
            DecoyStrategy::Reverse,
            DecoyStrategy::Shuffle,
            DecoyStrategy::SequonReverse,
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
        // GLVLG is a palindrome too; the redraw keeps the R.
        assert_ne!(&d[4..10], b"GLVLGR");
        assert_eq!(d[9], b'R');
    }

    #[test]
    fn decoy_segment_never_spells_another_target_segment() {
        // Reversing GASV gives VSAG, P2's segment; I/L folded: LIAGK reversed
        // is GAILK, which reads as P3's GALLK.
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

    /// Random proteins over the 20 residues (deterministic xorshift).
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
            },
        ];
        for proteases in &digests {
            for strategy in [
                DecoyStrategy::Reverse,
                DecoyStrategy::Shuffle,
                DecoyStrategy::SequonReverse,
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
                // Short segments cannot change; most others do.
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
    fn whole_sequence_sequon_reverse_is_a_permutation_and_recovers_sequons() {
        let db = make_db(&[("P1", b"NISCDENLTKRSTNVTAAK")]);
        let orig = &db.proteins[0].sequence;
        let d = &decoys(
            &db,
            DecoyStrategy::SequonReverse,
            &enzyme(Enzyme::NonSpecific),
        )[0];
        let (mut a, mut b) = (orig.clone(), d.clone());
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "composition preserved");
        let plain: Vec<u8> = orig.iter().rev().copied().collect();
        let (t, s, p) = (
            count_sequon_starts(orig),
            count_sequon_starts(d),
            count_sequon_starts(&plain),
        );
        assert!(s <= t, "no overshoot ({s} <= {t})");
        assert!(s >= p, "at least plain reversal ({s} >= {p})");
    }

    #[test]
    fn peptide_level_sequon_reverse_keeps_segment_sequon_counts_within_target() {
        let db = random_db(200);
        let p = Proteases::default();
        let plain = reversed(&db);
        let sequon = decoys(&db, DecoyStrategy::SequonReverse, &p);
        let (mut t, mut r, mut s) = (0, 0, 0);
        for ((target, rev), seq) in db.proteins.iter().zip(&plain).zip(&sequon) {
            let cuts = compute_cleavage_positions(&target.sequence, p.primary, &p.extra);
            for w in cuts.windows(2) {
                let (a, b) = (w[0] as usize, w[1] as usize);
                let (ct, cs) = (
                    count_sequon_starts(&target.sequence[a..b]),
                    count_sequon_starts(&seq[a..b]),
                );
                assert!(cs <= ct.max(count_sequon_starts(&rev[a..b])));
                t += ct;
                r += count_sequon_starts(&rev[a..b]);
                s += cs;
            }
        }
        assert!(s > r, "target {t}, reverse {r}, sequon-reverse {s}");
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
        // Decoys are built as "<prefix>_<orig>", so detection requires the "_".
        assert_eq!(decoy_accession_needle("XXX"), "XXX_");
        assert_eq!(decoy_accession_needle("XXX_"), "XXX_"); // normalize is idempotent
        assert_eq!(decoy_accession_needle(""), "XXX_"); // empty → default prefix
        assert!(is_decoy_accession("XXX_P12345", "XXX"));
        assert!(
            !is_decoy_accession("XXXP12345", "XXX"),
            "no '_' delimiter ⇒ target, not decoy"
        );
        assert!(!is_decoy_accession("P12345", "XXX"));
        // Short custom prefix: a real "rev..." target must NOT be misread as decoy.
        assert!(is_decoy_accession("rev_sp|P1", "rev"));
        assert!(!is_decoy_accession("reverse_kinase", "rev"));
    }

    #[test]
    fn affix_matches_prefix_or_suffix() {
        // Prefix (andes-native) still works with no suffix configured.
        assert!(is_decoy_accession_affix("XXX_P12345", "XXX", None));
        assert!(!is_decoy_accession_affix("P12345", "XXX", None));
        // Suffix form recognizes quantms "<orig>_rev" decoys that the
        // prefix test alone misses.
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
        // Either affix qualifies; an empty suffix is ignored (no spurious match).
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
        assert_eq!(DecoyStrategy::default(), DecoyStrategy::Reverse);
    }
}
