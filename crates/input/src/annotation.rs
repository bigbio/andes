//! Entries built from [`Annotation`]s of a reference protein, and knockoff
//! annotations, whatever database the annotations come from.

use std::sync::Arc;

use model::{
    AminoAcid, Annotation, Edit, EntryKind, EntryOrigin, Enzyme, ModSite, Protein,
    ProteinAnnotation, Substitution, VariantRegion,
};

/// How knockoff variant entries are made. A knockoff is a real variant's
/// substitution moved to another position of the same sequence that holds the
/// same residue, in a matched context and far from every known variant: an
/// almost certainly false target that a mass coincidence or an isotope error
/// fits as well as it fits the real variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KnockoffOptions {
    /// Knockoffs per substitution, at most.
    pub per_variant: usize,
    /// The search's enzyme: a knockoff sits in another enzymatic span of a
    /// similar length.
    pub enzyme: Enzyme,
}

/// The entry `annotation` makes from `reference`, or `None` for a kind that
/// makes no entry (`Reference`). The entry carries the reference's sites that
/// its edits keep, plus the annotation's own; its origin names the
/// reference's family (for a reference with no origin, its primary
/// accession).
pub fn build_entry(reference: &Protein, annotation: &Annotation) -> Option<Protein> {
    let base = &reference.sequence;
    let family: Arc<str> = reference.origin().map_or_else(
        || Arc::from(primary_accession(&reference.accession)),
        |o| o.family.clone(),
    );
    let origin = |variant_delta_da: f64| EntryOrigin {
        kind: annotation.kind,
        family: family.clone(),
        feature: annotation.id.clone(),
        variant_delta_da,
        knockoff_of: annotation.knockoff_of.clone(),
    };
    let edits = triples(&annotation.edits);
    let (sequence, info) = match annotation.kind {
        EntryKind::Reference => return None,
        EntryKind::Variant => {
            let sequence = apply(base, &annotation.edits);
            let mut mod_sites = sites_through(reference.mod_sites(), &edits);
            mod_sites.extend(annotation.mods.iter().copied());
            let info = ProteinAnnotation::new(mod_sites, None)
                .with_variant(variant_regions(base, &annotation.edits))
                .with_origin(Some(origin(delta_da(base, &annotation.edits))))
                .with_substitutions(substitutions(base, &sequence, &annotation.edits));
            (sequence, info)
        }
        EntryKind::Isoform => {
            let sequence = apply(base, &annotation.edits);
            let mut mod_sites = sites_through(reference.mod_sites(), &edits);
            mod_sites.extend(annotation.mods.iter().copied());
            let info = ProteinAnnotation::new(mod_sites, None).with_origin(Some(origin(0.0)));
            (sequence, info)
        }
        EntryKind::Chain => {
            let c = annotation.processing?;
            let (start, end) = (c.start as usize, c.end as usize);
            let mut mod_sites: Vec<ModSite> = reference
                .mod_sites()
                .iter()
                .filter(|s| (start..end).contains(&(s.position as usize)))
                .map(|s| ModSite {
                    position: s.position - c.start,
                    ..*s
                })
                .collect();
            mod_sites.extend(annotation.mods.iter().copied());
            let info = ProteinAnnotation::new(mod_sites, None)
                .with_chain(Some(c.ends))
                .with_origin(Some(origin(0.0)));
            (base[start..end].to_vec(), info)
        }
    };
    let accession = annotation
        .accession
        .clone()
        .unwrap_or_else(|| suffixed_accession(&reference.accession, &annotation.id));
    let description = annotation.description.clone().unwrap_or_else(|| {
        let label = match (annotation.kind, annotation.processing) {
            (EntryKind::Chain, Some(c)) => {
                format!("CHAIN={}:{}-{}", annotation.id, c.start + 1, c.end)
            }
            (EntryKind::Isoform, _) => format!("ISOFORM={}", annotation.id),
            _ => format!("VARIANT={}", annotation.id),
        };
        format!("{} {label}", reference.description)
    });
    Some(Protein {
        accession,
        description,
        sequence,
        annotation: Some(Arc::new(info)),
    })
}

/// The primary accession of a database accession: the middle field of
/// `sp|P04637|P53_HUMAN`, the part after the database prefix of `sp:P04637`
/// (PEFF), else the whole accession.
pub fn primary_accession(accession: &str) -> &str {
    let parts: Vec<&str> = accession.split('|').collect();
    match parts.as_slice() {
        [_, acc, _] => acc,
        [one] => one.split_once(':').map_or(one, |(_, acc)| acc),
        _ => accession,
    }
}

/// `accession` with `_<id>` added to its middle field (`sp|P04637|P53_HUMAN`
/// → `sp|P04637_<id>|P53_HUMAN`), or appended when it has no such field.
pub fn suffixed_accession(accession: &str, id: &str) -> String {
    let parts: Vec<&str> = accession.split('|').collect();
    match parts.as_slice() {
        [db, acc, name] => format!("{db}|{acc}_{id}|{name}"),
        _ => format!("{accession}_{id}"),
    }
}

/// `base` with `edits` applied.
fn apply(base: &[u8], edits: &[Edit]) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(base.len());
    let mut at = 0;
    for e in edits {
        sequence.extend_from_slice(&base[at..e.start as usize]);
        sequence.extend_from_slice(&e.alt);
        at = e.end as usize;
    }
    sequence.extend_from_slice(&base[at..]);
    sequence
}

/// `(start, end, replacement length)` of each edit.
fn triples(edits: &[Edit]) -> Vec<(usize, usize, usize)> {
    edits
        .iter()
        .map(|e| (e.start as usize, e.end as usize, e.alt.len()))
        .collect()
}

/// The parts of the edited sequence that read differently from `base`, in its
/// coordinates. Each edit's parts come from the sequences, not the edit's
/// coordinates: an indel in a repeat can slide along it, spans past the shared
/// prefix or suffix would read like reference peptides, and so would one
/// inside a run the replacement keeps.
fn variant_regions(base: &[u8], edits: &[Edit]) -> Vec<VariantRegion> {
    let mut current = base.to_vec();
    let mut shift = 0i64;
    let mut regions = Vec::new();
    for e in edits {
        let start = (e.start as i64 + shift) as usize;
        let end = (e.end as i64 + shift) as usize;
        regions.extend(VariantRegion::of_edit(&current, start, end, &e.alt));
        current.splice(start..end, e.alt.iter().copied());
        shift += e.alt.len() as i64 - (e.end - e.start) as i64;
    }
    regions
}

/// Mass of the replacement residues minus the mass of those they replace.
fn delta_da(base: &[u8], edits: &[Edit]) -> f64 {
    let mass = |r: &[u8]| -> f64 {
        r.iter()
            .filter_map(|&a| AminoAcid::standard(a))
            .map(|a| a.mass)
            .sum()
    };
    edits
        .iter()
        .map(|e| mass(&e.alt) - mass(&base[e.start as usize..e.end as usize]))
        .sum()
}

/// For equal-length edits, each residue that differs and the one it replaced,
/// in the edited sequence's coordinates, so a form a searched mod makes
/// canonical is skipped. Empty when any edit changes the length.
fn substitutions(base: &[u8], sequence: &[u8], edits: &[Edit]) -> Vec<Substitution> {
    if edits
        .iter()
        .any(|e| e.alt.len() != (e.end - e.start) as usize)
    {
        return Vec::new();
    }
    edits
        .iter()
        .flat_map(|e| e.start as usize..e.end as usize)
        .filter(|&i| sequence[i] != base[i])
        .map(|i| Substitution {
            position: i as u32,
            original: base[i],
        })
        .collect()
}

/// Where residue `p` of the edited sequence's source ends up, `None` when an
/// edit replaced it.
fn map_through(edits: &[(usize, usize, usize)], p: usize) -> Option<usize> {
    let mut shift = 0i64;
    for &(start, end, len) in edits {
        if p < start {
            break;
        }
        if p < end {
            return None;
        }
        shift += len as i64 - (end - start) as i64;
    }
    Some((p as i64 + shift) as usize)
}

/// `sites` carried through `edits`: a site keeps its residue, and is dropped
/// when an edit replaced that residue.
fn sites_through(sites: &[ModSite], edits: &[(usize, usize, usize)]) -> Vec<ModSite> {
    sites
        .iter()
        .filter_map(|s| {
            map_through(edits, s.position as usize).map(|p| ModSite {
                position: p as u32,
                ..*s
            })
        })
        .collect()
}

/// The change an edit makes, 1-based: `D7H`, `KR120-121EE`, `F508del`.
pub(crate) fn change_label(original: &[u8], begin: u32, end: u32, replacement: &[u8]) -> String {
    let original = String::from_utf8_lossy(original);
    let range = if begin == end {
        begin.to_string()
    } else {
        format!("{begin}-{end}")
    };
    if replacement.is_empty() {
        format!("{original}{range}del")
    } else {
        format!("{original}{range}{}", String::from_utf8_lossy(replacement))
    }
}

/// Up to `ko.per_variant` knockoff annotations of each single-residue
/// substitution in `real`, annotations of `sequence`. Their positions keep
/// clear of every position in `listed` (every variant the source lists on the
/// sequence, whatever its evidence). `seed_key` (the accession) and the
/// variant's id seed the choice. A knockoff of `<id>` is named
/// `KO<i>_<change>_<id>`. Returns the knockoffs and how many substitutions had
/// no position for one.
pub fn knockoff_annotations(
    sequence: &[u8],
    seed_key: &str,
    real: &[Annotation],
    listed: &[usize],
    ko: KnockoffOptions,
) -> (Vec<Annotation>, usize) {
    let mut out = Vec::new();
    let mut untested = 0;
    for a in real
        .iter()
        .filter(|a| a.kind == EntryKind::Variant && a.is_single_substitution())
    {
        let e = &a.edits[0];
        let seed = fnv1a(format!("{seed_key}:{}", a.id).as_bytes());
        let positions = knockoff_positions(sequence, e.start as usize, listed, ko, seed);
        if positions.is_empty() {
            untested += 1;
            continue;
        }
        for (i, q) in positions.into_iter().enumerate() {
            let at = q as u32 + 1;
            let change = change_label(&sequence[q..q + 1], at, at, &e.alt);
            let mut k = Annotation::new(
                &format!("KO{}_{change}_{}", i + 1, a.id),
                EntryKind::Variant,
            );
            k.edits = vec![Edit {
                start: q as u32,
                end: q as u32 + 1,
                alt: e.alt.clone(),
            }];
            k.knockoff_of = Some(a.id.clone());
            out.push(k);
        }
    }
    (out, untested)
}

/// The class a knockoff site's next residue must share with the real site's:
/// hydrophobic, aromatic, polar, basic, acidic, or the residue itself.
fn residue_class(r: u8) -> u8 {
    match r {
        b'A' | b'V' | b'L' | b'I' | b'M' => b'h',
        b'F' | b'W' | b'Y' => b'a',
        b'S' | b'T' | b'N' | b'Q' => b'p',
        b'K' | b'R' | b'H' => b'+',
        b'D' | b'E' => b'-',
        other => other,
    }
}

/// The `[start, end)` span between the enzyme's cleavage sites around residue
/// `i` (no missed cleavage).
pub(crate) fn enzymatic_span(seq: &[u8], i: usize, enzyme: Enzyme) -> (usize, usize) {
    let cut = |j: usize| {
        j == 0
            || j == seq.len()
            || if enzyme.is_n_term() {
                enzyme.is_cleavable_before(seq[j])
            } else {
                enzyme.is_cleavable_after(seq[j - 1])
            }
    };
    let mut a = i;
    while !cut(a) {
        a -= 1;
    }
    let mut b = i + 1;
    while !cut(b) {
        b += 1;
    }
    (a, b)
}

/// Whether residue `i` is the N of an N-x-S/T sequon (x not P).
fn sequon(seq: &[u8], i: usize) -> bool {
    seq[i] == b'N' && i + 2 < seq.len() && seq[i + 1] != b'P' && matches!(seq[i + 2], b'S' | b'T')
}

/// Knockoff positions for a substitution of residue `p` of `seq`: positions
/// `q` holding the same residue, whose next residue matches (the same residue
/// when the site is N, Q, K or R or either next residue is P or G, else the
/// same [`residue_class`]), with the same sequon status, in another enzymatic
/// span whose length is within 3 of `p`'s and 6..=50, and at least 10 residues
/// from every position in `listed`. Up to `ko.per_variant` of them, drawn
/// without replacement by a generator seeded with `seed`.
pub(crate) fn knockoff_positions(
    seq: &[u8],
    p: usize,
    listed: &[usize],
    ko: KnockoffOptions,
    seed: u64,
) -> Vec<usize> {
    if p + 1 >= seq.len() {
        return Vec::new();
    }
    let o = seq[p];
    let (pa, pb) = enzymatic_span(seq, p, ko.enzyme);
    let next_ok = |q: usize| {
        let (a, b) = (seq[p + 1], seq[q + 1]);
        if b"NQKR".contains(&o) || matches!(a, b'P' | b'G') || matches!(b, b'P' | b'G') {
            a == b
        } else {
            residue_class(a) == residue_class(b)
        }
    };
    let mut eligible: Vec<usize> = (0..seq.len() - 1)
        .filter(|&q| {
            seq[q] == o
                && next_ok(q)
                && sequon(seq, q) == sequon(seq, p)
                && listed.iter().all(|&r| r.abs_diff(q) >= 10)
                && {
                    let (qa, qb) = enzymatic_span(seq, q, ko.enzyme);
                    qa != pa && (qb - qa).abs_diff(pb - pa) <= 3 && (6..=50).contains(&(qb - qa))
                }
        })
        .collect();
    let mut state = seed;
    let take = ko.per_variant.min(eligible.len());
    for i in 0..take {
        let j = i + (splitmix64(&mut state) % (eligible.len() - i) as u64) as usize;
        eligible.swap(i, j);
    }
    eligible.truncate(take);
    eligible
}

/// 64-bit FNV-1a hash.
fn fnv1a(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, &b| {
        (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// The next value of a SplitMix64 generator.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}
