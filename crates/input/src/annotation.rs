//! Entries built from [`Annotation`]s of a reference protein, and knockoff
//! annotations, whatever database the annotations come from.

use std::sync::Arc;

use model::{
    AminoAcid, Annotation, Edit, EntryKind, EntryOrigin, Enzyme, ModSite, Protein,
    ProteinAnnotation, ProteinDb, SiteAnchor, SiteClaim, Substitution, VariantRegion,
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
        source: annotation.source.clone(),
        tier: annotation.tier,
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
            k.source = a.source.clone();
            k.tier = a.tier;
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

/// Counts from [`add_site_knockoffs`], for the run log.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SiteKnockoffStats {
    /// Real site claims considered.
    pub sites: usize,
    /// Knockoff sites added.
    pub knockoffs: usize,
    /// Real sites for which no position qualified (or terminal sites).
    pub untested: usize,
}

const PHOSPHO_DA: f64 = 79.966331;
const METHYL_DA: f64 = 14.01565;
const DIMETHYL_DA: f64 = 28.0313;
const DEAMIDATION_DA: f64 = 0.984016;
const OXIDATION_DA: f64 = 15.994915;

/// The motif status a knockoff site must share with its real site, for the
/// modification `mass` on residue `seq[q]`: phospho S/T (proline at +1, a
/// basic residue at −3, an acidic one at +3), phospho Y (an acidic residue at
/// −1..−4), methyl or dimethyl R, deamidated N and hydroxy P (glycine at +1),
/// and N's sequon. Any other modification has no motif status.
fn motif_status(seq: &[u8], q: usize, mass: f64) -> [bool; 3] {
    let at = |i: isize| -> u8 {
        let j = q as isize + i;
        if j >= 0 && (j as usize) < seq.len() {
            seq[j as usize]
        } else {
            0
        }
    };
    let is = |m: f64| (mass - m).abs() <= model::protein::SITE_MASS_TOLERANCE;
    match seq[q] {
        b'S' | b'T' if is(PHOSPHO_DA) => [
            at(1) == b'P',
            matches!(at(-3), b'R' | b'K'),
            matches!(at(3), b'D' | b'E'),
        ],
        b'Y' if is(PHOSPHO_DA) => [(1..=4).any(|i| matches!(at(-i), b'D' | b'E')), false, false],
        b'R' if is(METHYL_DA) || is(DIMETHYL_DA) => [at(1) == b'G', false, false],
        b'N' if is(DEAMIDATION_DA) => [at(1) == b'G', sequon(seq, q), false],
        b'N' => [sequon(seq, q), false, false],
        b'P' if is(OXIDATION_DA) => [at(1) == b'G', false, false],
        _ => [false; 3],
    }
}

/// The span a peptide modified at `q` lies in: the enzymatic span, merged
/// with the next one for a modified K or R, after which the enzyme does not
/// cut.
fn modified_span(seq: &[u8], q: usize, enzyme: Enzyme) -> (usize, usize) {
    let (a, b) = enzymatic_span(seq, q, enzyme);
    if matches!(seq[q], b'K' | b'R') && b == q + 1 && b < seq.len() {
        (a, enzymatic_span(seq, b, enzyme).1)
    } else {
        (a, b)
    }
}

/// Add up to `ko.per_variant` knockoff sites for each real site claim of every
/// reference protein in `db`. A knockoff site holds the same residue in a
/// matched context: the same digestion (a modified K or R blocks cleavage,
/// so the merged span is compared), the next-residue rule of variant
/// knockoffs, the same motif status ([`motif_status`]), another span whose
/// length is within 3 of the real one's and 6..=50, and at least 5 residues
/// from every real site of the same modification. Terminal sites get none.
/// A knockoff site is also a `mod_sites` entry, so a `sites=annotated` mod may
/// be placed there, and a site claim naming the real site it copies.
pub fn add_site_knockoffs(db: &mut ProteinDb, ko: KnockoffOptions) -> SiteKnockoffStats {
    let mut stats = SiteKnockoffStats::default();
    if ko.per_variant == 0 {
        return stats;
    }
    for p in db.proteins.iter_mut() {
        if p.entry_kind() != EntryKind::Reference {
            continue;
        }
        let Some(info) = p.annotation.as_deref() else {
            continue;
        };
        let real: Vec<SiteClaim> = info
            .site_claims
            .iter()
            .filter(|c| c.knockoff_of.is_none())
            .copied()
            .collect();
        if real.is_empty() {
            continue;
        }
        let seq = &p.sequence;
        let same = |a: f64, b: f64| (a - b).abs() <= model::protein::SITE_MASS_TOLERANCE;
        let mut sites = info.mod_sites.clone();
        let mut claims = info.site_claims.clone();
        for c in &real {
            stats.sites += 1;
            let pos = c.position as usize;
            let terminal = info.mod_sites.iter().any(|s| {
                s.position == c.position
                    && same(s.mass_delta, c.mass_delta)
                    && s.anchor != SiteAnchor::Residue
            });
            if terminal || pos + 1 >= seq.len() {
                stats.untested += 1;
                continue;
            }
            let (pa, pb) = modified_span(seq, pos, ko.enzyme);
            let status = motif_status(seq, pos, c.mass_delta);
            let next_ok = |q: usize| {
                let (a, b) = (seq[pos + 1], seq[q + 1]);
                if matches!(a, b'P' | b'G') || matches!(b, b'P' | b'G') {
                    a == b
                } else {
                    residue_class(a) == residue_class(b)
                }
            };
            let mut eligible: Vec<usize> = (0..seq.len() - 1)
                .filter(|&q| {
                    seq[q] == seq[pos]
                        && next_ok(q)
                        && motif_status(seq, q, c.mass_delta) == status
                        && real
                            .iter()
                            .filter(|r| same(r.mass_delta, c.mass_delta))
                            .all(|r| (r.position as usize).abs_diff(q) >= 5)
                        && {
                            let (qa, qb) = modified_span(seq, q, ko.enzyme);
                            qa != pa
                                && (qb - qa).abs_diff(pb - pa) <= 3
                                && (6..=50).contains(&(qb - qa))
                        }
                })
                .collect();
            let seed = fnv1a(
                format!(
                    "{}:{}:{}",
                    primary_accession(&p.accession),
                    pos,
                    c.mass_delta.to_bits()
                )
                .as_bytes(),
            );
            let mut state = seed;
            let take = ko.per_variant.min(eligible.len());
            for i in 0..take {
                let j = i + (splitmix64(&mut state) % (eligible.len() - i) as u64) as usize;
                eligible.swap(i, j);
            }
            eligible.truncate(take);
            if eligible.is_empty() {
                stats.untested += 1;
                continue;
            }
            for q in eligible {
                sites.push(ModSite {
                    position: q as u32,
                    mass_delta: c.mass_delta,
                    anchor: SiteAnchor::Residue,
                });
                claims.push(SiteClaim {
                    position: q as u32,
                    mass_delta: c.mass_delta,
                    tier: c.tier,
                    knockoff_of: Some(c.position),
                });
                stats.knockoffs += 1;
            }
        }
        let mut updated = info.clone();
        updated.mod_sites = ProteinAnnotation::new(sites, None).mod_sites;
        updated = updated.with_site_claims(claims);
        p.annotation = Some(Arc::new(updated));
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::PriorTier;

    const KO: KnockoffOptions = KnockoffOptions {
        per_variant: 4,
        enzyme: Enzyme::Trypsin,
    };

    /// A reference protein whose phospho site (S at 0-based `site`) is
    /// annotated, with `seq`.
    fn reference(seq: &[u8], site: u32) -> ProteinDb {
        let sites = vec![ModSite {
            position: site,
            mass_delta: PHOSPHO_DA,
            anchor: SiteAnchor::Residue,
        }];
        let claims = vec![SiteClaim {
            position: site,
            mass_delta: PHOSPHO_DA,
            tier: PriorTier::Curated,
            knockoff_of: None,
        }];
        ProteinDb {
            proteins: vec![Protein {
                accession: "sp|P1|A".into(),
                description: String::new(),
                sequence: seq.to_vec(),
                annotation: Some(Arc::new(
                    ProteinAnnotation::new(sites, None).with_site_claims(claims),
                )),
            }],
        }
    }

    #[test]
    fn knockoff_sites_match_the_motif_and_stay_clear_of_real_sites() {
        // Spans GGGGSPGGGK (S then P: proline-directed) alternate with
        // GGGGSAGGGK (S then A). The real site is the S at 4 (S-P).
        let mut seq = Vec::new();
        for _ in 0..4 {
            seq.extend(b"GGGGSPGGGK");
            seq.extend(b"GGGGSAGGGK");
        }
        let mut db = reference(&seq, 4);
        let stats = add_site_knockoffs(&mut db, KO);
        assert_eq!((stats.sites, stats.knockoffs, stats.untested), (1, 3, 0));
        let info = db.proteins[0].annotation.as_deref().unwrap();
        let ko: Vec<&SiteClaim> = info
            .site_claims
            .iter()
            .filter(|c| c.knockoff_of.is_some())
            .collect();
        assert_eq!(ko.len(), 3);
        for c in ko {
            let q = c.position as usize;
            assert_eq!(&seq[q..q + 2], b"SP", "same motif: proline at +1");
            assert!(q.abs_diff(4) >= 5);
            assert_eq!(c.knockoff_of, Some(4));
            assert!(
                model::protein::has_site(&info.mod_sites, c.position, PHOSPHO_DA),
                "placeable"
            );
        }
    }

    #[test]
    fn a_terminal_site_or_one_with_no_match_gets_no_knockoff() {
        let mut db = reference(b"GGGGSPGGGKGGGGSAGGGK", 4);
        let stats = add_site_knockoffs(&mut db, KO);
        assert_eq!((stats.knockoffs, stats.untested), (0, 1));
    }

    #[test]
    fn a_modified_lysine_is_compared_by_its_merged_span() {
        // Acetyl-K at 9: trypsin does not cut there, so the peptide is
        // 0..20. A K whose merged span is as long qualifies; one followed by
        // the protein end does not.
        let seq = b"GGGGGGGGGKGGGGGGGGGRGGGGGGGGGKGGGGGGGGGRGGGGGGGGGK".to_vec();
        let acetyl = 42.010565;
        let mut db = ProteinDb {
            proteins: vec![Protein {
                accession: "sp|P2|B".into(),
                description: String::new(),
                sequence: seq.clone(),
                annotation: Some(Arc::new(
                    ProteinAnnotation::new(
                        vec![ModSite {
                            position: 9,
                            mass_delta: acetyl,
                            anchor: SiteAnchor::Residue,
                        }],
                        None,
                    )
                    .with_site_claims(vec![SiteClaim {
                        position: 9,
                        mass_delta: acetyl,
                        tier: PriorTier::Curated,
                        knockoff_of: None,
                    }]),
                )),
            }],
        };
        add_site_knockoffs(&mut db, KO);
        let info = db.proteins[0].annotation.as_deref().unwrap();
        let ko: Vec<u32> = info
            .site_claims
            .iter()
            .filter(|c| c.knockoff_of.is_some())
            .map(|c| c.position)
            .collect();
        assert_eq!(ko, vec![29]);
    }
}
