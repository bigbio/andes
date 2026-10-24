//! Modification-site claims and their q-values, after rescoring.
//!
//! A reference protein's annotated sites are claims ("residue r carries
//! modification M"), and each has knockoff sites: the same modification on
//! matched residues no annotation names (see `input::add_site_knockoffs`). A
//! target PSM that places M on an annotated site is evidence for that claim;
//! one that places M on a knockoff site is a false claim. Within a stratum
//! (modification × residue × prior tier), the knockoff PSMs above a score
//! estimate the false claims above it, as for variants
//! ([`crate::knockoff_fdr`]). A PSM's own q-value is not changed: whether a
//! peptide is identified does not depend on which residue carries the
//! modification. Only the site claim gets a q-value.
//!
//! Real unannotated sites exist, so some knockoff sites are really modified,
//! which makes the estimate reject some true claims. In the other direction,
//! a knockoff site lies in another peptide, so the estimate measures whether
//! a residue carries the modification (occupancy), not which one it carries:
//! the right site with a near-isobaric wrong modification passes it. The
//! site posteriors give each accepted site its identity for that.

use std::collections::{BTreeMap, HashMap, HashSet};

use model::{EntryKind, PriorTier, Protein, SiteClaim};
use search::SearchIndex;

use crate::percolator::PercolatorPsm;

/// Tested claim PSMs at the cutoff a modification × residue needs within one
/// tier to get its own stratum (as for variants).
const MIN_STRATUM_TESTED: usize = 50;

/// One accepted site claim.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteRow {
    pub family: String,
    pub accession: String,
    /// 1-based.
    pub position: u32,
    pub residue: char,
    pub mass_delta: f64,
    pub tier: &'static str,
    /// PSMs at q ≤ the PSM cutoff placing the modification there.
    pub psms: usize,
    pub best_psm: String,
    pub site_q: f64,
    /// Share of the posterior on the peptides covering the site that places
    /// a modification there; `None` without row posteriors.
    pub occupancy: Option<f64>,
    /// Given a modification there, the share of this one.
    pub identity: Option<f64>,
    /// The other modifications placed there, as `mass:share`, largest first.
    pub alternatives: String,
    /// Posterior mass placing a modification there, in spectra.
    pub identity_spectra: f64,
}

/// One stratum, as counted for the run log.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteCount {
    pub stratum: String,
    /// Claim PSMs at the PSM cutoff, and of them the tested ones.
    pub claims: usize,
    pub tested: usize,
    /// Knockoff-site PSMs at the PSM cutoff, each counting 1/k.
    pub knockoff_weighted: f64,
    /// Claim PSMs whose site q-value passes the claim cutoff.
    pub accepted: usize,
}

/// What [`site_claims`] found.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SiteOutcome {
    pub rows: Vec<SiteRow>,
    pub counts: Vec<SiteCount>,
    /// Rows that place a `sites=annotated` modification on a knockoff site:
    /// a placement only the knockoffs allow, never to be reported.
    pub knockoff_rows: HashSet<String>,
}

/// The residues of a peptide string (`K.PEPS+79.96633TIDE.K`) and the
/// modifications on them, as (0-based residue index, mass).
pub fn peptide_mods(peptide: &str) -> (Vec<u8>, Vec<(usize, f64)>) {
    let core = match (peptide.find('.'), peptide.rfind('.')) {
        (Some(a), Some(b)) if b > a && peptide[..a].len() <= 2 => &peptide[a + 1..b],
        _ => peptide,
    };
    let bytes = core.as_bytes();
    let (mut residues, mut mods) = (Vec::new(), Vec::new());
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        if c.is_ascii_uppercase() {
            residues.push(c);
            i += 1;
        } else if c == b'+' || c == b'-' {
            let start = i;
            i += 1;
            while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
                i += 1;
            }
            if let Ok(m) = core[start..i].parse::<f64>() {
                mods.push((residues.len().saturating_sub(1), m));
            }
        } else {
            i += 1;
        }
    }
    (residues, mods)
}

/// Start of every occurrence of `pattern` in `seq`, overlapping ones included.
fn occurrences<'a>(seq: &'a [u8], pattern: &'a [u8]) -> impl Iterator<Item = usize> + 'a {
    let mut from = 0;
    std::iter::from_fn(move || {
        if pattern.is_empty() {
            return None;
        }
        let start = from + memchr::memmem::find(seq.get(from..)?, pattern)?;
        from = start + 1;
        Some(start)
    })
}

/// The first mass of each class of claim masses within
/// [`model::protein::SITE_MASS_TOLERANCE`] of each other, ascending.
fn mass_classes(proteins: &[&Protein]) -> Vec<f64> {
    let mut masses: Vec<f64> = proteins
        .iter()
        .filter_map(|p| p.annotation.as_deref())
        .flat_map(|a| a.site_claims.iter().map(|c| c.mass_delta))
        .collect();
    masses.sort_by(f64::total_cmp);
    let mut reps: Vec<f64> = Vec::new();
    for m in masses {
        if reps
            .last()
            .is_none_or(|r| (m - r).abs() > model::protein::SITE_MASS_TOLERANCE)
        {
            reps.push(m);
        }
    }
    reps
}

/// The bits of the class representative `mass` belongs to: one key for every
/// claim mass of one modification, whatever their last bits.
fn canonical_mass_bits(reps: &[f64], mass: f64) -> u64 {
    reps.iter()
        .find(|r| (mass - **r).abs() <= model::protein::SITE_MASS_TOLERANCE)
        .copied()
        .unwrap_or(mass)
        .to_bits()
}

/// A real site: (protein, position, mass bits).
type SiteKey = (usize, u32, u64);
/// A site's best site q-value, its PSM, its PSMs at the cutoff, its tier.
type SiteBest = (f64, String, usize, PriorTier);

/// One row's part in a stratum.
enum Part {
    Claim {
        site: (usize, u32, u64),
        tested: bool,
        /// Knockoff sites of the claimed site (`usize::MAX` untested).
        k: usize,
    },
    Knockoff {
        weight: f64,
        /// Knockoff sites of the real site this one copies.
        k: usize,
    },
}

/// Site q-values for every claim of `targets`' PSMs, and the accepted site
/// claims: a site whose best PSM passes `fdr` and whose site q-value passes
/// `claim_fdr`. `restricted` lists the masses of the `sites=annotated` mods,
/// whose knockoff-site placements are not reported (their rows get q = 1).
/// `posteriors` (each target row's posterior among its spectrum's rows, see
/// [`crate::site_identity`]) give each site its occupancy and identity.
pub fn site_claims(
    targets: &mut HashMap<String, PercolatorPsm>,
    idx: &SearchIndex,
    fdr: f64,
    claim_fdr: f64,
    restricted: &[f64],
    posteriors: Option<&HashMap<String, f64>>,
) -> SiteOutcome {
    let same = |a: f64, b: f64| (a - b).abs() <= model::protein::SITE_MASS_TOLERANCE;
    let proteins: Vec<&Protein> = idx
        .iter_target_proteins()
        .filter(|p| p.entry_kind() == EntryKind::Reference)
        .filter(|p| {
            p.annotation
                .as_deref()
                .is_some_and(|a| !a.site_claims.is_empty())
        })
        .collect();
    if proteins.is_empty() {
        return SiteOutcome::default();
    }
    let by_acc: HashMap<&str, usize> = proteins
        .iter()
        .enumerate()
        .map(|(i, p)| (p.accession.as_str(), i))
        .collect();
    let claims_of = |i: usize| -> &[SiteClaim] {
        &proteins[i]
            .annotation
            .as_deref()
            .expect("filtered on site claims")
            .site_claims
    };
    // Claim masses within tolerance of each other are one modification; every
    // key below uses the class's first mass.
    let reps = mass_classes(&proteins);
    let canon = |m: f64| canonical_mass_bits(&reps, m);
    // Knockoffs per real site, keyed by (protein, position, mass bits).
    let mut per_site: HashMap<(usize, u32, u64), usize> = HashMap::new();
    for i in 0..proteins.len() {
        for c in claims_of(i) {
            if let Some(of) = c.knockoff_of {
                let real = claims_of(i).iter().find(|r| {
                    r.knockoff_of.is_none() && r.position == of && same(r.mass_delta, c.mass_delta)
                });
                if let Some(r) = real {
                    *per_site
                        .entry((i, r.position, canon(r.mass_delta)))
                        .or_insert(0) += 1;
                }
            }
        }
    }

    let rank = |p: &PercolatorPsm| p.score.unwrap_or(-p.q_value);
    // Per (mod, residue, tier): (rank, row id, part).
    type Key = (u64, u8, PriorTier);
    let mut by_key: HashMap<Key, Vec<(f64, String, Part)>> = HashMap::new();
    let mut knockoff_rows = HashSet::new();
    for (id, psm) in targets.iter() {
        let (residues, mods) = peptide_mods(&psm.peptide);
        if mods.is_empty() {
            continue;
        }
        // Each with the modification's offset in the peptide: (stratum, site,
        // tested, offset) and (stratum, weight, restricted, offset).
        type RealPlacement = (Key, (usize, u32, u64), Option<usize>, usize);
        let mut real: Vec<RealPlacement> = Vec::new();
        let mut knock: Vec<(Key, usize, bool, usize)> = Vec::new();
        for acc in psm.proteins.split(';') {
            let Some(&i) = by_acc.get(acc) else {
                continue;
            };
            let seq = &proteins[i].sequence;
            // Every occurrence of the peptide in the protein places its mods.
            for (start, &(at, mass)) in occurrences(seq, &residues)
                .flat_map(|start| mods.iter().map(move |m| (start, m)))
            {
                let pos = (start + at) as u32;
                let Some(c) = claims_of(i)
                    .iter()
                    .find(|c| c.position == pos && same(c.mass_delta, mass))
                else {
                    continue;
                };
                let key = (canon(c.mass_delta), seq[pos as usize], c.tier);
                match c.knockoff_of {
                    Some(of) => {
                        let real_bits = claims_of(i)
                            .iter()
                            .find(|r| {
                                r.knockoff_of.is_none()
                                    && r.position == of
                                    && same(r.mass_delta, c.mass_delta)
                            })
                            .map(|r| canon(r.mass_delta))
                            .unwrap_or(canon(c.mass_delta));
                        let k = per_site.get(&(i, of, real_bits)).copied().unwrap_or(1);
                        let restricted = restricted.iter().any(|&m| same(m, mass));
                        knock.push((key, k, restricted, at));
                    }
                    None => {
                        let site = (i, pos, canon(c.mass_delta));
                        real.push((key, site, per_site.get(&site).copied(), at));
                    }
                }
            }
        }
        // A placement that is a knockoff site in one protein and a real claim
        // in another (a paralog annotated differently) supports that claim.
        knock.retain(|k| !real.iter().any(|r| r.3 == k.3 && r.0 .0 == k.0 .0));
        if let Some(&(key, k, _, _)) = knock.first() {
            by_key.entry(key).or_default().push((
                rank(psm),
                id.clone(),
                Part::Knockoff {
                    weight: 1.0 / k as f64,
                    k,
                },
            ));
            if knock.iter().any(|k| k.2) {
                knockoff_rows.insert(id.clone());
            }
        } else {
            for (key, site, k, _) in real {
                by_key.entry(key).or_default().push((
                    rank(psm),
                    id.clone(),
                    Part::Claim {
                        site,
                        tested: k.is_some(),
                        k: k.unwrap_or(usize::MAX),
                    },
                ));
            }
        }
    }

    // Strata: mod × residue within a tier when it has enough tested claims at
    // the cutoff, else across its small tiers; mods and residues never pool.
    let tiers: std::collections::BTreeSet<PriorTier> = by_key.keys().map(|k| k.2).collect();
    let mut strata: HashMap<String, Vec<(f64, String, Part)>> = HashMap::new();
    for ((mass_bits, residue, tier), parts) in by_key {
        let tested = parts
            .iter()
            .filter(|(_, id, p)| {
                matches!(p, Part::Claim { tested: true, .. }) && targets[id].q_value <= fdr
            })
            .count();
        let base = format!("{:+.4}@{}", f64::from_bits(mass_bits), residue as char);
        let name = if tiers.len() < 2 {
            base
        } else if tested >= MIN_STRATUM_TESTED {
            format!("{base}/{}", tier.name())
        } else {
            format!("{base}/pooled")
        };
        strata.entry(name).or_default().extend(parts);
    }

    // Best site q-value per real site, with its PSM, and PSMs per site.
    let mut best: BTreeMap<SiteKey, SiteBest> = BTreeMap::new();
    let mut counts = Vec::new();
    for (name, mut parts) in strata {
        let mut count = SiteCount {
            stratum: name,
            claims: 0,
            tested: 0,
            knockoff_weighted: 0.0,
            accepted: 0,
        };
        // A PSM counts once per stratum on the claim side, however many
        // proteins hold its claim (a peptide shared across genes), as a
        // knockoff PSM does on the knockoff side.
        let (mut seen_claim, mut seen_tested) = (HashSet::new(), HashSet::new());
        for (_, id, part) in &parts {
            if targets[id].q_value <= fdr {
                match part {
                    Part::Knockoff { weight, .. } => count.knockoff_weighted += weight,
                    Part::Claim { tested, .. } => {
                        count.claims += seen_claim.insert(id.as_str()) as usize;
                        count.tested += (*tested && seen_tested.insert(id.as_str())) as usize;
                    }
                }
            }
        }
        parts.sort_by(|a, b| {
            b.0.total_cmp(&a.0).then(
                matches!(b.2, Part::Knockoff { .. }).cmp(&matches!(a.2, Part::Knockoff { .. })),
            )
        });
        let offset =
            crate::knockoff_fdr::finite_sample_term(parts.iter().map(|(_, _, part)| match *part {
                Part::Claim { k, .. } | Part::Knockoff { k, .. } => k,
            }));
        let mut fdp = Vec::with_capacity(parts.len());
        let (mut w, mut t) = (0.0f64, 0usize);
        let mut counted = HashSet::new();
        for (_, id, part) in &parts {
            // Only PSMs that pass the PSM cutoff count, on both sides; a PSM
            // counts once however many proteins hold its claim.
            if targets[id].q_value <= fdr {
                match part {
                    Part::Knockoff { weight, .. } => w += weight,
                    Part::Claim { tested, .. } => {
                        t += (*tested && counted.insert(id.as_str())) as usize
                    }
                }
            }
            fdp.push(((offset + w) / t.max(1) as f64).min(1.0));
        }
        for i in (0..fdp.len().saturating_sub(1)).rev() {
            fdp[i] = fdp[i].min(fdp[i + 1]);
        }
        for ((_, id, part), q) in parts.iter().zip(fdp) {
            let Part::Claim { site, tested, .. } = part else {
                continue;
            };
            let psm_q = targets[id].q_value;
            if psm_q > fdr {
                continue;
            }
            let site_q = if *tested { q } else { 1.0 };
            count.accepted += (site_q <= claim_fdr) as usize;
            let tier = claims_of(site.0)
                .iter()
                .find(|c| {
                    c.knockoff_of.is_none()
                        && c.position == site.1
                        && canon(c.mass_delta) == site.2
                })
                .map_or(PriorTier::Catalogue, |c| c.tier);
            let e = best
                .entry(*site)
                .or_insert((f64::INFINITY, String::new(), 0, tier));
            e.2 += 1;
            if site_q < e.0 || (site_q == e.0 && *id < e.1) {
                e.0 = site_q;
                e.1 = id.clone();
            }
        }
        counts.push(count);
    }
    counts.sort_by(|a, b| a.stratum.cmp(&b.stratum));

    for id in &knockoff_rows {
        if let Some(p) = targets.get_mut(id) {
            p.q_value = 1.0;
            p.pep = 1.0;
        }
    }
    let identity = posteriors.map(|post| site_posteriors(targets, &proteins, &by_acc, post));
    let rows = best
        .into_iter()
        .filter(|(_, (q, ..))| *q <= claim_fdr)
        .map(|((i, pos, bits), (q, psm, n, tier))| {
            let p = proteins[i];
            let (occupancy, share, alternatives, spectra) = identity
                .as_ref()
                .and_then(|m| m.get(&(i, pos)))
                .map_or((None, None, String::new(), 0.0), |s| {
                    s.describe(f64::from_bits(bits), &same)
                });
            SiteRow {
                family: p.origin().map_or_else(
                    || primary(&p.accession).to_string(),
                    |o| o.family.to_string(),
                ),
                accession: p.accession.clone(),
                position: pos + 1,
                residue: p.sequence[pos as usize] as char,
                mass_delta: f64::from_bits(bits),
                tier: tier.name(),
                psms: n,
                best_psm: psm,
                site_q: q,
                occupancy,
                identity: share,
                alternatives,
                identity_spectra: spectra,
            }
        })
        .collect();
    SiteOutcome {
        rows,
        counts,
        knockoff_rows,
    }
}

/// Posterior mass at one site: per modification mass placed there, and on
/// every row whose peptide covers it.
#[derive(Debug, Default)]
struct SitePosterior {
    mods: Vec<(f64, f64)>,
    covering: f64,
}

impl SitePosterior {
    /// (occupancy, identity of `mass`, the other mods as `mass:share`,
    /// posterior mass placing a mod).
    fn describe(
        &self,
        mass: f64,
        same: &dyn Fn(f64, f64) -> bool,
    ) -> (Option<f64>, Option<f64>, String, f64) {
        let placed: f64 = self.mods.iter().map(|m| m.1).sum();
        if placed <= 0.0 {
            return (None, None, String::new(), 0.0);
        }
        let mut merged: Vec<(f64, f64)> = Vec::new();
        for &(m, w) in &self.mods {
            match merged.iter_mut().find(|x| same(x.0, m)) {
                Some(x) => x.1 += w,
                None => merged.push((m, w)),
            }
        }
        merged.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.total_cmp(&b.0)));
        let own = merged.iter().find(|x| same(x.0, mass)).map_or(0.0, |x| x.1);
        let alternatives = merged
            .iter()
            .filter(|x| !same(x.0, mass) && x.1 / placed >= 0.01)
            .map(|x| format!("{:+.4}:{:.2}", x.0, x.1 / placed))
            .collect::<Vec<_>>()
            .join(";");
        (
            Some((placed / self.covering.max(placed)).min(1.0)),
            Some(own / placed),
            alternatives,
            placed,
        )
    }
}

/// Posterior mass at every real claim site of `proteins`.
fn site_posteriors(
    targets: &HashMap<String, PercolatorPsm>,
    proteins: &[&Protein],
    by_acc: &HashMap<&str, usize>,
    posteriors: &HashMap<String, f64>,
) -> HashMap<(usize, u32), SitePosterior> {
    let sites: Vec<Vec<u32>> = proteins
        .iter()
        .map(|p| {
            let mut v: Vec<u32> = p
                .annotation
                .as_deref()
                .map(|a| {
                    a.site_claims
                        .iter()
                        .filter(|c| c.knockoff_of.is_none())
                        .map(|c| c.position)
                        .collect()
                })
                .unwrap_or_default();
            v.sort_unstable();
            v.dedup();
            v
        })
        .collect();
    let mut out: HashMap<(usize, u32), SitePosterior> = HashMap::new();
    for (id, psm) in targets {
        let Some(&w) = posteriors.get(id) else {
            continue;
        };
        let (residues, mods) = peptide_mods(&psm.peptide);
        for acc in psm.proteins.split(';') {
            let Some(&i) = by_acc.get(acc) else {
                continue;
            };
            for start in occurrences(&proteins[i].sequence, &residues) {
                let (lo, hi) = (start as u32, (start + residues.len()) as u32);
                let from = sites[i].partition_point(|&p| p < lo);
                for &pos in sites[i][from..].iter().take_while(|&&p| p < hi) {
                    let s = out.entry((i, pos)).or_default();
                    s.covering += w;
                    for &(at, mass) in &mods {
                        if start as u32 + at as u32 == pos {
                            s.mods.push((mass, w));
                        }
                    }
                }
            }
        }
    }
    out
}

fn primary(accession: &str) -> &str {
    let parts: Vec<&str> = accession.split('|').collect();
    match parts.as_slice() {
        [_, acc, _] => acc,
        [one] => one.split_once(':').map_or(one, |(_, acc)| acc),
        _ => accession,
    }
}

/// `field` with the characters that would break the table's rows or columns
/// (tab, line breaks) as spaces.
fn tsv_field(field: &str) -> std::borrow::Cow<'_, str> {
    if field.contains(['\t', '\n', '\r']) {
        field.replace(['\t', '\n', '\r'], " ").into()
    } else {
        field.into()
    }
}

/// Write the site table: one row per accepted site claim.
pub fn write_sites_tsv(path: &std::path::Path, rows: &[SiteRow]) -> std::io::Result<()> {
    use std::io::Write;
    crate::atomic::write_atomic(path, |path| {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "family\taccession\tposition\tresidue\tmodification_mass\tprior_tier\tpsms\tbest_psm\tsite_q\tconfidence\toccupancy\tidentity\tidentity_alternatives\tidentity_spectra"
        )?;
        let opt = |v: Option<f64>| v.map_or_else(|| "NA".to_string(), |v| format!("{v:.4}"));
        for r in rows {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{:.2}",
                tsv_field(&r.family),
                tsv_field(&r.accession),
                r.position,
                r.residue,
                r.mass_delta,
                r.tier,
                r.psms,
                tsv_field(&r.best_psm),
                r.site_q,
                1.0 - r.site_q,
                opt(r.occupancy),
                opt(r.identity),
                tsv_field(&r.alternatives),
                r.identity_spectra
            )?;
        }
        w.flush()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::{ModSite, ProteinAnnotation, ProteinDb, SiteAnchor};
    use std::sync::Arc;

    const PHOSPHO: f64 = 79.966331;

    #[test]
    fn peptide_strings_give_residues_and_mods() {
        let (r, m) = peptide_mods("K.PEPS+79.96633TIDE.K");
        assert_eq!(r, b"PEPSTIDE");
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].0, 3);
        assert!((m[0].1 - 79.96633).abs() < 1e-9);
        let (r, m) = peptide_mods("R.A+229.16293VFVK+229.16293.D");
        assert_eq!(r, b"AVFVK");
        assert_eq!(m.iter().map(|x| x.0).collect::<Vec<_>>(), [0, 4]);
        assert!(peptide_mods("-.PEPTIDEK.-").1.is_empty());
    }

    /// Spans `<X>GGGSPGGGK` with a different first residue each; the S at 4
    /// is an annotated phospho site, the S at 14 and 24 its knockoffs.
    fn index() -> SearchIndex {
        SearchIndex::from_target_db(
            &ProteinDb {
                proteins: vec![protein(
                    "sp|P1|A",
                    &[(4, None), (14, Some(4)), (24, Some(4))],
                )],
            },
            "XXX",
        )
    }

    /// A protein of the [`index`] sequence with phospho claims at the given
    /// positions, each a real claim or a knockoff of another.
    fn protein(accession: &str, claims: &[(u32, Option<u32>)]) -> Protein {
        let seq: Vec<u8> = b"ACDEF"
            .iter()
            .flat_map(|&x| {
                let mut s = vec![x];
                s.extend(b"GGGSPGGGK");
                s
            })
            .collect();
        let site = |p: u32| ModSite {
            position: p,
            mass_delta: PHOSPHO,
            anchor: SiteAnchor::Residue,
        };
        let claim = |&(p, of): &(u32, Option<u32>)| SiteClaim {
            position: p,
            mass_delta: PHOSPHO,
            tier: PriorTier::Curated,
            knockoff_of: of,
        };
        Protein {
            accession: accession.into(),
            description: String::new(),
            sequence: seq,
            annotation: Some(Arc::new(
                ProteinAnnotation::new(claims.iter().map(|c| site(c.0)).collect(), None)
                    .with_site_claims(claims.iter().map(claim).collect()),
            )),
        }
    }

    fn psm(id: &str, score: f64, peptide: &str) -> (String, PercolatorPsm) {
        (
            id.to_string(),
            PercolatorPsm {
                psm_id: id.to_string(),
                q_value: 0.001,
                pep: 0.0,
                score: Some(score),
                peptide: peptide.into(),
                proteins: "sp|P1|A".into(),
            },
        )
    }

    fn rows() -> HashMap<String, PercolatorPsm> {
        (0..60)
            .map(|i| {
                psm(
                    &format!("c{i}"),
                    100.0 - i as f64,
                    "K.AGGGS+79.96633PGGGK.C",
                )
            })
            .chain([
                psm("k1", 99.5, "K.CGGGS+79.96633PGGGK.D"),
                psm("k2", 50.0, "K.DGGGS+79.96633PGGGK.E"),
                psm("plain", 98.0, "K.AGGGSPGGGK.C"),
            ])
            .collect()
    }

    #[test]
    fn a_site_claim_gets_its_knockoff_q_value() {
        let idx = index();
        let mut t = rows();
        let out = site_claims(&mut t, &idx, 0.01, 0.05, &[PHOSPHO], None);
        // Two knockoff PSMs of a site with two knockoffs weigh 1/2 each. The
        // threshold just above the second one (50 claims, one knockoff) gives
        // (1/2 + 1/2) / 50 = 0.02, the least FDP any claim reaches.
        assert_eq!(out.rows.len(), 1);
        let r = &out.rows[0];
        assert_eq!(
            (r.position, r.residue, r.psms, r.tier),
            (5, 'S', 60, "curated")
        );
        assert!((r.site_q - 0.02).abs() < 1e-9, "{}", r.site_q);
        assert_eq!(
            out.knockoff_rows,
            HashSet::from(["k1".to_string(), "k2".to_string()])
        );
        assert_eq!(
            t["k1"].q_value, 1.0,
            "a restricted mod's knockoff placement is not reported"
        );
        assert_eq!(t["c0"].q_value, 0.001, "a PSM's own q-value is unchanged");
        // A stricter claim cutoff accepts no site.
        let mut t = rows();
        assert!(site_claims(&mut t, &idx, 0.01, 0.015, &[PHOSPHO], None)
            .rows
            .is_empty());
    }

    /// The finite-sample term is the stratum's own: 1/k for the fewest
    /// knockoffs of a site it holds (here P1's two), not 1/k of the site with
    /// the most knockoffs in the database (P2's four).
    #[test]
    fn the_finite_sample_term_uses_the_strata_fewest_knockoffs() {
        let idx = SearchIndex::from_target_db(
            &ProteinDb {
                proteins: vec![
                    protein("sp|P1|A", &[(4, None), (14, Some(4)), (24, Some(4))]),
                    protein(
                        "sp|P2|A",
                        &[
                            (4, None),
                            (14, Some(4)),
                            (24, Some(4)),
                            (34, Some(4)),
                            (44, Some(4)),
                        ],
                    ),
                ],
            },
            "XXX",
        );
        let mut t: HashMap<_, _> = (0..100)
            .map(|i| {
                psm(
                    &format!("c{i}"),
                    100.0 - i as f64,
                    "K.AGGGS+79.96633PGGGK.C",
                )
            })
            .collect();
        let (id, mut p2) = psm("p2", 1.0, "K.AGGGS+79.96633PGGGK.C");
        p2.proteins = "sp|P2|A".into();
        t.insert(id, p2);
        let out = site_claims(&mut t, &idx, 0.01, 0.05, &[PHOSPHO], None);
        let p1 = out.rows.iter().find(|r| r.accession == "sp|P1|A").unwrap();
        assert!((p1.site_q - 0.5 / 101.0).abs() < 1e-12, "{}", p1.site_q);
    }

    #[test]
    fn posteriors_give_a_site_its_occupancy_and_identity() {
        let idx = index();
        let mut t = rows();
        // Three rows place another mass on the site, and the unmodified row
        // covers it: identity 60/63, occupancy 63/64.
        for i in 0..3 {
            let (id, p) = psm(&format!("x{i}"), 1.0, "K.AGGGS+42.01056PGGGK.C");
            t.insert(id, p);
        }
        let post: HashMap<String, f64> = t
            .keys()
            .filter(|k| !k.starts_with('k'))
            .map(|k| (k.clone(), 1.0))
            .collect();
        let out = site_claims(&mut t, &idx, 0.01, 0.05, &[PHOSPHO], Some(&post));
        let r = &out.rows[0];
        assert!((r.identity.unwrap() - 60.0 / 63.0).abs() < 1e-9);
        assert!((r.occupancy.unwrap() - 63.0 / 64.0).abs() < 1e-9);
        assert_eq!(r.alternatives, "+42.0106:0.05");
        assert!((r.identity_spectra - 63.0).abs() < 1e-9);
        let out = site_claims(&mut rows(), &idx, 0.01, 0.05, &[PHOSPHO], None);
        assert_eq!((out.rows[0].occupancy, out.rows[0].identity), (None, None));
    }

    #[test]
    fn an_unrestricted_mod_keeps_its_knockoff_site_psms() {
        let idx = index();
        let mut t = rows();
        let out = site_claims(&mut t, &idx, 0.01, 0.05, &[], None);
        assert!(out.knockoff_rows.is_empty());
        assert_eq!(t["k1"].q_value, 0.001);
        assert_eq!(out.rows.len(), 1, "the site claim is still controlled");
    }

    #[test]
    fn a_psm_shared_across_genes_counts_once() {
        // P2 is a second gene with the same sequence and claims: each claim
        // PSM names both, and still counts once against the knockoffs.
        let idx = SearchIndex::from_target_db(
            &ProteinDb {
                proteins: vec![
                    protein("sp|P1|A", &[(4, None), (14, Some(4)), (24, Some(4))]),
                    protein("sp|P2|B", &[(4, None), (14, Some(4)), (24, Some(4))]),
                ],
            },
            "XXX",
        );
        let mut t = rows();
        for p in t.values_mut() {
            p.proteins = "sp|P1|A;sp|P2|B".into();
        }
        let out = site_claims(&mut t, &idx, 0.01, 0.05, &[PHOSPHO], None);
        assert_eq!(out.rows.len(), 2);
        for r in &out.rows {
            assert!((r.site_q - 0.02).abs() < 1e-9, "{}", r.site_q);
        }
        assert_eq!(out.counts[0].claims, 60);
    }

    #[test]
    fn a_paralogs_real_claim_wins_over_a_knockoff_site() {
        // Position 14 is a knockoff site of P1 but a real claim of P2: the
        // PSM placing phospho there supports P2's claim and stays reported.
        let idx = SearchIndex::from_target_db(
            &ProteinDb {
                proteins: vec![
                    protein("sp|P1|A", &[(4, None), (14, Some(4)), (24, Some(4))]),
                    protein("sp|P2|B", &[(14, None)]),
                ],
            },
            "XXX",
        );
        let mut t = rows();
        t.get_mut("k1").unwrap().proteins = "sp|P1|A;sp|P2|B".into();
        let out = site_claims(&mut t, &idx, 0.01, 0.05, &[PHOSPHO], None);
        assert!(!out.knockoff_rows.contains("k1"));
        assert_eq!(t["k1"].q_value, 0.001);
        assert!(out.knockoff_rows.contains("k2"));
    }

    /// A protein of `chunks`, each `<X>GGGSPGGGK`, with the given claims.
    fn protein_of(accession: &str, chunks: &[u8], claims: &[(u32, Option<u32>)]) -> Protein {
        let mut p = protein(accession, claims);
        p.sequence = chunks
            .iter()
            .flat_map(|&x| {
                let mut s = vec![x];
                s.extend(b"GGGSPGGGK");
                s
            })
            .collect();
        p
    }

    #[test]
    fn a_peptide_repeated_in_one_protein_reports_every_site() {
        // The chunk `AGGGSPGGGK` occurs twice; the S of each is a claimed site
        // (positions 4 and 14) with two knockoff sites of its own.
        let idx = SearchIndex::from_target_db(
            &ProteinDb {
                proteins: vec![protein_of(
                    "sp|P1|A",
                    b"AACDEF",
                    &[
                        (4, None),
                        (14, None),
                        (24, Some(4)),
                        (34, Some(4)),
                        (44, Some(14)),
                        (54, Some(14)),
                    ],
                )],
            },
            "XXX",
        );
        let mut t = rows();
        for (id, score, pep) in [
            ("k3", 1.0, "K.EGGGS+79.96633PGGGK.F"),
            ("k4", 0.5, "K.FGGGS+79.96633PGGGK.G"),
        ] {
            let (id, p) = psm(id, score, pep);
            t.insert(id, p);
        }
        let out = site_claims(&mut t, &idx, 0.01, 0.05, &[PHOSPHO], None);
        let sites: Vec<(u32, usize)> = out.rows.iter().map(|r| (r.position, r.psms)).collect();
        assert_eq!(sites, [(5, 60), (15, 60)]);
    }

    #[test]
    fn occurrences_include_overlapping_ones() {
        assert_eq!(occurrences(b"AAAA", b"AA").collect::<Vec<_>>(), [0, 1, 2]);
        assert_eq!(occurrences(b"PEPXPEP", b"PEP").collect::<Vec<_>>(), [0, 4]);
        assert!(occurrences(b"PEP", b"K").next().is_none());
        assert!(occurrences(b"PEP", b"").next().is_none());
    }

    #[test]
    fn masses_differing_in_the_last_bits_share_one_key() {
        let near = f64::from_bits(PHOSPHO.to_bits() + 3);
        let acetyl = 42.010565;
        let p = |m: f64| {
            let mut p = protein("sp|P|A", &[(4, None)]);
            let a = p.annotation.take().unwrap();
            let mut claims = a.site_claims.clone();
            claims[0].mass_delta = m;
            p.annotation = Some(Arc::new(
                ProteinAnnotation::new(a.mod_sites.clone(), None).with_site_claims(claims),
            ));
            p
        };
        let (a, b, c) = (p(PHOSPHO), p(near), p(acetyl));
        let reps = mass_classes(&[&a, &b, &c]);
        assert_eq!(reps, [acetyl, PHOSPHO]);
        assert_eq!(
            canonical_mass_bits(&reps, near),
            canonical_mass_bits(&reps, PHOSPHO)
        );
        assert_ne!(
            canonical_mass_bits(&reps, acetyl),
            canonical_mass_bits(&reps, PHOSPHO)
        );
    }

    #[test]
    fn knockoff_claims_with_nearly_equal_masses_leave_the_result_unchanged() {
        let make = |knock_mass: f64| {
            let mut p = protein("sp|P1|A", &[(4, None), (14, Some(4)), (24, Some(4))]);
            let a = p.annotation.take().unwrap();
            let mut claims = a.site_claims.clone();
            for c in claims.iter_mut().filter(|c| c.knockoff_of.is_some()) {
                c.mass_delta = knock_mass;
            }
            p.annotation = Some(Arc::new(
                ProteinAnnotation::new(a.mod_sites.clone(), None).with_site_claims(claims),
            ));
            SearchIndex::from_target_db(&ProteinDb { proteins: vec![p] }, "XXX")
        };
        let run = |idx: &SearchIndex| {
            let mut t = rows();
            site_claims(&mut t, idx, 0.01, 0.05, &[PHOSPHO], None)
        };
        let exact = run(&make(PHOSPHO));
        let near = run(&make(f64::from_bits(PHOSPHO.to_bits() + 3)));
        assert_eq!(exact.rows.len(), 1);
        assert_eq!(exact.rows, near.rows);
        assert_eq!(exact.counts, near.counts);
    }

    #[test]
    fn the_site_table_escapes_tabs_and_line_breaks() {
        let row = SiteRow {
            family: "fam\tily".into(),
            accession: "sp|P1|A\nB".into(),
            position: 5,
            residue: 'S',
            mass_delta: PHOSPHO,
            tier: "curated",
            psms: 3,
            best_psm: "run/t\r1_1_1".into(),
            site_q: 0.01,
            occupancy: None,
            identity: None,
            alternatives: String::new(),
            identity_spectra: 0.0,
        };
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sites.tsv");
        write_sites_tsv(&path, &[row]).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        for l in &lines {
            assert_eq!(l.split('\t').count(), 14, "{l}");
        }
        assert!(lines[1].starts_with("fam ily\tsp|P1|A B\t5\tS\t"));
    }
}
