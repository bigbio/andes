//! The proteoform table: per protein family, the canonical entry and each
//! isoform, variant and chain entry with accepted evidence, its confidence and
//! its quantity.
//!
//! A peptide is *specific* to a member when, within its family, no other
//! member explains it, and, for an isoform, variant or chain member, no
//! reference entry of any family holds it (a paralog, or a reference sharing it
//! up to I/L); otherwise it is *shared*. An isoform, variant or chain
//! entry is listed only with specific evidence. A member's confidence is
//! `1 - q` of its best specific PSM (the class-aware q-value of
//! [`crate::class_fdr`]), or of its best shared PSM when it has no specific one
//! (then `evidence` is `shared`). Its quantity is the summed abundance of its
//! specific peptides in each run or channel.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use model::{EntryKind, Protein};
use search::SearchIndex;

use crate::class_fdr::Cutoffs;
use crate::percolator::PercolatorPsm;

/// Quantities of identified peptides, from the quantification step.
#[derive(Debug, Clone, Default)]
pub struct QuantSummary {
    /// Run names (label-free) or channel names (isobaric), in column order.
    pub columns: Vec<String>,
    /// (bare peptide sequence, its proteins, column, abundance).
    pub values: Vec<(String, Vec<String>, usize, f64)>,
}

/// One family member of the proteoform table.
#[derive(Debug, Clone, PartialEq)]
pub struct ProteoformRow {
    pub family: String,
    pub member: String,
    pub kind: EntryKind,
    pub feature: String,
    pub specific_peptides: usize,
    pub shared_peptides: usize,
    pub psms: usize,
    /// `specific` or `shared`: which PSMs the confidence comes from.
    pub evidence: &'static str,
    pub best_psm: String,
    pub q_value: f64,
    /// For a variant: its substitution's mass coincidence with a common
    /// modification, and whether the canonical counterpart peptide was
    /// identified. Empty otherwise.
    pub alternative: String,
    /// The database the member's claim comes from (empty for a FASTA entry).
    pub source: String,
    /// How likely claims like the member's are to be real (empty for a
    /// FASTA entry).
    pub tier: &'static str,
    /// Summed abundance of the member's specific peptides, per quant column.
    pub quantity: Vec<f64>,
}

impl ProteoformRow {
    pub fn confidence(&self) -> f64 {
        1.0 - self.q_value
    }
}

/// The residues of a Percolator peptide string (`K.PEPT+79.97IDEK.A`).
pub fn bare_sequence(peptide: &str) -> String {
    let core = match (peptide.find('.'), peptide.rfind('.')) {
        (Some(a), Some(b)) if b > a => &peptide[a + 1..b],
        _ => peptide,
    };
    core.chars().filter(|c| c.is_ascii_uppercase()).collect()
}

/// A substitution mass that equals a common modification on the canonical
/// residue, so the variant could be the canonical peptide modified.
pub fn mass_coincidence(delta_da: f64) -> Option<&'static str> {
    [
        ("deamidation", 0.984016),
        ("methylation", 14.01565),
        ("oxidation", 15.994915),
    ]
    .into_iter()
    .find(|(_, m)| (delta_da.abs() - m).abs() < 0.02)
    .map(|(name, _)| name)
}

/// Build the table from the rescored targets passing their `cutoffs`.
pub fn proteoform_table(
    targets: &HashMap<String, PercolatorPsm>,
    idx: &SearchIndex,
    cutoffs: &Cutoffs,
    quant: Option<&QuantSummary>,
) -> Vec<ProteoformRow> {
    let proteins: HashMap<&str, &Protein> = idx
        .iter_target_proteins()
        .map(|p| (p.accession.as_str(), p))
        .collect();
    let family_of = |p: &Protein| -> String {
        p.origin()
            .map_or_else(|| p.accession.clone(), |o| o.family.to_string())
    };
    #[derive(Default)]
    struct Acc {
        specific: BTreeSet<String>,
        shared: BTreeSet<String>,
        psms: usize,
        best_specific: Option<(f64, String, String)>,
        best_shared: Option<(f64, String)>,
    }
    let better = |q: f64, id: &str, cur: &Option<(f64, String)>| -> bool {
        cur.as_ref()
            .is_none_or(|(cq, cid)| q < *cq || (q == *cq && id < cid.as_str()))
    };
    let mut members: BTreeMap<(String, String), Acc> = BTreeMap::new();
    for (id, t) in targets.iter().filter(|(id, t)| t.q_value <= cutoffs.of(id)) {
        let seq = bare_sequence(&t.peptide);
        let mut by_family: BTreeMap<String, BTreeSet<&str>> = BTreeMap::new();
        let mut has_reference = false;
        for acc in t.proteins.split(';').filter(|a| !a.is_empty()) {
            if let Some(p) = proteins.get(acc) {
                by_family.entry(family_of(p)).or_default().insert(acc);
                has_reference |= p.entry_kind() == EntryKind::Reference;
            }
        }
        for (family, accs) in by_family {
            for acc in accs.iter().copied() {
                // A reference of another family explains the peptide as well,
                // so it is no evidence for an alternative member.
                let specific = accs.len() == 1
                    && (proteins[acc].entry_kind() == EntryKind::Reference || !has_reference);
                let a = members
                    .entry((family.clone(), acc.to_string()))
                    .or_default();
                if specific {
                    a.specific.insert(seq.clone());
                    a.psms += 1;
                    let cur = a.best_specific.as_ref().map(|(q, i, _)| (*q, i.clone()));
                    if better(t.q_value, id, &cur) {
                        a.best_specific = Some((t.q_value, id.clone(), seq.clone()));
                    }
                } else {
                    a.shared.insert(seq.clone());
                    if better(t.q_value, id, &a.best_shared) {
                        a.best_shared = Some((t.q_value, id.clone()));
                    }
                }
            }
        }
    }
    // Peptides each family's canonical and isoform entries were identified by.
    let mut non_variant: HashMap<String, BTreeSet<String>> = HashMap::new();
    for ((family, acc), a) in &members {
        if proteins[acc.as_str()].entry_kind() != EntryKind::Variant {
            let set = non_variant.entry(family.clone()).or_default();
            set.extend(a.specific.iter().cloned());
            set.extend(a.shared.iter().cloned());
        }
    }
    let mut rows: Vec<ProteoformRow> = members
        .iter()
        // A non-canonical member is listed only with a peptide no other member
        // of its family explains; shared peptides show the family, not it.
        .filter(|((_, acc), a)| {
            proteins[acc.as_str()].entry_kind() == EntryKind::Reference || !a.specific.is_empty()
        })
        .map(|((family, acc), a)| {
            let p = proteins[acc.as_str()];
            let (evidence, q_value, best_psm) = match (&a.best_specific, &a.best_shared) {
                (Some((q, id, _)), _) => ("specific", *q, id.clone()),
                (None, Some((q, id))) => ("shared", *q, id.clone()),
                (None, None) => ("shared", 1.0, String::new()),
            };
            let alternative = match (p.entry_kind(), &a.best_specific) {
                (EntryKind::Variant, Some((_, _, seq))) => {
                    variant_alternative(p, seq, non_variant.get(family))
                }
                _ => String::new(),
            };
            let quantity = quant.map_or_else(Vec::new, |qs| {
                let mut v = vec![0.0; qs.columns.len()];
                for (seq, prots, col, value) in &qs.values {
                    if a.specific.contains(seq) && prots.iter().any(|x| x == acc) {
                        if let Some(slot) = v.get_mut(*col) {
                            *slot += value;
                        }
                    }
                }
                v
            });
            ProteoformRow {
                family: family.clone(),
                member: acc.clone(),
                kind: p.entry_kind(),
                feature: p
                    .origin()
                    .map_or_else(String::new, |o| o.feature.to_string()),
                source: p
                    .origin()
                    .and_then(|o| o.source.as_deref())
                    .unwrap_or_default()
                    .to_string(),
                tier: p.origin().map_or("", |o| o.tier.name()),
                specific_peptides: a.specific.len(),
                shared_peptides: a.shared.len(),
                psms: a.psms,
                evidence,
                best_psm,
                q_value,
                alternative,
                quantity,
            }
        })
        .collect();
    let order = |k: EntryKind| match k {
        EntryKind::Reference => 0,
        EntryKind::Isoform => 1,
        EntryKind::Variant => 2,
        EntryKind::Chain => 3,
    };
    rows.sort_by(|a, b| {
        a.family
            .cmp(&b.family)
            .then(order(a.kind).cmp(&order(b.kind)))
            .then(a.member.cmp(&b.member))
    });
    rows
}

/// Why a variant's best specific peptide may be the canonical peptide after
/// all: its substitution weighs the same as a common modification, and/or the
/// canonical counterpart peptide (the substitutions reverted) was identified.
fn variant_alternative(p: &Protein, seq: &str, identified: Option<&BTreeSet<String>>) -> String {
    let mut parts: Vec<String> = Vec::new();
    if let Some(name) = p
        .origin()
        .and_then(|o| mass_coincidence(o.variant_delta_da))
    {
        parts.push(format!("{name}-like mass"));
    }
    let haystack = String::from_utf8_lossy(&p.sequence);
    if let Some(start) = haystack.find(seq) {
        let mut counterpart: Vec<u8> = seq.as_bytes().to_vec();
        let mut changed = false;
        for s in p.substitutions() {
            let at = s.position as usize;
            if at >= start && at < start + seq.len() {
                counterpart[at - start] = s.original;
                changed = true;
            }
        }
        let counterpart = String::from_utf8_lossy(&counterpart).into_owned();
        if changed && identified.is_some_and(|set| set.contains(&counterpart)) {
            parts.push(format!("canonical counterpart {counterpart} identified"));
        }
    }
    parts.join("; ")
}

/// Write the table as TSV.
pub fn write_proteoform_tsv(
    path: &std::path::Path,
    rows: &[ProteoformRow],
    quant_columns: &[String],
) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    write!(
        w,
        "family\tmember\tkind\tfeature\tconfidence\tq_value\tevidence\tspecific_peptides\tshared_peptides\tpsms\tbest_psm\talternative_explanation\tsource\tprior_tier"
    )?;
    for c in quant_columns {
        write!(w, "\tquantity_{c}")?;
    }
    writeln!(w)?;
    for r in rows {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.family,
            r.member,
            r.kind.name(),
            r.feature,
            r.confidence(),
            r.q_value,
            r.evidence,
            r.specific_peptides,
            r.shared_peptides,
            r.psms,
            r.best_psm,
            r.alternative,
            r.source,
            r.tier
        )?;
        for v in &r.quantity {
            write!(w, "\t{v}")?;
        }
        writeln!(w)?;
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::{EntryOrigin, ProteinAnnotation, ProteinDb, Substitution};
    use std::sync::Arc;

    fn entry(
        acc: &str,
        seq: &[u8],
        kind: EntryKind,
        feature: &str,
        delta: f64,
        subs: Vec<Substitution>,
    ) -> Protein {
        Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: seq.to_vec(),
            annotation: Some(Arc::new(
                ProteinAnnotation::new(Vec::new(), None)
                    .with_origin(Some(EntryOrigin {
                        kind,
                        family: Arc::from("P1"),
                        feature: Arc::from(feature),
                        variant_delta_da: delta,
                        knockoff_of: None,
                        source: None,
                        tier: model::PriorTier::Catalogue,
                    }))
                    .with_substitutions(subs),
            )),
        }
    }

    /// Canonical P1, isoform P1-2 (shares AAAK, has its own CCCK), and
    /// variant N10D of P1 (VLDGSTR from canonical VLNGSTR).
    fn index() -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![
                entry(
                    "sp|P1|A",
                    b"AAAKEEEKVLNGSTR",
                    EntryKind::Reference,
                    "",
                    0.0,
                    vec![],
                ),
                entry(
                    "sp|P1-2|A",
                    b"AAAKCCCK",
                    EntryKind::Isoform,
                    "P1-2",
                    0.0,
                    vec![],
                ),
                entry(
                    "sp|P1_VAR_1_N10D|A",
                    b"AAAKEEEKVLDGSTR",
                    EntryKind::Variant,
                    "VAR_1_N10D",
                    0.984016,
                    vec![Substitution {
                        position: 10,
                        original: b'N',
                    }],
                ),
            ],
        };
        SearchIndex::from_target_db(&target, "XXX")
    }

    fn psm(id: &str, q: f64, pep: &str, proteins: &str) -> (String, PercolatorPsm) {
        (
            id.into(),
            PercolatorPsm {
                psm_id: id.into(),
                q_value: q,
                pep: 0.0,
                score: Some(1.0),
                peptide: format!("K.{pep}.A"),
                proteins: proteins.into(),
            },
        )
    }

    fn targets() -> HashMap<String, PercolatorPsm> {
        [
            psm("s1", 0.001, "AAAK", "sp|P1|A;sp|P1-2|A"),
            psm("s2", 0.002, "EEEK", "sp|P1|A"),
            psm("s3", 0.004, "CCCK", "sp|P1-2|A"),
            psm("s4", 0.003, "VLDGSTR", "sp|P1_VAR_1_N10D|A"),
            psm("s5", 0.005, "VLNGSTR", "sp|P1|A"),
            psm("s6", 0.5, "CCCK", "sp|P1-2|A"),
        ]
        .into_iter()
        .collect()
    }

    #[test]
    fn members_get_specific_and_shared_evidence() {
        let rows = proteoform_table(&targets(), &index(), &Cutoffs::uniform(0.01), None);
        let m: HashMap<&str, &ProteoformRow> =
            rows.iter().map(|r| (r.member.as_str(), r)).collect();
        assert_eq!(rows[0].member, "sp|P1|A", "canonical first");
        let can = m["sp|P1|A"];
        assert_eq!(
            (can.specific_peptides, can.shared_peptides, can.psms),
            (2, 1, 2)
        );
        assert_eq!((can.evidence, can.best_psm.as_str()), ("specific", "s2"));
        assert!((can.confidence() - 0.998).abs() < 1e-12);
        let iso = m["sp|P1-2|A"];
        assert_eq!(
            (iso.specific_peptides, iso.shared_peptides, iso.psms),
            (1, 1, 1)
        );
        assert!(
            (iso.q_value - 0.004).abs() < 1e-12,
            "the q 0.5 PSM is not accepted"
        );
        assert!(rows.iter().all(|r| r.family == "P1"));
    }

    #[test]
    fn a_variant_peptide_another_familys_reference_holds_is_shared() {
        // VLDGSTR is the P1 variant's peptide, and the reference of family Q9
        // holds it too: the reference explains it, so the variant is not listed.
        let mut proteins = index().db.proteins.clone();
        proteins.retain(|p| !p.accession.starts_with("XXX"));
        proteins.push(Protein {
            accession: "sp|Q9|B".into(),
            description: String::new(),
            sequence: b"VLDGSTRK".to_vec(),
            annotation: None,
        });
        let idx = SearchIndex::from_target_db(&ProteinDb { proteins }, "XXX");
        let t: HashMap<_, _> = [psm("v", 0.001, "VLDGSTR", "sp|P1_VAR_1_N10D|A;sp|Q9|B")]
            .into_iter()
            .collect();
        let rows = proteoform_table(&t, &idx, &Cutoffs::uniform(0.01), None);
        let members: Vec<&str> = rows.iter().map(|r| r.member.as_str()).collect();
        assert_eq!(members, ["sp|Q9|B"]);
        assert_eq!(rows[0].evidence, "specific");
    }

    #[test]
    fn an_alternative_member_needs_specific_evidence() {
        // Only the shared peptide AAAK is accepted: the family is present,
        // shown through its canonical entry, but nothing supports the isoform.
        let only_shared: HashMap<_, _> = [psm("s1", 0.001, "AAAK", "sp|P1|A;sp|P1-2|A")]
            .into_iter()
            .collect();
        let rows = proteoform_table(&only_shared, &index(), &Cutoffs::uniform(0.01), None);
        let members: Vec<&str> = rows.iter().map(|r| r.member.as_str()).collect();
        assert_eq!(members, ["sp|P1|A"]);
        assert_eq!(rows[0].evidence, "shared");
    }

    #[test]
    fn a_variant_reports_its_alternative_explanation() {
        let rows = proteoform_table(&targets(), &index(), &Cutoffs::uniform(0.01), None);
        let var = rows.iter().find(|r| r.kind == EntryKind::Variant).unwrap();
        assert_eq!(var.feature, "VAR_1_N10D");
        assert_eq!(
            var.alternative,
            "deamidation-like mass; canonical counterpart VLNGSTR identified"
        );
    }

    #[test]
    fn quantity_sums_specific_peptides_per_column() {
        let quant = QuantSummary {
            columns: vec!["run1".into(), "run2".into()],
            values: vec![
                ("EEEK".into(), vec!["sp|P1|A".into()], 0, 10.0),
                ("VLNGSTR".into(), vec!["sp|P1|A".into()], 0, 5.0),
                ("EEEK".into(), vec!["sp|P1|A".into()], 1, 7.0),
                (
                    "AAAK".into(),
                    vec!["sp|P1|A".into(), "sp|P1-2|A".into()],
                    0,
                    100.0,
                ),
                ("CCCK".into(), vec!["sp|P1-2|A".into()], 1, 3.0),
            ],
        };
        let rows = proteoform_table(&targets(), &index(), &Cutoffs::uniform(0.01), Some(&quant));
        let m: HashMap<&str, &ProteoformRow> =
            rows.iter().map(|r| (r.member.as_str(), r)).collect();
        assert_eq!(
            m["sp|P1|A"].quantity,
            vec![15.0, 7.0],
            "shared AAAK is not counted"
        );
        assert_eq!(m["sp|P1-2|A"].quantity, vec![0.0, 3.0]);
    }
}
