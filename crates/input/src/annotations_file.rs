//! Annotation files (`--annotations <name>.annotations.parquet`): claims about
//! the reference proteins of any target database, one row per claim, read
//! into [`Annotation`]s and built into entries like UniProt's own.
//!
//! Columns (key-value metadata `andes.annotations.schema = 1`):
//!
//! | column | type | |
//! |---|---|---|
//! | `annotation_id` | utf8 | stable id, naming one claim of its reference across files |
//! | `reference_accession` | utf8 | accession (first header token) of the reference |
//! | `reference_md5` | utf8 | md5 (hex) of the reference sequence it was written against |
//! | `source_db` | utf8 | where the claim comes from (`COSMIC`, `dbSNP`, `sample_vcf`, …) |
//! | `edits` | list<struct<start, end, ref: utf8, alt: utf8>> | 0-based half-open, sorted, disjoint |
//! | `mods` | list<struct<position, residue: utf8, mass: f64, anchor: utf8>> | in the entry's coordinates |
//! | `processing` | struct<start, end> | a mature chain's span |
//!
//! Integer fields may be any integer type. A row needs edits, mods or
//! processing: edits make a variant entry, processing a chain entry, and mods
//! alone add sites to the reference itself. Unknown columns are ignored.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, AsArray, ListArray, StructArray};
use arrow::compute::cast;
use arrow::datatypes::{DataType, Float64Type, UInt32Type};
use md5::{Digest, Md5};
use model::{
    Annotation, ChainEnds, Edit, EntryKind, ModSite, PriorTier, Processing, Protein,
    ProteinAnnotation, ProteinDb, SiteAnchor,
};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::annotation::{build_entry, knockoff_annotations, primary_accession, KnockoffOptions};

/// The schema version this reader understands.
pub const ANNOTATIONS_SCHEMA: &str = "1";
const SCHEMA_KEY: &str = "andes.annotations.schema";

#[derive(thiserror::Error, Debug)]
pub enum AnnotationsError {
    #[error("{path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("{path}: {message}")]
    Format { path: String, message: String },
    #[error("annotation id {id} is given to different claims ({first} and {second})")]
    DuplicateId {
        id: String,
        first: String,
        second: String,
    },
}

/// Counts from loading annotation files, for the run log.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AnnotationLoadStats {
    pub files: usize,
    pub rows: usize,
    pub variant_entries: usize,
    pub chain_entries: usize,
    pub knockoff_entries: usize,
    /// Sites added to reference proteins (rows with mods only).
    pub sites: usize,
    /// Rows whose claim another row already made (merged into it).
    pub merged: usize,
    /// Rows whose claim an entry of the database already makes (a UniProt
    /// XML or PEFF entry): merged into that entry, which keeps the better tier.
    pub merged_into_database: usize,
    /// Single-residue variants for which no position qualified for a knockoff.
    pub variants_without_knockoff: usize,
    /// Rows left out, by reason.
    pub dropped: BTreeMap<&'static str, usize>,
}

/// One claim as read (from an annotation file or a PEFF header), before it
/// is checked against its reference.
pub(crate) struct Row {
    pub(crate) id: String,
    pub(crate) reference: String,
    pub(crate) md5: String,
    /// `(start, end, ref, alt)`, 0-based half-open.
    pub(crate) edits: Vec<(u32, u32, String, String)>,
    /// `(position, residue, mass, anchor)`.
    pub(crate) mods: Vec<(u32, String, f64, String)>,
    pub(crate) processing: Option<(u32, u32)>,
    pub(crate) source: String,
    pub(crate) tier: PriorTier,
}

/// The tier of a claim from what its source says: the sample's own call
/// (`genotype` not reference), then a common allele (`population_af` ≥ 0.01),
/// then curated or experimental evidence or a benign clinical class, else a
/// catalogue claim.
pub(crate) fn tier_of(
    genotype: &str,
    population_af: Option<f64>,
    evidence: &str,
    clinical_significance: &str,
) -> PriorTier {
    let g = genotype.trim();
    if !g.is_empty() && !matches!(g, "0/0" | "0|0" | "./." | ".|." | ".") {
        PriorTier::Sample
    } else if population_af.is_some_and(|af| af >= 0.01) {
        PriorTier::Common
    } else if matches!(evidence.trim(), "experimental" | "curated")
        || clinical_significance
            .to_ascii_lowercase()
            .contains("benign")
    {
        PriorTier::Curated
    } else {
        PriorTier::Catalogue
    }
}

/// Read `paths` and add their claims to `db` (see [`add_claims`]).
pub fn load_annotations(
    paths: &[impl AsRef<Path>],
    db: &mut ProteinDb,
    knockoffs: Option<KnockoffOptions>,
) -> Result<AnnotationLoadStats, AnnotationsError> {
    let mut stats = AnnotationLoadStats::default();
    let mut rows = Vec::new();
    for path in paths {
        stats.files += 1;
        rows.extend(read_rows(path.as_ref())?);
    }
    add_claims(rows, db, knockoffs, &mut stats)?;
    Ok(stats)
}

/// Check `rows` against their reference proteins in `db` and add them: sites
/// to the references, then (in reference order, after every existing protein)
/// each reference's variant entries, their knockoffs (with `knockoffs`) and
/// its chain entries. Only `Reference` proteins can be annotated.
pub(crate) fn add_claims(
    rows: Vec<Row>,
    db: &mut ProteinDb,
    knockoffs: Option<KnockoffOptions>,
    stats: &mut AnnotationLoadStats,
) -> Result<(), AnnotationsError> {
    let index: HashMap<String, usize> = db
        .proteins
        .iter()
        .enumerate()
        .filter(|(_, p)| p.entry_kind() == EntryKind::Reference)
        .map(|(i, p)| (p.accession.clone(), i))
        .collect();
    let md5s: HashMap<usize, String> = index
        .values()
        .map(|&i| (i, sequence_md5(&db.proteins[i].sequence)))
        .collect();

    // Claims per reference in row order, each claim's content key by
    // (reference, id), and the content keys seen.
    let mut claims: BTreeMap<usize, Vec<Annotation>> = BTreeMap::new();
    let mut sites: BTreeMap<usize, Vec<ModSite>> = BTreeMap::new();
    let mut listed: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut ids: HashMap<String, String> = HashMap::new();
    // Where each content key's claim went, so a better tier can replace a
    // merged claim's.
    let mut seen: HashMap<String, (usize, usize)> = HashMap::new();
    for row in rows {
        stats.rows += 1;
        let skip = |stats: &mut AnnotationLoadStats, reason: &'static str| {
            *stats.dropped.entry(reason).or_insert(0) += 1;
        };
        let Some(&r) = index.get(&row.reference) else {
            skip(stats, "reference not in the database");
            continue;
        };
        if !md5s[&r].eq_ignore_ascii_case(&row.md5) {
            skip(stats, "reference sequence differs (md5)");
            continue;
        }
        let reference = &db.proteins[r];
        let claim = match claim_of(&row, reference) {
            Ok(c) => c,
            Err(reason) => {
                skip(stats, reason);
                continue;
            }
        };
        listed.entry(r).or_default().extend(
            claim
                .edits
                .iter()
                .flat_map(|e| e.start as usize..(e.end as usize).max(e.start as usize + 1)),
        );
        let key = format!(
            "{}\t{:?}\t{:?}\t{:?}",
            r, claim.edits, claim.mods, claim.processing
        );
        // An id names one claim of its reference.
        let id = format!("{}\t{}", row.reference, row.id);
        match ids.get(&id) {
            Some(k) if *k != key => {
                return Err(AnnotationsError::DuplicateId {
                    id: row.id.clone(),
                    first: k.clone(),
                    second: key,
                })
            }
            _ => {
                ids.insert(id, key.clone());
            }
        }
        if let Some(&(cr, ci)) = seen.get(&key) {
            stats.merged += 1;
            if let Some(kept) = claims.get_mut(&cr).and_then(|c| c.get_mut(ci)) {
                if claim.tier < kept.tier {
                    kept.tier = claim.tier;
                    kept.source = claim.source;
                }
            }
            continue;
        }
        if claim.edits.is_empty() && claim.processing.is_none() {
            seen.insert(key, (usize::MAX, 0));
            sites.entry(r).or_default().extend(claim.mods);
        } else {
            let list = claims.entry(r).or_default();
            seen.insert(key, (r, list.len()));
            list.push(claim);
        }
    }

    for (r, extra) in sites {
        stats.sites += extra.len();
        let p = &mut db.proteins[r];
        let mut info = p
            .annotation
            .as_deref()
            .cloned()
            .unwrap_or_else(|| ProteinAnnotation::new(Vec::new(), None));
        let mut all = info.mod_sites.clone();
        all.extend(extra);
        info.mod_sites = ProteinAnnotation::new(all, None).mod_sites;
        p.annotation = Some(Arc::new(info));
    }

    // The database's own alternative entries (not knockoffs), by what they
    // read: a claim equal to one of them joins it.
    let mut existing: HashMap<(Arc<str>, EntryKind, Vec<u8>), usize> = HashMap::new();
    // Each such entry with its knockoffs, by (family, feature).
    let mut members: HashMap<(Arc<str>, Arc<str>), Vec<usize>> = HashMap::new();
    for (i, p) in db.proteins.iter().enumerate() {
        let Some(o) = p.origin() else {
            continue;
        };
        if o.kind == EntryKind::Reference {
            continue;
        }
        match &o.knockoff_of {
            Some(of) => members
                .entry((o.family.clone(), of.clone()))
                .or_default()
                .push(i),
            None => {
                existing.insert((o.family.clone(), o.kind, p.sequence.clone()), i);
                members
                    .entry((o.family.clone(), o.feature.clone()))
                    .or_default()
                    .push(i);
            }
        }
    }

    for (r, list) in claims {
        let reference = db.proteins[r].clone();
        let mut new = Vec::with_capacity(list.len());
        for a in list {
            let Some(entry) = build_entry(&reference, &a) else {
                continue;
            };
            let o = entry.origin().expect("a built entry has an origin");
            match existing.get(&(o.family.clone(), o.kind, entry.sequence.clone())) {
                Some(&i) => {
                    stats.merged_into_database += 1;
                    let key = {
                        let e = db.proteins[i].origin().expect("an entry with an origin");
                        (e.family.clone(), e.feature.clone())
                    };
                    promote(db, &members[&key], a.tier, a.source.clone());
                }
                None => new.push(a),
            }
        }
        let (variants, chains): (Vec<Annotation>, Vec<Annotation>) =
            new.into_iter().partition(|a| a.kind == EntryKind::Variant);
        let ko = match knockoffs.filter(|k| k.per_variant > 0) {
            Some(ko) => {
                let (ko, untested) = knockoff_annotations(
                    &reference.sequence,
                    primary_accession(&reference.accession),
                    &variants,
                    &listed[&r],
                    ko,
                );
                stats.variants_without_knockoff += untested;
                ko
            }
            None => Vec::new(),
        };
        for (group, count) in [
            (&variants, &mut stats.variant_entries),
            (&ko, &mut stats.knockoff_entries),
            (&chains, &mut stats.chain_entries),
        ] {
            for a in group {
                if let Some(p) = build_entry(&reference, a) {
                    db.proteins.push(p);
                    *count += 1;
                }
            }
        }
    }
    Ok(())
}

/// Give the entries at `indices` (an entry and its knockoffs) `tier` and
/// `source` when `tier` is better than the first one's own.
fn promote(db: &mut ProteinDb, indices: &[usize], tier: PriorTier, source: Option<Arc<str>>) {
    let better = indices
        .first()
        .and_then(|&i| db.proteins[i].origin())
        .is_some_and(|o| tier < o.tier);
    if !better {
        return;
    }
    for &i in indices {
        let p = &mut db.proteins[i];
        let mut info = p
            .annotation
            .as_deref()
            .cloned()
            .expect("an entry with an origin");
        if let Some(po) = info.origin.as_mut() {
            po.tier = tier;
            po.source = source.clone();
        }
        p.annotation = Some(Arc::new(info));
    }
}

/// The claim a row makes about `reference`, or why it is left out.
fn claim_of(row: &Row, reference: &Protein) -> Result<Annotation, &'static str> {
    let seq = &reference.sequence;
    let mut edits = Vec::with_capacity(row.edits.len());
    for (start, end, rf, alt) in &row.edits {
        let (s, e) = (*start as usize, *end as usize);
        if s > e || e > seq.len() {
            return Err("edit outside the reference");
        }
        if rf.as_bytes() != &seq[s..e] {
            return Err("edit ref does not match the reference");
        }
        if !alt.bytes().all(|c| c.is_ascii_uppercase()) {
            return Err("edit alt is not a residue string");
        }
        edits.push(Edit {
            start: *start,
            end: *end,
            alt: alt.as_bytes().to_vec(),
        });
    }
    if edits.windows(2).any(|w| w[1].start < w[0].end) {
        return Err("edits not sorted and disjoint");
    }
    let il = |c: &u8| if *c == b'I' { b'L' } else { *c };
    if !edits.is_empty()
        && edits.iter().all(|e| {
            e.alt.len() == (e.end - e.start) as usize
                && e.alt
                    .iter()
                    .map(il)
                    .eq(seq[e.start as usize..e.end as usize].iter().map(il))
        })
    {
        return Err("edits only swap I and L");
    }
    if !edits.is_empty() && row.processing.is_some() {
        return Err("edits on a mature chain are not supported");
    }
    let processing = match row.processing {
        Some((start, end)) => {
            let (s, e) = (start as usize, end as usize);
            if s >= e || e > seq.len() {
                return Err("processing span outside the reference");
            }
            let ends = ChainEnds {
                n_term: s > 0 && !(s == 1 && seq[0] == b'M'),
                c_term: e < seq.len(),
            };
            if !ends.n_term && !ends.c_term {
                return Err("processing span has no new terminus");
            }
            Some(Processing { start, end, ends })
        }
        None => None,
    };
    if edits.is_empty() && processing.is_none() && row.mods.is_empty() {
        return Err("no edits, mods or processing");
    }
    // The sequence mods refer to: the reference with the edits, or the chain.
    let entry_seq: Vec<u8> = match processing {
        Some(p) => seq[p.start as usize..p.end as usize].to_vec(),
        None => {
            let mut out = Vec::with_capacity(seq.len());
            let mut at = 0;
            for e in &edits {
                out.extend_from_slice(&seq[at..e.start as usize]);
                out.extend_from_slice(&e.alt);
                at = e.end as usize;
            }
            out.extend_from_slice(&seq[at..]);
            out
        }
    };
    let mut mods = Vec::with_capacity(row.mods.len());
    for (position, residue, mass, anchor) in &row.mods {
        let p = *position as usize;
        if p >= entry_seq.len() || residue.as_bytes() != [entry_seq[p]] {
            return Err("mod residue does not match");
        }
        let anchor = match anchor.as_str() {
            "" | "residue" => SiteAnchor::Residue,
            "protein_n_term" => SiteAnchor::NTerminus,
            "protein_c_term" => SiteAnchor::CTerminus,
            _ => return Err("unknown mod anchor"),
        };
        mods.push(ModSite {
            position: *position,
            mass_delta: *mass,
            anchor,
        });
    }
    let kind = if !edits.is_empty() {
        EntryKind::Variant
    } else if processing.is_some() {
        EntryKind::Chain
    } else {
        EntryKind::Reference
    };
    let mut a = Annotation::new(&row.id, kind);
    a.edits = edits;
    a.mods = mods;
    a.processing = processing;
    a.source = (!row.source.is_empty()).then(|| Arc::from(row.source.as_str()));
    a.tier = row.tier;
    Ok(a)
}

/// md5 (lowercase hex) of a sequence's uppercase residues: the
/// `reference_md5` an annotation file gives for its reference.
pub fn sequence_md5(seq: &[u8]) -> String {
    let upper: Vec<u8> = seq.iter().map(u8::to_ascii_uppercase).collect();
    Md5::digest(&upper)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn read_rows(path: &Path) -> Result<Vec<Row>, AnnotationsError> {
    let name = path.display().to_string();
    let bad = |message: String| AnnotationsError::Format {
        path: name.clone(),
        message,
    };
    let file = File::open(path).map_err(|source| AnnotationsError::Io {
        path: name.clone(),
        source,
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| bad(e.to_string()))?;
    let version = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|k| k.key == SCHEMA_KEY))
        .and_then(|k| k.value.clone());
    if version.as_deref() != Some(ANNOTATIONS_SCHEMA) {
        return Err(bad(format!(
            "not an annotations file of schema {ANNOTATIONS_SCHEMA} (metadata {SCHEMA_KEY} = {})",
            version.as_deref().unwrap_or("missing")
        )));
    }
    let reader = builder.build().map_err(|e| bad(e.to_string()))?;
    let mut rows = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| bad(e.to_string()))?;
        let column = |n: &str| batch.column_by_name(n).cloned();
        let required = |n: &str| column(n).ok_or_else(|| bad(format!("missing column {n}")));
        let id = utf8(&required("annotation_id")?).map_err(bad)?;
        let reference = utf8(&required("reference_accession")?).map_err(bad)?;
        let md5 = utf8(&required("reference_md5")?).map_err(bad)?;
        let source = utf8(&required("source_db")?).map_err(bad)?;
        let optional_utf8 = |n: &str| -> Result<Vec<Option<String>>, AnnotationsError> {
            match column(n) {
                Some(c) => utf8(&c).map_err(bad),
                None => Ok(vec![None; batch.num_rows()]),
            }
        };
        let genotype = optional_utf8("genotype")?;
        let evidence = optional_utf8("evidence")?;
        let clinical = optional_utf8("clinical_significance")?;
        let population_af = match column("population_af") {
            Some(c) => {
                let c = cast(&c, &DataType::Float64).map_err(|e| bad(e.to_string()))?;
                let c = c.as_primitive::<Float64Type>();
                (0..c.len())
                    .map(|i| (!c.is_null(i)).then(|| c.value(i)))
                    .collect()
            }
            None => vec![None; batch.num_rows()],
        };
        let edits = match column("edits") {
            Some(c) => list_of_structs(&c, &["start", "end", "ref", "alt"]).map_err(bad)?,
            None => vec![None; batch.num_rows()],
        };
        let mods = match column("mods") {
            Some(c) => {
                list_of_structs(&c, &["position", "residue", "mass", "anchor"]).map_err(bad)?
            }
            None => vec![None; batch.num_rows()],
        };
        let processing = match column("processing") {
            Some(c) => processing_spans(&c).map_err(bad)?,
            None => vec![None; batch.num_rows()],
        };
        for i in 0..batch.num_rows() {
            let (Some(id), Some(reference), Some(md5)) = (&id[i], &reference[i], &md5[i]) else {
                return Err(bad(format!("row {i}: a required value is null")));
            };
            let edits = match &edits[i] {
                Some(s) => (0..s.len())
                    .map(|j| -> Result<_, String> {
                        Ok((
                            u32_at(&s.fields[0], j)?,
                            u32_at(&s.fields[1], j)?,
                            str_at(&s.fields[2], j),
                            str_at(&s.fields[3], j),
                        ))
                    })
                    .collect::<Result<Vec<_>, String>>()
                    .map_err(bad)?,
                None => Vec::new(),
            };
            let mods = match &mods[i] {
                Some(s) => (0..s.len())
                    .map(|j| -> Result<_, String> {
                        Ok((
                            u32_at(&s.fields[0], j)?,
                            str_at(&s.fields[1], j),
                            f64_at(&s.fields[2], j)?,
                            str_at(&s.fields[3], j),
                        ))
                    })
                    .collect::<Result<Vec<_>, String>>()
                    .map_err(bad)?,
                None => Vec::new(),
            };
            rows.push(Row {
                id: id.clone(),
                reference: reference.clone(),
                md5: md5.clone(),
                edits,
                mods,
                processing: processing[i],
                source: source[i].clone().unwrap_or_default(),
                tier: tier_of(
                    genotype[i].as_deref().unwrap_or(""),
                    population_af[i],
                    evidence[i].as_deref().unwrap_or(""),
                    clinical[i].as_deref().unwrap_or(""),
                ),
            });
        }
    }
    Ok(rows)
}

/// The fields of one row's list of structs, cast to the types the reader
/// uses, each of the same length.
#[derive(Clone)]
struct Fields {
    fields: Vec<ArrayRef>,
}

impl Fields {
    fn len(&self) -> usize {
        self.fields.first().map_or(0, |f| f.len())
    }
}

fn utf8(col: &ArrayRef) -> Result<Vec<Option<String>>, String> {
    let c = cast(col, &DataType::Utf8).map_err(|e| e.to_string())?;
    let s = c.as_string::<i32>();
    Ok((0..s.len())
        .map(|i| (!s.is_null(i)).then(|| s.value(i).to_string()))
        .collect())
}

/// Per row, the struct fields `names` of a list<struct> column (absent fields
/// read as null), or `None` for a null or empty list.
fn list_of_structs(col: &ArrayRef, names: &[&str]) -> Result<Vec<Option<Fields>>, String> {
    let list = cast(
        col,
        &DataType::List(Arc::new(arrow::datatypes::Field::new(
            "item",
            match col.data_type() {
                DataType::List(f) | DataType::LargeList(f) => f.data_type().clone(),
                other => return Err(format!("expected a list column, found {other}")),
            },
            true,
        ))),
    )
    .map_err(|e| e.to_string())?;
    let list: &ListArray = list.as_list::<i32>();
    let mut out = Vec::with_capacity(list.len());
    for i in 0..list.len() {
        if list.is_null(i) || list.value_length(i) == 0 {
            out.push(None);
            continue;
        }
        let item = list.value(i);
        let st: &StructArray = item
            .as_any()
            .downcast_ref::<StructArray>()
            .ok_or("expected a list of structs")?;
        let fields = names
            .iter()
            .map(|n| {
                st.column_by_name(n)
                    .cloned()
                    .unwrap_or_else(|| arrow::array::new_null_array(&DataType::Utf8, st.len()))
            })
            .collect();
        out.push(Some(Fields { fields }));
    }
    Ok(out)
}

fn processing_spans(col: &ArrayRef) -> Result<Vec<Option<(u32, u32)>>, String> {
    let st: &StructArray = col
        .as_any()
        .downcast_ref::<StructArray>()
        .ok_or("processing: expected a struct column")?;
    let (Some(s), Some(e)) = (st.column_by_name("start"), st.column_by_name("end")) else {
        return Err("processing: needs start and end".into());
    };
    (0..st.len())
        .map(|i| {
            if st.is_null(i) {
                Ok(None)
            } else {
                Ok(Some((u32_at(s, i)?, u32_at(e, i)?)))
            }
        })
        .collect()
}

fn u32_at(a: &ArrayRef, i: usize) -> Result<u32, String> {
    let c = cast(a, &DataType::UInt32).map_err(|e| e.to_string())?;
    let c = c.as_primitive::<UInt32Type>();
    if c.is_null(i) {
        return Err("a required integer is null or negative".into());
    }
    Ok(c.value(i))
}

fn f64_at(a: &ArrayRef, i: usize) -> Result<f64, String> {
    let c = cast(a, &DataType::Float64).map_err(|e| e.to_string())?;
    let c = c.as_primitive::<Float64Type>();
    if c.is_null(i) {
        return Err("a required mass is null".into());
    }
    Ok(c.value(i))
}

fn str_at(a: &ArrayRef, i: usize) -> String {
    match cast(a, &DataType::Utf8) {
        Ok(c) => {
            let c = c.as_string::<i32>();
            if c.is_null(i) {
                String::new()
            } else {
                c.value(i).to_string()
            }
        }
        Err(_) => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Float64Array, Int64Array, RecordBatch, StringArray};
    use arrow::buffer::{NullBuffer, OffsetBuffer};
    use arrow::datatypes::{Field, Fields as ArrowFields, Schema};
    use model::Enzyme;
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;
    use parquet::format::KeyValue;

    type EditRow<'a> = (i64, i64, &'a str, &'a str);
    type ModRow<'a> = (i64, &'a str, f64, &'a str);

    /// One annotations row: id, reference, md5, edits, mods, processing.
    struct In<'a> {
        id: &'a str,
        reference: &'a str,
        md5: String,
        edits: Vec<EditRow<'a>>,
        mods: Vec<ModRow<'a>>,
        processing: Option<(i64, i64)>,
        genotype: Option<&'a str>,
    }

    fn list_of(fields: Vec<(Arc<Field>, ArrayRef)>, lengths: Vec<usize>) -> ArrayRef {
        let st = StructArray::from(fields);
        let item = Arc::new(Field::new("item", st.data_type().clone(), true));
        Arc::new(ListArray::new(
            item,
            OffsetBuffer::from_lengths(lengths),
            Arc::new(st),
            None,
        ))
    }

    fn write(path: &Path, rows: &[In], version: &str) {
        let s = |f: &dyn Fn(&In) -> String| -> ArrayRef {
            Arc::new(StringArray::from(rows.iter().map(f).collect::<Vec<_>>()))
        };
        let e: Vec<&EditRow> = rows.iter().flat_map(|r| &r.edits).collect();
        let edits = list_of(
            vec![
                (
                    Arc::new(Field::new("start", DataType::Int64, false)),
                    Arc::new(Int64Array::from(e.iter().map(|x| x.0).collect::<Vec<_>>()))
                        as ArrayRef,
                ),
                (
                    Arc::new(Field::new("end", DataType::Int64, false)),
                    Arc::new(Int64Array::from(e.iter().map(|x| x.1).collect::<Vec<_>>())),
                ),
                (
                    Arc::new(Field::new("ref", DataType::Utf8, false)),
                    Arc::new(StringArray::from(e.iter().map(|x| x.2).collect::<Vec<_>>())),
                ),
                (
                    Arc::new(Field::new("alt", DataType::Utf8, false)),
                    Arc::new(StringArray::from(e.iter().map(|x| x.3).collect::<Vec<_>>())),
                ),
            ],
            rows.iter().map(|r| r.edits.len()).collect(),
        );
        let m: Vec<&ModRow> = rows.iter().flat_map(|r| &r.mods).collect();
        let mods = list_of(
            vec![
                (
                    Arc::new(Field::new("position", DataType::Int64, false)),
                    Arc::new(Int64Array::from(m.iter().map(|x| x.0).collect::<Vec<_>>()))
                        as ArrayRef,
                ),
                (
                    Arc::new(Field::new("residue", DataType::Utf8, false)),
                    Arc::new(StringArray::from(m.iter().map(|x| x.1).collect::<Vec<_>>())),
                ),
                (
                    Arc::new(Field::new("mass", DataType::Float64, false)),
                    Arc::new(Float64Array::from(
                        m.iter().map(|x| x.2).collect::<Vec<_>>(),
                    )),
                ),
                (
                    Arc::new(Field::new("anchor", DataType::Utf8, false)),
                    Arc::new(StringArray::from(m.iter().map(|x| x.3).collect::<Vec<_>>())),
                ),
            ],
            rows.iter().map(|r| r.mods.len()).collect(),
        );
        let proc_fields = ArrowFields::from(vec![
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]);
        let processing: ArrayRef = Arc::new(
            StructArray::try_new(
                proc_fields,
                vec![
                    Arc::new(Int64Array::from(
                        rows.iter()
                            .map(|r| r.processing.map_or(0, |p| p.0))
                            .collect::<Vec<_>>(),
                    )),
                    Arc::new(Int64Array::from(
                        rows.iter()
                            .map(|r| r.processing.map_or(0, |p| p.1))
                            .collect::<Vec<_>>(),
                    )),
                ],
                Some(NullBuffer::from(
                    rows.iter()
                        .map(|r| r.processing.is_some())
                        .collect::<Vec<_>>(),
                )),
            )
            .unwrap(),
        );
        let columns: Vec<(&str, ArrayRef)> = vec![
            ("annotation_id", s(&|r| r.id.to_string())),
            ("reference_accession", s(&|r| r.reference.to_string())),
            ("reference_md5", s(&|r| r.md5.clone())),
            ("source_db", s(&|_| "custom".to_string())),
            ("edits", edits),
            ("mods", mods),
            ("processing", processing),
            (
                "genotype",
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.genotype).collect::<Vec<_>>(),
                )),
            ),
        ];
        let schema = Arc::new(Schema::new(
            columns
                .iter()
                .map(|(n, a)| Field::new(*n, a.data_type().clone(), true))
                .collect::<Vec<_>>(),
        ));
        let batch = RecordBatch::try_new(
            schema.clone(),
            columns.into_iter().map(|(_, a)| a).collect(),
        )
        .unwrap();
        let props = WriterProperties::builder()
            .set_key_value_metadata(Some(vec![KeyValue::new(
                SCHEMA_KEY.to_string(),
                version.to_string(),
            )]))
            .build();
        let mut w = ArrowWriter::try_new(File::create(path).unwrap(), schema, Some(props)).unwrap();
        w.write(&batch).unwrap();
        w.close().unwrap();
    }

    /// Eight tryptic spans `GGGGDLGGGK`, D at 0-based 4, 14, …, 74.
    fn db() -> ProteinDb {
        ProteinDb {
            proteins: vec![Protein {
                accession: "sp|P00009|KO_HUMAN".into(),
                description: "Knockoff test".into(),
                sequence: b"GGGGDLGGGK".repeat(8),
                annotation: None,
            }],
        }
    }

    fn md5() -> String {
        sequence_md5(&b"GGGGDLGGGK".repeat(8))
    }

    const KO: KnockoffOptions = KnockoffOptions {
        per_variant: 4,
        enzyme: Enzyme::Trypsin,
    };

    fn row<'a>(id: &'a str) -> In<'a> {
        In {
            id,
            reference: "sp|P00009|KO_HUMAN",
            md5: md5(),
            edits: Vec::new(),
            mods: Vec::new(),
            processing: None,
            genotype: None,
        }
    }

    #[test]
    fn rows_become_variants_knockoffs_chains_and_sites() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.annotations.parquet");
        write(
            &path,
            &[
                In {
                    edits: vec![(4, 5, "D", "N")],
                    ..row("V1")
                },
                In {
                    processing: Some((10, 80)),
                    ..row("C1")
                },
                In {
                    mods: vec![(9, "K", 114.042927, "residue")],
                    ..row("S1")
                },
            ],
            ANNOTATIONS_SCHEMA,
        );
        let mut db = db();
        let stats = load_annotations(&[&path], &mut db, Some(KO)).unwrap();
        assert_eq!(
            (
                stats.variant_entries,
                stats.knockoff_entries,
                stats.chain_entries,
                stats.sites
            ),
            (1, 4, 1, 1)
        );
        assert_eq!(
            db.proteins[0].mod_sites().len(),
            1,
            "the site is on the reference"
        );
        let v = &db.proteins[1];
        assert_eq!(v.accession, "sp|P00009_V1|KO_HUMAN");
        assert_eq!(v.sequence[4], b'N');
        assert_eq!(
            v.mod_sites().len(),
            1,
            "the variant inherits the reference's site"
        );
        let o = v.origin().unwrap();
        assert_eq!((o.kind, &*o.family), (EntryKind::Variant, "P00009"));
        assert!(db.proteins[2..6].iter().all(|p| p.is_knockoff()));
        let c = &db.proteins[6];
        assert_eq!(c.entry_kind(), EntryKind::Chain);
        assert_eq!(c.sequence.len(), 70);
        assert!(c.description.ends_with("CHAIN=C1:11-80"));
    }

    #[test]
    fn rows_that_do_not_fit_are_counted_not_loaded() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.annotations.parquet");
        write(
            &path,
            &[
                In {
                    md5: "0".repeat(32),
                    edits: vec![(4, 5, "D", "N")],
                    ..row("bad_md5")
                },
                In {
                    edits: vec![(4, 5, "E", "N")],
                    ..row("bad_ref")
                },
                In {
                    reference: "sp|Q0|X",
                    edits: vec![(4, 5, "D", "N")],
                    ..row("no_ref")
                },
                In {
                    edits: vec![(5, 6, "L", "I")],
                    ..row("il")
                },
                In {
                    mods: vec![(4, "S", 79.966331, "residue")],
                    ..row("bad_site")
                },
                In { ..row("empty") },
            ],
            ANNOTATIONS_SCHEMA,
        );
        let mut db = db();
        let stats = load_annotations(&[&path], &mut db, None).unwrap();
        assert_eq!(db.proteins.len(), 1);
        assert_eq!(
            stats.dropped.values().sum::<usize>(),
            6,
            "{:?}",
            stats.dropped
        );
        assert_eq!(stats.dropped["reference sequence differs (md5)"], 1);
        assert_eq!(stats.dropped["edits only swap I and L"], 1);
    }

    #[test]
    fn files_merge_by_content_and_ids_stay_unique() {
        let dir = tempfile::tempdir().unwrap();
        let (a, b, c) = (
            dir.path().join("a.parquet"),
            dir.path().join("b.parquet"),
            dir.path().join("c.parquet"),
        );
        write(
            &a,
            &[In {
                edits: vec![(4, 5, "D", "N")],
                ..row("V1")
            }],
            ANNOTATIONS_SCHEMA,
        );
        write(
            &b,
            &[In {
                edits: vec![(4, 5, "D", "N")],
                ..row("rs1")
            }],
            ANNOTATIONS_SCHEMA,
        );
        write(
            &c,
            &[In {
                edits: vec![(14, 15, "D", "N")],
                ..row("V1")
            }],
            ANNOTATIONS_SCHEMA,
        );
        let mut d = db();
        let stats = load_annotations(&[&a, &b], &mut d, None).unwrap();
        assert_eq!((stats.variant_entries, stats.merged), (1, 1));
        let mut d = db();
        assert!(matches!(
            load_annotations(&[&a, &c], &mut d, None),
            Err(AnnotationsError::DuplicateId { .. })
        ));
    }

    #[test]
    fn a_file_without_the_schema_key_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.parquet");
        write(
            &path,
            &[In {
                edits: vec![(4, 5, "D", "N")],
                ..row("V1")
            }],
            "0",
        );
        let mut d = db();
        assert!(matches!(
            load_annotations(&[&path], &mut d, None),
            Err(AnnotationsError::Format { .. })
        ));
    }

    /// A UniProt XML entry with the `GGGGDLGGGK` x8 sequence and one curated
    /// variant, D5N (`VAR_000009`).
    fn ko_xml() -> String {
        let seq = String::from_utf8(b"GGGGDLGGGK".repeat(8)).unwrap();
        format!(
            r#"<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot">
<entry dataset="Swiss-Prot" version="1">
  <accession>P00009</accession>
  <name>KO_HUMAN</name>
  <protein><recommendedName><fullName>Knockoff test</fullName></recommendedName></protein>
  <organism><name type="scientific">Homo sapiens</name></organism>
  <proteinExistence type="evidence at protein level"/>
  <feature type="sequence variant" id="VAR_000009" evidence="1">
    <original>D</original><variation>N</variation>
    <location><position position="5"/></location>
  </feature>
  <evidence key="1" type="ECO:0000269"/>
  <sequence length="80" mass="1" checksum="0" modified="2000-01-01" version="1">{seq}</sequence>
</entry>
</uniprot>"#
        )
    }

    #[test]
    fn a_variant_from_a_file_equals_the_same_variant_from_uniprot_xml() {
        let xml = ko_xml();
        let opts = crate::UniProtOptions {
            variants: true,
            knockoffs: Some(KO),
            ..crate::UniProtOptions::default()
        };
        let (xml_db, _) = crate::load_uniprot_xml(xml.as_bytes(), &opts).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.annotations.parquet");
        write(
            &path,
            &[In {
                edits: vec![(4, 5, "D", "N")],
                ..row("VAR_000009_D5N")
            }],
            ANNOTATIONS_SCHEMA,
        );
        let mut file_db = db();
        load_annotations(&[&path], &mut file_db, Some(KO)).unwrap();
        assert_eq!(xml_db.proteins.len(), file_db.proteins.len());
        for (x, f) in xml_db.proteins.iter().zip(&file_db.proteins).skip(1) {
            assert_eq!(x.accession, f.accession);
            assert_eq!(x.sequence, f.sequence);
            // The same claim, apart from where it is said to come from.
            let strip = |p: &Protein| {
                let mut a = p.annotation.as_deref().cloned().unwrap();
                if let Some(o) = a.origin.as_mut() {
                    o.source = None;
                    o.tier = PriorTier::Catalogue;
                }
                a
            };
            assert_eq!(strip(x), strip(f), "{}", x.accession);
        }
    }

    #[test]
    fn tiers_follow_the_first_rule_that_holds() {
        assert_eq!(tier_of("0/1", Some(0.5), "curated", ""), PriorTier::Sample);
        assert_eq!(tier_of("0/0", Some(0.02), "", ""), PriorTier::Common);
        assert_eq!(
            tier_of("", Some(0.001), "", "Likely benign"),
            PriorTier::Curated
        );
        assert_eq!(tier_of("./.", None, "experimental", ""), PriorTier::Curated);
        assert_eq!(tier_of("", None, "any", "pathogenic"), PriorTier::Catalogue);
    }

    #[test]
    fn a_merged_claim_keeps_the_best_tier() {
        let dir = tempfile::tempdir().unwrap();
        let (a, b) = (dir.path().join("a.parquet"), dir.path().join("b.parquet"));
        write(
            &a,
            &[In {
                edits: vec![(4, 5, "D", "N")],
                ..row("cat")
            }],
            ANNOTATIONS_SCHEMA,
        );
        write(
            &b,
            &[In {
                edits: vec![(4, 5, "D", "N")],
                genotype: Some("0/1"),
                ..row("own")
            }],
            ANNOTATIONS_SCHEMA,
        );
        let mut d = db();
        load_annotations(&[&a, &b], &mut d, None).unwrap();
        let v = d.proteins[1].origin().unwrap();
        assert_eq!((v.tier, &*v.feature), (PriorTier::Sample, "cat"));
    }

    #[test]
    fn a_claim_the_database_already_makes_joins_its_entry() {
        let opts = crate::UniProtOptions {
            variants: true,
            knockoffs: Some(KO),
            ..crate::UniProtOptions::default()
        };
        let (mut xml_db, _) = crate::load_uniprot_xml(ko_xml().as_bytes(), &opts).unwrap();
        let before = xml_db.proteins.len();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.parquet");
        let reference = xml_db.proteins[0].accession.clone();
        write(
            &path,
            &[In {
                reference: &reference,
                edits: vec![(4, 5, "D", "N")],
                genotype: Some("0/1"),
                ..row("own")
            }],
            ANNOTATIONS_SCHEMA,
        );
        let stats = load_annotations(&[&path], &mut xml_db, Some(KO)).unwrap();
        assert_eq!(
            (
                stats.merged_into_database,
                stats.variant_entries,
                stats.knockoff_entries
            ),
            (1, 0, 0)
        );
        assert_eq!(
            xml_db.proteins.len(),
            before,
            "no second entry, no second knockoffs"
        );
        for p in &xml_db.proteins[1..] {
            let o = p.origin().unwrap();
            assert_eq!(o.tier, PriorTier::Sample, "{}", p.accession);
            assert_eq!(o.source.as_deref(), Some("custom"));
        }
        // A weaker claim leaves the entry as it is.
        let (mut xml_db, _) = crate::load_uniprot_xml(ko_xml().as_bytes(), &opts).unwrap();
        write(
            &path,
            &[In {
                reference: &reference,
                edits: vec![(4, 5, "D", "N")],
                ..row("cat")
            }],
            ANNOTATIONS_SCHEMA,
        );
        load_annotations(&[&path], &mut xml_db, Some(KO)).unwrap();
        assert_eq!(
            xml_db.proteins[1].origin().unwrap().tier,
            PriorTier::Curated
        );
    }
}
