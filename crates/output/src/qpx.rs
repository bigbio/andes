use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, Float64Builder, Int32Builder, ListBuilder, StringBuilder,
    StructBuilder,
};
use arrow::datatypes::{DataType, Field, Fields, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;

use model::mass::{ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::tolerance::{PrecursorTolerance, Tolerance};
use search::candidate_gen::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::search_index::SearchIndex;
use search::search_params::SearchParams;

use crate::percolator::PercolatorPsm;
use crate::pin::{
    exp_mass_text, feature_value_type, format_feature_value, format_spec_id, psm_feature_values,
    spec_id_run,
};
use crate::quant_out::ModRecord;
use crate::row_context::{emitted_rows, resolve_accession};

const QPX_VERSION: &str = "1.0";
const CREATOR: &str = "andes";

pub fn omitted_psm_rows(
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    rescore_decoys: Option<&HashSet<String>>,
    excluded: &HashSet<String>,
) -> std::io::Result<HashSet<String>> {
    let Some(m) = rescore else {
        return Ok(excluded.clone());
    };
    let mut rows = rescore_rows(spectra, queues, candidates);
    rows.retain(|r| {
        !excluded.contains(&r.spec_id)
            || m.contains_key(&r.spec_id)
            || rescore_decoys.is_some_and(|d| d.contains(&r.spec_id))
    });
    let mut omitted = audit_rescore_join(&rows, m, rescore_decoys)?;
    omitted.extend(excluded.iter().cloned());
    Ok(omitted)
}

#[allow(clippy::too_many_arguments)]
pub fn write_qpx(
    out_dir: &Path,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    params: &SearchParams,
    search_index: &SearchIndex,
    fragment_tol: &Tolerance,
    run_id: &str,
    primary_ms_run_paths: &[String],
    run_ends: &[usize],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    rescore_decoys: Option<&HashSet<String>>,
    excluded: &HashSet<String>,
    ipf: Option<&HashMap<String, crate::ipf::IpfQpx>>,
) -> std::io::Result<()> {
    if run_ends.len() != primary_ms_run_paths.len()
        || run_ends.windows(2).any(|w| w[0] > w[1])
        || run_ends.last().copied().unwrap_or(0) != spectra.len()
    {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "QPX: run_ends {:?} do not partition {} spectra over {} inputs",
                run_ends,
                spectra.len(),
                primary_ms_run_paths.len()
            ),
        ));
    }
    std::fs::create_dir_all(out_dir)?;

    let competed_away = omitted_psm_rows(spectra, queues, candidates, rescore, rescore_decoys, excluded)?;

    let mut row_offset = 0usize;
    let batches = spectra
        .chunks(1024)
        .zip(queues.chunks(1024))
        .enumerate()
        .map(|(chunk, (spectra, queues))| {
            let first = chunk * 1024;
            let reference_files: Vec<&str> = (first..first + spectra.len())
                .map(|i| primary_ms_run_paths[run_ends.partition_point(|&end| end <= i)].as_str())
                .collect();
            let batch = build_psms_batch(
                spectra,
                queues,
                candidates,
                search_index,
                run_id,
                &reference_files,
                rescore,
                &competed_away,
                row_offset,
                ipf,
            )?;
            row_offset += batch.num_rows();
            Ok(batch)
        });
    write_parquet_batches(
        &out_dir.join("psms.parquet"),
        "psms",
        psms_schema(),
        batches,
    )?;

    let prot_batch = build_proteins_batch(queues, candidates, search_index, run_id)?;
    write_parquet(&out_dir.join("proteins.parquet"), "proteins", prot_batch)?;

    let sp_batch = build_search_params_batch(params, fragment_tol, run_id, primary_ms_run_paths)?;
    write_parquet(
        &out_dir.join("search_params.parquet"),
        "search_params",
        sp_batch,
    )?;

    let pg_batch = build_protein_groups_batch();
    write_parquet(
        &out_dir.join("protein_groups.parquet"),
        "protein_groups",
        pg_batch,
    )?;

    Ok(())
}

fn write_parquet(path: &Path, file_type: &str, batch: RecordBatch) -> std::io::Result<()> {
    write_parquet_batches(path, file_type, batch.schema(), std::iter::once(Ok(batch)))
}

fn write_parquet_batches(
    path: &Path,
    file_type: &str,
    schema: Arc<Schema>,
    batches: impl IntoIterator<Item = std::io::Result<RecordBatch>>,
) -> std::io::Result<()> {
    let schema_with_meta = Arc::new(Schema::new_with_metadata(
        schema.fields().clone(),
        qpx_metadata(file_type),
    ));
    let props = WriterProperties::builder().build();
    crate::atomic::write_atomic(path, |path| {
        let file = std::fs::File::create(path)?;
        let mut writer =
            ArrowWriter::try_new(file, schema_with_meta.clone(), Some(props)).map_err(to_io)?;
        for batch in batches {
            let batch = batch?;
            let batch = RecordBatch::try_new(schema_with_meta.clone(), batch.columns().to_vec())
                .map_err(to_io)?;
            writer.write(&batch).map_err(to_io)?;
            writer.flush().map_err(to_io)?;
        }
        writer.close().map_err(to_io)?;
        Ok(())
    })
}

fn qpx_metadata(file_type: &str) -> std::collections::HashMap<String, String> {
    let now = chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string();
    let mut m = std::collections::HashMap::new();
    m.insert("qpx_version".to_string(), QPX_VERSION.to_string());
    m.insert("creator".to_string(), CREATOR.to_string());
    m.insert("software_provider".to_string(), CREATOR.to_string());
    m.insert("file_type".to_string(), file_type.to_string());
    m.insert("creation_date".to_string(), now);
    m.insert("uuid".to_string(), uuid_v4());
    m
}

fn to_io<E: std::fmt::Display>(e: E) -> std::io::Error {
    std::io::Error::other(e.to_string())
}

pub(crate) fn uuid_v4() -> String {
    use std::hash::{BuildHasher, Hasher};
    let mk = || {
        let mut h = std::collections::hash_map::RandomState::new().build_hasher();
        h.write_u128(
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0),
        );
        h.write_usize(&h as *const _ as usize);
        h.finish()
    };
    let mut bytes = [0u8; 16];
    bytes[..8].copy_from_slice(&mk().to_le_bytes());
    bytes[8..].copy_from_slice(&mk().to_le_bytes());
    bytes[6] = (bytes[6] & 0x0f) | 0x40;
    bytes[8] = (bytes[8] & 0x3f) | 0x80;
    format!(
        "{:02x}{:02x}{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}",
        bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        bytes[8], bytes[9], bytes[10], bytes[11], bytes[12], bytes[13], bytes[14], bytes[15],
    )
}

fn modification_fields() -> Fields {
    vec![
        Field::new("name", DataType::Utf8, true),
        Field::new("accession", DataType::Utf8, true),
        Field::new(
            "positions",
            DataType::List(Arc::new(Field::new(
                "element",
                DataType::Struct(crate::feature_parquet::position_fields()),
                true,
            ))),
            true,
        ),
    ]
    .into()
}

fn additional_score_fields() -> Fields {
    vec![
        Field::new("score_name", DataType::Utf8, true),
        Field::new("score_value", DataType::Float64, true),
        Field::new("higher_better", DataType::Boolean, true),
    ]
    .into()
}

fn protein_accession_fields() -> Fields {
    vec![
        Field::new("accession", DataType::Utf8, false),
        Field::new("aa_before", DataType::Utf8, true),
        Field::new("aa_after", DataType::Utf8, true),
        Field::new("start", DataType::Int32, true),
        Field::new("end", DataType::Int32, true),
    ]
    .into()
}

fn metavalue_fields() -> Fields {
    vec![
        Field::new("name", DataType::Utf8, true),
        Field::new("value", DataType::Utf8, true),
        Field::new("value_type", DataType::Utf8, true),
    ]
    .into()
}

fn list_of(item: DataType) -> DataType {
    DataType::List(Arc::new(Field::new("element", item, true)))
}

fn psms_schema() -> Arc<Schema> {
    let fields = vec![
        Field::new("sequence", DataType::Utf8, true),
        Field::new("peptidoform", DataType::Utf8, true),
        Field::new(
            "modifications",
            list_of(DataType::Struct(modification_fields())),
            true,
        ),
        Field::new("precursor_charge", DataType::Int32, true),
        Field::new("posterior_error_probability", DataType::Float64, true),
        Field::new("is_decoy", DataType::Boolean, true),
        Field::new("calculated_mz", DataType::Float64, true),
        Field::new("observed_mz", DataType::Float64, true),
        Field::new(
            "additional_scores",
            list_of(DataType::Struct(additional_score_fields())),
            true,
        ),
        Field::new(
            "protein_accessions",
            list_of(DataType::Struct(protein_accession_fields())),
            true,
        ),
        Field::new("predicted_rt", DataType::Float64, true),
        Field::new("reference_file_name", DataType::Utf8, true),
        Field::new("cv_params", DataType::Utf8, true),
        Field::new("scan", DataType::Int32, true),
        Field::new("rt", DataType::Float64, true),
        Field::new("ion_mobility", DataType::Float64, true),
        Field::new("spectrum_reference", DataType::Utf8, true),
        Field::new("score", DataType::Float64, true),
        Field::new("score_type", DataType::Utf8, true),
        Field::new("higher_score_better", DataType::Boolean, true),
        Field::new("hit_index", DataType::Int32, true),
        Field::new("peptide_identification_index", DataType::Int32, true),
        Field::new(
            "psm_metavalues",
            list_of(DataType::Struct(metavalue_fields())),
            true,
        ),
        Field::new(
            "spectrum_metavalues",
            list_of(DataType::Struct(metavalue_fields())),
            true,
        ),
        Field::new("run_identifier", DataType::Utf8, true),
        Field::new("mz_array", list_of(DataType::Float32), true),
        Field::new("intensity_array", list_of(DataType::Float32), true),
        Field::new("charge_array", list_of(DataType::Int32), true),
        Field::new("ion_type_array", list_of(DataType::Utf8), true),
    ];
    Arc::new(Schema::new(fields))
}

struct RescoreRow {
    spec_id: String,
    is_decoy: bool,
    group: (String, i32, String),
}

fn audit_rescore_join(
    rows: &[RescoreRow],
    rescore: &HashMap<String, PercolatorPsm>,
    decoys: Option<&HashSet<String>>,
) -> std::io::Result<HashSet<String>> {
    let present = |r: &RescoreRow| {
        rescore.contains_key(&r.spec_id) || decoys.is_some_and(|d| d.contains(&r.spec_id))
    };
    let mut winners: HashMap<&(String, i32, String), usize> = HashMap::new();
    if decoys.is_some() {
        for r in rows.iter().filter(|r| present(r)) {
            *winners.entry(&r.group).or_insert(0) += 1;
        }
    }
    let mut targets = 0usize;
    let mut missing = 0usize;
    let mut consumed: HashSet<&str> = HashSet::new();
    let mut competed_away: HashSet<String> = HashSet::new();
    for r in rows.iter().filter(|r| !r.is_decoy) {
        targets += 1;
        if rescore.contains_key(&r.spec_id) {
            consumed.insert(r.spec_id.as_str());
        } else if winners.get(&r.group) == Some(&1) {
            competed_away.insert(r.spec_id.clone());
        } else {
            missing += 1;
        }
    }
    if let Some(d) = decoys {
        for r in rows
            .iter()
            .filter(|r| r.is_decoy && !d.contains(&r.spec_id))
        {
            if winners.get(&r.group) == Some(&1) {
                competed_away.insert(r.spec_id.clone());
            }
        }
    }
    if missing > 0 {
        let extra = rescore.len().saturating_sub(consumed.len());
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Percolator (--rescore) join incomplete: {missing}/{targets} non-decoy PSM \
                 rows had no matching Percolator result ({extra} Percolator keys unused). \
                 This corrupts FDR — the idparquet SpecId derivation has drifted from the \
                 PIN written for Percolator. Refusing to emit silently-null q-values."
            ),
        ));
    }
    Ok(competed_away)
}

fn rescore_rows(
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
) -> Vec<RescoreRow> {
    let mut rows = Vec::new();
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        let sorted = queue.clone().into_rank_sorted_vec();
        let emitted = emitted_rows(spec, &sorted, candidates);
        let ranks: Vec<u32> = emitted.iter().map(|(r, _)| *r).collect();
        let psms: Vec<&PsmMatch> = emitted.iter().map(|(_, p)| &**p).collect();
        let multi_row = psms.len() > 1;
        let spec_id_str = if spec.title.is_empty() {
            format!("scan={}", spec.scan.unwrap_or(0))
        } else {
            spec.title.clone()
        };
        let scan_num = spec.scan.unwrap_or(0);
        for (hit, psm) in psms.iter().enumerate() {
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            let precursor_mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
            let spec_id = format_spec_id(&spec_id_str, scan_num, ranks[hit], hit, multi_row);
            let run = spec_id_run(&spec_id).to_string();
            rows.push(RescoreRow {
                spec_id,
                is_decoy: cand.is_decoy,
                group: (run, scan_num, exp_mass_text(precursor_mz, psm.charge_used)),
            });
        }
    }
    rows
}

#[allow(clippy::too_many_arguments)]
fn build_psms_batch(
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    search_index: &SearchIndex,
    run_id: &str,
    reference_files: &[&str],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    competed_away: &HashSet<String>,
    row_offset: usize,
    ipf: Option<&HashMap<String, crate::ipf::IpfQpx>>,
) -> std::io::Result<RecordBatch> {
    let schema = psms_schema();

    let mut sequence = StringBuilder::new();
    let mut peptidoform = StringBuilder::new();
    let mut modifications = list_struct_builder(modification_fields());
    let mut precursor_charge = Int32Builder::new();
    let mut pep = Float64Builder::new();
    let mut is_decoy = BooleanBuilder::new();
    let mut calculated_mz = Float64Builder::new();
    let mut observed_mz = Float64Builder::new();
    let mut additional_scores = list_struct_builder(additional_score_fields());
    let mut protein_accessions = list_struct_builder(protein_accession_fields());
    let mut predicted_rt = Float64Builder::new();
    let mut reference_file_name = StringBuilder::new();
    let mut cv_params = StringBuilder::new();
    let mut scan = Int32Builder::new();
    let mut rt = Float64Builder::new();
    let mut ion_mobility = Float64Builder::new();
    let mut spectrum_reference = StringBuilder::new();
    let mut score = Float64Builder::new();
    let mut score_type = StringBuilder::new();
    let mut higher_score_better = BooleanBuilder::new();
    let mut hit_index = Int32Builder::new();
    let mut peptide_identification_index = Int32Builder::new();
    let mut psm_metavalues = list_struct_builder(metavalue_fields());
    let mut spectrum_metavalues = list_struct_builder(metavalue_fields());
    let mut run_identifier = StringBuilder::new();
    let mut mz_array =
        primitive_list_builder(arrow::array::Float32Builder::new(), DataType::Float32);
    let mut intensity_array =
        primitive_list_builder(arrow::array::Float32Builder::new(), DataType::Float32);
    let mut charge_array = primitive_list_builder(Int32Builder::new(), DataType::Int32);
    let mut ion_type_array = primitive_list_builder(StringBuilder::new(), DataType::Utf8);

    let mut pep_id_index = row_offset;

    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        let sorted = queue.clone().into_rank_sorted_vec();
        let rows = emitted_rows(spec, &sorted, candidates);
        let ranks: Vec<u32> = rows.iter().map(|(r, _)| *r).collect();
        let psms: Vec<&PsmMatch> = rows.iter().map(|(_, p)| &**p).collect();
        let multi_row = psms.len() > 1;
        let spec_id_str = if spec.title.is_empty() {
            format!("scan={}", spec.scan.unwrap_or(0))
        } else {
            spec.title.clone()
        };
        for (hit, &psm) in psms.iter().enumerate() {
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            let z = psm.charge_used as f64;
            let precursor_mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
            let calc_mz = cand.peptide.mass() / z + PROTON;
            let adj_observed_mz = precursor_mz - ISOTOPE * (psm.isotope_offset as f64) / z;
            let scan_num = spec.scan.unwrap_or(0);

            let spec_id = format_spec_id(&spec_id_str, scan_num, ranks[hit], hit, multi_row);
            if competed_away.contains(&spec_id) {
                continue;
            }
            let matched = rescore.and_then(|m| m.get(&spec_id));

            sequence.append_value(bare_sequence(cand));
            peptidoform.append_value(peptidoform_string(cand));
            let localized = ipf.and_then(|m| {
                m.get(&format_spec_id(&spec_id_str, scan_num, ranks[hit], hit, multi_row))
            });
            append_modifications_scored(
                &mut modifications,
                cand,
                localized.map(|l| &l.site_probs),
            );
            precursor_charge.append_value(psm.charge_used as i32);
            match matched {
                Some(p) => pep.append_value(p.pep),
                None => pep.append_null(),
            }
            is_decoy.append_value(cand.is_decoy);
            calculated_mz.append_value(calc_mz);
            observed_mz.append_value(adj_observed_mz);
            append_additional_scores(
                &mut additional_scores,
                psm,
                matched.map(|p| p.q_value),
                localized,
            );
            append_protein_accessions(&mut protein_accessions, psm, candidates, search_index);
            if psm.features.predicted_rt_min != 0.0 {
                predicted_rt.append_value(psm.features.predicted_rt_min as f64);
            } else {
                predicted_rt.append_null();
            }
            reference_file_name.append_value(reference_files[spec_idx]);
            cv_params.append_null();
            scan.append_value(scan_num);
            match spec.rt_seconds {
                Some(v) => rt.append_value(v),
                None => rt.append_null(),
            }
            ion_mobility.append_null();
            spectrum_reference.append_value(spectrum_ref(
                spec,
                scan_num,
                reference_files[spec_idx],
            ));
            score.append_value(psm.score as f64);
            score_type.append_value("andes:RawScore");
            higher_score_better.append_value(true);
            hit_index.append_value(hit as i32);
            peptide_identification_index.append_value(i32::try_from(pep_id_index).map_err(to_io)?);
            {
                let mvb = psm_metavalues.values();
                for (name, value, fmt) in psm_feature_values(psm, ranks[hit]) {
                    let value_str = format_feature_value(value, fmt);
                    let prefixed = format!("andes:{name}");
                    append_metavalue(mvb, &prefixed, &value_str, feature_value_type(fmt));
                }
                if let Some(l) = localized {
                    append_metavalue(mvb, "andes:ipf_best_sites", &l.best_sites, "string");
                    append_metavalue(
                        mvb,
                        "andes:ipf_site_probabilities",
                        &l.site_probs_label,
                        "string",
                    );
                    append_metavalue(
                        mvb,
                        "andes:ipf_truncated",
                        if l.truncated { "1" } else { "0" },
                        "int",
                    );
                }
                psm_metavalues.append(true);
            }
            {
                let sb = spectrum_metavalues.values();
                append_metavalue(sb, "MS:1001115", &scan_num.to_string(), "string");
                spectrum_metavalues.append(true);
            }
            run_identifier.append_value(run_id);

            let has_peaks = !spec.peaks.is_empty();
            for &(mz, inten) in &spec.peaks {
                mz_array.values().append_value(mz as f32);
                intensity_array.values().append_value(inten);
            }
            mz_array.append(has_peaks);
            intensity_array.append(has_peaks);
            charge_array.append(false);
            ion_type_array.append(false);

            pep_id_index += 1;
        }
    }

    let columns: Vec<ArrayRef> = vec![
        Arc::new(sequence.finish()),
        Arc::new(peptidoform.finish()),
        Arc::new(modifications.finish()),
        Arc::new(precursor_charge.finish()),
        Arc::new(pep.finish()),
        Arc::new(is_decoy.finish()),
        Arc::new(calculated_mz.finish()),
        Arc::new(observed_mz.finish()),
        Arc::new(additional_scores.finish()),
        Arc::new(protein_accessions.finish()),
        Arc::new(predicted_rt.finish()),
        Arc::new(reference_file_name.finish()),
        Arc::new(cv_params.finish()),
        Arc::new(scan.finish()),
        Arc::new(rt.finish()),
        Arc::new(ion_mobility.finish()),
        Arc::new(spectrum_reference.finish()),
        Arc::new(score.finish()),
        Arc::new(score_type.finish()),
        Arc::new(higher_score_better.finish()),
        Arc::new(hit_index.finish()),
        Arc::new(peptide_identification_index.finish()),
        Arc::new(psm_metavalues.finish()),
        Arc::new(spectrum_metavalues.finish()),
        Arc::new(run_identifier.finish()),
        Arc::new(mz_array.finish()),
        Arc::new(intensity_array.finish()),
        Arc::new(charge_array.finish()),
        Arc::new(ion_type_array.finish()),
    ];
    RecordBatch::try_new(schema, columns).map_err(to_io)
}

pub fn bare_sequence(cand: &Candidate) -> String {
    cand.peptide
        .residues
        .iter()
        .map(|aa| aa.residue as char)
        .collect()
}

/// ProForma 2.0 text of a candidate: a residue's modification follows it in brackets, an
/// N-terminal one precedes the sequence as `[mod]-` and a C-terminal one follows it as `-[mod]`.
/// The tag is the modification's accession, or its signed mass delta when it has none.
pub fn peptidoform_string(cand: &Candidate) -> String {
    use model::modification::ModLocation::{CTerm, NTerm, ProtCTerm, ProtNTerm};
    let residues = &cand.peptide.residues;
    let last = residues.len().saturating_sub(1);
    let mut n_term = String::new();
    let mut body = String::new();
    let mut c_term = String::new();
    for (i, aa) in residues.iter().enumerate() {
        body.push(aa.residue as char);
        for m in aa.mod_.iter().flat_map(|m| m.components()) {
            let out = match m.location {
                NTerm | ProtNTerm if i == 0 => &mut n_term,
                CTerm | ProtCTerm if i == last => &mut c_term,
                _ => &mut body,
            };
            push_mod_tag(out, m);
        }
    }
    let mut s = n_term;
    if !s.is_empty() {
        s.push('-');
    }
    s.push_str(&body);
    if !c_term.is_empty() {
        s.push('-');
        s.push_str(&c_term);
    }
    s
}

/// Reads a ProForma peptidoform as written by `peptidoform_string`: `[mod]-` before the sequence is
/// N-terminal (position 0) and `-[mod]` after it C-terminal (position length + 1).
pub fn parse_peptidoform(peptidoform: &str) -> (String, Vec<ModRecord>) {
    let (mut seq, mut mods, mut c_term) = (String::new(), Vec::new(), false);
    let mut chars = peptidoform.chars();
    while let Some(c) = chars.next() {
        match c {
            '-' => c_term = !seq.is_empty(),
            '[' => {
                let tag: String = chars.by_ref().take_while(|&c| c != ']').collect();
                let n = seq.chars().count() as i32;
                let (position, amino_acid) = match (n, c_term) {
                    (0, _) => (0, None),
                    (_, true) => (n + 1, None),
                    _ => (n, seq.chars().last()),
                };
                mods.push(ModRecord {
                    accession: tag.starts_with("UNIMOD:").then(|| tag.clone()),
                    name: tag,
                    position,
                    amino_acid,
                });
            }
            _ => seq.push(c),
        }
    }
    (seq, mods)
}

fn push_mod_tag(s: &mut String, m: &model::modification::Modification) {
    match m.accession.as_deref() {
        Some(acc) if !acc.is_empty() => {
            s.push('[');
            s.push_str(acc);
            s.push(']');
        }
        _ => s.push_str(&format!("[{:+.5}]", m.mass_delta)),
    }
}

#[cfg(test)]
fn append_modifications(builder: &mut ListBuilder<StructBuilder>, cand: &Candidate) {
    append_modifications_scored(builder, cand, None)
}

fn append_modifications_scored(
    builder: &mut ListBuilder<StructBuilder>,
    cand: &Candidate,
    site_probs: Option<&HashMap<(i32, String), f64>>,
) {
    let sb = builder.values();
    let len = cand.peptide.residues.len();
    for (i, aa) in cand.peptide.residues.iter().enumerate() {
        let Some(m) = ModRecord::of_residue(i, len, aa) else {
            continue;
        };
        sb.field_builder::<StringBuilder>(0)
            .unwrap()
            .append_value(&m.name);
        sb.field_builder::<StringBuilder>(1)
            .unwrap()
            .append_option(m.accession.as_deref());
        let pos_list = sb.field_builder::<ListBuilder<StructBuilder>>(2).unwrap();
        let psb = pos_list.values();
        psb.field_builder::<Int32Builder>(0)
            .unwrap()
            .append_value(m.position);
        psb.field_builder::<StringBuilder>(1)
            .unwrap()
            .append_option(m.amino_acid.map(|c| c.to_string()));
        let scores = psb.field_builder::<ListBuilder<StructBuilder>>(2).unwrap();
        match site_probs.and_then(|s| s.get(&(m.position, m.name.clone()))) {
            Some(&p) => {
                let s = scores.values();
                s.field_builder::<StringBuilder>(0)
                    .unwrap()
                    .append_value("ipf_site_probability");
                s.field_builder::<Float64Builder>(1)
                    .unwrap()
                    .append_value(p);
                s.field_builder::<BooleanBuilder>(2)
                    .unwrap()
                    .append_value(true);
                s.append(true);
                scores.append(true);
            }
            None => scores.append_null(),
        }
        psb.append(true);
        pos_list.append(true);
        sb.append(true);
    }
    builder.append(true);
}

fn append_additional_scores(
    builder: &mut ListBuilder<StructBuilder>,
    psm: &search::psm::PsmMatch,
    q_value: Option<f64>,
    ipf: Option<&crate::ipf::IpfQpx>,
) {
    let f = &psm.features;
    let entries: [(&str, f64); 9] = [
        ("RankScore", psm.rank_score as f64),
        ("RawScore", f.strong_score as f64),
        ("NumMatchedMainIons", f.num_matched_main_ions as f64),
        ("TailorScore", f.tailor_score as f64),
        ("DeltaRankScore", f.delta_raw_score as f64),
        ("EdgeScore", f.edge_score as f64),
        (
            "ExplainedIonCurrentRatio",
            f.explained_ion_current_ratio as f64,
        ),
        ("MeanErrorTop7", f.mean_error_top7 as f64),
        ("RichIonLLR", f.rich_ion_llr as f64),
    ];
    let sb = builder.values();
    for (name, value) in entries {
        sb.field_builder::<StringBuilder>(0)
            .unwrap()
            .append_value(name);
        sb.field_builder::<Float64Builder>(1)
            .unwrap()
            .append_value(value);
        sb.field_builder::<BooleanBuilder>(2)
            .unwrap()
            .append_value(true);
        sb.append(true);
    }
    if let Some(q) = q_value {
        sb.field_builder::<StringBuilder>(0)
            .unwrap()
            .append_value("q-value");
        sb.field_builder::<Float64Builder>(1)
            .unwrap()
            .append_value(q);
        sb.field_builder::<BooleanBuilder>(2)
            .unwrap()
            .append_value(false);
        sb.append(true);
    }
    if let Some(l) = ipf {
        for (name, value, higher) in [
            ("ipf_best_prob", l.best_prob, true),
            ("ipf_pep", l.pep, false),
            ("ipf_q", l.q, false),
        ] {
            sb.field_builder::<StringBuilder>(0)
                .unwrap()
                .append_value(name);
            sb.field_builder::<Float64Builder>(1)
                .unwrap()
                .append_value(value);
            sb.field_builder::<BooleanBuilder>(2)
                .unwrap()
                .append_value(higher);
            sb.append(true);
        }
    }
    builder.append(true);
}

fn append_protein_accessions(
    builder: &mut ListBuilder<StructBuilder>,
    psm: &search::psm::PsmMatch,
    candidates: &[Candidate],
    search_index: &SearchIndex,
) {
    let sb = builder.values();
    for &cidx in &psm.candidate_idxs {
        let cand = &candidates[cidx as usize];
        let accession = resolve_accession(cand, search_index);
        let aa_before = (cand.peptide.pre as char).to_string();
        let aa_after = (cand.peptide.post as char).to_string();
        let start = cand.start_offset_in_protein as i32;
        let end = start + cand.peptide.length() as i32;
        sb.field_builder::<StringBuilder>(0)
            .unwrap()
            .append_value(&accession);
        sb.field_builder::<StringBuilder>(1)
            .unwrap()
            .append_value(&aa_before);
        sb.field_builder::<StringBuilder>(2)
            .unwrap()
            .append_value(&aa_after);
        sb.field_builder::<Int32Builder>(3)
            .unwrap()
            .append_value(start);
        sb.field_builder::<Int32Builder>(4)
            .unwrap()
            .append_value(end);
        sb.append(true);
    }
    builder.append(true);
}

fn spectrum_ref(spec: &Spectrum, scan: i32, reference_file: &str) -> String {
    if spec.title.is_empty() {
        return format!("scan={scan}");
    }
    let stem = std::path::Path::new(reference_file)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    spec.title
        .strip_prefix(stem)
        .and_then(|rest| rest.strip_prefix('/'))
        .unwrap_or(&spec.title)
        .to_string()
}

fn append_metavalue(sb: &mut StructBuilder, name: &str, value: &str, value_type: &str) {
    sb.field_builder::<StringBuilder>(0)
        .unwrap()
        .append_value(name);
    sb.field_builder::<StringBuilder>(1)
        .unwrap()
        .append_value(value);
    sb.field_builder::<StringBuilder>(2)
        .unwrap()
        .append_value(value_type);
    sb.append(true);
}

fn list_struct_builder(fields: Fields) -> ListBuilder<StructBuilder> {
    let element_field = Arc::new(Field::new(
        "element",
        DataType::Struct(fields.clone()),
        true,
    ));
    let child_builders: Vec<Box<dyn arrow::array::ArrayBuilder>> = fields
        .iter()
        .map(|f| make_builder_for(f.data_type()))
        .collect();
    ListBuilder::new(StructBuilder::new(fields, child_builders)).with_field(element_field)
}

fn primitive_list_builder<T: arrow::array::ArrayBuilder>(
    values_builder: T,
    inner: DataType,
) -> ListBuilder<T> {
    let element_field = Arc::new(Field::new("element", inner, true));
    ListBuilder::new(values_builder).with_field(element_field)
}

fn make_builder_for(dt: &DataType) -> Box<dyn arrow::array::ArrayBuilder> {
    match dt {
        DataType::Utf8 => Box::new(StringBuilder::new()),
        DataType::Float64 => Box::new(Float64Builder::new()),
        DataType::Int32 => Box::new(Int32Builder::new()),
        DataType::Boolean => Box::new(BooleanBuilder::new()),
        DataType::List(item) => match item.data_type() {
            DataType::Struct(inner) => Box::new(list_struct_builder(inner.clone())),
            other => panic!("unsupported QPX nested list item type: {other:?}"),
        },
        other => panic!("unsupported QPX struct field type: {other:?}"),
    }
}

fn build_proteins_batch(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    search_index: &SearchIndex,
    run_id: &str,
) -> std::io::Result<RecordBatch> {
    let fields = vec![
        Field::new("accession", DataType::Utf8, false),
        Field::new("score", DataType::Float64, false),
        Field::new("rank", DataType::Int32, false),
        Field::new("coverage", DataType::Float64, true),
        Field::new("sequence", DataType::Utf8, true),
        Field::new("description", DataType::Utf8, true),
        Field::new("is_decoy", DataType::Boolean, true),
        Field::new("run_identifier", DataType::Utf8, false),
    ];
    let schema = Arc::new(Schema::new(fields));

    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut accession = StringBuilder::new();
    let mut score = Float64Builder::new();
    let mut rank = Int32Builder::new();
    let mut coverage = Float64Builder::new();
    let mut sequence = StringBuilder::new();
    let mut description = StringBuilder::new();
    let mut is_decoy = BooleanBuilder::new();
    let mut run_identifier = StringBuilder::new();

    for queue in queues {
        if queue.is_empty() {
            continue;
        }
        for psm in queue.clone().into_rank_sorted_vec() {
            for &cidx in &psm.candidate_idxs {
                let cand = &candidates[cidx as usize];
                let acc = resolve_accession(cand, search_index);
                if !seen.insert(acc.clone()) {
                    continue;
                }
                accession.append_value(&acc);
                score.append_value(0.0);
                rank.append_value(0);
                coverage.append_null();
                sequence.append_null();
                description.append_null();
                is_decoy.append_value(cand.is_decoy);
                run_identifier.append_value(run_id);
            }
        }
    }

    let columns: Vec<ArrayRef> = vec![
        Arc::new(accession.finish()),
        Arc::new(score.finish()),
        Arc::new(rank.finish()),
        Arc::new(coverage.finish()),
        Arc::new(sequence.finish()),
        Arc::new(description.finish()),
        Arc::new(is_decoy.finish()),
        Arc::new(run_identifier.finish()),
    ];
    RecordBatch::try_new(schema, columns).map_err(to_io)
}

fn build_protein_groups_batch() -> RecordBatch {
    let data_struct = |inner: DataType| -> DataType {
        DataType::Struct(Fields::from(vec![
            Field::new("name", DataType::Utf8, true),
            Field::new("values", list_of(inner), true),
        ]))
    };
    let fields = vec![
        Field::new("group_type", DataType::Utf8, true),
        Field::new("probability", DataType::Float64, true),
        Field::new("accessions", list_of(DataType::Utf8), true),
        Field::new("run_identifier", DataType::Utf8, true),
        Field::new("group_index", DataType::Int32, true),
        Field::new("float_data", list_of(data_struct(DataType::Float64)), true),
        Field::new("string_data", list_of(data_struct(DataType::Utf8)), true),
        Field::new("integer_data", list_of(data_struct(DataType::Int64)), true),
    ];
    RecordBatch::new_empty(Arc::new(Schema::new(fields)))
}

fn build_search_params_batch(
    params: &SearchParams,
    fragment_tol: &Tolerance,
    run_id: &str,
    primary_ms_run_paths: &[String],
) -> std::io::Result<RecordBatch> {
    let fields = vec![
        Field::new("run_identifier", DataType::Utf8, false),
        Field::new("search_engine", DataType::Utf8, false),
        Field::new("search_engine_version", DataType::Utf8, true),
        Field::new("score_type", DataType::Utf8, false),
        Field::new("higher_score_better", DataType::Boolean, false),
        Field::new("mass_type", DataType::Utf8, false),
        Field::new("precursor_mass_tolerance", DataType::Float64, false),
        Field::new("precursor_mass_tolerance_ppm", DataType::Boolean, false),
        Field::new("fragment_mass_tolerance", DataType::Float64, false),
        Field::new("fragment_mass_tolerance_ppm", DataType::Boolean, false),
        Field::new("digestion_enzyme", DataType::Utf8, true),
        Field::new("enzyme_term_specificity", DataType::Utf8, true),
        Field::new("missed_cleavages", DataType::Int32, false),
        Field::new("fixed_modifications", list_of(DataType::Utf8), false),
        Field::new("variable_modifications", list_of(DataType::Utf8), false),
        Field::new("primary_ms_run_paths", list_of(DataType::Utf8), false),
    ];
    let schema = Arc::new(Schema::new(fields));

    let (prec_tol, prec_ppm) = precursor_tol_scalar(&params.precursor_tolerance);
    let (frag_tol, frag_ppm) = match fragment_tol {
        Tolerance::Ppm(v) => (*v, true),
        Tolerance::Da(v) => (*v, false),
    };
    let term_spec = match params.num_tolerable_termini {
        2 => "FULL",
        1 => "SEMI",
        _ => "NONE",
    };

    let mut fixed_mods: Vec<String> = Vec::new();
    let mut variable_mods: Vec<String> = Vec::new();
    for m in params.aa_set.distinct_mods() {
        let residue = match m.residue {
            model::modification::ResidueSpec::Specific(r) => (r as char).to_string(),
            model::modification::ResidueSpec::Wildcard => "*".to_string(),
        };
        let label = format!("{} ({})", m.name, residue);
        if m.fixed {
            fixed_mods.push(label);
        } else {
            variable_mods.push(label);
        }
    }
    fixed_mods.sort();
    fixed_mods.dedup();
    variable_mods.sort();
    variable_mods.dedup();

    let columns: Vec<ArrayRef> = vec![
        Arc::new(arrow::array::StringArray::from(vec![run_id])),
        Arc::new(arrow::array::StringArray::from(vec!["andes"])),
        Arc::new(arrow::array::StringArray::from(vec![env!(
            "CARGO_PKG_VERSION"
        )])),
        Arc::new(arrow::array::StringArray::from(vec!["andes:RawScore"])),
        Arc::new(arrow::array::BooleanArray::from(vec![true])),
        Arc::new(arrow::array::StringArray::from(vec!["MONOISOTOPIC"])),
        Arc::new(arrow::array::Float64Array::from(vec![prec_tol])),
        Arc::new(arrow::array::BooleanArray::from(vec![prec_ppm])),
        Arc::new(arrow::array::Float64Array::from(vec![frag_tol])),
        Arc::new(arrow::array::BooleanArray::from(vec![frag_ppm])),
        Arc::new(arrow::array::StringArray::from(vec![params.enzyme.name()])),
        Arc::new(arrow::array::StringArray::from(vec![term_spec])),
        Arc::new(arrow::array::Int32Array::from(vec![
            params.max_missed_cleavages as i32,
        ])),
        Arc::new(string_list_singleton(fixed_mods)),
        Arc::new(string_list_singleton(variable_mods)),
        Arc::new(string_list_singleton(primary_ms_run_paths.to_vec())),
    ];
    RecordBatch::try_new(schema, columns).map_err(to_io)
}

fn string_list_singleton(values: Vec<String>) -> arrow::array::ListArray {
    let mut b = primitive_list_builder(StringBuilder::new(), DataType::Utf8);
    for v in values {
        b.values().append_value(v);
    }
    b.append(true);
    b.finish()
}

fn precursor_tol_scalar(t: &PrecursorTolerance) -> (f64, bool) {
    let one = |tol: &Tolerance| match tol {
        Tolerance::Ppm(v) => (*v, true),
        Tolerance::Da(v) => (*v, false),
    };
    let (lo, lo_ppm) = one(&t.left);
    let (hi, hi_ppm) = one(&t.right);
    if hi >= lo {
        (hi, hi_ppm)
    } else {
        (lo, lo_ppm)
    }
}

#[cfg(test)]
mod rescore_join_tests {
    use super::*;

    fn psm(id: &str) -> PercolatorPsm {
        PercolatorPsm {
            pooled_q_value: None,
            psm_id: id.to_string(),
            q_value: 0.0,
            pep: 0.0,
            score: None,
            peptide: "PEPTIDE".to_string(),
            proteins: String::new(),
        }
    }

    fn row(id: &str, is_decoy: bool, scan: i32, exp_mass: &str) -> RescoreRow {
        RescoreRow {
            spec_id: id.to_string(),
            is_decoy,
            group: (String::new(), scan, exp_mass.to_string()),
        }
    }

    #[test]
    fn rescore_join_errors_on_unmatched_target() {
        let mut map = HashMap::new();
        map.insert("scan=1_1_1".to_string(), psm("scan=1_1_1"));
        let emitted = vec![row("scan=2_2_1", false, 2, "1000")];
        let err = audit_rescore_join(&emitted, &map, None).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let err = audit_rescore_join(&emitted, &map, Some(&HashSet::new())).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }

    #[test]
    fn rescore_join_reads_groups_per_input() {
        let in_run = |run: &str, id: &str, is_decoy: bool| RescoreRow {
            spec_id: format!("{run}/{id}"),
            is_decoy,
            group: (run.to_string(), 7, "900.25".to_string()),
        };
        let emitted = vec![
            in_run("a", "s7_7_1", false),
            in_run("a", "s7_7_2", false),
            in_run("b", "s7_7_1", true),
            in_run("b", "s7_7_2", false),
        ];
        let mut map = HashMap::new();
        map.insert("a/s7_7_1".to_string(), psm("a/s7_7_1"));
        let decoys: HashSet<String> = ["b/s7_7_1".to_string()].into();
        let away = audit_rescore_join(&emitted, &map, Some(&decoys)).unwrap();
        let mut away: Vec<&str> = away.iter().map(String::as_str).collect();
        away.sort();
        assert_eq!(away, ["a/s7_7_2", "b/s7_7_2"]);

        let merged: Vec<RescoreRow> = emitted
            .iter()
            .map(|r| RescoreRow {
                spec_id: r.spec_id.clone(),
                is_decoy: r.is_decoy,
                group: (String::new(), 7, "900.25".to_string()),
            })
            .collect();
        assert!(audit_rescore_join(&merged, &map, Some(&decoys)).is_err());
    }

    #[test]
    fn rescore_join_ok_when_targets_matched_and_decoys_exempt() {
        let mut map = HashMap::new();
        map.insert("scan=1_1_1".to_string(), psm("scan=1_1_1"));
        let emitted = vec![
            row("scan=1_1_1", false, 1, "1000"),
            row("scan=9_9_1", true, 9, "1000"),
        ];
        assert!(audit_rescore_join(&emitted, &map, None).unwrap().is_empty());
    }

    #[test]
    fn rescore_join_accepts_targets_removed_by_target_decoy_competition() {
        let emitted = vec![
            row("s1_1_0", false, 1, "1000.5"),
            row("s1_2_1", false, 1, "1000.5"),
            row("s1_1_2", false, 1, "1500.25"),
            row("s2_1_0", true, 2, "900"),
            row("s2_2_1", false, 2, "900"),
            row("s2_3_2", false, 2, "900"),
        ];
        let mut map = HashMap::new();
        for id in ["s1_1_0", "s1_1_2"] {
            map.insert(id.to_string(), psm(id));
        }
        let decoys: HashSet<String> = ["s2_1_0".to_string()].into();
        let away = audit_rescore_join(&emitted, &map, Some(&decoys)).unwrap();
        let mut away: Vec<&str> = away.iter().map(String::as_str).collect();
        away.sort();
        assert_eq!(away, ["s1_2_1", "s2_2_1", "s2_3_2"]);

        assert!(audit_rescore_join(&emitted, &map, None).is_err());

        let lost = vec![
            row("s3_1_0", false, 3, "800"),
            row("s3_2_1", false, 3, "800"),
        ];
        assert!(audit_rescore_join(&lost, &map, Some(&decoys)).is_err());

        let mut both = map.clone();
        both.insert("s1_2_1".to_string(), psm("s1_2_1"));
        let three = vec![
            row("s1_1_0", false, 1, "1000.5"),
            row("s1_2_1", false, 1, "1000.5"),
            row("s1_3_2", false, 1, "1000.5"),
        ];
        assert!(audit_rescore_join(&three, &both, Some(&decoys)).is_err());
    }

    #[test]
    fn a_modification_without_an_accession_is_null() {
        use arrow::array::{Array, ListArray, StructArray};
        use model::amino_acid::AminoAcid;
        use model::modification::{leak_mod, ModLocation, Modification, ResidueSpec};
        let m = leak_mod(Modification {
            name: "Custom".into(),
            mass_delta: 14.0,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        });
        let residues = vec![
            AminoAcid::standard(b'P').unwrap(),
            AminoAcid::standard(b'K').unwrap().with_mod(m),
        ];
        let cand = Candidate {
            peptide: model::peptide::Peptide::new(residues, b'K', b'R'),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        };
        let mut b = list_struct_builder(modification_fields());
        append_modifications(&mut b, &cand);
        b.append(true);
        let list: ListArray = b.finish();
        let structs = list.value(0);
        let structs = structs.as_any().downcast_ref::<StructArray>().unwrap();
        assert_eq!(structs.len(), 1);
        assert!(structs.column(1).is_null(0));
    }

    #[test]
    fn peptidoform_tags_become_modifications() {
        let (seq, mods) = parse_peptidoform("LC[UNIMOD:4]AK[+8.01420]");
        assert_eq!(seq, "LCAK");
        assert_eq!(mods.len(), 2);
        assert_eq!(
            (
                mods[0].position,
                mods[0].amino_acid,
                mods[0].accession.as_deref()
            ),
            (2, Some('C'), Some("UNIMOD:4"))
        );
        assert_eq!((mods[1].position, mods[1].accession.clone()), (4, None));
        let (seq, mods) = parse_peptidoform("[UNIMOD:1]-AM[UNIMOD:35]K-[UNIMOD:2]");
        let at: Vec<_> = mods.iter().map(|m| (m.position, m.amino_acid)).collect();
        assert_eq!(seq, "AMK");
        assert_eq!(at, vec![(0, None), (2, Some('M')), (4, None)]);
    }

    #[test]
    fn terminal_modifications_are_written_as_proforma_terminal_tags() {
        use model::amino_acid::AminoAcid;
        use model::modification::{leak_mod, ModLocation, Modification, ResidueSpec};
        let mk = |acc: &str, mass: f64, residue: ResidueSpec, location: ModLocation| {
            let mut m = Modification::new(acc, mass, residue, location, false);
            m.accession = Some(acc.to_string());
            m
        };
        let acetyl = mk(
            "UNIMOD:1",
            42.010565,
            ResidueSpec::Wildcard,
            ModLocation::ProtNTerm,
        );
        let oxidation = mk(
            "UNIMOD:35",
            15.994915,
            ResidueSpec::Specific(b'M'),
            ModLocation::Anywhere,
        );
        let amidated = mk(
            "UNIMOD:2",
            -0.984016,
            ResidueSpec::Wildcard,
            ModLocation::CTerm,
        );
        let cand = |residues: Vec<AminoAcid>| Candidate {
            peptide: model::peptide::Peptide::new(residues, b'-', b'R'),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: true,
            is_protein_c_term: false,
        };
        let aa = |r: u8| AminoAcid::standard(r).unwrap();

        let n_acetyl = cand(vec![
            aa(b'A').with_mod(leak_mod(acetyl.clone())),
            aa(b'P'),
            aa(b'S'),
            aa(b'K'),
        ]);
        assert_eq!(peptidoform_string(&n_acetyl), "[UNIMOD:1]-APSK");

        let mut stacked = Modification::new(
            "UNIMOD:1+UNIMOD:35",
            58.00548,
            ResidueSpec::Specific(b'M'),
            ModLocation::ProtNTerm,
            false,
        );
        stacked.parts = vec![leak_mod(acetyl), leak_mod(oxidation)];
        let both = cand(vec![aa(b'M').with_mod(leak_mod(stacked)), aa(b'K')]);
        assert_eq!(peptidoform_string(&both), "[UNIMOD:1]-M[UNIMOD:35]K");

        let c_amide = cand(vec![aa(b'P'), aa(b'K').with_mod(leak_mod(amidated))]);
        assert_eq!(peptidoform_string(&c_amide), "PK-[UNIMOD:2]");
    }

    #[test]
    fn modification_positions_follow_the_qpx_layout() {
        use arrow::array::{Array, Int32Array, ListArray, StringArray, StructArray};
        use model::amino_acid::AminoAcid;
        use model::modification::{leak_mod, ModLocation, Modification, ResidueSpec};
        let mk = |name: &str, residue: ResidueSpec, location: ModLocation| {
            leak_mod(Modification {
                name: name.into(),
                mass_delta: 14.0,
                residue,
                location,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            })
        };
        let n_term = mk("Acetyl", ResidueSpec::Wildcard, ModLocation::NTerm);
        let side = mk(
            "Phospho",
            ResidueSpec::Specific(b'S'),
            ModLocation::Anywhere,
        );
        let c_term = mk("Amidated", ResidueSpec::Wildcard, ModLocation::CTerm);
        let residues = vec![
            AminoAcid::standard(b'A').unwrap().with_mod(n_term),
            AminoAcid::standard(b'P').unwrap(),
            AminoAcid::standard(b'S').unwrap().with_mod(side),
            AminoAcid::standard(b'K').unwrap().with_mod(c_term),
        ];
        let cand = Candidate {
            peptide: model::peptide::Peptide::new(residues, b'K', b'R'),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        };
        let mut b = list_struct_builder(modification_fields());
        append_modifications(&mut b, &cand);
        b.append(true);
        let list: ListArray = b.finish();
        let structs = list.value(0);
        let structs = structs.as_any().downcast_ref::<StructArray>().unwrap();
        let mut got = Vec::new();
        for i in 0..structs.len() {
            let positions = structs
                .column(2)
                .as_any()
                .downcast_ref::<ListArray>()
                .unwrap();
            let entry = positions.value(i);
            let entry = entry.as_any().downcast_ref::<StructArray>().unwrap();
            let pos = entry
                .column(0)
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap();
            let aa = entry
                .column(1)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            got.push((
                pos.value(0),
                (!aa.is_null(0)).then(|| aa.value(0).to_string()),
            ));
        }
        assert_eq!(got, vec![(0, None), (3, Some("S".to_string())), (5, None)]);
        assert_eq!(
            peptidoform_string(&cand),
            "[+14.00000]-APS[+14.00000]K-[+14.00000]"
        );

        let mut b = list_struct_builder(modification_fields());
        let probs: HashMap<(i32, String), f64> =
            [((3, "Phospho".to_string()), 0.97)].into_iter().collect();
        append_modifications_scored(&mut b, &cand, Some(&probs));
        b.append(true);
        let list: ListArray = b.finish();
        let structs = list.value(0);
        let structs = structs.as_any().downcast_ref::<StructArray>().unwrap();
        let positions = structs.column(2).as_any().downcast_ref::<ListArray>().unwrap();
        let score_of = |i: usize| {
            let entry = positions.value(i);
            let entry = entry.as_any().downcast_ref::<StructArray>().unwrap();
            let scores = entry.column(2).as_any().downcast_ref::<ListArray>().unwrap();
            if scores.is_null(0) {
                return None;
            }
            let s = scores.value(0);
            let s = s.as_any().downcast_ref::<StructArray>().unwrap();
            let name = s.column(0).as_any().downcast_ref::<StringArray>().unwrap();
            let v = s
                .column(1)
                .as_any()
                .downcast_ref::<arrow::array::Float64Array>()
                .unwrap();
            Some((name.value(0).to_string(), v.value(0)))
        };
        assert_eq!(score_of(0), None);
        assert_eq!(score_of(1), Some(("ipf_site_probability".to_string(), 0.97)));
        assert_eq!(score_of(2), None);
    }
}
