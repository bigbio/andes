//! QPX `.idparquet` output writer (OpenMS QPXFile 1.0 compatible).
//!
//! Writes a parquet PSM bundle DIRECTORY alongside the existing `.pin`/`.tsv`
//! so that OpenMS / quantms tooling can consume andes results natively. The
//! directory mirrors the layout produced by OpenMS `QPXFile::store`:
//!
//! ```text
//! <name>.idparquet/
//!   psms.parquet
//!   proteins.parquet
//!   search_params.parquet
//! ```
//!
//! Each file carries the QPX schema metadata key-value pairs
//! (`qpx_version`, `creator`, `software_provider`, `file_type`,
//! `creation_date`, `uuid`) on its Arrow schema, exactly as OpenMS writes them
//! (verified against a real OpenMS-emitted bundle). `creator`/`software_provider`
//! are `"andes"` (this engine) rather than `"OpenMS"`.
//!
//! # Column semantics (psms.parquet)
//!
//! The schema (column names + Arrow types) is copied verbatim from the OpenMS
//! `QPXPSMSchema` so a downstream reader can use one code path for OpenMS-native and
//! andes output. andes sources every field it can from the same per-PSM data
//! the PIN writer uses ([`crate::pin`]); fields andes does not compute (PEP,
//! predicted RT, ion mobility, cv params, metavalues) are written null/empty —
//! see the module-level field-mapping note in the closing report.
//!
//! Rows are emitted rank-1-first per spectrum (same ordering as the PIN: each
//! queue is `into_rank_sorted_vec()`), one row per retained PSM.

use std::collections::HashMap;
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
use search::psm::TopNQueue;
use search::search_index::SearchIndex;
use search::search_params::SearchParams;

use crate::percolator::PercolatorPsm;
use crate::pin::{
    feature_value_type, format_feature_value, format_spec_id, psm_feature_values,
};
use crate::row_context::{iter_ranked_by_rank_score, resolve_accession};

/// QPX format version string written into every file's schema metadata.
const QPX_VERSION: &str = "1.0";
/// Engine name written into `creator` / `software_provider`.
const CREATOR: &str = "andes";

// ── public API ───────────────────────────────────────────────────────────────

/// Write a QPX `.idparquet/` bundle at `out_dir`.
///
/// `out_dir` is created (recursively) if missing; the three parquet members are
/// written into it. `spectra`/`queues` are parallel slices (`queues[i]` holds
/// the top-N PSMs for `spectra[i]`); `candidates` is the per-search candidate
/// pool; `search_index` resolves protein accessions; `params`/`fragment_tol`
/// populate `search_params.parquet` (`fragment_tol` is the resolved per-model
/// `param.mme`). `run_id` is a stable identifier for this run, written into the
/// `run_identifier` column of every row. `primary_ms_run_paths` are the input
/// spectrum file paths.
///
/// `rescore` is an optional `SpecId → PercolatorPsm` map (from a `--rescore`
/// Percolator run). When `Some`, each row's `posterior_error_probability` is
/// filled from the matched PEP and a `("q-value", q, false)` entry is appended
/// to `additional_scores`. When `None` (rescore off) the bundle is byte-schema-
/// identical to the pre-rescore output (PEP null, no q-value entry). The join key
/// is reconstructed via [`crate::pin::format_spec_id`] so it matches the PIN's
/// `SpecId` exactly.
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
    rescore: Option<&HashMap<String, PercolatorPsm>>,
) -> std::io::Result<()> {
    std::fs::create_dir_all(out_dir)?;

    let reference_file = primary_ms_run_paths.first().cloned().unwrap_or_default();

    let psm_batch = build_psms_batch(
        spectra,
        queues,
        candidates,
        search_index,
        run_id,
        &reference_file,
        rescore,
    )?;
    write_parquet(&out_dir.join("psms.parquet"), "psms", psm_batch)?;

    let prot_batch = build_proteins_batch(queues, candidates, search_index, run_id)?;
    write_parquet(&out_dir.join("proteins.parquet"), "proteins", prot_batch)?;

    let sp_batch = build_search_params_batch(params, fragment_tol, run_id, primary_ms_run_paths)?;
    write_parquet(&out_dir.join("search_params.parquet"), "search_params", sp_batch)?;

    // OpenMS' QPX reader (PSMArrowIO::importFromParquet) requires a
    // protein_groups member to load the bundle, even though andes does not run
    // protein inference (that is the downstream ProteinInference/Percolator
    // step). Emit a valid zero-row file so OpenMS tools can read the idparquet.
    let pg_batch = build_protein_groups_batch();
    write_parquet(&out_dir.join("protein_groups.parquet"), "protein_groups", pg_batch)?;

    Ok(())
}

// ── parquet file writer (metadata + single batch) ─────────────────────────────

/// Write one RecordBatch to `path`, stamping the QPX schema metadata
/// (`file_type` = `file_type`) onto the schema before writing.
fn write_parquet(path: &Path, file_type: &str, batch: RecordBatch) -> std::io::Result<()> {
    // Re-key the batch's schema with the QPX metadata so the written file
    // carries it (OpenMS readers key on these). UUID is per-file.
    let meta = qpx_metadata(file_type);
    let schema_with_meta =
        Arc::new(Schema::new_with_metadata(batch.schema().fields().clone(), meta));
    // Rebuild the batch against the metadata-bearing schema (same columns).
    let batch = RecordBatch::try_new(schema_with_meta.clone(), batch.columns().to_vec())
        .map_err(to_io)?;

    let props = WriterProperties::builder().build();
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema_with_meta, Some(props)).map_err(to_io)?;
    writer.write(&batch).map_err(to_io)?;
    writer.close().map_err(to_io)?;
    Ok(())
}

/// Build the QPX schema-metadata map for a file of the given `file_type`.
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

/// Generate a random UUID v4 string without pulling in the `uuid` crate (which
/// is not in the workspace lockfile). 122 random bits + the version/variant
/// nibbles, sourced from `RandomState`-seeded hashing of high-resolution
/// clocks — sufficient for a per-file identifier (not a security token).
fn uuid_v4() -> String {
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
    // Set version (4) and variant (RFC 4122) nibbles.
    bytes[6] = (bytes[6] & 0x0f) | 0x40;
    bytes[8] = (bytes[8] & 0x3f) | 0x80;
    format!(
        "{:02x}{:02x}{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}",
        bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        bytes[8], bytes[9], bytes[10], bytes[11], bytes[12], bytes[13], bytes[14], bytes[15],
    )
}

// ── shared nested-type field definitions ──────────────────────────────────────

/// `positions` inner list element: `struct<position: string, scores: double>`.
fn position_fields() -> Fields {
    vec![
        Field::new("position", DataType::Utf8, true),
        Field::new("scores", DataType::Float64, true),
    ]
    .into()
}

/// `modifications` element:
/// `struct<name: string, accession: string, positions: list<position>>`.
fn modification_fields() -> Fields {
    vec![
        Field::new("name", DataType::Utf8, true),
        Field::new("accession", DataType::Utf8, true),
        Field::new(
            "positions",
            DataType::List(Arc::new(Field::new(
                "element",
                DataType::Struct(position_fields()),
                true,
            ))),
            true,
        ),
    ]
    .into()
}

/// `additional_scores` element:
/// `struct<score_name: string, score_value: double, higher_better: bool>`.
fn additional_score_fields() -> Fields {
    vec![
        Field::new("score_name", DataType::Utf8, true),
        Field::new("score_value", DataType::Float64, true),
        Field::new("higher_better", DataType::Boolean, true),
    ]
    .into()
}

/// `protein_accessions` element:
/// `struct<accession: string, aa_before: string, aa_after: string, start: int32, end: int32>`.
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

/// `*_metavalues` element: `struct<name: string, value: string, value_type: string>`.
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

/// Full `psms.parquet` Arrow schema (column order + types verbatim from QPX).
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

// ── psms.parquet builder ──────────────────────────────────────────────────────

/// Validate the Percolator (`--rescore`) join. When rescoring is on, every
/// emitted **non-decoy** PSM row must have a matching Percolator result; a
/// missing match means the `SpecId` derivation drifted between the PIN and the
/// idparquet, which silently corrupts downstream FDR (those target rows would
/// otherwise be written with a null q-value). Decoys are exempt — Percolator
/// writes them to a separate file the caller does not load. On failure the
/// error reports how many target rows were unmatched and how many Percolator
/// keys went unused.
fn audit_rescore_join(
    emitted: &[(String, bool)],
    rescore: &HashMap<String, PercolatorPsm>,
) -> std::io::Result<()> {
    let mut targets = 0usize;
    let mut missing = 0usize;
    let mut consumed: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for (spec_id, is_decoy) in emitted {
        if *is_decoy {
            continue;
        }
        targets += 1;
        if rescore.contains_key(spec_id) {
            consumed.insert(spec_id.as_str());
        } else {
            missing += 1;
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
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_psms_batch(
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    search_index: &SearchIndex,
    run_id: &str,
    reference_file: &str,
    rescore: Option<&HashMap<String, PercolatorPsm>>,
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
    let mut mz_array = primitive_list_builder(
        arrow::array::Float32Builder::new(),
        DataType::Float32,
    );
    let mut intensity_array = primitive_list_builder(
        arrow::array::Float32Builder::new(),
        DataType::Float32,
    );
    let mut charge_array = primitive_list_builder(Int32Builder::new(), DataType::Int32);
    let mut ion_type_array = primitive_list_builder(StringBuilder::new(), DataType::Utf8);

    // Per-spectrum identification index (matches OpenMS, which gives each
    // spectrum's hit-list a stable index). Incremented once per emitted PSM.
    let mut pep_id_index: i32 = 0;

    // When rescoring is on, record each emitted row's (SpecId, is_decoy) so the
    // Percolator join can be audited for completeness after the loop. Empty (no
    // allocation) and unused when rescore is off, keeping that path identical.
    let mut rescore_emitted: Vec<(String, bool)> = Vec::new();

    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        let psms = queue.clone().into_rank_sorted_vec();
        // Reproduce the PIN's SpecId derivation EXACTLY so the Percolator join
        // key matches: `rank` from iter_ranked_by_rank_score, `row_idx` = the
        // per-scan emission index (== `hit`), `multi_row` = scan emits >1 row.
        let ranks: Vec<u32> = iter_ranked_by_rank_score(&psms).map(|(r, _)| r).collect();
        let multi_row = psms.len() > 1;
        let spec_id_str = if spec.title.is_empty() {
            format!("scan={}", spec.scan.unwrap_or(0))
        } else {
            spec.title.clone()
        };
        for (hit, psm) in psms.iter().enumerate() {
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            let z = psm.charge_used as f64;
            let precursor_mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
            let calc_mz = cand.peptide.mass() / z + PROTON;
            // Apply the matched isotope offset so observed/calc agree (PIN parity).
            let adj_observed_mz = precursor_mz - ISOTOPE * (psm.isotope_offset as f64) / z;
            let scan_num = spec.scan.unwrap_or(0);

            // Percolator join: reconstruct this row's SpecId and look up the
            // matched target-PSM result (PEP + q-value). `None` when rescore is
            // off, or when a row had no matching Percolator output (e.g. a decoy,
            // which Percolator writes to the decoy file we do not load here).
            let matched = if let Some(m) = rescore {
                let spec_id = format_spec_id(&spec_id_str, scan_num, ranks[hit], hit, multi_row);
                let hit_result = m.get(&spec_id);
                rescore_emitted.push((spec_id, cand.is_decoy));
                hit_result
            } else {
                None
            };

            sequence.append_value(bare_sequence(cand));
            peptidoform.append_value(peptidoform_string(cand));
            append_modifications(&mut modifications, cand);
            precursor_charge.append_value(psm.charge_used as i32);
            match matched {
                Some(p) => pep.append_value(p.pep),
                None => pep.append_null(),
            }
            is_decoy.append_value(cand.is_decoy);
            calculated_mz.append_value(calc_mz);
            observed_mz.append_value(adj_observed_mz);
            append_additional_scores(&mut additional_scores, psm, matched.map(|p| p.q_value));
            append_protein_accessions(&mut protein_accessions, psm, candidates, search_index);
            predicted_rt.append_null();
            reference_file_name.append_value(reference_file);
            cv_params.append_null();
            scan.append_value(scan_num);
            match spec.rt_seconds {
                Some(v) => rt.append_value(v),
                None => rt.append_null(),
            }
            ion_mobility.append_null();
            spectrum_reference.append_value(spectrum_ref(spec, scan_num));
            // Headline score = the rank/RawScore the PIN emits (psm.score, the
            // ranking signal; rounded in the PIN, kept full-precision here).
            score.append_value(psm.score as f64);
            score_type.append_value("andes:RawScore");
            higher_score_better.append_value(true);
            hit_index.append_value(hit as i32);
            peptide_identification_index.append_value(pep_id_index);
            // Per-PSM discriminative features (the column OpenMS'
            // PercolatorAdapter reads). Sourced from the SHARED pin.rs feature
            // vector so the idparquet carries the SAME features as the PIN.
            // `rank` (1-based) gates DeltaRankScore exactly as the PIN does.
            //
            // Every feature name is written with the `andes:` prefix (mirroring
            // comet's `COMET:` convention) so OpenMS' PercolatorAdapter collects
            // them as andes-specific features. The prefix lives ONLY here — the
            // PIN writer's own column names (in `crate::pin`) are unprefixed.
            {
                let mvb = psm_metavalues.values();
                for (name, value, fmt) in psm_feature_values(psm, ranks[hit]) {
                    let value_str = format_feature_value(value, fmt);
                    let prefixed = format!("andes:{name}");
                    append_metavalue(mvb, &prefixed, &value_str, feature_value_type(fmt));
                }
                psm_metavalues.append(true);
            }
            // One spectrum metavalue: the MS:1001115 scan number (as OpenMS does).
            {
                let sb = spectrum_metavalues.values();
                append_metavalue(sb, "MS:1001115", &scan_num.to_string(), "string");
                spectrum_metavalues.append(true);
            }
            run_identifier.append_value(run_id);

            // Spectrum peak arrays (charge/ion-type unknown → null lists).
            for &(mz, inten) in &spec.peaks {
                mz_array.values().append_value(mz as f32);
                intensity_array.values().append_value(inten);
            }
            mz_array.append(true);
            intensity_array.append(true);
            charge_array.append(false);
            ion_type_array.append(false);

            pep_id_index += 1;
        }
    }

    // Fail loud if the Percolator join was incomplete (target rows with no
    // matching result), rather than silently emitting null q-values.
    if let Some(m) = rescore {
        audit_rescore_join(&rescore_emitted, m)?;
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

// ── per-PSM field derivations ─────────────────────────────────────────────────

/// Bare residue sequence (no mods, no flanks).
fn bare_sequence(cand: &Candidate) -> String {
    cand.peptide
        .residues
        .iter()
        .map(|aa| aa.residue as char)
        .collect()
}

/// ProForma-ish peptidoform: residues with `[UNIMOD:NN]` (or `[+delta]`) after
/// each modified residue. Falls back to the accession when known, else the mass
/// delta. No flanks (OpenMS peptidoform convention).
fn peptidoform_string(cand: &Candidate) -> String {
    let mut s = String::new();
    for aa in &cand.peptide.residues {
        s.push(aa.residue as char);
        if let Some(m) = aa.mod_.as_ref() {
            match m.accession.as_deref() {
                Some(acc) if !acc.is_empty() => {
                    s.push('[');
                    s.push_str(acc);
                    s.push(']');
                }
                _ => s.push_str(&format!("[{:+.5}]", m.mass_delta)),
            }
        }
    }
    s
}

/// Append the `modifications` list cell for one candidate: one struct per
/// modified residue with `name`, `accession`, and a single `positions` entry
/// (`position` = `"<residue><1-based index>"`, `scores` = null — andes does not
/// localize).
fn append_modifications(builder: &mut ListBuilder<StructBuilder>, cand: &Candidate) {
    let sb = builder.values();
    for (i, aa) in cand.peptide.residues.iter().enumerate() {
        let Some(m) = aa.mod_.as_ref() else { continue };
        sb.field_builder::<StringBuilder>(0).unwrap().append_value(&m.name);
        sb.field_builder::<StringBuilder>(1)
            .unwrap()
            .append_value(m.accession.as_deref().unwrap_or(""));
        // positions: list<struct<position: string, scores: double>>
        let pos_list = sb.field_builder::<ListBuilder<StructBuilder>>(2).unwrap();
        let psb = pos_list.values();
        let position = format!("{}.{}", aa.residue as char, i + 1);
        psb.field_builder::<StringBuilder>(0).unwrap().append_value(&position);
        psb.field_builder::<Float64Builder>(1).unwrap().append_null();
        psb.append(true);
        pos_list.append(true);
        sb.append(true);
    }
    builder.append(true);
}

/// Append the `additional_scores` list cell: the informative per-PSM andes
/// scores beyond the headline RawScore. All andes scores are higher-is-better.
///
/// When `q_value` is `Some` (a `--rescore` Percolator result matched this row), a
/// trailing `("q-value", q, higher_better=false)` entry is appended. When `None`
/// (rescore off, or no match) no such entry is written, so the default bundle
/// stays byte-identical to the pre-rescore output.
fn append_additional_scores(
    builder: &mut ListBuilder<StructBuilder>,
    psm: &search::psm::PsmMatch,
    q_value: Option<f64>,
) {
    let f = &psm.features;
    let entries: [(&str, f64); 9] = [
        ("RankScore", psm.rank_score as f64),
        ("RawScore", f.strong_score as f64),
        ("NumMatchedMainIons", f.num_matched_main_ions as f64),
        ("TailorScore", f.tailor_score as f64),
        ("DeltaRankScore", f.delta_raw_score as f64),
        ("EdgeScore", f.edge_score as f64),
        ("ExplainedIonCurrentRatio", f.explained_ion_current_ratio as f64),
        ("MeanErrorTop7", f.mean_error_top7 as f64),
        ("RichIonLLR", f.rich_ion_llr as f64),
    ];
    let sb = builder.values();
    for (name, value) in entries {
        sb.field_builder::<StringBuilder>(0).unwrap().append_value(name);
        sb.field_builder::<Float64Builder>(1).unwrap().append_value(value);
        sb.field_builder::<BooleanBuilder>(2).unwrap().append_value(true);
        sb.append(true);
    }
    // Percolator q-value (lower is better) — only when rescore matched this row.
    if let Some(q) = q_value {
        sb.field_builder::<StringBuilder>(0).unwrap().append_value("q-value");
        sb.field_builder::<Float64Builder>(1).unwrap().append_value(q);
        sb.field_builder::<BooleanBuilder>(2).unwrap().append_value(false);
        sb.append(true);
    }
    builder.append(true);
}

/// Append the `protein_accessions` list cell: one struct per protein the
/// peptide maps to. `aa_before`/`aa_after` = the peptide's `pre`/`post` flanks;
/// `start`/`end` = the peptide's 0-based offset in the source protein (from the
/// candidate) and `start + length`.
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
        sb.field_builder::<StringBuilder>(0).unwrap().append_value(&accession);
        sb.field_builder::<StringBuilder>(1).unwrap().append_value(&aa_before);
        sb.field_builder::<StringBuilder>(2).unwrap().append_value(&aa_after);
        sb.field_builder::<Int32Builder>(3).unwrap().append_value(start);
        sb.field_builder::<Int32Builder>(4).unwrap().append_value(end);
        sb.append(true);
    }
    builder.append(true);
}

/// OpenMS-style `spectrum_reference` native id; `spec.title` if non-empty,
/// else a `scan=N` placeholder.
fn spectrum_ref(spec: &Spectrum, scan: i32) -> String {
    if spec.title.is_empty() {
        format!("scan={scan}")
    } else {
        spec.title.clone()
    }
}

fn append_metavalue(sb: &mut StructBuilder, name: &str, value: &str, value_type: &str) {
    sb.field_builder::<StringBuilder>(0).unwrap().append_value(name);
    sb.field_builder::<StringBuilder>(1).unwrap().append_value(value);
    sb.field_builder::<StringBuilder>(2).unwrap().append_value(value_type);
    sb.append(true);
}

/// Build a `ListBuilder<StructBuilder>` for a `list<struct<...>>` column from
/// the struct's field definitions.
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

/// Build a `ListBuilder` for a `list<T>` primitive column whose element field is
/// named `"element"` (QPX/Parquet LIST convention), matching the schema produced
/// by [`list_of`]. `inner` is the leaf `DataType` (nullable, as `list_of` emits).
fn primitive_list_builder<T: arrow::array::ArrayBuilder>(
    values_builder: T,
    inner: DataType,
) -> ListBuilder<T> {
    let element_field = Arc::new(Field::new("element", inner, true));
    ListBuilder::new(values_builder).with_field(element_field)
}

/// Construct an Arrow `ArrayBuilder` for the given field type. Handles the
/// nested `positions` list inside `modifications` (the only list-of-struct
/// nested inside a struct in the QPX PSM schema).
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

// ── proteins.parquet builder ──────────────────────────────────────────────────

/// Minimal `proteins.parquet`: one row per distinct accession seen in the
/// emitted PSMs (rank-1-first ordering of first appearance). `score`/`rank`/
/// `coverage`/`sequence`/`description`/`modifications` are left at OpenMS-
/// compatible defaults (andes does not run protein inference — that is the
/// downstream ProteinInference step).
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

// ── protein_groups.parquet builder ────────────────────────────────────────────

/// Empty `protein_groups.parquet` with the OpenMS QPX protein-group schema.
///
/// andes does not perform protein inference, so there are no groups to emit, but
/// OpenMS' `PSMArrowIO::importFromParquet` requires this bundle member to be
/// present to read the `.idparquet` directory at all (e.g. PercolatorAdapter).
/// We therefore write a schema-correct, zero-row file. Downstream OpenMS
/// ProteinInference / Percolator populate the actual groups.
fn build_protein_groups_batch() -> RecordBatch {
    // struct<name: utf8, values: list<inner>> — the OpenMS metavalue container
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

// ── search_params.parquet builder ─────────────────────────────────────────────

/// Single-row `search_params.parquet` describing the resolved search.
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

    // Partition the declared mods into fixed/variable name lists.
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
        Arc::new(arrow::array::StringArray::from(vec![env!("CARGO_PKG_VERSION")])),
        Arc::new(arrow::array::StringArray::from(vec!["andes:RawScore"])),
        Arc::new(arrow::array::BooleanArray::from(vec![true])),
        Arc::new(arrow::array::StringArray::from(vec!["MONOISOTOPIC"])),
        Arc::new(arrow::array::Float64Array::from(vec![prec_tol])),
        Arc::new(arrow::array::BooleanArray::from(vec![prec_ppm])),
        Arc::new(arrow::array::Float64Array::from(vec![frag_tol])),
        Arc::new(arrow::array::BooleanArray::from(vec![frag_ppm])),
        Arc::new(arrow::array::StringArray::from(vec![params.enzyme.name()])),
        Arc::new(arrow::array::StringArray::from(vec![term_spec])),
        Arc::new(arrow::array::Int32Array::from(vec![params.max_missed_cleavages as i32])),
        Arc::new(string_list_singleton(fixed_mods)),
        Arc::new(string_list_singleton(variable_mods)),
        Arc::new(string_list_singleton(primary_ms_run_paths.to_vec())),
    ];
    RecordBatch::try_new(schema, columns).map_err(to_io)
}

/// Build a single-row `List<Utf8>` array whose one cell holds `values`.
fn string_list_singleton(values: Vec<String>) -> arrow::array::ListArray {
    let mut b = primitive_list_builder(StringBuilder::new(), DataType::Utf8);
    for v in values {
        b.values().append_value(v);
    }
    b.append(true);
    b.finish()
}

/// Collapse a (possibly asymmetric) precursor tolerance to a single
/// `(value, is_ppm)` for the QPX scalar columns. Uses the wider side.
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
            psm_id: id.to_string(),
            q_value: 0.0,
            pep: 0.0,
            peptide: "PEPTIDE".to_string(),
            proteins: String::new(),
        }
    }

    #[test]
    fn rescore_join_errors_on_unmatched_target() {
        let mut map = HashMap::new();
        map.insert("scan=1_1_1".to_string(), psm("scan=1_1_1"));
        // A target (non-decoy) row whose SpecId is absent from Percolator output.
        let emitted = vec![("scan=2_2_1".to_string(), false)];
        let err = audit_rescore_join(&emitted, &map).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }

    #[test]
    fn rescore_join_ok_when_targets_matched_and_decoys_exempt() {
        let mut map = HashMap::new();
        map.insert("scan=1_1_1".to_string(), psm("scan=1_1_1"));
        let emitted = vec![
            ("scan=1_1_1".to_string(), false), // target, matched
            ("scan=9_9_1".to_string(), true),  // decoy, no match is fine
        ];
        assert!(audit_rescore_join(&emitted, &map).is_ok());
    }
}
