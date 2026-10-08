//! `quantms.feature.parquet`: the QPX feature view (OpenMS
//! `QPXFeatureSchema::schema()`, as `IsobaricWorkflow -out_qpx` and
//! `ConsensusMapArrowExport` write it).
//!
//! One row is "a quantified peptidoform in one run": a label-free feature
//! carries a single `{label: "LFQ", intensity}`; an isobaric PSM carries one
//! entry per channel (`TMT126`, `TMT127N`, ...). Normalized or protein-level
//! values are never written here; they belong to the quantms / OpenMS layer.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayBuilder, ArrayRef, BooleanBuilder, Float32Builder, Float64Builder, Int16Builder,
    Int32Builder, Int64Builder, ListBuilder, StringBuilder, StructBuilder,
};
use arrow::datatypes::{DataType, Field, Fields, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::quant_out::ModRecord;

/// QPX version of the feature view written here (the psm/feature views are
/// unchanged between 1.0 and 1.1; OpenMS stamps 1.1).
const FEATURE_QPX_VERSION: &str = "1.1";
/// The columns a reader re-derives a feature's identity from (OpenMS
/// `QPXIdentity::FEATURE_COMPOSITE`).
const FEATURE_IDENTITY_COMPOSITE: &str = "run_file_name,peptidoform,charge,rt,scan,observed_mz";

/// A protein the peptide maps to: `(accession, start, end, pre, post)` with
/// 1-based inclusive positions and flanking residues when known.
pub type ProteinPosition = (
    String,
    Option<i32>,
    Option<i32>,
    Option<String>,
    Option<String>,
);

/// One quantified feature row.
#[derive(Debug, Clone, Default)]
pub struct FeatureRecord {
    pub feature_id: i64,
    pub sequence: String,
    pub peptidoform: String,
    pub modifications: Vec<ModRecord>,
    pub charge: i16,
    pub pep: Option<f64>,
    pub is_decoy: bool,
    pub calculated_mz: f32,
    pub observed_mz: f32,
    pub mass_error_ppm: Option<f32>,
    /// `(score_name, score_value, higher_better)`.
    pub additional_scores: Vec<(String, f64, bool)>,
    pub run_file_name: String,
    /// `(cv_name, cv_value)`.
    pub cv_params: Vec<(String, String)>,
    /// Scan number components of the identifying PSM.
    pub scan: Vec<i32>,
    pub rt: Option<f32>,
    pub missed_cleavages: Option<i16>,
    /// `(label, intensity)`.
    pub intensities: Vec<(String, f32)>,
    /// `(label, [(intensity_name, value)])`.
    pub additional_intensities: Vec<(String, Vec<(String, f32)>)>,
    /// Protein accessions with optional positions `(accession, start, end, pre, post)`.
    pub proteins: Vec<ProteinPosition>,
    pub id_run_file_name: Option<String>,
    pub rt_start: Option<f32>,
    pub rt_stop: Option<f32>,
    /// Row indices into `psms.parquet`.
    pub psm_ids: Vec<i64>,
}

// ── nested types (QPXPSMSchema / QPXFeatureSchema) ────────────────────────────

fn score_fields() -> Fields {
    Fields::from(vec![
        Field::new("score_name", DataType::Utf8, false),
        Field::new("score_value", DataType::Float64, false),
        Field::new("higher_better", DataType::Boolean, true),
    ])
}

fn list_type(item: DataType) -> DataType {
    DataType::List(Arc::new(Field::new("element", item, true)))
}

fn position_fields() -> Fields {
    Fields::from(vec![
        Field::new("position", DataType::Int32, false),
        Field::new("amino_acid", DataType::Utf8, true),
        Field::new("scores", list_type(DataType::Struct(score_fields())), true),
    ])
}

fn modification_fields() -> Fields {
    Fields::from(vec![
        Field::new("name", DataType::Utf8, false),
        Field::new("accession", DataType::Utf8, true),
        Field::new(
            "positions",
            list_type(DataType::Struct(position_fields())),
            false,
        ),
    ])
}

fn cv_fields() -> Fields {
    Fields::from(vec![
        Field::new("cv_name", DataType::Utf8, false),
        Field::new("cv_value", DataType::Utf8, false),
    ])
}

fn intensity_fields() -> Fields {
    Fields::from(vec![
        Field::new("label", DataType::Utf8, false),
        Field::new("intensity", DataType::Float32, false),
    ])
}

fn ai_entry_fields() -> Fields {
    Fields::from(vec![
        Field::new("intensity_name", DataType::Utf8, false),
        Field::new("intensity_value", DataType::Float32, false),
    ])
}

fn ai_fields() -> Fields {
    Fields::from(vec![
        Field::new("label", DataType::Utf8, false),
        Field::new(
            "intensities",
            list_type(DataType::Struct(ai_entry_fields())),
            false,
        ),
    ])
}

fn pg_accession_fields() -> Fields {
    Fields::from(vec![
        Field::new("accession", DataType::Utf8, false),
        Field::new("start", DataType::Int32, true),
        Field::new("end", DataType::Int32, true),
        Field::new("pre", DataType::Utf8, true),
        Field::new("post", DataType::Utf8, true),
    ])
}

fn pg_position_fields() -> Fields {
    Fields::from(vec![
        Field::new("protein_accession", DataType::Utf8, false),
        Field::new("start", DataType::Int32, false),
        Field::new("end", DataType::Int32, false),
    ])
}

/// The feature schema, in `QPXFeatureSchema::schema()` order.
pub fn feature_schema() -> Schema {
    Schema::new(vec![
        Field::new("feature_id", DataType::Int64, false),
        Field::new("sequence", DataType::Utf8, false),
        Field::new("peptidoform", DataType::Utf8, false),
        Field::new(
            "modifications",
            list_type(DataType::Struct(modification_fields())),
            true,
        ),
        Field::new("charge", DataType::Int16, false),
        Field::new("posterior_error_probability", DataType::Float64, true),
        Field::new("is_decoy", DataType::Boolean, false),
        Field::new("calculated_mz", DataType::Float32, false),
        Field::new("observed_mz", DataType::Float32, false),
        Field::new("mass_error_ppm", DataType::Float32, true),
        Field::new(
            "additional_scores",
            list_type(DataType::Struct(score_fields())),
            true,
        ),
        Field::new("predicted_rt", DataType::Float32, true),
        Field::new("run_file_name", DataType::Utf8, false),
        Field::new("cv_params", list_type(DataType::Struct(cv_fields())), true),
        Field::new("scan", list_type(DataType::Int32), false),
        Field::new("rt", DataType::Float32, true),
        Field::new("ion_mobility", DataType::Float32, true),
        Field::new("missed_cleavages", DataType::Int16, true),
        Field::new(
            "intensities",
            list_type(DataType::Struct(intensity_fields())),
            true,
        ),
        Field::new(
            "additional_intensities",
            list_type(DataType::Struct(ai_fields())),
            true,
        ),
        Field::new(
            "pg_accessions",
            list_type(DataType::Struct(pg_accession_fields())),
            true,
        ),
        Field::new("anchor_protein", DataType::Utf8, true),
        Field::new("unique", DataType::Boolean, true),
        Field::new("pg_global_qvalue", DataType::Float64, true),
        Field::new(
            "pg_positions",
            list_type(DataType::Struct(pg_position_fields())),
            true,
        ),
        Field::new("ion_mobility_start", DataType::Float32, true),
        Field::new("ion_mobility_stop", DataType::Float32, true),
        Field::new("gg_accessions", list_type(DataType::Utf8), true),
        Field::new("gg_names", list_type(DataType::Utf8), true),
        Field::new("id_run_file_name", DataType::Utf8, true),
        Field::new("rt_start", DataType::Float32, true),
        Field::new("rt_stop", DataType::Float32, true),
        Field::new("psm_ids", list_type(DataType::Int64), true),
    ])
}

// ── builders ──────────────────────────────────────────────────────────────────

fn struct_list(fields: Fields, children: Vec<Box<dyn ArrayBuilder>>) -> ListBuilder<StructBuilder> {
    let element = Arc::new(Field::new(
        "element",
        DataType::Struct(fields.clone()),
        true,
    ));
    ListBuilder::new(StructBuilder::new(fields, children)).with_field(element)
}

fn score_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        score_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(Float64Builder::new()),
            Box::new(BooleanBuilder::new()),
        ],
    )
}

fn position_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        position_fields(),
        vec![
            Box::new(Int32Builder::new()),
            Box::new(StringBuilder::new()),
            Box::new(score_list_builder()),
        ],
    )
}

fn modification_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        modification_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(StringBuilder::new()),
            Box::new(position_list_builder()),
        ],
    )
}

fn cv_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        cv_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(StringBuilder::new()),
        ],
    )
}

fn intensity_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        intensity_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(Float32Builder::new()),
        ],
    )
}

fn ai_entry_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        ai_entry_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(Float32Builder::new()),
        ],
    )
}

fn ai_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        ai_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(ai_entry_list_builder()),
        ],
    )
}

fn pg_accession_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        pg_accession_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(Int32Builder::new()),
            Box::new(Int32Builder::new()),
            Box::new(StringBuilder::new()),
            Box::new(StringBuilder::new()),
        ],
    )
}

fn pg_position_list_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        pg_position_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(Int32Builder::new()),
            Box::new(Int32Builder::new()),
        ],
    )
}

fn prim_list<T: ArrayBuilder>(values: T, inner: DataType) -> ListBuilder<T> {
    ListBuilder::new(values).with_field(Arc::new(Field::new("element", inner, true)))
}

/// Build the Arrow batch for `records` (the whole table in one batch).
pub fn build_feature_batch(records: &[FeatureRecord]) -> std::io::Result<RecordBatch> {
    let schema = Arc::new(feature_schema());

    let mut feature_id = Int64Builder::new();
    let mut sequence = StringBuilder::new();
    let mut peptidoform = StringBuilder::new();
    let mut modifications = modification_list_builder();
    let mut charge = Int16Builder::new();
    let mut pep = Float64Builder::new();
    let mut is_decoy = BooleanBuilder::new();
    let mut calculated_mz = Float32Builder::new();
    let mut observed_mz = Float32Builder::new();
    let mut mass_error_ppm = Float32Builder::new();
    let mut additional_scores = score_list_builder();
    let mut predicted_rt = Float32Builder::new();
    let mut run_file_name = StringBuilder::new();
    let mut cv_params = cv_list_builder();
    let mut scan = prim_list(Int32Builder::new(), DataType::Int32);
    let mut rt = Float32Builder::new();
    let mut ion_mobility = Float32Builder::new();
    let mut missed_cleavages = Int16Builder::new();
    let mut intensities = intensity_list_builder();
    let mut additional_intensities = ai_list_builder();
    let mut pg_accessions = pg_accession_list_builder();
    let mut anchor_protein = StringBuilder::new();
    let mut unique = BooleanBuilder::new();
    let mut pg_global_qvalue = Float64Builder::new();
    let mut pg_positions = pg_position_list_builder();
    let mut im_start = Float32Builder::new();
    let mut im_stop = Float32Builder::new();
    let mut gg_accessions = prim_list(StringBuilder::new(), DataType::Utf8);
    let mut gg_names = prim_list(StringBuilder::new(), DataType::Utf8);
    let mut id_run_file_name = StringBuilder::new();
    let mut rt_start = Float32Builder::new();
    let mut rt_stop = Float32Builder::new();
    let mut psm_ids = prim_list(Int64Builder::new(), DataType::Int64);

    for r in records {
        feature_id.append_value(r.feature_id);
        sequence.append_value(&r.sequence);
        peptidoform.append_value(&r.peptidoform);
        {
            let sb = modifications.values();
            for m in &r.modifications {
                sb.field_builder::<StringBuilder>(0)
                    .unwrap()
                    .append_value(&m.name);
                match &m.accession {
                    Some(a) => sb
                        .field_builder::<StringBuilder>(1)
                        .unwrap()
                        .append_value(a),
                    None => sb.field_builder::<StringBuilder>(1).unwrap().append_null(),
                }
                let positions = sb.field_builder::<ListBuilder<StructBuilder>>(2).unwrap();
                let pb = positions.values();
                pb.field_builder::<Int32Builder>(0)
                    .unwrap()
                    .append_value(m.position);
                pb.field_builder::<StringBuilder>(1)
                    .unwrap()
                    .append_value(m.amino_acid.to_string());
                // scores: null (andes does not localize)
                pb.field_builder::<ListBuilder<StructBuilder>>(2)
                    .unwrap()
                    .append_null();
                pb.append(true);
                positions.append(true);
                sb.append(true);
            }
            modifications.append(true);
        }
        charge.append_value(r.charge);
        match r.pep {
            Some(v) => pep.append_value(v),
            None => pep.append_null(),
        }
        is_decoy.append_value(r.is_decoy);
        calculated_mz.append_value(r.calculated_mz);
        observed_mz.append_value(r.observed_mz);
        match r.mass_error_ppm {
            Some(v) => mass_error_ppm.append_value(v),
            None => mass_error_ppm.append_null(),
        }
        {
            let sb = additional_scores.values();
            for (name, value, higher) in &r.additional_scores {
                sb.field_builder::<StringBuilder>(0)
                    .unwrap()
                    .append_value(name);
                sb.field_builder::<Float64Builder>(1)
                    .unwrap()
                    .append_value(*value);
                sb.field_builder::<BooleanBuilder>(2)
                    .unwrap()
                    .append_value(*higher);
                sb.append(true);
            }
            additional_scores.append(true);
        }
        predicted_rt.append_null();
        run_file_name.append_value(&r.run_file_name);
        {
            let sb = cv_params.values();
            for (name, value) in &r.cv_params {
                sb.field_builder::<StringBuilder>(0)
                    .unwrap()
                    .append_value(name);
                sb.field_builder::<StringBuilder>(1)
                    .unwrap()
                    .append_value(value);
                sb.append(true);
            }
            cv_params.append(true);
        }
        for s in &r.scan {
            scan.values().append_value(*s);
        }
        scan.append(true);
        match r.rt {
            Some(v) => rt.append_value(v),
            None => rt.append_null(),
        }
        ion_mobility.append_null();
        match r.missed_cleavages {
            Some(v) => missed_cleavages.append_value(v),
            None => missed_cleavages.append_null(),
        }
        {
            let sb = intensities.values();
            for (label, value) in &r.intensities {
                sb.field_builder::<StringBuilder>(0)
                    .unwrap()
                    .append_value(label);
                sb.field_builder::<Float32Builder>(1)
                    .unwrap()
                    .append_value(*value);
                sb.append(true);
            }
            intensities.append(true);
        }
        {
            let sb = additional_intensities.values();
            for (label, entries) in &r.additional_intensities {
                sb.field_builder::<StringBuilder>(0)
                    .unwrap()
                    .append_value(label);
                let inner = sb.field_builder::<ListBuilder<StructBuilder>>(1).unwrap();
                let ib = inner.values();
                for (name, value) in entries {
                    ib.field_builder::<StringBuilder>(0)
                        .unwrap()
                        .append_value(name);
                    ib.field_builder::<Float32Builder>(1)
                        .unwrap()
                        .append_value(*value);
                    ib.append(true);
                }
                inner.append(true);
                sb.append(true);
            }
            additional_intensities.append(true);
        }
        {
            let sb = pg_accessions.values();
            for (acc, start, end, pre, post) in &r.proteins {
                sb.field_builder::<StringBuilder>(0)
                    .unwrap()
                    .append_value(acc);
                sb.field_builder::<Int32Builder>(1)
                    .unwrap()
                    .append_option(*start);
                sb.field_builder::<Int32Builder>(2)
                    .unwrap()
                    .append_option(*end);
                sb.field_builder::<StringBuilder>(3)
                    .unwrap()
                    .append_option(pre.as_deref());
                sb.field_builder::<StringBuilder>(4)
                    .unwrap()
                    .append_option(post.as_deref());
                sb.append(true);
            }
            pg_accessions.append(true);
        }
        match r.proteins.first() {
            Some((acc, ..)) => anchor_protein.append_value(acc),
            None => anchor_protein.append_null(),
        }
        if r.proteins.is_empty() {
            unique.append_null();
        } else {
            unique.append_value(r.proteins.len() == 1);
        }
        pg_global_qvalue.append_null();
        {
            let sb = pg_positions.values();
            for (acc, start, end, _, _) in &r.proteins {
                if let (Some(s), Some(e)) = (start, end) {
                    sb.field_builder::<StringBuilder>(0)
                        .unwrap()
                        .append_value(acc);
                    sb.field_builder::<Int32Builder>(1)
                        .unwrap()
                        .append_value(*s);
                    sb.field_builder::<Int32Builder>(2)
                        .unwrap()
                        .append_value(*e);
                    sb.append(true);
                }
            }
            pg_positions.append(true);
        }
        im_start.append_null();
        im_stop.append_null();
        gg_accessions.append_null();
        gg_names.append_null();
        match &r.id_run_file_name {
            Some(v) => id_run_file_name.append_value(v),
            None => id_run_file_name.append_null(),
        }
        match r.rt_start {
            Some(v) => rt_start.append_value(v),
            None => rt_start.append_null(),
        }
        match r.rt_stop {
            Some(v) => rt_stop.append_value(v),
            None => rt_stop.append_null(),
        }
        if r.psm_ids.is_empty() {
            psm_ids.append_null();
        } else {
            for p in &r.psm_ids {
                psm_ids.values().append_value(*p);
            }
            psm_ids.append(true);
        }
    }

    let columns: Vec<ArrayRef> = vec![
        Arc::new(feature_id.finish()),
        Arc::new(sequence.finish()),
        Arc::new(peptidoform.finish()),
        Arc::new(modifications.finish()),
        Arc::new(charge.finish()),
        Arc::new(pep.finish()),
        Arc::new(is_decoy.finish()),
        Arc::new(calculated_mz.finish()),
        Arc::new(observed_mz.finish()),
        Arc::new(mass_error_ppm.finish()),
        Arc::new(additional_scores.finish()),
        Arc::new(predicted_rt.finish()),
        Arc::new(run_file_name.finish()),
        Arc::new(cv_params.finish()),
        Arc::new(scan.finish()),
        Arc::new(rt.finish()),
        Arc::new(ion_mobility.finish()),
        Arc::new(missed_cleavages.finish()),
        Arc::new(intensities.finish()),
        Arc::new(additional_intensities.finish()),
        Arc::new(pg_accessions.finish()),
        Arc::new(anchor_protein.finish()),
        Arc::new(unique.finish()),
        Arc::new(pg_global_qvalue.finish()),
        Arc::new(pg_positions.finish()),
        Arc::new(im_start.finish()),
        Arc::new(im_stop.finish()),
        Arc::new(gg_accessions.finish()),
        Arc::new(gg_names.finish()),
        Arc::new(id_run_file_name.finish()),
        Arc::new(rt_start.finish()),
        Arc::new(rt_stop.finish()),
        Arc::new(psm_ids.finish()),
    ];
    RecordBatch::try_new(schema, columns).map_err(|e| std::io::Error::other(e.to_string()))
}

/// Write `records` to `path` as a QPX feature file with the schema metadata
/// OpenMS readers key on (`file_type = feature_file`).
pub fn write_feature_parquet(path: &Path, records: &[FeatureRecord]) -> std::io::Result<()> {
    let now = chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string();
    let mut meta = std::collections::HashMap::new();
    meta.insert("qpx_version".to_string(), FEATURE_QPX_VERSION.to_string());
    meta.insert("file_type".to_string(), "feature_file".to_string());
    meta.insert("creator".to_string(), "andes".to_string());
    meta.insert("software_provider".to_string(), "andes".to_string());
    meta.insert("creation_date".to_string(), now);
    meta.insert("compression_format".to_string(), "SNAPPY".to_string());
    meta.insert("uuid".to_string(), crate::qpx::uuid_v4());
    meta.insert("primary_key".to_string(), "feature_id".to_string());
    meta.insert(
        "identity_composite".to_string(),
        FEATURE_IDENTITY_COMPOSITE.to_string(),
    );
    let schema = Arc::new(Schema::new_with_metadata(
        feature_schema().fields().clone(),
        meta,
    ));
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| std::io::Error::other(e.to_string()))?;
    for chunk in records.chunks(4096) {
        let batch = build_feature_batch(chunk)?;
        let batch = RecordBatch::try_new(schema.clone(), batch.columns().to_vec())
            .map_err(|e| std::io::Error::other(e.to_string()))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::other(e.to_string()))?;
        writer
            .flush()
            .map_err(|e| std::io::Error::other(e.to_string()))?;
    }
    writer
        .close()
        .map_err(|e| std::io::Error::other(e.to_string()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    fn record(label: &str) -> FeatureRecord {
        FeatureRecord {
            feature_id: 7,
            sequence: "PEPTIDEK".into(),
            peptidoform: "PEPTIDEK[UNIMOD:737]".into(),
            modifications: vec![ModRecord {
                name: "TMT10plex".into(),
                accession: Some("UNIMOD:737".into()),
                position: 8,
                amino_acid: 'K',
            }],
            charge: 2,
            pep: Some(0.01),
            is_decoy: false,
            calculated_mz: 500.25,
            observed_mz: 500.2512,
            mass_error_ppm: Some(2.4),
            additional_scores: vec![("precursor_purity".into(), 0.9, true)],
            run_file_name: "run1".into(),
            cv_params: vec![],
            scan: vec![1234],
            rt: Some(600.0),
            missed_cleavages: Some(0),
            intensities: vec![(label.to_string(), 1e5), ("TMT127N".into(), 2e5)],
            additional_intensities: vec![("TMT126".into(), vec![("raw".into(), 1.1e5)])],
            proteins: vec![(
                "P12345".into(),
                Some(10),
                Some(17),
                Some("K".into()),
                Some("A".into()),
            )],
            id_run_file_name: Some("run1".into()),
            rt_start: None,
            rt_stop: None,
            psm_ids: vec![42],
        }
    }

    #[test]
    fn writes_the_qpx_feature_schema_and_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("quantms.feature.parquet");
        write_feature_parquet(&p, &[record("TMT126"), record("TMT126")]).unwrap();
        let file = std::fs::File::open(&p).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        let schema = builder.schema().clone();
        let names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
        assert_eq!(names[0], "feature_id");
        assert_eq!(names[18], "intensities");
        assert_eq!(names[names.len() - 1], "psm_ids");
        assert_eq!(names.len(), 33);
        let md = schema.metadata();
        assert_eq!(
            md.get("file_type").map(String::as_str),
            Some("feature_file")
        );
        assert_eq!(
            md.get("primary_key").map(String::as_str),
            Some("feature_id")
        );
        assert_eq!(md.get("qpx_version").map(String::as_str), Some("1.1"));
        let mut rows = 0;
        for b in builder.build().unwrap() {
            rows += b.unwrap().num_rows();
        }
        assert_eq!(rows, 2);
    }

    #[test]
    fn empty_table_is_valid() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("f.parquet");
        write_feature_parquet(&p, &[]).unwrap();
        let file = std::fs::File::open(&p).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        assert_eq!(builder.schema().fields().len(), 33);
    }

    #[test]
    fn multiple_batches_preserve_feature_ids_and_psm_links() {
        use arrow::array::{Int64Array, ListArray};
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("features.parquet");
        let records: Vec<_> = (0..4097)
            .map(|i| {
                let mut r = record("LFQ");
                r.feature_id = i;
                r.psm_ids = vec![i + 10];
                r
            })
            .collect();
        write_feature_parquet(&path, &records).unwrap();
        let builder =
            ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(path).unwrap()).unwrap();
        assert_eq!(builder.metadata().num_row_groups(), 2);
        let mut row = 0i64;
        for batch in builder.build().unwrap() {
            let batch = batch.unwrap();
            let ids = batch
                .column_by_name("feature_id")
                .unwrap()
                .as_any()
                .downcast_ref::<Int64Array>()
                .unwrap();
            let links = batch
                .column_by_name("psm_ids")
                .unwrap()
                .as_any()
                .downcast_ref::<ListArray>()
                .unwrap();
            for i in 0..batch.num_rows() {
                assert_eq!(ids.value(i), row);
                let link = links.value(i);
                assert_eq!(
                    link.as_any().downcast_ref::<Int64Array>().unwrap().value(0),
                    row + 10
                );
                row += 1;
            }
        }
        assert_eq!(row, 4097);
    }
}
