//! The QPX protein-group view (`<prefix>.pg.parquet`) from protein inference.
//!
//! One row per group, targets and the decoy groups that won their family's
//! competition (`is_decoy`), so a reader can recompute the protein FDR. The
//! view is identification-only here: `label` and `intensity` are null.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    new_null_array, ArrayRef, BooleanBuilder, Float64Builder, Int32Builder, Int64Builder,
    ListBuilder, StringBuilder, StructBuilder,
};
use arrow::datatypes::{DataType, Field, Fields, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::feature_parquet::{
    cv_fields, cv_list_builder, list_type, score_fields, score_list_builder, struct_list,
};
use crate::protein_inference::{Inference, ProteinGroup};

/// QPX version of the pg view written here.
const PG_QPX_VERSION: &str = "1.1";
/// The columns a reader re-derives a group's identity from.
const PG_IDENTITY_COMPOSITE: &str = "pg_accessions,grouped_runs,label";

fn peptide_per_protein_fields() -> Fields {
    Fields::from(vec![
        Field::new("protein_name", DataType::Utf8, false),
        Field::new("peptide_count", DataType::Int32, false),
    ])
}

fn peptide_counts_fields() -> Fields {
    Fields::from(vec![
        Field::new("unique_sequences", DataType::Int32, false),
        Field::new("total_sequences", DataType::Int32, false),
    ])
}

fn feature_counts_fields() -> Fields {
    Fields::from(vec![
        Field::new("unique_features", DataType::Int32, false),
        Field::new("total_features", DataType::Int32, false),
    ])
}

fn additional_intensity_type() -> DataType {
    let pair = Fields::from(vec![
        Field::new("intensity_name", DataType::Utf8, false),
        Field::new("intensity_value", DataType::Float32, false),
    ]);
    list_type(DataType::Struct(Fields::from(vec![
        Field::new("label", DataType::Utf8, false),
        Field::new("intensities", list_type(DataType::Struct(pair)), true),
    ])))
}

/// The pg view's columns, in the QPX schema's order.
pub fn pg_schema() -> Schema {
    Schema::new(vec![
        Field::new("pg_id", DataType::Int64, false),
        Field::new("pg_accessions", list_type(DataType::Utf8), false),
        Field::new("pg_names", list_type(DataType::Utf8), true),
        Field::new("gg_accessions", list_type(DataType::Utf8), true),
        Field::new("gg_names", list_type(DataType::Utf8), true),
        Field::new("gg_qvalue", DataType::Float64, true),
        Field::new("anchor_protein", DataType::Utf8, false),
        Field::new("grouped_runs", list_type(DataType::Utf8), false),
        Field::new("global_qvalue", DataType::Float64, true),
        Field::new("pg_qvalue", DataType::Float64, true),
        Field::new("label", DataType::Utf8, true),
        Field::new("intensity", DataType::Float32, true),
        Field::new("additional_intensities", additional_intensity_type(), true),
        Field::new("is_decoy", DataType::Boolean, false),
        Field::new("contaminant", DataType::Boolean, true),
        Field::new(
            "peptides",
            list_type(DataType::Struct(peptide_per_protein_fields())),
            false,
        ),
        Field::new(
            "peptide_counts",
            DataType::Struct(peptide_counts_fields()),
            true,
        ),
        Field::new(
            "feature_counts",
            DataType::Struct(feature_counts_fields()),
            true,
        ),
        Field::new("sequence_coverage", DataType::Float32, true),
        Field::new("molecular_weight", DataType::Float32, true),
        Field::new(
            "additional_scores",
            list_type(DataType::Struct(score_fields())),
            true,
        ),
        Field::new("cv_params", list_type(DataType::Struct(cv_fields())), true),
    ])
}

fn build_batch(
    schema: &Arc<Schema>,
    rows: &[&ProteinGroup],
    runs: &[String],
    run_level: bool,
) -> std::io::Result<RecordBatch> {
    let n = rows.len();
    let mut pg_id = Int64Builder::new();
    let mut accessions = ListBuilder::new(StringBuilder::new()).with_field(Arc::new(Field::new(
        "element",
        DataType::Utf8,
        true,
    )));
    let mut anchor = StringBuilder::new();
    let mut grouped_runs = ListBuilder::new(StringBuilder::new()).with_field(Arc::new(Field::new(
        "element",
        DataType::Utf8,
        true,
    )));
    let mut global_q = Float64Builder::new();
    let mut pg_q = Float64Builder::new();
    let mut is_decoy = BooleanBuilder::new();
    let mut peptides = struct_list(
        peptide_per_protein_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(Int32Builder::new()),
        ],
    );
    let mut counts = StructBuilder::new(
        peptide_counts_fields(),
        vec![Box::new(Int32Builder::new()), Box::new(Int32Builder::new())],
    );
    let mut scores = score_list_builder();
    let mut cvs = cv_list_builder();
    for g in rows {
        pg_id.append_value(g.id as i64);
        accessions
            .values()
            .extend(g.accessions.iter().map(|a| Some(a.as_str())));
        accessions.append(true);
        anchor.append_value(&g.lead);
        grouped_runs
            .values()
            .extend(runs.iter().map(|r| Some(r.as_str())));
        grouped_runs.append(true);
        global_q.append_value(g.q_value);
        if run_level {
            pg_q.append_value(g.q_value);
        } else {
            pg_q.append_null();
        }
        is_decoy.append_value(g.is_decoy);
        {
            let sb = peptides.values();
            for (acc, count) in &g.peptides_per_accession {
                sb.field_builder::<StringBuilder>(0)
                    .expect("protein_name")
                    .append_value(acc);
                sb.field_builder::<Int32Builder>(1)
                    .expect("peptide_count")
                    .append_value(*count as i32);
                sb.append(true);
            }
            peptides.append(true);
        }
        counts
            .field_builder::<Int32Builder>(0)
            .expect("unique_sequences")
            .append_value(g.unique as i32);
        counts
            .field_builder::<Int32Builder>(1)
            .expect("total_sequences")
            .append_value(g.peptides as i32);
        counts.append(true);
        {
            let sb = scores.values();
            sb.field_builder::<StringBuilder>(0)
                .expect("score_name")
                .append_value("andes:best_peptide_score");
            sb.field_builder::<Float64Builder>(1)
                .expect("score_value")
                .append_value(g.score);
            sb.field_builder::<BooleanBuilder>(2)
                .expect("higher_better")
                .append_value(true);
            sb.append(true);
            scores.append(true);
        }
        {
            let sb = cvs.values();
            let mut cv = |name: &str, value: String| {
                sb.field_builder::<StringBuilder>(0)
                    .expect("cv_name")
                    .append_value(name);
                sb.field_builder::<StringBuilder>(1)
                    .expect("cv_value")
                    .append_value(value);
                sb.append(true);
            };
            cv("andes:families", g.families.join(";"));
            cv("andes:psms", g.psms.to_string());
            cv("andes:razor_peptides", g.razor.to_string());
            if !g.subsumed.is_empty() {
                cv("andes:subsumed", g.subsumed.join(";"));
            }
            cvs.append(true);
        }
    }
    let null = |name: &str| -> ArrayRef {
        let field = schema.field_with_name(name).expect("pg column");
        new_null_array(field.data_type(), n)
    };
    let columns: Vec<ArrayRef> = vec![
        Arc::new(pg_id.finish()),
        Arc::new(accessions.finish()),
        null("pg_names"),
        null("gg_accessions"),
        null("gg_names"),
        null("gg_qvalue"),
        Arc::new(anchor.finish()),
        Arc::new(grouped_runs.finish()),
        Arc::new(global_q.finish()),
        Arc::new(pg_q.finish()),
        null("label"),
        null("intensity"),
        null("additional_intensities"),
        Arc::new(is_decoy.finish()),
        null("contaminant"),
        Arc::new(peptides.finish()),
        Arc::new(counts.finish()),
        null("feature_counts"),
        null("sequence_coverage"),
        null("molecular_weight"),
        Arc::new(scores.finish()),
        Arc::new(cvs.finish()),
    ];
    RecordBatch::try_new(schema.clone(), columns).map_err(|e| std::io::Error::other(e.to_string()))
}

/// Write `inf`'s groups as the QPX pg view. `runs` are the run names the
/// groups were inferred over (`grouped_runs`). `global_qvalue` is the group's
/// q-value; `pg_qvalue` (run level) is set only when `run_level`, i.e. when
/// the inference covered one run.
pub fn write_pg_parquet(
    path: &Path,
    inf: &Inference,
    runs: &[String],
    run_level: bool,
) -> std::io::Result<()> {
    let now = chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string();
    let mut meta = std::collections::HashMap::new();
    meta.insert("qpx_version".to_string(), PG_QPX_VERSION.to_string());
    meta.insert("file_type".to_string(), "pg_file".to_string());
    meta.insert("creator".to_string(), "andes".to_string());
    meta.insert("software_provider".to_string(), "andes".to_string());
    meta.insert("creation_date".to_string(), now);
    meta.insert("compression_format".to_string(), "SNAPPY".to_string());
    meta.insert("uuid".to_string(), crate::qpx::uuid_v4());
    meta.insert("primary_key".to_string(), "pg_id".to_string());
    meta.insert(
        "identity_composite".to_string(),
        PG_IDENTITY_COMPOSITE.to_string(),
    );
    let schema = Arc::new(Schema::new_with_metadata(
        pg_schema().fields().clone(),
        meta,
    ));
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();
    crate::atomic::write_atomic(path, |path| {
        let file = std::fs::File::create(path)?;
        let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
            .map_err(|e| std::io::Error::other(e.to_string()))?;
        let rows: Vec<&ProteinGroup> = inf.groups.iter().chain(&inf.decoy_groups).collect();
        for chunk in rows.chunks(4096) {
            let batch = build_batch(&schema, chunk, runs, run_level)?;
            writer
                .write(&batch)
                .map_err(|e| std::io::Error::other(e.to_string()))?;
        }
        writer
            .close()
            .map_err(|e| std::io::Error::other(e.to_string()))?;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Array, BooleanArray, Float64Array, Int64Array, ListArray, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    fn group(id: usize, lead: &str, q: f64, decoy: bool) -> ProteinGroup {
        ProteinGroup {
            id,
            lead: lead.into(),
            families: vec![lead.into()],
            accessions: vec![lead.into(), format!("{lead}-2")],
            peptides: 3,
            unique: 2,
            razor: 1,
            psms: 7,
            score: 4.5,
            q_value: q,
            subsumed: vec!["S".into()],
            peptides_per_accession: vec![(lead.into(), 3), (format!("{lead}-2"), 1)],
            is_decoy: decoy,
        }
    }

    #[test]
    fn writes_the_qpx_pg_view_with_targets_and_decoys() {
        let inf = Inference {
            groups: vec![group(1, "P1", 0.0, false), group(2, "P2", 0.02, false)],
            peptide_group: Default::default(),
            decoy_groups: vec![group(3, "XXX_P9", 0.02, true)],
        };
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("exp.pg.parquet");
        write_pg_parquet(&path, &inf, &["run_a".into(), "run_b".into()], false).unwrap();

        let file = std::fs::File::open(&path).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        let meta = builder.schema().metadata().clone();
        assert_eq!(meta["file_type"], "pg_file");
        assert_eq!(meta["primary_key"], "pg_id");
        let batch = builder.build().unwrap().next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 3);
        assert_eq!(batch.schema().fields().len(), pg_schema().fields().len());
        let col = |n: &str| batch.column(batch.schema().index_of(n).unwrap()).clone();
        let ids = col("pg_id");
        let ids = ids.as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(ids.values(), &[1, 2, 3]);
        let decoy = col("is_decoy");
        let decoy = decoy.as_any().downcast_ref::<BooleanArray>().unwrap();
        assert!(!decoy.value(0) && decoy.value(2));
        let anchor = col("anchor_protein");
        assert_eq!(
            anchor
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .value(1),
            "P2"
        );
        let gq = col("global_qvalue");
        assert_eq!(
            gq.as_any().downcast_ref::<Float64Array>().unwrap().value(1),
            0.02
        );
        assert_eq!(
            col("pg_qvalue").null_count(),
            3,
            "dataset level: no run-level q"
        );
        let runs = col("grouped_runs");
        let runs = runs.as_any().downcast_ref::<ListArray>().unwrap();
        assert_eq!(runs.value(0).len(), 2);
        assert_eq!(col("label").null_count(), 3);
    }
}
