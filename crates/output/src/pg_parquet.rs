use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    new_null_array, ArrayRef, BooleanBuilder, Float32Builder, Float64Builder, Int32Builder,
    Int64Builder, ListBuilder, StringBuilder, StructBuilder,
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
use crate::protein_quant::{
    LabelQuantity, ProteinQuant, QuantUnit, PROTEIN_QUANT_VERSION, TOP3, TOTAL_INTENSITY,
};

const CONTAMINANT_PREFIX: &str = "CONTAMINANT_";

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FeatureCounts {
    pub unique: usize,
    pub total: usize,
}

const PG_QPX_VERSION: &str = "1.1";
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

fn intensity_pair_fields() -> Fields {
    Fields::from(vec![
        Field::new("intensity_name", DataType::Utf8, false),
        Field::new("intensity_value", DataType::Float32, false),
    ])
}

fn labelled_intensity_fields() -> Fields {
    Fields::from(vec![
        Field::new("label", DataType::Utf8, false),
        Field::new(
            "intensities",
            list_type(DataType::Struct(intensity_pair_fields())),
            true,
        ),
    ])
}

fn additional_intensity_type() -> DataType {
    list_type(DataType::Struct(labelled_intensity_fields()))
}

fn additional_intensity_builder() -> ListBuilder<StructBuilder> {
    struct_list(
        labelled_intensity_fields(),
        vec![
            Box::new(StringBuilder::new()),
            Box::new(struct_list(
                intensity_pair_fields(),
                vec![
                    Box::new(StringBuilder::new()),
                    Box::new(Float32Builder::new()),
                ],
            )),
        ],
    )
}

/// One row of the view: a group with identification evidence only (no unit), or a group's
/// quantity in one unit and label.
struct PgRow<'a> {
    id: i64,
    group: &'a ProteinGroup,
    quant: Option<(&'a QuantUnit, &'a str, LabelQuantity)>,
}

/// The rows of the view in output order: groups at the cut, targets then decoys; a group with
/// quantities gets one row per unit (in unit-name order) and label (in channel order), any
/// other group one identification row. Without quantities the row id is the group id, as
/// before quantities existed; with them it numbers the rows from 1.
fn pg_rows<'a>(
    groups: &[&'a ProteinGroup],
    units: &'a [QuantUnit],
    quant: Option<&ProteinQuant>,
) -> Vec<PgRow<'a>> {
    let mut order: Vec<usize> = (0..units.len()).collect();
    order.sort_by(|&a, &b| units[a].name.cmp(&units[b].name));
    let mut rows = Vec::new();
    for &g in groups {
        let cells = quant.map(|q| q.group(&g.lead)).unwrap_or_default();
        if cells.is_empty() {
            let id = if quant.is_some() {
                rows.len() as i64 + 1
            } else {
                g.id as i64
            };
            rows.push(PgRow {
                id,
                group: g,
                quant: None,
            });
            continue;
        }
        for &u in &order {
            let Some((_, labels)) = cells.iter().find(|(cu, _)| *cu == u) else {
                continue;
            };
            for (label, value) in units[u].labels.iter().zip(labels) {
                rows.push(PgRow {
                    id: rows.len() as i64 + 1,
                    group: g,
                    quant: Some((&units[u], label.as_str(), *value)),
                });
            }
        }
    }
    rows
}

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
    rows: &[PgRow<'_>],
    runs: &[String],
    run_level: bool,
    features: Option<&HashMap<usize, FeatureCounts>>,
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
    let mut contaminant = BooleanBuilder::new();
    let mut feature_counts = StructBuilder::new(
        feature_counts_fields(),
        vec![Box::new(Int32Builder::new()), Box::new(Int32Builder::new())],
    );
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
    let mut label = StringBuilder::new();
    let mut intensity = Float32Builder::new();
    let mut additional = additional_intensity_builder();
    for row in rows {
        let g = row.group;
        pg_id.append_value(row.id);
        accessions
            .values()
            .extend(g.accessions.iter().map(|a| Some(a.as_str())));
        accessions.append(true);
        anchor.append_value(&g.lead);
        let unit_runs = row.quant.map_or(runs, |(u, _, _)| u.runs.as_slice());
        grouped_runs
            .values()
            .extend(unit_runs.iter().map(|r| Some(r.as_str())));
        grouped_runs.append(true);
        match row.quant {
            Some((_, l, value)) => {
                label.append_value(l);
                intensity.append_option(value.total_intensity.map(|v| v as f32));
                let named: Vec<(&str, f64)> = [
                    (TOTAL_INTENSITY, value.total_intensity),
                    (TOP3, value.top3),
                ]
                .into_iter()
                .filter_map(|(n, v)| v.map(|v| (n, v)))
                .collect();
                if named.is_empty() {
                    additional.append(false);
                } else {
                    let sb = additional.values();
                    sb.field_builder::<StringBuilder>(0)
                        .expect("label")
                        .append_value(l);
                    let inner = sb
                        .field_builder::<ListBuilder<StructBuilder>>(1)
                        .expect("intensities");
                    let ib = inner.values();
                    for (n, v) in named {
                        ib.field_builder::<StringBuilder>(0)
                            .expect("intensity_name")
                            .append_value(n);
                        ib.field_builder::<Float32Builder>(1)
                            .expect("intensity_value")
                            .append_value(v as f32);
                        ib.append(true);
                    }
                    inner.append(true);
                    sb.append(true);
                    additional.append(true);
                }
            }
            None => {
                label.append_null();
                intensity.append_null();
                additional.append(false);
            }
        }
        global_q.append_value(g.q_value);
        if run_level {
            pg_q.append_value(g.q_value);
        } else {
            pg_q.append_null();
        }
        is_decoy.append_value(g.is_decoy);
        contaminant.append_value(
            g.accessions
                .iter()
                .any(|a| a.starts_with(CONTAMINANT_PREFIX)),
        );
        match features {
            Some(f) => {
                let c = f.get(&g.id).copied().unwrap_or_default();
                feature_counts
                    .field_builder::<Int32Builder>(0)
                    .expect("unique_features")
                    .append_value(c.unique as i32);
                feature_counts
                    .field_builder::<Int32Builder>(1)
                    .expect("total_features")
                    .append_value(c.total as i32);
                feature_counts.append(true);
            }
            None => {
                feature_counts
                    .field_builder::<Int32Builder>(0)
                    .expect("unique_features")
                    .append_null();
                feature_counts
                    .field_builder::<Int32Builder>(1)
                    .expect("total_features")
                    .append_null();
                feature_counts.append(false);
            }
        }
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
            if row.quant.is_some() {
                cv("quantification_method", TOTAL_INTENSITY.to_string());
                cv("andes:group_id", g.id.to_string());
            }
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
        Arc::new(label.finish()),
        Arc::new(intensity.finish()),
        Arc::new(additional.finish()),
        Arc::new(is_decoy.finish()),
        Arc::new(contaminant.finish()),
        Arc::new(peptides.finish()),
        Arc::new(counts.finish()),
        Arc::new(feature_counts.finish()),
        null("sequence_coverage"),
        null("molecular_weight"),
        Arc::new(scores.finish()),
        Arc::new(cvs.finish()),
    ];
    RecordBatch::try_new(schema.clone(), columns).map_err(|e| std::io::Error::other(e.to_string()))
}

pub fn write_pg_parquet(
    path: &Path,
    inf: &Inference,
    runs: &[String],
    run_level: bool,
    max_q: f64,
    features: Option<&HashMap<usize, FeatureCounts>>,
    quant: Option<&ProteinQuant>,
) -> std::io::Result<()> {
    let now = chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string();
    let mut meta = HashMap::new();
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
    if quant.is_some() {
        meta.insert(
            "andes_protein_quant_version".to_string(),
            PROTEIN_QUANT_VERSION.to_string(),
        );
        meta.insert(
            "andes_protein_quant".to_string(),
            format!(
                "intensity = {TOTAL_INTENSITY}: sum over the group's unique peptides; {TOP3}: mean of its \
                 3 most intense unique peptides in the unit (fewer when it has fewer), ranked by the \
                 intensity summed over the unit's labels; peptide = peptidoform summed over charges, \
                 fractions and rows; unit = run (label-free) or plex (isobaric); raw, not normalised"
            ),
        );
    }
    let units = quant.map(ProteinQuant::units).unwrap_or_default();
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
        let groups: Vec<&ProteinGroup> = inf
            .groups
            .iter()
            .chain(&inf.decoy_groups)
            .filter(|g| g.q_value <= max_q)
            .collect();
        let rows = pg_rows(&groups, &units, quant);
        for chunk in rows.chunks(4096) {
            let batch = build_batch(&schema, chunk, runs, run_level, features)?;
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
    use arrow::array::{
        Array, BooleanArray, Float64Array, Int32Array, Int64Array, ListArray, StringArray,
        StructArray,
    };
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
            shared_peptides: Default::default(),
            decoy_groups: vec![group(3, "XXX_P9", 0.02, true)],
        };
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("exp.pg.parquet");
        write_pg_parquet(
            &path,
            &inf,
            &["run_a".into(), "run_b".into()],
            false,
            1.0,
            None,
            None,
        )
        .unwrap();

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
        assert_eq!(col("feature_counts").null_count(), 3);
    }

    #[test]
    fn keeps_the_groups_at_the_cut_with_their_feature_counts() {
        let mut con = group(2, "CONTAMINANT_P7", 0.005, false);
        con.accessions = vec!["CONTAMINANT_P7".into()];
        let inf = Inference {
            groups: vec![group(1, "P1", 0.0, false), con, group(3, "P3", 0.2, false)],
            peptide_group: Default::default(),
            shared_peptides: Default::default(),
            decoy_groups: vec![
                group(4, "XXX_P8", 0.005, true),
                group(5, "XXX_P9", 0.3, true),
            ],
        };
        let counts: HashMap<usize, FeatureCounts> = [(
            1,
            FeatureCounts {
                unique: 4,
                total: 6,
            },
        )]
        .into();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("exp.pg.parquet");
        write_pg_parquet(&path, &inf, &["r".into()], false, 0.01, Some(&counts), None).unwrap();
        let file = std::fs::File::open(&path).unwrap();
        let batch = ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap()
            .build()
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        let col = |n: &str| batch.column(batch.schema().index_of(n).unwrap()).clone();
        let ids = col("pg_id");
        let ids = ids.as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(ids.values(), &[1, 2, 4], "groups above q 0.01 are left out");
        let con = col("contaminant");
        let con = con.as_any().downcast_ref::<BooleanArray>().unwrap();
        assert_eq!((con.value(0), con.value(1)), (false, true));
        let fc = col("feature_counts");
        let fc = fc.as_any().downcast_ref::<StructArray>().unwrap();
        let total = fc.column(1).as_any().downcast_ref::<Int32Array>().unwrap();
        assert_eq!((total.value(0), total.value(1)), (6, 0));
    }

    #[test]
    fn a_quantified_group_gets_one_row_per_unit_and_label_with_total_and_top3() {
        let inf = Inference {
            groups: vec![group(7, "P1", 0.0, false), group(8, "P2", 0.0, false)],
            peptide_group: Default::default(),
            shared_peptides: Default::default(),
            decoy_groups: vec![group(9, "XXX_P9", 0.0, true)],
        };
        let mut q = ProteinQuant::new();
        let l: Vec<String> = vec!["TMT126".into(), "TMT127N".into()];
        for (pep, v) in [("A", [40.0, 4.0]), ("B", [30.0, 0.0]), ("C", [20.0, 2.0]), ("D", [10.0, 1.0])] {
            q.add("P1", "plex2", "p2_f1", pep, &l, &v).unwrap();
        }
        q.add("P1", "plex1", "p1_f2", "A", &l, &[5.0, 6.0]).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("exp.pg.parquet");
        write_pg_parquet(&path, &inf, &["p1_f2".into(), "p2_f1".into()], false, 0.01, None, Some(&q)).unwrap();
        let file = std::fs::File::open(&path).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        assert_eq!(builder.schema().metadata()["andes_protein_quant_version"], PROTEIN_QUANT_VERSION);
        let batch = builder.build().unwrap().next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 6, "P1: 2 plexes x 2 channels; P2 and the decoy: one row each");
        let col = |n: &str| batch.column(batch.schema().index_of(n).unwrap()).clone();
        let ids = col("pg_id");
        assert_eq!(ids.as_any().downcast_ref::<Int64Array>().unwrap().values(), &[1, 2, 3, 4, 5, 6]);
        let label = col("label");
        let label = label.as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(
            (0..4).map(|i| label.value(i)).collect::<Vec<_>>(),
            ["TMT126", "TMT127N", "TMT126", "TMT127N"]
        );
        assert!(label.is_null(4) && label.is_null(5));
        let runs = col("grouped_runs");
        let runs = runs.as_any().downcast_ref::<ListArray>().unwrap();
        let first = runs.value(0);
        assert_eq!(first.as_any().downcast_ref::<StringArray>().unwrap().value(0), "p1_f2");
        assert_eq!(runs.value(4).len(), 2, "an identification row keeps every run");
        let it = col("intensity");
        let it = it.as_any().downcast_ref::<arrow::array::Float32Array>().unwrap();
        assert_eq!((it.value(0), it.value(1)), (5.0, 6.0));
        assert_eq!((it.value(2), it.value(3)), (100.0, 7.0));
        let ai = col("additional_intensities");
        let ai = ai.as_any().downcast_ref::<ListArray>().unwrap();
        let entry = ai.value(3);
        let entry = entry.as_any().downcast_ref::<StructArray>().unwrap();
        let pairs = entry.column(1).as_any().downcast_ref::<ListArray>().unwrap().value(0);
        let pairs = pairs.as_any().downcast_ref::<StructArray>().unwrap();
        let names = pairs.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        let values = pairs.column(1).as_any().downcast_ref::<arrow::array::Float32Array>().unwrap();
        assert_eq!((names.value(0), values.value(0)), (TOTAL_INTENSITY, 7.0));
        assert_eq!((names.value(1), values.value(1)), (TOP3, 3.0), "TMT127N averages A and C of the top three A, B, C; B has no value");
        assert!(ai.is_null(4));
    }
}
