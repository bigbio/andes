//! Per-project error of a pooled crosslink FDR: at each q cut of the pooled list, the target and
//! decoy counts of every project's share, and the decoy estimate of that share. The table of
//! quantified cross-link features is in `features`.

pub mod features;

use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float64Array, StringArray, UInt32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::versioned::{read_table, write_table, ReadTable, TableVersion, BATCH_ROWS};

const TABLE: TableVersion = TableVersion {
    key: "andes.xl_projects.schema",
    version: "1.0",
    major: 1,
    what: "per-project crosslink error table",
};

/// One project's share of the CSMs a pooled q cut accepts, for one link class.
#[derive(Debug, Clone, PartialEq)]
pub struct XlProjectRow {
    pub project: String,
    pub link: String,
    pub q_cut: f64,
    pub tt: usize,
    pub td: usize,
    pub dd: usize,
    /// The decoy estimate of the false target-target share: (TD - DD) / TT, or for
    /// homotypic links (TD + DD) / TT; 0 when the share has no target-target CSM.
    pub decoy_fdr: f64,
}

fn schema() -> Schema {
    Schema::new(vec![
        Field::new("project", DataType::Utf8, false),
        Field::new("link", DataType::Utf8, false),
        Field::new("q_cut", DataType::Float64, false),
        Field::new("tt", DataType::UInt32, false),
        Field::new("td", DataType::UInt32, false),
        Field::new("dd", DataType::UInt32, false),
        Field::new("decoy_fdr", DataType::Float64, false),
    ])
}

fn batch(schema: &Arc<Schema>, rows: &[XlProjectRow]) -> Result<RecordBatch, String> {
    let u32s = |f: &dyn Fn(&XlProjectRow) -> usize| -> ArrayRef {
        Arc::new(UInt32Array::from_iter_values(
            rows.iter().map(|r| u32::try_from(f(r)).unwrap_or(u32::MAX)),
        ))
    };
    let cols: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from_iter_values(
            rows.iter().map(|r| r.project.as_str()),
        )),
        Arc::new(StringArray::from_iter_values(
            rows.iter().map(|r| r.link.as_str()),
        )),
        Arc::new(Float64Array::from_iter_values(rows.iter().map(|r| r.q_cut))),
        u32s(&|r| r.tt),
        u32s(&|r| r.td),
        u32s(&|r| r.dd),
        Arc::new(Float64Array::from_iter_values(
            rows.iter().map(|r| r.decoy_fdr),
        )),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

pub fn write_xl_projects(path: &Path, rows: &[XlProjectRow]) -> Result<(), String> {
    let schema = Arc::new(schema());
    let batches = rows.chunks(BATCH_ROWS).map(|c| batch(&schema, c));
    write_table(path, &TABLE, TABLE.version, schema.clone(), batches, &[])
}

pub fn read_xl_projects(path: &Path) -> Result<ReadTable, String> {
    read_table(path, &TABLE)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::versioned::column;

    #[test]
    fn rows_round_trip_with_the_schema_version() {
        let t = tempfile::tempdir().unwrap();
        let p = t.path().join("c.xl.projects.parquet");
        let rows = vec![
            XlProjectRow {
                project: "PXD1".into(),
                link: "intra".into(),
                q_cut: 0.01,
                tt: 200,
                td: 5,
                dd: 1,
                decoy_fdr: 0.02,
            },
            XlProjectRow {
                project: "PXD2".into(),
                link: "homotypic".into(),
                q_cut: 0.05,
                tt: 0,
                td: 0,
                dd: 0,
                decoy_fdr: 0.0,
            },
        ];
        write_xl_projects(&p, &rows).unwrap();
        let back = read_xl_projects(&p).unwrap();
        assert_eq!(back.meta.get(TABLE.key).map(String::as_str), Some("1.0"));
        let b = &back.batches[0];
        assert_eq!(b.num_rows(), 2);
        let project = column::<StringArray>(b, "project", &p).unwrap();
        let tt = column::<UInt32Array>(b, "tt", &p).unwrap();
        let fdr = column::<Float64Array>(b, "decoy_fdr", &p).unwrap();
        assert_eq!(
            (project.value(0), tt.value(0), fdr.value(0)),
            ("PXD1", 200, 0.02)
        );
        assert_eq!(project.value(1), "PXD2");
    }
}
