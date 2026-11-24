//! The finalize output rows: schema, evidence codes, chunked writer and reader, and the merge
//! into the final MBR table.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, Float64Array, Float64Builder, UInt16Array, UInt16Builder, UInt32Array,
    UInt32Builder, UInt8Array, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use output::mbr::{EvidenceType, MbrClaim, MbrParquetWriter, MbrRow};
use parquet::arrow::arrow_reader::ParquetRecordBatchReader;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use quant::claims::Claim;
use quant::library::Library;

use super::scores::{column, temp_reader, temp_writer, ROWS_KEY};
use super::{MbrFeature, MAX_OPEN_PARTS, OUT_GROUP_ROWS};

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct OutRow {
    pub(super) precursor: u32,
    pub(super) run: u16,
    pub(super) evidence_type: EvidenceType,
    pub(super) quantity: Option<f64>,
    pub(super) quantity_se: Option<f64>,
    pub(super) quantity_normalized: Option<f64>,
    pub(super) transfer_q: Option<f64>,
    pub(super) quant_q: Option<f64>,
    pub(super) apex_rt_s: Option<f64>,
    pub(super) measurement_q: Option<f64>,
    pub(super) library_q: f64,
    pub(super) mass_err_ppm: Option<f64>,
    pub(super) claim: Claim,
    pub(super) transfer_pep: Option<f64>,
}

impl OutRow {
    fn feature(&self) -> MbrFeature {
        MbrFeature {
            precursor: self.precursor,
            run: self.run,
            evidence_type: self.evidence_type,
            quantity: self.quantity,
            quantity_se: self.quantity_se,
            quantity_normalized: self.quantity_normalized,
            transfer_q: self.transfer_q,
            quant_q: self.quant_q,
            apex_rt_s: self.apex_rt_s,
            library_q: self.library_q,
            mass_err_ppm: self.mass_err_ppm,
            claim: self.claim,
        }
    }
}

pub(super) const EVIDENCE_CODES: [EvidenceType; 6] = [
    EvidenceType::Identified,
    EvidenceType::IdentifiedBelowThreshold,
    EvidenceType::TransferMs2,
    EvidenceType::TransferMs1,
    EvidenceType::BelowThreshold,
    EvidenceType::SharedFeature,
];

fn evidence_code(e: EvidenceType) -> u8 {
    EVIDENCE_CODES
        .iter()
        .position(|&x| x == e)
        .expect("every evidence type has a code") as u8
}

const OUT_VALUES: [&str; 7] = [
    "quantity",
    "quantity_se",
    "quantity_normalized",
    "transfer_q",
    "quant_q",
    "apex_rt_s",
    "measurement_q",
];

const OUT_CLAIM: [&str; 7] = [
    "claimant",
    "claim_dt_s",
    "claim_identified",
    "claim_co_apex",
    "claim_n",
    "claim_n_co_apex",
    "claim_n_identified_co_apex",
];

fn out_schema() -> Arc<Schema> {
    let mut fields = vec![
        Field::new("precursor", DataType::UInt32, false),
        Field::new("run", DataType::UInt16, false),
        Field::new("evidence_type", DataType::UInt8, false),
    ];
    fields.extend(
        OUT_VALUES
            .iter()
            .map(|n| Field::new(*n, DataType::Float64, true)),
    );
    fields.push(Field::new("library_q", DataType::Float64, false));
    fields.push(Field::new("mass_err_ppm", DataType::Float64, true));
    fields.extend(
        OUT_CLAIM
            .iter()
            .map(|n| Field::new(*n, DataType::Float64, true)),
    );
    fields.push(Field::new("transfer_pep", DataType::Float64, true));
    Arc::new(Schema::new(fields))
}

pub(super) struct OutWriter {
    pub(super) path: PathBuf,
    pub(super) schema: Arc<Schema>,
    writer: ArrowWriter<File>,
}

impl OutWriter {
    pub(super) fn create(path: &Path) -> Result<Self, String> {
        let schema = out_schema();
        let writer = temp_writer(
            path,
            ROWS_KEY,
            schema.clone(),
            Compression::ZSTD(ZstdLevel::try_new(1).map_err(|e| e.to_string())?),
            OUT_GROUP_ROWS,
        )?;
        Ok(Self {
            path: path.to_path_buf(),
            schema,
            writer,
        })
    }

    pub(super) fn write(&mut self, rows: &[OutRow]) -> Result<(), String> {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", self.path.display());
        for c in rows.chunks(OUT_GROUP_ROWS) {
            let mut precursor = UInt32Builder::with_capacity(c.len());
            let mut run = UInt16Builder::with_capacity(c.len());
            let mut ev = UInt8Builder::with_capacity(c.len());
            let mut values: Vec<Float64Builder> = (0..OUT_VALUES.len() + 3 + OUT_CLAIM.len())
                .map(|_| Float64Builder::with_capacity(c.len()))
                .collect();
            for r in c {
                precursor.append_value(r.precursor);
                run.append_value(r.run);
                ev.append_value(evidence_code(r.evidence_type));
                let v = [
                    r.quantity,
                    r.quantity_se,
                    r.quantity_normalized,
                    r.transfer_q,
                    r.quant_q,
                    r.apex_rt_s,
                    r.measurement_q,
                    Some(r.library_q),
                    r.mass_err_ppm,
                    r.claim.claimant.map(f64::from),
                    r.claim.dt_s.map(f64::from),
                    Some(f64::from(u8::from(r.claim.identified))),
                    Some(f64::from(u8::from(r.claim.co_apex))),
                    Some(f64::from(r.claim.n)),
                    Some(f64::from(r.claim.n_co_apex)),
                    Some(f64::from(r.claim.n_identified_co_apex)),
                    r.transfer_pep,
                ];
                for (b, x) in values.iter_mut().zip(v) {
                    b.append_option(x);
                }
            }
            let mut cols: Vec<ArrayRef> = vec![
                Arc::new(precursor.finish()),
                Arc::new(run.finish()),
                Arc::new(ev.finish()),
            ];
            cols.extend(values.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
            let b = RecordBatch::try_new(self.schema.clone(), cols).map_err(|e| err(&e))?;
            self.writer.write(&b).map_err(|e| err(&e))?;
        }
        Ok(())
    }

    pub(super) fn finish(self) -> Result<(), String> {
        let path = self.path;
        self.writer
            .close()
            .map(|_| ())
            .map_err(|e| format!("{}: {e}", path.display()))
    }
}

pub(super) struct OutReader {
    pub(super) path: PathBuf,
    pub(super) reader: ParquetRecordBatchReader,
    pub(super) rows: std::vec::IntoIter<OutRow>,
}

impl OutReader {
    pub(super) fn open(path: &Path) -> Result<Self, String> {
        Ok(Self {
            path: path.to_path_buf(),
            reader: temp_reader(path, ROWS_KEY, OUT_GROUP_ROWS)?,
            rows: Vec::new().into_iter(),
        })
    }

    fn next_row(&mut self) -> Result<Option<OutRow>, String> {
        loop {
            if let Some(r) = self.rows.next() {
                return Ok(Some(r));
            }
            let Some(b) = self.reader.next() else {
                return Ok(None);
            };
            let b = b.map_err(|e| format!("{}: {e}", self.path.display()))?;
            self.rows = self.decode(&b)?.into_iter();
        }
    }

    fn decode(&self, b: &RecordBatch) -> Result<Vec<OutRow>, String> {
        let path = self.path.as_path();
        let expected = OUT_VALUES.len() + 6 + OUT_CLAIM.len();
        if b.num_columns() != expected {
            return Err(format!(
                "{}: {} columns where {} were expected",
                path.display(),
                b.num_columns(),
                expected
            ));
        }
        let precursor = column::<UInt32Array>(b, 0, path)?;
        let run = column::<UInt16Array>(b, 1, path)?;
        let ev = column::<UInt8Array>(b, 2, path)?;
        let values: Vec<&Float64Array> = (3..b.num_columns())
            .map(|i| column::<Float64Array>(b, i, path))
            .collect::<Result<_, _>>()?;
        let opt = |k: usize, i: usize| (!values[k].is_null(i)).then(|| values[k].value(i));
        (0..b.num_rows())
            .map(|i| {
                let evidence_type = *EVIDENCE_CODES
                    .get(ev.value(i) as usize)
                    .ok_or_else(|| format!("{}: unknown evidence code", path.display()))?;
                Ok(OutRow {
                    precursor: precursor.value(i),
                    run: run.value(i),
                    evidence_type,
                    quantity: opt(0, i),
                    quantity_se: opt(1, i),
                    quantity_normalized: opt(2, i),
                    transfer_q: opt(3, i),
                    quant_q: opt(4, i),
                    apex_rt_s: opt(5, i),
                    measurement_q: opt(6, i),
                    library_q: values[7].value(i),
                    mass_err_ppm: opt(8, i),
                    claim: Claim {
                        claimant: opt(9, i).map(|x| x as u32),
                        dt_s: opt(10, i).map(|x| x as f32),
                        identified: opt(11, i).is_some_and(|x| x != 0.0),
                        co_apex: opt(12, i).is_some_and(|x| x != 0.0),
                        n: opt(13, i).map_or(0, |x| x as u16),
                        n_co_apex: opt(14, i).map_or(0, |x| x as u16),
                        n_identified_co_apex: opt(15, i).map_or(0, |x| x as u16),
                        ..Claim::default()
                    },
                    transfer_pep: opt(16, i),
                })
            })
            .collect()
    }
}

pub(super) fn merge_out(
    paths: &[PathBuf],
    dir: &Path,
    level: usize,
    sink: &mut dyn FnMut(&[OutRow]) -> Result<(), String>,
) -> Result<(), String> {
    if paths.len() > MAX_OPEN_PARTS {
        let mut merged = Vec::new();
        for (g, group) in paths.chunks(MAX_OPEN_PARTS).enumerate() {
            let p = dir.join(format!("rows-merged-{level}-{g:05}.parquet"));
            let mut w = OutWriter::create(&p)?;
            merge_out(group, dir, level + 1, &mut |rows| w.write(rows))?;
            w.finish()?;
            for q in group {
                let _ = std::fs::remove_file(q);
            }
            merged.push(p);
        }
        return merge_out(&merged, dir, level + 1, sink);
    }
    let mut readers: Vec<OutReader> = paths
        .iter()
        .map(|p| OutReader::open(p))
        .collect::<Result<_, _>>()?;
    let mut heads: Vec<Option<OutRow>> = Vec::with_capacity(readers.len());
    let mut heap = std::collections::BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        let head = r.next_row()?;
        if let Some(h) = &head {
            heap.push(std::cmp::Reverse((h.precursor, h.run, i)));
        }
        heads.push(head);
    }
    let mut out = Vec::with_capacity(OUT_GROUP_ROWS);
    while let Some(std::cmp::Reverse((_, _, i))) = heap.pop() {
        out.push(heads[i].take().expect("a queued table has a head row"));
        if let Some(h) = readers[i].next_row()? {
            heap.push(std::cmp::Reverse((h.precursor, h.run, i)));
            heads[i] = Some(h);
        }
        if out.len() == OUT_GROUP_ROWS {
            sink(&out)?;
            out.clear();
        }
    }
    if !out.is_empty() {
        sink(&out)?;
    }
    Ok(())
}

pub(super) fn write_merged(
    paths: &[PathBuf],
    dir: &Path,
    lib: &Library,
    stems: &[String],
    mbr_path: &Path,
    on_features: &mut dyn FnMut(Vec<MbrFeature>) -> Result<(), String>,
) -> Result<usize, String> {
    use std::sync::mpsc::sync_channel;
    let stopped = || "the output merge stopped early".to_string();
    let (to_table, table_rows) = sync_channel::<Arc<Vec<OutRow>>>(4);
    let (to_features, feature_rows) = sync_channel::<Arc<Vec<OutRow>>>(4);
    std::thread::scope(|scope| {
        let merger = scope.spawn(move || {
            merge_out(paths, dir, 0, &mut |rows| {
                let b = Arc::new(rows.to_vec());
                to_table.send(b.clone()).map_err(|_| stopped())?;
                to_features.send(b).map_err(|_| stopped())
            })
        });
        let table = scope.spawn(move || -> Result<(MbrParquetWriter, usize), String> {
            let mut writer = MbrParquetWriter::create(mbr_path)?;
            let mut n = 0;
            for rows in table_rows {
                for r in rows.iter() {
                    writer.push(&MbrRow {
                        precursor_key: &lib.meta[r.precursor as usize].key,
                        run: &stems[r.run as usize],
                        evidence_type: r.evidence_type,
                        quantity: r.quantity,
                        quantity_se: r.quantity_se,
                        quantity_normalized: r.quantity_normalized,
                        transfer_q: r.transfer_q,
                        quant_q: r.quant_q,
                        apex_rt_s: r.apex_rt_s,
                        library_q: r.library_q,
                        claim: r.evidence_type.is_transfer().then(|| MbrClaim {
                            claimed_by: r.claim.claimant.map(|p| lib.meta[p as usize].key.as_str()),
                            dt_s: r.claim.dt_s.map(f64::from),
                            identified: r.claim.identified,
                            co_apex: r.claim.co_apex,
                            n_co_apex: u32::from(r.claim.n_co_apex),
                            n_identified_co_apex: u32::from(r.claim.n_identified_co_apex),
                        }),
                        measurement_q: r.measurement_q,
                        transfer_pep: r.transfer_pep,
                    })?;
                }
                n += rows.len();
            }
            Ok((writer, n))
        });
        let mut handed = Ok(());
        for rows in feature_rows.iter() {
            handed = on_features(rows.iter().map(OutRow::feature).collect());
            if handed.is_err() {
                break;
            }
        }
        drop(feature_rows);
        let merged = merger.join().expect("output merge thread");
        let written = table.join().expect("output table thread");
        let (writer, n) = written?;
        handed?;
        merged?;
        writer.finish()?;
        Ok(n)
    })
}
