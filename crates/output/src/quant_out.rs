//! Quantification tables: the flat `tmt.tsv` / `lfq.tsv` files and the record
//! type the QPX feature parquet is written from.
//!
//! The TSVs keep a stable wide layout so notebooks written for `tmt.tsv` /
//! `lfq.tsv` keep working: one row per quantified PSM with one column per
//! channel, and one row per peptidoform-charge with one column per run.

use std::io::{self, BufWriter, Write};
use std::path::Path;

use quant::isobaric::Plex;
use quant::lfq::FeatureQuant;

/// One modification on a quantified peptidoform (for the parquet).
#[derive(Debug, Clone, PartialEq)]
pub struct ModRecord {
    pub name: String,
    pub accession: Option<String>,
    /// 1-based residue position.
    pub position: i32,
    pub amino_acid: char,
}

/// Identification side of a quantified row, shared by both modules.
#[derive(Debug, Clone)]
pub struct QuantId {
    /// Index of the input file (`--spectrum` order).
    pub file_idx: usize,
    /// File stem of the run (`run_file_name`).
    pub run: String,
    pub spec_id: String,
    pub scan: i32,
    pub rt_seconds: Option<f64>,
    pub sequence: String,
    pub peptidoform: String,
    pub modifications: Vec<ModRecord>,
    pub charge: u8,
    pub proteins: Vec<String>,
    pub is_decoy: bool,
    pub calculated_mz: f64,
    pub observed_mz: f64,
    pub q_value: Option<f64>,
    pub pep: Option<f64>,
    pub missed_cleavages: Option<i16>,
    /// 0-based index of the PSM row in `psms.parquet` (the feature's `psm_ids`).
    pub psm_index: i64,
}

/// One quantified PSM of an isobaric run.
#[derive(Debug, Clone)]
pub struct TmtRow {
    pub id: QuantId,
    /// Native id of the scan the reporters were read from (the MS2 itself, or
    /// its MS3).
    pub quant_scan_id: String,
    pub quant_ms_level: u8,
    pub purity: Option<f32>,
    pub raw: Vec<f32>,
    pub corrected: Vec<f32>,
    /// Per-channel signal-to-noise of the corrected intensities (the raw ones
    /// without a correction matrix), when the reporter scan records noise.
    pub sn: Option<Vec<f32>>,
    /// Suggested weight of this PSM in a protein or peptide roll-up (0 for a
    /// PSM the quantification gate removed).
    pub quant_weight: f32,
}

impl TmtRow {
    /// Summed raw reporter intensity.
    pub fn intensity_sum(&self) -> f32 {
        self.raw.iter().map(|v| v.max(0.0)).sum()
    }

    /// Summed reporter signal-to-noise, when noise is known.
    pub fn sn_sum(&self) -> Option<f32> {
        self.sn.as_ref().map(|sn| sn.iter().sum())
    }
}

/// One quantified precursor feature of a label-free run.
#[derive(Debug, Clone)]
pub struct LfqRow {
    pub id: QuantId,
    pub feature: FeatureQuant,
    pub score: f32,
    pub feature_q_value: Option<f64>,
    /// The decoy twin's score (None when the twin found no peak).
    pub decoy_score: Option<f32>,
}

fn fmt_opt<T: std::fmt::Display>(v: Option<T>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

fn fmt_f32(v: f32) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!("{v:.6e}")
    }
}

/// Write `tmt.tsv`: one row per quantified PSM, the corrected channel
/// intensities in kit order (`tmt_126 … tmt_131`), the raw values in a
/// `raw_` block when a correction matrix was applied, then the summed reporter
/// intensity, summed signal-to-noise and roll-up weight, and a per-channel
/// `sn_` block (empty cells when the reporter scan carries no noise).
pub fn write_tmt_tsv(path: &Path, plex: Plex, rows: &[TmtRow], corrected: bool) -> io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    let prefix = if matches!(plex, Plex::Itraq4 | Plex::Itraq8) {
        "itraq"
    } else {
        "tmt"
    };
    let mut header = vec![
        "filename",
        "scannr",
        "spec_id",
        "quant_scan",
        "quant_ms_level",
        "rt",
        "peptide",
        "charge",
        "proteins",
        "is_decoy",
        "psm_q_value",
        "psm_pep",
        "purity",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    for c in plex.channels() {
        header.push(format!("{prefix}_{}", c.name));
    }
    if corrected {
        for c in plex.channels() {
            header.push(format!("raw_{prefix}_{}", c.name));
        }
    }
    for name in ["reporter_intensity_sum", "reporter_sn_sum", "quant_weight"] {
        header.push(name.to_string());
    }
    for c in plex.channels() {
        header.push(format!("sn_{prefix}_{}", c.name));
    }
    writeln!(w, "{}", header.join("\t"))?;
    for r in rows {
        let id = &r.id;
        let mut cols: Vec<String> = vec![
            id.run.clone(),
            id.scan.to_string(),
            id.spec_id.clone(),
            r.quant_scan_id.clone(),
            r.quant_ms_level.to_string(),
            fmt_opt(id.rt_seconds.map(|t| format!("{t:.3}"))),
            id.peptidoform.clone(),
            id.charge.to_string(),
            id.proteins.join(";"),
            (id.is_decoy as u8).to_string(),
            fmt_opt(id.q_value),
            fmt_opt(id.pep),
            fmt_opt(r.purity.map(|p| format!("{p:.4}"))),
        ];
        let shown = if corrected { &r.corrected } else { &r.raw };
        cols.extend(shown.iter().map(|&v| fmt_f32(v)));
        if corrected {
            cols.extend(r.raw.iter().map(|&v| fmt_f32(v)));
        }
        cols.push(fmt_f32(r.intensity_sum()));
        cols.push(fmt_opt(r.sn_sum().map(|v| format!("{v:.2}"))));
        cols.push(format!("{:.4}", r.quant_weight));
        match &r.sn {
            Some(sn) => cols.extend(sn.iter().map(|v| format!("{v:.2}"))),
            None => cols.extend(plex.channels().iter().map(|_| String::new())),
        }
        writeln!(w, "{}", cols.join("\t"))?;
    }
    w.flush()
}

/// Write `lfq.tsv` (wide): one row per (peptidoform, charge) with the feature
/// area in every run (`runs` gives the column order), plus the best feature
/// q-value and score across runs.
pub fn write_lfq_tsv(path: &Path, runs: &[String], rows: &[LfqRow]) -> io::Result<()> {
    use std::collections::BTreeMap;
    // key → (proteins, per-run row)
    // (peptidoform, charge) → (proteins, one slot per run)
    type RunSlots<'a> = (Vec<String>, Vec<Option<&'a LfqRow>>);
    let mut groups: BTreeMap<(String, u8), RunSlots<'_>> = BTreeMap::new();
    for r in rows {
        let e = groups
            .entry((r.id.peptidoform.clone(), r.id.charge))
            .or_insert_with(|| (r.id.proteins.clone(), vec![None; runs.len()]));
        if r.id.file_idx < runs.len() {
            // Keep the best-scoring feature per run (two PSM anchors can map to
            // the same peptidoform-charge in one run only through the
            // best-PSM selection upstream; this is a safety net).
            let slot = &mut e.1[r.id.file_idx];
            if slot.is_none_or(|s| r.score > s.score) {
                *slot = Some(r);
            }
        }
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    let mut header = vec![
        "peptide",
        "charge",
        "proteins",
        "feature_q_value",
        "score",
        "cosine",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    header.extend(runs.iter().cloned());
    writeln!(w, "{}", header.join("\t"))?;
    for ((pep, z), (proteins, per_run)) in groups {
        let best = per_run.iter().flatten().max_by(|a, b| {
            a.score
                .partial_cmp(&b.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let (q, score, cos) = match best {
            Some(b) => (
                fmt_opt(b.feature_q_value.map(|q| format!("{q:.5}"))),
                format!("{:.4}", b.score),
                format!("{:.4}", b.feature.cosine),
            ),
            None => (String::new(), String::new(), String::new()),
        };
        let mut cols = vec![pep, z.to_string(), proteins.join(";"), q, score, cos];
        for slot in per_run {
            cols.push(match slot {
                Some(r) => format!("{:.6e}", r.feature.area),
                None => String::new(),
            });
        }
        writeln!(w, "{}", cols.join("\t"))?;
    }
    w.flush()
}

/// Write `lfq_features.tsv` (long): one row per quantified feature in one run
/// with its apex, boundaries, areas and the decoy-competition result.
pub fn write_lfq_features_tsv(path: &Path, rows: &[LfqRow]) -> io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "filename\tpeptide\tcharge\tproteins\tscannr\tspec_id\tpsm_q_value\tcalc_mz\t\
         rt_anchor\trt_apex\trt_start\trt_stop\tn_scans\tarea\tapex_intensity\t\
         isotope_areas\tcosine\tn_isotopes\tscore\tdecoy_score\tfeature_q_value"
    )?;
    for r in rows {
        let f = &r.feature;
        let areas: Vec<String> = f.isotope_areas.iter().map(|a| format!("{a:.4e}")).collect();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{:.6e}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{}",
            r.id.run,
            r.id.peptidoform,
            r.id.charge,
            r.id.proteins.join(";"),
            r.id.scan,
            r.id.spec_id,
            fmt_opt(r.id.q_value),
            r.id.calculated_mz,
            fmt_opt(r.id.rt_seconds.map(|t| format!("{t:.3}"))),
            f.apex_rt,
            f.rt_start,
            f.rt_stop,
            f.n_scans,
            f.area,
            fmt_f32(f.apex_intensity),
            areas.join(";"),
            f.cosine,
            f.n_isotopes_found,
            r.score,
            fmt_opt(r.decoy_score.map(|s| format!("{s:.4}"))),
            fmt_opt(r.feature_q_value.map(|q| format!("{q:.5}"))),
        )?;
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn id(run: &str, file_idx: usize, pep: &str, z: u8) -> QuantId {
        QuantId {
            file_idx,
            run: run.to_string(),
            spec_id: format!("{run}.scan=5"),
            scan: 5,
            rt_seconds: Some(100.0),
            sequence: pep.to_string(),
            peptidoform: pep.to_string(),
            modifications: vec![],
            charge: z,
            proteins: vec!["P1".into()],
            is_decoy: false,
            calculated_mz: 500.0,
            observed_mz: 500.001,
            q_value: Some(0.001),
            pep: None,
            missed_cleavages: Some(0),
            psm_index: 0,
        }
    }

    fn feature(area: f64) -> FeatureQuant {
        FeatureQuant {
            apex_rt: 101.0,
            rt_start: 95.0,
            rt_stop: 107.0,
            area,
            apex_intensity: 1000.0,
            isotope_areas: vec![area * 0.6, area * 0.4],
            cosine: 0.98,
            n_isotopes_found: 2,
            rt_delta_s: 1.0,
            n_scans: 13,
        }
    }

    #[test]
    fn tmt_tsv_has_one_column_per_channel_and_raw_block() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("tmt.tsv");
        let rows = vec![TmtRow {
            id: id("run1", 0, "PEPTIDEK", 2),
            quant_scan_id: "scan=5".into(),
            quant_ms_level: 2,
            purity: Some(0.9),
            raw: vec![1.0, 2.0, 3.0, 4.0],
            corrected: vec![1.1, 2.1, 3.1, 4.1],
            sn: Some(vec![11.0, 21.0, 31.0, 41.0]),
            quant_weight: 0.81,
        }];
        write_tmt_tsv(&p, Plex::Itraq4, &rows, true).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(header.len(), 13 + 8 + 3 + 4);
        assert_eq!(header[13], "itraq_114");
        assert_eq!(header[17], "raw_itraq_114");
        assert_eq!(header[21..24], ["reporter_intensity_sum", "reporter_sn_sum", "quant_weight"]);
        assert_eq!(header[24], "sn_itraq_114");
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(row[6], "PEPTIDEK");
        assert_eq!(row[12], "0.9000");
        assert_eq!(row[13], "1.100000e0");
        assert_eq!(row[17], "1.000000e0");
        assert_eq!(row[21], "1.000000e1");
        assert_eq!(row[22], "104.00");
        assert_eq!(row[23], "0.8100");
        assert_eq!(row[27], "41.00");
        // Without noise the S/N cells stay empty and the layout is unchanged.
        let mut no_noise = rows[0].clone();
        no_noise.sn = None;
        write_tmt_tsv(&p, Plex::Itraq4, &[no_noise], false).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let mut lines = text.lines();
        assert_eq!(lines.next().unwrap().split('\t').count(), 13 + 4 + 3 + 4);
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(row[18], "");
        assert_eq!(row[20..24], ["", "", "", ""]);
    }

    /// Assert every non-empty cell outside `text_cols` parses as `f64`
    /// (`;`-joined lists cell by cell).
    fn assert_numeric_cells_parse(path: &Path, text_cols: &[&str]) {
        let text = std::fs::read_to_string(path).unwrap();
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let mut n = 0;
        for line in lines {
            for (name, cell) in header.iter().zip(line.split('\t')) {
                if text_cols.contains(name) || cell.is_empty() {
                    continue;
                }
                for v in cell.split(';') {
                    assert!(
                        v.parse::<f64>().is_ok(),
                        "{}: column {name} holds {v:?}",
                        path.display()
                    );
                    n += 1;
                }
            }
        }
        assert!(n > 0, "no numeric cells in {}", path.display());
    }

    /// Values in [1, 10) have exponent 0; every numeric cell of every quant
    /// table must still be a valid float.
    #[test]
    fn every_numeric_cell_parses_as_a_float() {
        let dir = tempfile::tempdir().unwrap();
        let tmt = dir.path().join("tmt.tsv");
        let tmt_rows = vec![TmtRow {
            id: id("run1", 0, "PEPTIDEK", 2),
            quant_scan_id: "scan=5".into(),
            quant_ms_level: 2,
            purity: Some(0.9),
            raw: vec![3.4, 0.0, 12.5, 4.2e6],
            corrected: vec![3.3, 0.0, 9.99, 4.1e6],
            sn: Some(vec![0.5, 0.0, 2.0, 1.0e4]),
            quant_weight: 1.0,
        }];
        write_tmt_tsv(&tmt, Plex::Itraq4, &tmt_rows, true).unwrap();
        assert_numeric_cells_parse(
            &tmt,
            &["filename", "spec_id", "quant_scan", "peptide", "proteins"],
        );

        let mut f = feature(3.4);
        f.apex_intensity = 5.5;
        let lfq_rows = vec![LfqRow {
            id: id("a", 0, "PEPTIDEK", 2),
            feature: f,
            score: 0.9,
            feature_q_value: Some(0.001),
            decoy_score: Some(0.1),
        }];
        let lfq = dir.path().join("lfq.tsv");
        write_lfq_tsv(&lfq, &["a".into()], &lfq_rows).unwrap();
        assert_numeric_cells_parse(&lfq, &["peptide", "proteins"]);
        let features = dir.path().join("lfq_features.tsv");
        write_lfq_features_tsv(&features, &lfq_rows).unwrap();
        assert_numeric_cells_parse(&features, &["filename", "peptide", "proteins", "spec_id"]);
    }

    #[test]
    fn lfq_tsv_is_wide_over_runs() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("lfq.tsv");
        let rows = vec![
            LfqRow {
                id: id("a", 0, "PEPTIDEK", 2),
                feature: feature(1e6),
                score: 0.9,
                feature_q_value: Some(0.001),
                decoy_score: None,
            },
            LfqRow {
                id: id("b", 1, "PEPTIDEK", 2),
                feature: feature(2e6),
                score: 0.8,
                feature_q_value: Some(0.002),
                decoy_score: Some(0.1),
            },
            LfqRow {
                id: id("b", 1, "OTHERK", 3),
                feature: feature(5e5),
                score: 0.5,
                feature_q_value: None,
                decoy_score: Some(0.6),
            },
        ];
        write_lfq_tsv(&p, &["a".into(), "b".into()], &rows).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "peptide\tcharge\tproteins\tfeature_q_value\tscore\tcosine\ta\tb"
        );
        assert_eq!(lines.len(), 3);
        let other: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(other[0], "OTHERK");
        assert_eq!(other[6], "");
        let pep: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(pep[3], "0.00100");
        assert!(pep[6].starts_with("1.000000e"), "{}", pep[6]);
        assert!(pep[7].starts_with("2.000000e"), "{}", pep[7]);
        let fp = dir.path().join("lfq_features.tsv");
        write_lfq_features_tsv(&fp, &rows).unwrap();
        assert_eq!(std::fs::read_to_string(&fp).unwrap().lines().count(), 4);
    }
}
