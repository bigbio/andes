use crate::pin::RowMeta;
use crate::text::Texts;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn non_finite(v: f64) -> Option<String> {
    if v.is_nan() {
        Some(if v.is_sign_negative() { "-nan" } else { "nan" }.to_string())
    } else if v.is_infinite() {
        Some(if v < 0.0 { "-inf" } else { "inf" }.to_string())
    } else {
        None
    }
}

fn strip_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

pub fn fmt_g(v: f64, precision: usize) -> String {
    if let Some(s) = non_finite(v) {
        return s;
    }
    let p = precision.max(1);
    if v == 0.0 {
        return if v.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    let e = format!("{:.*e}", p - 1, v);
    let (mantissa, exp) = e.split_once('e').unwrap_or((e.as_str(), "0"));
    let x: i32 = exp.parse().unwrap_or(0);
    if x < -4 || x >= p as i32 {
        let sign = if x < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", strip_zeros(mantissa), sign, x.unsigned_abs())
    } else {
        let fixed = format!("{:.*}", (p as i32 - 1 - x) as usize, v);
        strip_zeros(&fixed).to_string()
    }
}

pub fn fmt_fixed(v: f64, precision: usize) -> String {
    match non_finite(v) {
        Some(s) => s,
        None => format!("{:.*}", precision, v),
    }
}

fn io_err(path: &Path, e: std::io::Error) -> String {
    format!("write {}: {e}", path.display())
}

#[allow(clippy::too_many_arguments)]
fn write_one(
    path: &Path,
    label: i8,
    meta: &[RowMeta],
    texts: &Texts,
    file_names: &[String],
    order: &[usize],
    rows: &[usize],
    cols: [&[f64]; 3],
) -> Result<(), String> {
    let file = File::create(path).map_err(|e| io_err(path, e))?;
    let mut out = BufWriter::with_capacity(1 << 20, file);
    let mut header = String::from("PSMId\t");
    if !file_names.is_empty() {
        header.push_str("filename\t");
    }
    header.push_str("score\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n");
    out.write_all(header.as_bytes())
        .map_err(|e| io_err(path, e))?;
    let mut line = String::new();
    for &k in order {
        let m = &meta[rows[k]];
        if m.label != label {
            continue;
        }
        let t = texts.get(rows[k])?;
        line.clear();
        line.push_str(&t.psm_id);
        line.push('\t');
        if !file_names.is_empty() {
            line.push_str(file_names.get(m.file as usize).map_or("", |s| s.as_str()));
            line.push('\t');
        }
        for col in cols {
            line.push_str(&fmt_g(col[k], 6));
            line.push('\t');
        }
        line.push_str(&t.peptide);
        line.push('\t');
        let mut first = true;
        for p in t.proteins.split('\t').filter(|p| !p.is_empty()) {
            if !first {
                line.push('\t');
            }
            line.push_str(p);
            first = false;
        }
        line.push('\n');
        out.write_all(line.as_bytes())
            .map_err(|e| io_err(path, e))?;
    }
    out.flush().map_err(|e| io_err(path, e))
}

pub struct Results<'a> {
    pub rows: &'a [usize],
    pub scores: &'a [f64],
    pub q: &'a [f64],
    pub pep: &'a [f64],
}

pub fn write_results(
    target: &Path,
    decoy: &Path,
    meta: &[RowMeta],
    texts: &Texts,
    file_names: &[String],
    res: &Results,
) -> Result<(), String> {
    let (rows, scores) = (res.rows, res.scores);
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        let (ma, mb) = (&meta[rows[a]], &meta[rows[b]]);
        scores[b]
            .total_cmp(&scores[a])
            .then(mb.scan.cmp(&ma.scan))
            .then(mb.exp_mass.total_cmp(&ma.exp_mass))
            .then(mb.label.cmp(&ma.label))
            .then(rows[a].cmp(&rows[b]))
    });
    let cols = [scores, res.q, res.pep];
    write_one(target, 1, meta, texts, file_names, &order, rows, cols)?;
    write_one(decoy, -1, meta, texts, file_names, &order, rows, cols)
}

pub fn write_weights(
    path: &Path,
    names: &[String],
    w_norm: &[Vec<f64>],
    w_raw: &[Vec<f64>],
) -> Result<(), String> {
    let mut s = String::new();
    s.push_str("# This file contains the weights from each cross validation bin from percolator training\n");
    s.push_str("# First line is the feature names, followed by normalized weights, and the raw weights of bin 1\n");
    s.push_str("# This is repeated for the other bins\n");
    for (f, (norm, raw)) in w_norm.iter().zip(w_raw).enumerate() {
        s.push_str(&names.join("\t"));
        s.push_str("\tm0\n");
        for (i, &v) in norm.iter().enumerate() {
            if i > 0 {
                s.push('\t');
                s.push_str(&fmt_fixed(v, 4));
            } else if f == 0 {
                s.push_str(&fmt_g(v, 3));
            } else {
                s.push_str(&fmt_fixed(v, 3));
            }
        }
        s.push('\n');
        let cells: Vec<String> = raw.iter().map(|&v| fmt_fixed(v, 4)).collect();
        s.push_str(&cells.join("\t"));
        s.push('\n');
    }
    std::fs::write(path, s).map_err(|e| io_err(path, e))
}
