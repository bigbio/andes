use crate::budget::Reservation;
use crate::pin::RowMeta;
use crate::scores::psm_order;
use crate::text::{RowText, Texts};
use rayon::prelude::*;
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
        None => format!("{v:.precision$}"),
    }
}

fn io_err(path: &Path, e: std::io::Error) -> String {
    format!("write {}: {e}", path.display())
}

const LINES_PER_TASK: usize = 4096;

const LINES_PER_WINDOW: usize = 16 * LINES_PER_TASK;

struct Lines<'a> {
    meta: &'a [RowMeta],
    file_names: &'a [String],
    rows: &'a [usize],
    cols: [&'a [f64]; 3],
}

impl Lines<'_> {
    fn push(&self, k: usize, t: &RowText, line: &mut String) {
        let m = &self.meta[self.rows[k]];
        line.push_str(&t.psm_id);
        line.push('\t');
        if !self.file_names.is_empty() {
            line.push_str(
                self.file_names
                    .get(m.file as usize)
                    .map_or("", |s| s.as_str()),
            );
            line.push('\t');
        }
        for col in self.cols {
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
    }
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
    let lines = Lines {
        meta,
        file_names,
        rows,
        cols,
    };
    let spilled = texts.is_spilled();
    let (window_lines, _held) = if spilled {
        texts.batch_rows(LINES_PER_WINDOW)
    } else {
        (LINES_PER_WINDOW, Reservation::none())
    };
    for window in order.chunks(window_lines) {
        let read: Vec<RowText> = if spilled {
            let ids: Vec<usize> = window
                .iter()
                .filter(|&&k| meta[rows[k]].label == label)
                .map(|&k| rows[k])
                .collect();
            texts.get_many(&ids)?
        } else {
            Vec::new()
        };
        let starts: Vec<usize> = window
            .chunks(LINES_PER_TASK)
            .scan(0usize, |n, ks| {
                let at = *n;
                *n += ks.iter().filter(|&&k| meta[rows[k]].label == label).count();
                Some(at)
            })
            .collect();
        let parts: Vec<Result<String, String>> = window
            .par_chunks(LINES_PER_TASK)
            .zip(&starts)
            .map(|(ks, &at)| {
                let mut buf = String::new();
                let mut j = at;
                for &k in ks {
                    if meta[rows[k]].label == label {
                        if spilled {
                            lines.push(k, &read[j], &mut buf);
                            j += 1;
                        } else {
                            lines.push(k, &texts.get(rows[k])?, &mut buf);
                        }
                    }
                }
                Ok(buf)
            })
            .collect();
        for part in parts {
            out.write_all(part?.as_bytes())
                .map_err(|e| io_err(path, e))?;
        }
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
    let order = psm_order(meta, res.rows, res.scores, |i| texts.get(i).map(|t| t.psm_id))?;
    write_ordered(target, decoy, meta, texts, file_names, res, &order)
}

/// `write_results` with the rows already in output order (`order` from `psm_order` on the
/// same rows and scores).
pub(crate) fn write_ordered(
    target: &Path,
    decoy: &Path,
    meta: &[RowMeta],
    texts: &Texts,
    file_names: &[String],
    res: &Results,
    order: &[usize],
) -> Result<(), String> {
    let (rows, scores) = (res.rows, res.scores);
    let cols = [scores, res.q, res.pep];
    write_one(target, 1, meta, texts, file_names, order, rows, cols)?;
    write_one(decoy, -1, meta, texts, file_names, order, rows, cols)
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
