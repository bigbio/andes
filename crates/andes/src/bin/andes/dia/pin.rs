//! Reading a PIN into the rows the DIA steps index: each row's spectrum, peptidoform, charge,
//! label and score. Lines are parsed in parallel blocks; peptidoforms are numbered in the order
//! they first appear, I/L twins as one peptidoform (the first one seen).

use std::path::Path;

use model::Spectrum;
use onnx_infer::predictors::form::{parse_peptide, Form, Mass};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::dia::rt::{locate, row_charge, spectrum_bases, PinRow, PinTable, NONE};

pub(crate) fn read_pin(
    pin: &Path,
    spectra: &[Spectrum],
    new_cols: &[String],
) -> Result<PinTable, String> {
    let layout: std::cell::OnceCell<PinLayout> = std::cell::OnceCell::new();
    let bases = spectrum_bases(spectra);
    let mut forms: Vec<(Form, Mass)> = Vec::new();
    let mut form_ix: FxHashMap<Form, u32> = FxHashMap::default();
    let mut by_peptide: FxHashMap<String, u32> = FxHashMap::default();
    let mut rows = Vec::new();
    let mut unlocated = 0usize;
    let mut unencoded: Vec<String> = Vec::new();
    crate::dia::blocks::read_blocks(
        pin,
        |header| {
            let _ = layout.set(PinLayout::of(pin, header, new_cols)?);
            Ok(())
        },
        |first, lines| {
            let Some(l) = layout.get() else {
                return Err(format!("{}: empty PIN", pin.display()));
            };
            let parsed: Vec<Result<ParsedRow<'_>, String>> = lines
                .par_iter()
                .enumerate()
                .map(|(i, line)| l.parse(pin, first + i, line, &bases, spectra))
                .collect();
            let mut block: Vec<ParsedRow<'_>> = Vec::with_capacity(parsed.len());
            for r in parsed {
                block.push(r?);
            }
            let known: Vec<bool> = block
                .par_iter()
                .map(|r| by_peptide.contains_key(r.peptide))
                .collect();
            let mut fresh: Vec<&str> = Vec::new();
            let mut seen: rustc_hash::FxHashSet<&str> = rustc_hash::FxHashSet::default();
            for (r, _) in block.iter().zip(&known).filter(|(_, &k)| !k) {
                if seen.insert(r.peptide) {
                    fresh.push(r.peptide);
                }
            }
            let parsed_fresh: Vec<Option<(Form, Mass)>> =
                fresh.par_iter().map(|p| parse_peptide(p)).collect();
            for (p, fm) in fresh.iter().zip(parsed_fresh) {
                if fm.is_none() {
                    unencoded.push((*p).to_string());
                }
                let f = fm.map_or(NONE, |(f, m)| {
                    let key = Form {
                        seq: model::isomers::fold_isoleucine_text(&f.seq),
                        mods: f.mods.clone(),
                    };
                    *form_ix.entry(key).or_insert_with(|| {
                        forms.push((f, m));
                        (forms.len() - 1) as u32
                    })
                });
                by_peptide.insert((*p).to_string(), f);
            }
            unlocated += block.iter().filter(|r| r.row.spec == NONE).count();
            let by_peptide = &by_peptide;
            let located: Vec<PinRow> = block
                .into_par_iter()
                .map(|r| {
                    let mut row = r.row;
                    row.form = by_peptide.get(r.peptide).copied().unwrap_or(NONE);
                    row
                })
                .collect();
            rows.extend(located);
            Ok(())
        },
    )?;
    let Some(layout) = layout.into_inner() else {
        return Err(format!("{}: empty PIN", pin.display()));
    };
    if unlocated > 0 {
        log::warn!(
            target: "predict",
            "{unlocated} PIN rows match no searched spectrum; they are marked missing"
        );
    }
    report_unencoded(&unencoded, &rows, forms.len());
    Ok(PinTable {
        rows,
        forms,
        peptide_col: layout.c_pep,
    })
}

/// Examples of unencodable peptidoforms named in the [`report_unencoded`] warning.
const UNENCODED_EXAMPLES: usize = 3;

/// One warning per run for the peptidoforms the predictors cannot encode (a modification missing
/// from the Unimod table of `onnx_infer::rt::mods`, or a residue outside the 20), with how
/// many PIN rows they cover and a few of them.
fn report_unencoded(unencoded: &[String], rows: &[PinRow], encoded: usize) {
    static ONCE: std::sync::Once = std::sync::Once::new();
    if unencoded.is_empty() || ONCE.is_completed() {
        return;
    }
    ONCE.call_once(|| ());
    let n_rows = rows.iter().filter(|r| r.form == NONE).count();
    let mut examples: Vec<&str> = unencoded.iter().map(String::as_str).collect();
    examples.sort_unstable();
    examples.truncate(UNENCODED_EXAMPLES);
    log::warn!(
        target: "predict",
        "{} of {} peptidoforms ({n_rows} PIN rows) carry a modification or residue the \
         predictors do not encode, e.g. {}; those rows get no predicted RT or MS2 features",
        unencoded.len(),
        unencoded.len() + encoded,
        examples.join(", ")
    );
}

/// The columns of a PIN header that [`read_pin`] reads.
struct PinLayout {
    c_id: usize,
    c_label: usize,
    c_scan: usize,
    c_pep: usize,
    c_score: usize,
    c_exp_mass: Option<usize>,
    charge_cols: Vec<(usize, u8)>,
    charge_edges: (u8, u8),
    need: usize,
}

/// A parsed PIN row before its peptidoform is looked up.
struct ParsedRow<'a> {
    row: PinRow,
    peptide: &'a str,
}

impl PinLayout {
    fn of(pin: &Path, header: &str, new_cols: &[String]) -> Result<Self, String> {
        let h: Vec<&str> = header.split('\t').collect();
        let col = |name: &str| {
            h.iter()
                .position(|c| *c == name)
                .ok_or_else(|| format!("{}: PIN has no {name} column", pin.display()))
        };
        let (c_id, c_label, c_scan, c_pep) = (
            col("SpecId")?,
            col("Label")?,
            col("ScanNr")?,
            col("Peptide")?,
        );
        let c_score = col("RawScore").or_else(|_| col("RankScore"))?;
        for c in new_cols {
            if h.contains(&c.as_str()) {
                return Err(format!("{}: PIN already has a {c} column", pin.display()));
            }
        }
        let charge_cols: Vec<(usize, u8)> = h
            .iter()
            .enumerate()
            .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
            .collect();
        let charge_edges = (
            charge_cols.iter().map(|c| c.1).min().unwrap_or(0),
            charge_cols.iter().map(|c| c.1).max().unwrap_or(0),
        );
        let c_exp_mass = h.iter().position(|c| *c == "ExpMass");
        let need = charge_cols
            .iter()
            .map(|(i, _)| *i)
            .chain([c_id, c_label, c_scan, c_pep, c_score])
            .chain(c_exp_mass)
            .max()
            .unwrap_or(0);
        Ok(Self {
            c_id,
            c_label,
            c_scan,
            c_pep,
            c_score,
            c_exp_mass,
            charge_cols,
            charge_edges,
            need,
        })
    }

    /// Parses data line `row` (0-based, after the header) of the PIN.
    fn parse<'a>(
        &self,
        pin: &Path,
        row: usize,
        line: &'a str,
        bases: &FxHashMap<String, u32>,
        spectra: &[Spectrum],
    ) -> Result<ParsedRow<'a>, String> {
        let x: Vec<&str> = line.split('\t').collect();
        if x.len() <= self.need {
            return Err(format!(
                "{}: line {} has {} columns, expected more than {}",
                pin.display(),
                row + 2,
                x.len(),
                self.need
            ));
        }
        let scan: i32 = x[self.c_scan].trim().parse().map_err(|_| {
            format!(
                "{}: line {} ScanNr {:?} is not an integer",
                pin.display(),
                row + 2,
                x[self.c_scan]
            )
        })?;
        let spec = locate(bases, spectra, x[self.c_id], scan);
        let hot = self
            .charge_cols
            .iter()
            .find(|(i, _)| x[*i] == "1")
            .map(|(_, z)| *z);
        let charge = row_charge(
            hot,
            self.charge_edges,
            self.c_exp_mass.and_then(|i| x[i].trim().parse().ok()),
            (spec != NONE).then(|| spectra[spec as usize].precursor_mz),
        );
        Ok(ParsedRow {
            row: PinRow {
                spec,
                form: NONE,
                charge,
                decoy: x[self.c_label].trim() == "-1",
                score: x[self.c_score].trim().parse().unwrap_or(f64::NEG_INFINITY),
            },
            peptide: x[self.c_pep],
        })
    }
}
