//! Reading a PIN into the rows the DIA steps index: each row's spectrum, peptidoform, charge,
//! label and score. Lines are parsed in parallel blocks; peptidoforms are numbered in the order
//! they first appear.

use std::collections::HashMap;
use std::path::Path;

use model::Spectrum;
use onnx_infer::predictors::form::{parse_peptide, Form, Mass};
use rayon::prelude::*;

use crate::dia::rt::{locate, row_charge, spectrum_bases, PinRow, PinTable, NONE};

pub(crate) fn read_pin(
    pin: &Path,
    spectra: &[Spectrum],
    new_cols: &[String],
) -> Result<PinTable, String> {
    let layout: std::cell::OnceCell<PinLayout> = std::cell::OnceCell::new();
    let bases = spectrum_bases(spectra);
    let mut forms: Vec<(Form, Mass)> = Vec::new();
    let mut form_ix: HashMap<Form, u32> = HashMap::new();
    let mut by_peptide: HashMap<String, u32> = HashMap::new();
    let mut rows = Vec::new();
    let mut unlocated = 0usize;
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
            let mut fresh: Vec<&str> = Vec::new();
            let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
            for r in &block {
                if !by_peptide.contains_key(r.peptide) && seen.insert(r.peptide) {
                    fresh.push(r.peptide);
                }
            }
            let parsed_fresh: Vec<Option<(Form, Mass)>> =
                fresh.par_iter().map(|p| parse_peptide(p)).collect();
            for (p, fm) in fresh.iter().zip(parsed_fresh) {
                let f = fm.map_or(NONE, |(f, m)| {
                    *form_ix.entry(f.clone()).or_insert_with(|| {
                        forms.push((f, m));
                        (forms.len() - 1) as u32
                    })
                });
                by_peptide.insert((*p).to_string(), f);
            }
            rows.reserve(block.len());
            for r in block {
                unlocated += usize::from(r.row.spec == NONE);
                let mut row = r.row;
                row.form = by_peptide.get(r.peptide).copied().unwrap_or(NONE);
                rows.push(row);
            }
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
    Ok(PinTable {
        rows,
        forms,
        peptide_col: layout.c_pep,
    })
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
        bases: &HashMap<String, u32>,
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
