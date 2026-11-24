//! Training inputs: the spectra and labelled PSMs of each run, and the search parameters used to
//! label them.

use std::path::Path;

use input::{MgfReader, MzMLReader};
use model::Spectrum;
use search::SearchParams;

use super::args::TrainFromSearchArgs;
use crate::models::build_aa_set;
use crate::spectra::spectrum_ext_lower;

/// Spectra a reader returned as errors: how many, and the first error text.
#[derive(Default)]
pub(super) struct ParseFailures {
    pub(super) count: usize,
    pub(super) first: Option<String>,
}

impl ParseFailures {
    pub(super) fn add(&mut self, e: impl std::fmt::Display) {
        self.count += 1;
        if self.first.is_none() {
            self.first = Some(e.to_string());
        }
    }

    /// Adds spectra a reader skipped on its own, with the errors it kept.
    fn add_skipped(&mut self, count: usize, errors: &[String]) {
        self.count += count;
        if self.first.is_none() {
            self.first = errors.first().cloned();
        }
    }

    /// One warning for all failed spectra of `path`, none when every spectrum parsed.
    fn report(&self, path: &Path) {
        if self.count == 0 {
            return;
        }
        let first = self.first.as_deref().unwrap_or("no message");
        log::warn!(
            target: "train",
            "{} spectra of {} could not be parsed and were skipped (first: {first})",
            self.count,
            path.display()
        );
    }
}

pub(crate) fn load_spectra_for_train(
    path: &Path,
) -> Result<Vec<Spectrum>, Box<dyn std::error::Error>> {
    let ext_lower = spectrum_ext_lower(path);
    let mut spectra = Vec::new();
    let mut failed = ParseFailures::default();
    match ext_lower.as_deref() {
        Some("mzml") => {
            let mut reader =
                MzMLReader::new(input::open_buf_maybe_gz(path)?).with_ms_level_range(2, 2);
            for item in reader.by_ref() {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => failed.add(e),
                }
            }
            failed.add_skipped(reader.skipped_count(), reader.skipped_errors());
        }
        Some("raw") => {
            #[cfg(feature = "thermo")]
            {
                let reader = input::ThermoRawReader::open(path)
                    .map_err(|e| format!("open Thermo .raw {}: {e}", path.display()))?
                    .with_ms_level(Some(2));
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => failed.add(e),
                    }
                }
            }
            #[cfg(not(feature = "thermo"))]
            {
                return Err(format!(
                    "native Thermo `.raw` training input requires building with \
                     `--features thermo` (and the .NET 8 runtime): {}",
                    path.display()
                )
                .into());
            }
        }
        Some("d") => {
            #[cfg(feature = "timstof")]
            {
                let reader = input::TimsTofReader::open(path)
                    .map_err(|e| format!("open Bruker .d {}: {e}", path.display()))?;
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => failed.add(e),
                    }
                }
            }
            #[cfg(not(feature = "timstof"))]
            {
                return Err(format!(
                    "native Bruker `.d` training input requires building with \
                     `--features timstof`: {}",
                    path.display()
                )
                .into());
            }
        }
        _ => {
            let reader = MgfReader::new(input::open_buf_maybe_gz(path)?);
            for item in reader {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => failed.add(e),
                }
            }
        }
    }
    failed.report(path);
    Ok(spectra)
}

pub(crate) fn build_train_search_params(
    args: &TrainFromSearchArgs,
) -> Result<SearchParams, Box<dyn std::error::Error>> {
    let aa = build_aa_set(args.fixed_mods.as_ref(), args.variable_mods.as_ref())?;
    let mut params = SearchParams::default_tryptic(aa);
    params.charge_range = 2..=5;
    Ok(params)
}

pub(crate) fn first_nglyco_site(seq: &[u8]) -> Option<usize> {
    (0..seq.len()).find(|&i| {
        seq[i] == b'N'
            && i + 2 < seq.len()
            && seq[i + 1] != b'P'
            && (seq[i + 2] == b'S' || seq[i + 2] == b'T')
    })
}

pub(crate) fn load_labels_from_tsv(
    path: &std::path::Path,
    spectra: &[model::spectrum::Spectrum],
    aa_set: &model::amino_acid::set::AminoAcidSet,
) -> Result<Vec<model_train::labeled::LabeledMatch>, Box<dyn std::error::Error>> {
    use model::peptide::Peptide;
    use std::collections::HashMap;

    let mut scan_to_idx: HashMap<i32, usize> = HashMap::new();
    for (i, s) in spectra.iter().enumerate() {
        if let Some(sc) = s.scan {
            scan_to_idx.entry(sc).or_insert(i);
        }
    }

    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("reading labels {}: {e}", path.display()))?;
    let mut lines = text.lines();
    let header = lines.next().ok_or("labels file is empty")?;
    let cols: Vec<String> = header
        .split('\t')
        .map(|c| c.trim().to_ascii_lowercase())
        .collect();
    let col = |names: &[&str]| -> Option<usize> {
        cols.iter().position(|c| names.iter().any(|n| c == n))
    };
    let scan_c = col(&["scan", "scannr", "scan_number"]).ok_or("labels: no scan column")?;
    let pep_c = col(&["peptide", "backbone", "sequence", "peptide_sequence"])
        .ok_or("labels: no peptide column")?;
    let chg_c = col(&["charge", "z", "precursor_charge"]).ok_or("labels: no charge column")?;
    let gly_c = col(&["glycan_mass", "glycan", "glycan_neutral"]);
    let site_c = col(&["glycosite", "site", "glyco_site"]);
    let max_c = scan_c.max(pep_c).max(chg_c);

    let anywhere: std::collections::HashMap<u8, f64> =
        aa_set.fixed_anywhere_deltas().into_iter().collect();
    let terminal = aa_set.fixed_terminal_deltas();
    let at = |loc: model::modification::ModLocation| -> f64 {
        terminal
            .iter()
            .filter(|(l, _)| *l == loc)
            .map(|(_, d)| d)
            .sum()
    };
    let (n_term, c_term) = (
        at(model::modification::ModLocation::NTerm),
        at(model::modification::ModLocation::CTerm),
    );
    let decorate = |seq: &str| -> String {
        let mut d = String::with_capacity(seq.len() + 4);
        d.push_str("-.");
        let n = seq.len();
        for (i, &b) in seq.as_bytes().iter().enumerate() {
            d.push(b as char);
            let mut delta = anywhere.get(&b).copied().unwrap_or(0.0);
            if i == 0 {
                delta += n_term;
            }
            if i + 1 == n {
                delta += c_term;
            }
            if delta != 0.0 {
                d.push_str(&format!("+{:.5}", delta));
            }
        }
        d.push_str(".-");
        d
    };

    let mut labels = Vec::new();
    let mut seen_scans: std::collections::HashSet<i32> = std::collections::HashSet::new();
    let (mut miss_scan, mut miss_pep, mut miss_other) = (0usize, 0usize, 0usize);
    let (mut dup_scan, mut charge_mismatch) = (0usize, 0usize);
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= max_c {
            miss_other += 1;
            continue;
        }
        let scan: i32 = match f[scan_c].trim().parse() {
            Ok(v) => v,
            Err(_) => {
                miss_other += 1;
                continue;
            }
        };
        let charge: u8 = match f[chg_c].trim().parse() {
            Ok(v) => v,
            Err(_) => {
                miss_other += 1;
                continue;
            }
        };
        let idx = match scan_to_idx.get(&scan) {
            Some(&v) => v,
            None => {
                miss_scan += 1;
                continue;
            }
        };
        let seq = f[pep_c].trim();
        let glyco: Option<(usize, f64)> = gly_c.and_then(|gc| {
            let gmass: f64 = f.get(gc)?.trim().parse().ok()?;
            if gmass <= 0.0 {
                return None;
            }
            let site = site_c
                .and_then(|sc| f.get(sc)?.trim().parse::<usize>().ok())
                .or_else(|| first_nglyco_site(seq.as_bytes()))?;
            (site < seq.len()).then_some((site, gmass))
        });
        let mut peptide = match Peptide::from_str(&decorate(seq), aa_set) {
            Ok(p) => p,
            Err(_) => {
                miss_pep += 1;
                continue;
            }
        };
        if let Some((site, gmass)) = glyco {
            if let Some(res) = peptide.residues.get_mut(site) {
                let base = res.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
                res.mod_ = Some(model::modification::leak_mod(
                    model::modification::Modification {
                        name: "Glycan".to_string(),
                        mass_delta: base + gmass,
                        residue: model::modification::ResidueSpec::Specific(res.residue),
                        location: model::modification::ModLocation::Anywhere,
                        fixed: false,
                        accession: None,
                        neutral_losses: Vec::new(),
                        loss_class: 1,
                        annotated_site_mass: None,
                        blocks_cleavage: false,
                        parts: Vec::new(),
                        diagnostic_ions: Vec::new(),
                    },
                ));
                let charge = peptide.charge;
                let residues = std::mem::take(&mut peptide.residues);
                peptide = Peptide::new(residues, peptide.pre, peptide.post);
                peptide.charge = charge;
            }
        }
        if !seen_scans.insert(scan) {
            dup_scan += 1;
            continue;
        }
        if let Some(spec_z) = spectra[idx].precursor_charge {
            if spec_z != charge as i32 {
                charge_mismatch += 1;
            }
        }
        labels.push(model_train::labeled::LabeledMatch {
            spectrum_index: idx,
            peptide,
            charge,
            confidence: 0.0,
        });
    }
    log::info!(
        "loaded {} labels from {} ({} skipped: {} scan-not-found, {} dup-scan, {} unparseable-peptide, {} malformed; {} charge-mismatch vs mzML)",
        labels.len(),
        path.display(),
        miss_scan + dup_scan + miss_pep + miss_other,
        miss_scan,
        dup_scan,
        miss_pep,
        miss_other,
        charge_mismatch,
    );
    if charge_mismatch * 5 > labels.len().max(1) {
        log::warn!(
            "{charge_mismatch} labels ({}%) disagree with the mzML precursor \
             charge; the scan->spectrum mapping or the label charges may be wrong.",
            charge_mismatch * 100 / labels.len().max(1),
        );
    }
    Ok(labels)
}

pub(super) struct TrainRun {
    pub(super) dataset: String,
    pub(super) stem: String,
    pub(super) spectra: Vec<Spectrum>,
    pub(super) labels: Vec<model_train::labeled::LabeledMatch>,
}

impl TrainRun {
    pub(super) fn keep_labelled(
        path: &Path,
        spectra: Vec<Spectrum>,
        mut labels: Vec<model_train::labeled::LabeledMatch>,
    ) -> Self {
        let mut kept = Vec::with_capacity(labels.len());
        let mut new_index = std::collections::HashMap::new();
        for label in &mut labels {
            let idx = *new_index.entry(label.spectrum_index).or_insert_with(|| {
                kept.push(spectra[label.spectrum_index].clone());
                kept.len() - 1
            });
            label.spectrum_index = idx;
        }
        let dataset = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());
        let stem = run_stem(path);
        TrainRun {
            dataset,
            stem,
            spectra: kept,
            labels,
        }
    }
}

pub(super) fn run_stem(path: &Path) -> String {
    let mut name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string());
    for ext in [".gz", ".mzml", ".mgf", ".d", ".raw"] {
        if name.len() > ext.len() && name.to_ascii_lowercase().ends_with(ext) {
            name.truncate(name.len() - ext.len());
        }
    }
    name
}
