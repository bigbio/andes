//! `quant extract`: the donor PSMs of each library entry and the per-run extraction it drives.

use std::collections::HashMap;
use std::path::Path;

use model::Tolerance;
use output::extract::write_extract_table_tagged;
use output::library::{StoredLibrary, FINGERPRINT_KEY};
use output::ms2::read_ms2_table;
use output::quantities::runs::read_run_quant_mods;
use quant::evidence::EvidenceParams;
use quant::library::{Library, PsmIn};
use quant::ms2::{fragment_ladder, library_spectrum, residue_masses, LibrarySpectrum, Ms2RunIndex};
use rayon::prelude::*;

use super::{
    artifact, EXTRACT_SUFFIX, MS1_SUFFIX, MS2_SUFFIX, RUN_SCALES_KEY, SELECTION_KEY, TABLE_SUFFIX,
    TOL_KEY,
};
use crate::quant::mbr::{
    attach_ms2, attach_shared, extract_run, ExtractCounts, ExtractRow, Ms2Counts, TransferSelection,
};

#[derive(Debug, Default, PartialEq)]
pub(crate) struct SpectrumCounts {
    pub with_spectrum: usize,
    pub no_scan: usize,
    pub no_masses: usize,
    pub too_few: usize,
}

pub(crate) type EntrySpectra = (
    Vec<Option<LibrarySpectrum>>,
    Vec<Vec<(u16, f64)>>,
    SpectrumCounts,
);

pub(crate) fn scan_of(by_key: &HashMap<&str, usize>, spec_id: &str) -> Option<usize> {
    let mut id = spec_id;
    for _ in 0..2 {
        id = id.rsplit_once('_')?.0;
        if let Some(&i) = by_key.get(id) {
            return Some(i);
        }
    }
    None
}

#[derive(Clone, Copy)]
pub(super) struct BestForm {
    pub(super) donor: bool,
    pub(super) score: f64,
    pub(super) rt: f64,
    pub(super) run: u16,
    pub(super) site: u32,
}

impl BestForm {
    pub(super) fn beats(&self, other: &Self) -> bool {
        if self.donor != other.donor {
            return self.donor;
        }
        self.score > other.score
            || (self.score == other.score
                && (self.rt.total_cmp(&other.rt).then(self.run.cmp(&other.run))).is_lt())
    }
}

pub(crate) struct DonorPsms {
    slot: Vec<u32>,
    pub(super) scans: Vec<Vec<(u16, f64)>>,
    best: Vec<(BestDonor, Box<str>)>,
}

#[derive(Clone, Copy)]
pub(super) struct BestDonor {
    pub(super) score: f64,
    pub(super) rt: f64,
    pub(super) run: u16,
}

impl DonorPsms {
    pub(crate) fn new(n_precursors: usize) -> Self {
        Self {
            slot: vec![u32::MAX; n_precursors],
            scans: Vec::new(),
            best: Vec::new(),
        }
    }

    pub(crate) fn add(&mut self, p: &PsmIn, spec_id: &str) {
        let id = p.precursor as usize;
        if id >= self.slot.len() {
            return;
        }
        let new = BestDonor {
            score: p.score,
            rt: p.rt,
            run: p.run,
        };
        let i = match self.slot[id] {
            u32::MAX => {
                self.slot[id] = u32::try_from(self.scans.len())
                    .ok()
                    .filter(|&i| i != u32::MAX)
                    .expect("fewer than 2^32 - 1 donor precursors");
                self.scans.push(Vec::new());
                self.best.push((new, spec_id.into()));
                self.scans.len() - 1
            }
            i => {
                let i = i as usize;
                let (b, s) = &mut self.best[i];
                if p.score > b.score
                    || (p.score == b.score
                        && (p.rt.total_cmp(&b.rt).then(p.run.cmp(&b.run))).is_lt())
                {
                    *b = new;
                    *s = spec_id.into();
                }
                i
            }
        };
        self.scans[i].push((p.run, p.rt));
    }

    pub(crate) fn finish(
        self,
        dir: &Path,
        stems: &[String],
        lib: &Library,
    ) -> Result<EntrySpectra, String> {
        let n = lib.meta.len();
        let Self {
            slot,
            mut scans,
            best: donor_best,
        } = self;
        let mut entry_scans: Vec<Vec<(u16, f64)>> = vec![Vec::new(); n];
        let mut best: Vec<Option<(u16, Box<str>)>> = vec![None; n];
        let mut donor_best: Vec<Option<(BestDonor, Box<str>)>> =
            donor_best.into_iter().map(Some).collect();
        for e in &lib.entries {
            let id = e.precursor as usize;
            let Some(&i) = slot.get(id).filter(|&&i| i != u32::MAX) else {
                continue;
            };
            let mut d = std::mem::take(&mut scans[i as usize]);
            d.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
            entry_scans[id] = d;
            best[id] = donor_best[i as usize].take().map(|(b, s)| (b.run, s));
        }
        drop(scans);
        drop(donor_best);
        drop(slot);
        let mut spectra: Vec<Option<LibrarySpectrum>> = vec![None; n];
        let mut c = SpectrumCounts::default();
        let mut wanted: Vec<Vec<(usize, &str)>> = vec![Vec::new(); stems.len()];
        for (id, b) in best.iter().enumerate() {
            if let Some((run, spec_id)) = b {
                if let Some(w) = wanted.get_mut(*run as usize) {
                    w.push((id, &**spec_id));
                }
            }
        }
        type Made = Vec<(usize, Result<LibrarySpectrum, u8>)>;
        let per_run: Vec<Result<Made, String>> = wanted
            .par_iter()
            .zip(stems.par_iter())
            .map(|(wanted, stem)| {
                if wanted.is_empty() {
                    return Ok(Vec::new());
                }
                let ms2_path = artifact(dir, stem, MS2_SUFFIX);
                let table = read_ms2_table(&ms2_path)?;
                if table.run != *stem {
                    return Err(format!(
                        "{}: the table is of run {:?}, not {stem:?}",
                        ms2_path.display(),
                        table.run
                    ));
                }
                let mods = read_run_quant_mods(&artifact(dir, stem, TABLE_SUFFIX))?;
                let by_key: HashMap<&str, usize> = table
                    .scans
                    .iter()
                    .enumerate()
                    .map(|(i, s)| (s.spec_key.as_str(), i))
                    .collect();
                Ok(wanted
                    .iter()
                    .map(|&(id, spec_id)| {
                        let m = &lib.meta[id];
                        let Some(scan) = scan_of(&by_key, spec_id) else {
                            return (id, Err(0));
                        };
                        let Ok(masses) = residue_masses(split_key(&m.key), &mods) else {
                            return (id, Err(1));
                        };
                        let ladder = fragment_ladder(&masses, m.charge.max(1));
                        let peaks = &table.scans[scan].peaks;
                        (
                            id,
                            library_spectrum(&ladder, peaks, table.fragment_tol).ok_or(2),
                        )
                    })
                    .collect())
            })
            .collect();
        for made in per_run {
            for (id, s) in made? {
                match s {
                    Ok(s) => {
                        c.with_spectrum += 1;
                        spectra[id] = Some(s);
                    }
                    Err(0) => c.no_scan += 1,
                    Err(1) => c.no_masses += 1,
                    Err(_) => c.too_few += 1,
                }
            }
        }
        Ok((spectra, entry_scans, c))
    }
}

fn log_extract(stem: &str, rows: usize, c: &ExtractCounts) {
    log::info!(
        target: "extract",
        "{stem}: {rows} rows{}; skipped: {} entries at a FAIMS voltage the run lacks, {} identified pairs with a non-finite best PSM RT, {} identified pairs without a prediction from other runs (no pseudo-transfer), {} transfers without a prediction, {} pseudo-transfers and {} transfers whose donor runs are all in other alignment components, {} transfers left out by the transfer selection; {} identified pairs and {} of their decoy rows measured on the larger peak at the position the other runs predict; MS1 scan spacing {}, width scale {:.3} s ({} fits, {} rounds{}); width prior from {} fits{}",
        if c.unaligned {
            " (unaligned run: identified rows only)"
        } else {
            ""
        },
        c.foreign_cv,
        c.nonfinite_best_rt,
        c.no_pseudo_prediction,
        c.no_prediction,
        c.pseudo_other_component,
        c.other_component,
        c.not_selected,
        c.consensus_peak,
        c.consensus_peak_decoys,
        c.scales
            .spacing_s
            .map_or_else(|| "unknown".to_string(), |s| format!("{s:.3} s")),
        c.scales.width_s,
        c.scales.width_fits,
        c.scales.rounds,
        if c.scales.fallback {
            "; too few fits: the seed width is used"
        } else {
            ""
        },
        c.width_fits,
        if c.width_prior_fallback {
            " (too few: the width scale is used)"
        } else {
            ""
        }
    );
}

fn log_ms2(stem: &str, c: &Ms2Counts, secs: f64) {
    let part = |role: usize| {
        ["target", "mass-shift decoy", "RT decoy"]
            .iter()
            .enumerate()
            .filter(|&(f, _)| c.rows[role][f] > 0)
            .map(|(f, name)| {
                let n = c.rows[role][f] as f64;
                format!(
                    "{name} {} (scan in window {:.3}, supported {:.4})",
                    c.rows[role][f],
                    c.with_scan[role][f] as f64 / n,
                    c.supported[role][f] as f64 / n
                )
            })
            .collect::<Vec<_>>()
            .join(", ")
    };
    log::debug!(
        target: "extract",
        "{stem}: MS2 evidence ({secs:.1} s): transfers: {}; pseudo-transfers: {}; {} rows of entries without a library spectrum",
        part(0),
        part(1),
        c.no_spectrum
    );
}

pub(crate) fn extract_and_write(
    dir: &Path,
    stored: &StoredLibrary,
    fingerprint: &str,
    run: u16,
    ep: &EvidenceParams,
    selection: &TransferSelection,
) -> Result<Vec<ExtractRow>, String> {
    let stem = &stored.runs[run as usize];
    let ms1 = artifact(dir, stem, MS1_SUFFIX);
    let mut ex = extract_run(run, &ms1, &stored.library, ep, selection)?;
    log_extract(stem, ex.rows.len(), &ex.counts);
    let t = std::time::Instant::now();
    let ms2_path = artifact(dir, stem, MS2_SUFFIX);
    let table = read_ms2_table(&ms2_path)?;
    if table.run != *stem {
        return Err(format!(
            "{}: the table is of run {:?}, not {stem:?}",
            ms2_path.display(),
            table.run
        ));
    }
    let tol = table.fragment_tol;
    let index = Ms2RunIndex::new(table.scans);
    let c = attach_ms2(
        &mut ex.rows,
        &stored.library,
        &stored.spectra,
        &stored.donor_scans,
        &index,
        tol,
    );
    log_ms2(stem, &c, t.elapsed().as_secs_f64());
    if quant::shared::Model::shipped().needs_ms2_evidence() {
        let t = std::time::Instant::now();
        let mods = read_run_quant_mods(&artifact(dir, stem, TABLE_SUFFIX))?;
        let sc = attach_shared(
            &mut ex.rows,
            &stored.library,
            &stored.donor_scans,
            &mods,
            &index,
            tol,
        );
        log::debug!(
            target: "extract",
            "{stem}: shared-feature evidence ({:.1} s): {} transfers with a claimant, {} with a covering scan, {} scored against a claimant's distinguishing ions",
            t.elapsed().as_secs_f64(),
            sc.claimed,
            sc.with_scan,
            sc.scored
        );
    }
    drop(index);
    let path = artifact(dir, stem, EXTRACT_SUFFIX);
    let mut tags = vec![
        (FINGERPRINT_KEY.to_string(), fingerprint.to_string()),
        (TOL_KEY.to_string(), ep.tol_ppm.to_string()),
        (RUN_SCALES_KEY.to_string(), ex.counts.scales.describe()),
    ];
    if !selection.is_off() {
        let aligned = stored.library.alignment.maps.iter().flatten().count();
        tags.push((SELECTION_KEY.to_string(), selection.describe(aligned)));
    }
    write_extract_table_tagged(&path, &ex.rows, &tags)?;
    Ok(ex.rows)
}

pub(crate) fn evidence_params(tol: Option<Tolerance>) -> Result<EvidenceParams, String> {
    match tol {
        None => Ok(EvidenceParams::default()),
        Some(Tolerance::Ppm(v)) => Ok(EvidenceParams {
            tol_ppm: v,
            ..EvidenceParams::default()
        }),
        Some(Tolerance::Da(v)) => Err(format!(
            "match-between-runs extraction needs a ppm --lfq-tol, got {v} Da"
        )),
    }
}

pub(crate) fn split_key(key: &str) -> &str {
    key.rsplit_once('/').map_or(key, |(p, _)| p)
}
