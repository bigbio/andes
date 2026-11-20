//! Cross-run candidate transfer of a DIA study (`andes quant extract` on a DIA quant directory).
//!
//! In study mode each per-run `andes --dia` process stores its nominees (the pass-2 inputs:
//! precursor, predicted RT, nominating scan RT, pool `dia_llr`, top predicted fragments) as
//! `<quant-dir>/dia_pass2/nominees/<stem>.parquet`. The transfer step of run r reads every
//! run's nominee keys and gives r, besides its own nominees, every (peptide string, charge,
//! decoy) key nominated in another run but not in r, with the descriptors of the other run
//! that nominated it with the highest `dia_llr` (the earlier run in stem order at ties).
//! Targets and decoys follow the same rule, and a decoy key only ever matches the same decoy.
//!
//! Pass 2 then runs again over r's spectra on the own and the transferred nominees. A
//! transferred nominee has `dia_llr` [`crate::dia_pass2::TRANSFER_DIA_LLR`] (a flag value) and
//! no nominating scan in r, and it is left out of every fit set of the run (RT ridge, presence
//! rates, joint-model fits), so the run's own nominees alone shape the calibration and the
//! null. The resulting pass-2 table (`dia_pass2/xfer/<stem>.parquet`, column `xfer` 1 for a
//! transferred row) replaces the run's own table in `andes quant finalize`.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use model::scan::Ms1Scan;
use model::Spectrum;
use output::dia_nominee_table::{
    read_nominee_keys, read_nominee_meta, read_nominees, write_nominees, NomineeRow,
    NomineeTable,
};
use output::dia_pass2_table::{read_pass2_meta, write_pass2_with};
use quant::Ms1RunIndex;
use rayon::prelude::*;

use crate::dia_pass2::{Nominee, TRANSFER_DIA_LLR};
use crate::dia_xrun::{precursor_key, PASS2_DIR};

/// Subdirectory of `dia_pass2/` with each run's nominee table.
const NOMINEE_DIR: &str = "nominees";
/// Subdirectory of `dia_pass2/` with each run's pass-2 table over own and transferred nominees.
const XFER_DIR: &str = "xfer";
/// Nominee-table metadata: the run's spectrum file as searched.
const SPECTRUM_KEY: &str = "andes.dia.spectrum";
/// Nominee-table metadata: the run's learned fragment m/z shift (ppm).
const SHIFT_KEY: &str = "andes.dia.fragment_shift_ppm";
/// Transfer-table metadata: the study's run stems (one per line) the transfer was built from.
pub(crate) const RUNS_KEY: &str = "andes.dia.xfer.runs";

pub(crate) fn nominee_path(quant_dir: &Path, stem: &str) -> PathBuf {
    quant_dir
        .join(PASS2_DIR)
        .join(NOMINEE_DIR)
        .join(format!("{stem}.parquet"))
}

pub(crate) fn xfer_path(quant_dir: &Path, stem: &str) -> PathBuf {
    quant_dir
        .join(PASS2_DIR)
        .join(XFER_DIR)
        .join(format!("{stem}.parquet"))
}

/// Stores the nominees of a run with its spectrum file and fragment m/z shift.
pub(crate) fn write_nominees_table(
    path: &Path,
    noms: &[Nominee],
    spectrum: &Path,
    shift_ppm: f64,
) -> Result<(), String> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
    }
    let spectrum = std::fs::canonicalize(spectrum).unwrap_or_else(|_| spectrum.to_path_buf());
    let meta = [
        (SPECTRUM_KEY, spectrum.display().to_string()),
        (SHIFT_KEY, format!("{shift_ppm}")),
    ];
    write_nominees(
        path,
        noms.len(),
        |i| {
            let n = &noms[i];
            NomineeRow {
                spec_id: &n.spec_id,
                decoy: n.decoy,
                charge: n.charge,
                peptide: &n.peptide,
                proteins: &n.proteins,
                seq: &n.seq,
                prec_mz: n.prec_mz,
                pred_rt: n.pred_rt,
                apex_rt: n.apex_rt,
                dia_llr: n.dia_llr,
                mobility: n.mobility,
                ions: &n.ions,
            }
        },
        &meta,
    )
}

/// The nominees of a table as run 0 of the extraction. Transferred rows (`source` = the stem
/// of the run they come from) get the flag `dia_llr`, an unknown nominating RT and a SpecId
/// that names their source run.
fn nominees_of(t: &NomineeTable, source: Option<&str>) -> Vec<Nominee> {
    (0..t.len())
        .map(|i| Nominee {
            run: 0,
            charge: t.charge[i],
            decoy: t.decoy[i],
            prec_mz: t.prec_mz[i],
            pred_rt: t.pred_rt[i],
            rt_missing: false,
            apex_rt: if source.is_some() { f64::NAN } else { t.apex_rt[i] },
            dia_llr: if source.is_some() { TRANSFER_DIA_LLR } else { t.dia_llr[i] },
            mobility: t.mobility[i],
            ions: t.ions(i).collect(),
            seq: t.seq[i].as_str().into(),
            peptide: t.peptide[i].as_str().into(),
            proteins: t.proteins[i].as_str().into(),
            spec_id: match source {
                Some(s) => format!("xfer_{s}_{}", t.spec_id[i]).into(),
                None => t.spec_id[i].as_str().into(),
            },
            xfer: source.is_some(),
        })
        .collect()
}

/// Best source of every key nominated in another run but not in the receiving run: the run
/// with the highest `dia_llr`, the earlier source at ties (sources are offered in study order),
/// the earlier row within a source. `best` maps key -> (source, row, dia_llr).
fn offer(
    best: &mut HashMap<u128, (u32, u32, f64)>,
    own: &HashSet<u128>,
    source: u32,
    keys: &[u128],
    llr: &[f64],
) {
    for (row, (&k, &l)) in keys.iter().zip(llr).enumerate() {
        if own.contains(&k) {
            continue;
        }
        match best.get_mut(&k) {
            Some(b) if l > b.2 => *b = (source, row as u32, l),
            Some(_) => {}
            None => {
                best.insert(k, (source, row as u32, l));
            }
        }
    }
}

/// The chosen rows of each source run, sorted.
fn rows_by_source(best: HashMap<u128, (u32, u32, f64)>, n_runs: usize) -> Vec<Vec<u32>> {
    let mut rows = vec![Vec::new(); n_runs];
    for (s, r, _) in best.into_values() {
        rows[s as usize].push(r);
    }
    for r in &mut rows {
        r.sort_unstable();
    }
    rows
}

fn keys_of(peptide: &[String], charge: &[u8], decoy: &[bool]) -> Vec<u128> {
    (0..peptide.len())
        .into_par_iter()
        .map(|i| precursor_key(&peptide[i], charge[i], decoy[i]))
        .collect()
}

/// The MS2 spectra and MS1 scans of one DIA run, read as the search read them.
fn read_run_spectra(path: &Path) -> Result<(Vec<Spectrum>, Vec<Ms1Scan>), String> {
    let (is_mzml, is_raw, is_d, _) = crate::spectra::input_format_flags(path);
    let mut errors = 0usize;
    let mut spectra = Vec::new();
    let mut take = |r: Result<Spectrum, String>| match r {
        Ok(s) => spectra.push(s),
        Err(_) => errors += 1,
    };
    let mut ms1 = if is_mzml {
        let mut reader = input::MzMLReader::new(
            input::open_buf_maybe_gz(path).map_err(|e| format!("open {}: {e}", path.display()))?,
        )
        .with_ms_level_range(2, 2)
        .with_run_ms1_capture(true);
        for r in reader.by_ref() {
            take(r.map_err(|e| e.to_string()));
        }
        reader.take_run_ms1()
    } else if is_raw {
        #[cfg(feature = "thermo")]
        {
            let mut reader = input::ThermoRawReader::open(path)
                .map_err(|e| format!("open Thermo .raw {}: {e}", path.display()))?
                .with_ms_level(Some(2))
                .with_run_ms1_capture(true);
            for r in reader.by_ref() {
                take(r.map_err(|e| e.to_string()));
            }
            reader.take_run_ms1()
        }
        #[cfg(not(feature = "thermo"))]
        {
            return Err(format!(
                "{}: this andes build has no Thermo .raw support; rebuild with `--features thermo`",
                path.display()
            ));
        }
    } else if is_d {
        #[cfg(feature = "timstof")]
        {
            let mut reader = input::TimsTofReader::open(path)
                .map_err(|e| format!("open Bruker .d {}: {e}", path.display()))?
                .with_run_ms1_capture(true)
                .with_ms1_mobility(true);
            for r in reader.by_ref() {
                take(r.map_err(|e| e.to_string()));
            }
            reader.take_run_ms1()
        }
        #[cfg(not(feature = "timstof"))]
        {
            return Err(format!(
                "{}: this andes build has no Bruker .d support; rebuild with `--features timstof`",
                path.display()
            ));
        }
    } else {
        return Err(format!(
            "{}: a DIA run needs MS1 survey scans (mzML, Thermo .raw or Bruker .d)",
            path.display()
        ));
    };
    if errors > 0 {
        log::warn!(
            "{}: {errors} spectra could not be read and are left out of the transfer",
            path.display()
        );
    }
    ms1.retain(|s| s.rt.is_finite());
    ms1.sort_by(|a, b| a.rt.total_cmp(&b.rt));
    Ok((spectra, ms1))
}

/// `andes quant extract --run <stem>` on a DIA study (module doc): writes the run's pass-2
/// table over its own and its transferred nominees.
pub(crate) fn extract_run(quant_dir: &Path, stem: &str) -> Result<(), String> {
    let t0 = std::time::Instant::now();
    let stems = crate::dia_xrun::study_runs(quant_dir)?;
    let me = stems.iter().position(|s| s == stem).ok_or_else(|| {
        format!(
            "run {stem:?} is not a DIA run of {} (runs: {})",
            quant_dir.display(),
            stems.join(", ")
        )
    })?;
    for s in &stems {
        let p = nominee_path(quant_dir, s);
        if !p.is_file() {
            return Err(format!(
                "{}: run {s} has no nominee table; search it again with --dia --quant-dir",
                p.display()
            ));
        }
    }
    let own_path = nominee_path(quant_dir, stem);
    let meta = read_nominee_meta(&own_path)?;
    let get = |k: &str| {
        meta.iter()
            .find(|(m, _)| m == k)
            .map(|(_, v)| v.clone())
            .ok_or_else(|| format!("{}: no {k} metadata", own_path.display()))
    };
    let spectrum = PathBuf::from(get(SPECTRUM_KEY)?);
    let shift: f64 = get(SHIFT_KEY)?
        .parse()
        .map_err(|_| format!("{}: unreadable {SHIFT_KEY}", own_path.display()))?;
    let own = read_nominees(&own_path, None)?;
    let own_keys: HashSet<u128> = keys_of(&own.peptide, &own.charge, &own.decoy)
        .into_iter()
        .collect();
    let mut noms = nominees_of(&own, None);
    drop(own);
    let mut best: HashMap<u128, (u32, u32, f64)> = HashMap::new();
    for (o, s) in stems.iter().enumerate().filter(|(o, _)| *o != me) {
        let k = read_nominee_keys(&nominee_path(quant_dir, s))?;
        let keys = keys_of(&k.peptide, &k.charge, &k.decoy);
        offer(&mut best, &own_keys, o as u32, &keys, &k.dia_llr);
    }
    drop(own_keys);
    let n_own = noms.len();
    let mut per_source: Vec<String> = Vec::new();
    for (o, rows) in rows_by_source(best, stems.len()).into_iter().enumerate() {
        if rows.is_empty() {
            continue;
        }
        let t = read_nominees(&nominee_path(quant_dir, &stems[o]), Some(&rows))?;
        per_source.push(format!("{} from {}", t.len(), stems[o]));
        noms.extend(nominees_of(&t, Some(&stems[o])));
    }
    let count = |r: std::ops::Range<usize>, decoy: bool| {
        noms[r].iter().filter(|n| n.decoy == decoy).count()
    };
    log::info!(
        target: "dia",
        "run {stem}: {n_own} own nominees ({} targets, {} decoys), {} transferred ({} targets, \
         {} decoys; {}) [{:.2}s]",
        count(0..n_own, false),
        count(0..n_own, true),
        noms.len() - n_own,
        count(n_own..noms.len(), false),
        count(n_own..noms.len(), true),
        if per_source.is_empty() { "none".to_string() } else { per_source.join(", ") },
        t0.elapsed().as_secs_f64()
    );
    let t_read = std::time::Instant::now();
    let (spectra, ms1) = read_run_spectra(&spectrum)?;
    if ms1.is_empty() {
        return Err(format!("{}: no MS1 survey scans", spectrum.display()));
    }
    log::info!(
        target: "dia",
        "run {stem}: read {} MS2 spectra and {} MS1 scans from {} [{:.2}s]",
        spectra.len(),
        ms1.len(),
        spectrum.display(),
        t_read.elapsed().as_secs_f64()
    );
    let t_p2 = std::time::Instant::now();
    let idx = vec![Some(Ms1RunIndex::new(ms1))];
    let p2 = crate::dia_pass2::extract(&noms, &spectra, &[spectra.len()], &idx, &[shift]);
    drop(idx);
    drop(spectra);
    let table = crate::dia_pass2::to_table(&noms, &p2);
    drop(p2);
    let xfer_kept = table
        .feature(crate::dia_pass2::XFER_COLUMN)
        .map_or(0, |v| v.iter().filter(|&&x| x > 0.5).count());
    let mut out_meta: Vec<(&str, String)> = read_pass2_meta(
        &quant_dir.join(PASS2_DIR).join(format!("{stem}.parquet")),
    )?
    .into_iter()
    .filter(|(k, _)| k == crate::dia_xrun::RESCORE_KEY)
    .map(|(_, v)| (crate::dia_xrun::RESCORE_KEY, v))
    .collect();
    out_meta.push((RUNS_KEY, stems.join("\n")));
    let path = xfer_path(quant_dir, stem);
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
    }
    write_pass2_with(&path, &table, &out_meta)?;
    log::info!(
        target: "dia",
        "run {stem}: pass 2 over {} nominees, {} inside a window ({xfer_kept} transferred), \
         wrote {} [{:.2}s; total {:.2}s]",
        noms.len(),
        table.len(),
        path.display(),
        t_p2.elapsed().as_secs_f64(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(p: &str, z: u8, d: bool) -> u128 {
        precursor_key(p, z, d)
    }

    #[test]
    fn transfer_matches_keys_across_runs_and_takes_the_best_source() {
        let own: HashSet<u128> = [key("PEPTIDEK", 2, false)].into_iter().collect();
        let mut best = HashMap::new();
        let src1 = [key("PEPTIDEK", 2, false), key("PEPTIDEK", 3, false), key("AAAK", 2, false)];
        offer(&mut best, &own, 1, &src1, &[50.0, 5.0, 7.0]);
        let src2 = [key("PEPTIDEK", 3, false), key("AAAK", 2, false), key("PEPTIDEK", 2, true)];
        offer(&mut best, &own, 2, &src2, &[9.0, 7.0, -3.0]);
        assert!(!best.contains_key(&key("PEPTIDEK", 2, false)));
        assert_eq!(best[&key("PEPTIDEK", 3, false)], (2, 0, 9.0));
        assert_eq!(best[&key("AAAK", 2, false)], (1, 2, 7.0));
        assert_eq!(best[&key("PEPTIDEK", 2, true)], (2, 2, -3.0));
        let rows = rows_by_source(best, 3);
        assert_eq!(rows, vec![vec![], vec![2], vec![0, 2]]);
    }

    #[test]
    fn transfer_treats_decoys_by_the_target_rule() {
        let own: HashSet<u128> = [key("PEPTIDEK", 2, false), key("EDITPEPK", 2, true)]
            .into_iter()
            .collect();
        let mut best = HashMap::new();
        let other = [
            key("PEPTIDEK", 2, true),
            key("EDITPEPK", 2, false),
            key("PEPTIDEK", 2, false),
            key("EDITPEPK", 2, true),
            key("KPEPTIDE", 3, false),
            key("EDIKTPEP", 3, true),
        ];
        offer(&mut best, &own, 1, &other, &[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        let picked: Vec<u32> = rows_by_source(best, 2).remove(1);
        assert_eq!(picked, vec![0, 1, 4, 5]);
    }

    #[test]
    fn transferred_nominees_are_flagged_and_keep_their_source() {
        let t = NomineeTable {
            spec_id: vec!["s_9_1".into()],
            decoy: vec![true],
            charge: vec![2],
            peptide: vec!["-.EDITPEPK.-".into()],
            proteins: vec!["XXX_P1".into()],
            seq: vec!["EDITPEPK".into()],
            prec_mz: vec![450.0],
            pred_rt: vec![12.0],
            apex_rt: vec![12.5],
            dia_llr: vec![30.0],
            mobility: vec![f64::NAN],
            ion_off: vec![0, 1],
            ion_mz: vec![300.0],
            ion_int: vec![1.0],
        };
        let own = nominees_of(&t, None);
        assert!(!own[0].xfer && own[0].dia_llr == 30.0 && own[0].apex_rt == 12.5);
        let x = nominees_of(&t, Some("run02"));
        assert!(x[0].xfer && x[0].decoy);
        assert_eq!(x[0].dia_llr, TRANSFER_DIA_LLR);
        assert!(x[0].apex_rt.is_nan());
        assert_eq!(&*x[0].spec_id, "xfer_run02_s_9_1");
        assert_eq!(x[0].ions, vec![(300.0, 1.0)]);
        assert!(!crate::dia_pass2::confident(&x[0]) && !crate::dia_pass2::null_member(&x[0]));
        assert!(crate::dia_pass2::null_member(&own[0]));
    }
}
