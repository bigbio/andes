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
//! transferred nominee has `dia_llr` [`crate::dia::pass2::TRANSFER_DIA_LLR`] (a flag value) and
//! no nominating scan in r; its predicted RT is its raw andes-rt prediction through r's own
//! andes-rt calibration (stored with r's nominees), so it sits on r's time axis exactly as r's
//! own PIN rows do, and it is left out of every fit set of the run (RT ridge, presence
//! rates, joint-model fits), so the run's own nominees alone shape the calibration and the
//! null. The resulting pass-2 table (`dia_pass2/xfer/<stem>.parquet`, column `xfer` 1 for a
//! transferred row) replaces the run's own table in `andes quant finalize`.
//!
//! The transfer table records a fingerprint of every run's stored inputs (nominee table,
//! pass-2 table, pass-1 PIN); `andes quant finalize` refuses a transfer table whose inputs
//! were stored again after it was built.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use model::hash::{fnv1a64_extend, FNV1A64_OFFSET};
use model::scan::Ms1Scan;
use model::Spectrum;
use onnx_infer::predictors::calibration::RtCalibration;
use output::dia::nominees::{
    read_nominee_keys, read_nominee_meta, read_nominees, write_nominees, NomineeRow,
    NomineeTable,
};
use output::dia::pass2::{read_pass2_meta, write_pass2_with, Pass2Table};
use quant::Ms1RunIndex;
use rayon::prelude::*;

use crate::dia::pass2::{Nominee, TRANSFER_DIA_LLR};
use crate::dia::study::{precursor_key, PASS2_DIR};

/// Subdirectory of `dia_pass2/` with each run's nominee table.
const NOMINEE_DIR: &str = "nominees";
/// Subdirectory of `dia_pass2/` with each run's pass-2 table over own and transferred nominees.
const XFER_DIR: &str = "xfer";
/// Nominee-table metadata: the run's spectrum file as searched.
const SPECTRUM_KEY: &str = "andes.dia.spectrum";
/// Nominee-table metadata: the run's learned fragment m/z shift (ppm).
const SHIFT_KEY: &str = "andes.dia.fragment_shift_ppm";
/// Nominee-table metadata: the run's andes-rt calibration (raw prediction to minutes).
const RT_CAL_KEY: &str = "andes.dia.rt_calibration";
/// Nominee-table metadata: `--dia-mobility-slices` of the run's search.
const MOBILITY_SLICES_KEY: &str = "andes.dia.mobility_slices";
/// Nominee-table and pass-2-table metadata: a content hash of the table's rows.
pub(crate) const STAMP_KEY: &str = "andes.dia.stamp";
/// Transfer-table metadata: the study's run stems (one per line) the transfer was built from.
pub(crate) const RUNS_KEY: &str = "andes.dia.xfer.runs";
/// Transfer-table metadata: the fingerprint of every run's stored inputs ([`run_inputs`]), one
/// line per run in study order.
pub(crate) const INPUTS_KEY: &str = "andes.dia.xfer.inputs";

/// Where and how a run's nominees are stored (study mode).
pub(crate) struct NomineesOut<'a> {
    pub(crate) path: &'a Path,
    /// The run's spectrum file as searched.
    pub(crate) spectrum: &'a Path,
    /// `--dia-mobility-slices` of the search, applied again when the run is read for transfer.
    pub(crate) mobility_slices: usize,
}

/// Content hash of a run's nominees, with their count.
fn nominee_stamp(noms: &[Nominee]) -> String {
    let mut h = FNV1A64_OFFSET;
    for n in noms {
        h = fnv1a64_extend(h, n.spec_id.bytes().chain(n.peptide.bytes()));
        h = fnv1a64_extend(h, [u8::from(n.decoy), n.charge, u8::from(n.rt_missing)]);
        for v in [n.prec_mz, n.pred_rt, n.raw_rt, n.apex_rt, n.dia_llr, n.mobility] {
            h = fnv1a64_extend(h, v.to_bits().to_le_bytes());
        }
        for &(mz, p) in &n.ions {
            h = fnv1a64_extend(h, mz.to_bits().to_le_bytes());
            h = fnv1a64_extend(h, p.to_bits().to_le_bytes());
        }
    }
    format!("{h:016x}/{}", noms.len())
}

/// Content hash of a pass-2 table, with its row count.
pub(crate) fn pass2_stamp(t: &Pass2Table) -> String {
    let mut h = FNV1A64_OFFSET;
    for i in 0..t.len() {
        h = fnv1a64_extend(h, t.spec_id[i].bytes().chain(t.peptide[i].bytes()));
        h = fnv1a64_extend(h, t.run[i].to_le_bytes());
        h = fnv1a64_extend(h, [u8::from(t.decoy[i]), t.charge[i]]);
        h = fnv1a64_extend(h, t.prec_mz[i].to_bits().to_le_bytes());
    }
    for (name, v) in &t.features {
        h = fnv1a64_extend(h, name.bytes());
        for x in v {
            h = fnv1a64_extend(h, x.to_bits().to_le_bytes());
        }
    }
    format!("{h:016x}/{}", t.len())
}

/// A run's andes-rt calibration as nominee-table metadata (exact: shortest round-trip floats).
fn encode_rt_cal(c: Option<&RtCalibration>) -> String {
    match c {
        None => "none".to_string(),
        Some(RtCalibration::Linear { slope, intercept }) => format!("linear {slope} {intercept}"),
        Some(RtCalibration::Isotonic { x, y }) => {
            let mut out = String::from("isotonic");
            for (a, b) in x.iter().zip(y) {
                out.push_str(&format!(" {a} {b}"));
            }
            out
        }
    }
}

fn decode_rt_cal(s: &str) -> Result<Option<RtCalibration>, String> {
    let mut it = s.split(' ');
    let kind = it.next().unwrap_or_default();
    let v: Vec<f64> = it
        .map(|x| x.parse::<f64>().map_err(|_| format!("unreadable {RT_CAL_KEY} value {x:?}")))
        .collect::<Result<_, _>>()?;
    match (kind, v.len()) {
        ("none", 0) => Ok(None),
        ("linear", 2) => Ok(Some(RtCalibration::Linear {
            slope: v[0],
            intercept: v[1],
        })),
        ("isotonic", n) if n >= 2 && n % 2 == 0 => Ok(Some(RtCalibration::Isotonic {
            x: v.iter().step_by(2).copied().collect(),
            y: v.iter().skip(1).step_by(2).copied().collect(),
        })),
        _ => Err(format!("unreadable {RT_CAL_KEY} metadata {kind:?} with {} values", v.len())),
    }
}

/// How a run is read again for transfer, from its nominee-table metadata.
struct RunRead {
    spectrum: PathBuf,
    shift_ppm: f64,
    rt_cal: Option<RtCalibration>,
    mobility_slices: usize,
}

fn run_read(path: &Path) -> Result<RunRead, String> {
    let meta = read_nominee_meta(path)?;
    let get = |k: &str| {
        meta.iter()
            .find(|(m, _)| m == k)
            .map(|(_, v)| v.as_str())
            .ok_or_else(|| format!("{}: no {k} metadata", path.display()))
    };
    Ok(RunRead {
        spectrum: PathBuf::from(get(SPECTRUM_KEY)?),
        shift_ppm: get(SHIFT_KEY)?
            .parse()
            .map_err(|_| format!("{}: unreadable {SHIFT_KEY}", path.display()))?,
        rt_cal: decode_rt_cal(get(RT_CAL_KEY)?).map_err(|e| format!("{}: {e}", path.display()))?,
        mobility_slices: get(MOBILITY_SLICES_KEY)?
            .parse()
            .map_err(|_| format!("{}: unreadable {MOBILITY_SLICES_KEY}", path.display()))?,
    })
}

/// Fingerprint of one run's stored inputs: the content stamps of its nominee table and its
/// pass-2 table, and the byte size of its pass-1 PIN.
pub(crate) fn run_inputs(quant_dir: &Path, stem: &str) -> Result<String, String> {
    let stamp = |meta: Vec<(String, String)>, p: &Path| {
        meta.into_iter()
            .find(|(k, _)| k == STAMP_KEY)
            .map(|(_, v)| v)
            .ok_or_else(|| {
                format!(
                    "{}: no {STAMP_KEY} metadata; search run {stem} again with --dia --quant-dir",
                    p.display()
                )
            })
    };
    let nom = nominee_path(quant_dir, stem);
    let p2 = quant_dir.join(PASS2_DIR).join(format!("{stem}.parquet"));
    let pin = quant_dir.join(PASS2_DIR).join(format!("{stem}.pin"));
    let pin_bytes = std::fs::metadata(&pin)
        .map_err(|e| format!("{}: {e}", pin.display()))?
        .len();
    Ok(format!(
        "{stem}\t{}\t{}\t{pin_bytes}",
        stamp(read_nominee_meta(&nom)?, &nom)?,
        stamp(read_pass2_meta(&p2)?, &p2)?
    ))
}

/// Checks that the transfer table at `table` (its [`INPUTS_KEY`] value `built`) was built from
/// the inputs the study's runs hold now (`current`, [`run_inputs`] of every stem in order).
pub(crate) fn check_inputs(
    quant_dir: &Path,
    table: &Path,
    stem: &str,
    built: Option<&str>,
    current: &[String],
) -> Result<(), String> {
    let again = format!(
        "run `andes quant extract --quant-dir {} --run {stem}` again",
        quant_dir.display()
    );
    let Some(built) = built else {
        return Err(format!(
            "{}: no record of the inputs it was built from (stale transfer); {again}",
            table.display()
        ));
    };
    let built: Vec<&str> = built.lines().collect();
    let stale: Vec<&str> = current
        .iter()
        .enumerate()
        .filter(|(i, c)| built.get(*i) != Some(&c.as_str()))
        .map(|(_, c)| c.split('\t').next().unwrap_or_default())
        .collect();
    if stale.is_empty() && built.len() == current.len() {
        return Ok(());
    }
    Err(format!(
        "{}: built before run(s) {} were stored again (stale transfer); {again}",
        table.display(),
        if stale.is_empty() { "of another study".to_string() } else { stale.join(", ") }
    ))
}

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

/// Stores the nominees of a run with its spectrum file, its mobility slicing, its fragment m/z
/// shift, its andes-rt calibration and a content stamp.
pub(crate) fn write_nominees_table(
    out: &NomineesOut,
    noms: &[Nominee],
    shift_ppm: f64,
    rt_cal: Option<&RtCalibration>,
) -> Result<(), String> {
    let path = out.path;
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
    }
    let spectrum =
        std::fs::canonicalize(out.spectrum).unwrap_or_else(|_| out.spectrum.to_path_buf());
    let meta = [
        (SPECTRUM_KEY, spectrum.display().to_string()),
        (SHIFT_KEY, format!("{shift_ppm}")),
        (RT_CAL_KEY, encode_rt_cal(rt_cal)),
        (MOBILITY_SLICES_KEY, out.mobility_slices.to_string()),
        (STAMP_KEY, nominee_stamp(noms)),
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
                raw_pred_rt: n.raw_rt,
                rt_missing: n.rt_missing,
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
/// of the run they come from) get the flag `dia_llr`, an unknown nominating RT, a SpecId that
/// names their source run, and their raw andes-rt prediction through the receiving run's
/// calibration `rt_cal`; without a raw prediction or a calibration they keep the source run's
/// predicted RT and count as missing (they never enter a fit either way).
fn nominees_of(
    t: &NomineeTable,
    source: Option<&str>,
    rt_cal: Option<&RtCalibration>,
) -> Vec<Nominee> {
    let on_receiver = |i: usize| {
        rt_cal
            .filter(|_| source.is_some() && t.raw_pred_rt[i].is_finite())
            .map(|c| c.apply(t.raw_pred_rt[i]))
            .filter(|v| v.is_finite())
    };
    (0..t.len())
        .map(|i| Nominee {
            run: 0,
            charge: t.charge[i],
            decoy: t.decoy[i],
            prec_mz: t.prec_mz[i],
            pred_rt: match source {
                Some(_) => on_receiver(i).unwrap_or(t.pred_rt[i]),
                None => t.pred_rt[i],
            },
            rt_missing: match source {
                Some(_) => on_receiver(i).is_none(),
                None => t.rt_missing[i],
            },
            raw_rt: t.raw_pred_rt[i],
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
fn read_run_spectra(
    path: &Path,
    mobility_slices: usize,
) -> Result<(Vec<Spectrum>, Vec<Ms1Scan>), String> {
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
                .with_dia_mobility_slices(mobility_slices)
                .with_run_ms1_capture(true)
                .with_ms1_mobility(true);
            for r in reader.by_ref() {
                take(r.map_err(|e| e.to_string()));
            }
            reader.take_run_ms1()
        }
        #[cfg(not(feature = "timstof"))]
        {
            let _ = mobility_slices;
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
    let stems = crate::dia::study::study_runs(quant_dir)?;
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
    let inputs: Vec<String> = stems
        .iter()
        .map(|s| run_inputs(quant_dir, s))
        .collect::<Result<_, _>>()?;
    let read = run_read(&own_path)?;
    let own = read_nominees(&own_path, None)?;
    let own_keys: HashSet<u128> = keys_of(&own.peptide, &own.charge, &own.decoy)
        .into_iter()
        .collect();
    let mut noms = nominees_of(&own, None, None);
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
        noms.extend(nominees_of(&t, Some(&stems[o]), read.rt_cal.as_ref()));
    }
    let count = |r: std::ops::Range<usize>, decoy: bool| {
        noms[r].iter().filter(|n| n.decoy == decoy).count()
    };
    let unmapped = noms[n_own..].iter().filter(|n| n.rt_missing).count();
    if unmapped > 0 {
        log::warn!(
            "run {stem}: {unmapped} of {} transferred nominees have no raw andes-rt prediction \
             or the run has no andes-rt calibration; they keep their source run's predicted RT",
            noms.len() - n_own
        );
    }
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
    let spectrum = &read.spectrum;
    let (spectra, ms1) = read_run_spectra(spectrum, read.mobility_slices)?;
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
    let p2 =
        crate::dia::pass2::extract(&noms, &spectra, &[spectra.len()], &idx, &[read.shift_ppm]);
    drop(idx);
    drop(spectra);
    let table = crate::dia::pass2::to_table(&noms, &p2);
    drop(p2);
    let xfer_kept = table
        .feature(crate::dia::pass2::XFER_COLUMN)
        .map_or(0, |v| v.iter().filter(|&&x| x > 0.5).count());
    let mut out_meta: Vec<(&str, String)> = read_pass2_meta(
        &quant_dir.join(PASS2_DIR).join(format!("{stem}.parquet")),
    )?
    .into_iter()
    .filter(|(k, _)| k == crate::dia::study::RESCORE_KEY)
    .map(|(_, v)| (crate::dia::study::RESCORE_KEY, v))
    .collect();
    out_meta.push((RUNS_KEY, stems.join("\n")));
    out_meta.push((INPUTS_KEY, inputs.join("\n")));
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
            raw_pred_rt: vec![0.25],
            rt_missing: vec![false],
            apex_rt: vec![12.5],
            dia_llr: vec![30.0],
            mobility: vec![f64::NAN],
            ion_off: vec![0, 1],
            ion_mz: vec![300.0],
            ion_int: vec![1.0],
        };
        let own = nominees_of(&t, None, None);
        assert!(!own[0].xfer && own[0].dia_llr == 30.0 && own[0].apex_rt == 12.5);
        let x = nominees_of(&t, Some("run02"), None);
        assert!(x[0].xfer && x[0].decoy);
        assert_eq!(x[0].dia_llr, TRANSFER_DIA_LLR);
        assert!(x[0].apex_rt.is_nan());
        assert_eq!(&*x[0].spec_id, "xfer_run02_s_9_1");
        assert_eq!(x[0].ions, vec![(300.0, 1.0)]);
        assert!(!crate::dia::pass2::confident(&x[0]) && !crate::dia::pass2::null_member(&x[0]));
        assert!(crate::dia::pass2::null_member(&own[0]));
    }

    fn one_nominee_table(raw: f64, pred: f64) -> NomineeTable {
        NomineeTable {
            spec_id: vec!["s_9_1".into()],
            decoy: vec![false],
            charge: vec![2],
            peptide: vec!["-.PEPTIDEK.-".into()],
            proteins: vec!["P1".into()],
            seq: vec!["PEPTIDEK".into()],
            prec_mz: vec![450.0],
            pred_rt: vec![pred],
            raw_pred_rt: vec![raw],
            rt_missing: vec![!raw.is_finite()],
            apex_rt: vec![pred],
            dia_llr: vec![30.0],
            mobility: vec![f64::NAN],
            ion_off: vec![0, 1],
            ion_mz: vec![300.0],
            ion_int: vec![1.0],
        }
    }

    #[test]
    fn a_transferred_nominee_takes_the_receiving_runs_rt_calibration() {
        let receiver = RtCalibration::Linear {
            slope: 2.0,
            intercept: 1.0,
        };
        let t = one_nominee_table(10.0, 33.0);
        let x = nominees_of(&t, Some("run02"), Some(&receiver));
        assert_eq!(x[0].pred_rt, 21.0, "raw 10 through the receiving run's map");
        assert!(!x[0].rt_missing);
        let own = nominees_of(&t, None, Some(&receiver));
        assert_eq!(own[0].pred_rt, 33.0, "own nominees keep their own calibrated RT");
        let iso = RtCalibration::Isotonic {
            x: vec![0.0, 20.0],
            y: vec![5.0, 45.0],
        };
        assert_eq!(nominees_of(&t, Some("run02"), Some(&iso))[0].pred_rt, 25.0);
        let no_raw = one_nominee_table(f64::NAN, 33.0);
        let x = nominees_of(&no_raw, Some("run02"), Some(&receiver));
        assert_eq!(x[0].pred_rt, 33.0, "no raw prediction: the source run's RT");
        assert!(x[0].rt_missing);
        let x = nominees_of(&t, Some("run02"), None);
        assert_eq!(x[0].pred_rt, 33.0, "no calibration in the receiving run");
        assert!(x[0].rt_missing);
    }

    #[test]
    fn rt_calibration_metadata_round_trips_exactly() {
        for c in [
            None,
            Some(RtCalibration::Linear {
                slope: 0.1 + 0.2,
                intercept: -1e-7,
            }),
            Some(RtCalibration::Isotonic {
                x: vec![-3.5, 0.1 + 0.2, 1e300],
                y: vec![1.0 / 3.0, 2.0, 2.0],
            }),
        ] {
            assert_eq!(decode_rt_cal(&encode_rt_cal(c.as_ref())).unwrap(), c);
        }
        assert!(decode_rt_cal("isotonic 1").is_err());
        assert!(decode_rt_cal("cubic 1 2").is_err());
    }

    fn nominee(rt: f64) -> Nominee {
        Nominee {
            run: 0,
            charge: 2,
            decoy: false,
            prec_mz: 450.0,
            pred_rt: rt,
            rt_missing: false,
            raw_rt: rt / 10.0,
            apex_rt: rt,
            dia_llr: 30.0,
            mobility: f64::NAN,
            ions: vec![(300.0, 1.0)],
            seq: "PEPTIDEK".into(),
            peptide: "-.PEPTIDEK.-".into(),
            proteins: "P1".into(),
            spec_id: "s_9_1".into(),
            xfer: false,
        }
    }

    fn store(quant_dir: &Path, stem: &str, rt: f64) {
        let spectrum = quant_dir.join(format!("{stem}.mzML"));
        let path = nominee_path(quant_dir, stem);
        let out = NomineesOut {
            path: &path,
            spectrum: &spectrum,
            mobility_slices: 3,
        };
        write_nominees_table(&out, &[nominee(rt)], 1.5, None).unwrap();
        let pin = quant_dir.join(format!("{stem}.pin"));
        std::fs::write(&pin, format!("SpecId\tLabel\n{stem}/s_9_1\t1\n")).unwrap();
        let table = Pass2Table {
            spec_id: vec!["s_9_1".into()],
            run: vec![0],
            decoy: vec![false],
            charge: vec![2],
            peptide: vec!["-.PEPTIDEK.-".into()],
            proteins: vec!["P1".into()],
            seq: vec!["PEPTIDEK".into()],
            prec_mz: vec![450.0],
            features: vec![("p2_llr".into(), vec![rt as f32])],
            ..Default::default()
        };
        crate::dia::study::store_run(quant_dir, &spectrum, &pin, &table, None, &Default::default()).unwrap();
    }

    #[test]
    fn nominee_metadata_keeps_the_reading_options() {
        let dir = tempfile::tempdir().unwrap();
        store(dir.path(), "run01", 12.0);
        let r = run_read(&nominee_path(dir.path(), "run01")).unwrap();
        assert_eq!(r.mobility_slices, 3);
        assert_eq!(r.shift_ppm, 1.5);
        assert!(r.rt_cal.is_none());
        assert!(r.spectrum.ends_with("run01.mzML"));
    }

    #[test]
    fn storing_a_run_again_after_extract_makes_the_transfer_stale() {
        let dir = tempfile::tempdir().unwrap();
        let qd = dir.path();
        let stems = ["run01".to_string(), "run02".to_string()];
        store(qd, "run01", 12.0);
        store(qd, "run02", 14.0);
        let inputs = |qd: &Path| -> Vec<String> {
            stems.iter().map(|s| run_inputs(qd, s).unwrap()).collect()
        };
        let built = inputs(qd).join("\n");
        let table = xfer_path(qd, "run01");
        check_inputs(qd, &table, "run01", Some(&built), &inputs(qd)).unwrap();
        store(qd, "run02", 14.0);
        check_inputs(qd, &table, "run01", Some(&built), &inputs(qd)).unwrap();
        store(qd, "run02", 15.0);
        let e = check_inputs(qd, &table, "run01", Some(&built), &inputs(qd)).unwrap_err();
        assert!(e.contains("run(s) run02 were stored again"), "{e}");
        assert!(e.contains("--run run01"), "{e}");
        let e = check_inputs(qd, &table, "run01", None, &inputs(qd)).unwrap_err();
        assert!(e.contains("stale transfer"), "{e}");
    }
}
