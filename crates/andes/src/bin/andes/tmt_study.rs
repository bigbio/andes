//! Study mode of `andes quant finalize` for isobaric (TMT/iTRAQ) runs.
//!
//! Each run is searched on its own with `--tmt` and without rescoring, so its
//! `<stem>.tmt.parquet` holds every rank-1 target PSM with its reporter ions and the run's
//! quantification gate on purity and summed S/N.
//! The runs' PINs are rescored together (one Percolator run over every run gives one score
//! scale and one FDR), and `andes quant finalize --quant-dir <dir>` on the directory of the
//! per-run tables then:
//! 1. joins the pooled q, PEP and score onto every run's rows on (run, SpecId) and keeps the
//!    rows at pooled q <= `--psm-fdr` (and PEP <= `--pep` when given), as every quantification;
//! 2. keeps the run's gate (purity and summed S/N, quantification-quality criteria only) and
//!    weights the rows that pass by purity squared times 1 - pooled PEP (`quant_weight`);
//! 3. links each row's peptide to the study's protein group at `--protein-fdr` (the inference
//!    of `andes inference` over the same pooled results);
//! 4. adds the run's plex, fraction, instrument and the sample and role (sample, reference or
//!    blank) of every channel from `--sdrf`;
//! 5. writes `<base>.tmt.parquet` (`output::tmt_study_table`) and with `--output-parquet` the
//!    QPX views `quantms.feature.parquet` (the rows in a group) and `quantms.pg.parquet`.
//!
//! The reporter values are those of the run tables (raw or impurity-corrected, as the runs
//! were searched); nothing is normalised here.

use std::collections::HashMap;
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash, Hasher};
use std::path::{Path, PathBuf};

use arrow::array::{
    Array, BooleanArray, Float32Array, Float64Array, Int32Array, StringArray, UInt8Array,
};
use arrow::record_batch::RecordBatch;
use output::feature_parquet::{FeatureGroup, FeatureRecord};
use output::quant_out::{TMT_CORRECTED_KEY, TMT_MBR_KEY, TMT_PLEX_KEY, TMT_UNMIX_KEY};
use output::tmt_study_table::{
    is_tmt_study_table, RunDesignColumns, StudyGroup, TmtStudyLayout, TmtStudyRow, TmtStudyWriter,
    TMT_STUDY_PEP_KEY, TMT_STUDY_PROTEIN_FDR_KEY, TMT_STUDY_PSM_FDR_KEY, TMT_STUDY_RUNS_KEY,
    TMT_STUDY_SDRF_KEY,
};
use output::FeatureParquetWriter;
use quant::isobaric::Plex;

use crate::quant_cmd::{GroupTally, ProteinGroups};
use crate::sdrf::{ChannelSample, IsobaricDesign};

/// File-name suffix of a per-run isobaric table.
pub(crate) const TMT_SUFFIX: &str = ".tmt.parquet";

/// Per-run isobaric tables of a quant directory (study tables left out), sorted; empty when
/// it holds none.
pub(crate) fn study_tables(quant_dir: &Path, output: &Path) -> Result<Vec<PathBuf>, String> {
    let rd = std::fs::read_dir(quant_dir).map_err(|e| format!("{}: {e}", quant_dir.display()))?;
    let out = with_suffix(output, TMT_SUFFIX);
    let mut paths = Vec::new();
    for e in rd.filter_map(|e| e.ok()) {
        let p = e.path();
        let is_table = p
            .file_name()
            .and_then(|n| n.to_str())
            .is_some_and(|n| n.ends_with(TMT_SUFFIX));
        if !is_table || p == out || !p.is_file() || is_tmt_study_table(&p)? {
            continue;
        }
        paths.push(p);
    }
    paths.sort();
    Ok(paths)
}

pub(crate) fn with_suffix(base: &Path, suffix: &str) -> PathBuf {
    let mut p = base.as_os_str().to_owned();
    p.push(suffix);
    PathBuf::from(p)
}

/// One row of a per-run table.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RunRow {
    pub filename: String,
    pub scannr: i32,
    pub spec_id: String,
    pub quant_scan: String,
    pub quant_ms_level: u8,
    pub rt: Option<f64>,
    pub peptide: String,
    pub charge: u8,
    pub proteins: String,
    pub q_value: Option<f64>,
    pub purity: Option<f32>,
    pub intensities: Vec<f32>,
    pub raw: Option<Vec<f32>>,
    pub intensity_sum: f32,
    pub sn_sum: Option<f32>,
    pub sn: Vec<Option<f32>>,
    pub interference: Option<f32>,
    pub unmixed: Vec<Option<f32>>,
    pub gate_passed: bool,
}

/// A per-run table with its layout.
pub(crate) struct RunTable {
    pub layout: TmtStudyLayout,
    pub unmix_mode: String,
    pub rows: Vec<RunRow>,
}

fn col<'a, T: 'static>(b: &'a RecordBatch, name: &str, path: &Path) -> Result<&'a T, String> {
    let i = b
        .schema()
        .index_of(name)
        .map_err(|_| format!("{}: no column {name}", path.display()))?;
    b.column(i).as_any().downcast_ref::<T>().ok_or_else(|| {
        format!(
            "{}: column {name} has type {}",
            path.display(),
            b.column(i).data_type()
        )
    })
}

fn opt_f32(a: &Float32Array, i: usize) -> Option<f32> {
    a.is_valid(i).then(|| a.value(i))
}

fn opt_f64(a: &Float64Array, i: usize) -> Option<f64> {
    a.is_valid(i).then(|| a.value(i))
}

/// Reads a per-run `<stem>.tmt.parquet` of a search without rescoring.
pub(crate) fn read_run_table(path: &Path) -> Result<RunTable, String> {
    let t = output::read_tmt_parquet(path)?;
    let meta = |k: &str| t.meta.get(k).cloned();
    if meta(TMT_MBR_KEY).is_some() {
        return Err(format!(
            "{}: written by a search with --mbr; a TMT study joins the pooled rescoring onto \
             runs searched without rescoring",
            path.display()
        ));
    }
    let plex = Plex::parse(
        &meta(TMT_PLEX_KEY)
            .ok_or_else(|| format!("{}: no {TMT_PLEX_KEY} metadata", path.display()))?,
    )
    .map_err(|e| format!("{}: {e}", path.display()))?;
    let corrected = meta(TMT_CORRECTED_KEY).as_deref() == Some("true");
    let unmix_mode = meta(TMT_UNMIX_KEY).unwrap_or_else(|| "off".to_string());
    let layout = TmtStudyLayout {
        plex,
        corrected,
        unmixed: unmix_mode != "off",
    };
    let prefix = if matches!(plex, Plex::Itraq4 | Plex::Itraq8) {
        "itraq"
    } else {
        "tmt"
    };
    let names: Vec<&str> = plex.channels().iter().map(|c| c.name).collect();
    let mut rows = Vec::new();
    for b in &t.batches {
        let block = |head: &str| -> Result<Vec<&Float32Array>, String> {
            names
                .iter()
                .map(|c| col::<Float32Array>(b, &format!("{head}{prefix}_{c}"), path))
                .collect()
        };
        let filename = col::<StringArray>(b, "filename", path)?;
        let scannr = col::<Int32Array>(b, "scannr", path)?;
        let spec_id = col::<StringArray>(b, "spec_id", path)?;
        let quant_scan = col::<StringArray>(b, "quant_scan", path)?;
        let level = col::<UInt8Array>(b, "quant_ms_level", path)?;
        let rt = col::<Float64Array>(b, "rt", path)?;
        let peptide = col::<StringArray>(b, "peptide", path)?;
        let charge = col::<UInt8Array>(b, "charge", path)?;
        let proteins = col::<StringArray>(b, "proteins", path)?;
        let q = col::<Float64Array>(b, "psm_q_value", path)?;
        let purity = col::<Float32Array>(b, "purity", path)?;
        let shown = block("")?;
        let raw = if corrected {
            block("raw_")?
        } else {
            Vec::new()
        };
        let sum = col::<Float32Array>(b, "reporter_intensity_sum", path)?;
        let sn_sum = col::<Float32Array>(b, "reporter_sn_sum", path)?;
        let sn = block("sn_")?;
        let (interference, unmix) = if layout.unmixed {
            (
                Some(col::<Float32Array>(b, "interference", path)?),
                block("unmix_")?,
            )
        } else {
            (None, Vec::new())
        };
        let gate = col::<BooleanArray>(b, "quant_gate_passed", path)?;
        for i in 0..b.num_rows() {
            rows.push(RunRow {
                filename: filename.value(i).to_string(),
                scannr: scannr.value(i),
                spec_id: spec_id.value(i).to_string(),
                quant_scan: quant_scan.value(i).to_string(),
                quant_ms_level: level.value(i),
                rt: opt_f64(rt, i),
                peptide: peptide.value(i).to_string(),
                charge: charge.value(i),
                proteins: proteins.value(i).to_string(),
                q_value: opt_f64(q, i),
                purity: opt_f32(purity, i),
                intensities: shown.iter().map(|a| a.value(i)).collect(),
                raw: corrected.then(|| raw.iter().map(|a| a.value(i)).collect()),
                intensity_sum: sum.value(i),
                sn_sum: opt_f32(sn_sum, i),
                sn: sn.iter().map(|a| opt_f32(a, i)).collect(),
                interference: interference.and_then(|a| opt_f32(a, i)),
                unmixed: unmix.iter().map(|a| opt_f32(a, i)).collect(),
                gate_passed: gate.value(i),
            });
        }
    }
    if let Some(r) = rows.iter().find(|r| r.q_value.is_some()) {
        return Err(format!(
            "{}: {} carries a per-run q-value: the run was rescored on its own and its table \
             holds only the PSMs at that run's threshold; search the runs of a TMT study \
             without --rescore and rescore their PINs together",
            path.display(),
            r.spec_id
        ));
    }
    Ok(RunTable {
        layout,
        unmix_mode,
        rows,
    })
}

/// Join key of a PSM: its run and its SpecId.
pub(crate) fn psm_key(run: &str, spec_id: &str) -> u128 {
    let half = |salt: u8| {
        let mut h = BuildHasherDefault::<DefaultHasher>::default().build_hasher();
        salt.hash(&mut h);
        run.hash(&mut h);
        spec_id.hash(&mut h);
        h.finish()
    };
    (u128::from(half(0)) << 64) | u128::from(half(1))
}

/// A pooled target PSM at the PSM threshold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Pooled {
    pub q: f64,
    pub pep: f64,
    pub score: Option<f64>,
    pub joined: bool,
}

/// The pooled target PSMs at q <= `max_q` and, when given, PEP <= `max_pep`, keyed by (run,
/// SpecId); the run is the SpecId's `<run>/` prefix. Returns them with the number of target
/// PSMs read.
pub(crate) fn pooled_targets(
    paths: &[PathBuf],
    max_q: f64,
    max_pep: Option<f64>,
) -> Result<(HashMap<u128, Pooled>, usize), String> {
    let mut map = HashMap::new();
    let mut read = 0usize;
    for p in paths {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", p.display());
        let file = std::fs::File::open(p).map_err(|e| err(&e))?;
        let mut reader =
            output::PsmResultReader::new(std::io::BufReader::with_capacity(1 << 20, file))
                .map_err(|e| err(&e))?;
        while let Some(row) = reader.next_row().map_err(|e| err(&e))? {
            read += 1;
            if row.q_value > max_q || max_pep.is_some_and(|c| row.pep > c) {
                continue;
            }
            let key = psm_key(output::spec_id_run(row.psm_id), row.psm_id);
            let v = Pooled {
                q: row.q_value,
                pep: row.pep,
                score: row.score,
                joined: false,
            };
            if map.insert(key, v).is_some() {
                return Err(err(&format!(
                    "duplicate PSMId `{}` in the pooled results",
                    row.psm_id
                )));
            }
        }
    }
    Ok((map, read))
}

/// The SDRF columns of one run: the samples of the plex channels in kit order.
pub(crate) fn run_design(
    design: Option<&IsobaricDesign>,
    run: &str,
    plex: Plex,
) -> Result<RunDesignColumns, String> {
    let channels = plex.channels();
    let Some(d) = design else {
        return Ok(RunDesignColumns {
            channels: vec![None; channels.len()],
            ..Default::default()
        });
    };
    let r = d.runs.get(run).ok_or_else(|| {
        format!(
            "--sdrf {}: run {run} has no `comment[data file]` row",
            d.path.display()
        )
    })?;
    let mut out: Vec<Option<(String, &'static str)>> = vec![None; channels.len()];
    for (name, ChannelSample { sample, role }) in &r.channels {
        let i = channels
            .iter()
            .position(|c| c.name == name.as_str())
            .or_else(|| {
                channels
                    .iter()
                    .position(|c| c.name == format!("{name}N").as_str())
            })
            .ok_or_else(|| {
                format!(
                    "--sdrf {}: run {run} names channel {name}, which {} does not have",
                    d.path.display(),
                    plex.name()
                )
            })?;
        if out[i].is_some() {
            return Err(format!(
                "--sdrf {}: run {run} names channel {} twice",
                d.path.display(),
                channels[i].name
            ));
        }
        out[i] = Some((sample.clone(), role.as_str()));
    }
    Ok(RunDesignColumns {
        plex: Some(r.plex),
        fraction: r.fraction.clone(),
        instrument: r.instrument.clone(),
        channels: out,
    })
}

/// Counts of the join and the gate over the study.
#[derive(Debug, Default, Clone, PartialEq)]
pub(crate) struct StudyCounts {
    pub run_rows: usize,
    pub kept: usize,
    pub not_at_q: usize,
    pub gate_passed: usize,
    pub grouped: usize,
    pub unique: usize,
    /// SpecId prefix (`<spectrum file stem>/`) to the run whose rows carry it.
    pub spec_prefixes: HashMap<String, String>,
}

/// Joins one run's rows to the pooled results, applies the gate and links the groups.
pub(crate) fn study_rows(
    rows: Vec<RunRow>,
    pooled: &mut HashMap<u128, Pooled>,
    groups: &ProteinGroups<'_>,
    tally: &mut GroupTally,
    counts: &mut StudyCounts,
) -> Result<Vec<TmtStudyRow>, String> {
    let mut out = Vec::with_capacity(rows.len());
    for r in rows {
        counts.run_rows += 1;
        let prefix = output::spec_id_run(&r.spec_id);
        match counts.spec_prefixes.get(prefix) {
            Some(run) if *run != r.filename => {
                return Err(format!(
                    "runs {run} and {} have SpecIds with the same prefix `{prefix}/` (e.g. {}); \
                     the pooled results cannot be told apart by run",
                    r.filename, r.spec_id
                ));
            }
            Some(_) => {}
            None => {
                counts
                    .spec_prefixes
                    .insert(prefix.to_string(), r.filename.clone());
            }
        }
        let Some(p) = pooled.get_mut(&psm_key(prefix, &r.spec_id)) else {
            counts.not_at_q += 1;
            continue;
        };
        p.joined = true;
        let p = *p;
        let passed = r.gate_passed;
        counts.kept += 1;
        counts.gate_passed += usize::from(passed);
        let (sequence, _) = crate::quant_cmd::parse_peptidoform(&r.peptide);
        let group = tally.assign(groups, &sequence).map(|(accessions, g)| {
            counts.grouped += 1;
            counts.unique += usize::from(g.unique);
            StudyGroup {
                accessions: accessions
                    .iter()
                    .map(|a| a.0.as_str())
                    .collect::<Vec<_>>()
                    .join(";"),
                anchor: g.anchor,
                unique: g.unique,
                q_value: g.global_q,
            }
        });
        out.push(TmtStudyRow {
            filename: r.filename,
            scannr: r.scannr,
            spec_id: r.spec_id,
            quant_scan: r.quant_scan,
            quant_ms_level: r.quant_ms_level,
            rt: r.rt,
            peptide: r.peptide,
            sequence,
            charge: r.charge,
            proteins: r.proteins,
            psm_q_value: p.q,
            psm_pep: p.pep,
            psm_score: p.score,
            purity: r.purity,
            intensities: r.intensities,
            raw: r.raw,
            intensity_sum: r.intensity_sum,
            sn_sum: r.sn_sum,
            sn: r.sn,
            interference: r.interference,
            unmixed: r.unmixed,
            quant_gate_passed: passed,
            quant_weight: if passed {
                crate::quant_run::quant_weight(r.purity, Some(p.pep))
            } else {
                0.0
            },
            group,
        });
    }
    Ok(out)
}

/// Theoretical m/z of a peptidoform (residues, `[UNIMOD:n]` or `[+mass]` tags) at `charge`;
/// `None` when a tag has no bundled Unimod mass.
pub(crate) fn peptidoform_mz(peptidoform: &str, charge: u8) -> Option<f64> {
    let mut mods: HashMap<String, f64> = HashMap::new();
    let mut rest = peptidoform;
    while let Some(start) = rest.find('[') {
        let end = rest[start..].find(']')? + start;
        let tag = &rest[start + 1..end];
        if !(tag.starts_with('+') || tag.starts_with('-')) {
            mods.insert(
                tag.to_string(),
                model::isomer_class::bundled_unimod_mass(tag)?,
            );
        }
        rest = &rest[end + 1..];
    }
    let masses = quant::ms2::residue_masses(peptidoform, &mods).ok()?;
    let z = f64::from(charge.max(1));
    Some((masses.iter().sum::<f64>() + model::mass::H2O) / z + model::mass::PROTON)
}

/// The QPX feature of a study row in a group; `None` when the m/z cannot be computed.
pub(crate) fn feature_record(
    r: &TmtStudyRow,
    plex: Plex,
    feature_id: i64,
) -> Option<FeatureRecord> {
    let g = r.group.as_ref()?;
    let mz = peptidoform_mz(&r.peptide, r.charge)? as f32;
    let (_, modifications) = crate::quant_cmd::parse_peptidoform(&r.peptide);
    let labels: Vec<String> = plex
        .channels()
        .iter()
        .map(|c| format!("{}{}", plex.label_prefix(), c.name))
        .collect();
    Some(FeatureRecord {
        feature_id,
        sequence: r.sequence.clone(),
        peptidoform: r.peptide.clone(),
        modifications,
        charge: i16::from(r.charge),
        pep: Some(r.psm_pep),
        is_decoy: false,
        calculated_mz: mz,
        observed_mz: mz,
        mass_error_ppm: None,
        additional_scores: vec![
            ("psm_q_value".into(), r.psm_q_value, false),
            ("quant_weight".into(), f64::from(r.quant_weight), false),
        ],
        run_file_name: r.filename.clone(),
        cv_params: Vec::new(),
        scan: if r.scannr > 0 {
            vec![r.scannr]
        } else {
            Vec::new()
        },
        rt: r.rt.map(|t| t as f32),
        missed_cleavages: None,
        intensities: labels
            .into_iter()
            .zip(r.intensities.iter().copied())
            .collect(),
        additional_intensities: Vec::new(),
        proteins: g
            .accessions
            .split(';')
            .map(|a| (a.to_string(), None, None, None, None))
            .collect(),
        id_run_file_name: Some(r.filename.clone()),
        rt_start: None,
        rt_stop: None,
        psm_ids: Vec::new(),
        mbr: None,
        group: Some(FeatureGroup {
            anchor: g.anchor.clone(),
            unique: g.unique,
            global_q: g.q_value,
        }),
    })
}

/// `andes quant finalize` on a TMT study (module doc).
pub(crate) fn finalize(
    a: &crate::cli::QuantFinalizeArgs,
    tables: &[PathBuf],
) -> Result<(), String> {
    let t0 = std::time::Instant::now();
    let database = a.database.as_ref().ok_or(
        "a TMT study needs --database (with the search's database options) for the protein \
         groups",
    )?;
    if a.psms.is_empty() || a.decoy_psms.is_empty() {
        return Err(
            "a TMT study needs --psms and --decoy-psms: the pooled Percolator results of every \
             run's PIN"
                .into(),
        );
    }
    let psm_fdr = a.proteins.psm_fdr;
    let design = a
        .sdrf
        .as_deref()
        .map(crate::sdrf::isobaric_design)
        .transpose()?;
    log::info!(
        target: "tmt",
        "TMT study: {} run tables in {}, pooled results {}, {}",
        tables.len(),
        a.quant_dir.display(),
        a.psms.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", "),
        match &design {
            Some(d) => format!("SDRF {} ({} runs, {} plexes)", d.path.display(), d.runs.len(), d.plexes),
            None => "no SDRF (plex, fraction, instrument and channel samples are empty)".into(),
        }
    );
    if design.is_none() {
        log::warn!(
            "a TMT study without --sdrf: the table has no plex, fraction, instrument or channel \
             samples"
        );
    }
    let inf = crate::inference_cmd::infer_from_files(database, &a.proteins, &a.psms, &a.decoy_psms)
        .map_err(|e| e.to_string())?;
    let groups = ProteinGroups {
        inf: &inf,
        max_q: a.proteins.protein_fdr,
    };
    let (mut pooled, read) = pooled_targets(&a.psms, psm_fdr, a.pep)?;
    let pooled_at_q = pooled.len();
    log::info!(
        target: "tmt",
        "pooled results: {read} target PSMs, {pooled_at_q} at q <= {psm_fdr}{} [{:.2}s]",
        a.pep.map_or_else(String::new, |p| format!(" and PEP <= {p}")),
        t0.elapsed().as_secs_f64()
    );

    let out_path = with_suffix(&a.output, TMT_SUFFIX);
    if let Some(parent) = out_path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent).map_err(|e| format!("create {}: {e}", parent.display()))?;
    }
    let mut features = match &a.output_parquet {
        Some(dir) => {
            std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
            let path = dir.join("quantms.feature.parquet");
            Some((FeatureParquetWriter::new(&path), path))
        }
        None => None,
    };
    let mut writer: Option<(TmtStudyWriter, TmtStudyLayout, String)> = None;
    let mut tally = GroupTally::default();
    let mut counts = StudyCounts::default();
    let mut runs: Vec<String> = Vec::new();
    let mut no_mz = 0usize;
    let mut feature_id = 0i64;
    for path in tables {
        let table = read_run_table(path)?;
        match &writer {
            None => {
                let meta = vec![
                    (TMT_PLEX_KEY, table.layout.plex.name().to_string()),
                    (TMT_CORRECTED_KEY, table.layout.corrected.to_string()),
                    (TMT_UNMIX_KEY, table.unmix_mode.clone()),
                    (TMT_STUDY_PSM_FDR_KEY, psm_fdr.to_string()),
                    (
                        TMT_STUDY_PROTEIN_FDR_KEY,
                        a.proteins.protein_fdr.to_string(),
                    ),
                    (
                        TMT_STUDY_PEP_KEY,
                        a.pep.map_or_else(|| "none".to_string(), |p| p.to_string()),
                    ),
                    (TMT_STUDY_RUNS_KEY, tables.len().to_string()),
                    (
                        TMT_STUDY_SDRF_KEY,
                        design
                            .as_ref()
                            .and_then(|d| d.path.file_name())
                            .map(|n| n.to_string_lossy().into_owned())
                            .unwrap_or_default(),
                    ),
                ];
                let w = TmtStudyWriter::create(&out_path, table.layout, &meta)?;
                writer = Some((w, table.layout, table.unmix_mode.clone()));
            }
            Some((_, layout, mode)) if *layout != table.layout || *mode != table.unmix_mode => {
                return Err(format!(
                    "{}: {} (correction {}, unmixing {}) differs from the study's first table \
                     ({}, correction {}, unmixing {}); search every run with the same --tmt, \
                     --tmt-correction and --tmt-unmix",
                    path.display(),
                    table.layout.plex.name(),
                    table.layout.corrected,
                    table.unmix_mode,
                    layout.plex.name(),
                    layout.corrected,
                    mode
                ));
            }
            Some(_) => {}
        }
        let (w, layout, _) = writer.as_mut().expect("opened above");
        let mut by_run: Vec<(String, Vec<RunRow>)> = Vec::new();
        for r in table.rows {
            match by_run.iter_mut().find(|(n, _)| *n == r.filename) {
                Some((_, v)) => v.push(r),
                None => by_run.push((r.filename.clone(), vec![r])),
            }
        }
        for (run, rows) in by_run {
            if runs.contains(&run) {
                return Err(format!(
                    "{}: run {run} is also in another table of {}",
                    path.display(),
                    a.quant_dir.display()
                ));
            }
            let d = run_design(design.as_ref(), &run, layout.plex)?;
            let before = counts.clone();
            let rows = study_rows(rows, &mut pooled, &groups, &mut tally, &mut counts)?;
            if rows.is_empty() && counts.run_rows > before.run_rows {
                log::warn!(
                    "run {run}: none of its {} PSMs is in the pooled results at q <= {psm_fdr}; \
                     were its PIN rows rescored with the others?",
                    counts.run_rows - before.run_rows
                );
            }
            log::info!(
                target: "tmt",
                "run {run}: {} PSMs with reporters, {} at pooled q <= {psm_fdr}, {} passing the \
                 gate, {} in a protein group at q <= {}{}",
                counts.run_rows - before.run_rows,
                counts.kept - before.kept,
                counts.gate_passed - before.gate_passed,
                counts.grouped - before.grouped,
                a.proteins.protein_fdr,
                d.plex.map_or(String::new(), |p| format!(
                    " (plex {p}, fraction {})",
                    d.fraction.as_deref().unwrap_or("none")
                ))
            );
            w.push_run(&rows, &d)?;
            if let Some((fw, fpath)) = features.as_mut() {
                let mut recs = Vec::new();
                for r in rows.iter().filter(|r| r.group.is_some()) {
                    feature_id += 1;
                    match feature_record(r, layout.plex, feature_id) {
                        Some(f) => recs.push(f),
                        None => {
                            feature_id -= 1;
                            no_mz += 1;
                        }
                    }
                }
                fw.push(recs)
                    .map_err(|e| format!("write {}: {e}", fpath.display()))?;
            }
            runs.push(run);
        }
    }
    let (w, _, _) =
        writer.ok_or_else(|| format!("{}: no TMT run tables", a.quant_dir.display()))?;
    let written = w.rows();
    w.finish()?;
    let unjoined = pooled.values().filter(|p| !p.joined).count();
    if let Some(d) = &design {
        let missing = d.runs.keys().filter(|r| !runs.contains(r)).count();
        if missing > 0 {
            log::warn!(
                "--sdrf {}: {missing} of its {} runs have no table in {}",
                d.path.display(),
                d.runs.len(),
                a.quant_dir.display()
            );
        }
    }
    if no_mz > 0 {
        log::warn!(
            "{no_mz} grouped rows left out of quantms.feature.parquet: a modification has no \
             bundled Unimod mass, so the m/z is unknown (they are in {})",
            out_path.display()
        );
    }
    if let (Some((fw, fpath)), Some(dir)) = (features, &a.output_parquet) {
        fw.finish()
            .map_err(|e| format!("write {}: {e}", fpath.display()))?;
        log::info!(
            target: "tmt",
            "wrote {} ({} features: the rows in a protein group)",
            fpath.display(),
            feature_id
        );
        tally.write_pg(dir, &groups, &runs)?;
    }
    let pct = |n: usize, d: usize| {
        if d == 0 {
            0.0
        } else {
            100.0 * n as f64 / d as f64
        }
    };
    let s = |label: &str, n: usize, note: String| {
        log::info!(target: "summary", "{label:<44}{n:>12}   {note}");
    };
    log::info!(target: "summary", "TMT study  {} runs  {}", runs.len(), out_path.display());
    s(
        "pooled target PSMs at q <= --psm-fdr",
        pooled_at_q,
        format!("of {read}"),
    );
    s(
        "  joined to a run's reporter row",
        pooled_at_q - unjoined,
        format!("{:.1} %", pct(pooled_at_q - unjoined, pooled_at_q)),
    );
    s(
        "  without a reporter row",
        unjoined,
        "no quantification scan, scan range short of the reporters, or run not in --quant-dir"
            .into(),
    );
    s("run PSMs with reporters", counts.run_rows, String::new());
    s("  above q or not rescored", counts.not_at_q, String::new());
    s("rows written", written, String::new());
    s(
        "  passing the gate",
        counts.gate_passed,
        format!(
            "{:.1} %; failed the run's gate (purity, S/N) {}; weight purity^2 x (1 - PEP)",
            pct(counts.gate_passed, written),
            written - counts.gate_passed
        ),
    );
    s(
        "  in a protein group at --protein-fdr",
        counts.grouped,
        format!(
            "{} unique; no group {}, group above the cut {}",
            counts.unique, tally.no_group, tally.above_cut
        ),
    );
    s(
        "target protein groups at --protein-fdr",
        inf.groups
            .iter()
            .filter(|g| g.q_value <= a.proteins.protein_fdr)
            .count(),
        format!("{} with rows", tally.counts.len()),
    );
    log::info!(
        target: "summary",
        "wall {:.1} s",
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_row(run: &str, scan: i32, peptide: &str, gate: bool) -> RunRow {
        RunRow {
            filename: run.into(),
            scannr: scan,
            spec_id: format!("{run}/scan={scan}_{scan}_1"),
            quant_scan: format!("ms3_{scan}"),
            quant_ms_level: 3,
            rt: Some(1.0),
            peptide: peptide.into(),
            charge: 2,
            proteins: "sp|P1|A".into(),
            q_value: None,
            purity: Some(0.8),
            intensities: vec![1.0; 6],
            raw: None,
            intensity_sum: 6.0,
            sn_sum: None,
            sn: vec![None; 6],
            interference: None,
            unmixed: Vec::new(),
            gate_passed: gate,
        }
    }

    fn pooled(entries: &[(&str, f64, f64)]) -> HashMap<u128, Pooled> {
        entries
            .iter()
            .map(|&(id, q, pep)| {
                (
                    psm_key(output::spec_id_run(id), id),
                    Pooled {
                        q,
                        pep,
                        score: Some(1.0),
                        joined: false,
                    },
                )
            })
            .collect()
    }

    fn inference(peptides: &[(&str, &str)]) -> output::Inference {
        let mut inf = output::Inference::default();
        for (i, (pep, prot)) in peptides.iter().enumerate() {
            inf.groups.push(output::protein_inference::ProteinGroup {
                id: i,
                lead: prot.to_string(),
                families: vec![prot.to_string()],
                accessions: vec![prot.to_string()],
                peptides: 1,
                unique: 1,
                razor: 0,
                psms: 1,
                score: 1.0,
                q_value: 0.0,
                subsumed: Vec::new(),
                peptides_per_accession: Vec::new(),
                is_decoy: false,
            });
            inf.peptide_group.insert(pep.to_string(), i);
        }
        inf
    }

    #[test]
    fn rows_join_on_run_and_spec_id_and_the_weight_uses_the_pooled_pep() {
        let inf = inference(&[("PEPTIDEK", "sp|P1|A")]);
        let groups = ProteinGroups {
            inf: &inf,
            max_q: 0.01,
        };
        let mut map = pooled(&[
            ("a/scan=1_1_1", 0.001, 0.01),
            ("a/scan=2_2_1", 0.002, 0.2),
            ("b/scan=3_3_1", 0.003, 0.01),
        ]);
        let rows = vec![
            run_row("a", 1, "PEPTIDEK", true),
            run_row("a", 2, "PEPTIDEK", true),
            run_row("a", 3, "PEPTIDEK", true),
            run_row("a", 4, "PEPTIDEK", false),
        ];
        let mut tally = GroupTally::default();
        let mut c = StudyCounts::default();
        let out = study_rows(rows, &mut map, &groups, &mut tally, &mut c).unwrap();
        assert_eq!(out.len(), 2, "scan 3 is run b's, scan 4 is not rescored");
        assert_eq!(c.not_at_q, 2);
        assert!(out[0].quant_gate_passed && out[1].quant_gate_passed);
        assert!((out[0].quant_weight - 0.64 * 0.99).abs() < 1e-6);
        assert!((out[1].quant_weight - 0.64 * 0.8).abs() < 1e-6);
        assert_eq!(c.gate_passed, 2);
        let joined: Vec<bool> = {
            let mut v: Vec<(f64, bool)> = map.values().map(|p| (p.q, p.joined)).collect();
            v.sort_by(|x, y| x.0.total_cmp(&y.0));
            v.into_iter().map(|x| x.1).collect()
        };
        assert_eq!(joined, vec![true, true, false]);
        assert_eq!(out[0].psm_q_value, 0.001);
    }

    #[test]
    fn runs_sharing_a_spec_id_prefix_are_refused() {
        let inf = inference(&[("PEPTIDEK", "sp|P1|A")]);
        let groups = ProteinGroups {
            inf: &inf,
            max_q: 0.01,
        };
        let mut c = StudyCounts::default();
        let mut tally = GroupTally::default();
        let mut map = HashMap::new();
        study_rows(
            vec![run_row("a", 1, "PEPTIDEK", true)],
            &mut map,
            &groups,
            &mut tally,
            &mut c,
        )
        .unwrap();
        let mut r = run_row("b", 1, "PEPTIDEK", true);
        r.spec_id = "a/scan=2_2_1".into();
        let e = study_rows(vec![r], &mut map, &groups, &mut tally, &mut c).unwrap_err();
        assert!(e.contains("same prefix `a/`"), "{e}");
    }

    #[test]
    fn the_pep_threshold_applies_with_the_q_threshold() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("t.psms");
        std::fs::write(
            &p,
            "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
             a/s_1_1\t3\t0.001\t0.01\tK.PEPK.A\tP1\n\
             a/s_2_1\t2\t0.002\t0.3\tK.PEPK.A\tP1\n\
             a/s_3_1\t1\t0.2\t0.5\tK.PEPK.A\tP1\n",
        )
        .unwrap();
        let paths = [p];
        assert_eq!(pooled_targets(&paths, 0.01, None).unwrap().0.len(), 2);
        let (kept, read) = pooled_targets(&paths, 0.01, Some(0.1)).unwrap();
        assert_eq!((kept.len(), read), (1, 3));
    }

    #[test]
    fn rows_link_to_the_group_of_their_bare_sequence() {
        let inf = inference(&[("PEPTIDEK", "sp|P1|A"), ("OTHERK", "sp|P2|B")]);
        let groups = ProteinGroups {
            inf: &inf,
            max_q: 1.0,
        };
        let mut map = pooled(&[("a/scan=1_1_1", 0.001, 0.01), ("a/scan=2_2_1", 0.001, 0.01)]);
        let rows = vec![
            run_row("a", 1, "PEPTIDEK[UNIMOD:737]", true),
            run_row("a", 2, "UNSEENK", true),
        ];
        let mut tally = GroupTally::default();
        let mut c = StudyCounts::default();
        let out = study_rows(rows, &mut map, &groups, &mut tally, &mut c).unwrap();
        let g = out[0].group.as_ref().expect("grouped");
        assert_eq!(out[0].sequence, "PEPTIDEK");
        assert_eq!((g.anchor.as_str(), g.unique), ("sp|P1|A", true));
        assert!(out[1].group.is_none());
        assert_eq!((c.grouped, tally.no_group), (1, 1));
    }

    #[test]
    fn sdrf_channels_map_to_kit_order_with_reference_and_blank() {
        let text = "source name\tcharacteristics[biological replicate]\tcomment[label]\tcomment[data file]\n\
                    GIS\tpooled\tTMT126\tr1.raw\n\
                    S1\t1\tTMT127\tr1.raw\n\
                    x-blank\t1\tTMT131\tr1.raw\n";
        let mut d = crate::sdrf::parse_isobaric_design(text).unwrap();
        d.path = PathBuf::from("d.sdrf.tsv");
        let cols = run_design(Some(&d), "r1", Plex::Tmt6).unwrap();
        assert_eq!(cols.plex, Some(1));
        assert_eq!(cols.channels[0], Some(("GIS".into(), "reference")));
        assert_eq!(cols.channels[1], Some(("S1".into(), "sample")));
        assert_eq!(cols.channels[2], None);
        assert_eq!(cols.channels[5], Some(("x-blank".into(), "blank")));
        let e = run_design(Some(&d), "r2", Plex::Tmt6).unwrap_err();
        assert!(e.contains("no `comment[data file]` row"), "{e}");
        let e = run_design(Some(&d), "r1", Plex::Itraq4).unwrap_err();
        assert!(e.contains("does not have"), "{e}");
        let tmt10 = "source name\tcomment[label]\tcomment[data file]\nS\tTMT131\tr.raw\n";
        let d = crate::sdrf::parse_isobaric_design(tmt10).unwrap();
        let cols = run_design(Some(&d), "r", Plex::Tmt11).unwrap();
        assert_eq!(cols.channels[9], Some(("S".into(), "sample")));
    }

    #[test]
    fn peptidoform_mz_adds_unimod_and_mass_tags() {
        let plain = peptidoform_mz("PEPTIDEK", 2).unwrap();
        let tmt = peptidoform_mz("PEPTIDEK[UNIMOD:737]", 2).unwrap();
        assert!(
            (tmt - plain - 229.162932 / 2.0).abs() < 1e-4,
            "{plain} {tmt}"
        );
        let tag = peptidoform_mz("PEPTIDEK[+229.16293]", 2).unwrap();
        assert!((tag - tmt).abs() < 1e-4);
        assert!(peptidoform_mz("PEPTIDEK[NoSuchMod]", 2).is_none());
    }
}
