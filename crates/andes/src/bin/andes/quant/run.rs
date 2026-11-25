//! Per-run quantification (`andes` search with quant on, and `quant extract`): collect the scans
//! of a run, pick the identifications, quantify them and write the run's quant tables.

mod artifacts;
mod collect;
mod hits;
pub(crate) mod labels;
mod lfq;
mod ms1;
mod scans;
mod settings;

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use model::mass::PROTON;
use model::spectrum::Spectrum;
use model::IsomerClass;
use output::features::{FeatureRecord, MbrColumns};
use output::{LfqRow, PercolatorPsm, QuantId, TmtRow};
use quant::isobaric::{reporter_noise, signal_to_noise, Plex};
use quant::lfq::{feature_score, LfqTarget};
use quant::tdc::picked_qvalues;
use quant::unmix::{unmix_run, ScanInput, Strength, UnmixRun};
use search::candidates::Candidate;
use search::index::SearchIndex;
use search::psm::TopNQueue;

pub(crate) use artifacts::quant_weight;
use artifacts::{write_feature_outputs, write_run_artifacts};
use collect::GateCounts;
pub(crate) use collect::{QuantCollector, RunMs1, RunScans, TmtMbrInput};
use hits::{
    base_record, is_confident, modifications_of, peptidoform_of, proteins_of, quant_id,
    select_hits, select_rank1, sole_protein, Hit,
};
pub(crate) use labels::LabelSettings;
use lfq::{competition_pairs, faims_groups, group_of, lfq_target, max_apex_of, quantify_targets};
use ms1::Ms1Cache;
pub(crate) use scans::{check_unique_stems, run_stem};
use scans::{covers_reporters, isolation_window_record, median_f32, native_id, survey_scan};
pub(crate) use settings::{label_spec, QuantSettings};

const PURITY_ISOTOPE_PPM: f64 = 10.0;

const MS1_RESOLUTION_CHECK_PPM: f64 = 10.0;
const MS1_MIN_ISOTOPE_PARTNER_FRACTION: f64 = 0.5;

const MS2_TOP_PEAKS: usize = 100;

const REPORTER_REGION_MARGIN: f64 = 1.0;

#[allow(clippy::too_many_arguments)]
fn quantify_label(
    collector: &QuantCollector,
    ls: &LabelSettings,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    index: &SearchIndex,
    hits: &[Hit<'_>],
    label_decoys: Option<&HashMap<String, PercolatorPsm>>,
    cutoffs: Option<&output::Cutoffs>,
    parquet_dir: Option<&Path>,
    groups: Option<&crate::quant::ProteinGroups<'_>>,
    next_feature_id: &mut i64,
    dir: &Path,
    stem: &str,
) -> Result<Vec<FeatureRecord>, String> {
    let decoy_hits: Vec<Hit<'_>> = match label_decoys {
        Some(d) => select_rank1(
            collector,
            spectra,
            queues,
            candidates,
            Some(d),
            None,
            cutoffs,
            true,
        )
        .into_iter()
        .filter(|h| h.cand.is_decoy)
        .collect(),
        None => Vec::new(),
    };
    let (rows, label_records, label_meta) = labels::quantify(
        collector,
        ls,
        spectra,
        hits,
        &decoy_hits,
        candidates,
        index,
        parquet_dir.is_some(),
        next_feature_id,
    )?;
    let runs: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
    let channels: Vec<String> = ls.spec.channels.iter().map(|c| c.name.clone()).collect();
    let path = dir.join(format!("{stem}.label.parquet"));
    output::write_label_parquet(
        &path,
        &rows,
        &output::LabelTableMeta {
            channels: &channels,
            reference: ls.spec.reference,
            runs: &runs,
            extra: &label_meta,
        },
    )?;
    log::info!(
        target: "quant",
        "wrote {} ({} channel rows)",
        path.display(),
        rows.len()
    );
    if parquet_dir.is_some() && ls.spec.channels.iter().any(|c| c.qpx_label.is_none()) {
        log::warn!(
            target: "quant",
            "--label: channels outside the SDRF/QPX label vocabulary (SILAC light / \
             medium / heavy, DIMETHYL0/4/6/8) are not written to quantms.feature.parquet"
        );
    }
    let ratios = labels::protein_ratios(
        &rows,
        |seq| {
            let pg = groups?;
            match pg.inf.group_at(seq, pg.max_q) {
                output::protein::inference::GroupAt::Group { group, .. } => {
                    Some((group.accessions.join(";"), group.lead.clone()))
                }
                _ => None,
            }
        },
        labels::MIN_RATIO_COUNT,
    );
    let ppath = dir.join(format!("{stem}.label_proteins.tsv"));
    std::fs::write(&ppath, labels::protein_ratio_tsv(&ratios))
        .map_err(|e| format!("write {}: {e}", ppath.display()))?;
    log::info!(
        target: "quant",
        "wrote {} ({} protein ratios from at least {} peptides)",
        ppath.display(),
        ratios.len(),
        labels::MIN_RATIO_COUNT
    );
    Ok(label_records)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_quant(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    index: &SearchIndex,
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    omitted: &HashSet<String>,
    cutoffs: Option<&output::Cutoffs>,
    report_base: &Path,
    parquet_dir: Option<&Path>,
    groups: Option<&crate::quant::ProteinGroups<'_>>,
    mbr: Option<&TmtMbrInput<'_>>,
    label_decoys: Option<&HashMap<String, PercolatorPsm>>,
) -> Result<output::QuantSummary, String> {
    let s = &collector.settings;
    let mut summary = output::QuantSummary::default();
    let t0 = std::time::Instant::now();
    if let Some(dir) = &s.quant_dir {
        write_run_artifacts(collector, spectra, queues, candidates, index, dir)?;
    }
    let hits = select_hits(
        collector, spectra, queues, candidates, rescore, omitted, cutoffs,
    );
    log::info!(
        target: "quant",
        "{} rank-1 target PSMs selected{} ({})",
        hits.len(),
        match rescore {
            Some(_) => format!(" at q<={}", s.psm_fdr),
            None => " (no rescoring: every rank-1 target PSM; filter downstream)".to_string(),
        },
        s.describe()
    );
    let stem = report_base
        .file_stem()
        .map(|x| x.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let dir: PathBuf = report_base
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| ".".into());
    let hits_by_psm: HashMap<i64, &Hit<'_>> = if parquet_dir.is_some() {
        hits.iter().map(|h| (h.psm_index, h)).collect()
    } else {
        HashMap::new()
    };
    let mut records: Vec<FeatureRecord> = Vec::new();
    let mut next_feature_id: i64 = 1;

    let q = QuantInputs {
        collector,
        spectra,
        candidates,
        index,
        rescore,
        parquet_dir,
        mbr,
        hits: &hits,
        hits_by_psm: &hits_by_psm,
        stem: &stem,
        dir: &dir,
    };
    if let Some(plex) = s.plex {
        quantify_tmt(&q, plex, &mut summary, &mut records, &mut next_feature_id)?;
    }

    if s.lfq {
        quantify_lfq(&q, &mut summary, &mut records, &mut next_feature_id)?;
    }

    if let Some(ls) = &s.label {
        records.extend(quantify_label(
            collector,
            ls,
            spectra,
            queues,
            candidates,
            index,
            &hits,
            label_decoys,
            cutoffs,
            parquet_dir,
            groups,
            &mut next_feature_id,
            &dir,
            &stem,
        )?);
    }

    let quantified = s.plex.is_some() || s.lfq || s.label.is_some() || !records.is_empty();
    if let (Some(pdir), true) = (parquet_dir, quantified) {
        write_feature_outputs(pdir, groups, records)?;
    }
    log::info!(target: "quant", "done in {:.2} s", t0.elapsed().as_secs_f64());
    Ok(summary)
}

/// What the per-mode quantification passes of `run_quant` read: the collector and the search
/// result behind it, the selected hits (and by PSM index), the rescoring targets, the MBR input
/// and where the reports go.
#[derive(Clone, Copy)]
struct QuantInputs<'a> {
    collector: &'a QuantCollector,
    spectra: &'a [Spectrum],
    candidates: &'a [Candidate],
    index: &'a SearchIndex,
    rescore: Option<&'a HashMap<String, PercolatorPsm>>,
    parquet_dir: Option<&'a Path>,
    mbr: Option<&'a TmtMbrInput<'a>>,
    hits: &'a [Hit<'a>],
    hits_by_psm: &'a HashMap<i64, &'a Hit<'a>>,
    stem: &'a str,
    dir: &'a Path,
}

/// TMT reporter quantification: reporter rows per run, the optional reporter transfers, the
/// TMT report and, with a Parquet directory, one feature record per row.
fn quantify_tmt(
    q: &QuantInputs<'_>,
    plex: Plex,
    summary: &mut output::QuantSummary,
    records: &mut Vec<FeatureRecord>,
    next_feature_id: &mut i64,
) -> Result<(), String> {
    let QuantInputs {
        collector,
        spectra,
        candidates,
        index,
        hits,
        stem,
        dir,
        ..
    } = *q;
    let s = &collector.settings;
    let mut rows: Vec<TmtRow> = Vec::new();
    let (mut no_scan, mut out_of_range) = (0usize, 0usize);
    let mut gated = GateCounts::default();
    let mut with_noise = 0usize;
    let mut unmix_strengths: Vec<(String, f64)> = Vec::new();
    for (run_idx, run) in collector.runs.iter().enumerate() {
        let ms1 = run
            .ms1
            .as_ref()
            .map(Ms1Cache::load)
            .transpose()
            .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
        let run_start = rows.len();
        let mut unmix_inputs: Vec<ScanInput> = Vec::new();
        for hit in hits
            .iter()
            .filter(|h| h.run_idx == run_idx && !h.co_isolated)
        {
            let spec = &spectra[hit.spec_idx];
            let context = run.contexts.get(native_id(run, &spec.title));
            let (quant_scan_id, raw, scan_window, noise) = if s.tmt_level >= 3 {
                match run.ms3_reporters.get(native_id(run, &spec.title)) {
                    Some(r) => (r.id.clone(), r.raw.clone(), r.scan_window, r.noise.clone()),
                    None => {
                        no_scan += 1;
                        continue;
                    }
                }
            } else {
                match collector.ms2_reporters.get(hit.spec_idx) {
                    Some(raw) => (
                        spec.title.clone(),
                        raw.clone(),
                        context.and_then(|c| c.scan_window),
                        context
                            .and_then(|c| c.noise.as_ref())
                            .map(|n| reporter_noise(n, plex)),
                    ),
                    None => {
                        no_scan += 1;
                        continue;
                    }
                }
            };
            if !covers_reporters(scan_window, plex, s.tmt_tol) {
                out_of_range += 1;
                continue;
            }
            let purity = ms1.as_ref().and_then(|ms1| {
                let rt = spec.rt_seconds?;
                let scan = survey_scan(ms1, context, rt)?;
                let lo = spec.isolation_lower_offset?;
                let hi = spec.isolation_upper_offset?;
                let pmz = hit.psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
                quant::precursor_purity(
                    ms1.peaks(scan),
                    pmz,
                    spec.isolation_target_mz,
                    hit.psm.charge_used,
                    lo,
                    hi,
                    PURITY_ISOTOPE_PPM,
                )
            });
            let corrected = match &s.correction {
                Some(m) => m.correct(&raw),
                None => raw.clone(),
            };
            let sn = noise.map(|n| signal_to_noise(&corrected, &n));
            with_noise += usize::from(sn.is_some());
            let sn_sum = sn.as_ref().map(|v| v.iter().sum::<f32>());
            let passes = gated.check(s, purity, sn_sum);
            if s.unmix.is_some() {
                let z = hit.psm.charge_used.max(1);
                let window = ms1.as_ref().and_then(|ms1| {
                    isolation_window_record(
                        ms1,
                        context,
                        spec,
                        hit.cand.peptide.mass() / z as f64 + PROTON,
                        z,
                    )
                });
                unmix_inputs.push(ScanInput {
                    charge: z,
                    reporters: corrected.clone(),
                    pic: collector
                        .ms2_pic
                        .get(hit.spec_idx)
                        .map_or(0.0, |&v| v as f64),
                    window,
                    protein: sole_protein(hit, candidates),
                });
            }
            rows.push(TmtRow {
                id: quant_id(hit, run, spec, candidates, index),
                quant_scan_id,
                quant_ms_level: s.tmt_level,
                purity,
                raw,
                corrected,
                sn,
                quant_weight: if passes {
                    quant_weight(purity, hit.pep)
                } else {
                    0.0
                },
                quant_gate_passed: passes,
                interference: None,
                unmixed: None,
                transfer: None,
            });
        }
        if let Some(strength) = s.unmix {
            let UnmixRun {
                model,
                strength: applied,
                scans: unmixed,
            } = unmix_run(&unmix_inputs, strength);
            unmix_strengths.push((run.stem.clone(), applied));
            match &model {
                Some(m) => {
                    let fractions: Vec<f32> =
                        unmixed.iter().flatten().map(|u| u.interference).collect();
                    log::info!(
                        target: "quant",
                        "{}: interference model on {} scans (precursor yield {:.3}/{:.3}/{:.3} \
                         at charge <=2/3/>=4, co-isolated {:.3}, hidden {:.3}, R2(log) {:.2}); \
                         {} scans unmixed, median interference {:.2}, strength {:.2}",
                        run.stem,
                        m.n,
                        m.precursor[0],
                        m.precursor[1],
                        m.precursor[2],
                        m.co_isolated,
                        m.hidden,
                        m.r2_log,
                        fractions.len(),
                        median_f32(&fractions),
                        applied
                    );
                }
                None => log::warn!(
                    target: "quant",
                    "{}: too few scans with an isolation-window record and \
                     reporter ions to fit the interference model ({} scans, {} needed); \
                     its unmix_ columns are empty",
                    run.stem,
                    unmix_inputs.len(),
                    quant::unmix::MIN_FIT_SCANS
                ),
            }
            for (row, u) in rows[run_start..].iter_mut().zip(unmixed) {
                if let Some(u) = u {
                    row.interference = Some(u.interference);
                    row.unmixed = Some(u.reporters);
                }
            }
        }
    }
    if s.tmt_mbr {
        tmt_transfers(q, plex, &mut rows)?;
    }
    summary.columns = plex.channels().iter().map(|c| c.name.to_string()).collect();
    for r in rows.iter().filter(|r| r.quant_gate_passed) {
        for (col, &v) in r.corrected.iter().enumerate() {
            summary
                .values
                .push((r.id.sequence.clone(), r.id.proteins.clone(), col, v as f64));
        }
    }
    let tmt_path = dir.join(format!("{stem}.tmt.parquet"));
    let unmix = s.unmix.map(|strength| output::UnmixMeta {
        mode: match strength {
            Strength::Fixed(_) => "on",
            Strength::Learned => "auto",
        },
        strengths: unmix_strengths.clone(),
    });
    if s.tmt_mbr {
        output::write_tmt_parquet_mbr(
            &tmt_path,
            plex,
            &rows,
            s.correction.is_some(),
            unmix.as_ref(),
        )?;
    } else {
        output::write_tmt_parquet(
            &tmt_path,
            plex,
            &rows,
            s.correction.is_some(),
            unmix.as_ref(),
        )?;
    }
    log::info!(
        target: "quant",
        "wrote {} ({} PSMs with reporter ions passing the quantification gate{}{}{})",
        tmt_path.display(),
        rows.iter().filter(|r| r.quant_gate_passed).count(),
        if no_scan > 0 {
            format!(", {no_scan} without a quantification scan")
        } else {
            String::new()
        },
        if out_of_range > 0 {
            format!(
                ", {out_of_range} whose scan range does not cover the reporter ions \
                 (not quantified)"
            )
        } else {
            String::new()
        },
        gated.describe(s)
    );
    if with_noise == 0 && s.min_sn > 0.0 && !rows.is_empty() {
        log::warn!(
            target: "quant",
            "no reporter scan carries noise values (mzML converted without \
             noise arrays, or ion-trap reporter scans); --tmt-min-sn was not applied and \
             the sn_ columns are empty"
        );
    }
    tmt_records(q, plex, &rows, records, next_feature_id);
    Ok(())
}

/// TMT reporter transfers between runs (`--mbr` with `--tmt`), appended to `rows`.
fn tmt_transfers(q: &QuantInputs<'_>, plex: Plex, rows: &mut Vec<TmtRow>) -> Result<(), String> {
    let QuantInputs {
        collector,
        rescore,
        mbr,
        ..
    } = *q;
    let s = &collector.settings;
    let input = mbr.ok_or("--mbr with --tmt: the rescoring result is missing")?;
    let targets = rescore.ok_or("--mbr with --tmt needs --rescore")?;
    let quant_dir = s
        .quant_dir
        .as_deref()
        .ok_or("--mbr with --tmt needs the per-run tables of --quant-dir")?;
    let params = crate::quant::tmt::mbr::Params {
        plex,
        correction: s.correction.as_ref(),
        library_fdr: input.library_fdr,
        donor_q: s.psm_fdr,
        donor_pep: s.pep_cap,
        accept_q: input.accept_q,
        wrong_log2: input.wrong_log2,
    };
    let stems: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
    let make = || crate::quant::rescore_scores(targets, input.decoys, input.source);
    let (transfers, outcome) = crate::quant::tmt::mbr::run(quant_dir, &stems, &make, &params)?;
    log::info!(
        target: "quant",
        "--mbr transferred {} reporter rows, {} at quant_q<={}{}",
        outcome.candidates,
        outcome.accepted,
        input.accept_q,
        if outcome.model_fitted {
            ""
        } else {
            " (no risk model could be fitted: none accepted)"
        }
    );
    for t in transfers {
        let (sequence, modifications) = crate::quant::parse_peptidoform(&t.peptidoform);
        let file_idx = collector
            .runs
            .iter()
            .position(|r| r.stem == t.run)
            .unwrap_or(0);
        rows.push(TmtRow {
            id: QuantId {
                file_idx,
                run: t.run.clone(),
                spec_id: t.spec_key.clone(),
                scan: t.scan,
                rt_seconds: Some(t.rt),
                sequence,
                peptidoform: t.peptidoform.clone(),
                modifications,
                charge: t.charge,
                proteins: t.proteins.clone(),
                is_decoy: false,
                calculated_mz: t.calc_mz,
                observed_mz: t.precursor_mz,
                q_value: None,
                pep: None,
                missed_cleavages: None,
                psm_index: -1,
            },
            quant_scan_id: t.native_id.clone(),
            quant_ms_level: s.tmt_level,
            purity: t.purity,
            raw: t.raw.clone(),
            corrected: t.corrected.clone(),
            sn: None,
            quant_weight: if t.accepted {
                quant_weight(t.purity, Some(f64::from(t.info.p_wrong)))
            } else {
                0.0
            },
            quant_gate_passed: t.accepted,
            interference: None,
            unmixed: None,
            transfer: Some(t.info.clone()),
        });
    }
    Ok(())
}

/// One feature record per TMT row that is written (transfers only when their gate passed).
fn tmt_records(
    q: &QuantInputs<'_>,
    plex: Plex,
    rows: &[TmtRow],
    records: &mut Vec<FeatureRecord>,
    next_feature_id: &mut i64,
) {
    let QuantInputs {
        collector,
        candidates,
        index,
        parquet_dir,
        hits_by_psm,
        ..
    } = *q;
    let s = &collector.settings;
    for r in rows.iter().filter(|_| parquet_dir.is_some()) {
        let proteins = match &r.transfer {
            Some(_) if !r.quant_gate_passed => continue,
            Some(_) => {
                r.id.proteins
                    .iter()
                    .map(|a| (a.clone(), None, None, None, None))
                    .collect()
            }
            None => proteins_of(&hits_by_psm[&r.id.psm_index].psm, candidates, index),
        };
        let mut rec = base_record(&r.id, proteins, *next_feature_id);
        *next_feature_id += 1;
        let labels: Vec<String> = plex
            .channels()
            .iter()
            .map(|c| format!("{}{}", plex.label_prefix(), c.name))
            .collect();
        rec.intensities = labels
            .iter()
            .cloned()
            .zip(r.corrected.iter().copied())
            .collect();
        if s.correction.is_some() || r.sn.is_some() || r.unmixed.is_some() {
            rec.additional_intensities = labels
                .iter()
                .enumerate()
                .map(|(i, l)| {
                    let mut extra = Vec::new();
                    if s.correction.is_some() {
                        extra.push(("raw".to_string(), r.raw[i]));
                    }
                    if let Some(sn) = &r.sn {
                        extra.push(("sn".to_string(), sn[i]));
                    }
                    if let Some(u) = &r.unmixed {
                        extra.push(("unmixed".to_string(), u[i]));
                    }
                    (l.clone(), extra)
                })
                .collect();
        }
        if let Some(p) = r.purity {
            rec.additional_scores
                .push(("precursor_purity".into(), p as f64, true));
        }
        if let Some(v) = r.sn_sum() {
            rec.additional_scores
                .push(("reporter_sn_sum".into(), v as f64, true));
        }
        if let Some(v) = r.interference {
            rec.additional_scores
                .push(("interference".into(), v as f64, false));
        }
        rec.additional_scores
            .push(("quant_weight".into(), r.quant_weight as f64, true));
        rec.cv_params
            .push(("quant_scan".into(), r.quant_scan_id.clone()));
        rec.cv_params
            .push(("quant_ms_level".into(), r.quant_ms_level.to_string()));
        if s.tmt_mbr {
            match &r.transfer {
                Some(t) => {
                    rec.psm_ids.clear();
                    rec.id_run_file_name = None;
                    rec.additional_scores
                        .push(("p_wrong".into(), f64::from(t.p_wrong), false));
                    rec.additional_scores.push((
                        "noise_floor".into(),
                        f64::from(t.noise_floor),
                        false,
                    ));
                    rec.additional_scores
                        .push(("transfer_angle".into(), f64::from(t.angle), true));
                    rec.additional_scores.push((
                        "donor_batch_cos".into(),
                        f64::from(t.batch_cos),
                        true,
                    ));
                    rec.mbr = Some(MbrColumns {
                        evidence_type: "transfer_ms2".into(),
                        transfer_q: Some(t.quant_q),
                        quant_q: Some(t.quant_q),
                        quantity_se: None,
                        ..MbrColumns::default()
                    });
                }
                None => {
                    rec.mbr = Some(MbrColumns {
                        evidence_type: "identified".into(),
                        ..MbrColumns::default()
                    });
                }
            }
        }
        records.push(rec);
    }
}

/// Label-free MS1 quantification of the selected hits, the LFQ report and its feature records.
fn quantify_lfq(
    q: &QuantInputs<'_>,
    summary: &mut output::QuantSummary,
    records: &mut Vec<FeatureRecord>,
    next_feature_id: &mut i64,
) -> Result<(), String> {
    let QuantInputs {
        collector,
        spectra,
        candidates,
        index,
        parquet_dir,
        hits,
        hits_by_psm,
        stem,
        dir,
        ..
    } = *q;
    let s = &collector.settings;
    let n_iso = s.lfq_params.n_isotopes;
    let better = |h: &Hit<'_>, other: &Hit<'_>| match (h.q_value, other.q_value) {
        (Some(a), Some(b)) if a != b => a < b,
        _ => h.psm.rank_score > other.psm.rank_score,
    };
    type Precursor = (IsomerClass, u8, bool);
    let mut best: HashMap<(usize, Precursor), usize> = HashMap::new();
    let mut form: HashMap<Precursor, usize> = HashMap::new();
    for (i, h) in hits.iter().enumerate() {
        let precursor = (
            IsomerClass::of_peptide(&h.cand.peptide),
            h.psm.charge_used,
            h.cand.is_decoy,
        );
        for slot in [
            best.entry((h.run_idx, precursor.clone())).or_insert(i),
            form.entry(precursor).or_insert(i),
        ] {
            if better(h, &hits[*slot]) {
                *slot = i;
            }
        }
    }
    let params = &s.lfq_params;
    let mut all_rows: Vec<LfqRow> = Vec::new();
    for (run_idx, run) in collector.runs.iter().enumerate() {
        let Some(cache) = run.ms1.as_ref() else {
            continue;
        };
        let ms1 = cache
            .load()
            .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
        let groups = faims_groups(ms1);
        let (mut no_series, mut no_rt) = (0usize, 0usize);
        let mut targets: Vec<(usize, usize, LfqTarget)> = Vec::new();
        let mut form_of: HashMap<usize, usize> = HashMap::new();
        for ((_, precursor), &hi) in best.iter().filter(|((r, _), _)| *r == run_idx) {
            form_of.insert(hi, form[precursor]);
            let h = &hits[hi];
            let spec = &spectra[h.spec_idx];
            let Some(group) = group_of(&groups, run, spec) else {
                no_series += 1;
                continue;
            };
            let Some(rt) = spec.rt_seconds else {
                no_rt += 1;
                continue;
            };
            targets.push((hi, group, lfq_target(h, rt, n_iso)));
        }
        targets.sort_unstable_by_key(|(hi, _, _)| *hi);
        let quantified = quantify_targets(&groups, &targets, params);
        let max_apex = max_apex_of(&quantified);
        let qvals = picked_qvalues(&competition_pairs(&quantified, params, max_apex));
        let mut n_found = 0usize;
        for (k, (hi, feature, decoy)) in quantified.into_iter().enumerate() {
            let Some(feature) = feature else { continue };
            n_found += 1;
            let hit = &hits[hi];
            let spec = &spectra[hit.spec_idx];
            let mut id = quant_id(hit, run, spec, candidates, index);
            let shown = hits[form_of[&hi]].cand;
            id.peptidoform = peptidoform_of(shown);
            id.modifications = modifications_of(shown);
            all_rows.push(LfqRow {
                id,
                score: feature_score(&feature, params.rt_window_s, max_apex),
                feature_q_value: qvals[k],
                decoy_score: decoy
                    .as_ref()
                    .map(|d| feature_score(d, params.rt_window_s, max_apex)),
                feature,
            });
        }
        log::info!(
            target: "quant",
            "{}: {} of {} precursors quantified, {} at feature q<={} and envelope cosine>={}{}{}",
            run.stem,
            n_found,
            targets.len(),
            all_rows
                .iter()
                .filter(|r| r.id.file_idx == run_idx && is_confident(r, s))
                .count(),
            s.feature_fdr,
            s.lfq_params.min_cosine,
            if no_series > 0 {
                format!(
                    ", {no_series} whose MS2 has no FAIMS voltage with an MS1 series \
                     (not quantified)"
                )
            } else {
                String::new()
            },
            if no_rt > 0 {
                format!(", {no_rt} whose MS2 has no retention time (not quantified)")
            } else {
                String::new()
            }
        );
    }
    let runs: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
    let confident: Vec<LfqRow> = all_rows
        .iter()
        .filter(|r| is_confident(r, s))
        .cloned()
        .collect();
    summary.columns = runs.clone();
    for r in &confident {
        summary.values.push((
            r.id.sequence.clone(),
            r.id.proteins.clone(),
            r.id.file_idx,
            r.feature.area,
        ));
    }
    let lfq_path = dir.join(format!("{stem}.lfq.parquet"));
    let accepted: Vec<bool> = all_rows.iter().map(|r| is_confident(r, s)).collect();
    output::write_lfq_parquet(
        &lfq_path,
        &all_rows,
        &accepted,
        &output::LfqTableMeta {
            runs: &runs,
            feature_fdr: s.feature_fdr,
            min_cosine: s.lfq_params.min_cosine,
        },
    )?;
    log::info!(
        target: "quant",
        "wrote {} ({} features, {} accepted at feature q<={}, cosine>={})",
        lfq_path.display(),
        all_rows.len(),
        confident.len(),
        s.feature_fdr,
        s.lfq_params.min_cosine,
    );
    for r in all_rows.iter().filter(|_| parquet_dir.is_some()) {
        let hit = hits_by_psm[&r.id.psm_index];
        let proteins = proteins_of(&hit.psm, candidates, index);
        let mut rec = base_record(&r.id, proteins, *next_feature_id);
        *next_feature_id += 1;
        let f = &r.feature;
        rec.rt = Some(f.apex_rt as f32);
        rec.rt_start = Some(f.rt_start as f32);
        rec.rt_stop = Some(f.rt_stop as f32);
        rec.intensities = vec![("LFQ".to_string(), f.area as f32)];
        let mut extra = vec![("apex".to_string(), f.apex_intensity)];
        for (k, a) in f.isotope_areas.iter().enumerate() {
            extra.push((format!("isotope_{k}_area"), *a as f32));
        }
        rec.additional_intensities = vec![("LFQ".to_string(), extra)];
        rec.additional_scores
            .push(("envelope_cosine".into(), f.cosine as f64, true));
        rec.additional_scores
            .push(("rt_delta_s".into(), f.rt_delta_s, false));
        rec.additional_scores
            .push(("n_isotopes".into(), f.n_isotopes_found as f64, true));
        rec.additional_scores
            .push(("feature_score".into(), r.score as f64, true));
        if let Some(q) = r.feature_q_value {
            rec.additional_scores
                .push(("feature_q_value".into(), q, false));
        }
        if let Some(d) = r.decoy_score {
            rec.additional_scores
                .push(("decoy_feature_score".into(), d as f64, true));
        }
        rec.cv_params.push((
            "psm_rt".into(),
            format!("{:.3}", r.id.rt_seconds.unwrap_or(0.0)),
        ));
        records.push(rec);
    }
    Ok(())
}

#[cfg(test)]
#[path = "run/tests/run.rs"]
mod tests;
