//! `andes train` and `andes train-from-search`: fit the scoring model's tables from labelled
//! PSMs (a TSV, a search, or MSnet Parquet) and update an existing model.

mod args;
mod inputs;
pub(crate) mod intensity;
mod msnet;

use model::activation::ActivationMethod;
use model::{InstrumentType, Tolerance};
use model_train::accumulate::{merge, StatsAccumulator};
use model_train::counts::CountStats;
use model_train::estimate::{Estimator, EstimatorConfig};
use model_train::gate::evaluate_candidate;
use model_train::geometry::{corpus_charge_masses, derive_geometry};
use model_train::labeled::bootstrap_labels;
use model_train::store::{
    commit_update, ensure_single_file_store, update_add, update_decay, update_remove,
    update_reweight, write_all_models_with_sources_and_gbdt_pub, write_all_models_with_sources_pub,
    SourceLedger,
};
use model_train::ModelStore;
use scoring_crate::{Param, RankScorer};

use crate::cli::{GeometryArgs, Protocol};
use crate::models::{build_aa_set, load_param_from_store, load_seed_param, ModelEntryOwned};
use crate::train::intensity::GbdtMode;
pub(crate) use args::{TrainArgs, TrainFromSearchArgs};
use inputs::{build_train_search_params, load_labels_from_tsv, load_spectra_for_train, TrainRun};
pub(crate) use msnet::{read_msnet_parquet, MsnetPsm};

fn inherit_seed_gbdts(trained: &mut Param, seed: &Param) {
    trained.frag_intensity_model = seed.frag_intensity_model.clone();
    trained.rich_ion_model = seed.rich_ion_model.clone();
}

pub(crate) fn run_train_from_search(
    args: TrainFromSearchArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    let t0 = std::time::Instant::now();

    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    if let Some(ref update_model_id) = args.update_model.clone() {
        return run_train_update(args, update_model_id, t0);
    }

    if args.spectra.is_empty() {
        return Err("--spectra is required for initial training".into());
    }
    if !args.labels.is_empty() && args.labels.len() != args.spectra.len() {
        return Err(format!(
            "--labels must be given once per --spectra ({} labels, {} spectra)",
            args.labels.len(),
            args.spectra.len()
        )
        .into());
    }

    let (seed_model_id, seed_param): (String, Param) = load_seed_param(&args.seed_model)?;
    log::info!("seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);
    let activation: Option<ActivationMethod> = args
        .activation
        .as_deref()
        .map(|a| {
            ActivationMethod::from_name(a).ok_or_else(|| {
                format!("unknown --activation '{a}' (expected CID/HCD/ETD/EThcD/UVPD/PQD)")
            })
        })
        .transpose()?;
    let protocol = model::protocol::Protocol::from_name(&args.protocol).ok_or_else(|| {
        format!(
            "unknown --protocol '{}' (expected Automatic/TMT/iTRAQ/iTRAQPhospho/Phosphorylation/Standard)",
            args.protocol
        )
    })?;

    let mut runs: Vec<TrainRun> = Vec::with_capacity(args.spectra.len());
    for (i, spectra_path) in args.spectra.iter().enumerate() {
        log::info!("loading spectra from {} ...", spectra_path.display());
        let mut spectra = load_spectra_for_train(spectra_path)?;
        log::info!("loaded {} spectra", spectra.len());
        if let Some(act) = activation {
            let n = spectra.len();
            spectra.retain(|s| s.activation_method.is_none_or(|a| a == act));
            log::info!(
                "kept {} of {n} spectra with activation {} (or none recorded)",
                spectra.len(),
                act.name()
            );
        }
        let labels = if let Some(labels_path) = args.labels.get(i) {
            let aa_set = build_aa_set(args.fixed_mods.as_ref(), args.variable_mods.as_ref())?;
            load_labels_from_tsv(labels_path, &spectra, &aa_set)?
        } else {
            let database = args
                .database
                .clone()
                .ok_or("--database is required for initial training (or pass --labels)")?;
            let search_params = build_train_search_params(&args)?;
            log::info!("running seed search (train-fdr = {}) ...", args.train_fdr);
            bootstrap_labels(
                &spectra,
                &database,
                &seed_scorer,
                &search_params,
                args.train_fdr,
                &crate::cli::dir_or_output_dir(args.index_cache_dir.as_deref(), &args.out_store),
            )
            .map_err(|e| format!("bootstrap_labels: {e}"))?
        };
        log::info!(
            "{} confident labels from {}",
            labels.len(),
            spectra_path.display()
        );
        runs.push(TrainRun::keep_labelled(spectra_path, spectra, labels));
    }
    let n_labels: usize = runs.iter().map(|r| r.labels.len()).sum();
    log::info!("{n_labels} confident labels from {} run(s)", runs.len());

    if n_labels == 0 {
        return Err(format!(
            "no confident labels found at train-fdr={} — try a higher --train-fdr",
            args.train_fdr
        )
        .into());
    }

    let use_seed_geometry = args.seed_geometry;
    let mut template: Param = if !use_seed_geometry {
        log::info!(
            "deriving own partition geometry from {n_labels} PSMs (pass --seed-geometry to reuse the seed geometry)"
        );
        let corpus: Vec<(i32, f32)> = runs
            .iter()
            .flat_map(|r| corpus_charge_masses(&r.labels))
            .collect();
        let geo_cfg = GeometryArgs::default().to_config();
        let mut base = seed_param.clone();
        if let Some(act) = activation {
            base.data_type.activation = act;
        }
        derive_geometry(&corpus, &base, &geo_cfg)
    } else {
        log::info!("--seed-geometry set — reusing seed partition geometry");
        seed_param.clone()
    };
    if let Some(ppm) = args.fragment_tol_ppm {
        template.mme = Tolerance::Ppm(ppm);
        log::info!("training fragment tolerance overridden to {ppm} ppm");
    } else if let Some(da) = args.fragment_tol_da {
        template.mme = Tolerance::Da(da);
        log::info!("training fragment tolerance overridden to {da} Da");
    } else {
        log::info!("training at the seed fragment tolerance {:?}", template.mme);
    }
    let accum_scorer = RankScorer::new(&template);

    log::info!("accumulating ion-match statistics ...");
    let accumulator = StatsAccumulator::new(&accum_scorer);
    let run_stats: Vec<CountStats> = runs
        .iter()
        .map(|run| {
            let mut stats = CountStats::new();
            for label in &run.labels {
                let spec = &run.spectra[label.spectrum_index];
                accumulator.accumulate(&mut stats, spec, &label.peptide, label.charge);
            }
            stats
        })
        .collect();
    let stats = merge(run_stats.clone());

    log::info!("estimating model parameters ...");
    let cfg = EstimatorConfig::default();
    let estimator = Estimator::new(cfg);
    let mut trained_param = estimator.estimate(&stats, &template);
    inherit_seed_gbdts(&mut trained_param, &seed_param);
    if let Some(act) = activation {
        trained_param.data_type.activation = act;
    }
    trained_param.data_type.protocol = protocol;
    let instrument = match &args.instrument {
        Some(name) => {
            trained_param.data_type.instrument =
                InstrumentType::from_name(name).ok_or_else(|| {
                    format!(
                        "unknown --instrument '{name}' \
                         (expected LowRes/HighRes/QExactive/OrbitrapAstral/TOF/TimsTOF)"
                    )
                })?;
            name.clone()
        }
        None => trained_param.data_type.instrument.name().to_string(),
    };
    log::info!(
        "model key = ({}, {}, {}, {})",
        trained_param.data_type.activation.name(),
        trained_param.data_type.instrument.name(),
        trained_param
            .data_type
            .enzyme
            .as_ref()
            .map_or("any", |e| e.name()),
        trained_param.data_type.protocol.name(),
    );
    let n_partitions = trained_param.partitions.len();
    log::info!("trained model has {} partitions", n_partitions);

    let model_id = args
        .model_id
        .clone()
        .unwrap_or_else(|| format!("trained_{}_{}", instrument, args.protocol));
    log::info!("model ID = {model_id}");

    let date_str = args.date.clone().unwrap_or_else(format_today_iso8601);
    let mut sources: Vec<(SourceLedger, CountStats)> = Vec::with_capacity(runs.len());
    let mut seen_ids = std::collections::HashSet::new();
    for (run, run_stats) in runs.iter().zip(run_stats) {
        let mut source_id = if runs.len() == 1 {
            format!("bootstrap_{model_id}")
        } else {
            run.stem.clone()
        };
        let mut k = 2;
        while !seen_ids.insert(source_id.clone()) {
            source_id = format!("{}_{k}", run.stem);
            k += 1;
        }
        let ledger = SourceLedger {
            source_id,
            dataset: run.dataset.clone(),
            n_psms: run.labels.len() as i64,
            date: date_str.clone(),
            weight: 1.0,
            train_fdr: args.train_fdr as f32,
            decay_half_life_days: None,
            instrument: instrument.clone(),
            experiment_class: protocol.name().to_string(),
        };
        sources.push((ledger, run_stats));
    }
    let source_list = sources
        .iter()
        .map(|(l, _)| l.source_id.as_str())
        .collect::<Vec<_>>()
        .join(", ");

    let store_path = &args.out_store;
    let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
    if store_path.exists() {
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
        for id in store.model_ids() {
            if id == model_id {
                log::info!("overwriting existing model '{id}' in store");
                continue;
            }
            let p = store
                .load_param(&id)
                .map_err(|e| format!("reading model '{id}': {e}"))?;
            let src_ledgers = store.load_sources(&id).unwrap_or_default();
            let mut src = Vec::new();
            for l in src_ledgers {
                if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                    src.push((l, s));
                }
            }
            existing_other.push((id, p, src));
        }
    }

    let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
    all_entries.push((model_id.clone(), trained_param.clone(), sources));
    for (id, p, src) in existing_other {
        all_entries.push((id, p, src));
    }

    write_all_models_with_sources_pub(
        store_path,
        &all_entries
            .iter()
            .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
            .collect::<Vec<_>>(),
    )
    .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;

    log::info!(
        "wrote model '{model_id}' to {} (sources: {source_list}) [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );
    log::info!(
        "summary — labels={}, partitions={}, store={}",
        n_labels,
        n_partitions,
        store_path.display(),
    );

    Ok(())
}

pub(crate) fn run_train(args: TrainArgs) -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;

    let t0 = std::time::Instant::now();

    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    let mut psms: Vec<MsnetPsm> = Vec::new();
    let mut rows_read = 0usize;
    for input in &args.inputs {
        log::info!("reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        rows_read += part.len();
        log::info!("{} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    log::info!(
        "{rows_read} total PSM rows across {} file(s)",
        args.inputs.len()
    );

    let (seed_model_id, mut seed_param): (String, Param) =
        load_seed_param(&Some(args.seed_model.clone()))?;
    log::info!("seed model = {seed_model_id}");
    if let Some(ppm) = args.fragment_tol_ppm {
        seed_param.mme = Tolerance::Ppm(ppm);
        log::info!("fragment tolerance overridden to {ppm} ppm");
    } else if let Some(da) = args.fragment_tol_da {
        seed_param.mme = Tolerance::Da(da);
        log::info!("fragment tolerance overridden to {da} Da");
    } else {
        log::info!("using seed fragment tolerance {:?}", seed_param.mme);
    }

    let use_seed_geometry = args.seed_geometry;
    let template: Param = if !use_seed_geometry {
        log::info!(
            "deriving own partition geometry from {} PSMs (pass --seed-geometry to reuse the seed geometry)",
            psms.len()
        );
        let corpus: Vec<(i32, f32)> = psms
            .iter()
            .map(|p| (p.charge as i32, p.peptide.mass() as f32))
            .collect();
        let geo_cfg = GeometryArgs::default().to_config();
        derive_geometry(&corpus, &seed_param, &geo_cfg)
    } else {
        log::info!("--seed-geometry set — reusing seed partition geometry");
        seed_param.clone()
    };

    let seed_scorer = RankScorer::new(&template);

    log::info!("accumulating ion-match statistics ...");
    let stats = psms
        .par_iter()
        .fold(CountStats::new, |mut acc, psm| {
            let accumulator = StatsAccumulator::new(&seed_scorer);
            accumulator.accumulate(&mut acc, &psm.spectrum, &psm.peptide, psm.charge);
            acc
        })
        .collect::<Vec<_>>();
    let stats = merge(stats);
    log::info!("accumulated {} PSMs", psms.len());

    log::info!("estimating model parameters ...");
    let cfg = EstimatorConfig {
        pseudo: args.train_pseudo,
        noise_pseudo: args.train_noise_pseudo,
        min_count: args.train_min_count,
        backoff_weight: args.train_backoff_weight,
        error_scaling_factor_override: None,
        rank_smoothing: args.rank_smoothing,
    };
    log::info!(
        "estimator pseudo={} noise_pseudo={} backoff_weight={} min_count={}",
        cfg.pseudo,
        cfg.noise_pseudo,
        cfg.backoff_weight,
        cfg.min_count
    );
    let estimator = Estimator::new(cfg);

    let prior_param: Option<Param> = match &args.prior_model_store {
        Some(store_path) => {
            let prior_id = args
                .prior_model
                .clone()
                .unwrap_or_else(|| args.model_id.clone());
            let (_pid, p) = load_param_from_store(
                seed_param.data_type.activation,
                Some(seed_param.data_type.instrument),
                Protocol::Auto,
                model::enzyme::Enzyme::Trypsin,
                Some(store_path.as_path()),
                Some(&prior_id),
            )
            .map_err(|e| format!("loading --prior-model '{prior_id}': {e}"))?;
            log::info!("prior model = {prior_id} (from {})", store_path.display());
            Some(p)
        }
        None => None,
    };

    let mut trained_param = estimator.estimate_with_prior(&stats, &template, prior_param.as_ref());
    inherit_seed_gbdts(&mut trained_param, &seed_param);
    let n_partitions = trained_param.partitions.len();
    log::info!("trained model has {n_partitions} partitions");

    if let Some(act) = &args.activation {
        trained_param.data_type.activation = ActivationMethod::from_name(act).ok_or_else(|| {
            format!("unknown --activation '{act}' (expected CID/HCD/ETD/UVPD/PQD)")
        })?;
    }
    if let Some(inst) = &args.instrument {
        trained_param.data_type.instrument = InstrumentType::from_name(inst).ok_or_else(|| {
            format!("unknown --instrument '{inst}' (expected LowRes/HighRes/QExactive/TOF)")
        })?;
    }
    if let Some(enz) = &args.enzyme {
        trained_param.data_type.enzyme =
            Some(model::enzyme::Enzyme::from_name(enz).ok_or_else(|| {
                format!("unknown --enzyme '{enz}' (e.g. Trypsin/LysC/LysN/AspN/GluC/ArgC)")
            })?);
    }
    if let Some(prot) = &args.protocol {
        trained_param.data_type.protocol = model::protocol::Protocol::from_name(prot)
            .ok_or_else(|| format!("unknown --protocol '{prot}' (expected Automatic/TMT/iTRAQ/iTRAQPhospho/Phosphorylation/Standard)"))?;
    }
    log::info!(
        "model data_type = {:?}/{:?}/{:?}/{:?}",
        trained_param.data_type.activation,
        trained_param.data_type.instrument,
        trained_param.data_type.enzyme,
        trained_param.data_type.protocol,
    );

    let dataset = args
        .inputs
        .first()
        .and_then(|p| p.file_name())
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "msnet".to_string());
    let ledger = SourceLedger {
        source_id: args.source.clone(),
        dataset,
        n_psms: psms.len() as i64,
        date: format_today_iso8601(),
        weight: 1.0,
        train_fdr: 1.0,
        decay_half_life_days: None,
        instrument: trained_param.data_type.instrument.name().to_string(),
        experiment_class: trained_param.data_type.protocol.name().to_string(),
    };

    let gbdt_blob: Option<Vec<u8>> = if args.gbdt != GbdtMode::Off {
        use model_train::gbdt::dataset::{build_dataset, PsmRow};
        use model_train::gbdt::train::{train_gbdt, TrainParams};

        log::info!("building GBDT dataset from {} PSMs ...", psms.len());
        let gbdt_rows: Vec<PsmRow<'_>> = psms
            .iter()
            .map(|psm| PsmRow {
                spectrum: &psm.spectrum,
                peptide: &psm.peptide,
                charge: psm.charge,
            })
            .collect();

        let gbdt_dataset = build_dataset(&gbdt_rows, &seed_scorer);
        log::info!(
            "GBDT dataset: {} peak rows, {} positives",
            gbdt_dataset.y.len(),
            gbdt_dataset.y.iter().filter(|&&l| l == 1).count(),
        );

        let gbdt_params = TrainParams {
            allow_degenerate: args.allow_degenerate_model,
            ..TrainParams::default()
        };
        let trained_gbdt = train_gbdt(&gbdt_dataset, &gbdt_params, 42)?;
        log::info!("GBDT trained: {} trees", trained_gbdt.trees.len(),);
        Some(trained_gbdt.to_bytes())
    } else {
        None
    };

    let store_path = &args.out_store;
    let model_id = args.model_id.clone();

    {
        let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
        let mut existing_blobs: Vec<Option<Vec<u8>>> = Vec::new();
        if store_path.exists() {
            let store = ModelStore::open(store_path)
                .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
            for id in store.model_ids() {
                if id == model_id {
                    log::info!("overwriting existing model '{id}' in store");
                    continue;
                }
                let p = store
                    .load_param(&id)
                    .map_err(|e| format!("reading model '{id}': {e}"))?;
                let blob = p.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
                let src_ledgers = store.load_sources(&id).unwrap_or_default();
                let mut src = Vec::new();
                for l in src_ledgers {
                    if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                        src.push((l, s));
                    }
                }
                existing_other.push((id, p, src));
                existing_blobs.push(blob);
            }
        }

        let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
        all_entries.push((model_id.clone(), trained_param, vec![(ledger, stats)]));
        for (id, p, src) in existing_other {
            all_entries.push((id, p, src));
        }

        let mut all_blobs: Vec<Option<Vec<u8>>> = vec![gbdt_blob];
        all_blobs.extend(existing_blobs);

        write_all_models_with_sources_and_gbdt_pub(
            store_path,
            &all_entries
                .iter()
                .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
                .collect::<Vec<_>>(),
            &all_blobs,
        )
        .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    }

    log::info!(
        "wrote model '{model_id}' to {} (source '{}', {} PSMs, {n_partitions} partitions) [{:.2}s]",
        store_path.display(),
        args.source,
        psms.len(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

pub(crate) fn run_train_update(
    args: TrainFromSearchArgs,
    model_id: &str,
    t0: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    if args.fragment_tol_ppm.is_some() || args.fragment_tol_da.is_some() {
        return Err(
            "--fragment-tol-ppm/--fragment-tol-da apply to initial training only; \
             an update keeps the stored model's tolerance"
                .into(),
        );
    }
    let store_path = &args.out_store;
    let cfg = EstimatorConfig::default();
    ensure_single_file_store(store_path)?;

    let (candidate, new_sources) = if args.add {
        let spectra_path = match args.spectra.as_slice() {
            [one] => one.clone(),
            [] => return Err("--spectra is required with --add".into()),
            _ => return Err("--update --add takes a single --spectra file".into()),
        };
        let database = args
            .database
            .clone()
            .ok_or("--database is required with --add")?;
        let source_id = args
            .source_id
            .clone()
            .ok_or("--source-id is required with --add")?;

        log::info!(
            target: "train-update",
            "loading spectra from {} ...",
            spectra_path.display()
        );
        let spectra = load_spectra_for_train(&spectra_path)?;
        log::info!(target: "train-update", "loaded {} spectra", spectra.len());

        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening store {}: {e}", store_path.display()))?;
        let current_param = store
            .load_param(model_id)
            .map_err(|e| format!("loading model '{model_id}': {e}"))?;
        let current_scorer = RankScorer::new(&current_param);

        let search_params = build_train_search_params(&args)?;

        log::info!(
            target: "train-update",
            "running seed search (train-fdr={}) ...",
            args.train_fdr
        );
        let labels = bootstrap_labels(
            &spectra,
            &database,
            &current_scorer,
            &search_params,
            args.train_fdr,
            &crate::cli::dir_or_output_dir(args.index_cache_dir.as_deref(), &args.out_store),
        )
        .map_err(|e| format!("bootstrap_labels: {e}"))?;
        log::info!(target: "train-update", "{} confident labels", labels.len());

        if labels.is_empty() {
            return Err(format!(
                "no confident labels at train-fdr={} — try a higher --train-fdr",
                args.train_fdr
            )
            .into());
        }

        let accumulator = StatsAccumulator::new(&current_scorer);
        let mut stats = CountStats::new();
        for label in &labels {
            accumulator.accumulate(
                &mut stats,
                &spectra[label.spectrum_index],
                &label.peptide,
                label.charge,
            );
        }
        let stats = merge(vec![stats]);

        let spectra_filename = spectra_path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| spectra_path.display().to_string());
        let date_str = args.date.clone().unwrap_or_else(format_today_iso8601);
        let ledger = SourceLedger {
            source_id: source_id.clone(),
            dataset: spectra_filename,
            n_psms: labels.len() as i64,
            date: date_str,
            weight: 1.0,
            train_fdr: args.train_fdr as f32,
            decay_half_life_days: None,
            instrument: args
                .instrument
                .clone()
                .unwrap_or_else(|| current_param.data_type.instrument.name().to_string()),
            experiment_class: current_param.data_type.protocol.name().to_string(),
        };

        update_add(store_path, model_id, ledger, stats, cfg)
            .map_err(|e| format!("update_add: {e}"))?
    } else if let Some(ref sid) = args.remove_source.clone() {
        update_remove(store_path, model_id, sid, cfg).map_err(|e| format!("update_remove: {e}"))?
    } else if let Some(ref spec) = args.reweight.clone() {
        let (sid, weight) = parse_reweight_spec(spec)?;
        update_reweight(store_path, model_id, &sid, weight, cfg)
            .map_err(|e| format!("update_reweight: {e}"))?
    } else if let Some(half_life) = args.decay {
        update_decay(store_path, model_id, half_life, cfg)
            .map_err(|e| format!("update_decay: {e}"))?
    } else {
        return Err(
            "update mode requires one of: --add, --remove-source, --reweight, --decay".into(),
        );
    };

    let commit = if let Some(ref validate_path) = args.validate.clone() {
        let database = args
            .database
            .clone()
            .ok_or("--database is required with --validate")?;

        log::info!(
            target: "train-update",
            "running acceptance gate on {} ...",
            validate_path.display()
        );
        let val_spectra = load_spectra_for_train(validate_path)?;

        let store =
            ModelStore::open(store_path).map_err(|e| format!("opening store for gate: {e}"))?;
        let current_param = store
            .load_param(model_id)
            .map_err(|e| format!("loading current model for gate: {e}"))?;
        let current_scorer = RankScorer::new(&current_param);
        let candidate_scorer = RankScorer::new(&candidate);

        let search_params = build_train_search_params(&args)?;

        let delta = evaluate_candidate(
            &val_spectra,
            &database,
            &current_scorer,
            &candidate_scorer,
            &search_params,
            args.train_fdr,
            &crate::cli::dir_or_output_dir(args.index_cache_dir.as_deref(), &args.out_store),
        )
        .map_err(|e| format!("evaluate_candidate: {e}"))?;

        log::info!(
            target: "train-update",
            "gate — current={} PSMs, candidate={} PSMs at FDR={}",
            delta.current_count, delta.candidate_count, args.train_fdr
        );

        if delta.is_accepted() {
            log::info!(target: "train-update", "ACCEPTED (candidate >= current)");
            true
        } else {
            log::warn!(target: "train-update", "REJECTED (candidate < current)");
            if args.force {
                log::warn!(target: "train-update", "--force set, committing anyway");
                true
            } else {
                log::info!(target: "train-update", "skipping commit (use --force to override)");
                false
            }
        }
    } else {
        log::warn!(target: "train-update", "no --validate dataset; skipping acceptance gate");
        if args.force {
            log::info!(target: "train-update", "--force set, committing unconditionally");
        }
        true
    };

    if commit {
        commit_update(store_path, model_id, &candidate, &new_sources)
            .map_err(|e| format!("commit_update: {e}"))?;
        log::info!(
            target: "train-update",
            "committed model '{model_id}' to {} [{:.2}s]",
            store_path.display(),
            t0.elapsed().as_secs_f64(),
        );
    }

    Ok(())
}

pub(crate) fn parse_reweight_spec(spec: &str) -> Result<(String, f32), Box<dyn std::error::Error>> {
    let pos = spec
        .rfind('=')
        .ok_or_else(|| format!("--reweight value must be <source-id>=<weight>, got '{spec}'"))?;
    let sid = spec[..pos].to_string();
    let weight: f32 = spec[pos + 1..]
        .parse()
        .map_err(|e| format!("invalid weight in --reweight '{spec}': {e}"))?;
    Ok((sid, weight))
}

pub(crate) fn format_today_iso8601() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = secs / 86400;
    unix_days_to_iso8601(days)
}

pub(crate) fn unix_days_to_iso8601(days: u64) -> String {
    let jdn = days as i64 + 2_440_588;
    let a = jdn + 32044;
    let b = (4 * a + 3) / 146097;
    let c = a - (146097 * b) / 4;
    let d = (4 * c + 3) / 1461;
    let e = c - (1461 * d) / 4;
    let m = (5 * e + 2) / 153;
    let day = e - (153 * m + 2) / 5 + 1;
    let month = m + 3 - 12 * (m / 10);
    let year = 100 * b + d - 4800 + m / 10;
    format!("{:04}-{:02}-{:02}", year, month, day)
}

#[cfg(test)]
#[path = "train/tests/msnet.rs"]
mod train_from_msnet_tests;
