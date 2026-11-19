use std::fs::File;
use std::path::{Path, PathBuf};

use crate::cli::{parse_unit_fraction, GeometryArgs, Protocol};
use crate::model_select::{build_aa_set, load_param_from_store, load_seed_param, ModelEntryOwned};
use crate::spectra::spectrum_ext_lower;
use crate::train_intensity::GbdtMode;
use clap::Args;
use input::{MgfReader, MzMLReader};
use model::{
    activation::ActivationMethod, InstrumentType, ModLocation, Modification,
    ResidueSpec, Spectrum, Tolerance,
};
use model_train::{
    accumulate::{merge, StatsAccumulator},
    counts::CountStats,
    estimate::{Estimator, EstimatorConfig},
    gate::evaluate_candidate,
    geometry::{corpus_charge_masses, derive_geometry},
    labeled::bootstrap_labels,
    store::{
        commit_update, ensure_single_file_store, update_add, update_decay, update_remove,
        update_reweight, write_all_models_with_sources_and_gbdt_pub,
        write_all_models_with_sources_pub, SourceLedger,
    },
    ModelStore,
};
use scoring_crate::{Param, RankScorer};
use search::SearchParams;

fn inherit_seed_gbdts(trained: &mut Param, seed: &Param) {
    trained.frag_intensity_model = seed.frag_intensity_model.clone();
    trained.rich_ion_model = seed.rich_ion_model.clone();
}

/// Training arguments for `andes train-from-search`.
#[derive(Args, Debug)]
pub(crate) struct TrainFromSearchArgs {
    /// Reuse the seed model's geometry instead of deriving one from the corpus.
    /// Deriving (the default) fits segments and mass tiers to the training data.
    #[arg(long = "seed-geometry", default_value_t = false)]
    pub(crate) seed_geometry: bool,
    /// Input spectrum file (training data). Same format dispatch as for search:
    /// `.mzML`/`.mzml` → mzML reader; anything else → MGF reader.
    ///
    /// Repeat the flag to train one model from several runs. Each run is
    /// searched and labelled on its own (its own target-decoy q-values), keeps
    /// only its labelled spectra in memory, and becomes its own source in the
    /// store ledger (source id = file stem), so it can later be removed or
    /// reweighted with `--update`. `--labels` and `--update --add` take a
    /// single file.
    ///
    /// Required for initial training.  In `--update` mode with `--remove-source`
    /// or `--reweight` / `--decay`, `--spectra` is only required when
    /// `--validate` is also given (to run the acceptance gate).
    #[arg(long)]
    pub(crate) spectra: Vec<PathBuf>,

    /// Input FASTA target database (decoys are generated automatically).
    ///
    /// Required for initial training and for `--update --add`.
    /// In `--update` mode without `--add`, only required when `--validate` is
    /// given.
    #[arg(long)]
    pub(crate) database: Option<PathBuf>,

    /// EXTERNAL LABELS (SP-B glyco training): a TSV with columns `scan`,
    /// `peptide`, `charge`. When given, labels come from this file instead of a
    /// seed search (`bootstrap_labels`), and `--database` is not required. The
    /// `peptide` is the BARE backbone sequence (the glycan is stripped — the
    /// rank model scores backbone b/y ions, which in HCD come from the
    /// deglycosylated backbone); Cam-C and any `--fixed-mods` / `--variable-mods` are applied via the
    /// AminoAcidSet. Rows whose scan is absent from `--spectra` or whose peptide
    /// fails to parse are skipped. Use to train a glyco-regime rank model from
    /// MSFragger/pGlyco backbone IDs.
    ///
    /// With several `--spectra`, give one `--labels` per run, in the same order.
    #[arg(long = "labels")]
    pub(crate) labels: Vec<PathBuf>,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// When omitted, the bundled `hcd_qexactive_tryp` model is used as the seed.
    #[arg(long = "seed-model")]
    pub(crate) seed_model: Option<String>,

    /// Target-decoy q-value threshold for accepting PSMs as confident training
    /// labels. Use a lenient value (e.g. 0.1 or 0.5) for small fixtures.
    #[arg(long = "train-fdr", default_value = "0.01", value_parser = parse_unit_fraction)]
    pub(crate) train_fdr: f64,

    /// Instrument of the trained model (`LowRes`, `HighRes`, `QExactive`,
    /// `OrbitrapAstral`, `TOF`, `TimsTOF`). It is part of the model's selection
    /// key, so a search on data from this instrument picks the model. Default:
    /// the seed model's instrument (`QExactive` for the default seed).
    #[arg(long)]
    pub(crate) instrument: Option<String>,

    /// Activation of the trained model (`CID`, `HCD`, `ETD`, `EThcD`, `UVPD`, `PQD`),
    /// part of its selection key. It also sets the ion series of the derived geometry
    /// (EThcD: b, y, c and z•), and only spectra recorded with this activation are
    /// searched and trained on (spectra without a recorded activation are kept).
    /// Default: the seed model's activation, all spectra.
    #[arg(long, conflicts_with = "seed_geometry")]
    pub(crate) activation: Option<String>,

    /// Experiment-class / protocol tag (e.g. `Automatic`, `TMT`). Default: `Automatic`.
    #[arg(long, default_value = "Automatic")]
    pub(crate) protocol: String,

    /// Path to the Parquet model store to write (created if absent, appended
    /// otherwise). REQUIRED.
    #[arg(long = "out-store")]
    pub(crate) out_store: PathBuf,

    /// Directory of the candidate-index cache of the seed searches, created if
    /// missing. Default: the directory of `--out-store`. A directory that
    /// cannot be written is an error.
    #[arg(long = "index-cache-dir", value_name = "DIR")]
    pub(crate) index_cache_dir: Option<PathBuf>,

    /// Model ID written into the store. Default: `trained_<instrument>_<protocol>`.
    #[arg(long = "model-id")]
    pub(crate) model_id: Option<String>,

    /// Fixed modifications, a list of Name (sites) as for the search's
    /// `--fixed-mods` [default: Carbamidomethyl (C)].
    #[arg(long, value_name = "LIST", value_parser = crate::cli::parse_fixed_mods)]
    pub(crate) fixed_mods: Option<crate::cli::ModList>,

    /// Variable modifications, a list of Name (sites) as for the search's
    /// `--variable-mods` [default: Oxidation (M)].
    #[arg(long, value_name = "LIST", value_parser = crate::cli::parse_variable_mods)]
    pub(crate) variable_mods: Option<crate::cli::ModList>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    pub(crate) threads: usize,

    /// ISO 8601 date string (e.g. `2026-01-01`) recorded in the source ledger.
    /// When omitted, the current date is used for initial training; empty string
    /// is stored when `--date ""` is explicitly passed.
    #[arg(long)]
    pub(crate) date: Option<String>,

    /// Fragment match tolerance in ppm for initial training: the window in
    /// which training matches fragments, stored as the trained model's `mme`.
    /// The seed search that labels the PSMs keeps the seed's `mme`. When
    /// neither this nor `--fragment-tol-da` is given, the seed's `mme` is
    /// kept. Training reconstructs fragment m/z from nominal masses, so a ppm
    /// window this narrow marks most fragments as missing; keep the seed `mme`
    /// unless the comparison is the point.
    #[arg(long = "fragment-tol-ppm", conflicts_with = "fragment_tol_da")]
    pub(crate) fragment_tol_ppm: Option<f64>,

    /// Fragment match tolerance in Da for initial training. Overwrites the seed
    /// model's `mme` (see `--fragment-tol-ppm`).
    #[arg(long = "fragment-tol-da")]
    pub(crate) fragment_tol_da: Option<f64>,

    /// Switch to incremental update mode for this model ID.
    /// When set, one of `--add`, `--remove-source`, `--reweight`, or `--decay`
    /// must be provided.
    #[arg(long = "update", value_name = "MODEL_ID")]
    pub(crate) update_model: Option<String>,

    /// (Update mode) Add a new source from `--spectra`.
    /// Requires `--source-id` and `--database`.
    #[arg(long, requires = "update_model")]
    pub(crate) add: bool,

    /// (Update mode) Source identifier for the new source being added
    /// (used with `--add`).
    #[arg(long = "source-id", requires = "add", value_name = "ID")]
    pub(crate) source_id: Option<String>,

    /// (Update mode) Remove the source with this ID from the model.
    #[arg(long = "remove-source", requires = "update_model", value_name = "ID")]
    pub(crate) remove_source: Option<String>,

    /// (Update mode) Set a source's weight.  Format: `<source-id>=<weight>`,
    /// e.g. `--reweight s0=0.5`.
    #[arg(long = "reweight", requires = "update_model", value_name = "ID=W")]
    pub(crate) reweight: Option<String>,

    /// (Update mode) Apply exponential age-decay to all sources with this
    /// half-life in days.
    #[arg(long = "decay", requires = "update_model", value_name = "DAYS")]
    pub(crate) decay: Option<f32>,

    /// (Update mode) Held-out validation spectra for the acceptance gate.
    /// When omitted the gate is skipped (a warning is printed).
    #[arg(long = "validate", requires = "update_model")]
    pub(crate) validate: Option<PathBuf>,

    /// (Update mode) Commit the update even if the acceptance gate fails.
    #[arg(long, requires = "update_model")]
    pub(crate) force: bool,
}

/// Training arguments for `andes train`.
///
/// Trains a scoring model directly from externally-labeled, high-confidence
/// PSMs supplied as a "flat training parquet" (one row per PSM, each carrying
/// the spectrum peaks + identified peptide + resolved mod mass-deltas). This
/// bypasses the bootstrap-search label step entirely: every input row is a
/// label. The seed model supplies only structural hyperparameters (`mme`,
/// deconvolution, segments, frag/precursor offset tables, max_rank); all
/// learned distributions come from the input data.
#[derive(Args, Debug)]
pub(crate) struct TrainArgs {
    /// Reuse the seed model's geometry instead of deriving one from the corpus.
    /// Deriving (the default) fits segments and mass tiers to the training data.
    #[arg(long = "seed-geometry", default_value_t = false)]
    pub(crate) seed_geometry: bool,
    /// Input flat training parquet(s). Repeatable; stats accumulate across all
    /// inputs into a single model.
    #[arg(long = "in", required = true)]
    pub(crate) inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store")]
    pub(crate) out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    pub(crate) model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// Supplies structural hyperparameters only.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    pub(crate) seed_model: String,

    /// Override the trained model's activation in the store `data_type`
    /// (e.g. `CID`, `HCD`, `ETD`, `UVPD`). Defaults to the seed's value.
    /// Together with `--instrument/--enzyme/--protocol` this lets a new slug
    /// carry the correct selection columns even when seeded from a related model.
    #[arg(long = "activation")]
    pub(crate) activation: Option<String>,

    /// Override the trained model's instrument in the store `data_type`
    /// (e.g. `LowRes`, `HighRes`, `QExactive`, `TOF`). Defaults to the seed's value.
    #[arg(long = "instrument")]
    pub(crate) instrument: Option<String>,

    /// Override the trained model's enzyme in the store `data_type`
    /// (e.g. `Trypsin`, `LysC`, `LysN`). Defaults to the seed's value.
    #[arg(long = "enzyme")]
    pub(crate) enzyme: Option<String>,

    /// Override the trained model's protocol in the store `data_type`
    /// (e.g. `TMT`, `iTRAQ`, `Phosphorylation`, `Automatic`). Drives
    /// `experiment_class` model selection. Defaults to the seed's value.
    #[arg(long = "protocol")]
    pub(crate) protocol: Option<String>,

    /// Fragment match tolerance in ppm. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-da`. When
    /// neither is given, the seed model's `mme` is kept.
    #[arg(long = "fragment-tol-ppm", conflicts_with = "fragment_tol_da")]
    pub(crate) fragment_tol_ppm: Option<f64>,

    /// Fragment match tolerance in Da. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da")]
    pub(crate) fragment_tol_da: Option<f64>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    pub(crate) threads: usize,

    /// Laplace pseudo-count for rank/error tables (lower = sharper; default 1.0).
    #[arg(long = "train-pseudo", default_value_t = 1.0)]
    pub(crate) train_pseudo: f32,

    /// Laplace pseudo-count for the NOISE rank distribution (lower = sharper).
    /// Noise is abundant and concentrated, so it needs far less smoothing than
    /// signal ions; the signal `--train-pseudo` over-flattens it. Default 0.05.
    #[arg(long = "train-noise-pseudo", default_value_t = 0.05)]
    pub(crate) train_noise_pseudo: f32,

    /// Partition backoff prior weight (lower = less smoothing toward parent; default 20).
    #[arg(long = "train-backoff-weight", default_value_t = 20.0)]
    pub(crate) train_backoff_weight: f32,

    /// Minimum partition count before backoff blending (default 50).
    #[arg(long = "train-min-count", default_value_t = 50)]
    pub(crate) train_min_count: u64,

    /// Optional path to an independent prior model store. Sparse partitions in
    /// the trained model shrink toward the matching prior model instead of the
    /// corpus-internal pool. Must be own-data (NOT a bundled seed model) to stay
    /// relicense-safe.
    #[arg(long)]
    pub(crate) prior_model_store: Option<PathBuf>,

    /// Model id to load from `--prior-model-store` (defaults to the trained
    /// model id when omitted).
    #[arg(long)]
    pub(crate) prior_model: Option<String>,

    /// Apply widening rank-window smoothing to signal rank distributions
    /// (Kim et al., Nat Commun 5:5277, 2014).
    #[arg(long)]
    pub(crate) rank_smoothing: bool,

    /// Source identifier for the source ledger. Defaults to "msnet".
    #[arg(long, default_value = "msnet")]
    pub(crate) source: String,

    /// Whether to also train and embed a GBDT peak model. `on` (default) trains
    /// GBDT and writes the blob; `off` writes rank-core only (byte-identical to
    /// the pre-GBDT path).
    #[arg(long, default_value = "on")]
    pub(crate) gbdt: GbdtMode,

    /// Opt-in fallback (finding 3.6): downgrade a failed GBDT quality gate to a
    /// warning and embed the degenerate model anyway. Default off.
    #[arg(long = "allow-degenerate-model", hide = true, default_value_t = false)]
    pub(crate) allow_degenerate_model: bool,
}

/// Spectra a reader returned as errors: how many, and the first error text.
#[derive(Default)]
struct ParseFailures {
    count: usize,
    first: Option<String>,
}

impl ParseFailures {
    fn add(&mut self, e: impl std::fmt::Display) {
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
    aa_set: &model::aa_set::AminoAcidSet,
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
        terminal.iter().filter(|(l, _)| *l == loc).map(|(_, d)| d).sum()
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

struct TrainRun {
    dataset: String,
    stem: String,
    spectra: Vec<Spectrum>,
    labels: Vec<model_train::labeled::LabeledMatch>,
}

impl TrainRun {
    fn keep_labelled(
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

fn run_stem(path: &Path) -> String {
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
            log::info!(
                "running seed search (train-fdr = {}) ...",
                args.train_fdr
            );
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
    log::info!(
        "{n_labels} confident labels from {} run(s)",
        runs.len()
    );

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
        log::info!(
            "training at the seed fragment tolerance {:?}",
            template.mme
        );
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

    let date_str = args.date.clone().unwrap_or_else(|| {
        format_today_iso8601()
    });
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
            experiment_class: args.protocol.clone(),
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

pub(crate) struct MsnetPsm {
    pub(crate) spectrum: Spectrum,
    pub(crate) peptide: model::peptide::Peptide,
    pub(crate) charge: u8,
}

pub(crate) fn build_msnet_peptide(
    seq: &str,
    res_mod_pos: &[i32],
    res_mod_delta: &[f64],
    nterm_delta: f64,
    cterm_delta: f64,
) -> Result<model::peptide::Peptide, Box<dyn std::error::Error>> {
    if res_mod_pos.len() != res_mod_delta.len() {
        return Err(format!(
            "res_mod_pos ({}) and res_mod_delta ({}) length mismatch",
            res_mod_pos.len(),
            res_mod_delta.len()
        )
        .into());
    }
    let bytes = seq.as_bytes();
    if bytes.is_empty() {
        return Err("empty peptide sequence".into());
    }
    let n = bytes.len();

    let mut residue_delta = vec![0.0f64; n];
    let mut residue_modded = vec![false; n];
    for (&pos1, &delta) in res_mod_pos.iter().zip(res_mod_delta.iter()) {
        if pos1 < 1 || (pos1 as usize) > n {
            return Err(
                format!("res_mod_pos {pos1} out of range for sequence of length {n}").into(),
            );
        }
        let idx = (pos1 - 1) as usize;
        residue_delta[idx] += delta;
        residue_modded[idx] = true;
    }
    if nterm_delta != 0.0 {
        residue_delta[0] += nterm_delta;
        residue_modded[0] = true;
    }
    if cterm_delta != 0.0 {
        residue_delta[n - 1] += cterm_delta;
        residue_modded[n - 1] = true;
    }

    let mut residues = Vec::with_capacity(n);
    for (i, &r) in bytes.iter().enumerate() {
        let aa = model::AminoAcid::standard(r)
            .ok_or_else(|| format!("non-standard residue {:?} at position {}", r as char, i + 1))?;
        if residue_modded[i] {
            let m = Modification {
                name: "msnet".to_string(),
                mass_delta: residue_delta[i],
                residue: ResidueSpec::Specific(r),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            };
            residues.push(aa.with_mod(m));
        } else {
            residues.push(aa);
        }
    }

    Ok(model::peptide::Peptide::new(residues, b'_', b'-'))
}

pub(crate) fn read_msnet_parquet(path: &Path) -> Result<Vec<MsnetPsm>, Box<dyn std::error::Error>> {
    use arrow::array::{Array, Float32Array, Float64Array, Int32Array, ListArray, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = File::open(path).map_err(|e| format!("opening {}: {e}", path.display()))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| format!("parquet reader for {}: {e}", path.display()))?
        .build()
        .map_err(|e| format!("building parquet reader for {}: {e}", path.display()))?;

    let mut out = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| format!("reading batch from {}: {e}", path.display()))?;

        let col = |name: &str| -> Result<&std::sync::Arc<dyn Array>, Box<dyn std::error::Error>> {
            batch
                .column_by_name(name)
                .ok_or_else(|| format!("missing column '{name}' in {}", path.display()).into())
        };

        let seq = col("seq")?
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or("column 'seq' is not a STRING column")?;
        let charge = col("charge")?
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or("column 'charge' is not an INT32 column")?;
        let prec_mz = col("prec_mz")?
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or("column 'prec_mz' is not a DOUBLE column")?;
        let res_mod_pos = col("res_mod_pos")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_pos' is not a LIST column")?;
        let res_mod_delta = col("res_mod_delta")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_delta' is not a LIST column")?;
        let nterm = col("nterm_delta")?
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or("column 'nterm_delta' is not a DOUBLE column")?;
        let cterm = col("cterm_delta")?
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or("column 'cterm_delta' is not a DOUBLE column")?;
        let mz = col("mz")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'mz' is not a LIST column")?;
        let intensity = col("intensity")?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or("column 'intensity' is not a LIST column")?;

        let list_i32 =
            |list: &ListArray, i: usize| -> Result<Vec<i32>, Box<dyn std::error::Error>> {
                if list.is_null(i) {
                    return Ok(Vec::new());
                }
                let v = list.value(i);
                let a = v
                    .as_any()
                    .downcast_ref::<Int32Array>()
                    .ok_or("list element is not INT32")?;
                Ok((0..a.len()).map(|j| a.value(j)).collect())
            };
        let list_f64 =
            |list: &ListArray, i: usize| -> Result<Vec<f64>, Box<dyn std::error::Error>> {
                if list.is_null(i) {
                    return Ok(Vec::new());
                }
                let v = list.value(i);
                let a = v
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .ok_or("list element is not DOUBLE")?;
                Ok((0..a.len()).map(|j| a.value(j)).collect())
            };
        let list_f32 =
            |list: &ListArray, i: usize| -> Result<Vec<f32>, Box<dyn std::error::Error>> {
                if list.is_null(i) {
                    return Ok(Vec::new());
                }
                let v = list.value(i);
                let a = v
                    .as_any()
                    .downcast_ref::<Float32Array>()
                    .ok_or("list element is not FLOAT")?;
                Ok((0..a.len()).map(|j| a.value(j)).collect())
            };

        for i in 0..batch.num_rows() {
            let seq_s = seq.value(i);
            let ch = charge.value(i);
            if !(1..=255).contains(&ch) {
                return Err(format!("invalid charge {ch} at row {i} of {}", path.display()).into());
            }
            let charge_u8 = ch as u8;

            let positions = list_i32(res_mod_pos, i)?;
            let deltas = list_f64(res_mod_delta, i)?;
            let peptide =
                build_msnet_peptide(seq_s, &positions, &deltas, nterm.value(i), cterm.value(i))?;

            let mzs = list_f32(mz, i)?;
            let ints = list_f32(intensity, i)?;
            if mzs.len() != ints.len() {
                return Err(format!(
                    "mz ({}) and intensity ({}) length mismatch at row {i} of {}",
                    mzs.len(),
                    ints.len(),
                    path.display()
                )
                .into());
            }
            let mut peaks: Vec<(f64, f32)> = mzs
                .iter()
                .zip(ints.iter())
                .map(|(&m, &it)| (m as f64, it))
                .collect();
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            let spectrum = Spectrum {
                title: format!("row{i}"),
                precursor_mz: prec_mz.value(i),
                precursor_intensity: None,
                precursor_charge: Some(charge_u8 as i32),
                rt_seconds: None,
                scan: None,
                peaks,
                activation_method: None,
                isolation_lower_offset: None,
                isolation_upper_offset: None,
                isolation_target_mz: None,
                mobility_window: None,
            };

            out.push(MsnetPsm {
                spectrum,
                peptide,
                charge: charge_u8,
            });
        }
    }
    Ok(out)
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
        cfg.pseudo, cfg.noise_pseudo, cfg.backoff_weight, cfg.min_count
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
            log::info!(
                "prior model = {prior_id} (from {})",
                store_path.display()
            );
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
            .map(|psm| {
                PsmRow {
                    spectrum: &psm.spectrum,
                    peptide: &psm.peptide,
                    charge: psm.charge,
                }
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
            experiment_class: args.protocol.clone(),
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
mod train_from_msnet_tests {
    use super::*;
    use model::mass::H2O;

    #[test]
    fn run_stem_strips_spectrum_extensions() {
        assert_eq!(run_stem(Path::new("/x/a.mzML.gz")), "a");
        assert_eq!(run_stem(Path::new("/x/PXD1_run.d")), "PXD1_run");
        assert_eq!(run_stem(Path::new("/x/PXD1_run.d/")), "PXD1_run");
        assert_eq!(run_stem(Path::new("b.MGF")), "b");
        assert_eq!(run_stem(Path::new("c.raw")), "c");
        assert_eq!(run_stem(Path::new(".d")), ".d");
    }

    fn unmod_mass(seq: &[u8]) -> f64 {
        let rsum: f64 = seq
            .iter()
            .map(|&r| model::AminoAcid::standard(r).unwrap().mass)
            .sum();
        rsum + H2O
    }

    #[test]
    fn unmodified_peptide_mass_matches_reference() {
        let p = build_msnet_peptide("PEPTIDEK", &[], &[], 0.0, 0.0).unwrap();
        let expected = unmod_mass(b"PEPTIDEK");
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues.iter().all(|aa| !aa.is_modified()));
    }

    #[test]
    fn residue_mod_mass_is_correct() {
        const OX: f64 = 15.994915;
        let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], 0.0, 0.0).unwrap();
        let expected = unmod_mass(b"MPEPTIDE") + OX;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(
            p.residues[0].is_modified(),
            "residue 1 should carry the mod"
        );
        assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX);
    }

    #[test]
    fn nterm_mod_mass_is_correct() {
        const ACETYL: f64 = 42.010565;
        let p = build_msnet_peptide("PEPTIDEK", &[], &[], ACETYL, 0.0).unwrap();
        let expected = unmod_mass(b"PEPTIDEK") + ACETYL;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues[0].is_modified());
    }

    #[test]
    fn overlapping_residue_and_nterm_mods_sum() {
        const OX: f64 = 15.994915;
        const ACETYL: f64 = 42.010565;
        let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], ACETYL, 0.0).unwrap();
        let expected = unmod_mass(b"MPEPTIDE") + OX + ACETYL;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX + ACETYL);
    }

    #[test]
    fn mismatched_mod_arrays_error() {
        let r = build_msnet_peptide("PEPTIDE", &[1, 2], &[1.0], 0.0, 0.0);
        assert!(r.is_err());
    }

    #[test]
    fn nonstandard_residue_errors() {
        let r = build_msnet_peptide("PEPTBDE", &[], &[], 0.0, 0.0);
        assert!(r.is_err());
    }
}
