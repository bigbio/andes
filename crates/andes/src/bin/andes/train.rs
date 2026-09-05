//! Search-based model training, MSNet parquet training and model-store updates.

use std::fs::File;
use std::path::{Path, PathBuf};

use crate::cli::{parse_unit_fraction, GeometryArgs, Protocol};
use crate::model_select::{build_aa_set, load_param_from_store, load_seed_param, ModelEntryOwned};
use crate::rescore;
use crate::spectra::spectrum_ext_lower;
use crate::train_intensity::GbdtMode;
use clap::Args;
use input::{MgfReader, MzMLReader};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, ModLocation, Modification,
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
        commit_update, update_add, update_decay, update_remove, update_reweight,
        write_all_models_with_sources_and_gbdt_pub, write_all_models_with_sources_pub,
        SourceLedger,
    },
    ModelStore,
};
use scoring_crate::{Param, RankScorer};
use search::SearchParams;

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
    /// Required for initial training.  In `--update` mode with `--remove-source`
    /// or `--reweight` / `--decay`, `--spectra` is only required when
    /// `--validate` is also given (to run the acceptance gate).
    #[arg(long)]
    pub(crate) spectra: Option<PathBuf>,

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
    /// deglycosylated backbone); Cam-C and any `--mods` are applied via the
    /// AminoAcidSet. Rows whose scan is absent from `--spectra` or whose peptide
    /// fails to parse are skipped. Use to train a glyco-regime rank model from
    /// the reference engine/a glyco search engine backbone IDs (see docs/plans/glyco/50-roadmap/spb-design.md).
    #[arg(long = "labels")]
    pub(crate) labels: Option<PathBuf>,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// When omitted, the bundled `hcd_qexactive_tryp` model is used as the seed.
    #[arg(long = "seed-model")]
    pub(crate) seed_model: Option<String>,

    /// Target-decoy q-value threshold for accepting PSMs as confident training
    /// labels. Use a lenient value (e.g. 0.1 or 0.5) for small fixtures.
    #[arg(long = "train-fdr", default_value = "0.01", value_parser = parse_unit_fraction)]
    pub(crate) train_fdr: f64,

    /// Instrument tag to embed in the trained model's metadata. Default: `QExactive`.
    #[arg(long, default_value = "QExactive")]
    pub(crate) instrument: String,

    /// Experiment-class / protocol tag (e.g. `Automatic`, `TMT`). Default: `Automatic`.
    #[arg(long, default_value = "Automatic")]
    pub(crate) protocol: String,

    /// Path to the Parquet model store to write (created if absent, appended
    /// otherwise). REQUIRED.
    #[arg(long = "out-store")]
    pub(crate) out_store: PathBuf,

    /// Model ID written into the store. Default: `trained_<instrument>_<protocol>`.
    #[arg(long = "model-id")]
    pub(crate) model_id: Option<String>,

    /// Path to a mods.txt file (same format as `--mods` for search). When
    /// omitted, uses built-in defaults (Carbamidomethyl-C fixed, Oxidation-M + protein-N-term-Acetyl
    /// variable).
    #[arg(long)]
    pub(crate) mods: Option<PathBuf>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    pub(crate) threads: usize,

    /// ISO 8601 date string (e.g. `2026-01-01`) recorded in the source ledger.
    /// When omitted, the current date is used for initial training; empty string
    /// is stored when `--date ""` is explicitly passed.
    #[arg(long)]
    pub(crate) date: Option<String>,

    // ── Update mode ──────────────────────────────────────────────────────────
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

/// Available subcommands.
#[derive(clap::Args, Debug)]
pub(crate) struct RescorePinArgs {
    /// Input PIN.
    #[arg(long = "in")]
    pub(crate) input: PathBuf,
    /// Output target PSMs (Percolator `.psms` shape: PSMId, score, q-value, PEP).
    #[arg(long = "out-psms")]
    pub(crate) out_psms: PathBuf,
    /// Output decoy PSMs.
    #[arg(long = "out-dpsms")]
    pub(crate) out_dpsms: PathBuf,
    /// Cross-validation seed.
    #[arg(long = "seed", default_value_t = 42u64)]
    pub(crate) seed: u64,
}

pub(crate) fn run_rescore_pin(args: RescorePinArgs) -> Result<(), Box<dyn std::error::Error>> {
    let pin_text = std::fs::read_to_string(&args.input)
        .map_err(|e| format!("reading {}: {e}", args.input.display()))?;
    let rows = rescore::native_rescore_qvalues(&pin_text, args.seed)?;
    let mut t = std::io::BufWriter::new(std::fs::File::create(&args.out_psms)?);
    let mut d = std::io::BufWriter::new(std::fs::File::create(&args.out_dpsms)?);
    use std::io::Write;
    for w in [&mut t, &mut d] {
        writeln!(
            w,
            "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"
        )?;
    }
    let (mut nt, mut nd) = (0usize, 0usize);
    for (id, is_decoy, q, score) in &rows {
        let w: &mut dyn Write = if *is_decoy { &mut d } else { &mut t };
        writeln!(w, "{id}\t{score}\t{q}\t{q}\t-\t-")?;
        if *is_decoy {
            nd += 1
        } else {
            nt += 1
        }
    }
    eprintln!("rescore-pin: {nt} target and {nd} decoy rows written");
    Ok(())
}

/// Load all MS2 spectra from a path using the same format-dispatch logic as
/// the search path (mzML by extension, otherwise MGF).
pub(crate) fn load_spectra_for_train(
    path: &Path,
) -> Result<Vec<Spectrum>, Box<dyn std::error::Error>> {
    let ext_lower = spectrum_ext_lower(path);
    let mut spectra = Vec::new();
    match ext_lower.as_deref() {
        Some("mzml") => {
            let reader = MzMLReader::new(input::open_buf_maybe_gz(path)?).with_ms_level_range(2, 2);
            for item in reader {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => eprintln!("WARN: mzML parse: {e}"),
                }
            }
        }
        // Native Thermo `.raw` — MS2 only, same Spectrum stream as the search
        // path. Requires building with `--features thermo`.
        Some("raw") => {
            #[cfg(feature = "thermo")]
            {
                let reader = input::ThermoRawReader::open(path)
                    .map_err(|e| format!("open Thermo .raw {}: {e}", path.display()))?
                    .with_ms_level(Some(2));
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => eprintln!("WARN: .raw parse: {e}"),
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
        // Native Bruker timsTOF `.d` (a directory). Requires `--features timstof`.
        Some("d") => {
            #[cfg(feature = "timstof")]
            {
                let reader = input::TimsTofReader::open(path)
                    .map_err(|e| format!("open Bruker .d {}: {e}", path.display()))?;
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => eprintln!("WARN: .d parse: {e}"),
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
        // MGF (default reader).
        _ => {
            let reader = MgfReader::new(input::open_buf_maybe_gz(path)?);
            for item in reader {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => eprintln!("WARN: MGF parse: {e}"),
                }
            }
        }
    }
    Ok(spectra)
}

/// Build the `SearchParams` used by every training mode (initial bootstrap,
/// `--add`, and the acceptance gate) in one place, so they stay consistent with
/// each other and with the production search:
///   - charge span `2..=5` (the search binary's default; `default_tryptic` alone
///     is the narrow `2..=3`, which drops z=4/5 labels common on Astral/timsTOF),
///   - the `NumMods=` variable-mod limit from the `--mods` file when present
///     (the search path applies it too).
pub(crate) fn build_train_search_params(
    mods: &Option<PathBuf>,
) -> Result<SearchParams, Box<dyn std::error::Error>> {
    let aa = build_aa_set(mods)?;
    let mut params = SearchParams::default_tryptic(aa);
    params.charge_range = 2..=5;
    if let Some(path) = mods {
        if let Some(n) = AminoAcidSetBuilder::parse_num_mods_from_file(path)
            .map_err(|e| format!("parsing NumMods= from {}: {e}", path.display()))?
        {
            params.max_variable_mods_per_peptide = n;
        }
    }
    Ok(params)
}

/// First N-X-S/T sequon position (0-based index of the N; X != P), or None.
/// Used to place the glycan mass for training an ETD c/z model: the glycan rides
/// on glycosite-spanning c/z fragments, so it must be baked onto the sequon N.
pub(crate) fn first_nglyco_site(seq: &[u8]) -> Option<usize> {
    (0..seq.len()).find(|&i| {
        seq[i] == b'N'
            && i + 2 < seq.len()
            && seq[i + 1] != b'P'
            && (seq[i + 2] == b'S' || seq[i + 2] == b'T')
    })
}

/// Load external training labels from a `scan\tpeptide\tcharge` TSV (SP-B glyco
/// training). Column order is discovered from the header (case-insensitive), so
/// extra columns are ignored. The peptide is parsed with `aa_set` (Cam-C + any
/// `--mods` applied). Rows with an unknown scan or an unparseable peptide/charge
/// are skipped and counted. `confidence` is set to 0.0 (labels are pre-filtered
/// by the external engine; the value is unused by the accumulator).
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
    // Optional glyco columns: when present, the glycan mass is baked onto the
    // glycosite N so glycosite-spanning c/z fragments carry the intact glycan
    // (required to train a correct ETD c/z model — without it the model learns
    // that glycosite c/z ions are "usually missing", which is false).
    let gly_c = col(&["glycan_mass", "glycan", "glycan_neutral"]);
    let site_c = col(&["glycosite", "site", "glyco_site"]);
    let max_c = scan_c.max(pep_c).max(chg_c);

    // Fixed-mod deltas (e.g. Cam-C 57.02146) from the aa_set: the label peptides
    // are BARE sequences, so we annotate each fixed-mod residue with `+delta` and
    // wrap in `-.SEQ.-` flanking so `Peptide::from_str` mass-matches the fixed
    // variant (a bare `C` would otherwise parse as UNMODIFIED — wrong b/y masses).
    // Variable mods (Ox-M) are left off: the label doesn't say which are modified,
    // and the unmodified form is the correct default for the rank corpus.
    let fixed_deltas: std::collections::HashMap<u8, f64> =
        aa_set.fixed_mod_deltas().into_iter().collect();
    let decorate = |seq: &str| -> String {
        let mut d = String::with_capacity(seq.len() + 4);
        d.push_str("-.");
        for &b in seq.as_bytes() {
            d.push(b as char);
            if let Some(delta) = fixed_deltas.get(&b) {
                d.push_str(&format!("+{:.5}", delta));
            }
        }
        d.push_str(".-");
        d
    };

    let mut labels = Vec::new();
    // ONE label per spectrum: the accumulator must not tally the same spectrum
    // multiple times (it would over-weight that spectrum's rank/edge/charge
    // histograms). External glyco exports can carry rank alternatives or
    // glycoform/site alternatives for one scan, so reject duplicate scans loudly
    // rather than silently double-counting (Codex adversarial-review finding).
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
        // Parse the peptide BEFORE the dedup guard: an unparseable/rank-alternative
        // row must not "claim" the scan and cause a later VALID row for the same
        // scan to be dropped as a duplicate (Codex + code-review finding — external
        // exports can order a weaker alternative first).
        let seq = f[pep_c].trim();
        // Resolve the glycan placement (if this is a glyco corpus): glycan mass
        // from the column, glycosite from an explicit column or the first sequon.
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
        // Bake the glycan directly onto the glycosite residue. We do NOT put it in
        // the peptide string: Peptide::from_str only accepts `+delta` mods that
        // match a REGISTERED variant, and per-PSM glycan masses aren't registered
        // (that dropped ~96% of glyco labels as "unparseable"). The glycan delta on
        // the N flows into every glycosite-spanning c/z fragment mass generically.
        if let Some((site, gmass)) = glyco {
            if let Some(res) = peptide.residues.get_mut(site) {
                let base = res.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
                res.mod_ = Some(std::sync::Arc::new(model::modification::Modification {
                    name: "Glycan".to_string(),
                    mass_delta: base + gmass,
                    residue: model::modification::ResidueSpec::Specific(res.residue),
                    location: model::modification::ModLocation::Anywhere,
                    fixed: false,
                    accession: None,
                    neutral_losses: Vec::new(),
                    loss_class: 1,
                }));
            }
        }
        if !seen_scans.insert(scan) {
            // A second VALID row for a scan already labeled: skip (keep the first).
            dup_scan += 1;
            continue;
        }
        // Charge cross-check: a label charge that disagrees with the spectrum's
        // own precursor charge signals a stale annotation or a scan-mapping error
        // — count it for visibility (the label charge is used for scoring, so a
        // systematic mismatch means the corpus is built on the wrong spectra).
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
    eprintln!(
        "train: loaded {} labels from {} ({} skipped: {} scan-not-found, {} dup-scan, {} unparseable-peptide, {} malformed; {} charge-mismatch vs mzML)",
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
        eprintln!(
            "train: WARNING — {charge_mismatch} labels ({}%) disagree with the mzML precursor \
             charge; the scan->spectrum mapping or the label charges may be wrong.",
            charge_mismatch * 100 / labels.len().max(1),
        );
    }
    Ok(labels)
}

pub(crate) fn run_train_from_search(
    args: TrainFromSearchArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
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

    // ── Standard training path ────────────────────────────────────────────────

    // ── 2. Load spectra ───────────────────────────────────────────────────────
    let spectra_path = args
        .spectra
        .clone()
        .ok_or("--spectra is required for initial training")?;
    eprintln!("train: loading spectra from {} ...", spectra_path.display());
    let spectra = load_spectra_for_train(&spectra_path)?;
    eprintln!("train: loaded {} spectra", spectra.len());

    // ── 3. Load seed Param + RankScorer ──────────────────────────────────────
    let (seed_model_id, seed_param): (String, Param) = load_seed_param(&args.seed_model)?;
    eprintln!("train: seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4-5. Labels: external TSV (SP-B) OR a seed search (bootstrap) ──────────
    let labels = if let Some(ref labels_path) = args.labels {
        // SP-B glyco path: labels come from an external engine's backbone IDs.
        // No seed search / database needed — the peptides are given directly.
        let aa_set = build_aa_set(&args.mods)?;
        load_labels_from_tsv(labels_path, &spectra, &aa_set)?
    } else {
        let database = args
            .database
            .clone()
            .ok_or("--database is required for initial training (or pass --labels)")?;
        let search_params = build_train_search_params(&args.mods)?;
        eprintln!(
            "train: running seed search (train-fdr = {}) ...",
            args.train_fdr
        );
        bootstrap_labels(
            &spectra,
            &database,
            &seed_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("bootstrap_labels: {e}"))?
    };
    eprintln!("train: {} confident labels", labels.len());

    if labels.is_empty() {
        return Err(format!(
            "no confident labels found at train-fdr={} — try a higher --train-fdr",
            args.train_fdr
        )
        .into());
    }

    // ── 5b. Geometry template: own-derived (DEFAULT) or seed (opt-out) ─────────
    // andes derives the partition/segment geometry from THIS corpus by default
    // (own geometry — no MS-GF+ partition structure); the seed supplies only
    // non-geometry metadata. Opt out with --seed-geometry to reuse the
    // seed's geometry (e.g. to reproduce a legacy model). Own geometry is
    // entrapment-FDP-validated to beat seed geometry on honest PSMs AND speed
    // across Astral (+57%), UPS1 (+15%) and TMT (+50%).
    let use_seed_geometry = args.seed_geometry;
    let template: Param = if !use_seed_geometry {
        eprintln!(
            "train: deriving own partition geometry from {} PSMs (pass --seed-geometry to reuse the seed geometry)",
            labels.len()
        );
        let corpus = corpus_charge_masses(&labels);
        let geo_cfg = GeometryArgs::default().to_config();
        derive_geometry(&corpus, &seed_param, &geo_cfg)
    } else {
        eprintln!("train: --seed-geometry set — reusing seed partition geometry");
        seed_param.clone()
    };
    let accum_scorer = RankScorer::new(&template);

    // ── 6. Accumulate stats ───────────────────────────────────────────────────
    eprintln!("train: accumulating ion-match statistics ...");
    let accumulator = StatsAccumulator::new(&accum_scorer);
    let mut stats = CountStats::new();
    for label in &labels {
        let spec = &spectra[label.spectrum_index];
        accumulator.accumulate(&mut stats, spec, &label.peptide, label.charge);
    }
    let stats = merge(vec![stats]);

    // ── 7. Estimate model ─────────────────────────────────────────────────────
    eprintln!("train: estimating model parameters ...");
    let cfg = EstimatorConfig::default();
    let estimator = Estimator::new(cfg);
    let trained_param = estimator.estimate(&stats, &template);
    let n_partitions = trained_param.partitions.len();
    eprintln!("train: trained model has {} partitions", n_partitions);

    // ── 8. Determine model ID ─────────────────────────────────────────────────
    let model_id = args
        .model_id
        .clone()
        .unwrap_or_else(|| format!("trained_{}_{}", args.instrument, args.protocol));
    eprintln!("train: model ID = {model_id}");

    // ── 9. Build source ledger (Part A) ───────────────────────────────────────
    // Determine date: use args.date if provided, else today's date in ISO 8601.
    let date_str = args.date.clone().unwrap_or_else(|| {
        // Format today as YYYY-MM-DD using std::time.
        format_today_iso8601()
    });
    let spectra_filename = spectra_path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| spectra_path.display().to_string());
    let source_id = format!("bootstrap_{model_id}");
    let ledger = SourceLedger {
        source_id: source_id.clone(),
        dataset: spectra_filename,
        n_psms: labels.len() as i64,
        date: date_str,
        weight: 1.0,
        train_fdr: args.train_fdr as f32,
        instrument: args.instrument.clone(),
        experiment_class: args.protocol.clone(),
    };

    // ── 10. Write to store with source tracking ────────────────────────────────
    // Read existing OTHER models from the store (preserve them on append).
    let store_path = &args.out_store;
    let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
    if store_path.exists() {
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
        for id in store.model_ids() {
            if id == model_id {
                eprintln!("train: overwriting existing model '{id}' in store");
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

    // Combine all models: write the trained model + existing others together.
    let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
    all_entries.push((
        model_id.clone(),
        trained_param.clone(),
        vec![(ledger, stats)],
    ));
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

    eprintln!(
        "train: wrote model '{model_id}' to {} (source '{source_id}') [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );
    eprintln!(
        "train: summary — labels={}, partitions={}, store={}",
        labels.len(),
        n_partitions,
        store_path.display(),
    );

    Ok(())
}

// ════════════════════════════════════════════════════════════════════════════
// train-from-msnet: supervised training from externally-labeled PSM parquets
// ════════════════════════════════════════════════════════════════════════════

/// One confident, externally-labeled PSM read from a flat training parquet.
///
/// This is the in-memory form of one parquet row (see the column contract in
/// the `train-from-msnet` CLI docs). The peaks are stored already sorted
/// ascending by m/z (the parquet stores acquisition order; the reader sorts).
pub(crate) struct MsnetPsm {
    pub(crate) spectrum: Spectrum,
    pub(crate) peptide: model::peptide::Peptide,
    pub(crate) charge: u8,
}

/// Build a [`model::peptide::Peptide`] from a bare uppercase sequence plus
/// resolved modification mass-deltas.
///
/// Modifications are applied two ways, matching how the scoring code computes
/// peptide mass (`Peptide::new` sums `aa.mass + mod_.mass_delta` per residue):
///
/// - **Residue mods** (`res_mod_pos` / `res_mod_delta`, 1-based positions):
///   attached to the residue at that position via `AminoAcid::with_mod` with a
///   `ModLocation::Anywhere`, `ResidueSpec::Specific(residue)` modification.
/// - **Terminal mods** (`nterm_delta` / `cterm_delta`): folded onto the first /
///   last residue's mass-delta. Because `Peptide::new` only sums per-residue
///   `mod_.mass_delta`, a terminal delta must be carried by a residue's `mod_`
///   to be counted. If a residue already carries a residue mod, the terminal
///   delta is *added* to that residue's existing delta (a single combined
///   `Modification`); otherwise a fresh terminal `Modification` is attached.
///   This keeps `peptide.mass()` correct regardless of overlap.
///
/// Returns an error if `seq` contains a non-standard residue or the mod arrays
/// are misaligned.
pub(crate) fn build_msnet_peptide(
    seq: &str,
    res_mod_pos: &[i32],
    res_mod_delta: &[f64],
    nterm_delta: f64,
    cterm_delta: f64,
) -> Result<model::peptide::Peptide, Box<dyn std::error::Error>> {
    use std::sync::Arc;

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

    // Accumulate the total mod delta to apply to each residue (1-based -> 0-based).
    // Residue mods first, then terminal deltas folded onto the end residues.
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
            };
            residues.push(aa.with_mod(Arc::new(m)));
        } else {
            residues.push(aa);
        }
    }

    // Flanking residues per the spec: pre=`_`, post=`-`.
    Ok(model::peptide::Peptide::new(residues, b'_', b'-'))
}

/// Read one flat training parquet into a vector of [`MsnetPsm`].
///
/// Reads via the workspace `parquet`/`arrow` crates in record-batch chunks.
/// List columns (`res_mod_pos`, `res_mod_delta`, `mz`, `intensity`) are
/// decoded per-row from their `ListArray` offsets. Peaks are sorted ascending
/// by m/z (the parquet stores acquisition order).
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

        // Helper to pull a Vec<i32> out of one ListArray row.
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
            // Input is acquisition order; the scoring path requires ascending m/z.
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

// train-intensity: merge partial intensity stats into a finalized model parquet
// ─────────────────────────────────────────────────────────────────────────────

/// `andes train`: train a scoring model directly from externally-labeled PSM
/// parquets, reusing the existing accumulate → estimate → store machinery but
/// bypassing the bootstrap search.
pub(crate) fn run_train(args: TrainArgs) -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;

    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    // ── 2. Read all input parquets ────────────────────────────────────────────
    let mut psms: Vec<MsnetPsm> = Vec::new();
    let mut rows_read = 0usize;
    for input in &args.inputs {
        eprintln!("train: reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        rows_read += part.len();
        eprintln!("train:   {} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    eprintln!(
        "train: {rows_read} total PSM rows across {} file(s)",
        args.inputs.len()
    );

    // ── 3. Load seed Param and apply the fragment-tolerance override ──────────
    let (seed_model_id, mut seed_param): (String, Param) =
        load_seed_param(&Some(args.seed_model.clone()))?;
    eprintln!("train: seed model = {seed_model_id}");
    if let Some(ppm) = args.fragment_tol_ppm {
        seed_param.mme = Tolerance::Ppm(ppm);
        eprintln!("train: fragment tolerance overridden to {ppm} ppm");
    } else if let Some(da) = args.fragment_tol_da {
        seed_param.mme = Tolerance::Da(da);
        eprintln!("train: fragment tolerance overridden to {da} Da");
    } else {
        eprintln!("train: using seed fragment tolerance {:?}", seed_param.mme);
    }

    // ── 3b. Geometry template: own-derived (DEFAULT) or seed (opt-out) ─────────
    // andes derives the partition/segment geometry from THIS corpus by default
    // (own geometry — no MS-GF+ partition structure); the seed supplies only
    // non-geometry metadata. Opt out with --seed-geometry to reuse the
    // seed's geometry (e.g. to reproduce a legacy model). Own geometry is
    // entrapment-FDP-validated to beat seed geometry on honest PSMs AND speed
    // across Astral (+57%), UPS1 (+15%) and TMT (+50%).
    let use_seed_geometry = args.seed_geometry;
    let template: Param = if !use_seed_geometry {
        eprintln!(
            "train: deriving own partition geometry from {} PSMs (pass --seed-geometry to reuse the seed geometry)",
            psms.len()
        );
        let corpus: Vec<(i32, f32)> = psms
            .iter()
            .map(|p| (p.charge as i32, p.peptide.mass() as f32))
            .collect();
        let geo_cfg = GeometryArgs::default().to_config();
        derive_geometry(&corpus, &seed_param, &geo_cfg)
    } else {
        eprintln!("train: --seed-geometry set — reusing seed partition geometry");
        seed_param.clone()
    };

    // Build the scorer AFTER the tolerance override so accumulation uses it.
    let seed_scorer = RankScorer::new(&template);

    // ── 4. Accumulate ion-match statistics (parallel; per-worker CountStats) ──
    eprintln!("train: accumulating ion-match statistics ...");
    let stats = psms
        .par_iter()
        .fold(CountStats::new, |mut acc, psm| {
            let accumulator = StatsAccumulator::new(&seed_scorer);
            accumulator.accumulate(&mut acc, &psm.spectrum, &psm.peptide, psm.charge);
            acc
        })
        .collect::<Vec<_>>();
    let stats = merge(stats);
    eprintln!("train: accumulated {} PSMs", psms.len());

    // ── 5. Estimate the model (replaces all learned tables in the seed) ───────
    eprintln!("train: estimating model parameters ...");
    let cfg = EstimatorConfig {
        pseudo: args.train_pseudo,
        noise_pseudo: args.train_noise_pseudo,
        min_count: args.train_min_count,
        backoff_weight: args.train_backoff_weight,
        error_scaling_factor_override: None,
        rank_smoothing: args.rank_smoothing,
    };
    eprintln!(
        "train: estimator pseudo={} noise_pseudo={} backoff_weight={} min_count={}",
        cfg.pseudo, cfg.noise_pseudo, cfg.backoff_weight, cfg.min_count
    );
    let estimator = Estimator::new(cfg);

    // Optional independent prior: sparse partitions shrink toward this model
    // (Level 0 of the backoff hierarchy) instead of the corpus-internal pool.
    // When `--prior-model` is omitted, default to the trained model id. The
    // selection columns passed to `load_param_from_store` are inert here because
    // `model_id_override` is `Some` (it loads that exact id).
    let prior_param: Option<Param> = match &args.prior_model_store {
        Some(store_path) => {
            let prior_id = args
                .prior_model
                .clone()
                .unwrap_or_else(|| args.model_id.clone());
            // `load_param_from_store`'s activation/instrument/protocol are only
            // consulted for automatic selection; passing an explicit
            // `model_id_override` makes them inert, so `Protocol::Auto` (the CLI
            // enum the signature expects) is a harmless placeholder.
            let (_pid, p) = load_param_from_store(
                seed_param.data_type.activation,
                Some(seed_param.data_type.instrument),
                Protocol::Auto,
                model::enzyme::Enzyme::Trypsin, // inert (model_id_override below makes selection columns unused)
                Some(store_path.as_path()),
                Some(&prior_id),
            )
            .map_err(|e| format!("loading --prior-model '{prior_id}': {e}"))?;
            eprintln!(
                "train: prior model = {prior_id} (from {})",
                store_path.display()
            );
            Some(p)
        }
        None => None,
    };

    let mut trained_param = estimator.estimate_with_prior(&stats, &template, prior_param.as_ref());
    let n_partitions = trained_param.partitions.len();
    eprintln!("train: trained model has {n_partitions} partitions");

    // ── 5b. Override the selection-relevant data_type from flags ──────────────
    // The trained model inherits the seed's data_type; minting a NEW slug whose
    // (activation, instrument, enzyme, protocol) differs from the seed requires
    // overriding those columns explicitly, otherwise model selection (which keys
    // on these columns, not the model_id string) would never route to it.
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
    eprintln!(
        "train: model data_type = {:?}/{:?}/{:?}/{:?}",
        trained_param.data_type.activation,
        trained_param.data_type.instrument,
        trained_param.data_type.enzyme,
        trained_param.data_type.protocol,
    );

    // ── 6. Build the source ledger (sentinel train_fdr; pre-labeled input) ────
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
        train_fdr: 1.0, // sentinel: input is pre-labeled, no q-value filtering here
        instrument: trained_param.data_type.instrument.name().to_string(),
        experiment_class: trained_param.data_type.protocol.name().to_string(),
    };

    // ── 6b. Optionally train a GBDT peak model ────────────────────────────────
    let gbdt_blob: Option<Vec<u8>> = if args.gbdt != GbdtMode::Off {
        use model_train::gbdt::dataset::{build_dataset, PsmRow};
        use model_train::gbdt::train::{train_gbdt, TrainParams};

        eprintln!("train: building GBDT dataset from {} PSMs ...", psms.len());
        let gbdt_rows: Vec<PsmRow<'_>> = psms
            .iter()
            .map(|psm| {
                // Pass the full mod-carrying peptide: labels are mod-aware so b/y
                // ions over Cam-C/TMT/Ox-M land at the correct (shifted) m/z.
                PsmRow {
                    spectrum: &psm.spectrum,
                    peptide: &psm.peptide,
                    charge: psm.charge,
                }
            })
            .collect();

        let gbdt_dataset = build_dataset(&gbdt_rows, &seed_scorer);
        eprintln!(
            "train: GBDT dataset: {} peak rows, {} positives",
            gbdt_dataset.y.len(),
            gbdt_dataset.y.iter().filter(|&&l| l == 1).count(),
        );

        let gbdt_params = TrainParams {
            allow_degenerate: args.allow_degenerate_model,
            ..TrainParams::default()
        };
        // Hard-error on quality-gate failure (finding 3.6) instead of writing a
        // degenerate model into the store (unless opted out).
        let trained_gbdt = train_gbdt(&gbdt_dataset, &gbdt_params, 42)?;
        eprintln!("train: GBDT trained: {} trees", trained_gbdt.trees.len(),);
        Some(trained_gbdt.to_bytes())
    } else {
        None
    };

    // ── 7. Write to store, preserving any other existing models ───────────────
    // Unified path: ALWAYS reads existing models from the store, builds the
    // combined entries list (new model + existing others), and writes once via
    // write_all_models_with_sources_and_gbdt_pub.  The GBDT blob (Some or None)
    // flows through for the new model; existing models' blobs are re-serialised
    // from their loaded gbdt_peak_model so no blob is ever lost on re-write.
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
                    eprintln!("train: overwriting existing model '{id}' in store");
                    continue;
                }
                let p = store
                    .load_param(&id)
                    .map_err(|e| format!("reading model '{id}': {e}"))?;
                // Preserve the GBDT blob (if any) for the existing model.
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

        // Build the combined entries list: new model first, then existing others.
        let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
        all_entries.push((model_id.clone(), trained_param, vec![(ledger, stats)]));
        for (id, p, src) in existing_other {
            all_entries.push((id, p, src));
        }

        // Parallel blobs: new model gets gbdt_blob (Some or None); existing models
        // get their preserved blobs.
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

    eprintln!(
        "train: wrote model '{model_id}' to {} (source '{}', {} PSMs, {n_partitions} partitions) [{:.2}s]",
        store_path.display(),
        args.source,
        psms.len(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

/// Incremental update mode (Part D): `--update <MODEL_ID>` plus one of
/// `--add`, `--remove-source`, `--reweight`, `--decay`.
pub(crate) fn run_train_update(
    args: TrainFromSearchArgs,
    model_id: &str,
    t0: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let store_path = &args.out_store;
    let cfg = EstimatorConfig::default();

    // ── Dispatch to the right update operation ────────────────────────────────
    let (candidate, new_sources) = if args.add {
        // --add mode: search spectra, accumulate stats, call update_add.
        let spectra_path = args
            .spectra
            .clone()
            .ok_or("--spectra is required with --add")?;
        let database = args
            .database
            .clone()
            .ok_or("--database is required with --add")?;
        let source_id = args
            .source_id
            .clone()
            .ok_or("--source-id is required with --add")?;

        eprintln!(
            "train update: loading spectra from {} ...",
            spectra_path.display()
        );
        let spectra = load_spectra_for_train(&spectra_path)?;
        eprintln!("train update: loaded {} spectra", spectra.len());

        // Load the current stored model as the seed.
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening store {}: {e}", store_path.display()))?;
        let current_param = store
            .load_param(model_id)
            .map_err(|e| format!("loading model '{model_id}': {e}"))?;
        let current_scorer = RankScorer::new(&current_param);

        let search_params = build_train_search_params(&args.mods)?;

        eprintln!(
            "train update: running seed search (train-fdr={}) ...",
            args.train_fdr
        );
        let labels = bootstrap_labels(
            &spectra,
            &database,
            &current_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("bootstrap_labels: {e}"))?;
        eprintln!("train update: {} confident labels", labels.len());

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
            instrument: args.instrument.clone(),
            experiment_class: args.protocol.clone(),
        };

        update_add(store_path, model_id, ledger, stats, cfg)
            .map_err(|e| format!("update_add: {e}"))?
    } else if let Some(ref sid) = args.remove_source.clone() {
        update_remove(store_path, model_id, sid, cfg).map_err(|e| format!("update_remove: {e}"))?
    } else if let Some(ref spec) = args.reweight.clone() {
        // Parse "source-id=weight"
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

    // ── Acceptance gate (Part D) ──────────────────────────────────────────────
    let commit = if let Some(ref validate_path) = args.validate.clone() {
        let database = args
            .database
            .clone()
            .ok_or("--database is required with --validate")?;

        eprintln!(
            "train update: running acceptance gate on {} ...",
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

        let search_params = build_train_search_params(&args.mods)?;

        let delta = evaluate_candidate(
            &val_spectra,
            &database,
            &current_scorer,
            &candidate_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("evaluate_candidate: {e}"))?;

        eprintln!(
            "train update: gate — current={} PSMs, candidate={} PSMs at FDR={}",
            delta.current_count, delta.candidate_count, args.train_fdr
        );

        if delta.is_accepted() {
            eprintln!("train update: ACCEPTED (candidate >= current)");
            true
        } else {
            eprintln!("train update: REJECTED (candidate < current)");
            if args.force {
                eprintln!("train update: --force set, committing anyway");
                true
            } else {
                eprintln!("train update: skipping commit (use --force to override)");
                false
            }
        }
    } else {
        eprintln!("train update: no --validate dataset; skipping acceptance gate");
        if args.force {
            eprintln!("train update: --force set, committing unconditionally");
        }
        // Without --validate, commit unless user explicitly uses --force to control.
        // Default: commit (no gate run = no evidence of regression).
        true
    };

    if commit {
        commit_update(store_path, model_id, &candidate, &new_sources)
            .map_err(|e| format!("commit_update: {e}"))?;
        eprintln!(
            "train update: committed model '{model_id}' to {} [{:.2}s]",
            store_path.display(),
            t0.elapsed().as_secs_f64(),
        );
    }

    Ok(())
}

/// Parse `"source-id=weight"` from a `--reweight` argument.
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

/// Format today's date as `YYYY-MM-DD` using `std::time::SystemTime`.
pub(crate) fn format_today_iso8601() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Simple Gregorian calendar conversion from Unix timestamp (days since epoch).
    let days = secs / 86400;
    unix_days_to_iso8601(days)
}

pub(crate) fn unix_days_to_iso8601(days: u64) -> String {
    // Algorithm: Gregorian calendar from Julian Day Number.
    // JDN for 1970-01-01 = 2440588.
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

    /// Reference: unmodified peptide mass = sum of residue masses + H2O.
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

    /// Oxidation on residue 1 of "MPEPTIDE" must add exactly +15.994915 Da to
    /// the unmodified mass — verifying mods are actually applied, not dropped.
    #[test]
    fn residue_mod_mass_is_correct() {
        const OX: f64 = 15.994915;
        let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], 0.0, 0.0).unwrap();
        let expected = unmod_mass(b"MPEPTIDE") + OX;
        // Exact: Peptide::new sums aa.mass + mod_.mass_delta, so the delta is
        // added with the same arithmetic as the reference.
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(
            p.residues[0].is_modified(),
            "residue 1 should carry the mod"
        );
        assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX);
    }

    /// N-terminal Acetyl folds onto the first residue's mass-delta.
    #[test]
    fn nterm_mod_mass_is_correct() {
        const ACETYL: f64 = 42.010565;
        let p = build_msnet_peptide("PEPTIDEK", &[], &[], ACETYL, 0.0).unwrap();
        let expected = unmod_mass(b"PEPTIDEK") + ACETYL;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues[0].is_modified());
    }

    /// A residue mod and a terminal mod on the SAME residue must sum (not
    /// clobber) so the total mass stays correct.
    #[test]
    fn overlapping_residue_and_nterm_mods_sum() {
        const OX: f64 = 15.994915;
        const ACETYL: f64 = 42.010565;
        // Oxidation on residue 1 (M) AND N-term acetyl on the same first residue.
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
