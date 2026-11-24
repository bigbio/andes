//! Arguments of the subcommands: rescore, database, inference, protein, cohort and quant.

use std::path::PathBuf;

use clap::Args;
use model::Tolerance;
use search::decoy::DecoyStrategy;

use super::flags::{
    parse_decoy_strategy, parse_enzyme_list, parse_transfer_models, EnzymeSpecificity, Listed,
};
use super::options::{proteoform_fdr, DatabaseArgs};
use super::{
    ADVANCED_FDR, DECOY_STRATEGIES, DEFAULT_DECOY_STRATEGY, ENZYMES, H_ADVANCED, H_FDR, H_IO,
    H_ISOBARIC,
};
use crate::cli::parse::{parse_precursor_tol, parse_unit_fraction};

#[derive(Args, Debug)]
pub(crate) struct RescoreArgs {
    /// Fixed storage budget of the engine in GB, not a process-memory cap: PIN features, PSM texts, normalised training rows and score buffers spill to --temp-dir once they would take more than their share of it; per-PSM metadata (about 160 bytes a PSM) is charged to it but always stays in memory [default: what the PIN needs to stay in memory, at most 16 GB and at most the memory free to the job or machine less a margin]
    #[arg(long, value_name = "GB",
          value_parser = clap::value_parser!(u64).range(1..))]
    pub(crate) rescore_memory_gb: Option<u64>,

    /// Directory of the spill files, created if missing
    #[arg(long, value_name = "DIR", required = true)]
    pub(crate) temp_dir: PathBuf,

    /// Percolator arguments as given to the percolator binary, the PIN last
    #[arg(value_name = "PERCOLATOR_ARGS", required = true, num_args = 1..,
          trailing_var_arg = true, allow_hyphen_values = true)]
    pub(crate) args: Vec<String>,
}

/// `andes database`: the target + decoy database of a study or cohort.
#[derive(Args, Debug)]
pub(crate) struct DatabaseCmdArgs {
    /// Target protein FASTA, may be gzipped
    #[arg(long, value_name = "FILE")]
    pub(crate) fasta: PathBuf,

    /// Target + decoy FASTA to write
    ///
    /// The settings table <stem>.database.parquet and the overlap report
    /// <stem>.overlap.parquet are written beside it. A search given this FASTA as --database
    /// reads the settings table and searches these decoys.
    #[arg(long, value_name = "FILE")]
    pub(crate) output: PathBuf,

    /// Protease of the searches, as the search's --enzyme
    #[arg(long, value_name = "ENZYME", default_value = "trypsin",
          value_parser = Listed::new(ENZYMES, parse_enzyme_list))]
    pub(crate) enzyme: String,

    /// Peptide ends that must be cleavage sites, as the search's --enzyme-specificity
    #[arg(
        long,
        value_name = "SPECIFICITY",
        value_enum,
        ignore_case = true,
        default_value = "fully"
    )]
    pub(crate) enzyme_specificity: EnzymeSpecificity,

    /// Shortest peptide, as the search's --min-length
    #[arg(long, value_name = "N", default_value = "6")]
    pub(crate) min_length: u32,

    /// Longest peptide, as the search's --max-length
    #[arg(long, value_name = "N", default_value = "50")]
    pub(crate) max_length: u32,

    /// Missed cleavages, as the search's --max-missed-cleavages (give 3 for --glyco searches)
    #[arg(long, value_name = "N", default_value = "2")]
    pub(crate) max_missed_cleavages: u32,

    /// Another digest the cohort searches these decoys with, as ENZYMES/SPECIFICITY/MIN-MAX/MISSED; repeatable
    ///
    /// For example lysc/fully/7-30/2 or trypsin,gluc/semi/6-50/2. The decoys are built for
    /// the digest of the options above and every decoy peptide that any listed digest
    /// enumerates is checked against that digest's target peptides; a search may use the
    /// database with any one of the digests.
    #[arg(long, value_name = "DIGEST", value_parser = parse_check_digest)]
    pub(crate) check_digest: Vec<CheckDigest>,

    /// How decoys are generated, as the search's --decoy-strategy
    #[arg(long, value_name = "STRATEGY", default_value = DEFAULT_DECOY_STRATEGY,
          value_parser = Listed::new(DECOY_STRATEGIES, parse_decoy_strategy))]
    pub(crate) decoy_strategy: DecoyStrategy,

    /// Accession prefix of the decoys, as the search's --decoy-prefix
    #[arg(
        long,
        value_name = "PREFIX",
        default_value = "XXX_",
        hide_short_help = true
    )]
    pub(crate) decoy_prefix: String,

    /// Accession suffix that also marks a decoy, as the search's --decoy-suffix
    #[arg(long, value_name = "SUFFIX", hide_short_help = true)]
    pub(crate) decoy_suffix: Option<String>,

    /// Seed of the decoy shuffles, as the search's --decoy-seed.
    #[arg(long, hide = true, default_value_t = search::decoy::DEFAULT_DECOY_SEED)]
    pub(crate) decoy_seed: u64,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N", default_value_t = num_cpus::get(), hide_default_value = true)]
    pub(crate) threads: usize,
}

/// One `--check-digest` of `andes database`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CheckDigest {
    pub(crate) enzyme: String,
    pub(crate) specificity: EnzymeSpecificity,
    pub(crate) min_length: u32,
    pub(crate) max_length: u32,
    pub(crate) max_missed_cleavages: u32,
}

pub(crate) fn parse_check_digest(s: &str) -> Result<CheckDigest, String> {
    let usage =
        || format!("'{s}': expected ENZYMES/SPECIFICITY/MIN-MAX/MISSED, e.g. lysc/fully/7-30/2");
    let parts: Vec<&str> = s.split('/').map(str::trim).collect();
    let [enzyme, specificity, lengths, missed] = parts.as_slice() else {
        return Err(usage());
    };
    let enzyme = parse_enzyme_list(enzyme)?;
    let specificity = <EnzymeSpecificity as clap::ValueEnum>::from_str(specificity, true)
        .map_err(|_| format!("{}: specificity is fully, semi or non-specific", usage()))?;
    let (min, max) = lengths
        .split_once('-')
        .and_then(|(a, b)| Some((a.trim().parse::<u32>().ok()?, b.trim().parse::<u32>().ok()?)))
        .filter(|(a, b)| a <= b)
        .ok_or_else(|| format!("{}: lengths are MIN-MAX with MIN <= MAX", usage()))?;
    let missed = missed
        .parse::<u32>()
        .map_err(|_| format!("{}: missed cleavages are a number", usage()))?;
    Ok(CheckDigest {
        enzyme,
        specificity,
        min_length: min,
        max_length: max,
        max_missed_cleavages: missed,
    })
}

/// `andes inference`: protein groups over the rescored PSMs of one or more runs.
#[derive(Args, Debug)]
pub(crate) struct InferenceArgs {
    /// Percolator target PSM results (--results-psms) of every run
    ///
    /// Their scores must be on one scale: one Percolator run over all the runs (one file,
    /// or one file per run from that run).
    #[arg(long = "psms", value_name = "TSV", num_args = 1.., required = true)]
    pub(crate) psms: Vec<PathBuf>,

    /// Percolator decoy PSM results (--decoy-results-psms) of the same Percolator run
    #[arg(long = "decoy-psms", value_name = "TSV", num_args = 1.., required = true)]
    pub(crate) decoy_psms: Vec<PathBuf>,

    /// The protein database of the search, with its options below
    ///
    /// With the search's --enzyme, --fixed-mods, --variable-mods, decoy affixes,
    /// --annotations, --uniprot-* and --knockoffs, protein families, knockoff entries, knockoff sites and decoy pairs are
    /// the search's.
    #[arg(long, value_name = "FILE")]
    pub(crate) database: PathBuf,

    #[command(flatten)]
    pub(crate) proteins: ProteinArgs,

    /// Run names of the groups' grouped_runs [default: the runs in the --psms files' FileName
    /// column; one per file, named after it, for a file without that column]
    #[arg(long = "runs", value_name = "NAME", num_args = 1..)]
    pub(crate) runs: Vec<String>,

    /// Output prefix: writes <PREFIX>.pg.parquet (the QPX protein-group view)
    #[arg(short = 'o', long = "output", value_name = "PREFIX")]
    pub(crate) output: PathBuf,
}

/// How protein inference reads the database and which PSMs count as evidence:
/// the search's database options and the PSM thresholds.
#[derive(Args, Debug, Clone)]
pub(crate) struct ProteinArgs {
    #[command(flatten)]
    pub(crate) db: DatabaseArgs,

    /// The search's --config file, for its database and modification keys
    ///
    /// Its enzyme, modification lists (with their config-only settings), decoy affixes,
    /// UniProt options, annotations and knockoffs apply here, below the flags given on the
    /// command line. Its other keys are not read.
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) config: Option<PathBuf>,

    /// PSM q-value threshold of the PSMs and peptides counted per group
    #[arg(long = "psm-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) psm_fdr: f64,

    /// PSM q-value threshold of proteoform claims (variant-class PSMs) [default: the larger of --psm-fdr and 0.05]
    #[arg(long = "proteoform-fdr", value_name = "FDR", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = ADVANCED_FDR)]
    pub(crate) proteoform_fdr: Option<f64>,

    /// Protein q-value threshold of the reported groups (1 keeps every group)
    ///
    /// The pg view keeps the target and decoy groups at or below it (the decoys show the
    /// FDR at the threshold), and `quant finalize`'s feature view the features of those
    /// target groups.
    #[arg(long = "protein-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) protein_fdr: f64,
}

/// `andes cohort`: plan, search and run the stages of several projects' SDRFs together.
#[derive(clap::Subcommand, Debug)]
pub(crate) enum CohortCommand {
    /// Check a cohort and write its plan: every run's search settings and the FDR groups
    ///
    /// Searches nothing. Prints the groups and every problem found, and writes the plan
    /// only when there is no error.
    Plan(CohortPlanArgs),
    /// Crosslink FDR over the CSM tables of several runs, as one search would compute it
    ///
    /// Reads each run's <base>.xl.tsv, rescores the intra-protein links with one Percolator
    /// model over all runs' <base>.xl.pin (the native engine; homotypic and inter-links keep
    /// the pair score), groups all CSMs into residue pairs and protein pairs, and gives every level
    /// its q-values once over all runs. Writes <output>.xl.tsv,
    /// <output>.xl.residue_pairs.tsv, <output>.xl.ppi.tsv and each project's share of the
    /// pooled list, <output>.xl.projects.parquet.
    #[command(name = "xl-fdr", hide = true)]
    XlFdr(CohortXlFdrArgs),
    /// Plan the cohort, search every run, then run each group's stages
    ///
    /// Each run's outputs go to <output_dir>/runs/<run>/, each group's to
    /// <output_dir>/groups/<group>/. A finished search is kept while its spectrum file, its
    /// settings, the files they name and its outputs are unchanged, so a rerun searches only
    /// new, failed or changed runs.
    Run(CohortRunArgs),
}

#[derive(clap::Args, Debug)]
pub(crate) struct CohortRunArgs {
    /// Config file with a `cohort` section
    #[arg(value_name = "CONFIG")]
    pub(crate) config: PathBuf,
    /// Worker threads of each search [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

#[derive(clap::Args, Debug)]
pub(crate) struct CohortXlFdrArgs {
    /// The runs' CSM tables (<base>.xl.tsv from `andes --xl`)
    #[arg(value_name = "CSM_TABLE", required = true)]
    pub(crate) tables: Vec<PathBuf>,
    /// Output base name: <output>.xl.tsv, <output>.xl.residue_pairs.tsv, <output>.xl.ppi.tsv
    #[arg(long, short, value_name = "BASE")]
    pub(crate) output: PathBuf,
    /// The cohort's plan.json, which gives each run's project for the per-project error
    /// table [default: the run name up to its first `.`]
    #[arg(long, value_name = "PLAN_JSON")]
    pub(crate) plan: Option<PathBuf>,
    /// The database the runs searched; with it, inter-links whose peptide has a near-copy in
    /// another protein are estimated apart (inter_ambiguous) and inter-links are rescored
    #[arg(long, value_name = "FASTA")]
    pub(crate) database: Option<PathBuf>,
}

#[derive(clap::Args, Debug)]
pub(crate) struct CohortPlanArgs {
    /// Config file with a `cohort` section; its other sections are the settings every
    /// run shares
    #[arg(value_name = "CONFIG")]
    pub(crate) config: PathBuf,
    /// Where to write the plan [default: <output_dir>/plan.json]
    #[arg(long, short, value_name = "FILE")]
    pub(crate) output: Option<PathBuf>,
}

/// The `default_value` of the flag `id` of the clap arguments `A`, so that code which
/// needs a flag's default reads it from the one place it is declared.
pub(crate) fn default_value<A: clap::Args>(id: &str) -> String {
    let cmd = A::augment_args(clap::Command::new("defaults"));
    let value = cmd
        .get_arguments()
        .find(|a| a.get_id() == id)
        .and_then(|a| a.get_default_values().first())
        .map(|v| v.to_string_lossy().into_owned());
    value.unwrap_or_else(|| panic!("flag {id} has no default value"))
}

/// `andes quant`: the match-between-runs stages over a quant directory
/// written by a search with `--quant-dir`.
#[derive(clap::Subcommand, Debug)]
pub(crate) enum QuantCommand {
    /// Build the precursor library and the joint retention-time alignment
    ///
    /// From every run's PSM table and the Percolator results of those PSMs. Writes
    /// library.parquet and alignment.parquet into the quant directory.
    Library(QuantLibraryArgs),
    /// Measure every library precursor in one run from its MS1 table
    ///
    /// Writes <stem>.extract.parquet. On a DIA study (runs searched with --dia --quant-dir),
    /// scores the precursors nominated in the other runs in this run's spectra instead and
    /// writes dia_pass2/xfer/<stem>.parquet for finalize.
    Extract(QuantExtractArgs),
    /// Transfer q-values and the match-between-runs tables from every run's extract table
    ///
    /// Writes <base>.mbr.parquet and <base>.mbr_validation.parquet. On a TMT study (a
    /// directory of per-run <stem>.tmt.parquet tables from searches without rescoring, and no
    /// library), joins the pooled rescoring (--psms, --decoy-psms) onto every run's reporter
    /// rows, applies the quantification gate again, links the rows to the study's protein
    /// groups (--database) and writes <base>.tmt.parquet.
    Finalize(Box<QuantFinalizeArgs>),
}

impl QuantCommand {
    pub(crate) fn threads(&self) -> Option<usize> {
        match self {
            QuantCommand::Library(a) => a.threads,
            QuantCommand::Extract(a) => a.threads,
            QuantCommand::Finalize(a) => a.threads,
        }
    }
}

#[derive(Args, Debug)]
pub(crate) struct QuantLibraryArgs {
    /// Directory with each run's <stem>.quant.parquet and <stem>.ms1.parquet
    #[arg(long = "quant-dir", value_name = "DIR")]
    pub(crate) quant_dir: PathBuf,

    /// Run stems to include [default: every <stem>.quant.parquet in the directory]
    ///
    /// A run whose PSM or MS1 table is missing is an error.
    #[arg(long = "runs", value_name = "STEM", num_args = 1..)]
    pub(crate) runs: Vec<String>,

    /// Percolator target PSM results (--results-psms) with their score column
    ///
    /// One file over every run (SpecIds distinct across runs), or <stem>=<file> once per
    /// run. PSMs are joined to the run tables on (run, SpecId).
    #[arg(long = "psms", value_name = "TSV|STEM=TSV", num_args = 1.., required = true)]
    pub(crate) psms: Vec<String>,

    /// Percolator decoy PSM results (--decoy-results-psms), in the form of --psms
    #[arg(long = "decoy-psms", value_name = "TSV|STEM=TSV", num_args = 1.., required = true)]
    pub(crate) decoy_psms: Vec<String>,

    /// PSM q-value threshold of a donor: the search's --psm-fdr gives the result of --mbr
    ///
    /// A target PSM at or below it identifies its precursor in its run: those runs are the
    /// precursor's donors and get identified rows; the other runs get transfers.
    #[arg(long = "psm-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) psm_fdr: f64,

    /// Precursor q-value threshold of a library entry [default: inherits --psm-fdr]
    ///
    /// The best-PSM target/decoy q-value of the precursor. The search's --mbr-donor-fdr
    /// gives the result of --mbr.
    #[arg(long = "mbr-donor-fdr", value_name = "FDR", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = ADVANCED_FDR)]
    pub(crate) mbr_donor_fdr: Option<f64>,

    /// PEP threshold of a donor PSM, applied with --psm-fdr [default: none]
    ///
    /// The Percolator posterior_error_prob column. The search's --pep gives the result of
    /// --mbr.
    #[arg(long = "pep", value_name = "PEP", value_parser = parse_unit_fraction,
          help_heading = H_FDR)]
    pub(crate) pep: Option<f64>,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

#[derive(Args, Debug)]
pub(crate) struct QuantExtractArgs {
    /// Directory with the library and the run's <stem>.ms1.parquet
    #[arg(long = "quant-dir", value_name = "DIR")]
    pub(crate) quant_dir: PathBuf,

    /// Stem of the run to measure
    #[arg(long = "run", value_name = "STEM")]
    pub(crate) run: String,

    /// MS1 m/z tolerance of the template fit in ppm, e.g. 10ppm [default: 10ppm]
    ///
    /// The search's --lfq-tol gives the result of --mbr.
    #[arg(long = "lfq-tol", value_name = "TOL", value_parser = parse_precursor_tol)]
    pub(crate) lfq_tol: Option<Tolerance>,

    /// Measure only the transfers that finalize can accept
    ///
    /// true measures a transfer only when the entry's best donor PSM passes the transfer
    /// donor threshold and it has at least --mbr-min-donor-fraction of the other aligned
    /// runs as donors; a skipped pair gets neither its target nor its decoy rows. false
    /// measures every transfer pair. The search's --mbr measures with true.
    #[arg(long = "mbr-preselect", value_name = "BOOL", default_value_t = true,
          action = clap::ArgAction::Set, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) mbr_preselect: bool,

    /// Fewest other aligned donor runs of a measured transfer, as a fraction of the aligned runs (0: no minimum)
    ///
    /// With --mbr-preselect true: the fraction of the library's aligned runs other than the
    /// acceptor, rounded up. A minimum leaves out some transfer fits (and their CPU time),
    /// including a few transfers that would have been accepted.
    #[arg(long = "mbr-min-donor-fraction", value_name = "FRACTION",
          default_value_t = crate::quant::mbr::DEFAULT_MIN_DONOR_FRACTION,
          value_parser = parse_unit_fraction, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) mbr_min_donor_fraction: f64,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

#[derive(Args, Debug)]
pub(crate) struct QuantFinalizeArgs {
    /// Directory with the library and every run's <stem>.extract.parquet (on a TMT study, every run's <stem>.tmt.parquet)
    #[arg(long = "quant-dir", value_name = "DIR")]
    pub(crate) quant_dir: PathBuf,

    /// Output base: writes <base>.mbr.parquet and <base>.mbr_validation.parquet (<base>.tmt.parquet on a TMT study)
    #[arg(short = 'o', long = "output", value_name = "BASE")]
    pub(crate) output: PathBuf,

    /// Quantification FDR: q-value threshold of the transfers
    ///
    /// A quantification-quality q-value, not an identification FDR. The search's
    /// --quant-fdr gives the result of --mbr. In a DIA study, the reported precursor q of the
    /// rows written to <base>.dia_quant.parquet.
    #[arg(long = "quant-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) quant_fdr: f64,

    /// q-value threshold of a transfer reported as transfer_ms1 [default: inherits --quant-fdr]
    #[arg(long = "quant-transfer-fdr", value_name = "FDR", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = ADVANCED_FDR)]
    pub(crate) quant_transfer_fdr: Option<f64>,

    /// Also write the QPX views into this directory (needs --database, --psms and --decoy-psms)
    ///
    /// The protein groups of every run's PSMs as quantms.pg.parquet, and the features of
    /// the groups at --protein-fdr as quantms.feature.parquet, each with its group.
    #[arg(long = "output-parquet", value_name = "DIR",
          requires_all = ["database", "inference_psms", "inference_decoy_psms"])]
    pub(crate) output_parquet: Option<PathBuf>,

    /// The protein database of the search, with its options (as for andes inference)
    #[arg(long, value_name = "FILE")]
    pub(crate) database: Option<PathBuf>,

    /// Percolator target PSM results of every run (as for andes inference)
    ///
    /// One file from one Percolator run over all the runs, or one per run.
    #[arg(id = "inference_psms", long = "psms", value_name = "TSV", num_args = 1..)]
    pub(crate) psms: Vec<PathBuf>,

    /// Percolator decoy PSM results of the same Percolator run
    #[arg(id = "inference_decoy_psms", long = "decoy-psms", value_name = "TSV", num_args = 1..)]
    pub(crate) decoy_psms: Vec<PathBuf>,

    #[command(flatten)]
    pub(crate) proteins: ProteinArgs,

    /// Training rows per transfer or identified-quality classifier
    ///
    /// Above it the classifier is trained on the rows of a hashed sample of (precursor,
    /// run) pairs (every row is still scored); below it every row trains.
    #[arg(long = "max-train-rows", value_name = "N",
          default_value_t = crate::quant::finalize::DEFAULT_MAX_TRAIN_ROWS,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_train_rows: usize,

    /// Members of the bagged transfer classifier (1 to 64)
    ///
    /// Each fold's transfer score is the mean of this many models, each trained on its own
    /// sample of the training rows with its own seed, which keeps the accepted transfers
    /// stable across seeds at hundreds of runs; training and scoring time grow with it.
    #[arg(long = "transfer-models", value_name = "N",
          default_value_t = crate::quant::finalize::DEFAULT_TRANSFER_MODELS,
          value_parser = parse_transfer_models, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) transfer_models: usize,

    /// Reference cells of the run normalisation (about 150 bytes each)
    ///
    /// Above it the fit uses a hashed sample of whole precursors, which coarsens the
    /// RT-local offsets, so keep it above the cells of the study when memory allows.
    #[arg(long = "max-norm-cells", value_name = "N",
          default_value_t = crate::quant::finalize::DEFAULT_MAX_NORM_CELLS,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_norm_cells: usize,

    /// Donor fits summarised in one pass over the extract tables (about 90 bytes each)
    ///
    /// With more, the precursors are summarised in groups, one pass each. The results do
    /// not depend on it.
    #[arg(long = "max-donor-fits", value_name = "N",
          default_value_t = crate::quant::finalize::DEFAULT_MAX_DONOR_FITS,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_donor_fits: usize,

    /// Held-out classifier scores kept in memory (4 bytes each)
    ///
    /// Above it they are sorted and written to temporary tables. The results do not depend
    /// on it.
    #[arg(long = "max-held-scores", value_name = "N",
          default_value_t = crate::quant::finalize::DEFAULT_MAX_HELD_SCORES,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_held_scores: usize,

    /// Directory of the temporary tables, e.g. a local disk [default: the directory of the output base]
    ///
    /// The tables are removed at the end.
    #[arg(long = "temp-dir", value_name = "DIR")]
    pub(crate) temp_dir: Option<PathBuf>,

    /// Seed of the transfer FDR folds and models and of the samples; on a DIA study, of the learned score and the rescoring
    #[arg(long = "seed", value_name = "N", default_value_t = crate::quant::MBR_SEED,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) seed: u64,

    /// Share of the precursors with an identified fit that report the trace area instead of the template area
    ///
    /// The trace area (trace_quantity of the extract) is the baseline-subtracted isotope
    /// traces integrated around the apex; the template area the Gaussian amplitude times
    /// the template width. The share, in [0, 1], goes to the brightest precursors by the
    /// median amplitude signal-to-noise of their identified fits, in every run: 0 reports
    /// the template area everywhere, 1 switches every such precursor. quant_q and the
    /// transfer q-values do not depend on it.
    #[arg(long = "quant-trace-fraction", value_name = "FRACTION",
          default_value_t = crate::quant::finalize::DEFAULT_TRACE_FRACTION,
          value_parser = parse_unit_fraction, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) quant_trace_fraction: f64,

    /// TMT study: the SDRF of the runs, for each run's plex, fraction, instrument and channel samples
    ///
    /// The runs are matched by data file name. A channel is a reference when
    /// characteristics[pooled sample] is pooled or SN=..., or characteristics[biological
    /// replicate] is pooled; blank when a characteristics cell, or the last word of the source
    /// name, is blank or empty. A plex is the files whose channels have the same samples and
    /// technical replicates.
    #[arg(long, value_name = "FILE", help_heading = H_ISOBARIC)]
    pub(crate) sdrf: Option<PathBuf>,

    /// TMT study: PEP threshold of a quantified PSM, applied with --psm-fdr [default: none]
    ///
    /// The pooled Percolator posterior_error_prob column, as the search's --pep. Every
    /// quantified row is weighted by 1 - PEP in quant_weight.
    #[arg(long = "pep", value_name = "PEP", value_parser = parse_unit_fraction,
          help_heading = H_FDR)]
    pub(crate) pep: Option<f64>,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

impl ProteinArgs {
    pub(crate) fn proteoform_fdr(&self) -> f64 {
        proteoform_fdr(self.proteoform_fdr, self.psm_fdr)
    }
}

impl QuantLibraryArgs {
    pub(crate) fn mbr_donor_fdr(&self) -> f64 {
        self.mbr_donor_fdr.unwrap_or(self.psm_fdr)
    }
}

impl QuantFinalizeArgs {
    pub(crate) fn quant_transfer_fdr(&self) -> f64 {
        self.quant_transfer_fdr.unwrap_or(self.quant_fdr)
    }
}
