//! Writers and tables for the results of a run: PIN files, Parquet and QPX views, atomic file
//! writes, FDR and protein-inference outputs and run provenance. It sits at the end of the
//! pipeline, turning search, rescoring and quantification results into files.

mod atomic;
pub use atomic::write_atomic;
pub mod tsv;
pub use tsv::{write_tsv, write_tsv_to};

pub mod pin;
pub use pin::{emitted_row_candidates, emitted_row_evidence, format_spec_id, spec_id_run, write_pin, write_pin_to};

pub mod statistics;
pub use statistics::{
    classify_mods, fact_rows, psm_facts, psm_facts_of_rows, write_statistics_log, FactRow,
    FactRows, FixedModSpec, IdentificationSummary, PeptideMods, PsmFacts, RunStatistics,
    SpectraTally, SummaryContext,
};

pub mod report;
pub use report::{
    CountRow, InputFile, ParamRow, Reference, RunReport, StageRow, RUN_REPORT_SCHEMA_MAJOR,
};

pub mod qpx;
pub use qpx::{omitted_psm_rows, write_qpx};

pub mod rt;
pub use rt::{populate_rt_features, populate_rt_features_per_run};

pub mod provenance;
pub use provenance::populate_provenance_features;
pub use fdr::class::{apply_class_fdr, ClassCount, Cutoffs};

pub mod isobaric;
pub use isobaric::share_isobaric_references;

pub mod ipf;
mod spectra;
pub use fdr::knockoff::{apply_knockoff_fdr, KnockoffCount, KnockoffOutcome};
pub use fdr::sites::{site_claims, write_sites_tsv, SiteCount, SiteOutcome, SiteRow};
pub use protein::groups::{write_pg_parquet, FeatureCounts};
pub use protein::quantities::ProteinQuant;
pub use protein::inference::{
    infer_protein_groups, write_protein_groups_tsv, GroupAt, Inference, ProteinGroup,
};
pub use protein::proteoforms::{proteoform_table, write_proteoform_tsv, QuantSummary};
pub use glycopeptides::rt::populate_glyco_rt_features;

pub mod offsets;
pub mod percolator;
pub use percolator::{
    apply_per_run_q, decoy_results_path, parse_psm_results, pep_is_degenerate,
    percolator_args_with_threads, resolve_backend, run_percolator, run_percolator_pooled,
    run_percolator_with_decoys, score_all_rows, PercolatorBackend, PercolatorError, PercolatorPsm,
    PsmResultReader, PsmResultRow, RowScores, read_psm_results, DEFAULT_PERCOLATOR_IMAGE,
};
pub use rescore::{run as run_native_percolator, MemLimit, Sizing};
pub(crate) mod rows;
pub use rows::emitted_rows;
pub use glycopeptides::pin::write_glyco_pin;

pub mod quantities;
pub use quantities::{
    read_label_parquet, read_lfq_parquet, read_tmt_parquet, write_label_parquet, write_lfq_parquet,
    write_tmt_parquet, write_tmt_parquet_mbr, LabelRow, LabelTableMeta, LfqRow, LfqTableMeta,
    ModRecord, QuantId, TmtRow, TransferInfo, UnmixMeta,
};

pub mod extract;
pub mod decoys;
pub mod tmt;
pub mod library;
pub mod ms1;
pub mod ms2;

pub mod features;
pub mod mbr;
pub mod versioned;
pub mod crosslinks;
pub mod fdr;
pub mod protein;
pub mod glycopeptides;
pub mod dia;
pub use features::{
    write_feature_parquet, FeatureGroup, FeatureParquetWriter, FeatureRecord, ProteinPosition,
};
pub use versioned::ReadTable;
