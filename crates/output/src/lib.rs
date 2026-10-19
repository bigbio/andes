//! Output writers: Percolator PIN, TSV.
//!
//! # Known column behaviors
//!
//! * **FragMethod**: emitted via `ActivationMethod::name()` (e.g. `"HCD"`,
//!   `"CID"`). Unknown activation is written as `"UNKNOWN"`.
//!
//! * **IsotopeError**: the precursor-matching loop tries multiple isotope
//!   offsets but does not record *which* offset produced the match. The TSV
//!   column is always written as `0`. Will be fixed once the winning
//!   isotope offset is threaded into `PsmMatch`.
//!
//! * **Decoy filtering**: this writer emits decoy PSMs with the decoy
//!   prefix preserved in the Protein column; downstream Percolator handles
//!   decoy labelling.
//!
//! * **QValue / PepQValue**: Not emitted; TDA columns are not currently
//!   produced.

pub mod tsv;
pub use tsv::{write_tsv, write_tsv_to};

pub mod pin;
pub use pin::{format_spec_id, write_pin, write_pin_to};

pub mod stats;
pub use stats::{write_statistics_log, RunStatistics};

pub mod qpx;
pub use qpx::write_qpx;

pub mod rt_wiring;
pub use rt_wiring::populate_rt_features;

pub mod provenance;
pub use provenance::populate_provenance_features;

pub mod class_fdr;
pub use class_fdr::{apply_class_fdr, ClassCount, Cutoffs};

pub mod isobaric_sharing;
pub use isobaric_sharing::share_isobaric_references;

pub mod knockoff_fdr;
pub use knockoff_fdr::{apply_knockoff_fdr, KnockoffCount, KnockoffOutcome};

pub mod proteoform_table;
pub use proteoform_table::{proteoform_table, write_proteoform_tsv, QuantSummary};

pub mod glyco_rt;
pub use glyco_rt::populate_glyco_rt_features;

pub mod percolator;
pub use percolator::{
    decoy_results_path, parse_psm_results, resolve_backend, run_percolator,
    run_percolator_with_decoys, PercolatorBackend, PercolatorError, PercolatorPsm,
    DEFAULT_PERCOLATOR_IMAGE,
};

pub(crate) mod percolator_enz;
pub(crate) mod row_context;
pub use row_context::emitted_rows;

pub mod glyco_pin;
pub use glyco_pin::write_glyco_pin;

pub mod quant_out;
pub use quant_out::{
    write_lfq_features_tsv, write_lfq_tsv, write_tmt_tsv, LfqRow, ModRecord, QuantId, TmtRow,
};

pub mod run_quant_table;

pub mod extract_table;
pub mod library_table;
pub mod ms1_table;
pub mod ms2_table;

pub mod feature_parquet;
pub mod mbr_out;
pub use feature_parquet::{write_feature_parquet, FeatureRecord, ProteinPosition};
