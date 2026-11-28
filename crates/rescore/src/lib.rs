//! Native rescoring engine: the Percolator algorithm (Käll et al., Nat. Methods 2007;4:923-925,
//! https://doi.org/10.1038/nmeth1113) run in process. A linear SVM is trained semi-supervised on
//! target and decoy PSMs with cross-validation, and its scores give q-values and PEPs. The engine
//! tracks upstream Percolator master (https://github.com/percolator/percolator, commit
//! a9783ee8) and writes Percolator's result files.

mod args;
mod budget;
mod codec;
mod cv;
mod direction;
mod folds;
mod merge;
mod normalize;
mod pep;
mod pin;
mod pread;
mod run;
mod scores;
mod store;
mod svm;
mod text;
mod write;

/// A PSM's row number in the PIN. A PIN holds fewer than 2^32 PSMs, so index lists that scale
/// with the PSM count store 4 bytes a row.
pub type RowId = u32;

pub use budget::{Budget, Reservation};
pub use cv::{train, CvConfig, CvResult};
pub use direction::init_direction;
pub use folds::{fold_groups, folds_by_precursor, folds_by_spectrum};
pub use merge::{merge_folds, normalize_fold_scores};
pub use normalize::{raw_weights, score, stdv, Norm};
pub use pep::{estimate_pi0, pep_from_q, too_good_separation};
pub use pin::{read_pin, sequence_key, MemLimit, Pin, RowMeta, Sizing};
pub use run::run;
pub use scores::{count_at, mix_max_q, tdc_q, tdc_winners, training_sets};
pub use store::Store;
pub use svm::{l2_svm_mfn, MfnOpts, Rows};
pub use text::{RowText, Texts};
pub use write::{fmt_fixed, fmt_g, write_results, write_weights, Results};
