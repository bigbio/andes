//! False discovery rate control beyond the PSM level: class-specific q-values, the
//! knockoff filter, site-level FDR and the identity of a site across PSMs.

pub mod class;
pub mod entrapment;
pub mod knockoff;
pub mod sites;
pub mod identity;
