//! Labelled-pair quantification (SILAC, dimethyl and other mass-shift labels): label
//! specifications, the joint fit of the channels' isotope envelopes, and channel classification.

mod classify;
mod fit;
mod spec;

pub use classify::{
    channel_qvalues, channel_score, classify_channel, classify_identified, classify_transferred,
    coelution_factor, decoy_min_sep, decoy_shifts, learn_rt_shifts, log2_ratio, mass_factor,
    median, weighted_qvalues, ChannelStatus, ClassifyParams,
};
pub use fit::{
    fit_group, ChannelEnvelope, ChannelFit, FitParams, GroupFit, GroupProbe, DEFAULT_WIDTHS,
    FREE_WIDTHS,
};
pub use spec::{Channel, LabelAtoms, LabelSite, LabelSpec, Site, PRESETS};

pub const D_H: f64 = model::mass::D_H;
pub const LOWEST_OFFSET: i32 = -2;
pub const HIGHEST_OFFSET: i32 = 5;
const N_OFFSETS: usize = (HIGHEST_OFFSET - LOWEST_OFFSET + 1) as usize;
const READ_OFFSETS: [i32; 5] = [-1, 0, 1, 2, 3];
pub const MAX_CHANNELS: usize = 5;
pub const MAX_COMPONENTS: usize = 2 * MAX_CHANNELS;

#[cfg(test)]
#[path = "label/tests/label.rs"]
mod tests;
