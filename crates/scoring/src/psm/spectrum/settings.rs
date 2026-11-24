//! Process-wide scoring settings: the peak filter, set once from the command line or the
//! environment of a test.

use std::sync::OnceLock;

pub(super) enum PeakFilterEnv {
    Disabled,
    Override(f64, usize),
    Unset,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ScoringSettings {
    pub peak_filter: Option<(f64, usize)>,
}

static SCORING_SETTINGS: OnceLock<ScoringSettings> = OnceLock::new();

pub fn init_scoring_settings(settings: ScoringSettings) {
    let _ = SCORING_SETTINGS.set(settings);
}

#[inline]
pub(crate) fn scoring_settings() -> ScoringSettings {
    SCORING_SETTINGS.get().copied().unwrap_or_default()
}

pub(super) fn peak_filter_env() -> &'static PeakFilterEnv {
    static CACHE: OnceLock<PeakFilterEnv> = OnceLock::new();
    CACHE.get_or_init(|| {
        let (w, k) = match scoring_settings().peak_filter {
            Some((w, k)) => (Some(w), Some(k)),
            None => (None, None),
        };
        match (w, k) {
            (Some(w), _) if w <= 0.0 => PeakFilterEnv::Disabled,
            (Some(w), Some(k)) => PeakFilterEnv::Override(w, k),
            _ => PeakFilterEnv::Unset,
        }
    })
}
