//! Activation methods: the EThcD/ETciD policy (label as EThcD, or as ETD with
//! `--ethcd-activation etd`) and the Thermo filter-string check for supplemental activation.

pub(crate) fn log_ethcd_once() {
    use std::sync::atomic::{AtomicBool, Ordering};
    static WARNED: AtomicBool = AtomicBool::new(false);
    if WARNED
        .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
        .is_ok()
    {
        if ethcd_as_etd() {
            log::info!(
                "EThcD/ETciD detected; --ethcd-activation etd is set, so these spectra \
                 stay classified as ETD (c/z stack active)."
            );
        } else {
            log::info!(
                "EThcD/ETciD detected (electron transfer + supplemental collisional \
                 activation); these spectra are labelled EThcD and scored with b/y and c/z \
                 ions. IN --glyco the c/z stack (c/z generation and scoring, the c/z \
                 truncation gate, the ETD DB fallback, --glyco-hcd-pair) runs only on ETD \
                 spectra; pass --ethcd-activation etd to label these spectra ETD."
            );
        }
    }
}

pub(super) fn activation_user_param_accession(name: &str) -> Option<&'static str> {
    match name.trim().to_ascii_lowercase().as_str() {
        "cid" | "collision-induced dissociation" => Some(CV_CID),
        "hcd" | "beam-type collision-induced dissociation" => Some(CV_HCD),
        "etd" | "electron transfer dissociation" => Some(CV_ETD),
        _ => None,
    }
}

pub fn filter_has_supplemental_collision(filter: &str) -> bool {
    let head = filter.split('[').next().unwrap_or("");
    let mut tokens = head.split_whitespace();
    if !tokens.any(|t| {
        let t = t.to_ascii_lowercase();
        t.len() > 2 && t.starts_with("ms") && t[2..].bytes().all(|b| b.is_ascii_digit())
    }) {
        return false;
    }
    let stages: Vec<&str> = tokens.filter(|t| t.contains('@')).collect();
    let Some(last) = stages.last() else {
        return false;
    };
    let precursor = |t: &str| t.split('@').next().unwrap_or("").to_string();
    let last_precursor = precursor(last);
    let methods: Vec<String> = stages
        .iter()
        .rev()
        .take_while(|t| precursor(t) == last_precursor)
        .flat_map(|t| t.split('@').skip(1))
        .map(|a| {
            a.chars()
                .take_while(char::is_ascii_alphabetic)
                .collect::<String>()
                .to_ascii_lowercase()
        })
        .collect();
    let electron = methods
        .iter()
        .any(|m| matches!(m.as_str(), "etd" | "ecd" | "netd"));
    let collisional = methods.iter().any(|m| matches!(m.as_str(), "hcd" | "cid"));
    electron && collisional
}

static ETHCD_AS_ETD: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

static ETHCD_POLICY_READ: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

pub fn init_ethcd_as_etd(as_etd: bool) {
    debug_assert!(
        !ETHCD_POLICY_READ.load(std::sync::atomic::Ordering::Relaxed) || !as_etd,
        "init_ethcd_as_etd(true) called after the policy was already read: every \
         spectrum read before this point was labelled under the default (EThcD)"
    );
    let _ = ETHCD_AS_ETD.set(as_etd);
}

#[inline]
pub fn ethcd_as_etd_policy() -> bool {
    ethcd_as_etd()
}

#[inline]
pub(super) fn ethcd_as_etd() -> bool {
    ETHCD_POLICY_READ.store(true, std::sync::atomic::Ordering::Relaxed);
    *ETHCD_AS_ETD.get().unwrap_or(&false)
}
