use std::collections::BTreeSet;

use model::enzyme::Enzyme;

#[derive(Debug, Clone, PartialEq)]
pub struct SelectionEntry {
    pub model_id: String,
    pub activation: String,
    pub instrument: String,
    pub enzyme: String,
    pub experiment_class: BTreeSet<String>,
}

#[derive(Debug, Clone)]
pub struct SelectionKey {
    pub activation: String,
    pub instrument: String,
    pub enzyme: String,
    pub experiment_class: BTreeSet<String>,
}

pub const RELAXED_ENZYME: &str = Enzyme::Trypsin.name();

const LABELING_TAGS: &[&str] = &["tmt", "itraq"];

fn is_labeling_tag(s: &str) -> bool {
    LABELING_TAGS.contains(&s)
}

fn contains_labeling_tag(set: &BTreeSet<String>) -> bool {
    set.iter().any(|s| is_labeling_tag(s.as_str()))
}

fn try_steps_1_to_3<'a>(
    entries: &'a [SelectionEntry],
    activation: &str,
    instrument: &str,
    enzyme: &str,
    key_class: &BTreeSet<String>,
) -> Option<&'a str> {
    let candidates: Vec<&SelectionEntry> = entries
        .iter()
        .filter(|e| e.activation == activation && e.instrument == instrument && e.enzyme == enzyme)
        .collect();

    if candidates.is_empty() {
        return None;
    }

    if let Some(hit) = candidates
        .iter()
        .filter(|e| &e.experiment_class == key_class)
        .min_by_key(|e| &e.model_id)
    {
        return Some(&hit.model_id);
    }

    {
        let subset_candidates: Vec<&&SelectionEntry> = candidates
            .iter()
            .filter(|e| {
                !e.experiment_class.is_empty()
                    && e.experiment_class != *key_class
                    && e.experiment_class.is_subset(key_class)
            })
            .collect();

        if !subset_candidates.is_empty() {
            let max_size = subset_candidates
                .iter()
                .map(|e| e.experiment_class.len())
                .max()
                .unwrap_or(0);

            let largest: Vec<&&&SelectionEntry> = subset_candidates
                .iter()
                .filter(|e| e.experiment_class.len() == max_size)
                .collect();

            let winner = largest
                .iter()
                .max_by_key(|e| {
                    let labeling_priority = if contains_labeling_tag(&e.experiment_class) {
                        1u8
                    } else {
                        0u8
                    };
                    (labeling_priority, std::cmp::Reverse(e.model_id.clone()))
                })
                .unwrap();

            return Some(&winner.model_id);
        }
    }

    for tag in LABELING_TAGS {
        if key_class.contains(*tag) {
            let tag_set: BTreeSet<String> = std::iter::once((*tag).to_string()).collect();
            if let Some(hit) = candidates
                .iter()
                .filter(|e| e.experiment_class == tag_set)
                .min_by_key(|e| &e.model_id)
            {
                return Some(&hit.model_id);
            }
        }
    }

    None
}

fn try_step_4<'a>(
    entries: &'a [SelectionEntry],
    activation: &str,
    instrument: &str,
    enzyme: &str,
) -> Option<&'a str> {
    entries
        .iter()
        .filter(|e| {
            e.activation == activation
                && e.instrument == instrument
                && e.enzyme == enzyme
                && e.experiment_class.is_empty()
        })
        .min_by_key(|e| &e.model_id)
        .map(|e| e.model_id.as_str())
}

pub fn select<'a>(
    entries: &'a [SelectionEntry],
    key: &SelectionKey,
    instrument_family: impl Fn(&str) -> String,
    generic_id: Option<&'a str>,
) -> Option<&'a str> {
    let family = instrument_family(&key.instrument);
    let has_family_fallback = family != key.instrument;

    if let Some(id) = try_steps_1_to_3(
        entries,
        &key.activation,
        &key.instrument,
        &key.enzyme,
        &key.experiment_class,
    ) {
        return Some(id);
    }

    if has_family_fallback {
        if let Some(id) = try_steps_1_to_3(
            entries,
            &key.activation,
            &family,
            &key.enzyme,
            &key.experiment_class,
        ) {
            return Some(id);
        }
    }

    if let Some(id) = try_step_4(entries, &key.activation, &key.instrument, &key.enzyme) {
        return Some(id);
    }

    if has_family_fallback {
        if let Some(id) = try_step_4(entries, &key.activation, &family, &key.enzyme) {
            return Some(id);
        }
    }

    generic_id
}

pub fn select_nearest<'a>(
    entries: &'a [SelectionEntry],
    key: &SelectionKey,
    instrument_family: impl Fn(&str) -> String + Copy,
    standard_id: &'a str,
) -> (&'a str, bool) {
    select_nearest_for(entries, key, None, instrument_family, standard_id)
}

pub fn select_nearest_for<'a>(
    entries: &'a [SelectionEntry],
    key: &SelectionKey,
    detected_instrument: Option<&str>,
    instrument_family: impl Fn(&str) -> String + Copy,
    standard_id: &'a str,
) -> (&'a str, bool) {
    let detected = detected_instrument.filter(|d| *d != key.instrument);
    if let Some(d) = detected {
        let k = SelectionKey {
            instrument: d.to_string(),
            ..key.clone()
        };
        if let Some(id) = select(entries, &k, instrument_family, None) {
            return (id, false);
        }
    }
    if let Some(id) = select(entries, key, instrument_family, None) {
        return (id, false);
    }
    let siblings = resolution_class_siblings(&key.instrument);
    for inst in siblings {
        let k = SelectionKey {
            instrument: (*inst).to_string(),
            ..key.clone()
        };
        if let Some(id) = select(entries, &k, instrument_family, None) {
            return (id, true);
        }
    }
    if key.enzyme != RELAXED_ENZYME {
        let tryptic = SelectionKey {
            enzyme: RELAXED_ENZYME.to_string(),
            ..key.clone()
        };
        if let Some(d) = detected {
            let k = SelectionKey {
                instrument: d.to_string(),
                ..tryptic.clone()
            };
            if let Some(id) = select(entries, &k, instrument_family, None) {
                return (id, true);
            }
        }
        if let Some(id) = select(entries, &tryptic, instrument_family, None) {
            return (id, true);
        }
        for inst in siblings {
            let k = SelectionKey {
                instrument: (*inst).to_string(),
                ..tryptic.clone()
            };
            if let Some(id) = select(entries, &k, instrument_family, None) {
                return (id, true);
            }
        }
    }
    if let Some(pick) = entries
        .iter()
        .filter(|e| e.activation == key.activation)
        .min_by(|a, b| {
            relaxed_instrument_rank(a, key)
                .cmp(&relaxed_instrument_rank(b, key))
                .then_with(|| a.model_id.cmp(&b.model_id))
        })
    {
        return (pick.model_id.as_str(), true);
    }
    (standard_id, true)
}

fn relaxed_instrument_rank(e: &SelectionEntry, key: &SelectionKey) -> (u8, bool, bool) {
    let labels = |set: &BTreeSet<String>| -> BTreeSet<String> {
        set.iter().filter(|s| is_labeling_tag(s)).cloned().collect()
    };
    let entry_labels = labels(&e.experiment_class);
    let labelling = if entry_labels == labels(&key.experiment_class) {
        0
    } else if entry_labels.is_empty() {
        1
    } else {
        2
    };
    (
        labelling,
        e.experiment_class != key.experiment_class,
        e.enzyme != key.enzyme,
    )
}

pub fn resolution_class_siblings(instrument: &str) -> &'static [&'static str] {
    match instrument {
        "QExactive" => &["HighRes"],
        "HighRes" => &["QExactive"],
        "TOF" | "TimsTOF" | "OrbitrapAstral" => &["QExactive", "HighRes"],
        _ => &[],
    }
}

pub fn parse_experiment_class(s: &str) -> BTreeSet<String> {
    let trimmed = s.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("standard") {
        return BTreeSet::new();
    }
    trimmed
        .split('+')
        .map(|tok| tok.trim().to_ascii_lowercase())
        .filter(|tok| !tok.is_empty())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn etd_query_prefers_an_etd_model_at_another_resolution_over_an_hcd_base() {
        let entries = vec![
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
            e("etd_highres_tryp", "ETD", "QExactive", "Trypsin", &[]),
        ];
        let k = key("ETD", "LowRes", "Trypsin", &[]);
        let (id, substituted) = select_nearest(&entries, &k, fam, "hcd_qexactive_tryp");
        assert_eq!(
            id, "etd_highres_tryp",
            "an ETD query must not land on a b/y model"
        );
        assert!(substituted, "an inexact match must still warn");
    }

    #[test]
    fn instrument_relaxation_keeps_the_labelling_and_class() {
        let entries = vec![
            e("etd_a_tmt", "ETD", "QExactive", "Trypsin", &["tmt"]),
            e("etd_b_phospho", "ETD", "QExactive", "Trypsin", &["phospho"]),
            e("etd_c", "ETD", "QExactive", "LysC", &[]),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
        ];
        let pick = |ec: &[&str]| {
            select_nearest(
                &entries,
                &key("ETD", "LowRes", "Trypsin", ec),
                fam,
                "hcd_qexactive_tryp",
            )
        };
        assert_eq!(
            pick(&[]),
            ("etd_c", true),
            "unlabelled, standard class first"
        );
        assert_eq!(pick(&["phospho"]), ("etd_b_phospho", true));
        assert_eq!(
            pick(&["tmt"]),
            ("etd_a_tmt", true),
            "a labelled query keeps its label"
        );
        assert_eq!(pick(&["itraq"]), ("etd_b_phospho", true));

        let only_tmt = vec![
            e("etd_a_tmt", "ETD", "QExactive", "Trypsin", &["tmt"]),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
        ];
        let k = key("ETD", "LowRes", "Trypsin", &[]);
        assert_eq!(
            select_nearest(&only_tmt, &k, fam, "hcd_qexactive_tryp"),
            ("etd_a_tmt", true)
        );

        let enz = vec![
            e("etd_a", "ETD", "QExactive", "LysC", &[]),
            e("etd_b", "ETD", "QExactive", "Trypsin", &[]),
        ];
        let k = key("ETD", "LowRes", "Trypsin", &[]);
        assert_eq!(select_nearest(&enz, &k, fam, "x"), ("etd_b", true));
    }

    #[test]
    fn instrument_relaxation_does_not_disturb_a_query_that_already_resolves() {
        let entries = vec![
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
            e("etd_highres_tryp", "ETD", "QExactive", "Trypsin", &[]),
        ];
        let k = key("HCD", "QExactive", "Trypsin", &[]);
        let (id, substituted) = select_nearest(&entries, &k, fam, "hcd_qexactive_tryp");
        assert_eq!(id, "hcd_qexactive_tryp");
        assert!(
            !substituted,
            "an exact hit must not be reported as substituted"
        );
    }

    #[test]
    fn falls_back_to_the_standard_base_when_no_model_shares_the_activation() {
        let entries = vec![e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[])];
        let k = key("UVPD", "LowRes", "Trypsin", &[]);
        let (id, substituted) = select_nearest(&entries, &k, fam, "hcd_qexactive_tryp");
        assert_eq!(id, "hcd_qexactive_tryp");
        assert!(substituted);
    }

    fn set(items: &[&str]) -> BTreeSet<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    fn e(id: &str, a: &str, i: &str, en: &str, ec: &[&str]) -> SelectionEntry {
        SelectionEntry {
            model_id: id.into(),
            activation: a.into(),
            instrument: i.into(),
            enzyme: en.into(),
            experiment_class: set(ec),
        }
    }

    fn key(a: &str, i: &str, en: &str, ec: &[&str]) -> SelectionKey {
        SelectionKey {
            activation: a.into(),
            instrument: i.into(),
            enzyme: en.into(),
            experiment_class: set(ec),
        }
    }

    fn fam(i: &str) -> String {
        match i {
            "OrbitrapAstral" => "QExactive",
            "TimsTOF" => "TOF",
            x => x,
        }
        .to_string()
    }

    fn manifest() -> Vec<SelectionEntry> {
        vec![
            e("astral", "HCD", "OrbitrapAstral", "Trypsin", &[]),
            e("qe_tmt", "HCD", "QExactive", "Trypsin", &["tmt"]),
            e("qe", "HCD", "QExactive", "Trypsin", &[]),
            e("generic", "HCD", "QExactive", "Trypsin", &[]),
        ]
    }

    #[test]
    fn exact_match_tie_break_is_deterministic() {
        let a = vec![
            e("zzz_model", "HCD", "QExactive", "Trypsin", &[]),
            e("aaa_model", "HCD", "QExactive", "Trypsin", &[]),
        ];
        let mut b = a.clone();
        b.reverse();
        let k = key("HCD", "QExactive", "Trypsin", &[]);
        assert_eq!(select(&a, &k, fam, None), Some("aaa_model"));
        assert_eq!(
            select(&a, &k, fam, None),
            select(&b, &k, fam, None),
            "selection must not depend on manifest row order"
        );
    }

    #[test]
    fn exact_instrument_hit() {
        assert_eq!(
            select(
                &manifest(),
                &key("HCD", "OrbitrapAstral", "Trypsin", &[]),
                fam,
                Some("generic")
            ),
            Some("astral")
        );
    }

    #[test]
    fn backs_off_instrument_family_then_experiment() {
        assert_eq!(
            select(
                &manifest(),
                &key("HCD", "OrbitrapAstral", "Trypsin", &["tmt"]),
                fam,
                Some("generic")
            ),
            Some("qe_tmt")
        );
    }

    #[test]
    fn falls_back_to_generic() {
        assert_eq!(
            select(
                &manifest(),
                &key("ETD", "LowRes", "Trypsin", &[]),
                fam,
                Some("generic")
            ),
            Some("generic")
        );
    }

    #[test]
    fn largest_subset_when_no_exact_combo() {
        let mut m = manifest();
        m.push(e("qe_phos", "HCD", "QExactive", "Trypsin", &["phospho"]));
        assert_eq!(
            select(
                &m,
                &key("HCD", "QExactive", "Trypsin", &["phospho", "tmt"]),
                fam,
                Some("generic")
            ),
            Some("qe_tmt")
        );
    }

    #[test]
    fn no_generic_id_returns_none() {
        assert_eq!(
            select(
                &manifest(),
                &key("ETD", "LowRes", "Trypsin", &[]),
                fam,
                None
            ),
            None
        );
    }

    #[test]
    fn family_fallback_empty_class() {
        assert_eq!(
            select(
                &manifest(),
                &key("HCD", "TimsTOF", "Trypsin", &[]),
                fam,
                Some("generic")
            ),
            Some("generic")
        );
    }

    #[test]
    fn step3_labeling_tag_alone() {
        let m = vec![
            e("qe_tmt", "HCD", "QExactive", "Trypsin", &["tmt"]),
            e("qe", "HCD", "QExactive", "Trypsin", &[]),
        ];
        assert_eq!(
            select(
                &m,
                &key("HCD", "QExactive", "Trypsin", &["phospho", "tmt"]),
                fam,
                Some("generic")
            ),
            Some("qe_tmt")
        );
    }

    #[test]
    fn step4_empty_class_fallback() {
        let m = vec![e("qe", "HCD", "QExactive", "Trypsin", &[])];
        assert_eq!(
            select(
                &m,
                &key("HCD", "QExactive", "Trypsin", &["phospho"]),
                fam,
                Some("generic")
            ),
            Some("qe")
        );
    }

    #[test]
    fn generic_highres_hcd_prefers_qexactive_over_alphabetical_astral() {
        let m = vec![
            e("hcd_astral_tryp", "HCD", "OrbitrapAstral", "Trypsin", &[]),
            e(
                "hcd_highres_nocleavage",
                "HCD",
                "HighRes",
                "NoCleavage",
                &[],
            ),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
            e("cid_lowres_tryp", "CID", "LowRes", "Trypsin", &[]),
        ];
        let (id, sub) = select_nearest(
            &m,
            &key("HCD", "HighRes", "Trypsin", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id, "hcd_qexactive_tryp");
        assert!(
            sub,
            "a resolution-class swap is still a substitution and must WARN"
        );
        let m2 = vec![e("hcd_highres_tryp_x", "HCD", "HighRes", "Trypsin", &[])];
        let (id2, _) = select_nearest(
            &m2,
            &key("HCD", "QExactive", "Trypsin", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id2, "hcd_highres_tryp_x");
        let (id3, sub3) = select_nearest(
            &m,
            &key("HCD", "QExactive", "Trypsin", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!((id3, sub3), ("hcd_qexactive_tryp", false));
    }

    #[test]
    fn dropped_cid_phospho_routes_to_cid_tryp() {
        let m = vec![
            e("cid_lowres_tryp", "CID", "LowRes", "Trypsin", &[]),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
        ];
        let (id, sub) = select_nearest(
            &m,
            &key("CID", "LowRes", "Trypsin", &["phospho"]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id, "cid_lowres_tryp");
        assert!(
            !sub,
            "empty-class fallback is the existing ladder, not a nearest substitution"
        );
    }

    #[test]
    fn enzyme_backoff_to_tryp() {
        let m = vec![e("cid_lowres_tryp", "CID", "LowRes", "Trypsin", &[])];
        let (id, sub) = select_nearest(
            &m,
            &key("CID", "LowRes", "GluC", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id, "cid_lowres_tryp");
        assert!(sub);
    }

    #[test]
    fn non_tryptic_query_relaxes_to_the_tryptic_model_of_its_activation_and_instrument() {
        let m = vec![
            e("cid_lowres_argc", "CID", "LowRes", "ArgC", &[]),
            e("cid_lowres_tryp", "CID", "LowRes", "Trypsin", &[]),
            e("hcd_astral_tryp", "HCD", "OrbitrapAstral", "Trypsin", &[]),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
            e(
                "hcd_qexactive_tryp_tmt",
                "HCD",
                "QExactive",
                "Trypsin",
                &["tmt"],
            ),
        ];
        let cases = [
            (key("CID", "LowRes", "Chymotrypsin", &[]), "cid_lowres_tryp"),
            (key("HCD", "QExactive", "LysC", &[]), "hcd_qexactive_tryp"),
            (
                key("HCD", "QExactive", "LysC", &["tmt"]),
                "hcd_qexactive_tryp_tmt",
            ),
            (key("HCD", "TOF", "NonSpecific", &[]), "hcd_qexactive_tryp"),
        ];
        for (k, want) in cases {
            let (id, sub) = select_nearest(&m, &k, fam, "hcd_qexactive_tryp");
            assert_eq!(id, want, "{k:?}");
            assert!(sub, "{k:?}");
        }
        assert_eq!(RELAXED_ENZYME, Enzyme::Trypsin.name());
    }

    #[test]
    fn detected_instrument_comes_before_its_family() {
        let m = vec![
            e("hcd_astral_tryp", "HCD", "OrbitrapAstral", "Trypsin", &[]),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
        ];
        let tryptic = key("HCD", "QExactive", "Trypsin", &[]);
        assert_eq!(
            select_nearest_for(
                &m,
                &tryptic,
                Some("OrbitrapAstral"),
                fam,
                "hcd_qexactive_tryp"
            ),
            ("hcd_astral_tryp", false)
        );
        assert_eq!(
            select_nearest(&m, &tryptic, fam, "hcd_qexactive_tryp"),
            ("hcd_qexactive_tryp", false)
        );
        let lysc = key("HCD", "QExactive", "LysC", &[]);
        assert_eq!(
            select_nearest_for(&m, &lysc, Some("OrbitrapAstral"), fam, "hcd_qexactive_tryp"),
            ("hcd_astral_tryp", true)
        );
        assert_eq!(
            select_nearest_for(&m, &lysc, Some("TimsTOF"), fam, "hcd_qexactive_tryp"),
            ("hcd_qexactive_tryp", true)
        );
    }

    #[test]
    fn same_enzyme_in_the_resolution_class_beats_the_tryptic_model() {
        let m = vec![
            e(
                "hcd_highres_nocleavage",
                "HCD",
                "HighRes",
                "NoCleavage",
                &[],
            ),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
        ];
        let (id, sub) = select_nearest(
            &m,
            &key("HCD", "QExactive", "NoCleavage", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id, "hcd_highres_nocleavage");
        assert!(sub);
    }

    #[test]
    fn no_enzyme_match_falls_back_to_standard_not_other_activation() {
        let m = vec![
            e("etd_highres_tryp", "ETD", "HighRes", "Trypsin", &[]),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
        ];
        let (id, sub) = select_nearest(
            &m,
            &key("HCD", "HighRes", "Trypsin", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id, "hcd_qexactive_tryp");
        assert!(sub);
    }

    #[test]
    fn falls_back_to_standard_base() {
        let m = vec![e("etd_highres_tryp", "ETD", "HighRes", "Trypsin", &[])];
        let (id, sub) = select_nearest(
            &m,
            &key("UVPD", "TimsTOF", "ArgC", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id, "hcd_qexactive_tryp");
        assert!(sub);
    }

    #[test]
    fn tof_without_a_tof_model_stays_high_res() {
        let m = vec![
            e("astral", "HCD", "OrbitrapAstral", "Trypsin", &[]),
            e("cid_lowres_tryp", "CID", "LowRes", "Trypsin", &[]),
            e("hcd_lowres_like", "HCD", "LowRes", "Trypsin", &[]),
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
            e(
                "hcd_qexactive_tryp_tmt",
                "HCD",
                "QExactive",
                "Trypsin",
                &["tmt"],
            ),
        ];
        for inst in ["TOF", "TimsTOF"] {
            let (id, sub) = select_nearest(
                &m,
                &key("HCD", inst, "Trypsin", &[]),
                fam,
                "hcd_qexactive_tryp",
            );
            assert_eq!(id, "hcd_qexactive_tryp", "{inst}");
            assert!(sub);
            let (id, sub) = select_nearest(
                &m,
                &key("HCD", inst, "Trypsin", &["tmt"]),
                fam,
                "hcd_qexactive_tryp",
            );
            assert_eq!(id, "hcd_qexactive_tryp_tmt", "{inst} TMT");
            assert!(sub);
        }
    }

    #[test]
    fn timstof_model_is_an_exact_hit() {
        let m = vec![
            e("hcd_qexactive_tryp", "HCD", "QExactive", "Trypsin", &[]),
            e("hcd_timstof_tryp", "HCD", "TimsTOF", "Trypsin", &[]),
        ];
        let (id, sub) = select_nearest(
            &m,
            &key("HCD", "TimsTOF", "Trypsin", &[]),
            fam,
            "hcd_qexactive_tryp",
        );
        assert_eq!(id, "hcd_timstof_tryp");
        assert!(!sub);
    }

    #[test]
    fn resolution_class_siblings_never_cross_to_low_res() {
        for inst in ["QExactive", "HighRes", "TOF", "TimsTOF", "OrbitrapAstral"] {
            assert!(
                !resolution_class_siblings(inst).contains(&"LowRes"),
                "{inst}"
            );
        }
        assert!(resolution_class_siblings("LowRes").is_empty());
    }

    #[test]
    fn exact_hit_not_substituted() {
        let m = manifest();
        let (id, sub) = select_nearest(
            &m,
            &key("HCD", "OrbitrapAstral", "Trypsin", &[]),
            fam,
            "generic",
        );
        assert_eq!(id, "astral");
        assert!(!sub);
    }

    #[test]
    fn parse_empty_string() {
        assert_eq!(parse_experiment_class(""), BTreeSet::new());
    }

    #[test]
    fn parse_standard_keyword() {
        assert_eq!(parse_experiment_class("standard"), BTreeSet::new());
        assert_eq!(parse_experiment_class("Standard"), BTreeSet::new());
    }

    #[test]
    fn parse_single_tag() {
        assert_eq!(parse_experiment_class("tmt"), set(&["tmt"]));
    }

    #[test]
    fn parse_compound_tag() {
        assert_eq!(
            parse_experiment_class("phospho+tmt"),
            set(&["phospho", "tmt"])
        );
    }

    #[test]
    fn parse_trims_whitespace() {
        assert_eq!(
            parse_experiment_class(" phospho + tmt "),
            set(&["phospho", "tmt"])
        );
    }

    #[test]
    fn parse_lowercases() {
        assert_eq!(
            parse_experiment_class("TMT+Phospho"),
            set(&["tmt", "phospho"])
        );
    }
}
