//! Selection key and backoff ladder for picking a scoring model from the manifest.
//!
//! The backoff ladder operates over a fixed `(activation, enzyme)` pair.
//! Ordering (first hit wins):
//!
//! Phase A — exact instrument (`key.instrument`):
//!   1. Exact `experiment_class` set match.
//!   2. The largest `experiment_class` subset of `key`'s set that any entry possesses.
//!      Tie-break: prefer the labeling tag (`"tmt"` or `"itraq"`) when present in the
//!      tied subsets, otherwise choose the lexicographically smallest model_id.
//!   3. The labeling tag alone if present in `key` (an entry whose set == `{"tmt"}` or
//!      `{"itraq"}` when the key's set contains it).
//!
//! Phase B — instrument family (`instrument_family(key.instrument)`), same steps 1–3.
//!
//! Phase C — empty/standard experiment_class fallback:
//!   4a. Exact instrument with empty `experiment_class`.
//!   4b. Instrument family with empty `experiment_class`.
//!
//! Last resort: `generic_id`.
//!
//! The split between phases A/B and phase C ensures that a more-specific experiment
//! match on a related instrument beats a less-specific empty-class match on the
//! exact instrument (e.g. `QExactive+tmt` wins over `OrbitrapAstral+standard` when
//! the query is `OrbitrapAstral+tmt`).

use std::collections::BTreeSet;

/// A row from the manifest representing one stored model.
#[derive(Debug, Clone, PartialEq)]
pub struct SelectionEntry {
    pub model_id: String,
    /// e.g. `"HCD"`, `"CID"`, `"ETD"`, `"UVPD"`
    pub activation: String,
    /// e.g. `"QExactive"`, `"OrbitrapAstral"`, `"TimsTOF"`, `"TOF"`, `"LowRes"`, `"HighRes"`
    pub instrument: String,
    /// e.g. `"Tryp"`
    pub enzyme: String,
    /// Parsed from the manifest's `"phospho+tmt"` string.
    pub experiment_class: BTreeSet<String>,
}

/// The query parameters for model selection.
#[derive(Debug, Clone)]
pub struct SelectionKey {
    pub activation: String,
    pub instrument: String,
    pub enzyme: String,
    pub experiment_class: BTreeSet<String>,
}

/// Labeling tags that receive priority in subset tie-breaking.
const LABELING_TAGS: &[&str] = &["tmt", "itraq"];

fn is_labeling_tag(s: &str) -> bool {
    LABELING_TAGS.contains(&s)
}

/// Returns `true` if `subset` contains a labeling tag.
fn contains_labeling_tag(set: &BTreeSet<String>) -> bool {
    set.iter().any(|s| is_labeling_tag(s.as_str()))
}

/// Try steps 1–3 (non-empty experiment_class matching) for a fixed
/// `(activation, instrument, enzyme)` triple.
///
/// Does NOT try the empty/standard fallback (step 4); that is handled separately
/// in `select` so that the family's specific match can beat the exact instrument's
/// empty fallback.
///
/// Returns `Some(model_id)` on the first step that finds a match, `None` otherwise.
fn try_steps_1_to_3<'a>(
    entries: &'a [SelectionEntry],
    activation: &str,
    instrument: &str,
    enzyme: &str,
    key_class: &BTreeSet<String>,
) -> Option<&'a str> {
    // Filter to candidates that match activation, instrument, and enzyme.
    let candidates: Vec<&SelectionEntry> = entries
        .iter()
        .filter(|e| {
            e.activation == activation && e.instrument == instrument && e.enzyme == enzyme
        })
        .collect();

    if candidates.is_empty() {
        return None;
    }

    // Step 1: exact experiment_class set match.
    if let Some(hit) = candidates
        .iter()
        .find(|e| &e.experiment_class == key_class)
    {
        return Some(&hit.model_id);
    }

    // Step 2: largest experiment_class subset of key's set.
    // A candidate entry qualifies if its experiment_class is a proper (or equal-smaller)
    // subset of key_class AND is not empty (empty is handled separately in step 4).
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
            // Find the maximum subset size.
            let max_size = subset_candidates
                .iter()
                .map(|e| e.experiment_class.len())
                .max()
                .unwrap_or(0);

            let largest: Vec<&&&SelectionEntry> = subset_candidates
                .iter()
                .filter(|e| e.experiment_class.len() == max_size)
                .collect();

            // Tie-break:
            // 1. Prefer entries whose experiment_class contains a labeling tag (tmt/itraq).
            // 2. Otherwise, prefer lexicographically smallest model_id.
            let winner = largest
                .iter()
                .max_by_key(|e| {
                    let labeling_priority = if contains_labeling_tag(&e.experiment_class) {
                        1u8
                    } else {
                        0u8
                    };
                    // To get lex-smallest with max_by_key, wrap in Reverse.
                    (labeling_priority, std::cmp::Reverse(e.model_id.clone()))
                })
                .unwrap();

            return Some(&winner.model_id);
        }
    }

    // Step 3: labeling tag alone — entry whose set == {"tmt"} or {"itraq"} when key contains it.
    for tag in LABELING_TAGS {
        if key_class.contains(*tag) {
            let tag_set: BTreeSet<String> = std::iter::once((*tag).to_string()).collect();
            if let Some(hit) = candidates.iter().find(|e| e.experiment_class == tag_set) {
                return Some(&hit.model_id);
            }
        }
    }

    None
}

/// Try step 4 (empty/standard fallback) for a fixed `(activation, instrument, enzyme)` triple.
fn try_step_4<'a>(
    entries: &'a [SelectionEntry],
    activation: &str,
    instrument: &str,
    enzyme: &str,
) -> Option<&'a str> {
    entries
        .iter()
        .find(|e| {
            e.activation == activation
                && e.instrument == instrument
                && e.enzyme == enzyme
                && e.experiment_class.is_empty()
        })
        .map(|e| e.model_id.as_str())
}

/// Pick the best `model_id` for `key` from `entries`, applying the backoff ladder.
///
/// `instrument_family` maps a specific instrument to its fallback family
/// (e.g. `"OrbitrapAstral"` → `"QExactive"`, `"TimsTOF"` → `"TOF"`); identity otherwise.
///
/// `generic_id` is the last-resort model (the historical default) returned when no
/// ladder step matches even after the family fallback.
pub fn select<'a>(
    entries: &'a [SelectionEntry],
    key: &SelectionKey,
    instrument_family: impl Fn(&str) -> String,
    generic_id: Option<&'a str>,
) -> Option<&'a str> {
    let family = instrument_family(&key.instrument);
    let has_family_fallback = family != key.instrument;

    // Phase A: steps 1–3 with exact instrument.
    if let Some(id) = try_steps_1_to_3(
        entries,
        &key.activation,
        &key.instrument,
        &key.enzyme,
        &key.experiment_class,
    ) {
        return Some(id);
    }

    // Phase B: steps 1–3 with instrument family.
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

    // Phase C4a: step 4 (empty class) with exact instrument.
    if let Some(id) = try_step_4(entries, &key.activation, &key.instrument, &key.enzyme) {
        return Some(id);
    }

    // Phase C4b: step 4 (empty class) with instrument family.
    if has_family_fallback {
        if let Some(id) = try_step_4(entries, &key.activation, &family, &key.enzyme) {
            return Some(id);
        }
    }

    // Last resort.
    generic_id
}

/// Parse a manifest experiment-class string (e.g. `"phospho+tmt"`) into a `BTreeSet<String>`.
///
/// An empty string or the literal `"standard"` maps to an empty set.
/// Otherwise the string is split on `'+'` and each token is lowercased and trimmed.
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
            e("astral", "HCD", "OrbitrapAstral", "Tryp", &[]),
            e("qe_tmt", "HCD", "QExactive", "Tryp", &["tmt"]),
            e("qe", "HCD", "QExactive", "Tryp", &[]),
            e("generic", "HCD", "QExactive", "Tryp", &[]),
        ]
    }

    #[test]
    fn exact_instrument_hit() {
        // OrbitrapAstral + empty class -> "astral" (step 1 exact match on Astral instrument)
        assert_eq!(
            select(
                &manifest(),
                &key("HCD", "OrbitrapAstral", "Tryp", &[]),
                fam,
                Some("generic")
            ),
            Some("astral")
        );
    }

    #[test]
    fn backs_off_instrument_family_then_experiment() {
        // No Astral+tmt entry -> family QExactive, exact tmt -> qe_tmt
        assert_eq!(
            select(
                &manifest(),
                &key("HCD", "OrbitrapAstral", "Tryp", &["tmt"]),
                fam,
                Some("generic")
            ),
            Some("qe_tmt")
        );
    }

    #[test]
    fn falls_back_to_generic() {
        // ETD / LowRes / Tryp has no match at all -> generic_id
        assert_eq!(
            select(
                &manifest(),
                &key("ETD", "LowRes", "Tryp", &[]),
                fam,
                Some("generic")
            ),
            Some("generic")
        );
    }

    #[test]
    fn largest_subset_when_no_exact_combo() {
        let mut m = manifest();
        m.push(e("qe_phos", "HCD", "QExactive", "Tryp", &["phospho"]));
        // key {phospho, tmt}: no combo model -> largest subset.
        // Both {phospho} (qe_phos) and {tmt} (qe_tmt) are size-1 subsets.
        // Tie-break: prefer the labeling tag (tmt) -> qe_tmt.
        assert_eq!(
            select(
                &m,
                &key("HCD", "QExactive", "Tryp", &["phospho", "tmt"]),
                fam,
                Some("generic")
            ),
            Some("qe_tmt")
        );
    }

    // ---- Additional coverage ----

    #[test]
    fn no_generic_id_returns_none() {
        assert_eq!(
            select(
                &manifest(),
                &key("ETD", "LowRes", "Tryp", &[]),
                fam,
                None
            ),
            None
        );
    }

    #[test]
    fn family_fallback_empty_class() {
        // TimsTOF has no entry; family TOF has none either; falls through to generic.
        assert_eq!(
            select(
                &manifest(),
                &key("HCD", "TimsTOF", "Tryp", &[]),
                fam,
                Some("generic")
            ),
            Some("generic")
        );
    }

    #[test]
    fn step3_labeling_tag_alone() {
        // key has {tmt,phospho}, no QExactive combo entry, no QExactive phospho+tmt entry,
        // and qe_phos is not in manifest -> step 3 finds qe_tmt (set == {"tmt"}).
        // (This mirrors the backs_off_instrument_family_then_experiment test but with
        //  an extra phospho class that has no dedicated model.)
        let m = vec![
            e("qe_tmt", "HCD", "QExactive", "Tryp", &["tmt"]),
            e("qe", "HCD", "QExactive", "Tryp", &[]),
        ];
        assert_eq!(
            select(
                &m,
                &key("HCD", "QExactive", "Tryp", &["phospho", "tmt"]),
                fam,
                Some("generic")
            ),
            Some("qe_tmt")
        );
    }

    #[test]
    fn step4_empty_class_fallback() {
        // key has {phospho} but no phospho entry; step 4 -> empty-class entry.
        let m = vec![e("qe", "HCD", "QExactive", "Tryp", &[])];
        assert_eq!(
            select(
                &m,
                &key("HCD", "QExactive", "Tryp", &["phospho"]),
                fam,
                Some("generic")
            ),
            Some("qe")
        );
    }

    // ---- parse_experiment_class ----

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
        assert_eq!(parse_experiment_class("TMT+Phospho"), set(&["tmt", "phospho"]));
    }
}
