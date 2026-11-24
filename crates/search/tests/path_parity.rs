fn standard_post_topn_fields(src: &str) -> Vec<String> {
    let mut out = Vec::new();
    let lines: Vec<&str> = src.lines().collect();
    for (i, line) in lines.iter().enumerate() {
        if !line.contains("fill_post_topn(|psm|") && !line.contains("let make_features = |") {
            continue;
        }
        let mut depth: i32 = 0;
        let mut opened = false;
        for l in &lines[i..] {
            depth += l.matches('{').count() as i32;
            depth -= l.matches('}').count() as i32;
            if l.contains('{') {
                opened = true;
            }
            for f in assigned_features(l) {
                out.push(f);
            }
            if opened && depth <= 0 {
                break;
            }
        }
    }
    out.sort();
    out.dedup();
    out
}

fn assigned_features(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = line;
    while let Some(pos) = rest.find("features.") {
        let after = &rest[pos + "features.".len()..];
        let name: String = after
            .chars()
            .take_while(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || *c == '_')
            .collect();
        let consumed = name.len().max(1);
        let tail = after[name.len()..].trim_start();
        if !name.is_empty() && tail.starts_with('=') && !tail.starts_with("==") {
            out.push(name);
        }
        rest = &after[consumed..];
    }
    out
}

const EXCEPTIONS: &[(&str, &str)] = &[
    (
        "strong_score_cal",
        "the calibrated strong score needs a per-spectrum population of strong \
         scores to z-score against. The glyco driver computes a strong score only \
         for the single collapse winner, so no such population exists. Filling it \
         requires scoring the accepted set, which is a real design change and not a \
         transcription fix — tracked separately.",
    ),
    (
        "ref_competitor",
        "the gap to the best Reference-entry candidate only matters when the \
         database holds UniProt isoform, variant or chain entries. The glyco driver \
         scores backbones per glycan composition and keeps no per-spectrum \
         Reference tracker; a glyco search over such entries reports 0 here \
         (documented in DOCS §1a).",
    ),
    ("ref_best_raw", "companion of ref_competitor, same reason."),
    (
        "open_site",
        "`--open-search` is refused with `--glyco` (search.rs), so the glyco driver \
         never sees an open-search PSM. The three Open* fields stay 0 there by design.",
    ),
    ("open_loc_score", "companion of open_site, same reason."),
    ("open_loc_margin", "companion of open_site, same reason."),
    ("open_delta", "companion of open_site, same reason."),
    ("open_used_shift", "companion of open_site, same reason."),
];

#[test]
fn glyco_path_fills_every_feature_the_standard_path_fills() {
    let standard = include_str!("../src/engine/run.rs");
    let glyco_driver = include_str!("../src/glycopeptides.rs");
    let glyco_rt = include_str!("../../output/src/glycopeptides/rt.rs");

    let expected = standard_post_topn_fields(standard);
    assert!(
        expected.len() >= 6,
        "parse failure: found only {} fill_post_topn assignments — the scan is broken, \
         not the code under test",
        expected.len()
    );

    let mut missing: Vec<String> = Vec::new();
    for field in &expected {
        let filled = glyco_driver.contains(&format!("features.{field} ="))
            || glyco_rt.contains(&format!("features.{field} ="));
        let excepted = EXCEPTIONS.iter().any(|(n, _)| n == field);
        if !filled && !excepted {
            missing.push(field.clone());
        }
    }

    assert!(
        missing.is_empty(),
        "these features are filled by the standard search's fill_post_topn but by \
         NEITHER glyco site, so they reach Percolator as constant defaults and silently \
         cost identifications: {missing:?}\n\
         Fix by mirroring the assignment in glycopeptides.rs (score-derived) or \
         glycopeptides/rt.rs (retention-time). If a field genuinely cannot be computed on the \
         glyco path, add it to `exceptions` WITH THE REASON."
    );
}

#[test]
fn every_documented_exception_is_still_a_real_gap() {
    let standard = include_str!("../src/engine/run.rs");
    let expected = standard_post_topn_fields(standard);
    for (name, _) in EXCEPTIONS {
        assert!(
            expected.iter().any(|f| f == name),
            "`{name}` is listed as an exception but the standard path no longer fills \
             it in fill_post_topn — drop the stale exception"
        );
    }
}
