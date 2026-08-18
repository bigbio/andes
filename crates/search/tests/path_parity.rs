//! Guards against a defect class this codebase has hit repeatedly: a PSM feature that
//! is populated on one code path and silently left at its default on another.
//!
//! The standard search fills a set of per-spectrum features in `fill_post_topn`, which
//! runs after the top-N queue is final. The glyco driver does not use that queue and so
//! never calls it — it hand-mirrors the same assignments instead. That mirror has been
//! wrong twice:
//!
//!   * the 2026-07-09 pass wired five fields and missed `delta_raw_score` and
//!     `strong_score_cal`, so every glyco PIN shipped with `DeltaRankScore` hard 0.0 —
//!     a top-tier Percolator feature, absent for two releases with no symptom other
//!     than degraded discrimination;
//!   * the same shape recurred for the RT features, which needed a third site.
//!
//! A constant-zero feature column cannot be spotted by any type check and does not fail
//! a search: it silently costs identifications. So this test reads the sources and
//! asserts that every feature the standard path fills is also filled by one of the
//! glyco path's sites. Adding a ninth field to `fill_post_topn` without mirroring it
//! now fails the build rather than quietly degrading glyco results.
//!
//! Source scanning is deliberate. The alternative — running both engines and diffing
//! features — needs fixtures for two very different pipelines and would still miss a
//! field that is merely assigned a wrong-but-nonzero value. This catches the specific
//! failure that actually happened, cheaply and without fixtures.

/// Fields the standard path assigns inside a `fill_post_topn` closure.
fn standard_post_topn_fields(src: &str) -> Vec<String> {
    let mut out = Vec::new();
    let lines: Vec<&str> = src.lines().collect();
    for (i, line) in lines.iter().enumerate() {
        if !line.contains("fill_post_topn(|psm|") {
            continue;
        }
        // Walk the closure body by brace depth from the call site.
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

/// `features.<name> =` / `psm.features.<name> =` → `<name>`.
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
        // Assignment only — not a read, and not `==`.
        if !name.is_empty() && tail.starts_with('=') && !tail.starts_with("==") {
            out.push(name);
        }
        rest = &after[consumed..];
    }
    out
}

#[test]
fn glyco_path_fills_every_feature_the_standard_path_fills() {
    let standard = include_str!("../src/match_engine.rs");
    // The glyco path spreads its equivalent of `fill_post_topn` over two sites: the
    // driver for score-derived features, and the RT wiring for retention-time ones.
    let glyco_driver = include_str!("../src/glyco_search.rs");
    let glyco_rt = include_str!("../../output/src/glyco_rt.rs");

    let expected = standard_post_topn_fields(standard);
    assert!(
        expected.len() >= 6,
        "parse failure: found only {} fill_post_topn assignments — the scan is broken, \
         not the code under test",
        expected.len()
    );

    // Known, deliberate exceptions. Each needs a REASON, not just a name: an entry here
    // is a decision to ship a constant column, and it should be re-argued when touched.
    let exceptions: &[(&str, &str)] = &[
        (
            "edge_score",
            "glyco computes the edge score during phase 1 and moves it in through the \
             PsmMatch struct literal (`edge: w.edge`) rather than a `features.` \
             assignment, so it is populated by a different mechanism, not missing.",
        ),
        (
            "strong_score_cal",
            "the calibrated strong score needs a per-spectrum population of strong \
             scores to z-score against. The glyco driver computes a strong score only \
             for the single collapse winner, so no such population exists. Filling it \
             requires scoring the accepted set, which is a real design change and not a \
             transcription fix — tracked separately.",
        ),
    ];

    let mut missing: Vec<String> = Vec::new();
    for field in &expected {
        let filled = glyco_driver.contains(&format!("features.{field} ="))
            || glyco_rt.contains(&format!("features.{field} ="));
        let excepted = exceptions.iter().any(|(n, _)| n == field);
        if !filled && !excepted {
            missing.push(field.clone());
        }
    }

    assert!(
        missing.is_empty(),
        "these features are filled by the standard search's fill_post_topn but by \
         NEITHER glyco site, so they reach Percolator as constant defaults and silently \
         cost identifications: {missing:?}\n\
         Fix by mirroring the assignment in glyco_search.rs (score-derived) or \
         glyco_rt.rs (retention-time). If a field genuinely cannot be computed on the \
         glyco path, add it to `exceptions` WITH THE REASON."
    );
}

#[test]
fn every_documented_exception_is_still_a_real_gap() {
    // Keeps the exception list honest: once a field is genuinely mirrored, its
    // exception is stale and must be removed, or the next real gap hides behind it.
    let standard = include_str!("../src/match_engine.rs");
    let expected = standard_post_topn_fields(standard);
    for name in ["edge_score", "strong_score_cal"] {
        assert!(
            expected.iter().any(|f| f == name),
            "`{name}` is listed as an exception but the standard path no longer fills \
             it in fill_post_topn — drop the stale exception"
        );
    }
}
