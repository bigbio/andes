use crate::backbone::{
    core_y_intensity, solve_backbone_min, BackboneCandidate, SpectrumStats, H2O,
};
use crate::glycan::db::GlycanComp;
use crate::glycan::first::GlycanCore;
use crate::oxonium::oxonium_gate;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Source {
    Db,
    DeNovo,
}

#[derive(Debug, Clone)]
pub struct BackboneHit {
    pub backbone_mass: f64,
    pub glycan: Option<GlycanComp>,
    pub source: Source,
    pub charge: u8,
    pub isotope_offset: i8,
    pub glycan_mass_residual: f64,
    pub peptide_first: bool,
}

pub fn db_branch(
    precursor_neutral: f64,
    glycans: &[GlycanComp],
    min_backbone: f64,
    charge: u8,
    isotope_offset: i8,
) -> Vec<BackboneHit> {
    let mut out: Vec<BackboneHit> = glycans
        .iter()
        .filter_map(|g| {
            let bb = precursor_neutral - g.mass;
            if bb >= min_backbone {
                Some(BackboneHit {
                    backbone_mass: bb,
                    glycan: Some(g.clone()),
                    source: Source::Db,
                    charge,
                    isotope_offset,
                    glycan_mass_residual: precursor_neutral - bb,
                    peptide_first: false,
                })
            } else {
                None
            }
        })
        .collect();

    out.sort_by(|a, b| {
        a.backbone_mass
            .to_bits()
            .cmp(&b.backbone_mass.to_bits())
            .then_with(|| {
                let ga = a.glycan.as_ref();
                let gb = b.glycan.as_ref();
                match (ga, gb) {
                    (Some(x), Some(y)) => x
                        .hexnac
                        .cmp(&y.hexnac)
                        .then(x.hex.cmp(&y.hex))
                        .then(x.fuc.cmp(&y.fuc))
                        .then(x.neuac.cmp(&y.neuac))
                        .then(x.neugc.cmp(&y.neugc))
                        .then(x.topology_id.cmp(&y.topology_id)),
                    _ => std::cmp::Ordering::Equal,
                }
            })
    });

    out
}

fn nearest_glycan(glycans: &[GlycanComp], target_mass: f64, tol_ppm: f64) -> Option<GlycanComp> {
    if target_mass <= 0.0 {
        return None;
    }
    let tol = (target_mass * tol_ppm * 1e-6_f64).max(0.02);
    let mut best: Option<(f64, &GlycanComp)> = None;
    for g in glycans {
        let d = (g.mass - target_mass).abs();
        if d <= tol && best.is_none_or(|(bd, _)| d < bd) {
            best = Some((d, g));
        }
    }
    best.map(|(_, g)| g.clone())
}

pub fn hybrid_candidates(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    glycans: &[GlycanComp],
    tol_ppm: f64,
    top_k: usize,
) -> Vec<BackboneHit> {
    hybrid_candidates_with_isotope(
        peaks,
        precursor_neutral,
        precursor_z,
        0,
        glycans,
        tol_ppm,
        top_k,
    )
}

pub fn hybrid_candidates_with_isotope(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    isotope_offset: i8,
    glycans: &[GlycanComp],
    tol_ppm: f64,
    top_k: usize,
) -> Vec<BackboneHit> {
    let presolved =
        solve_backbones_for_charge(peaks, precursor_neutral, precursor_z, tol_ppm, top_k);
    hybrid_candidates_presolved(
        presolved.as_deref(),
        peaks,
        precursor_neutral,
        precursor_z,
        isotope_offset,
        glycans,
        tol_ppm,
        top_k,
        false,
    )
}

pub fn solve_backbones_for_charge(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    tol_ppm: f64,
    top_k: usize,
) -> Option<Vec<BackboneCandidate>> {
    let ox = oxonium_gate(peaks, 0.10, tol_ppm);
    if !ox.fired {
        return None;
    }

    let mut dn = solve_backbone_min(peaks, precursor_neutral, precursor_z, tol_ppm, top_k, 2);
    if dn.is_empty() {
        dn = solve_backbone_min(peaks, precursor_neutral, precursor_z, tol_ppm, top_k, 1);
    }
    Some(dn)
}

#[allow(clippy::too_many_arguments)]
pub fn hybrid_candidates_presolved(
    presolved: Option<&[BackboneCandidate]>,
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    isotope_offset: i8,
    glycans: &[GlycanComp],
    tol_ppm: f64,
    top_k: usize,
    force_db_on_none: bool,
) -> Vec<BackboneHit> {
    const MIN_BACKBONE: f64 = 500.0;
    const MIN_GLYCAN: f64 = 406.0;

    let dn = match presolved {
        None => {
            if force_db_on_none {
                return db_branch(
                    precursor_neutral,
                    glycans,
                    MIN_BACKBONE,
                    precursor_z,
                    isotope_offset,
                );
            }
            return Vec::new();
        }
        Some(v) => v,
    };

    let mut combined: Vec<BackboneHit> = Vec::with_capacity(dn.len());
    for c in dn {
        let bb = c.backbone_mass - H2O;
        if bb >= precursor_neutral {
            continue;
        }
        if precursor_neutral - bb < MIN_GLYCAN {
            continue;
        }
        if bb < MIN_BACKBONE {
            continue;
        }
        let residual = precursor_neutral - bb;
        let glycan = nearest_glycan(glycans, residual, tol_ppm);
        let source = if glycan.is_some() {
            Source::Db
        } else {
            Source::DeNovo
        };
        combined.push(BackboneHit {
            backbone_mass: bb,
            glycan,
            source,
            charge: precursor_z,
            isotope_offset,
            glycan_mass_residual: residual,
            peptide_first: false,
        });
    }

    let y_first_has_evidence = !combined.is_empty();
    if y_first_has_evidence {
        combined.extend(db_branch(
            precursor_neutral,
            glycans,
            MIN_BACKBONE,
            precursor_z,
            isotope_offset,
        ));
    } else {
        combined = db_branch(
            precursor_neutral,
            glycans,
            MIN_BACKBONE,
            precursor_z,
            isotope_offset,
        );
    }

    combined.sort_by(|a, b| {
        a.backbone_mass
            .to_bits()
            .cmp(&b.backbone_mass.to_bits())
            .then_with(|| {
                let ord_a = if a.source == Source::Db { 0u8 } else { 1u8 };
                let ord_b = if b.source == Source::Db { 0u8 } else { 1u8 };
                ord_a.cmp(&ord_b)
            })
    });

    if combined.is_empty() {
        return combined;
    }

    let mut deduped: Vec<BackboneHit> = Vec::with_capacity(combined.len());
    let mut rep = combined.remove(0);

    for next in combined {
        let tol = (rep.backbone_mass * tol_ppm * 1e-6_f64).max(0.01);
        if (next.backbone_mass - rep.backbone_mass).abs() < tol {
            if rep.source == Source::DeNovo && next.source == Source::Db {
                rep = next;
            }
        } else {
            deduped.push(rep);
            rep = next;
        }
    }
    deduped.push(rep);

    if y_first_has_evidence && deduped.len() > top_k {
        let stats = SpectrumStats::new(peaks);
        let mut scored: Vec<(f64, BackboneHit)> = deduped
            .into_iter()
            .map(|h| {
                (
                    core_y_intensity(
                        peaks,
                        &stats,
                        h.backbone_mass + H2O,
                        tol_ppm,
                        precursor_z,
                        GlycanCore::NGlycan,
                    ),
                    h,
                )
            })
            .collect();
        scored.sort_by(|a, b| {
            b.0.total_cmp(&a.0).then_with(|| {
                a.1.backbone_mass
                    .to_bits()
                    .cmp(&b.1.backbone_mass.to_bits())
            })
        });
        scored.truncate(top_k);
        deduped = scored.into_iter().map(|(_, h)| h).collect();
    }

    deduped
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan::db::load_glycan_gdb;
    use crate::glycan::mass::{HEX, HEXNAC, PROTON};

    fn test_glycans() -> Vec<GlycanComp> {
        load_glycan_gdb(
            "H,N,A,G,F\n\
             (N(N(H(H)(H))))\n\
             (N(F)(N(H(H)(H))))\n\
             (N(N(H(H)(H(H)))))\n\
             (N(N(H(H)(H(N)))))\n",
        )
        .unwrap()
    }

    #[test]
    fn db_branch_recovers_backbone_for_known_glycan() {
        let glycans = test_glycans();
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        let hits = db_branch(precursor, &glycans, 500.0, 2, 0);
        assert!(!hits.is_empty(), "expected DB branch hits");

        let found = hits
            .iter()
            .any(|h| (h.backbone_mass - true_backbone).abs() < 0.01 && h.source == Source::Db);
        assert!(
            found,
            "did not find backbone at {:.4} in DB branch",
            true_backbone
        );
    }

    #[test]
    fn db_branch_filters_below_min() {
        let glycans = test_glycans();
        let precursor = 600.0;
        let hits = db_branch(precursor, &glycans, 500.0, 2, 0);
        for h in &hits {
            assert!(
                h.backbone_mass >= 500.0,
                "backbone below min: {}",
                h.backbone_mass
            );
        }
    }

    #[test]
    fn db_branch_is_sorted() {
        let glycans = test_glycans();
        let precursor = 4000.0;
        let hits = db_branch(precursor, &glycans, 500.0, 2, 0);
        for w in hits.windows(2) {
            assert!(
                w[0].backbone_mass <= w[1].backbone_mass + 1e-9,
                "not sorted: {} > {}",
                w[0].backbone_mass,
                w[1].backbone_mass
            );
        }
    }

    #[test]
    fn hybrid_union_contains_both_sources() {
        let glycans = test_glycans();
        let proton = PROTON;
        let steps = crate::glycan::mass::CORE_Y_STEPS;

        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;

        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (y0_neutral + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);
        assert!(!hits.is_empty(), "expected hybrid hits");

        let has_db = hits.iter().any(|h| h.source == Source::Db);
        assert!(
            has_db,
            "expected at least one annotated (Db) hit, got {:?}",
            hits
        );
    }

    #[test]
    fn db_union_stays_bounded_and_keeps_true_backbone() {
        let glycans = test_glycans();
        let proton = PROTON;
        let steps = crate::glycan::mass::CORE_Y_STEPS;
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;

        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (y0_neutral + proton, 1000.0),
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 600.0));
        }

        let top_k = 5;
        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, top_k);
        assert!(
            hits.len() <= top_k,
            "DB-union must stay bounded to top_k, got {}",
            hits.len()
        );
        let found = hits
            .iter()
            .any(|h| (h.backbone_mass - true_backbone_residue).abs() < 0.05);
        assert!(
            found,
            "true backbone (real core-Y ladder) must survive core-Y-ranked truncation"
        );
    }

    #[test]
    fn presolved_once_matches_per_isotope_union() {
        const ISOTOPE: f64 = 1.00335;
        let glycans = test_glycans();
        let proton = PROTON;
        let steps = crate::glycan::mass::CORE_Y_STEPS;
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone_residue = 1500.0_f64;
        let z = 2u8;
        let y0_neutral = true_backbone_residue + H2O;

        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (y0_neutral + proton, 1000.0),
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 600.0));
        }
        peaks.push((400.0, 300.0));
        peaks.push((true_backbone_residue + 20.025118 - 400.0, 300.0));
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));

        let tol = 20.0;
        let top_k = 50;
        let (iso_min, iso_max) = (-1i8, 2i8);
        let observed_neutral = true_backbone_residue + glycan_mass;

        let key = |h: &BackboneHit| {
            let src = match h.source {
                Source::Db => 0u8,
                Source::DeNovo => 1u8,
            };
            (
                (h.backbone_mass * 100.0).round() as i64,
                h.isotope_offset,
                src,
                h.glycan
                    .as_ref()
                    .map(|g| (g.hexnac, g.hex, g.fuc, g.neuac, g.neugc)),
            )
        };
        let sorted = |hits: &[BackboneHit]| {
            let mut k: Vec<_> = hits.iter().map(key).collect();
            k.sort();
            k
        };

        let mut old: Vec<BackboneHit> = Vec::new();
        for iso in iso_min..=iso_max {
            let pn = observed_neutral - (iso as f64) * ISOTOPE;
            if pn <= 0.0 {
                continue;
            }
            old.extend(hybrid_candidates_with_isotope(
                &peaks, pn, z, iso, &glycans, tol, top_k,
            ));
        }

        let pn_max = observed_neutral - (iso_min as f64) * ISOTOPE;
        let presolved = solve_backbones_for_charge(&peaks, pn_max, z, tol, top_k);
        let mut new: Vec<BackboneHit> = Vec::new();
        for iso in iso_min..=iso_max {
            let pn = observed_neutral - (iso as f64) * ISOTOPE;
            if pn <= 0.0 {
                continue;
            }
            new.extend(hybrid_candidates_presolved(
                presolved.as_deref(),
                &peaks,
                pn,
                z,
                iso,
                &glycans,
                tol,
                top_k,
                false,
            ));
        }

        assert!(
            !old.is_empty(),
            "expected non-empty union to make the test meaningful"
        );
        assert_eq!(
            sorted(&old),
            sorted(&new),
            "factored solve-once path must equal per-isotope union"
        );
    }

    #[test]
    fn hybrid_dedup_keeps_db_over_denovo() {
        let glycans = test_glycans();
        let proton = PROTON;
        let steps = crate::glycan::mass::CORE_Y_STEPS;

        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;

        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (y0_neutral + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);

        let near: Vec<&BackboneHit> = hits
            .iter()
            .filter(|h| (h.backbone_mass - true_backbone_residue).abs() < 0.02)
            .collect();
        assert_eq!(
            near.len(),
            1,
            "expected exactly one candidate near true backbone after dedup, got {}",
            near.len()
        );
        assert_eq!(
            near[0].source,
            Source::Db,
            "expected Db source to win dedup, got {:?}",
            near[0].source
        );
    }

    #[test]
    fn denovo_and_db_branches_agree_on_residue_mass_convention() {
        let proton = PROTON;
        let steps = crate::glycan::mass::CORE_Y_STEPS;

        let true_backbone_residue = 1500.0_f64;
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let precursor_neutral = true_backbone_residue + glycan_mass;

        let neutral_backbone = true_backbone_residue + crate::backbone::H2O;
        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (neutral_backbone + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((neutral_backbone + s + proton, 90.0));
        }

        let raw = crate::backbone::solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);
        assert!(!raw.is_empty(), "expected solve_backbone candidates");
        let raw_gap = raw[0].backbone_mass - true_backbone_residue;
        assert!(
            (raw_gap - crate::backbone::H2O).abs() < 0.01,
            "expected solve_backbone's raw output to be ~H2O ({:.4}) above the \
             residue-mass backbone, got gap={:.4} (raw={:.4}, residue={:.4})",
            crate::backbone::H2O,
            raw_gap,
            raw[0].backbone_mass,
            true_backbone_residue
        );

        let hits = hybrid_candidates(&peaks, precursor_neutral, 2, &[], 20.0, 5);
        let dn_hit = hits
            .iter()
            .filter(|h| h.source == Source::DeNovo)
            .min_by(|a, b| {
                (a.backbone_mass - true_backbone_residue)
                    .abs()
                    .total_cmp(&(b.backbone_mass - true_backbone_residue).abs())
            })
            .expect("expected at least one DeNovo BackboneHit");
        assert!(
            (dn_hit.backbone_mass - true_backbone_residue).abs() < 0.01,
            "DeNovo backbone_mass must match residue-mass convention: expected \
             ~{:.4}, got {:.4} (gap={:.4}); note the pre-fix bug would leave this \
             ~H2O ({:.4}) away instead",
            true_backbone_residue,
            dn_hit.backbone_mass,
            dn_hit.backbone_mass - true_backbone_residue,
            crate::backbone::H2O
        );
    }

    #[test]
    fn zero_core_y_known_glycan_is_recovered_by_db_fallback() {
        let glycans = test_glycans();
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;

        let peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (366.13950, 60.0),
            (900.0, 20.0),
        ];
        assert!(
            crate::backbone::solve_backbone_min(&peaks, precursor, 2, 20.0, 50, 1).is_empty(),
            "premise: even quorum-1 solver must be empty for a 0-core-Y spectrum"
        );

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 50);
        let found = hits.iter().any(|h| {
            (h.backbone_mass - true_backbone_residue).abs() < 0.02 && h.source == Source::Db
        });
        assert!(
            found,
            "DB fallback must recover the 0-core-Y known-glycan backbone; got {hits:?}"
        );
    }

    #[test]
    fn novel_glycan_hit_preserves_residual_mass() {
        let glycans: Vec<GlycanComp> = Vec::new();
        let proton = PROTON;
        let steps = crate::glycan::mass::CORE_Y_STEPS;

        let novel_glycan_mass = 1234.5678_f64;

        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + novel_glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;

        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (y0_neutral + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);
        let novel = hits
            .iter()
            .filter(|h| h.source == Source::DeNovo && h.glycan.is_none())
            .min_by(|a, b| {
                (a.backbone_mass - true_backbone_residue)
                    .abs()
                    .total_cmp(&(b.backbone_mass - true_backbone_residue).abs())
            })
            .expect("expected a novel (unannotated) DeNovo hit near the true backbone");

        assert!(
            (novel.glycan_mass_residual - (precursor - novel.backbone_mass)).abs() < 1e-6,
            "residual must be precursor − backbone; got {}",
            novel.glycan_mass_residual
        );
        assert!(
            (novel.glycan_mass_residual - novel_glycan_mass).abs() < 0.02,
            "residual must recover the novel glycan mass ~{novel_glycan_mass}, got {}",
            novel.glycan_mass_residual
        );
    }

    #[test]
    fn db_branch_records_charge_and_isotope_offset() {
        let glycans = test_glycans();
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        let hits = db_branch(precursor, &glycans, 500.0, 3, -1);
        assert!(!hits.is_empty());
        for h in &hits {
            assert_eq!(h.charge, 3, "charge must be threaded onto BackboneHit");
            assert_eq!(
                h.isotope_offset, -1,
                "isotope_offset must be threaded onto BackboneHit"
            );
        }
    }

    #[test]
    fn weak_ladder_spectrum_is_rescued_by_quorum1_cascade() {
        let glycans = test_glycans();
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;
        let proton = PROTON;

        let peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (y0_neutral + proton, 120.0),
        ];

        let primary = crate::backbone::solve_backbone(&peaks, precursor, 2, 20.0, 50);
        assert!(
            primary.is_empty(),
            "premise: quorum-2 solver should find nothing with a single core-Y rung, got {:?}",
            primary
        );

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 50);
        let found = hits.iter().any(|h| {
            (h.backbone_mass - true_backbone_residue).abs() < 0.05 && h.source == Source::Db
        });
        assert!(
            found,
            "quorum-1 cascade must recover the true backbone for a weak-ladder \
             known-glycan spectrum; got {:?}",
            hits
        );
    }

    #[test]
    fn hybrid_candidates_with_isotope_threads_offset_onto_all_hits() {
        let glycans = test_glycans();
        let proton = PROTON;
        let steps = crate::glycan::mass::CORE_Y_STEPS;

        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (true_backbone + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((true_backbone + s + proton, 90.0));
        }

        let hits = hybrid_candidates_with_isotope(&peaks, precursor, 2, 2, &glycans, 20.0, 10);
        assert!(!hits.is_empty(), "expected hybrid hits");
        for h in &hits {
            assert_eq!(
                h.isotope_offset, 2,
                "every hit must carry the caller's isotope offset"
            );
            assert_eq!(h.charge, 2, "every hit must carry the caller's charge");
        }
    }
}
