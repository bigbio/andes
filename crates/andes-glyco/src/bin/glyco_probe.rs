//! Phase-1 / SP-A hybrid glyco gate probe.
//!
//! Measures the "searchable-backbone rate": on real intact-N-glycopeptide MS2
//! spectra with the reference engine labile-NGLYCAN ground truth, does the hybrid backbone
//! candidate generator (DB-branch ∪ de-novo Y-ladder) land the true peptide-backbone
//! mass inside a ±20 ppm candidate window for ≥90% of spectra?
//!
//! Usage: glyco_probe <mzML> <truth.tsv>
//!
//! truth.tsv columns: scan, backbone_mass, precursor_mz, precursor_z, [..]
//! `backbone_mass` is the bare peptide backbone monoisotopic neutral mass
//! WITHOUT the glycan (the reference engine `calc_neutral_pep_mass`).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use andes_glyco::glycan_db::n_glycan_list;
use andes_glyco::glycan_mass::{CORE_Y_STEPS, PROTON};
use andes_glyco::hybrid::{db_branch, hybrid_candidates_with_isotope, BackboneHit, Source};
use andes_glyco::oxonium::oxonium_gate;
// Keep solve_backbone import for the de-novo-only comparison baseline.
use andes_glyco::backbone::solve_backbone;

/// Water mass (monoisotopic). The truth `backbone_mass` and `hybrid_candidates`
/// both use the RESIDUE convention (water NOT included), exactly like the
/// production driver (`glyco_search.rs` subtracts H2O from `precursor_neutral`).
/// Y-ion ladders, however, live at the NEUTRAL peptide mass (residue + H2O), so
/// rung checks add H2O back and the de-novo baseline (which returns NEUTRAL from
/// `solve_backbone`) subtracts it.
const H2O: f64 = 18.010565;

/// Neutron mass step for isotope-error handling. The production driver sweeps
/// `isotope_error_range` (default −1..=2); the probe mirrors that so its
/// find-rate/ceiling reflect the shipped search path (Codex re-review #4).
const ISOTOPE: f64 = 1.003355;
const ISO_MIN: i8 = -1;
const ISO_MAX: i8 = 2;

/// Symmetric ±20 ppm candidate-window check (floor 0.01 Da), per the gate spec.
fn in_window(solved: f64, truth: f64) -> bool {
    (solved - truth).abs() <= (truth * 20e-6).max(0.01)
}

struct Truth {
    backbone_mass: f64,
    precursor_mz: f64,
    precursor_z: u8,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: {} <mzML> <truth.tsv>", args[0]);
        std::process::exit(2);
    }
    let mzml_path = &args[1];
    let truth_path = &args[2];

    // --- build glycan DB once ---
    let glycans = n_glycan_list();
    eprintln!("glycan DB: {} compositions", glycans.len());

    // --- load truth ---
    let mut truth: HashMap<i32, Truth> = HashMap::new();
    {
        let f = File::open(truth_path).expect("open truth.tsv");
        for (i, line) in BufReader::new(f).lines().enumerate() {
            let line = line.unwrap();
            if i == 0 || line.trim().is_empty() {
                continue; // header
            }
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 4 {
                continue;
            }
            // Tolerant parse: precursor_z may be float-formatted ("2.0"); parse as
            // f64 and round. Skip malformed rows instead of panicking (CodeRabbit).
            let (Ok(scan), Ok(backbone_mass), Ok(precursor_mz), Ok(precursor_z)) = (
                c[0].parse::<i32>(),
                c[1].parse::<f64>(),
                c[2].parse::<f64>(),
                c[3].parse::<f64>().map(|z| z.round() as u8),
            ) else {
                continue;
            };
            truth.insert(scan, Truth { backbone_mass, precursor_mz, precursor_z });
        }
    }
    eprintln!("loaded {} truth scans", truth.len());

    // --- load MS2 spectra (only need MS2) ---
    let f = File::open(mzml_path).expect("open mzML");
    let (spectra, _ms1) = input::mzml::MzMLReader::new(BufReader::new(f))
        .with_ms_level_range(2, 2)
        .read_with_ms1()
        .expect("parse mzML");
    eprintln!("read {} MS2 spectra", spectra.len());

    // index spectra by scan number
    let mut by_scan: HashMap<i32, &model::spectrum::Spectrum> = HashMap::new();
    for s in &spectra {
        if let Some(sc) = s.scan {
            by_scan.insert(sc, s);
        }
    }

    let mut n_truth = 0usize;
    let mut n_found_spec = 0usize;
    let mut n_oxonium = 0usize;

    // De-novo-only baseline (for comparison with prior Phase-1 result).
    let mut n_denovo_searchable = 0usize;

    // DB-branch-only ceiling: precursor−glycan over the whole list, independent
    // of any core-Y ladder. This is the max find-rate the DB branch could add
    // (limited only by glycan-list coverage) — the sparse stratum the Y-first
    // cascade cannot anchor lives here.
    let mut n_db_only_searchable = 0usize;

    // Hybrid counters.
    let mut n_hybrid_searchable = 0usize;
    let mut n_hybrid_from_db = 0usize;       // first hit came from DB branch
    let mut n_hybrid_from_denovo = 0usize;   // first hit from de-novo only (DB miss)

    // Sparse stratum (≤2 true core-Y rungs present): hybrid tracking.
    let mut n_sparse = 0usize;
    let mut n_sparse_hybrid_searchable = 0usize;

    for (&scan, t) in &truth {
        n_truth += 1;
        let spec = match by_scan.get(&scan) {
            Some(s) => *s,
            None => continue,
        };
        n_found_spec += 1;

        let peaks = &spec.peaks;

        // oxonium gate (informational; hybrid runs DB branch even if it misses)
        let ox = oxonium_gate(peaks, 0.10, 20.0);
        if ox.fired {
            n_oxonium += 1;
        }

        let prec_z = t.precursor_z.max(1);
        // RESIDUE convention (H2O subtracted), matching the production driver
        // (glyco_search.rs) so `hybrid_candidates` returns backbone masses
        // directly comparable to the truth `backbone_mass` (also residue).
        let precursor_neutral = (t.precursor_mz - PROTON) * prec_z as f64 - H2O;

        // --- De-novo-only baseline (top_k=5, same as Phase-1) ---
        // Run over ALL matched truth scans (same denominator as hybrid) so the
        // de-novo baseline is directly comparable to the hybrid number.
        {
            let dn_cands = solve_backbone(peaks, precursor_neutral, prec_z, 20.0, 5);
            // solve_backbone returns NEUTRAL peptide mass (Y0-derived); truth is
            // residue → subtract H2O before comparing.
            if dn_cands.iter().any(|c| in_window(c.backbone_mass - H2O, t.backbone_mass)) {
                n_denovo_searchable += 1;
            }
        }

        // Isotope sweep mirroring the driver: try each offset's corrected
        // precursor neutral and union the candidates.
        let iso_neutrals: Vec<(i8, f64)> = (ISO_MIN..=ISO_MAX)
            .map(|iso| (iso, precursor_neutral - (iso as f64) * ISOTOPE))
            .filter(|&(_, pn)| pn > 0.0)
            .collect();

        // --- DB-branch-only ceiling (precursor − each known glycan), swept ---
        {
            let searchable = iso_neutrals.iter().any(|&(iso, pn)| {
                db_branch(pn, &glycans, 500.0, prec_z, iso, None)
                    .iter()
                    .any(|h| in_window(h.backbone_mass, t.backbone_mass))
            });
            if searchable {
                n_db_only_searchable += 1;
            }
        }

        // --- Hybrid candidates (isotope-swept, top_k=50 as in the driver) ---
        let mut hybrid: Vec<BackboneHit> = Vec::new();
        for &(iso, pn) in &iso_neutrals {
            hybrid.extend(hybrid_candidates_with_isotope(peaks, pn, prec_z, iso, &glycans, 20.0, 50));
        }

        // Attribute the hit to the CLOSEST candidate within the gate window (not
        // the first in mass-sorted order) so DB-vs-de-novo split reflects which
        // branch actually resolved the truth.
        let searchable_hit = hybrid
            .iter()
            .filter(|h| in_window(h.backbone_mass, t.backbone_mass))
            .min_by(|a, b| {
                let da = (a.backbone_mass - t.backbone_mass).abs();
                let db = (b.backbone_mass - t.backbone_mass).abs();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            });
        let hybrid_searchable = searchable_hit.is_some();
        if hybrid_searchable {
            n_hybrid_searchable += 1;
            match searchable_hit.unwrap().source {
                Source::Db | Source::Transferred => n_hybrid_from_db += 1,
                Source::DeNovo => n_hybrid_from_denovo += 1,
            }
        }

        // Sparse-stratum classification (identical to Phase-1 harness).
        let rung_present = |target_neutral: f64| -> bool {
            for &(p, _) in peaks {
                for zc in 1..=prec_z {
                    let pn = (p - PROTON) * zc as f64;
                    if (pn - target_neutral).abs() <= (target_neutral * 20e-6).max(0.01) {
                        return true;
                    }
                }
            }
            false
        };
        // Y-ion ladder lives at the NEUTRAL peptide mass (residue + H2O); truth
        // is residue, so add H2O back to locate the real Y-ions.
        let neutral_bb = t.backbone_mass + H2O;
        let mut rungs = if rung_present(neutral_bb) { 1 } else { 0 };
        for &s in CORE_Y_STEPS.iter() {
            if rung_present(neutral_bb + s) {
                rungs += 1;
            }
        }
        let is_sparse = rungs <= 2;
        if is_sparse {
            n_sparse += 1;
            if hybrid_searchable {
                n_sparse_hybrid_searchable += 1;
            }
        }
    }

    let pct = |a: usize, b: usize| if b == 0 { 0.0 } else { 100.0 * a as f64 / b as f64 };

    println!("=== SP-A Hybrid Glyco Backbone Gate ===");
    println!("glycan DB compositions:      {}", glycans.len());
    println!("truth scans (in tsv):        {}", n_truth);
    println!("matched in mzML:             {}", n_found_spec);
    println!(
        "oxonium-fired:               {} ({:.1}% of matched)",
        n_oxonium,
        pct(n_oxonium, n_found_spec)
    );
    println!();
    println!(
        "--- DB-branch-only ceiling (precursor − glycan, coverage-limited) ---"
    );
    println!(
        "searchable-backbone OVERALL: {} ({:.1}% of matched)",
        n_db_only_searchable,
        pct(n_db_only_searchable, n_found_spec)
    );
    println!();
    println!("--- De-novo-only baseline (solve_backbone top_k=5, all matched scans) ---");
    println!(
        "searchable-backbone OVERALL: {} ({:.1}% of matched)",
        n_denovo_searchable,
        pct(n_denovo_searchable, n_found_spec)
    );
    println!();
    println!("--- Hybrid (DB-branch ∪ de-novo, SP-A) ---");
    println!(
        "searchable-backbone OVERALL: {} ({:.1}% of matched)",
        n_hybrid_searchable,
        pct(n_hybrid_searchable, n_found_spec)
    );
    println!(
        "  source: DB branch:         {} ({:.1}% of searchable)",
        n_hybrid_from_db,
        pct(n_hybrid_from_db, n_hybrid_searchable)
    );
    println!(
        "  source: de-novo only:      {} ({:.1}% of searchable)",
        n_hybrid_from_denovo,
        pct(n_hybrid_from_denovo, n_hybrid_searchable)
    );
    println!(
        "sparse stratum (≤2 true core-Y rungs): n={}  searchable={} ({:.1}%)",
        n_sparse,
        n_sparse_hybrid_searchable,
        pct(n_sparse_hybrid_searchable, n_sparse)
    );

    // machine-readable summary line
    println!(
        "RESULT denovo_pct={:.1} hybrid_pct={:.1} db_pct={:.1} denovo_only_pct={:.1} sparse_hybrid_pct={:.1} n_matched={} n_db_hits={} n_denovo_hits={}",
        pct(n_denovo_searchable, n_found_spec),
        pct(n_hybrid_searchable, n_found_spec),
        pct(n_hybrid_from_db, n_hybrid_searchable),
        pct(n_hybrid_from_denovo, n_hybrid_searchable),
        pct(n_sparse_hybrid_searchable, n_sparse),
        n_found_spec,
        n_hybrid_from_db,
        n_hybrid_from_denovo,
    );

    // Exit with non-zero status if gate fails, so CI can catch it.
    let hybrid_pct = pct(n_hybrid_searchable, n_found_spec);
    if hybrid_pct < 90.0 {
        eprintln!(
            "GATE FAIL: hybrid searchable-backbone {:.1}% < 90% target",
            hybrid_pct
        );
        std::process::exit(1);
    } else {
        eprintln!(
            "GATE PASS: hybrid searchable-backbone {:.1}% >= 90% target",
            hybrid_pct
        );
    }
}
