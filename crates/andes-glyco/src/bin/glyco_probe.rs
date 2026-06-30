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
use andes_glyco::hybrid::{hybrid_candidates, Source};
use andes_glyco::oxonium::oxonium_gate;
// Keep solve_backbone import for the de-novo-only comparison baseline.
use andes_glyco::backbone::solve_backbone;

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
            let scan: i32 = c[0].parse().unwrap();
            truth.insert(
                scan,
                Truth {
                    backbone_mass: c[1].parse().unwrap(),
                    precursor_mz: c[2].parse().unwrap(),
                    precursor_z: c[3].parse().unwrap(),
                },
            );
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
        let precursor_neutral = (t.precursor_mz - PROTON) * prec_z as f64;

        // --- De-novo-only baseline (top_k=5, same as Phase-1) ---
        // Run over ALL matched truth scans (same denominator as hybrid) so the
        // de-novo baseline is directly comparable to the hybrid number.
        {
            let dn_cands = solve_backbone(peaks, precursor_neutral, prec_z, 20.0, 5);
            if dn_cands.iter().any(|c| in_window(c.backbone_mass, t.backbone_mass)) {
                n_denovo_searchable += 1;
            }
        }

        // --- Hybrid candidates ---
        let hybrid = hybrid_candidates(peaks, precursor_neutral, prec_z, &glycans, 20.0, 5);

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
                Source::Db => n_hybrid_from_db += 1,
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
        let mut rungs = if rung_present(t.backbone_mass) { 1 } else { 0 };
        for &s in CORE_Y_STEPS.iter() {
            if rung_present(t.backbone_mass + s) {
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
    println!("--- De-novo-only baseline (Phase-1 top_k=5, oxonium-gated) ---");
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
