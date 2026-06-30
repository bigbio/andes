//! Phase-1 glyco gate probe.
//!
//! Measures the "searchable-backbone rate": on real intact-N-glycopeptide MS2
//! spectra with the reference engine labile-NGLYCAN ground truth, does the backbone-mass
//! solver land the true peptide-backbone mass inside a ±20 ppm candidate window
//! for ≥70% of spectra — and does the rate hold in the sparse (≤1-core-Y)
//! stratum?
//!
//! Usage: glyco_probe <mzML> <truth.tsv>
//!
//! truth.tsv columns: scan, backbone_mass, precursor_mz, precursor_z, [..]
//! `backbone_mass` is the bare peptide backbone monoisotopic neutral mass
//! WITHOUT the glycan (the reference engine `calc_neutral_pep_mass`).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use andes_glyco::backbone::solve_backbone;
use andes_glyco::glycan_mass::{CORE_Y_STEPS, PROTON};
use andes_glyco::oxonium::oxonium_gate;
use input::mzml::MzMLReader;

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
    let (spectra, _ms1) = MzMLReader::new(BufReader::new(f))
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
    let mut n_searchable = 0usize;

    // sparse stratum: scans whose BEST candidate has core_y_hits <= 1.
    // (A candidate must have >=2 core-Y hits to survive the solver quorum, so
    //  "sparse" here means the spectrum yields little core-Y ladder support —
    //  we classify by the searchable candidate's core_y_hits when found, else by
    //  the top candidate's core_y_hits.)
    let mut n_sparse = 0usize;
    let mut n_sparse_searchable = 0usize;

    for (&scan, t) in &truth {
        n_truth += 1;
        let spec = match by_scan.get(&scan) {
            Some(s) => *s,
            None => continue,
        };
        n_found_spec += 1;

        let peaks = &spec.peaks;

        // oxonium gate
        let ox = oxonium_gate(peaks, 0.10, 20.0);
        if !ox.fired {
            continue;
        }
        n_oxonium += 1;

        let prec_z = t.precursor_z.max(1);
        let precursor_neutral = (t.precursor_mz - PROTON) * prec_z as f64;

        let cands = solve_backbone(peaks, precursor_neutral, prec_z, 20.0, 5);

        // best searchable candidate (if any)
        let searchable_cand = cands
            .iter()
            .find(|c| in_window(c.backbone_mass, t.backbone_mass));
        let searchable = searchable_cand.is_some();
        if searchable {
            n_searchable += 1;
        }

        // Sparse-stratum classification. The solver only emits candidates with
        // >=2 core-Y hits (quorum), so "best candidate core_y_hits <= 1" never
        // occurs; the meaningful sparse axis is how much real core-Y evidence the
        // SPECTRUM carries for the TRUE backbone. We count how many of {Y0, Y1..Y5}
        // for the true backbone are actually present (within 20 ppm, any charge
        // 1..=z); "sparse" = <=2 rungs present (Y0 + at most one core step).
        let prec_z_u = prec_z;
        let rung_present = |target_neutral: f64| -> bool {
            for &(p, _) in peaks {
                for zc in 1..=prec_z_u {
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
        let is_sparse = rungs <= 2; // <=2 true core-Y rungs present
        if is_sparse {
            n_sparse += 1;
            if searchable {
                n_sparse_searchable += 1;
            }
        }
    }

    let pct = |a: usize, b: usize| if b == 0 { 0.0 } else { 100.0 * a as f64 / b as f64 };

    println!("=== Phase-1 glyco backbone-solver gate ===");
    println!("truth scans (in tsv):        {}", n_truth);
    println!("matched in mzML:             {}", n_found_spec);
    println!(
        "oxonium-fired:               {} ({:.1}% of matched)",
        n_oxonium,
        pct(n_oxonium, n_found_spec)
    );
    println!(
        "searchable-backbone OVERALL: {} ({:.1}% of matched)",
        n_searchable,
        pct(n_searchable, n_found_spec)
    );
    println!(
        "  (as % of oxonium-fired):   {:.1}%",
        pct(n_searchable, n_oxonium)
    );
    println!(
        "sparse stratum (<=2 true core-Y rungs): n={}  searchable={} ({:.1}%)",
        n_sparse,
        n_sparse_searchable,
        pct(n_sparse_searchable, n_sparse)
    );

    // machine-readable summary line
    println!(
        "RESULT oxonium_pct={:.1} searchable_overall_pct={:.1} sparse_pct={:.1} n_matched={}",
        pct(n_oxonium, n_found_spec),
        pct(n_searchable, n_found_spec),
        pct(n_sparse_searchable, n_sparse),
        n_found_spec
    );
}
