//! Diagnostic companion to glyco_probe: for each truth scan, print the nearest
//! candidate backbone to truth (ppm error), candidate count, top candidate, and
//! whether a Y0/Y-ladder peak for the TRUE backbone actually exists in the
//! spectrum. Used to distinguish "solver genuinely can't find it" from a
//! harness/tolerance artifact.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use andes_glyco::backbone::solve_backbone;
use andes_glyco::glycan_mass::{CORE_Y_STEPS, PROTON};
use andes_glyco::oxonium::oxonium_gate;
use input::mzml::MzMLReader;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mzml_path = &args[1];
    let truth_path = &args[2];
    let limit: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(40);

    let mut truth: Vec<(i32, f64, f64, u8)> = Vec::new();
    {
        let f = File::open(truth_path).unwrap();
        for (i, line) in BufReader::new(f).lines().enumerate() {
            let line = line.unwrap();
            if i == 0 || line.trim().is_empty() {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            truth.push((
                c[0].parse().unwrap(),
                c[1].parse().unwrap(),
                c[2].parse().unwrap(),
                c[3].parse().unwrap(),
            ));
        }
    }

    let f = File::open(mzml_path).unwrap();
    let (spectra, _) = MzMLReader::new(BufReader::new(f))
        .with_ms_level_range(2, 2)
        .read_with_ms1()
        .unwrap();
    let mut by_scan: HashMap<i32, &model::spectrum::Spectrum> = HashMap::new();
    for s in &spectra {
        if let Some(sc) = s.scan {
            by_scan.insert(sc, s);
        }
    }

    // Tally: does a Y0 or any core-Y peak for the TRUE backbone exist in the
    // spectrum (within 20 ppm, at any charge 1..=z)?
    let mut have_true_y0 = 0usize;
    let mut have_true_yladder = 0usize; // >=2 core-Y rungs present
    let mut nearest_within_ppm = [0usize; 4]; // <20, <50, <200, <1000 ppm of nearest candidate
    let mut n = 0usize;
    let mut printed = 0usize;

    for &(scan, bb_true, pmz, z) in &truth {
        let spec = match by_scan.get(&scan) {
            Some(s) => *s,
            None => continue,
        };
        n += 1;
        let peaks = &spec.peaks;
        let ox = oxonium_gate(peaks, 0.10, 20.0);
        let z = z.max(1);
        let precursor_neutral = (pmz - PROTON) * z as f64;

        // Does the spectrum contain the true Y0 / Y-ladder peaks?
        let peak_exists = |target_neutral: f64| -> bool {
            for &(p, _) in peaks {
                for zc in 1..=z {
                    let pn = (p - PROTON) * zc as f64;
                    let tol = (target_neutral * 20e-6).max(0.01);
                    if (pn - target_neutral).abs() <= tol {
                        return true;
                    }
                }
            }
            false
        };
        let y0 = peak_exists(bb_true);
        if y0 {
            have_true_y0 += 1;
        }
        let rungs = CORE_Y_STEPS.iter().filter(|&&s| peak_exists(bb_true + s)).count()
            + if y0 { 1 } else { 0 };
        if rungs >= 2 {
            have_true_yladder += 1;
        }

        let cands = solve_backbone(peaks, precursor_neutral, z, 20.0, 5);
        let nearest = cands
            .iter()
            .map(|c| (c.backbone_mass - bb_true).abs() / bb_true * 1e6)
            .fold(f64::INFINITY, f64::min);
        if nearest < 20.0 {
            nearest_within_ppm[0] += 1;
        }
        if nearest < 50.0 {
            nearest_within_ppm[1] += 1;
        }
        if nearest < 200.0 {
            nearest_within_ppm[2] += 1;
        }
        if nearest < 1000.0 {
            nearest_within_ppm[3] += 1;
        }

        if printed < limit && nearest >= 20.0 {
            printed += 1;
            let top = cands
                .first()
                .map(|c| format!("{:.3}(h{},v{})", c.backbone_mass, c.core_y_hits, c.votes))
                .unwrap_or_else(|| "NONE".into());
            println!(
                "scan {:>6} z{} bb_true={:.3} ox={} ncand={} nearest_ppm={:.0} trueY0={} rungs={} top={}",
                scan, z, bb_true, ox.fired, cands.len(), nearest, y0, rungs, top
            );
        }
    }

    println!("--- summary over {} scans ---", n);
    println!("true Y0 peak present:        {} ({:.1}%)", have_true_y0, 100.0 * have_true_y0 as f64 / n as f64);
    println!("true Y-ladder (>=2 rungs):   {} ({:.1}%)", have_true_yladder, 100.0 * have_true_yladder as f64 / n as f64);
    println!("nearest candidate <20ppm:    {} ({:.1}%)", nearest_within_ppm[0], 100.0 * nearest_within_ppm[0] as f64 / n as f64);
    println!("nearest candidate <50ppm:    {} ({:.1}%)", nearest_within_ppm[1], 100.0 * nearest_within_ppm[1] as f64 / n as f64);
    println!("nearest candidate <200ppm:   {} ({:.1}%)", nearest_within_ppm[2], 100.0 * nearest_within_ppm[2] as f64 / n as f64);
    println!("nearest candidate <1000ppm:  {} ({:.1}%)", nearest_within_ppm[3], 100.0 * nearest_within_ppm[3] as f64 / n as f64);
}
