//! Manual regression check for the glycan-first pre-filter against real mouse
//! glycopeptide spectra. Gated on two env vars (the data is not vendored):
//!
//! ```text
//! ANDES_GLYCO_GDB=/path/to/pGlyco-N-Mouse.gdb \
//! ANDES_GLYCO_MGF=/path/to/glycopeptide-filtered.mgf \
//! cargo test -p andes-glyco --test glycan_first_mouse -- --nocapture
//! ```
//!
//! It loads the glycan DB, builds the Y-complementary ion index, stream-parses
//! the first `N` spectra from the MGF, and runs `search_glycans` on each. It is a
//! sanity check, not a gold-standard assertion: a glycopeptide spectrum should
//! fire the diagnostic gate and (when its core Y ions are present) retrieve
//! candidates.

use std::env;
use std::io::BufRead;

use andes_glyco::glycan_db::{load_glycan_gdb, GlycanComp};
use andes_glyco::glycan_first::{search_glycans, Glycan, GlycanConfig, GlycanCore, GlycanIonIndex};
use andes_glyco::glycan_mass::PROTON;

struct MgfSpectrum {
    pepmass: f64,
    charge: u8,
    peaks: Vec<(f64, f32)>,
}

/// Stream-parse at most `max_n` MGF spectra. Returns as soon as `max_n` complete
/// `BEGIN IONS … END IONS` blocks are seen (so a 500 MB filtered MGF is not
/// loaded whole).
fn parse_mgf_limited(path: &str, max_n: usize) -> Vec<MgfSpectrum> {
    let mut out = Vec::new();
    let mut cur: Option<MgfSpectrum> = None;
    let file = std::fs::File::open(path).expect("open MGF");
    for line in std::io::BufReader::new(file).lines() {
        let line = line.expect("read line");
        let line = line.trim();
        if line.starts_with("BEGIN IONS") {
            cur = Some(MgfSpectrum {
                pepmass: 0.0,
                charge: 0,
                peaks: Vec::new(),
            });
        } else if line.starts_with("END IONS") {
            if let Some(s) = cur.take() {
                out.push(s);
                if out.len() >= max_n {
                    break;
                }
            }
        } else if let Some(s) = cur.as_mut() {
            if let Some(v) = line.strip_prefix("PEPMASS=") {
                s.pepmass = v
                    .split_whitespace()
                    .next()
                    .and_then(|x| x.parse().ok())
                    .unwrap_or(0.0);
            } else if let Some(v) = line.strip_prefix("CHARGE=") {
                let v = v.trim().trim_end_matches('+').trim_end_matches('-');
                s.charge = v.parse().unwrap_or(0);
            } else if let Some((mz, inten)) = line.split_once(' ') {
                if let (Ok(mz), Ok(inten)) = (mz.trim().parse::<f64>(), inten.trim().parse::<f32>())
                {
                    s.peaks.push((mz, inten));
                }
            }
        }
    }
    out
}

#[test]
fn glycan_first_retrieves_candidates_on_mouse_spectra() {
    // Skip (not fail) when the data env vars are absent — the mouse data is not
    // vendored, so this test must be a no-op in normal CI.
    let (Ok(gdb), Ok(mgf)) = (env::var("ANDES_GLYCO_GDB"), env::var("ANDES_GLYCO_MGF")) else {
        eprintln!("skipping: ANDES_GLYCO_GDB / ANDES_GLYCO_MGF not set");
        return;
    };

    let comps: Vec<GlycanComp> = load_glycan_gdb(&std::fs::read_to_string(&gdb).unwrap()).unwrap();
    let glycans: Vec<Glycan> = comps
        .iter()
        .enumerate()
        .map(|(i, c)| Glycan {
            id: i as u32,
            composition: c.clone(),
        })
        .collect();

    let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
    let cfg = GlycanConfig {
        min_core_ions: 1, // permissive pre-filter, matches the driver
        top_k: 100,
        ..GlycanConfig::default()
    };

    let spectra = parse_mgf_limited(&mgf, 200);
    let mut fired = 0usize;
    let mut nonempty = 0usize;
    let mut total_candidates = 0usize;

    for spec in &spectra {
        let charge = if spec.charge > 0 { spec.charge } else { 2 };
        let neutral = (spec.pepmass - PROTON) * charge as f64;
        let mut peaks = spec.peaks.clone();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let res = search_glycans(&idx, &cfg, &comps, neutral, charge, &peaks);
        if res.is_glycopeptide {
            fired += 1;
            if !res.candidates.is_empty() {
                nonempty += 1;
                total_candidates += res.candidates.len();
            }
        }
        if res.is_glycopeptide && !res.candidates.is_empty() {
            let top = &res.candidates[0];
            let comp = &comps[top.glycan_id as usize];
            println!(
                "scan pepmass={} z={} neutral={:.3}: top glycan H{hex}N{hexnac}F{fuc}A{neuac}G{neugc} (mass {mass:.3}), score {score}, y={y} core={core}",
                spec.pepmass,
                charge,
                neutral,
                hex = comp.hex,
                hexnac = comp.hexnac,
                fuc = comp.fuc,
                neuac = comp.neuac,
                neugc = comp.neugc,
                mass = comp.mass,
                score = top.score,
                y = top.matched_y_ions,
                core = top.matched_core_ions,
            );
        }
    }

    println!(
        "glycan-first on mouse: {}/{} spectra passed diagnostic gate, {}/{} retrieved >=1 candidate, {} total candidates",
        fired,
        spectra.len(),
        nonempty,
        fired,
        total_candidates
    );
    // Sanity: on a glycopeptide-filtered run the gate should fire on the large
    // majority, and a substantial fraction should retrieve candidates.
    assert!(fired > 0, "diagnostic gate never fired — data or parser problem");
    assert!(
        nonempty > 0,
        "ion index retrieved nothing for every gated spectrum"
    );
}
