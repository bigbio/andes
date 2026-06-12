//! Diagnostic for prefix_score_cache[1087] == 0.0 bug.
//!
//! Loads HCD_QExactive_Tryp.param + PXD001819 mzML, finds scan=28787,
//! builds ScoredSpectrum at charge=2, then dumps:
//!  - per-segment ion-type list
//!  - per-prefix-ion theo_mz/segment for nominal in {974, 1087, 1216, 1345, 1561, 1920}
//!  - replicates `directional_node_score_inner` logic in user space (sums score per seg)
//!  - the production cached prefix score for cross-check
//!
//! Run:
//!   cargo run --release -p scoring --example dump_prefix_cache

use std::fs::File;
use std::io::BufReader;

use input::MzMLReader;
use model::tolerance::Tolerance;
use scoring::param_model::{IonType, Param};
use scoring::scoring::rank_scorer::RankScorer;
use scoring::scoring::scored_spectrum::ScoredSpectrum;

const PARAM_PATH: &str =
    "/Users/yperez/work/msgfplus-workspace/astral-speed-score-fix/resources/ionstat/CID_HighRes_Tryp.param";
const MZML_PATH: &str =
    "/Users/yperez/work/msgfplus-workspace/benchmark/data/PXD001819/UPS1_5000amol_R1.mzML";
const TARGET_SCAN: i32 = 28787;
const CHARGE: u8 = 2;

fn ion_label(ion: &IonType) -> String {
    match ion {
        IonType::Prefix { charge, offset_bits } => {
            format!("Prefix(c={},off={:.5})", charge, f32::from_bits(*offset_bits))
        }
        IonType::Suffix { charge, offset_bits } => {
            format!("Suffix(c={},off={:.5})", charge, f32::from_bits(*offset_bits))
        }
        IonType::Noise => "Noise".into(),
    }
}

fn mme_as_da(mme: &Tolerance, mz: f64) -> f64 {
    mme.as_da(mz)
}

fn main() {
    let param = Param::load_from_file(std::path::Path::new(PARAM_PATH)).expect("load param");
    let scorer = RankScorer::new(&param);
    println!("== Param ==");
    println!("num_segments     = {}", param.num_segments);
    println!("max_rank         = {}", param.max_rank);
    println!("mme              = {:?}", param.mme);

    println!("\n== ALL partitions (charge=2 only) ==");
    for p in &param.partitions {
        if p.charge != 2 {
            continue;
        }
        let logs = scorer.partition_ion_logs(p);
        let n_prefix = logs.iter().filter(|(ion, _)| ion.is_prefix()).count();
        let n_suffix = logs.iter().filter(|(ion, _)| ion.is_suffix()).count();
        println!(
            "  c={} pm={:.4} seg={} ions={} pfx={} sfx={}",
            p.charge, p.parent_mass, p.seg_num, logs.len(), n_prefix, n_suffix
        );
    }

    println!("\n== Reading mzML for scan={} ==", TARGET_SCAN);
    let f = File::open(MZML_PATH).expect("open mzML");
    let reader = MzMLReader::new(BufReader::new(f));
    let mut found = None;
    for spec_res in reader {
        let spec = spec_res.expect("parse spectrum");
        if spec.scan == Some(TARGET_SCAN) {
            found = Some(spec);
            break;
        }
    }
    let spec = found.expect("scan 28787 not found");
    let parent_mass = (spec.precursor_mz - 1.00727649) * (CHARGE as f64);
    println!("precursor_mz     = {:.5}", spec.precursor_mz);
    println!("parent_mass      = {:.5}", parent_mass);
    println!("peak_count       = {}", spec.peaks.len());

    let ss = ScoredSpectrum::new(&spec, &scorer, CHARGE);

    let num_segs = param.num_segments as usize;
    println!("\n== Per-segment partitions for THIS spectrum ==");
    let mut cached_ion_logs: Vec<Vec<(IonType, Vec<f32>)>> = Vec::with_capacity(num_segs);
    for seg in 0..num_segs {
        let p = param.partition_for(CHARGE, parent_mass, seg);
        let logs = scorer.partition_ion_logs(&p).to_vec();
        let n_prefix = logs.iter().filter(|(ion, _)| ion.is_prefix()).count();
        let n_suffix = logs.iter().filter(|(ion, _)| ion.is_suffix()).count();
        println!(
            "seg={} partition=(c={}, pm={:.3}, seg={}) total_ions={} prefix={} suffix={}",
            seg, p.charge, p.parent_mass, p.seg_num, logs.len(), n_prefix, n_suffix
        );
        for (ion, _logs) in &logs {
            println!("    {}", ion_label(ion));
        }
        cached_ion_logs.push(logs);
    }

    let max_rank = scorer.max_rank();
    let max_rank_idx = max_rank as usize;
    let mme = &param.mme;

    let targets = [974.0_f64, 1087.0, 1216.0, 1345.0, 1561.0, 1920.0];
    for nominal_mass in targets {
        println!("\n== nominal_mass = {:.1} (is_prefix=true) ==", nominal_mass);
        let mut total = 0.0_f32;
        let mut any_iter = false;
        for (seg, logs_slice) in cached_ion_logs.iter().enumerate().take(num_segs) {
            for (ion, logs) in logs_slice {
                if !ion.is_prefix() {
                    continue;
                }
                let theo_mz = ion.mz(nominal_mass);
                let seg_for_theo = param.segment_num(theo_mz, parent_mass);
                let in_segment = seg_for_theo == seg;
                let tol_da = mme_as_da(mme, theo_mz);
                let rank = ss.nearest_peak_rank(theo_mz, tol_da);
                let contribution_label;
                let contribution: f32 = if !in_segment {
                    contribution_label = "SKIP(seg mismatch)".to_string();
                    0.0
                } else {
                    any_iter = true;
                    match rank {
                        Some(r) => {
                            let idx = (r.min(max_rank).max(1) as usize) - 1;
                            if idx < logs.len() {
                                contribution_label = format!("matched rank={} idx={} score={:.4}", r, idx, logs[idx]);
                                logs[idx]
                            } else {
                                contribution_label = format!("matched rank={} but idx {} >= logs.len()={}", r, idx, logs.len());
                                0.0
                            }
                        }
                        None => {
                            if max_rank_idx < logs.len() {
                                contribution_label = format!("no peak; miss-slot[{}]={:.4}", max_rank_idx, logs[max_rank_idx]);
                                logs[max_rank_idx]
                            } else {
                                contribution_label = format!("no peak; miss-slot {} >= logs.len()={}", max_rank_idx, logs.len());
                                0.0
                            }
                        }
                    }
                };
                if in_segment {
                    total += contribution;
                }
                println!(
                    "  seg={} ion={} theo_mz={:.4} seg(theo)={} {} tol_da={:.4} | {}",
                    seg, ion_label(ion), theo_mz, seg_for_theo,
                    if in_segment { "(IN)" } else { "(OUT)" },
                    tol_da, contribution_label
                );
            }
        }
        let nominal_i32 = nominal_mass as i32;
        let cached = ss.cached_prefix_score(nominal_i32);
        println!(
            "  -> replicated_total={:.4} (any_in_segment_iter={}) cached_prefix_score({})={:?}",
            total, any_iter, nominal_i32, cached
        );
    }
}
