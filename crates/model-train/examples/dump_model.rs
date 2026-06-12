//! Diagnostic: dump the rank_dist_table Noise vs ion distributions for a model,
//! to understand the numeric form the Estimator must reproduce.
//! Run: cargo run -p model-train --example dump_model -- <path-to.param>
use scoring_crate::param_model::{IonType, Param};

fn summarize(label: &str, v: &[f32]) {
    let n = v.len();
    let sum: f32 = v.iter().sum();
    let head: Vec<String> = v.iter().take(6).map(|x| format!("{x:.5}")).collect();
    let missing = v.last().copied().unwrap_or(f32::NAN);
    println!(
        "  {label:<22} len={n} sum={sum:.4} head=[{}] missing(last)={missing:.5}",
        head.join(", ")
    );
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: dump_model <.param>");
    let p = Param::load_from_file(std::path::Path::new(&path)).expect("load param");
    println!("model: {path}");
    println!("num_segments={} max_rank={}", p.num_segments, p.max_rank);
    println!("partitions={}", p.partitions.len());

    // Pick a few partitions that have a Noise entry + a prefix ion.
    let mut shown = 0;
    for part in &p.partitions {
        let Some(table) = p.rank_dist_table.get(part) else { continue };
        let Some(noise) = table.get(&IonType::Noise) else { continue };
        // find a prefix ion
        let prefix = table.iter().find(|(it, _)| matches!(it, IonType::Prefix { .. }));
        let Some((ion, ion_v)) = prefix else { continue };
        println!(
            "\npartition charge={} parent_mass={:.1} seg={}",
            part.charge, part.parent_mass, part.seg_num
        );
        summarize("NOISE", noise);
        summarize(&format!("{ion:?}"), ion_v);
        // resulting log score at a few ranks: ln(ion / (noise * min(charge,segs)))
        let ion_charge: i32 = match ion {
            IonType::Prefix { charge, .. } | IonType::Suffix { charge, .. } => *charge,
            _ => 1,
        };
        let cos = ion_charge.min(p.num_segments) as f32;
        let n = ion_v.len().min(noise.len());
        let logs: Vec<String> = (0..n.min(6))
            .map(|i| format!("{:.2}", (ion_v[i] / (noise[i] * cos)).ln()))
            .collect();
        let last = n - 1;
        println!(
            "  log_scores head=[{}]  missing_log={:.2}  (charge_or_seg={cos})",
            logs.join(", "),
            (ion_v[last] / (noise[last] * cos)).ln()
        );
        shown += 1;
        if shown >= 3 { break; }
    }
}
