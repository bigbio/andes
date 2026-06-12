//! Diagnostic: dump main_ion picks per partition for a given param file.
//! Confirms whether the iter29 main_ion_from_param fix changes the dominant
//! ion for the dataset's bundled param.
use std::env;
use std::path::PathBuf;
use scoring::param_model::{Param, IonType};

fn main() {
    let path = env::args().nth(1).expect("usage: dump_main_ion <path/to/.param>");
    let param = Param::load_from_file(PathBuf::from(&path).as_path()).expect("load");
    println!("Param: {path}");
    println!("  num_segments={} num_partitions={}", param.num_segments, param.partitions.len());
    // Pick the (charge=2, seg=0) partition with the largest parent_mass
    // (representative of the bulk of the dataset).
    let mut seen: std::collections::BTreeMap<i32, Vec<f32>> = std::collections::BTreeMap::new();
    for p in &param.partitions {
        if p.seg_num != 0 { continue; }
        seen.entry(p.charge).or_default().push(p.parent_mass);
    }
    for (charge, mut masses) in seen {
        masses.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        // Print 3 representative masses: smallest, middle, largest
        let pick: Vec<f32> = vec![
            masses[0],
            masses[masses.len()/2],
            masses[masses.len()-1],
        ];
        for pm in pick {
            // The iter29 main_ion_from_param logic, replicated.
            let last_seg = (param.num_segments - 1).max(0) as usize;
            let part = param.partition_for(charge as u8, pm as f64, last_seg);
            // Aggregate frequencies across all segments for this (charge, parent_mass).
            let num_segs = param.num_segments.max(1) as usize;
            let mut ion_freq: std::collections::HashMap<IonType, f32> = std::collections::HashMap::new();
            for seg in 0..num_segs {
                let p = scoring::param_model::Partition { charge, parent_mass: part.parent_mass, seg_num: seg as i32 };
                if let Some(frags) = param.frag_off_table.get(&p) {
                    for f in frags {
                        if matches!(f.ion_type, IonType::Noise) { continue; }
                        *ion_freq.entry(f.ion_type).or_insert(0.0) += f.frequency;
                    }
                }
            }
            let mut entries: Vec<(IonType, f32)> = ion_freq.into_iter().collect();
            entries.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            print!("  charge={} pm={:.1} top_ions=", charge, pm);
            for (ion, freq) in entries.iter().take(3) {
                let kind = match ion {
                    IonType::Prefix { offset_bits, .. } => format!("b+{}", f32::from_bits(*offset_bits)),
                    IonType::Suffix { offset_bits, .. } => format!("y+{}", f32::from_bits(*offset_bits)),
                    IonType::Noise => "NOISE".to_string(),
                };
                print!("{}={:.4} ", kind, freq);
            }
            let main_kind = match entries.first().map(|(i, _)| i) {
                Some(IonType::Prefix { .. }) => "prefix (b-direction)",
                Some(IonType::Suffix { .. }) => "suffix (y-direction)",
                _ => "?",
            };
            println!("→ main_ion = {}", main_kind);
        }
    }
}
