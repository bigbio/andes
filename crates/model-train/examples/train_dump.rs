//! Diagnostic: train on a real dataset and dump the trained model's
//! distributions vs the seed for well-populated partitions, to localize the
//! calibration gap. Run:
//!   cargo run --release -p model-train --example train_dump -- <spectra> <fasta> [seed_slug]
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use input::{MgfReader, MzMLReader};
use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};
use scoring_crate::param_model::{IonType, Param};
use scoring_crate::RankScorer;
use search::SearchParams;

use model_train::{
    accumulate::{merge, StatsAccumulator},
    counts::CountStats,
    estimate::{Estimator, EstimatorConfig},
    labeled::bootstrap_labels,
    ModelStore,
};

fn aa_set() -> model::AminoAcidSet {
    let cam = Modification { name: "Carbamidomethyl".into(), mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'), location: ModLocation::Anywhere, fixed: true, accession: None,
        neutral_losses: Vec::new() };
    let ox = Modification { name: "Oxidation".into(), mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'), location: ModLocation::Anywhere, fixed: false, accession: None,
        neutral_losses: Vec::new() };
    AminoAcidSetBuilder::new_standard().add_fixed_mod(cam).add_variable_mod(ox).build().unwrap()
}

fn load(path: &Path) -> Vec<model::Spectrum> {
    let f = File::open(path).expect("open spectra");
    if path.extension().and_then(|e| e.to_str()).map(|e| e.eq_ignore_ascii_case("mzml")).unwrap_or(false) {
        MzMLReader::new(BufReader::new(f)).filter_map(|r| r.ok()).collect()
    } else {
        MgfReader::new(BufReader::new(f)).filter_map(|r| r.ok()).collect()
    }
}

fn row(label: &str, v: &[f32]) {
    let head: Vec<String> = v.iter().take(4).map(|x| format!("{x:.4}")).collect();
    let sum: f32 = v.iter().sum();
    println!("    {label:<8} sum={sum:.3} head=[{}] missing={:.4}", head.join(","), v.last().copied().unwrap_or(f32::NAN));
}

fn dump(tag: &str, p: &Param, part: &scoring_crate::param_model::Partition) {
    let Some(t) = p.rank_dist_table.get(part) else { println!("  [{tag}] no table"); return };
    let noise = t.get(&IonType::Noise);
    let ion = t.iter().find(|(it, _)| matches!(it, IonType::Prefix { .. }));
    println!("  [{tag}]");
    if let Some(n) = noise { row("noise", n); }
    if let Some((_, iv)) = ion { row("ion", iv); }
    if let (Some(n), Some((_, iv))) = (noise, ion) {
        let k = n.len().min(iv.len());
        let logs: Vec<String> = (0..k.min(4)).map(|i| format!("{:.2}", (iv[i]/(n[i].max(1e-9))).ln())).collect();
        println!("    log      head=[{}] missing={:.2}", logs.join(","), (iv[k-1]/(n[k-1].max(1e-9))).ln());
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!(
            "usage: {} <spectra.mzML|.mgf> <database.fasta> [seed_slug]",
            a.first().map(|s| s.as_str()).unwrap_or("train_dump")
        );
        std::process::exit(2);
    }
    let spectra_path = Path::new(&a[1]);
    let fasta = Path::new(&a[2]);
    let seed_slug = a.get(3).map(|s| s.as_str()).unwrap_or("hcd_qexactive_tryp");

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let store = ModelStore::open(&root.join("resources/ionstat/models.parquet")).expect("store");
    let seed_param: Param = store.load_param(seed_slug).expect("seed");
    let seed_scorer = RankScorer::new(&seed_param);

    let spectra = load(spectra_path);
    eprintln!("loaded {} spectra", spectra.len());
    let params = SearchParams::default_tryptic(aa_set());
    let labels = bootstrap_labels(&spectra, fasta, &seed_scorer, &params, 0.01).expect("bootstrap");
    eprintln!("{} labels", labels.len());

    let acc = StatsAccumulator::new(&seed_scorer);
    let mut stats = CountStats::new();
    for lm in &labels { acc.accumulate(&mut stats, &spectra[lm.spectrum_index], &lm.peptide, lm.charge); }
    let stats = merge(vec![stats]);
    let trained = Estimator::new(EstimatorConfig::default()).estimate(&stats, &seed_param);

    // Pick the 3 partitions with the most prefix-ion counts. Aggregate per
    // partition first (stats.rank is keyed by (Partition, IonType), so a
    // partition has several Prefix entries) so each partition appears once.
    let mut per_partition: std::collections::HashMap<scoring_crate::param_model::Partition, u64> =
        std::collections::HashMap::new();
    for ((p, it), v) in stats.rank.iter() {
        if matches!(it, IonType::Prefix { .. }) {
            *per_partition.entry(*p).or_default() += v.iter().sum::<u64>();
        }
    }
    let mut by_count: Vec<(scoring_crate::param_model::Partition, u64)> =
        per_partition.into_iter().collect();
    by_count.sort_by(|a, b| b.1.cmp(&a.1));
    for (part, cnt) in by_count.iter().take(3) {
        println!("\npartition charge={} mass={:.0} seg={}  ion_count={cnt}", part.charge, part.parent_mass, part.seg_num);
        dump("SEED", &seed_param, part);
        dump("TRAINED", &trained, part);
    }
}
