//! Local proof: train a model on the tiny BSA fixture and compare its 1% FDR
//! yield vs the seed (fallback), plus dump the trained Noise distribution shape.
//! Run: cargo run -p model-train --example local_yield
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use input::MgfReader;
use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};
use scoring_crate::param_model::IonType;
use scoring_crate::{Param, RankScorer};
use search::SearchParams;

use model_train::{
    accumulate::{merge, StatsAccumulator},
    counts::CountStats,
    estimate::{Estimator, EstimatorConfig},
    gate::evaluate_candidate,
    labeled::bootstrap_labels,
    ModelStore,
};

fn standard_aa_set() -> model::AminoAcidSet {
    let cam = Modification { name: "Carbamidomethyl".into(), mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'), location: ModLocation::Anywhere, fixed: true, accession: None,
        neutral_losses: Vec::new() };
    let ox = Modification { name: "Oxidation".into(), mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'), location: ModLocation::Anywhere, fixed: false, accession: None,
        neutral_losses: Vec::new() };
    AminoAcidSetBuilder::new_standard().add_fixed_mod(cam).add_variable_mod(ox).build().unwrap()
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mgf = root.join("test-fixtures/test.mgf");
    let fasta = root.join("test-fixtures/BSA.fasta");

    let f = File::open(&mgf).expect("open BSA.mgf");
    let spectra: Vec<model::Spectrum> =
        MgfReader::new(BufReader::new(f)).filter_map(|r| r.ok()).collect();
    eprintln!("loaded {} spectra", spectra.len());

    let store = ModelStore::open(&root.join("resources/ionstat/models.parquet")).expect("store");
    let seed_param: Param = store.load_param("hcd_qexactive_tryp").expect("seed");
    let seed_scorer = RankScorer::new(&seed_param);

    let params = SearchParams::default_tryptic(standard_aa_set());
    let train_fdr = 0.1;
    let labels = bootstrap_labels(&spectra, &fasta, &seed_scorer, &params, train_fdr).expect("bootstrap");
    eprintln!("{} confident labels at {train_fdr} FDR", labels.len());

    let acc = StatsAccumulator::new(&seed_scorer);
    let mut stats = CountStats::new();
    for lm in &labels {
        acc.accumulate(&mut stats, &spectra[lm.spectrum_index], &lm.peptide, lm.charge);
    }
    let stats = merge(vec![stats]);
    let trained = Estimator::new(EstimatorConfig::default()).estimate(&stats, &seed_param);

    // Dump the trained Noise vs ion shape for the first partition that has both.
    for part in &trained.partitions {
        if let Some(t) = trained.rank_dist_table.get(part) {
            let noise = t.get(&IonType::Noise);
            let ion = t.iter().find(|(it, _)| matches!(it, IonType::Prefix { .. }));
            if let (Some(noise), Some((_, iv))) = (noise, ion) {
                let nsum: f32 = noise.iter().sum();
                let head: Vec<String> = noise.iter().take(4).map(|x| format!("{x:.4}")).collect();
                eprintln!("TRAINED noise: missing(last)={:.4}  sum={:.4}  head=[{}]",
                    noise.last().unwrap(), nsum, head.join(","));
                eprintln!("TRAINED ion:   missing(last)={:.4}  head=[{}]",
                    iv.last().unwrap(),
                    iv.iter().take(4).map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(","));
                let cos = 1.0_f32; // charge-1 prefix, num_segments>=1
                let n = iv.len().min(noise.len());
                eprintln!("  log head=[{}]  missing_log={:.2}",
                    (0..n.min(4)).map(|i| format!("{:.2}", (iv[i]/(noise[i]*cos)).ln())).collect::<Vec<_>>().join(","),
                    (iv[n-1]/(noise[n-1]*cos)).ln());
                break;
            }
        }
    }

    let trained_scorer = RankScorer::new(&trained);
    let delta = evaluate_candidate(&spectra, &fasta, &seed_scorer, &trained_scorer, &params, train_fdr)
        .expect("evaluate");
    eprintln!("\nSEED @{train_fdr} FDR: {}   TRAINED @{train_fdr} FDR: {}",
        delta.current_count, delta.candidate_count);
}
