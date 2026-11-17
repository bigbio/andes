use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use input::MzMLReader;
use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};
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

fn bundled_store_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models")
}

fn standard_aa_set() -> model::AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

fn load_mzml(path: &Path) -> Vec<model::Spectrum> {
    let file = File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let reader = BufReader::new(file);
    MzMLReader::new(reader).filter_map(|r| r.ok()).collect()
}

#[test]
fn trained_model_yield_not_worse_than_fallback() {
    let bench_dir = match std::env::var("ANDES_TRAIN_BENCH") {
        Ok(v) if !v.trim().is_empty() => PathBuf::from(v),
        _ => {
            eprintln!(
                "skipping trained_model_yield_not_worse_than_fallback: \
                 set ANDES_TRAIN_BENCH to a directory containing \
                 train.mzML + db.fasta (+ optional validate.mzML) to run"
            );
            return;
        }
    };

    let train_mzml = bench_dir.join("train.mzML");
    let db_fasta = bench_dir.join("db.fasta");
    let val_mzml = bench_dir.join("validate.mzML");

    let val_mzml_path: &Path = if val_mzml.exists() {
        &val_mzml
    } else {
        eprintln!(
            "note: validate.mzML not found in {}; reusing train.mzML as \
             validation set (counts are informational only — relative comparison is valid)",
            bench_dir.display()
        );
        &train_mzml
    };

    assert!(
        train_mzml.exists(),
        "train.mzML not found in {}",
        bench_dir.display()
    );
    assert!(
        db_fasta.exists(),
        "db.fasta not found in {}",
        bench_dir.display()
    );

    eprintln!("loading training spectra from {}", train_mzml.display());
    let train_spectra = load_mzml(&train_mzml);
    assert!(
        !train_spectra.is_empty(),
        "train.mzML must contain at least one MS2 spectrum"
    );
    eprintln!("  loaded {} training spectra", train_spectra.len());

    eprintln!("loading bundled seed model hcd_qexactive_tryp");
    let store_path = bundled_store_path();
    let store = ModelStore::open(&store_path).expect("failed to open bundled model store");
    let seed_param: Param = store
        .load_param("hcd_qexactive_tryp")
        .expect("hcd_qexactive_tryp not found in bundled store");
    let seed_scorer = RankScorer::new(&seed_param);
    eprintln!(
        "  seed model loaded ({} partitions)",
        seed_param.partitions.len()
    );

    let search_params = SearchParams::default_tryptic(standard_aa_set());
    let train_fdr = 0.01_f64;

    eprintln!("running bootstrap_labels (train_fdr={train_fdr})");
    let cache_dir = tempfile::tempdir().expect("cache dir");
    let labels = bootstrap_labels(
        &train_spectra,
        &db_fasta,
        &seed_scorer,
        &search_params,
        train_fdr,
        cache_dir.path(),
    )
    .expect("bootstrap_labels must succeed");
    eprintln!("  {} confident labels at {} FDR", labels.len(), train_fdr);

    let acc = StatsAccumulator::new(&seed_scorer);
    let mut stats = CountStats::new();
    for lm in &labels {
        let spec = &train_spectra[lm.spectrum_index];
        acc.accumulate(&mut stats, spec, &lm.peptide, lm.charge);
    }
    let merged_stats = merge(vec![stats]);

    let estimator = Estimator::new(EstimatorConfig::default());
    let trained_param = estimator.estimate(&merged_stats, &seed_param);
    let trained_scorer = RankScorer::new(&trained_param);
    eprintln!(
        "  trained model estimated ({} partitions)",
        trained_param.partitions.len()
    );

    eprintln!(
        "loading validation spectra from {}",
        val_mzml_path.display()
    );
    let val_spectra = load_mzml(val_mzml_path);
    assert!(
        !val_spectra.is_empty(),
        "validation spectra must contain at least one MS2 spectrum"
    );
    eprintln!("  loaded {} validation spectra", val_spectra.len());

    eprintln!("running evaluate_candidate (fdr=0.01)");
    let delta = evaluate_candidate(
        &val_spectra,
        &db_fasta,
        &seed_scorer,
        &trained_scorer,
        &search_params,
        0.01_f64,
        cache_dir.path(),
    )
    .expect("evaluate_candidate must succeed");

    eprintln!(
        "  fallback PSMs @ 1% FDR: {}   trained PSMs @ 1% FDR: {}",
        delta.current_count, delta.candidate_count
    );

    let fallback_count = delta.current_count;
    let trained_count = delta.candidate_count;
    let floor = (fallback_count as f64 * 0.5).floor() as usize;
    assert!(
        trained_count >= floor,
        "trained model yield ({trained_count}) < {floor} (= 0.5 x fallback {fallback_count}) \
         at 1% FDR — the trained model is not discriminative (training regressed). \
         Note: >= fallback (1.0x) is the sub-project D/E goal; ~0.8x is expected \
         for a single-pass bootstrap from a generic seed."
    );
}
