//! Fragment-intensity model predictions shared by the PSMs of one scored chunk.
//!
//! The strong-mode PSM features predict the fragment intensities of every retained
//! candidate. In a chunk of spectra sorted by precursor mass the same peptidoform is a
//! candidate of many spectra (on DIA, every cycle of its window around its retention time),
//! so its prediction is computed once per chunk here and reused.

use std::hash::BuildHasher;
use std::sync::{Arc, Mutex};

use model::peptide::Peptide;
use rustc_hash::FxHashMap;
use scoring_crate::gbdt::GbdtPeakModel;
use scoring_crate::psm::predict_frag_intensities;

/// Shards of the cache; well above the thread count, so lock waits stay rare.
const SHARDS: usize = 64;

/// One cached prediction: the model's address, the precursor charge, the peptidoform (its
/// equality compares residues, residue masses, modifications, flanks and charge, a superset
/// of what the model reads) and the predicted intensities.
struct Entry {
    model: usize,
    charge: u8,
    peptide: Peptide,
    intensities: Arc<[f64]>,
}

/// Entries by the hash of their model, charge and peptidoform; a bucket holds every entry
/// whose hash collides.
type Shard = FxHashMap<u64, Vec<Entry>>;

pub struct FragPredCache {
    shards: Vec<Mutex<Shard>>,
}

impl Default for FragPredCache {
    fn default() -> Self {
        Self {
            shards: (0..SHARDS).map(|_| Mutex::new(Shard::default())).collect(),
        }
    }
}

fn find(
    bucket: Option<&Vec<Entry>>,
    model: usize,
    charge: u8,
    peptide: &Peptide,
) -> Option<Arc<[f64]>> {
    bucket?
        .iter()
        .find(|e| e.model == model && e.charge == charge && e.peptide == *peptide)
        .map(|e| Arc::clone(&e.intensities))
}

impl FragPredCache {
    /// The intensities of [`predict_frag_intensities`] for `peptide` at `charge`, predicted
    /// on the first request of the chunk.
    pub fn intensities(&self, model: &GbdtPeakModel, peptide: &Peptide, charge: u8) -> Arc<[f64]> {
        let addr = model as *const GbdtPeakModel as usize;
        let h = rustc_hash::FxBuildHasher.hash_one((addr, charge, peptide));
        let shard = &self.shards[(h >> 58) as usize % SHARDS];
        if let Some(hit) = find(
            shard.lock().expect("frag cache lock").get(&h),
            addr,
            charge,
            peptide,
        ) {
            return hit;
        }
        let intensities: Arc<[f64]> = predict_frag_intensities(model, peptide, charge).into();
        let mut map = shard.lock().expect("frag cache lock");
        let bucket = map.entry(h).or_default();
        if let Some(raced) = find(Some(bucket), addr, charge, peptide) {
            return raced;
        }
        bucket.push(Entry {
            model: addr,
            charge,
            peptide: peptide.clone(),
            intensities: Arc::clone(&intensities),
        });
        intensities
    }

    /// Drops every cached prediction.
    pub fn clear(&self) {
        for s in &self.shards {
            *s.lock().expect("frag cache lock") = Shard::default();
        }
    }
}
