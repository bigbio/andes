//! Fragment-intensity model predictions shared by the PSMs of one scored chunk.
//!
//! The strong-mode PSM features predict the fragment intensities of every retained
//! candidate. In a chunk of spectra sorted by precursor mass the same peptidoform is a
//! candidate of many spectra (on DIA, every cycle of its window around its retention time),
//! so its prediction is computed once per chunk here and reused.

use std::sync::{Arc, Mutex};

use model::peptide::Peptide;
use rustc_hash::FxHashMap;
use scoring_crate::gbdt_eval::GbdtPeakModel;
use scoring_crate::scoring::{frag_intensities_from_raw, predict_frag_raw};

/// Shards of the cache; well above the thread count, so lock waits stay rare.
const SHARDS: usize = 64;

/// Key: the model's address, the precursor charge and the peptidoform's text (residues,
/// modifications and flanks, a superset of what the model reads).
type Key = (usize, u8, String);

pub struct FragPredCache {
    shards: Vec<Mutex<FxHashMap<Key, Arc<[f32]>>>>,
}

impl Default for FragPredCache {
    fn default() -> Self {
        Self {
            shards: (0..SHARDS)
                .map(|_| Mutex::new(FxHashMap::default()))
                .collect(),
        }
    }
}

impl FragPredCache {
    /// The intensities of [`scoring_crate::scoring::predict_frag_intensities`] for `peptide`
    /// at `charge`, predicted on the first request of the chunk.
    pub fn intensities(&self, model: &GbdtPeakModel, peptide: &Peptide, charge: u8) -> Vec<f64> {
        use std::hash::BuildHasher;
        let key: Key = (
            model as *const GbdtPeakModel as usize,
            charge,
            peptide.to_string(),
        );
        let h = rustc_hash::FxBuildHasher.hash_one(&key);
        let shard = &self.shards[(h >> 58) as usize % SHARDS];
        let hit = shard.lock().expect("frag cache lock").get(&key).cloned();
        let raw = match hit {
            Some(raw) => raw,
            None => {
                let raw: Arc<[f32]> = predict_frag_raw(model, peptide, charge).into();
                shard
                    .lock()
                    .expect("frag cache lock")
                    .insert(key, Arc::clone(&raw));
                raw
            }
        };
        frag_intensities_from_raw(&raw)
    }

    /// Drops every cached prediction.
    pub fn clear(&self) {
        for s in &self.shards {
            *s.lock().expect("frag cache lock") = FxHashMap::default();
        }
    }
}
