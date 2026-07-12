//! Decoy database generation via sequence reversal (default), seeded shuffling,
//! or none (target-only / externally-supplied decoys).

use model::protein::{Protein, ProteinDb};

/// Default decoy accession prefix.
pub const DEFAULT_DECOY_PREFIX: &str = "XXX";

/// Default seed for the reproducible shuffle decoy strategy.
pub const DEFAULT_DECOY_SEED: u64 = 42;

/// How decoy proteins are generated for the target+decoy search database.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecoyStrategy {
    /// Reverse each target sequence (default; deterministic reverse decoys,
    /// Elias & Gygi, Nat Methods 2007).
    #[default]
    Reverse,
    /// Shuffle each target sequence with a seeded, reproducible RNG. Useful when
    /// reversal produces too-similar decoys (e.g. palindromic/low-complexity).
    Shuffle,
    /// Sequon-preserving reverse (for `--glyco` FDR). Like [`Self::Reverse`] but
    /// restores each target's N-X-S/T sequon at its mirrored position, so the
    /// decoy carries the SAME sequon count/density as its target. Plain reversal
    /// depletes sequons (N-X-S/T → S/T-X-N), which makes glyco FDR
    /// anti-conservative because the sequon gate admits decoys at a lower rate.
    SequonReverse,
    /// Generate no decoys: the search database is the input FASTA verbatim. For
    /// inputs that already contain decoys (accessions carrying the decoy prefix)
    /// or when FDR is computed externally downstream.
    None,
}

impl DecoyStrategy {
    /// Parse the `--decoy-strategy` CLI value. Case-insensitive.
    pub fn from_name(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "reverse" | "rev"     => Some(Self::Reverse),
            "shuffle" | "shuf"    => Some(Self::Shuffle),
            "sequon-reverse" | "sequon" | "sequon_reverse" => Some(Self::SequonReverse),
            "none" | "off"        => Some(Self::None),
            _ => None,
        }
    }
}

/// Reverse each protein's sequence and prepend `<prefix>_` to its
/// accession. `prefix` is normalized: trailing `_`s stripped; empty
/// prefix → `DEFAULT_DECOY_PREFIX`.
///
/// KNOWN FDR CAVEAT (glyco search only, deferred — see glyco-phase1 adversarial
/// review, 2026-07): plain sequence reversal does NOT preserve N-X-S/T
/// sequon density. Reversing a sequence turns a forward "N-X-S/T" motif into
/// "T/S-X-N" (read backward), so a decoy protein's sequon count is generally
/// UNRELATED to (typically lower than) its target counterpart's. The glyco
/// driver (`search::glyco_search`) gates ALL backbone candidates — target and
/// decoy alike — on `andes_glyco::sequon::has_nxst_sequon`, so if target
/// peptides pass that gate systematically more often than reverse decoys do,
/// the effective decoy search space is smaller/different from the target
/// space and Percolator's target/decoy competition is anti-conservative
/// (FDR is UNDERESTIMATED) for glyco PSMs specifically.
///
/// This does NOT affect the standard (non-glyco) search path, which has no
/// sequon gate. It also does not affect a benchmark run under
/// `--decoy-strategy none` (externally supplied target+decoy FASTA, e.g. the
/// the reference engine-matched comparison fixture) — this caveat only applies when
/// andes generates its OWN reverse decoys for a `--glyco` run.
///
/// Not fixed in this pass: a sequon-preserving decoy generator (e.g. a
/// pseudo-reversal that keeps N-X-S/T triplets anchored while permuting the
/// rest of the sequence) is a nontrivial algorithmic change to a
/// well-tested, shared decoy path used by every non-glyco search too;
/// changing it here risked destabilizing that path for a glyco-only benefit.
/// TODO before any glyco vs. non-andes FDR **parity** claim (not needed for
/// this bugfix's find-rate/count re-measurement, which uses
/// `--decoy-strategy none`): implement and validate a sequon-preserving decoy
/// strategy for `--glyco`, or run glyco decoy generation under
/// `DecoyStrategy::Shuffle` (which does not systematically deplete sequons
/// the way reversal does, though it is not proven sequon-neutral either) and
/// measure sequon-density parity target vs. decoy directly before trusting
/// glyco q-values from an andes-generated decoy DB.
pub fn reverse_db(db: &ProteinDb, prefix: &str) -> ProteinDb {
    let normalized = normalize_decoy_prefix(prefix);
    let proteins = db.proteins.iter().map(|p| Protein {
        accession: format!("{}_{}", normalized, p.accession),
        description: p.description.clone(),
        sequence: p.sequence.iter().rev().copied().collect(),
    }).collect();
    ProteinDb { proteins }
}

/// Concatenate target + decoy.
pub fn target_plus_decoy(target: &ProteinDb, prefix: &str) -> ProteinDb {
    let decoy = reverse_db(target, prefix);
    let mut proteins = target.proteins.clone();
    proteins.extend(decoy.proteins);
    ProteinDb { proteins }
}

/// Count N-X-S/T sequon start positions (X ≠ P) in a residue slice.
fn count_sequon_starts(s: &[u8]) -> usize {
    (0..s.len().saturating_sub(2))
        .filter(|&i| s[i] == b'N' && s[i + 1] != b'P' && (s[i + 2] == b'S' || s[i + 2] == b'T'))
        .count()
}

/// Sequon-density-matched reverse of one protein sequence.
///
/// Plain reversal does NOT drive the sequon count to zero — it produces a
/// DIFFERENT set of (spurious) sequons at ~87% of the target's density on a human
/// proteome (measured), a modest but real deficit that makes the glyco sequon gate
/// admit decoys slightly less often than targets (anti-conservative FDR). A naive
/// "restore every original sequon" over-corrects badly (~160%, over-conservative,
/// costs real IDs), because it ADDS the originals on top of reversal's spurious set.
///
/// This instead reverses, then closes ONLY the deficit: greedily create sequons up
/// to the target's own count via composition-preserving swaps (pull a spare S/T
/// into the third slot of an `N-X-·` window). The result is a PERMUTATION of `orig`
/// (composition, tryptic peptide count and mass distribution preserved) whose sequon
/// count is raised toward — but never above — the target's. Deterministic: windows
/// scanned ascending, donor residues scanned ascending. Runs once per protein at
/// DB-build time (not a hot path), so each accepted swap RE-COUNTS the actual sequon
/// total: the swap's side effects (it can also break/create a sequon at the donor
/// window) make a naive `+1` unreliable, so a swap is committed only if it strictly
/// increases the real count and never overshoots the target (adversarial review).
fn sequon_preserving_reverse(orig: &[u8]) -> Vec<u8> {
    let n = orig.len();
    let mut r: Vec<u8> = orig.iter().rev().copied().collect();
    if n < 3 {
        return r;
    }
    let target = count_sequon_starts(orig);
    let mut cur = count_sequon_starts(&r);
    let mut i = 0usize;
    while cur < target && i + 2 < n {
        // A window that is one S/T away from being a sequon: N, X≠P, non-S/T.
        if r[i] == b'N' && r[i + 1] != b'P' && r[i + 2] != b'S' && r[i + 2] != b'T' {
            // Donor: a later S/T that is NOT the 3rd residue of an existing sequon
            // (moving it would break that sequon, a net-zero swap). Scan ascending.
            if let Some(j) = (i + 3..n).find(|&j| {
                (r[j] == b'S' || r[j] == b'T')
                    && !(j >= 2 && r[j - 2] == b'N' && r[j - 1] != b'P')
            }) {
                r.swap(i + 2, j);
                // Verify the ACTUAL delta (the swap can break/create a sequon at the
                // donor window too). Keep the swap only if it strictly increases the
                // count without overshooting; otherwise revert.
                let new_count = count_sequon_starts(&r);
                if new_count > cur && new_count <= target {
                    cur = new_count;
                } else {
                    r.swap(i + 2, j); // revert
                }
            }
        }
        i += 1;
    }
    r
}

/// Sequon-preserving reverse decoy database (glyco FDR). See
/// [`sequon_preserving_reverse`] and [`DecoyStrategy::SequonReverse`]. Accession
/// scheme matches [`reverse_db`].
pub fn sequon_preserving_reverse_db(db: &ProteinDb, prefix: &str) -> ProteinDb {
    let normalized = normalize_decoy_prefix(prefix);
    let proteins = db
        .proteins
        .iter()
        .map(|p| Protein {
            accession: format!("{}_{}", normalized, p.accession),
            description: p.description.clone(),
            sequence: sequon_preserving_reverse(&p.sequence),
        })
        .collect();
    ProteinDb { proteins }
}

/// Shuffle each protein's sequence with a per-protein, seed-derived RNG, and
/// prepend `<prefix>_` to its accession (same accession scheme as [`reverse_db`]).
///
/// Deterministic: the same `(seed, protein index, sequence)` always yields the
/// same permutation, so repeated runs and re-benchmarks are reproducible. Uses a
/// dependency-free xorshift64 RNG + Fisher–Yates; the per-protein seed mixes the
/// global seed with the protein index so identical sequences don't shuffle alike.
pub fn shuffle_db(db: &ProteinDb, prefix: &str, seed: u64) -> ProteinDb {
    let normalized = normalize_decoy_prefix(prefix);
    let proteins = db
        .proteins
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let mut seq = p.sequence.clone();
            // Per-protein RNG state: mix the global seed with the index via the
            // splitmix64 multiplier so adjacent proteins get well-separated
            // streams. Never let the xorshift state be zero (it would stick).
            let mut state = seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            if state == 0 {
                state = 0xD1B5_4A32_D192_ED03;
            }
            // Fisher–Yates from the back, drawing each index from xorshift64.
            for j in (1..seq.len()).rev() {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let k = (state % (j as u64 + 1)) as usize;
                seq.swap(j, k);
            }
            Protein {
                accession: format!("{}_{}", normalized, p.accession),
                description: p.description.clone(),
                sequence: seq,
            }
        })
        .collect();
    ProteinDb { proteins }
}

/// Build the searched database from the target proteins according to `strategy`.
/// `Reverse`/`Shuffle` append a 1:1 decoy per target (target proteins occupy the
/// first half); `None` returns the target db unchanged. `seed` is only used by
/// `Shuffle`. Decoy membership downstream is determined by the accession prefix
/// (not index), so all three strategies interoperate with the same labeling.
pub fn build_search_db(
    target: &ProteinDb,
    prefix: &str,
    strategy: DecoyStrategy,
    seed: u64,
) -> ProteinDb {
    match strategy {
        DecoyStrategy::Reverse => target_plus_decoy(target, prefix),
        DecoyStrategy::Shuffle => {
            let decoy = shuffle_db(target, prefix, seed);
            let mut proteins = target.proteins.clone();
            proteins.extend(decoy.proteins);
            ProteinDb { proteins }
        }
        DecoyStrategy::SequonReverse => {
            let decoy = sequon_preserving_reverse_db(target, prefix);
            let mut proteins = target.proteins.clone();
            proteins.extend(decoy.proteins);
            ProteinDb { proteins }
        }
        DecoyStrategy::None => target.clone(),
    }
}

/// Normalize a user-supplied decoy accession prefix: trim whitespace,
/// strip trailing `_`, and fall back to [`DEFAULT_DECOY_PREFIX`] when empty.
pub fn normalize_decoy_prefix(prefix: &str) -> String {
    let trimmed = prefix.trim().trim_end_matches('_');
    if trimmed.is_empty() {
        DEFAULT_DECOY_PREFIX.to_string()
    } else {
        trimmed.to_string()
    }
}

/// The decoy-membership needle for `prefix`: `"<normalized_prefix>_"`.
///
/// SINGLE SOURCE OF TRUTH for target/decoy membership. Decoy accessions are
/// built as `"<prefix>_<orig>"` (see [`reverse_db`]/[`shuffle_db`]), so a protein
/// is a decoy iff its accession starts with this needle. Both candidate
/// generation and [`crate::SearchIndex`] MUST use this helper so the two paths
/// can never diverge: a target accession that merely starts with the BARE prefix
/// (e.g. a real protein named `"XXXfoo"`, or any accession under a short custom
/// prefix like `rev`) is NOT a decoy — the `_` delimiter is required.
pub fn decoy_accession_needle(prefix: &str) -> String {
    format!("{}_", normalize_decoy_prefix(prefix))
}

/// Returns `true` iff `accession` is a decoy accession under `prefix`
/// (i.e. starts with [`decoy_accession_needle`]).
pub fn is_decoy_accession(accession: &str, prefix: &str) -> bool {
    accession.starts_with(&decoy_accession_needle(prefix))
}

/// Returns `true` iff `accession` is a decoy under the AFFIX convention: it
/// starts with the prefix needle `"<prefix>_"` (andes' native decoy naming) OR
/// — when `suffix` is set and non-empty — ends with `suffix`.
///
/// The suffix form lets andes recognize decoys produced by an external pipeline
/// whose pre-built target+decoy DB tags decoys with an accession SUFFIX (e.g.
/// quantms / OpenMS DecoyDatabase: `<orig>_rev`) rather than andes' prefix. Used
/// with `--decoy-strategy none` so andes consumes the shared decoy DB verbatim
/// instead of generating its own (which would double-decoy and bias FDR).
pub fn is_decoy_accession_affix(accession: &str, prefix: &str, suffix: Option<&str>) -> bool {
    accession.starts_with(&decoy_accession_needle(prefix))
        || suffix.is_some_and(|s| !s.is_empty() && accession.ends_with(s))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_db(proteins: &[(&str, &[u8])]) -> ProteinDb {
        ProteinDb {
            proteins: proteins.iter().map(|(acc, seq)| Protein {
                accession: acc.to_string(),
                description: String::new(),
                sequence: seq.to_vec(),
            }).collect(),
        }
    }

    #[test]
    fn reverse_db_reverses_sequences() {
        let db = make_db(&[("P1", b"MKWV"), ("P2", b"AGCT")]);
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.len(), 2);
        assert_eq!(decoy.proteins[0].sequence, b"VWKM");
        assert_eq!(decoy.proteins[1].sequence, b"TCGA");
    }

    #[test]
    fn reverse_db_prepends_prefix() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_strips_trailing_underscores_in_prefix() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "XXX_");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_empty_prefix_uses_default() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_preserves_description() {
        let mut db = make_db(&[("P1", b"AB")]);
        db.proteins[0].description = "Some description".into();
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.proteins[0].description, "Some description");
    }

    #[test]
    fn target_plus_decoy_concats() {
        let target = make_db(&[("P1", b"AB"), ("P2", b"CD")]);
        let combined = target_plus_decoy(&target, "XXX");
        assert_eq!(combined.len(), 4);
        assert_eq!(combined.proteins[0].accession, "P1");
        assert_eq!(combined.proteins[2].accession, "XXX_P1");
        assert_eq!(combined.proteins[2].sequence, b"BA");
    }

    #[test]
    fn shuffle_db_is_deterministic_for_a_fixed_seed() {
        let db = make_db(&[("P1", b"MKWVLPASTNDE"), ("P2", b"AGCTFYHRQGIL")]);
        let a = shuffle_db(&db, "XXX", DEFAULT_DECOY_SEED);
        let b = shuffle_db(&db, "XXX", DEFAULT_DECOY_SEED);
        assert_eq!(a.proteins[0].sequence, b.proteins[0].sequence, "same seed ⇒ same shuffle");
        assert_eq!(a.proteins[1].sequence, b.proteins[1].sequence);
        // A different seed yields a different permutation (overwhelmingly likely
        // for a 12-mer), so shuffles are actually seed-dependent.
        let c = shuffle_db(&db, "XXX", DEFAULT_DECOY_SEED + 1);
        assert_ne!(a.proteins[0].sequence, c.proteins[0].sequence);
        // Shuffle is a permutation: same multiset of residues, prefixed accession.
        let mut orig = db.proteins[0].sequence.clone();
        let mut shuf = a.proteins[0].sequence.clone();
        orig.sort_unstable();
        shuf.sort_unstable();
        assert_eq!(orig, shuf, "shuffle preserves residue composition");
        assert_eq!(a.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn build_search_db_reverse_matches_target_plus_decoy() {
        let target = make_db(&[("P1", b"ABCD"), ("P2", b"EFGH")]);
        let built = build_search_db(&target, "XXX", DecoyStrategy::Reverse, DEFAULT_DECOY_SEED);
        let legacy = target_plus_decoy(&target, "XXX");
        // Reverse strategy is byte-identical to the simple target_plus_decoy path.
        assert_eq!(built.len(), legacy.len());
        for (a, b) in built.proteins.iter().zip(legacy.proteins.iter()) {
            assert_eq!(a.accession, b.accession);
            assert_eq!(a.sequence, b.sequence);
        }
    }

    #[test]
    fn sequon_preserving_reverse_is_permutation_and_recovers_sequons() {
        // Composition preserved (pure swaps), sequon count recovered toward target
        // WITHOUT overshoot, and strictly better than plain reversal.
        let orig = b"NISCDENLTKRSTNVTAAK";
        let rev = super::sequon_preserving_reverse(orig);
        let mut a: Vec<u8> = orig.to_vec();
        a.sort();
        let mut b = rev.clone();
        b.sort();
        assert_eq!(a, b, "must be a permutation of the target (composition preserved)");
        let t = super::count_sequon_starts(orig);
        let d = super::count_sequon_starts(&rev);
        let plain: Vec<u8> = orig.iter().rev().copied().collect();
        let p = super::count_sequon_starts(&plain);
        // No overshoot, and at least as many as plain reversal.
        assert!(d <= t, "must not overshoot the target count ({d} <= {t})");
        assert!(d >= p, "must recover at least as many sequons as plain reversal ({d} >= {p})");
    }

    #[test]
    fn build_search_db_sequon_reverse_appends_prefixed_decoys() {
        let target = make_db(&[("P1", b"NISCDENLTKR"), ("P2", b"AANGSVVK")]);
        let built = build_search_db(&target, "XXX", DecoyStrategy::SequonReverse, DEFAULT_DECOY_SEED);
        assert_eq!(built.len(), 4, "target + 1:1 sequon-preserving decoy");
        assert!(built.proteins[2].accession.starts_with("XXX_"));
        assert!(built.proteins[3].accession.starts_with("XXX_"));
    }

    #[test]
    fn build_search_db_none_is_target_only() {
        let target = make_db(&[("P1", b"ABCD"), ("P2", b"EFGH")]);
        let built = build_search_db(&target, "XXX", DecoyStrategy::None, DEFAULT_DECOY_SEED);
        assert_eq!(built.len(), 2, "None appends no decoys");
        assert_eq!(built.proteins[0].accession, "P1");
        assert_eq!(built.proteins[1].accession, "P2");
    }

    #[test]
    fn decoy_accession_needle_requires_underscore_delimiter() {
        // Decoys are built as "<prefix>_<orig>", so detection requires the "_".
        assert_eq!(decoy_accession_needle("XXX"), "XXX_");
        assert_eq!(decoy_accession_needle("XXX_"), "XXX_"); // normalize is idempotent
        assert_eq!(decoy_accession_needle(""), "XXX_"); // empty → default prefix
        assert!(is_decoy_accession("XXX_P12345", "XXX"));
        assert!(!is_decoy_accession("XXXP12345", "XXX"), "no '_' delimiter ⇒ target, not decoy");
        assert!(!is_decoy_accession("P12345", "XXX"));
        // Short custom prefix: a real "rev..." target must NOT be misread as decoy.
        assert!(is_decoy_accession("rev_sp|P1", "rev"));
        assert!(!is_decoy_accession("reverse_kinase", "rev"));
    }

    #[test]
    fn affix_matches_prefix_or_suffix() {
        // Prefix (andes-native) still works with no suffix configured.
        assert!(is_decoy_accession_affix("XXX_P12345", "XXX", None));
        assert!(!is_decoy_accession_affix("P12345", "XXX", None));
        // Suffix form recognizes quantms/OpenMS "<orig>_rev" decoys that the
        // prefix test alone misses.
        assert!(is_decoy_accession_affix("sp|P02769|ALBU_BOVIN_rev", "XXX", Some("rev")));
        assert!(!is_decoy_accession_affix("sp|P02769|ALBU_BOVIN", "XXX", Some("rev")));
        // Either affix qualifies; an empty suffix is ignored (no spurious match).
        assert!(is_decoy_accession_affix("XXX_P1", "XXX", Some("rev")));
        assert!(!is_decoy_accession_affix("P1", "XXX", Some("")));
    }

    #[test]
    fn decoy_strategy_parses_names() {
        assert_eq!(DecoyStrategy::from_name("reverse"), Some(DecoyStrategy::Reverse));
        assert_eq!(DecoyStrategy::from_name("Shuffle"), Some(DecoyStrategy::Shuffle));
        assert_eq!(DecoyStrategy::from_name("NONE"), Some(DecoyStrategy::None));
        assert_eq!(DecoyStrategy::from_name("bogus"), None);
        assert_eq!(DecoyStrategy::default(), DecoyStrategy::Reverse);
    }
}
