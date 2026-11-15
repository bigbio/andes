//! Crosslinking (XL-MS) search building blocks for andes, starting with the
//! MS-cleavable linker DSSO.
//!
//! This crate holds the pure, dependency-free parts of the crosslink search
//! (linker chemistry, the reporter-doublet detector, link-site rules, pair mass
//! lookup, decoy classes and pair FDR). Peptide scoring lives in
//! `search::xl_search`, which scores every peptide with the andes core scorer.
//!
//! # Design: how a crosslink maps onto the andes core scorer
//!
//! A crosslink-spectrum match (CSM) explains one MS2 spectrum as two peptides,
//! alpha and beta, joined through one residue each by the linker. Its neutral
//! mass is `M(alpha) + M(beta) + M(linker, intact)`.
//!
//! **Each peptide is scored alone, as a modified linear peptide.** The partner
//! (or what the cleaved linker leaves behind) is a mass on the link residue.
//! The peptide carries that mass as a modification and is scored with
//! `score_psm + psm_edge_score` (the rank-LLR node score plus the edge score)
//! on the spectrum's `ScoredSpectrum`, built once per spectrum with the run's
//! `RankScorer` at the precursor charge. There is no crosslink-specific
//! spectrum scorer: model selection, fragment tolerance, precursor-peak
//! filtering and partition choice are those of the standard search.
//!
//! The link-residue masses ([`crosslinker::LinkForm`]) cover what a DSSO
//! fragment carries:
//!
//! | form | mass on the link residue | which fragments |
//! |---|---|---|
//! | alkene stub | 54.0106 | peptide released by linker cleavage, and its b/y ions |
//! | thiol stub | 85.9826 | as alkene, the other half (sulfenic acid minus water) |
//! | sulfenic stub | 103.9932 | as thiol, before water loss |
//! | intact partner | partner peptide + 158.0038 | b/y ions of the uncleaved crosslink |
//!
//! Every allowed site ([`sites::link_sites`]: an unmodified Lys that is not the
//! peptide's C-terminal residue unless it is the protein's, or the protein
//! N-terminus) and every form is scored; the peptide's score is the best over
//! (site, form), and the per-form bests are kept as features.
//!
//! **Pair score.** The CSM score is `min(score_alpha, score_beta)`, so a strong
//! alpha cannot carry a junk beta. Both peptide scores, their sum and the
//! label-blind score gaps to the next distinct sequence travel with the CSM as
//! features.
//!
//! **Candidates from doublets.** DSSO cleaves in MS2 into the alkene and thiol
//! forms of each peptide, 31.9721 Da apart ([`doublet`]). A doublet fixes one
//! peptide's mass; the precursor fixes the other's as
//! `M(precursor) - M(peptide) - M(linker)` (isotope offsets 0..=1). A pair is
//! *full* when both masses carry a doublet, *half* when only one does and the
//! complement shows at least one stub ion. Peptides are looked up by mass in a
//! mass-sorted table of linkable candidates ([`pair::MassIndex`]).
//!
//! **Decoys.** Target and decoy peptides come from the same andes index (the
//! search's decoy strategy), are looked up and scored identically, and the
//! CSM class is TT, TD or DD ([`decoy::PairClass`]). Pair FDR is
//! `(TD - DD) / TT` ([`fdr`]), estimated separately for intra- and
//! inter-protein links, with the link class assigned by the same rule to
//! targets and decoys. CSMs roll up to residue pairs and protein pairs by the
//! best score of each group, and the same class-split estimate is made at each
//! level ([`rollup`]).
//!
//! **Rescoring.** Every CSM goes into a pair-aware PIN. Inter-links are
//! scored by the pair score: on the PXD042173 truth, Percolator models over
//! all CSMs or over the inter rows gave inter-link q-values below the actual
//! false share, while the pair score held it. With `--rescore` the intra-links
//! are scored by a Percolator model trained on the intra rows only (their own
//! PIN); the FDR stays split by link class, so the two score scales are never
//! compared.
//!
//! **Reported levels.** Intra-links are reported up to residue pairs (with
//! a q-value also among the residue pairs with at least two CSMs); inter-links
//! at CSM, residue-pair and protein-pair level. Intra protein pairs (a
//! protein with itself) stay in the protein-pair table without a q-value.
//!
//! **Correctness check.** [`pair::pairs_by_sum`] is the exhaustive pair
//! enumerator by precursor mass; it is tested against a brute-force all-pairs
//! loop, and the search driver can run it as an oracle (no doublet gate) on a
//! slice of spectra to measure what the doublet gate loses.

pub mod crosslinker;
pub mod decoy;
pub mod doublet;
pub mod fdr;
pub mod match_prob;
pub mod pair;
pub mod rollup;
pub mod sites;
