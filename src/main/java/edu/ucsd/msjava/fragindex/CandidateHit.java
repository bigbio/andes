package edu.ucsd.msjava.fragindex;

/**
 * A single Tier-1 survivor from {@link FragmentIndexCandidateGenerator}.
 * {@code slabId} + {@code localPeptideId} together identify the peptide
 * in the {@link FragmentIndex}: use
 * {@code index.peptideTable(slabId).sequence(localPeptideId)} to resolve
 * back to the peptide string and
 * {@code index.peptideTable(slabId).precursorMass(localPeptideId)} for the
 * neutral monoisotopic mass.
 *
 * <p>{@code newRankSum} is the accumulated {@code NewRankScorer} log-score
 * across matched fragment buckets for this peptide on this spectrum — the
 * ranking criterion used during top-K extraction.
 */
public record CandidateHit(int slabId, int localPeptideId, float newRankSum) {}
