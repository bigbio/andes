package edu.ucsd.msjava.msscorer;

/**
 * Package-internal helper that exposes {@link NewRankScorer#getPartition(int, float, int)}
 * and {@link NewRankScorer#getNumSegments()} to callers outside the {@code msscorer}
 * package without broadening the public surface of {@link NewRankScorer}.
 *
 * <p>Used by {@code edu.ucsd.msjava.fragindex.FragmentIndexCandidateGenerator} to
 * resolve the {@link Partition} corresponding to a spectrum's (charge, parentMass)
 * without constructing a {@link NewScoredSpectrum} (which would mutate the spectrum
 * via {@code filterPrecursorPeaks}).
 *
 * <p>Introduced in the speed-rewrite-v2 work. Stateless; safe for concurrent use.
 */
public final class ScorerPartitionAccess {

    private ScorerPartitionAccess() {}

    /**
     * Resolves the {@link Partition} used by {@link NewRankScorer#getNodeScore(Partition, edu.ucsd.msjava.msutil.IonType, int)}
     * for a spectrum at (charge, parentMass, lastSegment).
     */
    public static Partition lastSegmentPartition(NewRankScorer scorer, int charge, float parentMass) {
        int lastSeg = scorer.getNumSegments() - 1;
        return scorer.getPartition(charge, parentMass, lastSeg);
    }
}
