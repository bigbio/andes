package edu.ucsd.msjava.msgf;

/**
 * Streaming merger for PrimitiveGeneratingFunction score distributions
 * across isotope mass indices. Callers feed each GF via {@link #accept}
 * after constructing it; the group computes the GF, merges its
 * {@link ScoreDist} into a running aggregate, and releases the reference.
 * Peak memory is therefore one graph + one GF at a time, independent of
 * the number of mass indices.
 *
 * Math is identical to the previous register-all-then-merge approach
 * because ScoreDist.addProbDist with scoreDiff=0 and aaProb=1f is a
 * linear sum over the probability arrays.
 */
public class PrimitiveGeneratingFunctionGroup {
    private int minScore = Integer.MAX_VALUE;
    private int maxScore = Integer.MIN_VALUE;
    private ScoreDist mergedScoreDist = null;

    /**
     * Compute the supplied GF if needed and merge its distribution into
     * the running aggregate. The caller must drop its own reference to
     * {@code gf} after this call to allow its {@code distByNode} and
     * graph to be collected before the next mass index is built.
     */
    public void accept(PrimitiveGeneratingFunction gf) {
        if (!gf.isGFComputed()) {
            if (!gf.computeGeneratingFunction()) return;
        }
        ScoreDist dist = gf.getScoreDist();
        if (dist == null) return;

        int gfMin = gf.getMinScore();
        int gfMax = gf.getMaxScore();

        if (mergedScoreDist == null) {
            minScore = gfMin;
            maxScore = gfMax;
            mergedScoreDist = new ScoreDist(minScore, maxScore, false, true);
            mergedScoreDist.addProbDist(dist, 0, 1f);
            return;
        }

        int newMin = Math.min(minScore, gfMin);
        int newMax = Math.max(maxScore, gfMax);
        if (newMin != minScore || newMax != maxScore) {
            ScoreDist expanded = new ScoreDist(newMin, newMax, false, true);
            expanded.addProbDist(mergedScoreDist, 0, 1f);
            mergedScoreDist = expanded;
            minScore = newMin;
            maxScore = newMax;
        }
        mergedScoreDist.addProbDist(dist, 0, 1f);
    }

    public boolean isComputed() { return mergedScoreDist != null; }

    public double getSpectralProbability(int score) {
        return mergedScoreDist.getSpectralProbability(score);
    }

    public int getMaxScore() { return mergedScoreDist.getMaxScore(); }
    public ScoreDist getScoreDist() { return mergedScoreDist; }
}
