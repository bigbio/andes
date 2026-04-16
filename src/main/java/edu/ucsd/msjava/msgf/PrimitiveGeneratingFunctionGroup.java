package edu.ucsd.msjava.msgf;

/**
 * Groups multiple PrimitiveGeneratingFunction instances (one per isotope mass index)
 * and merges their score distributions.
 * Replaces GeneratingFunctionGroup<NominalMass> in the DB search hot path.
 */
public class PrimitiveGeneratingFunctionGroup {
    private PrimitiveGeneratingFunction[] gfs;
    private int count;
    private ScoreDist mergedScoreDist;

    public PrimitiveGeneratingFunctionGroup(int capacity) {
        this.gfs = new PrimitiveGeneratingFunction[capacity];
        this.count = 0;
    }

    public void register(PrimitiveGeneratingFunction gf) {
        if (count >= gfs.length) {
            PrimitiveGeneratingFunction[] newArr = new PrimitiveGeneratingFunction[gfs.length * 2];
            System.arraycopy(gfs, 0, newArr, 0, count);
            gfs = newArr;
        }
        gfs[count++] = gf;
    }

    public boolean computeGeneratingFunction() {
        int minScore = Integer.MAX_VALUE;
        int maxScore = Integer.MIN_VALUE;

        // Mirrors legacy GeneratingFunctionGroup: bounds are collected only for
        // GFs computed in this call. Safe under current usage in
        // DBScanner.computeSpecEValue, which constructs a fresh group + fresh
        // GFs per spectrum so isGFComputed() is always false here. If GFs are
        // ever cached/reused across calls, this loop must include their bounds
        // too or mergedScoreDist will be sized wrong.
        for (int i = 0; i < count; i++) {
            PrimitiveGeneratingFunction gf = gfs[i];
            if (!gf.isGFComputed() && gf.computeGeneratingFunction()) {
                if (gf.getMinScore() < minScore) minScore = gf.getMinScore();
                if (gf.getMaxScore() > maxScore) maxScore = gf.getMaxScore();
            }
        }

        if (minScore >= maxScore) return false;

        mergedScoreDist = new ScoreDist(minScore, maxScore, false, true);
        for (int i = 0; i < count; i++) {
            ScoreDist dist = gfs[i].getScoreDist();
            if (dist != null) {
                mergedScoreDist.addProbDist(dist, 0, 1f);
            }
        }
        return true;
    }

    public double getSpectralProbability(int score) {
        return mergedScoreDist.getSpectralProbability(score);
    }

    public int getMaxScore() {
        return mergedScoreDist.getMaxScore();
    }

    public ScoreDist getScoreDist() {
        return mergedScoreDist;
    }
}
