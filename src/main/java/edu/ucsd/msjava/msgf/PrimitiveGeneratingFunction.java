package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;

/**
 * Primitive-array–based generating function for computing spectral E-values.
 * Replaces GeneratingFunction<NominalMass> in the DB search hot path.
 *
 * All HashMaps are replaced with int[]/double[] arrays indexed by node index.
 * The inner DP loop operates on contiguous memory with zero object allocation.
 */
public class PrimitiveGeneratingFunction {
    // ---- GF score-distribution trace instrumentation (gated by system properties) ----
    // Mirrors the per-ion trace pattern in NewScoredSpectrum / DBScanner. JIT folds
    // these static-final booleans to a no-op when GFTRACE is false, so this is free
    // in production.
    public static final boolean GFTRACE = "true".equals(System.getProperty("msgfplus.gftrace"));
    public static final int GFTRACE_SCAN = Integer.parseInt(System.getProperty("msgfplus.trace.scan", "-1"));
    public static final String GFTRACE_PEP = System.getProperty("msgfplus.trace.pep", "");

    private final PrimitiveAminoAcidGraph graph;

    private ScoreDist distribution = null;
    private boolean isGFComputed = false;

    private int[] minScoreByNode;

    /**
     * Per-node score distributions retained only while {@link #GFTRACE} is enabled.
     * In production this stays null and the local {@code distByNode} inside
     * {@link #computeGeneratingFunction()} is GC'd as soon as that method returns.
     */
    private ScoreDist[] traceDistByNode = null;

    public PrimitiveGeneratingFunction(PrimitiveAminoAcidGraph graph) {
        this.graph = graph;
    }

    public boolean isGFComputed() { return isGFComputed; }
    public ScoreDist getScoreDist() { return distribution; }

    public int getMinScore() { return distribution.getMinScore(); }
    public int getMaxScore() { return distribution.getMaxScore(); }

    public double getSpectralProbability(int score) {
        if (distribution == null || !distribution.isProbSet()) return 1.0;
        return distribution.getSpectralProbability(score);
    }

    /** Graph accessor for trace consumers; package-private alternatives are awkward across packages. */
    public PrimitiveAminoAcidGraph getGraph() { return graph; }

    /** Retained only when {@link #GFTRACE} is true; null otherwise. */
    public ScoreDist[] getTraceDistByNode() { return traceDistByNode; }

    /**
     * Dump GF_NODE / GF_PROB / GF_TAIL lines for this single GF to {@code System.err},
     * matching the format emitted by Rust {@code msgf-trace --print-score-dist}.
     * No-op unless {@link #GFTRACE} is true and {@code traceDistByNode} was retained
     * during {@link #computeGeneratingFunction()}.
     *
     * @param scan         scan number label (for filtering downstream)
     * @param pepSeq       peptide label (for filtering downstream)
     * @param matchedScore matched PSM score used for the tail-sum / spec_prob line
     */
    public void dumpScoreDistTrace(int scan, String pepSeq, int matchedScore) {
        if (!GFTRACE || traceDistByNode == null || distribution == null) return;
        int[] activeNodes = graph.getActiveNodes();
        for (int ni = 0; ni < traceDistByNode.length; ni++) {
            ScoreDist d = traceDistByNode[ni];
            if (d == null) continue;
            int nodeMass = activeNodes[ni];
            System.err.printf("GF_NODE: scan=%d pep=%s node_idx=%d mass=%d min_score=%d max_score=%d%n",
                    scan, pepSeq, ni, nodeMass, d.getMinScore(), d.getMaxScore());
            // ScoreDist max-score is exclusive (probDistribution.length == maxScore - minScore).
            for (int s = d.getMinScore(); s < d.getMaxScore(); s++) {
                double p = d.getProbability(s);
                if (p == 0.0) continue;
                System.err.printf("GF_PROB: scan=%d pep=%s node_idx=%d score=%d prob=%.6e%n",
                        scan, pepSeq, ni, s, p);
            }
        }
        // Tail sum on this single GF's final distribution. Rust dumps over [matched_score, final_max).
        double tail = 0.0;
        int finalMin = distribution.getMinScore();
        int finalMax = distribution.getMaxScore();
        for (int s = Math.max(matchedScore, finalMin); s < finalMax; s++) {
            tail += distribution.getProbability(s);
        }
        double sp = distribution.getSpectralProbability(matchedScore);
        System.err.printf(
                "GF_TAIL: scan=%d pep=%s matched_score=%d spec_prob=%.6e tail_sum=%.6e final_min=%d final_max=%d%n",
                scan, pepSeq, matchedScore, sp, tail, finalMin, finalMax);
    }

    public void setUpScoreThreshold(int score) {
        int nodeCount = graph.getNodeCount();
        int[] activeNodes = graph.getActiveNodes();
        int[] edgeOffset = graph.getEdgeOffset();
        int[] edgePrevNode = graph.getEdgePrevNode();
        int[] edgeScoreArr = graph.getEdgeScore();
        int[] nodeScoresArr = graph.getNodeScores();
        int peptideMass = graph.getPeptideMass();
        int sourceIdx = graph.getSourceNodeIdx();

        int adjustedScore = score;
        Enzyme enzyme = graph.getEnzyme();
        if (enzyme != null) {
            adjustedScore -= graph.getAASet().getNeighboringAACleavageCredit();
        }

        minScoreByNode = new int[nodeCount];
        java.util.Arrays.fill(minScoreByNode, Integer.MAX_VALUE);

        int sinkIdx = graph.getSinkNodeIdx();
        minScoreByNode[sinkIdx] = adjustedScore;

        for (int e = edgeOffset[sinkIdx]; e < edgeOffset[sinkIdx + 1]; e++) {
            int prevMass = edgePrevNode[e];
            int prevIdx = graph.getNodeIndexForMass(prevMass);
            if (prevIdx < 0) continue;
            int newMin = adjustedScore - edgeScoreArr[e];
            if (newMin < minScoreByNode[prevIdx]) {
                minScoreByNode[prevIdx] = newMin;
            }
        }

        for (int ni = nodeCount - 1; ni >= 0; ni--) {
            if (ni == sourceIdx || ni == sinkIdx) {
                continue;
            }
            if (minScoreByNode[ni] == Integer.MAX_VALUE) continue;
            int curMass = activeNodes[ni];
            if (curMass == peptideMass) continue;
            int curNodeScore = nodeScoresArr[ni];

            for (int e = edgeOffset[ni]; e < edgeOffset[ni + 1]; e++) {
                int prevMass = edgePrevNode[e];
                int prevIdx = graph.getNodeIndexForMass(prevMass);
                if (prevIdx < 0) continue;
                int newMin = minScoreByNode[ni] - (curNodeScore + edgeScoreArr[e]);
                if (newMin < minScoreByNode[prevIdx]) {
                    minScoreByNode[prevIdx] = newMin;
                }
            }
        }
    }

    public boolean computeGeneratingFunction() {
        int nodeCount = graph.getNodeCount();
        int[] edgeOffset = graph.getEdgeOffset();
        int[] edgePrevNode = graph.getEdgePrevNode();
        float[] edgeProb = graph.getEdgeProb();
        int[] edgeScoreArr = graph.getEdgeScore();
        int[] nodeScoresArr = graph.getNodeScores();
        int sourceIdx = graph.getSourceNodeIdx();
        int sinkIdx = graph.getSinkNodeIdx();

        ScoreDist[] distByNode = new ScoreDist[nodeCount];

        ScoreDist sourceDist = new ScoreDist(0, 1, false, true);
        sourceDist.setProb(0, 1.0);
        distByNode[sourceIdx] = sourceDist;

        // Scratch buffer for valid edges.
        int maxEdgesPerNode = 0;
        for (int ni = 0; ni < nodeCount; ni++) {
            int count = edgeOffset[ni + 1] - edgeOffset[ni];
            if (count > maxEdgesPerNode) maxEdgesPerNode = count;
        }
        int[] validEdges = new int[maxEdgesPerNode];

        // DP over intermediate nodes (skip the explicit source node)
        for (int ni = 0; ni < nodeCount; ni++) {
            if (ni == sourceIdx) {
                continue;
            }
            int curNodeScore = nodeScoresArr[ni];

            if (minScoreByNode != null && minScoreByNode[ni] == Integer.MAX_VALUE) {
                continue;
            }

            int curMinScore;
            if (minScoreByNode != null) {
                curMinScore = minScoreByNode[ni];
            } else {
                curMinScore = Integer.MAX_VALUE;
            }
            int curMaxScore = Integer.MIN_VALUE;

            int validCount = 0;
            for (int e = edgeOffset[ni]; e < edgeOffset[ni + 1]; e++) {
                int prevMass = edgePrevNode[e];
                int prevIdx = graph.getNodeIndexForMass(prevMass);
                if (prevIdx < 0) continue;
                ScoreDist prevDist = distByNode[prevIdx];
                if (prevDist == null) continue;

                int combinedScore = curNodeScore + edgeScoreArr[e];
                int possibleMax = prevDist.getMaxScore() + combinedScore;
                if (possibleMax > curMaxScore) curMaxScore = possibleMax;

                if (minScoreByNode == null) {
                    int possibleMin = prevDist.getMinScore() + combinedScore;
                    if (possibleMin < curMinScore) curMinScore = possibleMin;
                }

                validEdges[validCount++] = e;
            }

            if (curMinScore >= curMaxScore || validCount == 0) {
                continue;
            }

            if (curMinScore < -10000 || curMaxScore > 10000) {
                continue;
            }

            ScoreDist curDist = new ScoreDist(curMinScore, curMaxScore, false, true);

            for (int vi = 0; vi < validCount; vi++) {
                int e = validEdges[vi];
                int prevMass = edgePrevNode[e];
                int prevIdx = graph.getNodeIndexForMass(prevMass);
                ScoreDist prevDist = distByNode[prevIdx];
                int combinedScore = curNodeScore + edgeScoreArr[e];
                curDist.addProbDist(prevDist, combinedScore, edgeProb[e]);
            }

            if (curDist.getProbability(curDist.getMaxScore() - 1) == 0) {
                curDist.setProb(curDist.getMaxScore() - 1, Float.MIN_VALUE);
            }

            distByNode[ni] = curDist;
        }

        // Process sink node — merge into final distribution
        ScoreDist sinkDist = distByNode[sinkIdx];
        if (sinkDist == null) return false;

        int minScore = sinkDist.getMinScore();
        int maxScore = sinkDist.getMaxScore();

        if (maxScore <= minScore) return false;

        // Apply neighboring AA adjustment
        Enzyme enzyme = graph.getEnzyme();
        AminoAcidSet aaSetLocal = graph.getAASet();
        ScoreDist finalDist;

        if (enzyme != null && enzyme.getResidues() != null) {
            int credit = aaSetLocal.getNeighboringAACleavageCredit();
            int penalty = aaSetLocal.getNeighboringAACleavagePenalty();
            finalDist = new ScoreDist(minScore + penalty, maxScore + credit, false, true);
            finalDist.addProbDist(sinkDist, credit, aaSetLocal.getProbCleavageSites());
            finalDist.addProbDist(sinkDist, penalty, 1 - aaSetLocal.getProbCleavageSites());
        } else {
            finalDist = sinkDist;
        }

        this.distribution = finalDist;
        this.isGFComputed = true;
        // Retain per-node distributions only when score-distribution trace is enabled
        // (system property -Dmsgfplus.gftrace=true). In production GFTRACE is false,
        // distByNode falls out of scope, and memory behavior matches pre-trace baseline.
        if (GFTRACE) {
            this.traceDistByNode = distByNode;
        }
        return true;
    }
}
