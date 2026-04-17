package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Modification.Location;

import java.util.ArrayList;

/**
 * Primitive-array–based amino acid graph for the generating function.
 * Replaces FlexAminoAcidGraph in the DB search hot path to eliminate
 * HashMap/ArrayList/NominalMass object overhead.
 *
 * Graph topology is stored in CSR (Compressed Sparse Row) format:
 *   edgeOffset[node+1] - edgeOffset[node] = number of incoming edges for node
 *   edgePrevNode[e], edgeProb[e], edgeMass[e], edgeScore[e] = edge data
 *
 * Node scores are stored in a flat int[] indexed by nominal mass.
 */
public class PrimitiveAminoAcidGraph {
    private final int peptideMass;
    private final AminoAcidSet aaSet;
    private final Enzyme enzyme;
    private final boolean direction;
    private final int minNodeMass;
    private final int massOffset;

    private int nodeCount;
    private int[] activeNodes;
    private int[] massToNodeIdx;

    private int totalEdges;
    private int[] edgeOffset;
    private int[] edgePrevNode;
    private float[] edgeProb;
    private float[] edgeMass;
    private int[] edgeScore;

    private int[] nodeScores;

    private int sourceNodeIdx;
    private int sinkNodeIdx;

    public PrimitiveAminoAcidGraph(
            AminoAcidSet aaSet,
            int peptideMass,
            Enzyme enzyme,
            ScoredSpectrum<NominalMass> scoredSpec,
            boolean useProteinNTerm,
            boolean useProteinCTerm
    ) {
        this.aaSet = aaSet;
        this.peptideMass = peptideMass;
        this.enzyme = enzyme;
        this.direction = scoredSpec.getMainIonDirection();

        Location sourceLocation;
        if (direction) {
            sourceLocation = useProteinNTerm ? Location.Protein_N_Term : Location.N_Term;
        } else {
            sourceLocation = useProteinCTerm ? Location.Protein_C_Term : Location.C_Term;
        }

        Location sinkLocation;
        if (direction) {
            sinkLocation = useProteinCTerm ? Location.Protein_C_Term : Location.C_Term;
        } else {
            sinkLocation = useProteinNTerm ? Location.Protein_N_Term : Location.N_Term;
        }

        ArrayList<AminoAcid> sourceAAs = aaSet.getAAList(sourceLocation);
        ArrayList<AminoAcid> anywhereAAs = aaSet.getAAList(Location.Anywhere);
        ArrayList<AminoAcid> sinkAAs = aaSet.getAAList(sinkLocation);

        int minMass = 0;
        for (AminoAcid aa : sourceAAs) {
            minMass = Math.min(minMass, aa.getNominalMass());
        }
        for (AminoAcid aa : anywhereAAs) {
            minMass = Math.min(minMass, 1 + aa.getNominalMass());
        }
        for (AminoAcid aa : sinkAAs) {
            minMass = Math.min(minMass, peptideMass - aa.getNominalMass());
        }
        this.minNodeMass = minMass;
        this.massOffset = -minMass;

        boolean[] reachable = new boolean[peptideMass - minNodeMass + 1];
        reachable[toDenseIndex(0)] = true;

        boolean addCleavageFromSource = enzyme != null && direction == enzyme.isNTerm();

        // Phase 1: discover reachable masses and count incoming edges per target mass.
        int[] inEdgeCountByMass = new int[peptideMass - minNodeMass + 1];

        // Forward edges from source (mass 0)
        for (AminoAcid aa : sourceAAs) {
            int nextMass = aa.getNominalMass();
            if (nextMass >= peptideMass || !isRepresentableMass(nextMass)) continue;
            reachable[toDenseIndex(nextMass)] = true;
            inEdgeCountByMass[toDenseIndex(nextMass)]++;
        }

        // Forward edges from intermediate nodes
        for (int curMass = 1; curMass < peptideMass; curMass++) {
            if (!reachable[toDenseIndex(curMass)]) continue;
            for (AminoAcid aa : anywhereAAs) {
                int nextMass = curMass + aa.getNominalMass();
                if (nextMass >= peptideMass || !isRepresentableMass(nextMass)) continue;
                reachable[toDenseIndex(nextMass)] = true;
                inEdgeCountByMass[toDenseIndex(nextMass)]++;
            }
        }

        // Backward edges to sink (peptideMass)
        boolean addCleavageToSink = enzyme != null && direction != enzyme.isNTerm();
        for (AminoAcid aa : sinkAAs) {
            int prevMass = peptideMass - aa.getNominalMass();
            if (!isRepresentableMass(prevMass) || !reachable[toDenseIndex(prevMass)]) continue;
            inEdgeCountByMass[toDenseIndex(peptideMass)]++;
        }
        reachable[toDenseIndex(peptideMass)] = true;

        // Phase 2: Count active nodes and build node index
        int count = 0;
        for (int m = minNodeMass; m <= peptideMass; m++) {
            if (reachable[toDenseIndex(m)]) count++;
        }
        this.nodeCount = count;
        this.activeNodes = new int[nodeCount];
        this.massToNodeIdx = new int[peptideMass - minNodeMass + 1];
        java.util.Arrays.fill(massToNodeIdx, -1);
        int idx = 0;
        activeNodes[idx] = 0;
        massToNodeIdx[toDenseIndex(0)] = idx;
        this.sourceNodeIdx = idx;
        idx++;
        for (int m = minNodeMass; m <= peptideMass; m++) {
            if (m == 0 || !reachable[toDenseIndex(m)]) {
                continue;
            }
            activeNodes[idx] = m;
            massToNodeIdx[toDenseIndex(m)] = idx;
            idx++;
        }
        this.sinkNodeIdx = getNodeIndexForMass(peptideMass);

        // Phase 3: Build CSR offsets from per-mass incoming edge counts.
        this.edgeOffset = new int[nodeCount + 1];
        for (int ni = 0; ni < nodeCount; ni++) {
            int mass = activeNodes[ni];
            edgeOffset[ni + 1] = edgeOffset[ni] + inEdgeCountByMass[toDenseIndex(mass)];
        }
        this.totalEdges = edgeOffset[nodeCount];

        this.edgePrevNode = new int[totalEdges];
        this.edgeProb = new float[totalEdges];
        this.edgeMass = new float[totalEdges];
        this.edgeScore = new int[totalEdges];

        // Phase 4: Fill CSR edges directly (same generation order as before).
        int[] writeCursor = java.util.Arrays.copyOf(edgeOffset, nodeCount);

        for (AminoAcid aa : sourceAAs) {
            int nextMass = aa.getNominalMass();
            if (nextMass >= peptideMass || !isRepresentableMass(nextMass)) continue;
            int cleavageScore = 0;
            if (addCleavageFromSource) {
                cleavageScore = enzyme.isCleavable(aa) ? aaSet.getPeptideCleavageCredit() : aaSet.getPeptideCleavagePenalty();
            }
            writeEdge(nextMass, 0, aa.getProbability(), aa.getMass(), cleavageScore, writeCursor);
        }

        for (int curMass = 1; curMass < peptideMass; curMass++) {
            if (!reachable[toDenseIndex(curMass)]) continue;
            for (AminoAcid aa : anywhereAAs) {
                int nextMass = curMass + aa.getNominalMass();
                if (nextMass >= peptideMass || !isRepresentableMass(nextMass)) continue;
                writeEdge(nextMass, curMass, aa.getProbability(), aa.getMass(), 0, writeCursor);
            }
        }

        for (AminoAcid aa : sinkAAs) {
            int prevMass = peptideMass - aa.getNominalMass();
            if (!isRepresentableMass(prevMass) || !reachable[toDenseIndex(prevMass)]) continue;
            int cleavageScore = 0;
            if (addCleavageToSink) {
                cleavageScore = enzyme.isCleavable(aa) ? aaSet.getPeptideCleavageCredit() : aaSet.getPeptideCleavagePenalty();
            }
            writeEdge(peptideMass, prevMass, aa.getProbability(), aa.getMass(), cleavageScore, writeCursor);
        }

        // Phase 5: Compute edge error scores and node scores.
        computeEdgeErrorScores(scoredSpec);
        this.edgeMass = null; // no longer needed after error scores computed
        computeNodeScores(scoredSpec);
    }

    private void writeEdge(int targetMass, int prevMass, float prob, float mass, int cleavageScore, int[] writeCursor) {
        int targetNodeIdx = getNodeIndexForMass(targetMass);
        if (targetNodeIdx < 0) {
            return;
        }
        int edgeIdx = writeCursor[targetNodeIdx]++;
        edgePrevNode[edgeIdx] = prevMass;
        edgeScore[edgeIdx] = cleavageScore;
        edgeProb[edgeIdx] = prob;
        edgeMass[edgeIdx] = mass;
    }

    private void computeEdgeErrorScores(ScoredSpectrum<NominalMass> scoredSpec) {
        // Cache one NominalMass per active node so per-edge prev-node lookup
        // is O(1) instead of allocating a fresh NominalMass on every edge.
        NominalMass[] nmByNode = new NominalMass[nodeCount];
        for (int ni = 0; ni < nodeCount; ni++) {
            nmByNode[ni] = new NominalMass(activeNodes[ni]);
        }

        for (int ni = 0; ni < nodeCount; ni++) {
            int curMass = activeNodes[ni];
            if (curMass == 0 || curMass == peptideMass) continue;

            NominalMass curNM = nmByNode[ni];
            for (int e = edgeOffset[ni]; e < edgeOffset[ni + 1]; e++) {
                int prevMass = edgePrevNode[e];
                int prevNodeIdx = getNodeIndexForMass(prevMass);
                NominalMass prevNM = (prevNodeIdx >= 0)
                        ? nmByNode[prevNodeIdx]
                        : new NominalMass(prevMass);
                int errorScore = scoredSpec.getEdgeScore(curNM, prevNM, edgeMass[e]);
                if (errorScore < -100 || errorScore > 100) {
                    errorScore = -4;
                }
                edgeScore[e] += errorScore;
            }
        }
    }

    private void computeNodeScores(ScoredSpectrum<NominalMass> scoredSpec) {
        this.nodeScores = new int[nodeCount];

        for (int ni = 1; ni < nodeCount; ni++) {
            int mass = activeNodes[ni];
            if (mass == peptideMass) {
                nodeScores[ni] = 0;
                continue;
            }
            int compMass = peptideMass - mass;
            NominalMass nodeNM = new NominalMass(mass);
            NominalMass compNM = new NominalMass(compMass);
            if (!direction) {
                nodeScores[ni] = scoredSpec.getNodeScore(compNM, nodeNM);
            } else {
                nodeScores[ni] = scoredSpec.getNodeScore(nodeNM, compNM);
            }
        }
    }

    // Accessors
    public int getPeptideMass() { return peptideMass; }
    public int getNodeCount() { return nodeCount; }
    public int[] getActiveNodes() { return activeNodes; }
    public int[] getMassToNodeIdx() { return massToNodeIdx; }
    public int getMassOffset() { return massOffset; }
    public int getSourceNodeIdx() { return sourceNodeIdx; }
    public int getSinkNodeIdx() { return sinkNodeIdx; }
    public int getTotalEdges() { return totalEdges; }
    public int[] getEdgeOffset() { return edgeOffset; }
    public int[] getEdgePrevNode() { return edgePrevNode; }
    public float[] getEdgeProb() { return edgeProb; }
    public int[] getEdgeScore() { return edgeScore; }
    public int getNodeScore(int nodeIdx) { return nodeScores[nodeIdx]; }
    public int[] getNodeScores() { return nodeScores; }
    public AminoAcidSet getAASet() { return aaSet; }
    public Enzyme getEnzyme() { return enzyme; }

    public int getNodeIndexForMass(int mass) {
        if (!isRepresentableMass(mass)) {
            return -1;
        }
        return massToNodeIdx[toDenseIndex(mass)];
    }

    private int toDenseIndex(int mass) {
        return mass + massOffset;
    }

    private boolean isRepresentableMass(int mass) {
        return mass >= minNodeMass && mass <= peptideMass;
    }
}
