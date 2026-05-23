package edu.ucsd.msjava.msscorer;

import edu.ucsd.msjava.msgf.NominalMass;

// Fast scorer for DB search, consider edges
public class DBScanScorer extends FastScorer {

    private float[] nodeMass = null;
    private NewRankScorer scorer = null;
    private Partition partition;
    private float probPeak;
    private boolean isNodeMassPRM;    // prefix: true, suffix: false

    public DBScanScorer(NewScoredSpectrum<NominalMass> scoredSpec, int peptideMass) {
        super(scoredSpec, peptideMass);
        this.scorer = scoredSpec.getScorer();

        nodeMass = new float[peptideMass];

        for (int i = 0; i < nodeMass.length; i++)
            nodeMass[i] = -1;

        isNodeMassPRM = scoredSpec.getMainIonDirection();
        // assign node mass
        nodeMass[0] = 0;
        for (int nominalMass = 1; nominalMass < nodeMass.length; nominalMass++) {
            nodeMass[nominalMass] = scoredSpec.getNodeMass(new NominalMass(nominalMass));
        }

        partition = scoredSpec.getPartition();
        probPeak = scoredSpec.getProbPeak();
    }

    // fromIndex: inclusive, toIndex: exclusive
    @Override
    public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr, int fromIndex, int toIndex, int numMods) {
        int nodeScore = super.getScore(prefixMassArr, nominalPrefixMassArr, fromIndex, toIndex, numMods);
        int edgeScore = 0;
        // iter28 per-edge trace gated by `msgfplus.trace=true` (cheap: skipped when off).
        boolean trace = "true".equals(System.getProperty("msgfplus.trace"));
        if (!isNodeMassPRM)    // reverse
        {
            int nominalPeptideMass = nominalPrefixMassArr[toIndex - 1];
            for (int i = toIndex - 2; i >= fromIndex; i--) {
                int cur = nominalPeptideMass - nominalPrefixMassArr[i];
                int prev = nominalPeptideMass - nominalPrefixMassArr[i + 1];
                float theoAa = (float) (prefixMassArr[i + 1] - prefixMassArr[i]);
                int es = getEdgeScoreInt(cur, prev, theoAa);
                if (trace) {
                    System.err.println("TRACE_JAVA_EDGE\tdir=reverse\ti=" + i + "\tcur=" + cur + "\tprev=" + prev + "\ttheoAa=" + theoAa + "\tedgeScore=" + es);
                }
                edgeScore += es;
            }
        } else                    // forward
        {
            for (int i = fromIndex; i <= toIndex - 2; i++) {
                int cur = nominalPrefixMassArr[i];
                int prev = nominalPrefixMassArr[i - 1];
                float theoAa = (float) (prefixMassArr[i] - prefixMassArr[i - 1]);
                int es = getEdgeScoreInt(cur, prev, theoAa);
                if (trace) {
                    System.err.println("TRACE_JAVA_EDGE\tdir=forward\ti=" + i + "\tcur=" + cur + "\tprev=" + prev + "\ttheoAa=" + theoAa + "\tedgeScore=" + es);
                }
                edgeScore += es;
            }
        }
        if (trace) {
            System.err.println("TRACE_JAVA_EDGE_TOTAL\tnodeScore=" + nodeScore + "\tedgeScore=" + edgeScore + "\tisPrefixMain=" + isNodeMassPRM);
        }
        return nodeScore + edgeScore;
    }

    @Override
    public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
        return getEdgeScoreInt(curNode.getNominalMass(), prevNode.getNominalMass(), theoMass);
    }

    private int getEdgeScoreInt(int curNominalMass, int prevNominalMass, float theoMass) {
        if (curNominalMass >= nodeMass.length || prevNominalMass >= nodeMass.length || curNominalMass < 0 || prevNominalMass < 0)
            return 0;
        int ionExistenceIndex = 0;
        float curMass = nodeMass[curNominalMass];
        if (curMass >= 0)
            ionExistenceIndex += 1;
        float prevMass = nodeMass[prevNominalMass];
        if (prevMass >= 0)
            ionExistenceIndex += 2;

        float edgeScore = scorer.getIonExistenceScore(partition, ionExistenceIndex, probPeak);
        if (ionExistenceIndex == 3) {
            edgeScore += scorer.getErrorScore(partition, curMass - prevMass - theoMass);
        }
        return Math.round(edgeScore);
    }
}
