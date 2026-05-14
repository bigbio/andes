package edu.ucsd.msjava.msscorer;

import edu.ucsd.msjava.msgf.FlexAminoAcidGraph;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peak;

// this does not use edge scores
public class FastScorer implements SimpleDBSearchScorer<NominalMass> {

    protected float[] prefixScore = null;
    protected float[] suffixScore = null;
    private boolean mainIonDirection;
    protected Peak precursor;
    protected ActivationMethod[] activationMethodArr;
    private int[] scanNumArr;

    public FastScorer(ScoredSpectrum<NominalMass> scoredSpec, int peptideMass) {
        prefixScore = new float[peptideMass];
        suffixScore = new float[peptideMass];
        for (int i = 0; i < prefixScore.length; i++)
            prefixScore[i] = Float.MIN_VALUE;
        for (int nominalMass = 1; nominalMass < peptideMass; nominalMass++) {
            NominalMass node = new NominalMass(nominalMass);
            prefixScore[nominalMass] = scoredSpec.getNodeScore(node, true);
            suffixScore[nominalMass] = scoredSpec.getNodeScore(node, false);
        }
        mainIonDirection = scoredSpec.getMainIonDirection();

        this.precursor = scoredSpec.getPrecursorPeak();
        this.activationMethodArr = scoredSpec.getActivationMethodArr();
        this.scanNumArr = scoredSpec.getScanNumArr();
    }

    public Peak getPrecursorPeak() {
        return precursor;
    }

    public ActivationMethod[] getActivationMethodArr() {
        return activationMethodArr;
    }

    public float getParentMass() {
        return precursor.getMass();
    }

    public float getPeptideMass() {
        return precursor.getMass() - (float) (Composition.H2O);
    }

    public int getCharge() {
        return precursor.getCharge();
    }


    // fromIndex: inclusive, toIndex: exclusive
    public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr, int fromIndex, int toIndex, int numMods) {
        int score = 0;
        int peptideMass = nominalPrefixMassArr[toIndex - 1];
        for (int i = fromIndex; i < toIndex - 1; i++) {
            int prefixMass = nominalPrefixMassArr[i];
            int suffixMass = peptideMass - prefixMass;
            int curScore;
            try {
                curScore = Math.round(prefixScore[prefixMass] + suffixScore[suffixMass]);
            } catch (ArrayIndexOutOfBoundsException e) {
                curScore = 0;
            }
            score += curScore;
        }

        score += FlexAminoAcidGraph.MODIFIED_EDGE_PENALTY * numMods;
        return score;
    }

    public int getNodeScore(NominalMass prefixMass, NominalMass suffixMass) {
        int preNormMass = prefixMass.getNominalMass();
        int sufNormMass = suffixMass.getNominalMass();
        if (preNormMass >= prefixScore.length || sufNormMass >= suffixScore.length || preNormMass < 0 || sufNormMass < 0)
            return 0;
        return Math.round(prefixScore[prefixMass.getNominalMass()] + suffixScore[suffixMass.getNominalMass()]);
    }

    public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
        return 0;
    }

    public boolean getMainIonDirection() {
        return mainIonDirection;
    }

    public float getNodeScore(NominalMass node, boolean isPrefix) {
        if (isPrefix)
            return prefixScore[node.getNominalMass()];
        else
            return suffixScore[node.getNominalMass()];
    }

    public int[] getScanNumArr() {
        return scanNumArr;
    }

    // -------------------------------------------------------------------------
    // Score-traceability instrumentation
    //
    // Mirrors `getScore` exactly (same loop, same formula). Emits per-split
    // tab-separated trace lines on `System.err` so the per-split contribution
    // can be diffed against the Rust port. Triggered from
    // `DBScanner.search(...)` when `-Dmsgfplus.trace=true` and both
    // `-Dmsgfplus.trace.scan=<scan>` and `-Dmsgfplus.trace.pep=<seq>` match.
    //
    // The line format is intentionally stable across the lifetime of the
    // parity work; the Rust port emits matching `TRACE_RUST*` lines.
    // -------------------------------------------------------------------------
    public int getScoreWithTrace(double[] prefixMassArr, int[] nominalPrefixMassArr,
                                 int fromIndex, int toIndex, int numMods,
                                 String pepSeq) {
        int scan = (scanNumArr != null && scanNumArr.length > 0) ? scanNumArr[0] : -1;
        int charge = (precursor != null) ? precursor.getCharge() : -1;
        // Header: scan, peptide, charge, range, mod count, table sizes.
        System.err.println(
                "TRACE_JAVA_HEADER" +
                        "\tscan=" + scan +
                        "\tpep=" + pepSeq +
                        "\tcharge=" + charge +
                        "\tfromIndex=" + fromIndex +
                        "\ttoIndex=" + toIndex +
                        "\tnumMods=" + numMods +
                        "\tprefScoreLen=" + prefixScore.length +
                        "\tsuffScoreLen=" + suffixScore.length);

        // FastScorer does not retain the per-segment ion-type list itself;
        // print the partition-cache shape we *do* have, so the per-ion dump
        // location is at least identified. The full per-(seg, ion) breakdown
        // lives in `NewScoredSpectrum` and is collapsed into the prefix/suffix
        // score arrays here, so we cannot re-emit it from this scorer.
        System.err.println(
                "TRACE_JAVA_IONS" +
                        "\tscan=" + scan +
                        "\tpep=" + pepSeq +
                        "\tnote=FastScorer-collapsed-prefSuff-tables-no-per-segment-detail");

        int score = 0;
        int peptideMass = nominalPrefixMassArr[toIndex - 1];
        for (int i = fromIndex; i < toIndex - 1; i++) {
            int prefixMass = nominalPrefixMassArr[i];
            int suffixMass = peptideMass - prefixMass;
            int curScore;
            float prefScore;
            float suffScore;
            boolean oob = false;
            try {
                prefScore = prefixScore[prefixMass];
                suffScore = suffixScore[suffixMass];
                curScore = Math.round(prefScore + suffScore);
            } catch (ArrayIndexOutOfBoundsException e) {
                prefScore = Float.NaN;
                suffScore = Float.NaN;
                curScore = 0;
                oob = true;
            }
            score += curScore;

            System.err.println(
                    "TRACE_JAVA" +
                            "\tpep=" + pepSeq +
                            "\tsplit=" + i +
                            "\tprefMass=" + prefixMass +
                            "\tsuffMass=" + suffixMass +
                            "\tprefScore=" + prefScore +
                            "\tsuffScore=" + suffScore +
                            "\tcontribution=" + curScore +
                            "\tcumulative=" + score +
                            (oob ? "\toob=true" : ""));
        }

        int modPenalty = edu.ucsd.msjava.msgf.FlexAminoAcidGraph.MODIFIED_EDGE_PENALTY * numMods;
        score += modPenalty;
        System.err.println(
                "TRACE_JAVA_FINAL" +
                        "\tpep=" + pepSeq +
                        "\tmodPenalty=" + modPenalty +
                        "\tnumMods=" + numMods +
                        "\trawScore=" + score);
        return score;
    }

}
