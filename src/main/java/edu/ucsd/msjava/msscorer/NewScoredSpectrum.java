package edu.ucsd.msjava.msscorer;

import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.*;

public class NewScoredSpectrum<T extends Matter> implements ScoredSpectrum<T> {

    private Spectrum spec;
    private NewRankScorer scorer;
    private Tolerance mme;

    private IonType[][] ionTypes;    // segmentNum, ionType
    private final int charge;
    private final float parentMass;
    private final Peak precursor;
    private final int[] scanNumArr;
    private ActivationMethod[] activationMethodArr;
    private IonType mainIon;
    private Partition partition;    // partition of the last segment
    private float probPeak;

    public NewScoredSpectrum(Spectrum spec, NewRankScorer scorer) {
        this.scorer = scorer;

        this.charge = spec.getCharge();
        this.parentMass = spec.getPrecursorMass();
        this.mme = scorer.mme;
        this.precursor = spec.getPrecursorPeak().clone();
        this.activationMethodArr = new ActivationMethod[1];
        if (spec.getActivationMethod() != null)
            activationMethodArr[0] = spec.getActivationMethod();
        else
            activationMethodArr[0] = scorer.getSpecDataType().getActivationMethod();
        this.scanNumArr = new int[1];
        scanNumArr[0] = spec.getScanNum();

        int numSegments = scorer.getNumSegments();
        ionTypes = new IonType[numSegments][];
        for (int seg = 0; seg < numSegments; seg++)
            ionTypes[seg] = scorer.getIonTypes(charge, parentMass, seg);

        // Diagnostic: dump partition triple per seg index for the trace target.
        if ("true".equals(System.getProperty("msgfplus.trace.getnode"))
                && matchesTargetTraceScan(scanNumArr)) {
            for (int seg = 0; seg < numSegments; seg++) {
                Partition pp = scorer.getPartition(charge, parentMass, seg);
                java.util.ArrayList<FragmentOffsetFrequency> fofs = scorer.getFragmentOFF(pp);
                int pfx = 0, sfx = 0;
                if (fofs != null) {
                    for (FragmentOffsetFrequency f : fofs) {
                        if (f.getIonType() instanceof IonType.PrefixIon) pfx++;
                        else if (f.getIonType() instanceof IonType.SuffixIon) sfx++;
                    }
                }
                System.err.println("TRACE_JAVA_GN_PART"
                        + "\tsegIndex=" + seg
                        + "\tpartition_charge=" + pp.getCharge()
                        + "\tpartition_pm=" + pp.getParentMass()
                        + "\tpartition_seg=" + pp.getSegNum()
                        + "\tionCount=" + ionTypes[seg].length
                        + "\tfofCount=" + (fofs == null ? -1 : fofs.size())
                        + "\tprefixFOF=" + pfx
                        + "\tsuffixFOF=" + sfx
                        + "\tquery_charge=" + charge
                        + "\tquery_pm=" + parentMass);
            }
        }

        // filter precursor peaks
        for (PrecursorOffsetFrequency off : scorer.getPrecursorOFF(spec.getCharge()))
            spec.filterPrecursorPeaks(mme, off.getReducedCharge(), off.getOffset());
        spec.setRanksOfPeaks();

        // deconvolute spectra
        if (scorer.applyDeconvolution())
            spec = spec.getDeconvolutedSpectrum(scorer.deconvolutionErrorTolerance());

        // for edge scoring
        partition = scorer.getPartition(spec.getCharge(), spec.getPrecursorMass(), scorer.getNumSegments() - 1);
        mainIon = scorer.getMainIonType(partition);

        float approxNumBins = spec.getPeptideMass() / (scorer.getMME().getValue() * 2);

        if (spec.size() == 0)
            probPeak = 1 / Math.max(approxNumBins, 1);
        else
            probPeak = spec.size() / Math.max(approxNumBins, 1);

        this.spec = spec;
    }

    public Peak getPrecursorPeak() {
        return precursor;
    }

    public ActivationMethod[] getActivationMethodArr() {
        return this.activationMethodArr;
    }

    public int getNodeScore(T prm, T srm) {
        float prefScore = getNodeScore(prm, true);
        float sufScore = getNodeScore(srm, false);
        return Math.round(prefScore + sufScore);
    }

    public int getEdgeScore(T curNode, T prevNode, float theoMass) {
        if (!scorer.supportEdgeScores())
            return 0;

        int ionExistenceIndex = 0;
        float curNodeMass = getNodeMass(curNode);
        if (curNodeMass >= 0)
            ionExistenceIndex += 1;
        Float prevNodeMass = getNodeMass(prevNode);
        if (prevNodeMass >= 0)
            ionExistenceIndex += 2;

        float edgeScore = scorer.getIonExistenceScore(partition, ionExistenceIndex, probPeak);
        if (ionExistenceIndex == 3)
            edgeScore += scorer.getErrorScore(partition, curNodeMass - prevNodeMass - theoMass);
        return Math.round(edgeScore);
    }

    public NewRankScorer getScorer() {
        return scorer;
    }

    public Partition getPartition() {
        return partition;
    }

    public float getProbPeak() {
        return probPeak;
    }

    public IonType getMainIon() {
        return mainIon;
    }

    public boolean getMainIonDirection() {
        return mainIon.isPrefixIon();
    }

    /** Returns the corrected m/z from the observed peak, or -1 if no peak was found. */
    public float getNodeMass(T node) {
        if (node.getNominalMass() == 0)
            return 0;
        float theoMass = mainIon.getMz(node.getMass());
        Peak p = spec.getPeakByMass(theoMass, scorer.getMME());
        if (p != null)
            return mainIon.getMass(p.getMz());
        else
            return -1;
    }

    public float getNodeScore(T node, boolean isPrefix) {
        return getNodeScore(node.getMass(), isPrefix);
    }

    public float getNodeScore(float nodeMass, boolean isPrefix) {
        // Diagnostic trace gated by -Dmsgfplus.trace.getnode=true, a target-mass filter
        // (fires only for nodeMass == 974 || 1087 || 1216 || 1561) AND a scan filter
        // (only fires when scanNumArr[0] matches -Dmsgfplus.trace.scan). The scan
        // filter is essential because getNodeScore is called per-spectrum to build
        // FastScorer's prefixScore[] table - without it the trace explodes to GB.
        // When the flag is off, behavior is bit-identical to the original.
        final boolean traceEnabled =
                "true".equals(System.getProperty("msgfplus.trace.getnode"))
                && isTargetTraceMass(nodeMass)
                && matchesTargetTraceScan(scanNumArr);
        if (traceEnabled) {
            System.err.println(
                    "TRACE_JAVA_GN_HEADER"
                            + "\tnodeMass=" + nodeMass
                            + "\tisPrefix=" + isPrefix
                            + "\tnumSegments=" + scorer.getNumSegments()
                            + "\tparentMass=" + parentMass
                            + "\tcharge=" + charge);
        }

        float score = 0;
        for (int segIndex = 0; segIndex < scorer.getNumSegments(); segIndex++) {
            if (traceEnabled) {
                System.err.println(
                        "TRACE_JAVA_GN_SEG"
                                + "\tnodeMass=" + nodeMass
                                + "\tsegIndex=" + segIndex
                                + "\tionsInSeg=" + ionTypes[segIndex].length);
            }
            for (IonType ion : ionTypes[segIndex]) {
                float theoMass;
                String ionClass;
                if (isPrefix)    // prefix
                {
                    if (ion instanceof IonType.PrefixIon) {
                        theoMass = ion.getMz(nodeMass);
                        ionClass = "PrefixIon";
                    } else {
                        if (traceEnabled) {
                            System.err.println(
                                    "TRACE_JAVA_GN_ION"
                                            + "\tnodeMass=" + nodeMass
                                            + "\tsegIndex=" + segIndex
                                            + "\tion=" + ion.getName()
                                            + "\tion_class=wrong-direction"
                                            + "\ttheoMass=NA"
                                            + "\tsegNum_of_theoMass=NA"
                                            + "\tmatch_segIdx=false"
                                            + "\tpeak_found=false"
                                            + "\tpeak_rank=-1"
                                            + "\tscored=0");
                        }
                        continue;
                    }
                } else {
                    if (ion instanceof IonType.SuffixIon) {
                        theoMass = ion.getMz(nodeMass);
                        ionClass = "SuffixIon";
                    } else {
                        if (traceEnabled) {
                            System.err.println(
                                    "TRACE_JAVA_GN_ION"
                                            + "\tnodeMass=" + nodeMass
                                            + "\tsegIndex=" + segIndex
                                            + "\tion=" + ion.getName()
                                            + "\tion_class=wrong-direction"
                                            + "\ttheoMass=NA"
                                            + "\tsegNum_of_theoMass=NA"
                                            + "\tmatch_segIdx=false"
                                            + "\tpeak_found=false"
                                            + "\tpeak_rank=-1"
                                            + "\tscored=0");
                        }
                        continue;
                    }
                }

                int segNum = scorer.getSegmentNum(theoMass, parentMass);
                if (segNum != segIndex) {
                    if (traceEnabled) {
                        System.err.println(
                                "TRACE_JAVA_GN_ION"
                                        + "\tnodeMass=" + nodeMass
                                        + "\tsegIndex=" + segIndex
                                        + "\tion=" + ion.getName()
                                        + "\tion_class=" + ionClass
                                        + "\ttheoMass=" + theoMass
                                        + "\tsegNum_of_theoMass=" + segNum
                                        + "\tmatch_segIdx=false"
                                        + "\tpeak_found=false"
                                        + "\tpeak_rank=-1"
                                        + "\tscored=0");
                    }
                    continue;
                }

                Peak p = spec.getPeakByMass(theoMass, mme);
                Partition part = scorer.getPartition(charge, parentMass, segNum);

                float scored;
                if (p != null)    // peak exists
                    scored = scorer.getNodeScore(part, ion, p.getRank());
                else    // missing peak
                    scored = scorer.getMissingIonScore(part, ion);
                score += scored;

                if (traceEnabled) {
                    System.err.println(
                            "TRACE_JAVA_GN_ION"
                                    + "\tnodeMass=" + nodeMass
                                    + "\tsegIndex=" + segIndex
                                    + "\tion=" + ion.getName()
                                    + "\tion_class=" + ionClass
                                    + "\ttheoMass=" + theoMass
                                    + "\tsegNum_of_theoMass=" + segNum
                                    + "\tmatch_segIdx=true"
                                    + "\tpeak_found=" + (p != null)
                                    + "\tpeak_rank=" + (p != null ? p.getRank() : -1)
                                    + "\tscored=" + scored);
                }
            }
        }

        if (traceEnabled) {
            System.err.println(
                    "TRACE_JAVA_GN_TOTAL"
                            + "\tnodeMass=" + nodeMass
                            + "\tisPrefix=" + isPrefix
                            + "\ttotal=" + score);
        }
        return score;
    }

    /**
     * Diagnostic helper: filter for {@link #getNodeScore(float, boolean)} tracing.
     * Returns true only for the four prefix-mass cliff points being investigated
     * (974, 1087, 1216, 1561). Comparison uses Math.round to tolerate the float
     * inputs (nominal masses are integers but they enter as float).
     */
    private static boolean isTargetTraceMass(float nodeMass) {
        // Use truncation (not round) because FastScorer's prefixScore[N] precompute
        // calls getNodeScore(new NominalMass(N).getMass(), true), where
        // NominalMass.getMass() = N / 0.999497 (INTEGER_MASS_SCALER). For N=1087
        // this yields 1087.547f which Math.round() pushes UP to 1088 -- missing the
        // target. (int) truncation maps 1087.547 -> 1087 as intended.
        int m = (int) nodeMass;
        return m == 974 || m == 1087 || m == 1216 || m == 1561;
    }

    /**
     * Diagnostic helper: returns true if -Dmsgfplus.trace.scan is unset (trace all)
     * OR if the spectrum's scan matches the configured target scan.
     */
    private static boolean matchesTargetTraceScan(int[] scanNumArr) {
        String s = System.getProperty("msgfplus.trace.scan");
        if (s == null || s.isEmpty()) return true;
        int target;
        try { target = Integer.parseInt(s); } catch (NumberFormatException e) { return true; }
        if (scanNumArr == null || scanNumArr.length == 0) return false;
        return scanNumArr[0] == target;
    }

    public float getExplainedIonCurrent(float residueMass, boolean isPrefix, Tolerance fragmentTolerance) {
        float explainedIonCurrent = 0;
        for (int segIndex = 0; segIndex < scorer.getNumSegments(); segIndex++) {
            for (IonType ion : ionTypes[segIndex]) {
                float theoMass;
                if (isPrefix)    // prefix
                {
                    if (ion instanceof IonType.PrefixIon)
                        theoMass = ion.getMz(residueMass);
                    else
                        continue;
                } else {
                    if (ion instanceof IonType.SuffixIon)
                        theoMass = ion.getMz(residueMass);
                    else
                        continue;
                }

                int segNum = scorer.getSegmentNum(theoMass, parentMass);
                if (segNum != segIndex)
                    continue;

                Peak p = spec.getPeakByMass(theoMass, fragmentTolerance);

                if (p != null)    // peak exists
                    explainedIonCurrent += p.getIntensity();
            }
        }
        return explainedIonCurrent;
    }

    public Pair<Float, Float> getMassErrorWithIntensity(float residueMass, boolean isPrefix, Tolerance fragmentTolerance) {
        Float error = null;
        float maxIntensity = 0;

        for (int segIndex = 0; segIndex < scorer.getNumSegments(); segIndex++) {
            for (IonType ion : ionTypes[segIndex]) {
                if (ion.getCharge() != 1)
                    continue;
                float theoMass;
                if (isPrefix)    // prefix
                {
                    if (ion instanceof IonType.PrefixIon)
                        theoMass = ion.getMz(residueMass);
                    else
                        continue;
                } else {
                    if (ion instanceof IonType.SuffixIon)
                        theoMass = ion.getMz(residueMass);
                    else
                        continue;
                }

                int segNum = scorer.getSegmentNum(theoMass, parentMass);
                if (segNum != segIndex)
                    continue;

                Peak p = spec.getPeakByMass(theoMass, fragmentTolerance);

                if (p != null)    // peak exists
                {
                    float err = (p.getMz() - theoMass) / theoMass * 1e6f;
                    float intensity = p.getIntensity();
                    if (intensity > maxIntensity) {
                        error = err;
                        maxIntensity = intensity;
                    }
                }
            }
        }
        if (error == null)
            return null;
        else {
            return new Pair<Float, Float>(error, maxIntensity);
        }
    }

    public Pair<Float, Float> getNodeMassAndScore(float residueMass, boolean isPrefix) {
        Float nodeMass = null;
        float nodeScore = 0;
        float curBestScore = 0;

        for (int segIndex = 0; segIndex < scorer.getNumSegments(); segIndex++) {
            for (IonType ion : ionTypes[segIndex]) {
                float theoMass;
                if (isPrefix)    // prefix
                {
                    if (ion instanceof IonType.PrefixIon)
                        theoMass = ion.getMz(residueMass);
                    else
                        continue;
                } else {
                    if (ion instanceof IonType.SuffixIon)
                        theoMass = ion.getMz(residueMass);
                    else
                        continue;
                }

                int segNum = scorer.getSegmentNum(theoMass, parentMass);
                if (segNum != segIndex)
                    continue;

                Peak p = spec.getPeakByMass(theoMass, mme);
                Partition part = scorer.getPartition(charge, parentMass, segNum);

                if (p != null)    // peak exists
                {
                    float score = scorer.getNodeScore(part, ion, p.getRank());
                    if (ion.getCharge() == 1 && score > curBestScore) {
                        nodeMass = ion.getMass(p.getMz());
                        curBestScore = score;
                    }
                    nodeScore += score;
                } else    // missing peak
                {
                    nodeScore += scorer.getMissingIonScore(part, ion);
                }
            }
        }
        return new Pair<Float, Float>(nodeMass, nodeScore);
    }

    public int[] getScanNumArr() {
        return scanNumArr;
    }
}
