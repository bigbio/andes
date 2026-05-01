package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.misc.ProgressData;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msgf.ScoredSpectrumSum;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.*;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.*;

import java.util.*;

public class ScoredSpectraMap {
    private final SpectraAccessor specAcc;
    private final List<SpecKey> specKeyList;
    private final Tolerance leftPrecursorMassTolerance;
    private final Tolerance rightPrecursorMassTolerance;
    private final int minIsotopeError;
    private final int maxIsotopeError;
    private final SpecDataType specDataType;
    /**
     * Achievement B (P2-cal) precursor mass shift in ppm. Applied to each
     * precursor mass when it first materialises from the spectrum. Zero means
     * no correction — the code path is bit-identical to a pre-calibration
     * build when this value is 0.0 (enforced by {@link #applyShift(float)}).
     */
    private final double precursorMassShiftPpm;

    private SortedMap<Double, SpecKey> pepMassSpecKeyMap;
    private Map<SpecKey, SimpleDBSearchScorer<NominalMass>> specKeyScorerMap;
    private Map<Pair<Integer, Integer>, SpecKey> specIndexChargeToSpecKeyMap;

    private Map<SpecKey, NewRankScorer> specKeyRankScorerMap;

    private boolean turnOffEdgeScoring = false;
    private boolean isolateSpectrumState = false;

    private ProgressData progress;

    public ScoredSpectraMap(
            SpectraAccessor specAcc,
            List<SpecKey> specKeyList,
            Tolerance leftPrecursorMassTolerance,
            Tolerance rightPrecursorMassTolerance,
            int minIsotopeError,
            int maxIsotopeError,
            SpecDataType specDataType,
            boolean storeRankScorer,
            boolean supportSpectrumSpecificErrorTolerance,
            double precursorMassShiftPpm
    ) {
        this.specAcc = specAcc;
        this.specKeyList = specKeyList;
        this.leftPrecursorMassTolerance = leftPrecursorMassTolerance;
        this.rightPrecursorMassTolerance = rightPrecursorMassTolerance;
        this.minIsotopeError = minIsotopeError;
        this.maxIsotopeError = maxIsotopeError;
        this.specDataType = specDataType;
        this.precursorMassShiftPpm = precursorMassShiftPpm;

        // Each ScoredSpectraMap is owned by exactly one RunMSGFPlus task (or the
        // MassCalibrator pre-pass, also single-threaded). The synchronized wrappers
        // these maps used to carry were defensive against a sharing pattern that
        // does not occur in production code paths. Plain Map/SortedMap is enough.
        pepMassSpecKeyMap = new TreeMap<>();
        specKeyScorerMap = new HashMap<>();
        specIndexChargeToSpecKeyMap = new HashMap<>();

        if (storeRankScorer)
            specKeyRankScorerMap = new HashMap<>();
        progress = null;
    }

    /**
     * Backwards-compatible ctor that defaults {@code precursorMassShiftPpm}
     * to 0.0. Existing callers that do not participate in calibration pick
     * up the no-op path and stay bit-identical.
     */
    public ScoredSpectraMap(
            SpectraAccessor specAcc,
            List<SpecKey> specKeyList,
            Tolerance leftPrecursorMassTolerance,
            Tolerance rightPrecursorMassTolerance,
            int minIsotopeError,
            int maxIsotopeError,
            SpecDataType specDataType,
            boolean storeRankScorer,
            boolean supportSpectrumSpecificErrorTolerance
    ) {
        this(specAcc, specKeyList, leftPrecursorMassTolerance, rightPrecursorMassTolerance,
                minIsotopeError, maxIsotopeError, specDataType,
                storeRankScorer, supportSpectrumSpecificErrorTolerance, 0.0);
    }

    public ScoredSpectraMap(
            SpectraAccessor specAcc,
            List<SpecKey> specKeyList,
            Tolerance leftPrecursorMassTolerance,
            Tolerance rightPrecursorMassTolerance,
            int maxNum13C,
            SpecDataType specDataType,
            boolean storeRankScorer,
            boolean supportSpectrumSpecificErrorTolerance
    ) {
        this(specAcc, specKeyList, leftPrecursorMassTolerance, rightPrecursorMassTolerance, 0, maxNum13C, specDataType, storeRankScorer, supportSpectrumSpecificErrorTolerance);
    }

    public ScoredSpectraMap(
            SpectraAccessor specAcc,
            List<SpecKey> specKeyList,
            Tolerance leftPrecursorMassTolerance,
            Tolerance rightPrecursorMassTolerance,
            int maxNum13C,
            SpecDataType specDataType,
            boolean storeRankScorer
    ) {
        this(specAcc, specKeyList, leftPrecursorMassTolerance, rightPrecursorMassTolerance, 0, maxNum13C, specDataType, storeRankScorer, false);
    }

    public ScoredSpectraMap turnOffEdgeScoring() {
        this.turnOffEdgeScoring = true;
        return this;
    }

    /**
     * Use cloned Spectrum snapshots while preprocessing so callers like the
     * calibration pre-pass do not mutate the shared SpectraAccessor cache.
     * The default remains false for the main search path to preserve current
     * behavior and allocation profile.
     */
    public ScoredSpectraMap isolateSpectrumState() {
        this.isolateSpectrumState = true;
        return this;
    }

    public SortedMap<Double, SpecKey> getPepMassSpecKeyMap() {
        return pepMassSpecKeyMap;
    }

    public Map<SpecKey, SimpleDBSearchScorer<NominalMass>> getSpecKeyScorerMap() {
        return specKeyScorerMap;
    }

    public SpectraAccessor getSpectraAccessor() {
        return specAcc;
    }

    public SpecDataType getSpecDataType() {
        return specDataType;
    }

    @Deprecated
    public Tolerance getLeftParentMassTolerance() {
        return getLeftPrecursorMassTolerance();
    }

    @Deprecated
    public Tolerance getRightParentMassTolerance() {
        return getRightPrecursorMassTolerance();
    }

    public Tolerance getLeftPrecursorMassTolerance() {
        return leftPrecursorMassTolerance;
    }

    public Tolerance getRightPrecursorMassTolerance() {
        return rightPrecursorMassTolerance;
    }

    public int getMaxIsotopeError() {
        return maxIsotopeError;
    }

    public int getMinIsotopeError() {
        return minIsotopeError;
    }

    public List<SpecKey> getSpecKeyList() {
        return specKeyList;
    }

    public SpecKey getSpecKey(int specIndex, int charge) {
        return specIndexChargeToSpecKeyMap.get(new Pair<Integer, Integer>(specIndex, charge));
    }

    public NewRankScorer getRankScorer(SpecKey specKey) {
        if (specKeyRankScorerMap == null)
            return null;
        else
            return this.specKeyRankScorerMap.get(specKey);
    }

    public ScoredSpectraMap makePepMassSpecKeyMap() {
        for (SpecKey specKey : specKeyList) {
            int specIndex = specKey.getSpecIndex();
            Spectrum spec = specAcc.getSpectrumBySpecIndex(specIndex);
            float peptideMass = (spec.getPrecursorPeak().getMz() - (float) Composition.ChargeCarrierMass()) * specKey.getCharge() - (float) Composition.H2O;
            peptideMass = applyShift(peptideMass);

            if (peptideMass > 0) {
                for (int delta = this.minIsotopeError; delta <= maxIsotopeError; delta++) {
                    float mass1 = peptideMass - delta * (float) Composition.ISOTOPE;
                    double mass1Key = (double) mass1;
                    while (pepMassSpecKeyMap.get(mass1Key) != null)
                        mass1Key = Math.nextUp(mass1Key);
                    pepMassSpecKeyMap.put(mass1Key, specKey);
                }
                specIndexChargeToSpecKeyMap.put(new Pair<Integer, Integer>(specIndex, specKey.getCharge()), specKey);

            } else {
                // Skip since precursor m/z is zero
            }
        }
        return this;
    }

    public void setProgressObj(ProgressData progObj) {
        progress = progObj;
    }

    public ProgressData getProgressObj() {
        return progress;
    }

    public void preProcessSpectra() {
        preProcessSpectra(0, specKeyList.size());
    }

    public void preProcessSpectra(int fromIndex, int toIndex) {
        if (progress == null) {
            progress = new ProgressData();
        }
        if (specDataType.getActivationMethod() != ActivationMethod.FUSION)
            preProcessIndividualSpectra(fromIndex, toIndex);
        else
            preProcessFusedSpectra(fromIndex, toIndex);
    }

    private void preProcessIndividualSpectra(int fromIndex, int toIndex) {
        NewRankScorer scorer = null;
        ActivationMethod activationMethod = specDataType.getActivationMethod();
        InstrumentType instType = specDataType.getInstrumentType();
        Enzyme enzyme = specDataType.getEnzyme();
        Protocol protocol = specDataType.getProtocol();

        if (activationMethod != ActivationMethod.ASWRITTEN && activationMethod != ActivationMethod.FUSION) {
            scorer = NewScorerFactory.get(activationMethod, instType, enzyme, protocol);
            if (this.turnOffEdgeScoring)
                scorer.doNotUseError();
        }
        int count = 0;
        int countIgnored = 0;
        int total = toIndex - fromIndex;
        for (SpecKey specKey : specKeyList.subList(fromIndex, toIndex)) {
            if (Thread.currentThread().isInterrupted()) {
                return;
            }

            int specIndex = specKey.getSpecIndex();
            Spectrum spec = specAcc.getSpectrumBySpecIndex(specIndex);
            if (activationMethod == ActivationMethod.ASWRITTEN || activationMethod == ActivationMethod.FUSION) {
                scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme, protocol);
                if (this.turnOffEdgeScoring)
                    scorer.doNotUseError();
            }
            int charge = specKey.getCharge();
            Spectrum scoringSpec = prepareSpectrumForScoring(spec, charge);

            NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(scoringSpec);

            float peptideMass = scoringSpec.getPrecursorMass() - (float) Composition.H2O;
            peptideMass = applyShift(peptideMass);
            float tolDaLeft = leftPrecursorMassTolerance.getToleranceAsDa(peptideMass);
            int maxNominalPeptideMass = NominalMass.toNominalMass(peptideMass) + Math.round(tolDaLeft - 0.4999f) - this.minIsotopeError;

            if (maxNominalPeptideMass > 0) {
                if (scorer.supportEdgeScores()) {
                    specKeyScorerMap.put(specKey, new DBScanScorer(scoredSpec, maxNominalPeptideMass));
                } else {
                    specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
                }

                if (specKeyRankScorerMap != null) {
                    specKeyRankScorerMap.put(specKey, scorer);
                }
            } else {
                countIgnored++;
                if (countIgnored <= 4) {
                    System.out.println("... ignoring spectrum at index " +
                            String.format("%1$5s", specKey.getSpecIndex()) +
                            " with invalid precursor ion of " + spec.getPrecursorMass() + " Da");
                }
            }

            count++;
            progress.report(count, total);
        }

        if (countIgnored > 1) {
            String threadName = Thread.currentThread().getName();
            System.out.println("Warning: Ignored " + countIgnored + " spectra with invalid precursor ions (" + threadName + ")");
        }
    }

    /**
     * Applies the learned precursor-mass calibration shift to a single mass.
     *
     * <p>When {@code precursorMassShiftPpm == 0.0} (the default and the
     * {@code -precursorCal off} path), this method returns the input
     * unchanged — the comparison is against the same {@code double} literal
     * that was stored in the field, so the check is exact and the code path
     * is bit-identical to a pre-calibration build. This is the non-negotiable
     * correctness gate for the feature.
     *
     * <p>When non-zero, applies {@code mass * (1 - shiftPpm * 1e-6)}, which
     * removes the positive bias learned by {@link MassCalibrator}.
     */
    private float applyShift(float peptideMass) {
        if (precursorMassShiftPpm == 0.0) {
            return peptideMass;
        }
        return peptideMass * (1.0f - (float) (precursorMassShiftPpm * 1e-6));
    }

    private void preProcessFusedSpectra(int fromIndex, int toIndex) {
        InstrumentType instType = specDataType.getInstrumentType();
        Enzyme enzyme = specDataType.getEnzyme();
        Protocol protocol = specDataType.getProtocol();

        for (SpecKey specKey : specKeyList.subList(fromIndex, toIndex)) {
            if (Thread.currentThread().isInterrupted()) {
                return;
            }

            ArrayList<Integer> specIndexList = specKey.getSpecIndexList();
            if (specIndexList == null) {
                specIndexList = new ArrayList<Integer>();
                specIndexList.add(specKey.getSpecIndex());
            }
            ArrayList<ScoredSpectrum<NominalMass>> scoredSpecList = new ArrayList<ScoredSpectrum<NominalMass>>();
            boolean supportEdgeScore = true;
            for (int specIndex : specIndexList) {
                if (Thread.currentThread().isInterrupted()) {
                    return;
                }

                Spectrum spec = specAcc.getSpectrumBySpecIndex(specIndex);

                NewRankScorer scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme, protocol);
                if (!scorer.supportEdgeScores())
                    supportEdgeScore = false;
                int charge = specKey.getCharge();
                Spectrum scoringSpec = prepareSpectrumForScoring(spec, charge);
                NewScoredSpectrum<NominalMass> sSpec = scorer.getScoredSpectrum(scoringSpec);
                scoredSpecList.add(sSpec);
            }

            if (scoredSpecList.size() == 0)
                continue;
            ScoredSpectrumSum<NominalMass> scoredSpec = new ScoredSpectrumSum<NominalMass>(scoredSpecList);
            float peptideMass = scoredSpec.getPrecursorPeak().getMass() - (float) Composition.H2O;
            float tolDaLeft = leftPrecursorMassTolerance.getToleranceAsDa(peptideMass);
            int maxNominalPeptideMass = NominalMass.toNominalMass(peptideMass) + Math.round(tolDaLeft - 0.4999f) + 1;
            if (supportEdgeScore)
                specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
            else
                specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
        }
    }

    Spectrum prepareSpectrumForScoring(Spectrum spec, int charge) {
        if (isolateSpectrumState) {
            Spectrum cloned = cloneSpectrum(spec);
            cloned.setCharge(charge);
            return cloned;
        }
        spec.setCharge(charge);
        return spec;
    }

    private static Spectrum cloneSpectrum(Spectrum spec) {
        Spectrum cloned = spec.getCloneWithoutPeakList();
        for (Peak peak : spec) {
            cloned.add(peak.clone());
        }
        return cloned;
    }
}
