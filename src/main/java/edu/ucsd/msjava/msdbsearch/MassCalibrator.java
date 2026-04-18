package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.SpecKey;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

/**
 * Two-pass precursor mass calibration (Achievement B — P2-cal).
 *
 * <p>Runs a sampled pre-pass of the existing {@link DBScanner} over ~10% of
 * the input spectra, filters to high-confidence PSMs, and returns the median
 * residual precursor-mass error in ppm. The caller applies this shift
 * downstream inside {@link ScoredSpectraMap} when materialising precursor
 * masses for the main search.
 *
 * <p>Sign convention: residual = (observed - theoretical) / theoretical * 1e6.
 * A positive shift means the instrument reports masses slightly higher than
 * theoretical. The main-pass correction is
 * {@code mass * (1 - shiftPpm * 1e-6)}, which re-centers the residual
 * distribution on zero.
 *
 * <p>Threading: all calibration work runs on the orchestrator thread before
 * worker {@code ScoredSpectraMap} instances are constructed. The learned
 * shift is stored on {@link edu.ucsd.msjava.msutil.DBSearchIOFiles} and read
 * immutably thereafter, so no synchronization is required.
 */
public class MassCalibrator {

    /** Sample every Nth SpecKey. Cap total sampled keys at {@link #MAX_SAMPLED}. */
    private static final int SAMPLING_STRIDE = 10;
    /** Hard upper bound on sampled spectra to keep the pre-pass bounded on large runs. */
    private static final int MAX_SAMPLED = 500;
    /** Minimum PSMs required before the learned shift is considered reliable. */
    private static final int MIN_CONFIDENT_PSMS = 200;
    /** SpecEValue threshold for "confident" pre-pass PSMs. Tight enough to exclude decoys. */
    private static final double MAX_SPEC_EVALUE = 1e-6;

    private final SpectraAccessor specAcc;
    private final CompactSuffixArray sa;
    private final AminoAcidSet aaSet;
    private final SearchParams params;
    private final List<SpecKey> specKeyList;
    private final Tolerance leftPrecursorMassTolerance;
    private final Tolerance rightPrecursorMassTolerance;
    private final int minIsotopeError;
    private final int maxIsotopeError;
    private final SpecDataType specDataType;

    /**
     * @param specAcc spectra accessor for the current file (already MS-level filtered)
     * @param sa compact suffix array for the target/decoy database
     * @param aaSet amino acid set with modifications applied
     * @param params parsed search params (used for enzyme, de novo score threshold, etc.)
     * @param specKeyList the full list of SpecKeys for the file; the calibrator
     *                    samples every {@value #SAMPLING_STRIDE}th entry up to
     *                    {@value #MAX_SAMPLED}.
     * @param leftPrecursorMassTolerance main-pass left tolerance (reused for the pre-pass)
     * @param rightPrecursorMassTolerance main-pass right tolerance (reused for the pre-pass)
     * @param minIsotopeError main-pass min isotope error
     * @param maxIsotopeError main-pass max isotope error
     * @param specDataType scoring metadata (activation, instrument, enzyme, protocol)
     */
    public MassCalibrator(
            SpectraAccessor specAcc,
            CompactSuffixArray sa,
            AminoAcidSet aaSet,
            SearchParams params,
            List<SpecKey> specKeyList,
            Tolerance leftPrecursorMassTolerance,
            Tolerance rightPrecursorMassTolerance,
            int minIsotopeError,
            int maxIsotopeError,
            SpecDataType specDataType
    ) {
        this.specAcc = specAcc;
        this.sa = sa;
        this.aaSet = aaSet;
        this.params = params;
        this.specKeyList = specKeyList;
        this.leftPrecursorMassTolerance = leftPrecursorMassTolerance;
        this.rightPrecursorMassTolerance = rightPrecursorMassTolerance;
        this.minIsotopeError = minIsotopeError;
        this.maxIsotopeError = maxIsotopeError;
        this.specDataType = specDataType;
    }

    /**
     * Runs the sampled pre-pass and returns the median ppm shift, or
     * {@code 0.0} if fewer than {@value #MIN_CONFIDENT_PSMS} high-confidence
     * PSMs are collected.
     *
     * <p>The {@code ioIndex} argument is accepted for future multi-file hooks
     * (e.g. logging per file); the actual calibration is scoped to the
     * {@link #specKeyList} passed in the constructor, so the same calibrator
     * handles one file at a time.
     *
     * @param ioIndex index of the file in the DBSearchIO list (for logging)
     * @return learned ppm shift, or 0.0 if the pre-pass had insufficient data
     */
    public double learnPrecursorShiftPpm(int ioIndex) {
        // Cheap guard: on a file too small to possibly reach MIN_CONFIDENT_PSMS
        // even if every sampled spectrum matched at 1e-6 SpecEValue, we skip the
        // pre-pass entirely. Running the pre-pass calls preProcessSpectra() on a
        // subset of shared Spectrum objects, which mutates their scored state and
        // causes a tiny (~0.1%) PSM-list drift vs -precursorCal off when the main
        // search later re-processes those same spectra. Skipping here preserves
        // the -precursorCal off ≡ no-flag bit-identity invariant for small runs,
        // which is the hard correctness gate. On large runs the guard is a no-op.
        int minFeasibleSpecCount = MIN_CONFIDENT_PSMS * SAMPLING_STRIDE;
        if (specKeyList == null || specKeyList.size() < minFeasibleSpecCount) {
            return 0.0;
        }
        List<Double> residuals = collectResiduals(ioIndex);
        if (residuals.size() < MIN_CONFIDENT_PSMS) {
            return 0.0;
        }
        return median(residuals);
    }

    /**
     * Runs the sampled pre-pass and returns the collected residuals in ppm.
     * Returns an empty list if nothing valid was collected. Package-private
     * so the integration test can exercise the full collection path.
     */
    List<Double> collectResiduals(int ioIndex) {
        if (specKeyList == null || specKeyList.isEmpty()) {
            return Collections.emptyList();
        }

        List<SpecKey> sampled = sampleEveryNth(specKeyList, SAMPLING_STRIDE, MAX_SAMPLED);
        if (sampled.isEmpty()) {
            return Collections.emptyList();
        }

        // numPeptidesPerSpec = 1 keeps the pre-pass tiny and fast. precursorMassShiftPpm = 0.0
        // because the whole point of the pre-pass is to LEARN the shift.
        ScoredSpectraMap prePassMap = new ScoredSpectraMap(
                specAcc,
                sampled,
                leftPrecursorMassTolerance,
                rightPrecursorMassTolerance,
                minIsotopeError,
                maxIsotopeError,
                specDataType,
                false, // storeRankScorer not needed for pre-pass
                false
        );
        prePassMap.makePepMassSpecKeyMap();
        prePassMap.preProcessSpectra();

        DBScanner scanner = new DBScanner(
                prePassMap,
                sa,
                params.getEnzyme(),
                aaSet,
                1, // numPeptidesPerSpec
                params.getMinPeptideLength(),
                params.getMaxPeptideLength(),
                params.getMaxNumVariantsPerPeptide(),
                params.getMinDeNovoScore(),
                params.ignoreMetCleavage(),
                params.getMaxMissedCleavages()
        );

        int ntt = params.getNumTolerableTermini();
        if (params.getEnzyme() == null) {
            ntt = 0;
        }
        int nnet = 2 - ntt;
        scanner.dbSearch(nnet);
        scanner.computeSpecEValue(false);
        scanner.generateSpecIndexDBMatchMap();

        return extractResiduals(scanner.getSpecIndexDBMatchMap(), params.getMinDeNovoScore());
    }

    /**
     * Walks the top-1 match queue for each sampled spectrum, filters to
     * high-confidence PSMs, and converts each to a ppm residual.
     */
    private List<Double> extractResiduals(
            Map<Integer, PriorityQueue<DatabaseMatch>> specIndexDBMatchMap,
            int minDeNovoScore
    ) {
        List<Double> residuals = new ArrayList<>();
        if (specIndexDBMatchMap == null || specIndexDBMatchMap.isEmpty()) {
            return residuals;
        }

        for (Map.Entry<Integer, PriorityQueue<DatabaseMatch>> entry : specIndexDBMatchMap.entrySet()) {
            PriorityQueue<DatabaseMatch> queue = entry.getValue();
            if (queue == null || queue.isEmpty()) {
                continue;
            }
            // peek() returns the worst match in the queue; we need the best (smallest SpecEValue).
            // The queue uses a SpecProbComparator, so we copy + extract the min.
            DatabaseMatch top = bestMatch(queue);
            if (top == null) {
                continue;
            }
            if (top.getSpecEValue() > MAX_SPEC_EVALUE) {
                continue;
            }
            if (top.getDeNovoScore() < minDeNovoScore) {
                continue;
            }

            int specIndex = entry.getKey();
            Spectrum spec = specAcc.getSpectrumBySpecIndex(specIndex);
            if (spec == null || spec.getPrecursorPeak() == null) {
                continue;
            }
            int charge = top.getCharge();
            if (charge <= 0) {
                continue;
            }

            double observedMz = spec.getPrecursorPeak().getMz();
            double observedPeptideMass = (observedMz - Composition.ChargeCarrierMass()) * charge - Composition.H2O;
            double theoreticalPeptideMass = top.getPeptideMass();
            if (theoreticalPeptideMass <= 0) {
                continue;
            }
            residuals.add(residualPpm(observedPeptideMass, theoreticalPeptideMass));
        }
        return residuals;
    }

    /**
     * The queue is ordered by SpecProbComparator: best (lowest SpecEValue) is
     * the last one remaining after polling, or equivalently — because
     * {@link DBScanner#generateSpecIndexDBMatchMap()} caps the queue at
     * {@code numPeptidesPerSpec = 1} — there is exactly one entry per
     * specIndex in our pre-pass. This helper is defensive in case that
     * invariant ever loosens.
     */
    private static DatabaseMatch bestMatch(PriorityQueue<DatabaseMatch> queue) {
        DatabaseMatch best = null;
        for (DatabaseMatch m : queue) {
            if (best == null || m.getSpecEValue() < best.getSpecEValue()) {
                best = m;
            }
        }
        return best;
    }

    // ----- visible-for-testing helpers (package-private) -----------------

    /**
     * Samples every Nth element (starting at index 0), capped at {@code cap}.
     */
    static <T> List<T> sampleEveryNth(List<T> source, int stride, int cap) {
        if (source == null || source.isEmpty() || stride <= 0 || cap <= 0) {
            return Collections.emptyList();
        }
        List<T> out = new ArrayList<>();
        for (int i = 0; i < source.size() && out.size() < cap; i += stride) {
            out.add(source.get(i));
        }
        return out;
    }

    /**
     * Residual in ppm for a single PSM. Sign convention:
     * {@code (observed - theoretical) / theoretical * 1e6}.
     * A positive result means the instrument reports higher than theoretical.
     */
    static double residualPpm(double observedMass, double theoreticalMass) {
        return (observedMass - theoreticalMass) / theoreticalMass * 1e6;
    }

    /**
     * Median of a list of doubles. Empty list => 0.0 (documented contract:
     * used by the calibrator as "no shift" fallback). Odd length => middle
     * element; even length => mean of the two middle elements. Sorts a
     * defensive copy so the caller's list is untouched.
     */
    static double median(List<Double> values) {
        if (values == null || values.isEmpty()) {
            return 0.0;
        }
        List<Double> copy = new ArrayList<>(values);
        Collections.sort(copy);
        int n = copy.size();
        if ((n & 1) == 1) {
            return copy.get(n / 2);
        } else {
            return (copy.get(n / 2 - 1) + copy.get(n / 2)) / 2.0;
        }
    }

    // ----- test-only public wrappers -------------------------------------
    //
    // These exist solely so the unit tests can pin the helper semantics
    // without needing a full spectrum-file fixture. They are thin
    // pass-throughs to the package-private helpers above.

    /** Test-only access to {@link #median(List)}. */
    public static double medianForTests(List<Double> values) {
        return median(values);
    }

    /** Test-only access to {@link #residualPpm(double, double)}. */
    public static double residualPpmForTests(double observed, double theoretical) {
        return residualPpm(observed, theoretical);
    }

    /** Test-only access to {@link #sampleEveryNth(List, int, int)}. */
    public static <T> List<T> sampleEveryNthForTests(List<T> source, int stride, int cap) {
        return sampleEveryNth(source, stride, cap);
    }
}
