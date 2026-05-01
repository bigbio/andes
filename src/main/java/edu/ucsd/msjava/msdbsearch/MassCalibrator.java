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
    /** Conservative lower bound for a tightened ppm half-window. */
    public static final float DEFAULT_TIGHTENED_WINDOW_FLOOR_PPM = 2.0f;
    /** Safety margin added after converting MAD to a Gaussian-equivalent sigma. */
    public static final float DEFAULT_TIGHTENED_WINDOW_MARGIN_PPM = 0.5f;
    /** Number of robust sigmas to keep when tightening precursor windows. */
    public static final float DEFAULT_TIGHTENED_WINDOW_SIGMA_MULTIPLIER = 3.0f;
    /** Gaussian-equivalent scale factor for MAD. */
    private static final double MAD_TO_SIGMA_SCALE = 1.4826;
    /**
     * Reject residuals whose magnitude exceeds this threshold. A genuine mass-accuracy
     * residual on any modern instrument is well under 50 ppm; values above this almost
     * always come from isotope-error matches (e.g. M+1 isotope at +1.003 Da on a 2 kDa
     * peptide = ~500 ppm residual) admitted by a wide {@code -ti} window. Filtering
     * before computing median + MAD prevents these outliers from contaminating the
     * robust spread estimate. Empirically the residual distribution drops off well
     * before this floor; isotope-shift contamination clusters near integer multiples
     * of (1.003 / mass) ppm.
     */
    static final double MAX_REASONABLE_RESIDUAL_PPM = 50.0;
    /** Sample every Nth SpecKey. Cap total sampled keys at {@link #maxSampled}. */
    private static final int SAMPLING_STRIDE = 10;
    /** Default upper bound on sampled spectra in the pre-pass. */
    public static final int DEFAULT_MAX_SAMPLED = 500;
    /** Default minimum PSMs required before the learned shift is considered reliable. */
    public static final int DEFAULT_MIN_CONFIDENT_PSMS = 200;
    /** System property to override {@link #DEFAULT_MAX_SAMPLED} at runtime. */
    public static final String MAX_SAMPLED_PROPERTY = "msgfplus.maxSampled";
    /** System property to override {@link #DEFAULT_MIN_CONFIDENT_PSMS} at runtime. */
    public static final String MIN_CONFIDENT_PSMS_PROPERTY = "msgfplus.minConfidentPsms";
    /** SpecEValue threshold for "confident" pre-pass PSMs. Tight enough to exclude decoys. */
    private static final double MAX_SPEC_EVALUE = 1e-6;
    /**
     * Size-guard threshold in SpecKeys. Below this, skip the pre-pass entirely.
     * SpecKey count is typically ~3× the spectrum count because charges 2-4 each get
     * their own SpecKey. The 10_000 threshold means "skip on anything smaller than a
     * ~3000-spectrum file" — too small to yield 200 confident PSMs reliably, and
     * small enough that the pre-pass's Spectrum-state mutation side-effect (which
     * would otherwise drift off-mode vs auto-mode results) is visible at unit-test
     * scale. Real datasets (PXD001819 ~66K SpecKeys, Astral ~75K, TMT ~40K) are
     * comfortably above this and run the calibrator as intended.
     */
    private static final int MIN_SPECKEYS_FOR_PREPASS = 10_000;

    private final SpectraAccessor specAcc;
    private final CompactSuffixArray sa;
    private final AminoAcidSet aaSet;
    private final SearchParams params;
    private final List<SpecKey> specKeyList;
    private final Tolerance leftPrecursorMassTolerance;
    private final Tolerance rightPrecursorMassTolerance;
    private final SpecDataType specDataType;
    /** Effective sampling cap; {@link #DEFAULT_MAX_SAMPLED} unless overridden via {@link #MAX_SAMPLED_PROPERTY}. */
    private final int maxSampled;
    /** Effective stratification floor; {@link #DEFAULT_MIN_CONFIDENT_PSMS} unless overridden via {@link #MIN_CONFIDENT_PSMS_PROPERTY}. */
    private final int minConfidentPsms;

    /** Immutable summary of the sampled calibration residuals for one file. */
    public static final class CalibrationStats {
        private final double shiftPpm;
        private final double robustSigmaPpm;
        private final int confidentPsmCount;

        public CalibrationStats(double shiftPpm, double robustSigmaPpm, int confidentPsmCount) {
            this.shiftPpm = shiftPpm;
            this.robustSigmaPpm = robustSigmaPpm;
            this.confidentPsmCount = confidentPsmCount;
        }

        public double getShiftPpm() {
            return shiftPpm;
        }

        public double getRobustSigmaPpm() {
            return robustSigmaPpm;
        }

        public int getConfidentPsmCount() {
            return confidentPsmCount;
        }

        public boolean hasReliableStats() {
            // The calibrator emits confidentPsmCount > 0 only when residuals
            // cleared the (configurable) minConfidentPsms threshold.
            return confidentPsmCount > 0;
        }
    }

    /**
     * @param specAcc spectra accessor for the current file (already MS-level filtered)
     * @param sa compact suffix array for the target/decoy database
     * @param aaSet amino acid set with modifications applied
     * @param params parsed search params (used for enzyme, de novo score threshold, etc.)
     * @param specKeyList the full list of SpecKeys for the file; the calibrator
     *                    samples every {@value #SAMPLING_STRIDE}th entry up to
     *                    {@value #DEFAULT_MAX_SAMPLED} (override via
     *                    system property {@code msgfplus.maxSampled}).
     * @param leftPrecursorMassTolerance main-pass left tolerance (reused for the pre-pass)
     * @param rightPrecursorMassTolerance main-pass right tolerance (reused for the pre-pass)
     * @param specDataType scoring metadata (activation, instrument, enzyme, protocol)
     *
     * Note: the user's {@code -ti} isotope-error window is intentionally NOT
     * propagated to the pre-pass. The pre-pass is fixed to isotope error 0 to
     * prevent isotope-shift contamination of the residual distribution.
     * See {@link #collectResiduals(int)}.
     */
    public MassCalibrator(
            SpectraAccessor specAcc,
            CompactSuffixArray sa,
            AminoAcidSet aaSet,
            SearchParams params,
            List<SpecKey> specKeyList,
            Tolerance leftPrecursorMassTolerance,
            Tolerance rightPrecursorMassTolerance,
            SpecDataType specDataType
    ) {
        this.specAcc = specAcc;
        this.sa = sa;
        this.aaSet = aaSet;
        this.params = params;
        this.specKeyList = specKeyList;
        this.leftPrecursorMassTolerance = leftPrecursorMassTolerance;
        this.rightPrecursorMassTolerance = rightPrecursorMassTolerance;
        this.specDataType = specDataType;
        this.maxSampled = readPositiveIntProperty(MAX_SAMPLED_PROPERTY, DEFAULT_MAX_SAMPLED);
        this.minConfidentPsms = readPositiveIntProperty(MIN_CONFIDENT_PSMS_PROPERTY, DEFAULT_MIN_CONFIDENT_PSMS);
    }

    /** Public accessor used by unit tests to exercise property parsing. */
    public static int readPositiveIntPropertyForTests(String name, int defaultValue) {
        return readPositiveIntProperty(name, defaultValue);
    }

    /**
     * Reads a positive-integer system property; falls back to {@code defaultValue}
     * for unset / non-numeric / non-positive values.
     */
    private static int readPositiveIntProperty(String name, int defaultValue) {
        String raw = System.getProperty(name);
        if (raw == null || raw.isEmpty()) return defaultValue;
        try {
            int parsed = Integer.parseInt(raw.trim());
            return parsed > 0 ? parsed : defaultValue;
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }

    /**
     * Runs the sampled pre-pass and returns the median ppm shift, or
     * {@code 0.0} if fewer than {@value #DEFAULT_MIN_CONFIDENT_PSMS} (override
     * via {@code msgfplus.minConfidentPsms}) high-confidence
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
        return learnCalibrationStats(ioIndex).getShiftPpm();
    }

    /**
     * Runs the sampled pre-pass and returns both the learned median shift and a
     * robust spread estimate for later tolerance tightening.
     */
    public CalibrationStats learnCalibrationStats(int ioIndex) {
        // Skip the pre-pass on small files where minConfidentPsms can't be reached.
        if (specKeyList == null || specKeyList.size() < MIN_SPECKEYS_FOR_PREPASS) {
            return new CalibrationStats(0.0, 0.0, 0);
        }
        List<Double> residuals = collectResiduals(ioIndex);
        if (residuals.size() < minConfidentPsms) {
            // count=0 is the "unreliable, do not apply" sentinel; CalibrationStats.hasReliableStats()
            // checks for count > 0.
            return new CalibrationStats(0.0, 0.0, 0);
        }
        double shiftPpm = median(residuals);
        double robustSigmaPpm = robustSigmaPpm(residuals, shiftPpm);
        return new CalibrationStats(shiftPpm, robustSigmaPpm, residuals.size());
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

        List<SpecKey> sampled = sampleEveryNth(specKeyList, SAMPLING_STRIDE, maxSampled);
        if (sampled.isEmpty()) {
            return Collections.emptyList();
        }

        // Force isotope error to 0 for the pre-pass: residuals are only meaningful
        // when the matched peptide's monoisotopic mass equals the observed precursor's
        // monoisotopic mass. With the user's wider -ti window (e.g. -1,2 on Astral),
        // PSMs whose precursor is the M+1 or M+2 isotope inject ~500 / ~1000 ppm
        // residuals into the pre-pass, contaminating median + MAD. Restricting the
        // pre-pass to isotope error 0 keeps the residual distribution clean.
        // numPeptidesPerSpec = 1 keeps the pre-pass tiny and fast. precursorMassShiftPpm = 0.0
        // because the whole point of the pre-pass is to LEARN the shift.
        ScoredSpectraMap prePassMap = new ScoredSpectraMap(
                specAcc,
                sampled,
                leftPrecursorMassTolerance,
                rightPrecursorMassTolerance,
                0,  // pre-pass minIsotopeError (overrides user's -ti to keep residuals clean)
                0,  // pre-pass maxIsotopeError
                specDataType,
                false, // storeRankScorer not needed for pre-pass
                false
        ).isolateSpectrumState();
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

        // Collect (residual, eValue) pairs so we can keep the cleanest subset
        // by spec_eValue. Stratification on a 393-PSM Astral pre-pass showed
        // sigma drops 4x (3.99 -> 0.99 ppm) when restricted to the top-200
        // most confident PSMs. Worst-half PSMs add residual scatter without
        // adding signal — they get filtered out post-collection.
        List<double[]> residualWithEval = new ArrayList<>();

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
            double residual = residualPpm(observedPeptideMass, theoreticalPeptideMass);
            // Reject isotope-error contamination before robust-stats aggregation.
            // See MAX_REASONABLE_RESIDUAL_PPM doc.
            if (Math.abs(residual) > MAX_REASONABLE_RESIDUAL_PPM) {
                continue;
            }
            residualWithEval.add(new double[]{residual, top.getSpecEValue()});
        }

        // Keep the top minConfidentPsms by spec_eValue (lowest eValue =
        // most confident). On Astral this drops sigma from ~4 ppm to ~1 ppm
        // because the worst-half PSMs (eValue near the 1e-6 threshold) are
        // dominated by residual scatter, not real instrument bias.
        residualWithEval.sort((a, b) -> Double.compare(a[1], b[1]));
        int keepN = Math.min(residualWithEval.size(), minConfidentPsms);
        for (int i = 0; i < keepN; i++) {
            residuals.add(residualWithEval.get(i)[0]);
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

    /**
     * Median absolute deviation around a known median. Empty list => 0.0.
     */
    static double medianAbsoluteDeviation(List<Double> values, double center) {
        if (values == null || values.isEmpty()) {
            return 0.0;
        }
        List<Double> deviations = new ArrayList<>(values.size());
        for (double value : values) {
            deviations.add(Math.abs(value - center));
        }
        return median(deviations);
    }

    /**
     * Robust Gaussian-equivalent sigma estimate derived from MAD.
     */
    static double robustSigmaPpm(List<Double> residuals, double center) {
        return MAD_TO_SIGMA_SCALE * medianAbsoluteDeviation(residuals, center);
    }

    /**
     * Conservative tightened ppm half-window for a calibrated main pass.
     */
    public static float tightenedTolerancePpm(float userPpm, double robustSigmaPpm, float sigmaMultiplier,
                                              float floorPpm, float marginPpm) {
        if (userPpm <= 0) {
            return userPpm;
        }
        double tightened = Math.max(floorPpm, sigmaMultiplier * robustSigmaPpm + marginPpm);
        return (float) Math.min(userPpm, tightened);
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

    /** Test-only access to {@link #medianAbsoluteDeviation(List, double)}. */
    public static double medianAbsoluteDeviationForTests(List<Double> values, double center) {
        return medianAbsoluteDeviation(values, center);
    }

    /** Test-only access to {@link #robustSigmaPpm(List, double)}. */
    public static double robustSigmaPpmForTests(List<Double> residuals, double center) {
        return robustSigmaPpm(residuals, center);
    }

    /** Test-only access to {@link #tightenedTolerancePpm(float, double, float, float, float)}. */
    public static float tightenedTolerancePpmForTests(float userPpm, double robustSigmaPpm,
                                                      float sigmaMultiplier, float floorPpm,
                                                      float marginPpm) {
        return tightenedTolerancePpm(userPpm, robustSigmaPpm, sigmaMultiplier, floorPpm, marginPpm);
    }
}
