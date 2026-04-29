package msgfplus;

import edu.ucsd.msjava.msdbsearch.MassCalibrator;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for {@link MassCalibrator} helpers.
 *
 * Pins the median + residual-ppm conventions that the rest of Achievement B
 * (two-pass precursor mass calibration) relies on. If these contracts move,
 * the whole calibration changes sign or starts drifting, so they are worth
 * nailing down explicitly.
 */
public class TestMassCalibrator {

    // ---- median() helper -------------------------------------------------

    @Test
    public void medianOdd() {
        Assert.assertEquals(3.0,
                MassCalibrator.medianForTests(new ArrayList<>(Arrays.asList(1.0, 3.0, 5.0))),
                1e-12);
    }

    @Test
    public void medianEven() {
        Assert.assertEquals(2.5,
                MassCalibrator.medianForTests(new ArrayList<>(Arrays.asList(1.0, 2.0, 3.0, 4.0))),
                1e-12);
    }

    @Test
    public void medianEmptyReturnsZero() {
        // Contract: an empty list returns 0.0 (no shift) rather than throwing,
        // so that the caller's "insufficient data" branch is trivially safe.
        Assert.assertEquals(0.0,
                MassCalibrator.medianForTests(Collections.emptyList()),
                0.0);
    }

    @Test
    public void medianUnsortedInput() {
        // Input is not required to be pre-sorted; helper sorts a defensive copy.
        Assert.assertEquals(3.0,
                MassCalibrator.medianForTests(new ArrayList<>(Arrays.asList(5.0, 1.0, 3.0))),
                1e-12);
    }

    @Test
    public void medianRobustToOutliers() {
        // This is why the calibrator uses the median, not the mean: a single
        // rogue match (e.g. a mis-assigned isotope peak) should not drag the
        // learned shift.
        Assert.assertEquals(3.0,
                MassCalibrator.medianForTests(new ArrayList<>(Arrays.asList(1.0, 2.0, 3.0, 4.0, 1000.0))),
                1e-12);
    }

    @Test
    public void medianSingleElement() {
        Assert.assertEquals(7.5,
                MassCalibrator.medianForTests(new ArrayList<>(Arrays.asList(7.5))),
                1e-12);
    }

    @Test
    public void medianAbsoluteDeviationUsesProvidedCenter() {
        List<Double> values = new ArrayList<>(Arrays.asList(1.0, 2.0, 4.0, 7.0));
        // Deviations from center=3 are [2,1,1,4] -> sorted [1,1,2,4] -> median 1.5
        Assert.assertEquals(1.5,
                MassCalibrator.medianAbsoluteDeviationForTests(values, 3.0),
                1e-12);
    }

    @Test
    public void robustSigmaPpmScalesMad() {
        List<Double> residuals = new ArrayList<>(Arrays.asList(9.0, 10.0, 11.0));
        // center=10, MAD=1 -> robust sigma = 1.4826
        Assert.assertEquals(1.4826,
                MassCalibrator.robustSigmaPpmForTests(residuals, 10.0),
                1e-6);
    }

    @Test
    public void tightenedTolerancePpmRespectsUserUpperBound() {
        float tightened = MassCalibrator.tightenedTolerancePpmForTests(
                10.0f, 0.2, 3.0f, 2.0f, 0.5f);
        // k*sigma + margin = 1.1, floor dominates -> 2.0 ppm
        Assert.assertEquals(2.0f, tightened, 1e-6f);
    }

    @Test
    public void tightenedTolerancePpmDoesNotExpandAlreadyTightWindow() {
        float tightened = MassCalibrator.tightenedTolerancePpmForTests(
                1.5f, 0.2, 3.0f, 2.0f, 0.5f);
        Assert.assertEquals(1.5f, tightened, 1e-6f);
    }

    @Test
    public void tightenedTolerancePpmTracksRobustSigmaWhenLargerThanFloor() {
        float tightened = MassCalibrator.tightenedTolerancePpmForTests(
                12.0f, 1.0, 3.0f, 2.0f, 0.5f);
        Assert.assertEquals(3.5f, tightened, 1e-6f);
    }

    @Test
    public void calibrationStatsCanBeReliableWithZeroShift() {
        MassCalibrator.CalibrationStats stats = new MassCalibrator.CalibrationStats(0.0, 0.8, 250);
        Assert.assertTrue(stats.hasReliableStats());
        Assert.assertEquals(0.0, stats.getShiftPpm(), 0.0);
        Assert.assertEquals(0.8, stats.getRobustSigmaPpm(), 1e-12);
        Assert.assertEquals(250, stats.getConfidentPsmCount());
    }

    // ---- residualPpm() sign convention ----------------------------------

    @Test
    public void residualPpmPositiveWhenObservedGreater() {
        // observed > theoretical => positive residual (instrument reports a
        // mass slightly HIGHER than theoretical; calibrator will apply
        // peptideMass * (1 - shiftPpm * 1e-6) to remove the bias).
        double residual = MassCalibrator.residualPpmForTests(1001.0, 1000.0);
        Assert.assertTrue("Expected positive residual, got " + residual, residual > 0);
        Assert.assertEquals(1000.0, residual, 0.5); // roughly 1000 ppm
    }

    @Test
    public void residualPpmNegativeWhenObservedSmaller() {
        double residual = MassCalibrator.residualPpmForTests(999.0, 1000.0);
        Assert.assertTrue("Expected negative residual, got " + residual, residual < 0);
        Assert.assertEquals(-1000.0, residual, 0.5);
    }

    @Test
    public void residualPpmZeroWhenEqual() {
        Assert.assertEquals(0.0,
                MassCalibrator.residualPpmForTests(1000.0, 1000.0),
                1e-12);
    }

    @Test
    public void residualPpmFivePpmShift() {
        // A 5 ppm shift on a 1000 Da peptide is 0.005 Da.
        double observed = 1000.0 + 1000.0 * 5e-6;
        double residual = MassCalibrator.residualPpmForTests(observed, 1000.0);
        Assert.assertEquals(5.0, residual, 1e-6);
    }

    // ---- sampleEveryNth cap ---------------------------------------------

    @Test
    public void sampleEveryNthReturnsExpectedCount() {
        List<Integer> source = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            source.add(i);
        }
        List<Integer> sampled = MassCalibrator.sampleEveryNthForTests(source, 10, 500);
        Assert.assertEquals(10, sampled.size());
        // Sanity: first element is index 0, last is index 90.
        Assert.assertEquals(Integer.valueOf(0), sampled.get(0));
        Assert.assertEquals(Integer.valueOf(90), sampled.get(9));
    }

    @Test
    public void sampleEveryNthRespectsCap() {
        List<Integer> source = new ArrayList<>();
        for (int i = 0; i < 10000; i++) {
            source.add(i);
        }
        // Every 10th of 10k = 1000 candidates; cap at 500.
        List<Integer> sampled = MassCalibrator.sampleEveryNthForTests(source, 10, 500);
        Assert.assertEquals(500, sampled.size());
    }

    @Test
    public void sampleEveryNthEmpty() {
        Assert.assertTrue(MassCalibrator.sampleEveryNthForTests(Collections.emptyList(), 10, 500).isEmpty());
    }

    @Test
    public void sampleEveryNthSmallerThanStride() {
        List<Integer> source = Arrays.asList(0, 1, 2);
        List<Integer> sampled = MassCalibrator.sampleEveryNthForTests(source, 10, 500);
        // Only index 0 hits the stride.
        Assert.assertEquals(1, sampled.size());
        Assert.assertEquals(Integer.valueOf(0), sampled.get(0));
    }
}
