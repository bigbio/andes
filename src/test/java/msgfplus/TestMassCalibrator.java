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
