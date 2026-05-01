package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.Spectrum;
import org.junit.Assert;
import org.junit.Test;

import java.util.Collections;

public class TestScoredSpectraMapIsolation {

    @Test
    public void defaultPathMutatesOriginalSpectrumCharge() {
        ScoredSpectraMap map = new ScoredSpectraMap(
                null,
                Collections.emptyList(),
                new Tolerance(10f, true),
                new Tolerance(10f, true),
                0,
                0,
                null,
                false,
                false);
        Spectrum original = new Spectrum(500f, 2, 100f);

        Spectrum prepared = map.prepareSpectrumForScoring(original, 3);

        Assert.assertSame(original, prepared);
        Assert.assertEquals(3, original.getCharge());
    }

    @Test
    public void isolatedPathClonesSpectrumBeforeChangingCharge() {
        ScoredSpectraMap map = new ScoredSpectraMap(
                null,
                Collections.emptyList(),
                new Tolerance(10f, true),
                new Tolerance(10f, true),
                0,
                0,
                null,
                false,
                false).isolateSpectrumState();
        Spectrum original = new Spectrum(500f, 2, 100f);

        Spectrum prepared = map.prepareSpectrumForScoring(original, 3);

        Assert.assertNotSame(original, prepared);
        Assert.assertEquals(2, original.getCharge());
        Assert.assertEquals(3, prepared.getCharge());
    }
}
