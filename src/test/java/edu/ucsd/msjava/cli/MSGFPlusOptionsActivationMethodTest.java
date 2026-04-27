package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.msutil.ActivationMethod;
import org.junit.Assert;
import org.junit.Test;

/**
 * Pins the {@code -m} ID -> {@link ActivationMethod} mapping. The legacy
 * dispatch went through the registry order (ASWRITTEN, CID, ETD, HCD, FUSION,
 * UVPD) with {@code FUSION} hidden by {@code addFragMethodParam(...,
 * doNotAddMergeMode=true)}, which shifted {@code UVPD} from registry slot 5
 * to the user-facing index 4. The Phase 4c rewrite originally hardcoded only
 * 0..3 and silently dropped UVPD; this test guards against regressing it
 * again.
 */
public class MSGFPlusOptionsActivationMethodTest {

    @Test
    public void defaultIsAsWritten() {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        Assert.assertSame(ActivationMethod.ASWRITTEN, opts.effectiveActivationMethod());
    }

    @Test
    public void mapsAllSupportedIndices() {
        Assert.assertSame(ActivationMethod.ASWRITTEN, withFragMethodId(0).effectiveActivationMethod());
        Assert.assertSame(ActivationMethod.CID,       withFragMethodId(1).effectiveActivationMethod());
        Assert.assertSame(ActivationMethod.ETD,       withFragMethodId(2).effectiveActivationMethod());
        Assert.assertSame(ActivationMethod.HCD,       withFragMethodId(3).effectiveActivationMethod());
        Assert.assertSame(ActivationMethod.UVPD,      withFragMethodId(4).effectiveActivationMethod());
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsOutOfRangeIndex() {
        withFragMethodId(5).effectiveActivationMethod();
    }

    private static MSGFPlusOptions withFragMethodId(int id) {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        opts.fragMethodId = id;
        return opts;
    }
}
