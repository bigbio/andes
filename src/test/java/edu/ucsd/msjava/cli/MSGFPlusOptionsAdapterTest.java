package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.params.ParamManager;
import org.junit.Assert;
import org.junit.Test;
import picocli.CommandLine;

/**
 * Phase 1 equivalence test: both the legacy
 * {@link ParamManager#parseParams(String[])} path and the new
 * picocli + {@link MSGFPlusOptionsAdapter} path must populate the
 * same {@link ParamManager} state for a representative CLI.
 *
 * If a future change drops a field from {@link MSGFPlusOptions} or the
 * adapter, this test catches the divergence before it reaches
 * {@code SearchParams.parse}.
 */
public class MSGFPlusOptionsAdapterTest {

    /** Canonical CLI a typical user passes to MS-GF+. */
    private static final String[] TYPICAL_CLI = {
            "-s", "src/test/resources/test.mgf",
            "-d", "src/test/resources/Tryp_Pig_Bov.fasta",
            "-t", "20ppm",
            "-ti", "-1,2",
            "-tda", "1",
            "-ntt", "2",
            "-thread", "4",
            "-minLength", "7",
            "-maxLength", "30",
            "-minCharge", "2",
            "-maxCharge", "4",
            "-n", "3",
            "-numMods", "2",
            "-addFeatures", "1",
            "-decoy", "XXX_",
    };

    @Test
    public void picocliPathPopulatesParamManagerEquivalentlyToLegacyPath() {
        ParamManager legacy = freshMSGFPlusParamManager();
        String legacyErr = legacy.parseParams(TYPICAL_CLI);
        Assert.assertNull("legacy parseParams returned error: " + legacyErr, legacyErr);

        ParamManager adapted = freshMSGFPlusParamManager();
        MSGFPlusOptions opts = new MSGFPlusOptions();
        new CommandLine(opts).parseArgs(TYPICAL_CLI);
        String adaptedErr = MSGFPlusOptionsAdapter.adapt(opts, adapted);
        Assert.assertNull("adapter returned error: " + adaptedErr, adaptedErr);

        // Compare every typed accessor that downstream SearchParams.parse reads.
        Assert.assertEquals(legacy.getDecoyProteinPrefix(), adapted.getDecoyProteinPrefix());
        Assert.assertEquals(legacy.getChargeCarrierMass(), adapted.getChargeCarrierMass(), 0.0);
        Assert.assertEquals(legacy.getNumTolerableTermini(), adapted.getNumTolerableTermini());
        Assert.assertEquals(legacy.getNumMatchesPerSpectrum(), adapted.getNumMatchesPerSpectrum());
        Assert.assertEquals(legacy.getTDA(), adapted.getTDA());
        Assert.assertEquals(legacy.getOutputAdditionalFeatures(), adapted.getOutputAdditionalFeatures());
        Assert.assertEquals(legacy.getMinPeptideLength(), adapted.getMinPeptideLength());
        Assert.assertEquals(legacy.getMaxPeptideLength(), adapted.getMaxPeptideLength());
        Assert.assertEquals(legacy.getMaxNumVariantsPerPeptide(), adapted.getMaxNumVariantsPerPeptide());
        Assert.assertEquals(legacy.getMinCharge(), adapted.getMinCharge());
        Assert.assertEquals(legacy.getMaxCharge(), adapted.getMaxCharge());
        Assert.assertEquals(legacy.getNumThreads(), adapted.getNumThreads());
        Assert.assertEquals(legacy.getOutputFormat(), adapted.getOutputFormat());
    }

    @Test
    public void picocliPathRejectsMissingRequiredFlags() {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        try {
            new CommandLine(opts).parseArgs(new String[] {"-t", "20ppm"});
            Assert.fail("expected picocli to reject CLI missing -s and -d");
        } catch (CommandLine.MissingParameterException expected) {
            // ok
        }
    }

    @Test
    public void picocliPathParsesAsymmetricTolerance() {
        ParamManager pm = freshMSGFPlusParamManager();
        String[] argv = {
                "-s", "src/test/resources/test.mgf",
                "-d", "src/test/resources/Tryp_Pig_Bov.fasta",
                "-t", "0.5Da,2.5Da",
        };
        MSGFPlusOptions opts = new MSGFPlusOptions();
        new CommandLine(opts).parseArgs(argv);
        String err = MSGFPlusOptionsAdapter.adapt(opts, pm);
        Assert.assertNull("adapter returned error on asymmetric tolerance: " + err, err);
        // Parity with legacy:
        ParamManager legacy = freshMSGFPlusParamManager();
        Assert.assertNull(legacy.parseParams(argv));
        Assert.assertEquals(legacy.getToleranceUnit(), pm.getToleranceUnit());
    }

    private static ParamManager freshMSGFPlusParamManager() {
        ParamManager pm = new ParamManager("MS-GF+", "test", "test", "test");
        pm.addMSGFPlusParams();
        return pm;
    }
}
