package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Regression tests for {@link MSGFPlusOptions#applyConfigFile} and the
 * downstream {@link SearchParams#parse} path.
 *
 * Pins the {@code CustomAA=} crash that was caught in code review: the
 * legacy hashtable-based config-file reader passed bare values to
 * {@code AminoAcidSet.parseConfigEntry}, but the modernized adapter
 * briefly re-prepended {@code "CustomAA="} which {@code parseConfigEntry}
 * does not strip — every {@code -conf} invocation containing a
 * {@code CustomAA=} line crashed via {@code System.exit(-1)}.
 */
public class MSGFPlusOptionsConfigFileTest {

    @Test
    public void configFileWithCustomAAParsesWithoutCrashing() throws IOException, URISyntaxException {
        // Build a minimal config file with the documented CustomAA= form.
        Path tmpDir = Files.createTempDirectory("msgfplus-customaa-");
        Path conf = tmpDir.resolve("with_custom_aa.txt");
        Files.write(conf, ("# Regression for the CustomAA= prefix bug\n"
                + "CustomAA=C3H5NO, U, custom, U, Selenocysteine\n"
                + "MinPepLength=7\n").getBytes(StandardCharsets.UTF_8));

        URI specUri = MSGFPlusOptionsConfigFileTest.class.getClassLoader()
                .getResource("test.mgf").toURI();
        URI dbUri = MSGFPlusOptionsConfigFileTest.class.getClassLoader()
                .getResource("Tryp_Pig_Bov.fasta").toURI();

        MSGFPlusOptions opts = new MSGFPlusOptions();
        opts.configFile = conf.toFile();
        opts.spectrumFile = new File(specUri);
        opts.databaseFile = new File(dbUri);

        SearchParams params = new SearchParams();
        String err = params.parse(opts);
        Assert.assertNull("SearchParams.parse must not crash on a config file with CustomAA= entries: " + err, err);

        // The custom AA list should reach opts.customAAs and be honored downstream.
        Assert.assertEquals(1, opts.customAAs.size());
        Assert.assertEquals("config-file MinPepLength=7 should win over the default of 6",
                7, opts.effectiveMinPeptideLength());

        // Cleanup.
        Files.deleteIfExists(conf);
        Files.deleteIfExists(tmpDir);
    }

    /**
     * Regression for the case-insensitive config-key match. The legacy
     * {@code ParamManager.parseConfigParamFile} matched names with
     * {@code equalsIgnoreCase}; the Phase 4c switch was exact-case so
     * {@code minCharge=} / {@code maxCharge=} from the test fixture
     * silently fell back to defaults instead of overriding them.
     */
    @Test
    public void configFileKeysAreMatchedCaseInsensitively() throws IOException {
        Path tmpDir = Files.createTempDirectory("msgfplus-caseinsens-");
        Path conf = tmpDir.resolve("mixed_case.txt");
        // Mix of canonical, lowercased-first-letter, and ALLCAPS forms.
        Files.write(conf, ("MinPepLength=8\n"
                + "maxpepLength=42\n"
                + "MINCHARGE=3\n"
                + "maxcharge=7\n"
                + "TDA=1\n").getBytes(StandardCharsets.UTF_8));

        MSGFPlusOptions opts = new MSGFPlusOptions();
        Assert.assertNull(opts.applyConfigFile(conf.toFile()));

        Assert.assertEquals(8,  opts.effectiveMinPeptideLength());
        Assert.assertEquals(42, opts.effectiveMaxPeptideLength());
        Assert.assertEquals(3,  opts.effectiveMinCharge());
        Assert.assertEquals(7,  opts.effectiveMaxCharge());
        Assert.assertEquals(1,  opts.effectiveTdaStrategy());

        Files.deleteIfExists(conf);
        Files.deleteIfExists(tmpDir);
    }

    /**
     * Pin the numeric/enum range validation that the legacy
     * {@code IntParameter.minValue}/{@code maxValue} machinery used to
     * enforce. After Phase 4c those checks initially disappeared; restoring
     * them ensures invalid CLI input produces a clean error string instead
     * of a stack trace from a downstream resolver.
     */
    @Test
    public void validateRejectsOutOfRangeFlags() {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        opts.spectrumFile = new File("anything.mgf");
        opts.databaseFile = new File("anything.fasta");

        opts.numThreads = 0;
        Assert.assertNotNull("numThreads=0 must be rejected", opts.validate());
        opts.numThreads = null;

        opts.fragMethodId = 99;
        Assert.assertNotNull("-m 99 must be rejected with a user-facing error", opts.validate());
        opts.fragMethodId = null;

        opts.numTolerableTermini = 5;
        Assert.assertNotNull("-ntt 5 must be rejected (valid 0..2)", opts.validate());
        opts.numTolerableTermini = null;

        opts.tdaStrategy = 2;
        Assert.assertNotNull("-tda 2 must be rejected (valid 0..1)", opts.validate());
        opts.tdaStrategy = null;

        // A clean invocation passes.
        Assert.assertNull(opts.validate());
    }
}
