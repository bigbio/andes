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
}
