package msgfplus;

import edu.ucsd.msjava.cli.MSGFPlusOptions;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParams.PrecursorCalMode;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.msutil.DBSearchIOFiles;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;

/**
 * Tests for the CLI scaffolding that Achievement B (two-pass precursor mass
 * calibration) layers on top of existing search parameters.
 * <p>
 * These tests pin:
 * <ol>
 *     <li>The {@code -precursorCal} flag parses cleanly with
 *         {@code auto}/{@code on}/{@code off} (case-insensitive) and defaults
 *         to {@code auto}.</li>
 *     <li>{@link DBSearchIOFiles#getPrecursorMassShiftPpm()} defaults to
 *         {@code 0.0} and survives a round-trip through its setter.</li>
 *     <li>Unknown values fall back to {@link PrecursorCalMode#AUTO} so that
 *         downstream code always has a sensible default.</li>
 * </ol>
 */
public class TestPrecursorCalScaffolding {

    private MSGFPlusOptions buildOpts() throws URISyntaxException {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        opts.configFile   = new File(SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI());
        opts.spectrumFile = new File(SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI());
        opts.databaseFile = new File(SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI());
        return opts;
    }

    @Test
    public void precursorCalDefaultIsAuto() throws URISyntaxException {
        MSGFPlusOptions opts = buildOpts();
        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(opts));
        Assert.assertEquals("Default -precursorCal should be AUTO",
                PrecursorCalMode.AUTO, params.getPrecursorCalMode());
    }

    @Test
    public void precursorCalOnIsParsed() throws URISyntaxException {
        MSGFPlusOptions opts = buildOpts();
        opts.precursorCalMode = "on";
        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(opts));
        Assert.assertEquals(PrecursorCalMode.ON, params.getPrecursorCalMode());
    }

    @Test
    public void precursorCalOffIsParsed() throws URISyntaxException {
        MSGFPlusOptions opts = buildOpts();
        opts.precursorCalMode = "off";
        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(opts));
        Assert.assertEquals(PrecursorCalMode.OFF, params.getPrecursorCalMode());
    }

    @Test
    public void precursorCalIsCaseInsensitive() throws URISyntaxException {
        MSGFPlusOptions opts = buildOpts();
        opts.precursorCalMode = "OFF";
        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(opts));
        Assert.assertEquals(PrecursorCalMode.OFF, params.getPrecursorCalMode());
    }

    @Test
    public void unknownPrecursorCalValueFallsBackToAuto() {
        // Unit-level contract: unknown strings must not crash the search path;
        // instead they silently fall back to AUTO.
        Assert.assertEquals(PrecursorCalMode.AUTO, PrecursorCalMode.fromString("bogus"));
        Assert.assertEquals(PrecursorCalMode.AUTO, PrecursorCalMode.fromString(null));
        Assert.assertEquals(PrecursorCalMode.AUTO, PrecursorCalMode.fromString(""));
    }

    @Test
    public void dbSearchIOFilesShiftDefaultsToZero() {
        DBSearchIOFiles ioFiles = new DBSearchIOFiles(
                new File("dummy.mzML"), SpecFileFormat.MZML, new File("dummy.mzid"));
        Assert.assertEquals("Default shift should be 0.0 ppm",
                0.0, ioFiles.getPrecursorMassShiftPpm(), 0.0);
    }

    @Test
    public void dbSearchIOFilesShiftRoundTrips() {
        DBSearchIOFiles ioFiles = new DBSearchIOFiles(
                new File("dummy.mzML"), SpecFileFormat.MZML, new File("dummy.mzid"));
        ioFiles.setPrecursorMassShiftPpm(4.2);
        Assert.assertEquals(4.2, ioFiles.getPrecursorMassShiftPpm(), 1e-12);
    }
}
