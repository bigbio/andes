package msgfplus;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParams.PrecursorCalMode;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.msutil.DBSearchIOFiles;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
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

    private ParamManager buildParamManager() throws URISyntaxException {
        ParamManager manager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE,
                "java -Xmx3500M -jar MSGFPlus.jar");
        manager.addMSGFPlusParams();

        URI paramUri = SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI();
        manager.getParameter("conf").parse(new File(paramUri).getAbsolutePath());

        URI specUri = SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI();
        manager.getParameter("s").parse(new File(specUri).getAbsolutePath());

        URI dbUri = SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI();
        manager.getParameter("d").parse(new File(dbUri).getAbsolutePath());
        return manager;
    }

    @Test
    public void precursorCalDefaultIsAuto() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals("Default -precursorCal should be AUTO",
                PrecursorCalMode.AUTO, params.getPrecursorCalMode());
    }

    @Test
    public void precursorCalOnIsParsed() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("precursorCal").parse("on"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals(PrecursorCalMode.ON, params.getPrecursorCalMode());
    }

    @Test
    public void precursorCalOffIsParsed() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("precursorCal").parse("off"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals(PrecursorCalMode.OFF, params.getPrecursorCalMode());
    }

    @Test
    public void precursorCalIsCaseInsensitive() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("precursorCal").parse("OFF"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
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
