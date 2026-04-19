package msgfplus;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParams.FragmentIndexMode;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;

/**
 * Tests for the CLI scaffolding that Phase 3 (fragment-index Tier-1 search)
 * layers on top of existing search parameters.
 * <p>
 * These tests pin:
 * <ol>
 *     <li>The {@code -useFragmentIndex} flag defaults to {@code off} when
 *         absent.</li>
 *     <li>{@code off}, {@code on}, and {@code compare} all parse (case-
 *         insensitively) into the corresponding
 *         {@link FragmentIndexMode}.</li>
 *     <li>Unknown values cause {@link SearchParams#parse(ParamManager)} to
 *         return a non-null error message.</li>
 * </ol>
 * <p>
 * Mirrors {@link TestPrecursorCalScaffolding} — this is pure param-system
 * plumbing; no search-time or scoring code is exercised.
 */
public class TestFragmentIndexModeScaffolding {

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
    public void fragmentIndexDefaultIsOff() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals("Default -useFragmentIndex should be OFF",
                FragmentIndexMode.OFF, params.getFragmentIndexMode());
    }

    @Test
    public void fragmentIndexOffIsParsed() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("useFragmentIndex").parse("off"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals(FragmentIndexMode.OFF, params.getFragmentIndexMode());
    }

    @Test
    public void fragmentIndexOnIsParsed() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("useFragmentIndex").parse("on"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals(FragmentIndexMode.ON, params.getFragmentIndexMode());
    }

    @Test
    public void fragmentIndexCompareIsParsed() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("useFragmentIndex").parse("compare"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals(FragmentIndexMode.COMPARE, params.getFragmentIndexMode());
    }

    @Test
    public void fragmentIndexIsCaseInsensitive() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("useFragmentIndex").parse("Compare"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));
        Assert.assertEquals(FragmentIndexMode.COMPARE, params.getFragmentIndexMode());
    }

    @Test
    public void invalidFragmentIndexValueReturnsError() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        // StringParameter.parse always accepts; validation happens in SearchParams.parse.
        Assert.assertNull(manager.getParameter("useFragmentIndex").parse("banana"));

        SearchParams params = new SearchParams();
        String error = params.parse(manager);
        Assert.assertNotNull("SearchParams.parse should reject unknown -useFragmentIndex values",
                error);
        Assert.assertTrue("Error message should mention useFragmentIndex",
                error.toLowerCase().contains("usefragmentindex"));
    }
}
