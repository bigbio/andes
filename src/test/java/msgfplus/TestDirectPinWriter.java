package msgfplus;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;

/**
 * Shape tests for the Percolator {@code .pin} output path (Q7).
 *
 * These exercise the CLI/flag plumbing and the header emitted by
 * {@link edu.ucsd.msjava.mzid.DirectPinWriter}. A full end-to-end
 * search-to-pin run is exercised by the integration tests under
 * {@code src/test/resources/} when external spectra are available;
 * here we focus on the parts we can verify without running the
 * search engine.
 */
public class TestDirectPinWriter {

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
    public void pinOutputFormatFlagIsAccepted() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        String err = manager.getParameter("outputFormat").parse("3");
        Assert.assertNull("parse('pin') should succeed but returned: " + err, err);
    }

    @Test
    public void writePinGetterReflectsOutputFormat() throws URISyntaxException {
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("outputFormat").parse("3"));

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(manager));

        Assert.assertTrue("writePin() should be true when outputFormat=pin", params.writePin());
        Assert.assertFalse("writeMzid() should be false when outputFormat=pin", params.writeMzid());
        Assert.assertFalse("writeTsv() should be false when outputFormat=pin", params.writeTsv());
    }

    @Test
    public void allOutputFormatEnumIndicesAreAccepted() throws URISyntaxException {
        // Sanity guard that adding the "pin" (index 3) entry didn't shift any
        // existing index. 0=mzid, 1=tsv, 2=both, 3=pin.
        for (String value : new String[]{"0", "1", "2", "3"}) {
            ParamManager manager = buildParamManager();
            String err = manager.getParameter("outputFormat").parse(value);
            Assert.assertNull("parse('" + value + "') should succeed but returned: " + err, err);
        }
    }

    @Test
    public void pinHeaderColumnsIncludeRequiredPercolatorFields() throws Exception {
        // Build a minimal result list so DirectPinWriter can emit a header.
        // We don't need real matches; an empty resultList still produces the
        // header row, which is what we're testing.
        ParamManager manager = buildParamManager();
        Assert.assertNull(manager.getParameter("outputFormat").parse("3"));

        SearchParams params = new SearchParams();
        Assert.assertNull(params.parse(manager));

        // DirectPinWriter needs a CompactSuffixArray and SpectraAccessor; we
        // can't construct those without running through BuildSA and loading
        // spectra. Instead, we verify the header shape indirectly by using
        // the Writer's internal header format via a small probe.
        //
        // Specifically: invoke DirectPinWriter via reflection on a stub output
        // stream. We assert the header line contains the Percolator-required
        // column names.
        java.lang.reflect.Method writeHeader =
                edu.ucsd.msjava.mzid.DirectPinWriter.class.getDeclaredMethod(
                        "writeHeader", java.io.PrintStream.class, int.class, int.class);
        writeHeader.setAccessible(true);

        // Build a DirectPinWriter with null sa/specAcc — header path doesn't
        // touch them. If the constructor starts using them, this test needs
        // to evolve; for now it's a pure shape check.
        java.lang.reflect.Constructor<?> ctor = edu.ucsd.msjava.mzid.DirectPinWriter.class
                .getDeclaredConstructor(
                        SearchParams.class,
                        edu.ucsd.msjava.msutil.AminoAcidSet.class,
                        edu.ucsd.msjava.msdbsearch.CompactSuffixArray.class,
                        edu.ucsd.msjava.msutil.SpectraAccessor.class,
                        int.class);
        Object writer = ctor.newInstance(params, params.getAASet(), null, null, 0);

        File tmp = File.createTempFile("msgfplus-pin-header-", ".pin");
        tmp.deleteOnExit();
        try (java.io.PrintStream ps = new java.io.PrintStream(new java.io.FileOutputStream(tmp))) {
            writeHeader.invoke(writer, ps, 2, 4); // minCharge=2, maxCharge=4
        }
        String header = new String(Files.readAllBytes(tmp.toPath()), java.nio.charset.StandardCharsets.UTF_8).trim();
        for (String required : new String[]{
                "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass",
                "RawScore", "DeNovoScore", "lnSpecEValue", "lnEValue", "IsotopeError",
                "PepLen", "dM", "absdM",
                "Charge2", "Charge3", "Charge4",
                "NumMatchedMainIons", "ExplainedIonCurrentRatio",
                "Peptide", "Proteins"}) {
            Assert.assertTrue("Pin header should contain " + required + ": " + header,
                    header.contains(required));
        }
        // Column separator should be tab.
        Assert.assertTrue("Header should be tab-separated", header.contains("\t"));
    }
}
