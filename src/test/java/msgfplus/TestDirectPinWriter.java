package msgfplus;

import edu.ucsd.msjava.msdbsearch.DatabaseMatch;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.mzid.DirectPinWriter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

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
                "lnDeltaSpecEValue", "matchedIonRatio",
                "Peptide", "Proteins"}) {
            Assert.assertTrue("Pin header should contain " + required + ": " + header,
                    header.contains(required));
        }
        // Column separator should be tab.
        Assert.assertTrue("Header should be tab-separated", header.contains("\t"));
        // The two extra features must come after the match-list features and before Peptide.
        int idxLast = header.indexOf("StdevRelErrorTop7");
        int idxLnDelta = header.indexOf("lnDeltaSpecEValue");
        int idxRatio = header.indexOf("matchedIonRatio");
        int idxPeptide = header.indexOf("Peptide");
        Assert.assertTrue("lnDeltaSpecEValue should come after StdevRelErrorTop7",
                idxLast > 0 && idxLnDelta > idxLast);
        Assert.assertTrue("matchedIonRatio should come after lnDeltaSpecEValue",
                idxRatio > idxLnDelta);
        Assert.assertTrue("Peptide should follow the extra features",
                idxPeptide > idxRatio);
    }

    // -----------------------------------------------------------------------
    // Helper tests for the two extra PSM-level features.
    // -----------------------------------------------------------------------

    @Test
    public void lnDeltaSpecEValueReturnsZeroForNonRank1() {
        Assert.assertEquals(0.0,
                DirectPinWriter.computeLnDeltaSpecEValue(2, 1e-10, 1e-5), 0.0);
        Assert.assertEquals(0.0,
                DirectPinWriter.computeLnDeltaSpecEValue(3, 1e-10, 1e-5), 0.0);
    }

    @Test
    public void lnDeltaSpecEValueReturnsLogRatioForRank1() {
        double rank1 = 1e-10;
        double rank2 = 1e-5;
        double expected = Math.log(rank1 / rank2); // negative: rank-1 more significant
        Assert.assertEquals(expected,
                DirectPinWriter.computeLnDeltaSpecEValue(1, rank1, rank2), 1e-12);
    }

    @Test
    public void lnDeltaSpecEValueIsZeroWhenRank2Missing() {
        Assert.assertEquals(0.0,
                DirectPinWriter.computeLnDeltaSpecEValue(1, 1e-10, Double.NaN), 0.0);
    }

    @Test
    public void lnDeltaSpecEValueIsZeroForNonPositiveInputs() {
        Assert.assertEquals(0.0,
                DirectPinWriter.computeLnDeltaSpecEValue(1, 0.0, 1e-5), 0.0);
        Assert.assertEquals(0.0,
                DirectPinWriter.computeLnDeltaSpecEValue(1, 1e-10, 0.0), 0.0);
        Assert.assertEquals(0.0,
                DirectPinWriter.computeLnDeltaSpecEValue(1, -1.0, 1e-5), 0.0);
    }

    @Test
    public void matchedIonRatioComputesNumMatchedOverPepLen() {
        Assert.assertEquals(0.5,
                DirectPinWriter.computeMatchedIonRatio("5", 10), 1e-12);
        Assert.assertEquals(1.0,
                DirectPinWriter.computeMatchedIonRatio("12", 12), 1e-12);
    }

    @Test
    public void matchedIonRatioHandlesMissingOrInvalidInput() {
        Assert.assertEquals(0.0,
                DirectPinWriter.computeMatchedIonRatio(null, 10), 0.0);
        Assert.assertEquals(0.0,
                DirectPinWriter.computeMatchedIonRatio("", 10), 0.0);
        Assert.assertEquals(0.0,
                DirectPinWriter.computeMatchedIonRatio("not-a-number", 10), 0.0);
    }

    @Test
    public void matchedIonRatioHandlesZeroOrNegativePepLen() {
        Assert.assertEquals(0.0,
                DirectPinWriter.computeMatchedIonRatio("5", 0), 0.0);
        Assert.assertEquals(0.0,
                DirectPinWriter.computeMatchedIonRatio("5", -1), 0.0);
    }

    @Test
    public void findRank2ReturnsDistinctNextBestSpecEValue() {
        // matchList is ordered worst-to-best: last element is rank-1.
        List<DatabaseMatch> matches = new ArrayList<>();
        matches.add(newMatch(1e-5));  // rank 3
        matches.add(newMatch(1e-8));  // rank 2
        matches.add(newMatch(1e-10)); // rank 1

        Assert.assertEquals(1e-8,
                DirectPinWriter.findRank2SpecEValue(matches, 0), 0.0);
    }

    @Test
    public void findRank2SkipsTiesWithRank1() {
        // Rank-1 and the next entry share a SpecEValue (tied rank-1 group);
        // rank-2 is the first *distinct* value below them.
        List<DatabaseMatch> matches = new ArrayList<>();
        matches.add(newMatch(1e-5));  // rank 3
        matches.add(newMatch(1e-10)); // rank 1 (tie)
        matches.add(newMatch(1e-10)); // rank 1 (tie)

        Assert.assertEquals(1e-5,
                DirectPinWriter.findRank2SpecEValue(matches, 0), 0.0);
    }

    @Test
    public void findRank2ReturnsNaNWhenOnlyOneRank() {
        List<DatabaseMatch> matches = new ArrayList<>();
        matches.add(newMatch(1e-10));
        Assert.assertTrue(
                Double.isNaN(DirectPinWriter.findRank2SpecEValue(matches, 0)));
    }

    @Test
    public void findRank2ReturnsNaNForEmptyList() {
        Assert.assertTrue(
                Double.isNaN(DirectPinWriter.findRank2SpecEValue(Collections.emptyList(), 0)));
    }

    private static DatabaseMatch newMatch(double specEValue) {
        DatabaseMatch m = new DatabaseMatch(0, (byte) 10, 100, 1000f, 1000, 2,
                "PEPTIDER", new ActivationMethod[]{ActivationMethod.CID});
        m.setSpecProb(specEValue);
        // DeNovoScore defaults to 0; test uses minDeNovoScore=0 so all matches qualify.
        return m;
    }
}
