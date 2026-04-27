package msgfplus;

import edu.ucsd.msjava.cli.MSGFPlusOptions;
import edu.ucsd.msjava.cli.OutputFormat;
import edu.ucsd.msjava.msdbsearch.DatabaseMatch;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.output.DirectPinWriter;
import picocli.CommandLine;
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
 * {@link edu.ucsd.msjava.output.DirectPinWriter}. A full end-to-end
 * search-to-pin run is exercised by the integration tests under
 * {@code src/test/resources/} when external spectra are available;
 * here we focus on the parts we can verify without running the
 * search engine.
 */
public class TestDirectPinWriter {

    private MSGFPlusOptions buildOpts() throws URISyntaxException {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        opts.configFile   = new File(SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI());
        opts.spectrumFile = new File(SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI());
        opts.databaseFile = new File(SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI());
        return opts;
    }

    @Test
    public void pinOutputFormatFlagIsAccepted() throws URISyntaxException {
        MSGFPlusOptions opts = buildOpts();
        opts.outputFormat = OutputFormat.PIN;
        Assert.assertEquals(OutputFormat.PIN, opts.effectiveOutputFormat());
    }

    @Test
    public void writePinGetterReflectsOutputFormat() throws URISyntaxException {
        MSGFPlusOptions opts = buildOpts();
        opts.outputFormat = OutputFormat.PIN;

        SearchParams params = new SearchParams();
        Assert.assertNull("SearchParams.parse should succeed", params.parse(opts));

        Assert.assertTrue("writePin() should be true when outputFormat=pin", params.writePin());
        Assert.assertFalse("writeTsv() should be false when outputFormat=pin", params.writeTsv());
    }

    @Test
    public void outputFormatAcceptsOnlyPinAndTsv() throws URISyntaxException {
        // Picocli matches enum values case-insensitively per the @Command setting.
        for (String value : new String[]{"pin", "PIN", "Pin", "tsv", "TSV", "Tsv"}) {
            MSGFPlusOptions opts = new MSGFPlusOptions();
            MSGFPlusOptions.commandLine(opts).parseArgs("-outputFormat", value);
            Assert.assertNotNull("'" + value + "' should parse to a valid OutputFormat", opts.outputFormat);
        }
        // Numeric forms (0/1) and removed legacy values (mzid, both, 2, 3) are
        // intentionally rejected -- the typed enum is part of the consistency
        // sweep called out in the parameter-modernization cleanup.
        for (String value : new String[]{"0", "1", "2", "3", "mzid", "both", ""}) {
            MSGFPlusOptions opts = new MSGFPlusOptions();
            try {
                MSGFPlusOptions.commandLine(opts).parseArgs("-outputFormat", value);
                Assert.fail("'" + value + "' should be rejected by picocli enum matching");
            } catch (CommandLine.ParameterException expected) {
                // ok
            }
        }
    }

    @Test
    public void pinHeaderColumnsIncludeRequiredPercolatorFields() throws Exception {
        MSGFPlusOptions opts = buildOpts();
        opts.outputFormat = OutputFormat.PIN;

        SearchParams params = new SearchParams();
        Assert.assertNull(params.parse(opts));

        // DirectPinWriter needs a CompactSuffixArray and SpectraAccessor; we
        // can't construct those without running through BuildSA and loading
        // spectra. Instead, we verify the header shape indirectly by using
        // the Writer's internal header format via a small probe.
        //
        // Specifically: invoke DirectPinWriter via reflection on a stub output
        // stream. We assert the header line contains the Percolator-required
        // column names.
        java.lang.reflect.Method writeHeader =
                edu.ucsd.msjava.output.DirectPinWriter.class.getDeclaredMethod(
                        "writeHeader", java.io.PrintStream.class, int.class, int.class);
        writeHeader.setAccessible(true);

        // Build a DirectPinWriter with null sa/specAcc — header path doesn't
        // touch them. If the constructor starts using them, this test needs
        // to evolve; for now it's a pure shape check.
        java.lang.reflect.Constructor<?> ctor = edu.ucsd.msjava.output.DirectPinWriter.class
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
                "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass",
                "RawScore", "DeNovoScore", "lnSpecEValue", "lnEValue", "isotope_error",
                "peplen", "dm", "absdm",
                "charge2", "charge3", "charge4",
                "enzN", "enzC", "enzInt",
                "NumMatchedMainIons", "longest_b", "longest_y", "longest_y_pct",
                "ExplainedIonCurrentRatio",
                "lnDeltaSpecEValue", "matchedIonRatio",
                "Peptide", "Proteins"}) {
            Assert.assertTrue("Pin header should contain " + required + ": " + header,
                    header.contains(required));
        }
        // Renamed columns must not appear under their legacy case-sensitive names.
        // We use tab-delimited matches to avoid accidental substring hits
        // (e.g., "dM" would otherwise trivially appear inside "ExpMass").
        for (String legacy : new String[]{"PepLen", "Charge2", "Charge3", "Charge4",
                "\tdM\t", "\tabsdM\t", "IsotopeError"}) {
            String probe = legacy.startsWith("\t") ? legacy : "\t" + legacy;
            Assert.assertFalse("Pin header should NOT contain legacy name " + legacy + ": " + header,
                    ("\t" + header + "\t").contains(probe));
        }
        // Column separator should be tab.
        Assert.assertTrue("Header should be tab-separated", header.contains("\t"));
        // mass must come right after CalcMass.
        Assert.assertTrue("mass should appear right after CalcMass: " + header,
                header.contains("\tCalcMass\tmass\t"));
        // enzN/enzC/enzInt must sit between the charge block and NumMatchedMainIons.
        int idxChargeLast = header.indexOf("charge4");
        int idxEnzN = header.indexOf("enzN");
        int idxEnzC = header.indexOf("enzC");
        int idxEnzInt = header.indexOf("enzInt");
        int idxNumMatched = header.indexOf("NumMatchedMainIons");
        Assert.assertTrue("enzN should come after the charge block",
                idxChargeLast > 0 && idxEnzN > idxChargeLast);
        Assert.assertTrue("enzC should come after enzN", idxEnzC > idxEnzN);
        Assert.assertTrue("enzInt should come after enzC", idxEnzInt > idxEnzC);
        Assert.assertTrue("NumMatchedMainIons should come after enzInt",
                idxNumMatched > idxEnzInt);
        // Ion-series run-length features must follow NumMatchedMainIons and precede
        // the ExplainedIonCurrent* ratios (they're part of the ion-structure block).
        int idxLongestB = header.indexOf("longest_b");
        int idxLongestY = header.indexOf("longest_y\t"); // tab-anchor to avoid matching longest_y_pct
        int idxLongestYPct = header.indexOf("longest_y_pct");
        int idxEIC = header.indexOf("ExplainedIonCurrentRatio");
        Assert.assertTrue("longest_b should come after NumMatchedMainIons",
                idxLongestB > idxNumMatched);
        Assert.assertTrue("longest_y should come after longest_b",
                idxLongestY > idxLongestB);
        Assert.assertTrue("longest_y_pct should come after longest_y",
                idxLongestYPct > idxLongestY);
        Assert.assertTrue("ExplainedIonCurrentRatio should come after longest_y_pct",
                idxEIC > idxLongestYPct);
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
    // Enzymatic-boundary helpers (mirror OpenMS PercolatorInfile::isEnz_).
    // -----------------------------------------------------------------------

    @Test
    public void enzymaticBoundaryTrypsinRulesMatchOpenMS() {
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('K', 'A', "trypsin"));
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('R', 'A', "trypsin"));
        Assert.assertFalse("KP is not a trypsin cleavage site",
                DirectPinWriter.isEnzymaticBoundary('K', 'P', "trypsin"));
        Assert.assertFalse("RP is not a trypsin cleavage site",
                DirectPinWriter.isEnzymaticBoundary('R', 'P', "trypsin"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('A', 'K', "trypsin"));
        Assert.assertTrue("N-terminal protein boundary is enzymatic",
                DirectPinWriter.isEnzymaticBoundary('-', 'A', "trypsin"));
        Assert.assertTrue("C-terminal protein boundary is enzymatic",
                DirectPinWriter.isEnzymaticBoundary('A', '-', "trypsin"));
    }

    @Test
    public void enzymaticBoundaryLysNLysCAspNGluCArgCMatchOpenMS() {
        // lys-c: cleave after K (unless c == P).
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('K', 'A', "lys-c"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('K', 'P', "lys-c"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('R', 'A', "lys-c"));
        // lys-n: cleave before K.
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('A', 'K', "lys-n"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('K', 'A', "lys-n"));
        // arg-c: cleave after R (unless c == P).
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('R', 'A', "arg-c"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('R', 'P', "arg-c"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('K', 'A', "arg-c"));
        // asp-n: cleave before D.
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('A', 'D', "asp-n"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('D', 'A', "asp-n"));
        // glu-c: cleave after E (unless c == P).
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('E', 'A', "glu-c"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('E', 'P', "glu-c"));
        Assert.assertFalse(DirectPinWriter.isEnzymaticBoundary('A', 'E', "glu-c"));
    }

    @Test
    public void enzymaticBoundaryUnknownEnzymeReturnsTrue() {
        // OpenMS default falls through to `true` when the enzyme name is unknown.
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('A', 'B', ""));
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('A', 'B', null));
        Assert.assertTrue(DirectPinWriter.isEnzymaticBoundary('A', 'B', "no-such-enzyme"));
    }

    @Test
    public void openMsEnzymeNameMapsKnownSingletons() {
        Assert.assertEquals("trypsin", DirectPinWriter.openMsEnzymeName(Enzyme.TRYPSIN));
        Assert.assertEquals("chymotrypsin", DirectPinWriter.openMsEnzymeName(Enzyme.CHYMOTRYPSIN));
        Assert.assertEquals("lys-c", DirectPinWriter.openMsEnzymeName(Enzyme.LysC));
        Assert.assertEquals("lys-n", DirectPinWriter.openMsEnzymeName(Enzyme.LysN));
        Assert.assertEquals("arg-c", DirectPinWriter.openMsEnzymeName(Enzyme.ArgC));
        Assert.assertEquals("asp-n", DirectPinWriter.openMsEnzymeName(Enzyme.AspN));
        Assert.assertEquals("glu-c", DirectPinWriter.openMsEnzymeName(Enzyme.GluC));
        Assert.assertEquals("", DirectPinWriter.openMsEnzymeName(null));
        Assert.assertEquals("", DirectPinWriter.openMsEnzymeName(Enzyme.UnspecificCleavage));
        Assert.assertEquals("", DirectPinWriter.openMsEnzymeName(Enzyme.NoCleavage));
        Assert.assertEquals("", DirectPinWriter.openMsEnzymeName(Enzyme.ALP));
        Assert.assertEquals("", DirectPinWriter.openMsEnzymeName(Enzyme.TrypsinPlusC));
    }

    @Test
    public void countInternalEnzymaticTrypsin() {
        // AKAKR, trypsin: i=1 (A,K)=false; i=2 (K,A)=true; i=3 (A,K)=false; i=4 (K,R)=true → 2.
        Assert.assertEquals(2, DirectPinWriter.countInternalEnzymatic("AKAKR", "trypsin"));
        // KP rule: RKPK → i=1 (R,K)=true; i=2 (K,P)=false (KP); i=3 (P,K)=false → 1.
        Assert.assertEquals(1, DirectPinWriter.countInternalEnzymatic("RKPK", "trypsin"));
    }

    @Test
    public void countInternalEnzymaticUnspecificEnzymeCountsEveryInterior() {
        // OpenMS default-true behavior: every interior boundary counts, giving peplen - 1.
        Assert.assertEquals(6, DirectPinWriter.countInternalEnzymatic("PEPTIDE", ""));
        Assert.assertEquals(6, DirectPinWriter.countInternalEnzymatic("PEPTIDE", null));
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
    public void sanitizeFeatureValueHandlesNaNAndInfinity() {
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue(null));
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue(""));
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue("NaN"));
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue("nan"));
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue("Infinity"));
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue("-Infinity"));
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue("Inf"));
        Assert.assertEquals("0", DirectPinWriter.sanitizeFeatureValue("-Inf"));
        Assert.assertEquals("1.5", DirectPinWriter.sanitizeFeatureValue("1.5"));
        Assert.assertEquals("-0.003", DirectPinWriter.sanitizeFeatureValue("-0.003"));
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
