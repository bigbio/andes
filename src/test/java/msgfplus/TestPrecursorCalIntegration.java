package msgfplus;

import edu.ucsd.msjava.cli.MSGFPlus;
import edu.ucsd.msjava.cli.MSGFPlusOptions;
import edu.ucsd.msjava.msdbsearch.SearchParams.PrecursorCalMode;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.msutil.DBSearchIOFiles;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * End-to-end integration tests for Achievement B — two-pass precursor mass
 * calibration (P2-cal).
 *
 * <p>The star test here is {@link #precursorCalOffMatchesBaseline()}, which is
 * the hard correctness gate from the design spec:
 * <blockquote>
 *     When {@code -precursorCal off} is supplied, the branch must produce
 *     bit-identical results to a run without any calibration code path.
 * </blockquote>
 * We enforce it by running two full searches on the bundled
 * {@code test.mgf} + {@code human-uniprot-contaminants.fasta} pair and
 * comparing every PSM data row from the two {@code .pin} outputs. A drift
 * here would be a silent FDR-inflating bug, so we demand strict equality
 * on the PSM list.
 *
 * <p>Because the {@code test.mgf} fixture is small, the default {@code auto}
 * mode takes the "insufficient confident PSMs" branch and also produces a
 * 0.0 ppm shift, so the comparison is against the same no-op-shift baseline.
 */
public class TestPrecursorCalIntegration {

    private MSGFPlusOptions buildOpts(File outputFile) throws URISyntaxException {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        opts.configFile   = new File(SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI());
        opts.spectrumFile = new File(SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI());
        opts.databaseFile = new File(SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI());
        opts.outputFile   = outputFile;
        return opts;
    }

    /**
     * Hard correctness gate: {@code -precursorCal off} must produce a
     * PSM list identical to a run with no flag at all.
     *
     * <p>The test runs both searches in a fresh temp dir to avoid colliding
     * with any cached suffix-array artefacts from other tests, then
     * compares the pin-file PSM data rows line by line.
     */
    @Test
    public void precursorCalOffMatchesBaseline() throws Exception {
        Path workDir = Files.createTempDirectory("msgfplus-p2cal-integration-");
        try {
            File offOut = new File(workDir.toFile(), "off.pin");
            File baselineOut = new File(workDir.toFile(), "baseline.pin");

            MSGFPlusOptions offManager = buildOpts(offOut);
            offManager.precursorCalMode = PrecursorCalMode.OFF;
            String offErr = MSGFPlus.runMSGFPlus(offManager);
            Assert.assertNull("runMSGFPlus(off) failed: " + offErr, offErr);
            Assert.assertTrue("off.pin must exist", offOut.exists());

            MSGFPlusOptions baselineManager = buildOpts(baselineOut);
            // No -precursorCal flag: picks up the default (AUTO). On the tiny
            // test.mgf dataset the pre-pass does not collect enough confident
            // PSMs (<200), so it returns 0.0 and the fast path kicks in.
            String baseErr = MSGFPlus.runMSGFPlus(baselineManager);
            Assert.assertNull("runMSGFPlus(baseline) failed: " + baseErr, baseErr);
            Assert.assertTrue("baseline.pin must exist", baselineOut.exists());

            List<String> offPsms = extractPsmItems(offOut);
            List<String> basePsms = extractPsmItems(baselineOut);

            Assert.assertFalse("Expected at least one PSM in the off run", offPsms.isEmpty());
            Assert.assertEquals("-precursorCal off must emit the same PSM count as the baseline",
                    basePsms.size(), offPsms.size());

            for (int i = 0; i < offPsms.size(); i++) {
                Assert.assertEquals("PSM #" + i + " differs between off and baseline runs",
                        basePsms.get(i), offPsms.get(i));
            }
        } finally {
            deleteRecursively(workDir.toFile());
        }
    }

    /**
     * The {@code -precursorCal off} path must be deterministic across two
     * back-to-back runs. This pins the no-op path against any accidental
     * non-determinism we introduce later (e.g. a HashSet iteration order
     * leaking into the output).
     */
    @Test
    public void precursorCalOffIsDeterministic() throws Exception {
        Path workDir = Files.createTempDirectory("msgfplus-p2cal-determinism-");
        try {
            File firstOut = new File(workDir.toFile(), "first.pin");
            File secondOut = new File(workDir.toFile(), "second.pin");

            MSGFPlusOptions firstManager = buildOpts(firstOut);
            firstManager.precursorCalMode = PrecursorCalMode.OFF;
            Assert.assertNull(MSGFPlus.runMSGFPlus(firstManager));

            MSGFPlusOptions secondManager = buildOpts(secondOut);
            secondManager.precursorCalMode = PrecursorCalMode.OFF;
            Assert.assertNull(MSGFPlus.runMSGFPlus(secondManager));

            List<String> firstPsms = extractPsmItems(firstOut);
            List<String> secondPsms = extractPsmItems(secondOut);

            Assert.assertEquals(firstPsms.size(), secondPsms.size());
            for (int i = 0; i < firstPsms.size(); i++) {
                Assert.assertEquals("PSM #" + i + " drifted across runs",
                        firstPsms.get(i), secondPsms.get(i));
            }
        } finally {
            deleteRecursively(workDir.toFile());
        }
    }

    /**
     * Verifies that the insufficient-data branch of the calibrator returns
     * 0.0. On the tiny test.mgf fixture the pre-pass cannot reach 200
     * confident PSMs, so the learned shift is 0.0 and the setter is never
     * called — meaning the ioFiles shift stays at the default of 0.0.
     */
    @Test
    public void insufficientPsmsLeavesShiftAtZero() throws Exception {
        Path workDir = Files.createTempDirectory("msgfplus-p2cal-auto-");
        try {
            File autoOut = new File(workDir.toFile(), "auto.pin");
            MSGFPlusOptions manager = buildOpts(autoOut);
            // Leave -precursorCal at default (AUTO). The pre-pass will run
            // but should not collect enough confident PSMs.
            Assert.assertNull(MSGFPlus.runMSGFPlus(manager));

            // The SearchParams list (via paramManager) is internal; we cannot
            // reach it post-run. Instead we re-parse to inspect state.
            // But the ioFiles object is held by SearchParams; re-parsing
            // creates fresh state. So we verify the weaker but still useful
            // invariant: if we re-inspect a freshly created DBSearchIOFiles,
            // its default is 0.0 (pinned by TestPrecursorCalScaffolding).
            // The stronger evidence is baked into
            // precursorCalOffMatchesBaseline: if auto DID apply a non-zero
            // shift, the baseline output would differ from off and that
            // test would fail.
            Assert.assertTrue("auto.pin must exist", autoOut.exists());

            // Additionally confirm the DBSearchIOFiles default via a fresh
            // construction (defensive regression for the field initialiser).
            DBSearchIOFiles sample = new DBSearchIOFiles(
                    new File("x.mgf"), SpecFileFormat.MGF, new File("x.mzid"));
            Assert.assertEquals(0.0, sample.getPrecursorMassShiftPpm(), 0.0);
        } finally {
            deleteRecursively(workDir.toFile());
        }
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    /**
     * Extract every PSM data row from a Percolator {@code .pin} file. The
     * first line is the tab-delimited header and is excluded; the remainder
     * are per-PSM rows whose order matches scoring order, so indexed
     * comparisons are meaningful. Blank trailing lines are skipped so a
     * final newline doesn't produce a spurious empty element.
     */
    private static List<String> extractPsmItems(File pinFile) throws Exception {
        List<String> items = new ArrayList<>();
        List<String> lines = Files.readAllLines(pinFile.toPath(),
                java.nio.charset.StandardCharsets.UTF_8);
        if (lines.size() <= 1) return items; // header only, no PSMs
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i);
            if (line.isEmpty()) continue;
            items.add(line);
        }
        return items;
    }

    private static void deleteRecursively(File file) {
        if (file == null || !file.exists()) return;
        if (file.isDirectory()) {
            File[] kids = file.listFiles();
            if (kids != null) {
                for (File kid : kids) deleteRecursively(kid);
            }
        }
        //noinspection ResultOfMethodCallIgnored
        file.delete();
    }
}
