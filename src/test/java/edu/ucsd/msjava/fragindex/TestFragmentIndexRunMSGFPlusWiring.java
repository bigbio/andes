package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * End-to-end wiring tests for Phase 3 commit 3.2 — the per-file
 * {@link FragmentIndex} build that sits between {@code SpecKey} list
 * construction and worker-task fan-out in
 * {@link MSGFPlus#runMSGFPlus(ParamManager)}.
 *
 * <p>This commit does not yet consume the index in the search loop, so
 * the tests pin the control-flow contract rather than any scoring
 * behaviour:
 *
 * <ul>
 *   <li>{@link #offModeDoesNotBuildFragmentIndex()} — the primary
 *       correctness gate. When {@code -useFragmentIndex off} is supplied,
 *       the build block must be skipped entirely; the
 *       {@code "Building fragment index"} log line must NOT appear. A
 *       regression here would mean the legacy path is paying an
 *       unnecessary fragment-index build cost, and we would no longer be
 *       bit-identical to HEAD on the off path.</li>
 * </ul>
 *
 * <p>The on-mode build is validated by end-to-end benchmarks and the
 * subsequent DBScanner wiring commit (3.4); we intentionally do NOT run
 * a full search here because the bundled
 * {@code human-uniprot-contaminants.fasta} (~82 MB) produced an OOM in
 * the earlier Phase 2b full-FASTA test, and full searches also take
 * several minutes per pass which is inappropriate for unit-test
 * latency.
 */
public class TestFragmentIndexRunMSGFPlusWiring {

    private ParamManager buildParamManager(File outputFile) throws URISyntaxException {
        ParamManager manager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE,
                "java -Xmx3500M -jar MSGFPlus.jar");
        manager.addMSGFPlusParams();

        URI paramUri = SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI();
        manager.getParameter("conf").parse(new File(paramUri).getAbsolutePath());

        URI specUri = SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI();
        manager.getParameter("s").parse(new File(specUri).getAbsolutePath());

        URI dbUri = SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI();
        manager.getParameter("d").parse(new File(dbUri).getAbsolutePath());

        manager.getParameter("o").parse(outputFile.getAbsolutePath());
        return manager;
    }

    /**
     * Primary correctness gate: with {@code -useFragmentIndex off} the
     * build block must never execute. We capture stdout, run a full
     * {@code MSGFPlus.runMSGFPlus} pass, and assert the telltale
     * {@code "Building fragment index"} log line is absent.
     */
    @Test
    public void offModeDoesNotBuildFragmentIndex() throws Exception {
        Path workDir = Files.createTempDirectory("msgfplus-fragindex-off-");
        PrintStream originalOut = System.out;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            File outFile = new File(workDir.toFile(), "off.pin");
            ParamManager manager = buildParamManager(outFile);
            Assert.assertNull(manager.getParameter("useFragmentIndex").parse("off"));

            System.setOut(new PrintStream(captured, true, StandardCharsets.UTF_8));
            try {
                String err = MSGFPlus.runMSGFPlus(manager);
                Assert.assertNull("runMSGFPlus(off) failed: " + err, err);
            } finally {
                System.setOut(originalOut);
            }

            String stdout = captured.toString(StandardCharsets.UTF_8);
            Assert.assertFalse(
                    "OFF mode must NOT emit the 'Building fragment index' log line; saw:\n" + stdout,
                    stdout.contains("Building fragment index"));
            Assert.assertFalse(
                    "OFF mode must NOT emit the 'Fragment index built' log line; saw:\n" + stdout,
                    stdout.contains("Fragment index built"));
            Assert.assertTrue("off.pin must exist", outFile.exists());
        } finally {
            System.setOut(originalOut);
            deleteRecursively(workDir.toFile());
        }
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
