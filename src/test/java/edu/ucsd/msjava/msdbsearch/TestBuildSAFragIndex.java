package edu.ucsd.msjava.msdbsearch;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.Assert;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Covers the {@code -buildFragIndex 0|1} CLI flag on {@link BuildSA}.
 *
 * <p>Uses a tiny synthetic FASTA written into a temp directory so each test
 * finishes in a few seconds — the bundled {@code human-uniprot-contaminants.fasta}
 * fixture is the full human proteome (82 MB) and triggers OOM in the
 * in-memory fragment-index builder. Phase 2b is only validating the CLI
 * plumbing and the logged "Fragment index built" line; a small hand-written
 * fixture is sufficient for that contract.
 *
 * <p>BuildSA's own stdout is redirected into a buffer so the test can
 * assert on the presence (or absence) of the "Fragment index built" log line.
 */
public class TestBuildSAFragIndex {
    private PrintStream originalOut;
    private ByteArrayOutputStream captured;
    private Path workDir;

    @Before
    public void setUp() throws Exception {
        originalOut = System.out;
        captured = new ByteArrayOutputStream();
        System.setOut(new PrintStream(captured));
        workDir = Files.createTempDirectory("msgfplus-buildsa-fragindex-");
    }

    @After
    public void tearDown() {
        System.setOut(originalOut);
        // Temp dir cleanup is best-effort; the OS cleans /tmp regularly.
    }

    /** Writes a small multi-protein FASTA that yields a handful of tryptic peptides. */
    private File writeTinyFasta() throws Exception {
        File fasta = new File(workDir.toFile(), "tiny.fasta");
        try (PrintStream ps = new PrintStream(fasta)) {
            ps.println(">sp|PROT1|PROT1_TEST test protein 1");
            ps.println("MAEKVLRKPIRAEDLGR");
            ps.println(">sp|PROT2|PROT2_TEST test protein 2");
            ps.println("SAMPLERSEQUENCEKPEPTIDER");
            ps.println(">sp|PROT3|PROT3_TEST test protein 3");
            ps.println("ACDEFGHIKLMNPQRSTVWYR");
        }
        return fasta;
    }

    @Test
    public void buildFragIndexFlagTriggersFragmentIndexBuild() throws Exception {
        File fasta = writeTinyFasta();

        BuildSA.main(new String[]{
                "-d", fasta.getAbsolutePath(),
                "-tda", "0",
                "-buildFragIndex", "1"
        });

        String out = captured.toString();
        Assert.assertTrue("expected 'Fragment index built' log line, got:\n" + out,
                out.contains("Fragment index built"));
        // Peptide count > 0 — this tiny fasta must yield at least some tryptic peptides.
        Assert.assertTrue("log should include a positive peptide count, got:\n" + out,
                out.matches("(?s).*Fragment index built for [^:]+: [1-9][0-9]* peptides.*"));
    }

    @Test
    public void defaultNoFlagDoesNotBuildFragmentIndex() throws Exception {
        File fasta = writeTinyFasta();

        BuildSA.main(new String[]{
                "-d", fasta.getAbsolutePath(),
                "-tda", "0"
        });

        String out = captured.toString();
        Assert.assertFalse("default path must NOT emit the 'Fragment index built' line, got:\n" + out,
                out.contains("Fragment index built"));
    }

    @Test
    public void buildFragIndexZeroIsSameAsDefault() throws Exception {
        File fasta = writeTinyFasta();

        BuildSA.main(new String[]{
                "-d", fasta.getAbsolutePath(),
                "-tda", "0",
                "-buildFragIndex", "0"
        });

        String out = captured.toString();
        Assert.assertFalse("-buildFragIndex 0 must NOT emit the 'Fragment index built' line, got:\n" + out,
                out.contains("Fragment index built"));
    }
}
