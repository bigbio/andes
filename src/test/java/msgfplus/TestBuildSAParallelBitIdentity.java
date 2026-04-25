package msgfplus;

import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

/**
 * BuildSA bit-identity test: the parallel temp-file path must produce
 * byte-identical .csarr/.cnlcp output to the single-thread direct-write path
 * past the 8-byte header (size + id; the id changes per-build via
 * non-deterministic UUID-hash in CompactFastaSequence).
 *
 * <p>Runs the same FASTA through both paths, captures the resulting bytes,
 * and verifies they match. Also verifies that the parallel path leaves no
 * temp files behind on success.
 */
public class TestBuildSAParallelBitIdentity {

    /** Mirror of {@code CompactSuffixArray.SA_BUILD_THREADS_PROPERTY} (package-private there). */
    private static final String SA_BUILD_THREADS_PROPERTY = "msgfplus.buildsa.threads";

    private static final String FIXTURE = "ecoli.fasta";

    @Test
    public void parallelMatchesSingleThreadByteForByte() throws Exception {
        File singleArtifacts = stageFastaIntoTempDir("buildsa-N1");
        File parallelArtifacts = stageFastaIntoTempDir("buildsa-N4");
        try {
            byte[] singleCsarr, singleCnlcp;
            byte[] parallelCsarr, parallelCnlcp;

            String prevThreads = System.getProperty(SA_BUILD_THREADS_PROPERTY);
            try {
                System.setProperty(SA_BUILD_THREADS_PROPERTY, "1");
                CompactFastaSequence seq1 = new CompactFastaSequence(singleArtifacts.getAbsolutePath());
                new CompactSuffixArray(seq1);
                singleCsarr = readBodyBytes(new File(stripExt(singleArtifacts.getAbsolutePath()) + ".csarr"));
                singleCnlcp = readBodyBytes(new File(stripExt(singleArtifacts.getAbsolutePath()) + ".cnlcp"));

                System.setProperty(SA_BUILD_THREADS_PROPERTY, "4");
                CompactFastaSequence seq4 = new CompactFastaSequence(parallelArtifacts.getAbsolutePath());
                new CompactSuffixArray(seq4);
                parallelCsarr = readBodyBytes(new File(stripExt(parallelArtifacts.getAbsolutePath()) + ".csarr"));
                parallelCnlcp = readBodyBytes(new File(stripExt(parallelArtifacts.getAbsolutePath()) + ".cnlcp"));
            } finally {
                if (prevThreads == null) {
                    System.clearProperty(SA_BUILD_THREADS_PROPERTY);
                } else {
                    System.setProperty(SA_BUILD_THREADS_PROPERTY, prevThreads);
                }
            }

            Assert.assertArrayEquals(".csarr post-header bytes must be identical between N=1 and N=4", singleCsarr, parallelCsarr);
            Assert.assertArrayEquals(".cnlcp post-header bytes must be identical between N=1 and N=4", singleCnlcp, parallelCnlcp);

            // No temp debris left behind in the parallel build's directory.
            File parentDir = parallelArtifacts.getAbsoluteFile().getParentFile();
            File[] debris = parentDir.listFiles((dir, name) -> name.contains(".buildsa-tmp."));
            Assert.assertNotNull(debris);
            Assert.assertEquals("BuildSA temp files must be cleaned up on success: " + java.util.Arrays.toString(debris),
                    0, debris.length);
        } finally {
            deleteDirRecursive(singleArtifacts.getParentFile());
            deleteDirRecursive(parallelArtifacts.getParentFile());
        }
    }

    /**
     * Copies the {@link #FIXTURE} into a fresh temp directory so we can build
     * {@code .canno / .cseq / .csarr / .cnlcp} alongside it without polluting
     * {@code src/test/resources}.
     */
    private static File stageFastaIntoTempDir(String prefix) throws Exception {
        Path tempDir = Files.createTempDirectory(prefix);
        File source = new File(TestBuildSAParallelBitIdentity.class.getClassLoader().getResource(FIXTURE).toURI());
        File dest = new File(tempDir.toFile(), source.getName());
        Files.copy(source.toPath(), dest.toPath(), StandardCopyOption.REPLACE_EXISTING);
        return dest;
    }

    /**
     * Read the file and skip both the 8-byte header (size int + id int — the id
     * is non-deterministic UUID-hash) and the 12-byte footer (lastModified long
     * + formatId int — the lastModified differs because the test stages the
     * fixture into a fresh temp dir per run, so the FASTA's mtime differs
     * between the two runs). The bytes between are the actual sort output we
     * want to compare bit-for-bit.
     */
    private static byte[] readBodyBytes(File f) throws IOException {
        byte[] all = Files.readAllBytes(f.toPath());
        int headerSize = 8;       // int size + int id
        int footerSize = 8 + 4;   // long lastModified + int formatId
        Assert.assertTrue("Output file too small: " + f, all.length >= headerSize + footerSize);
        int bodyLen = all.length - headerSize - footerSize;
        byte[] body = new byte[bodyLen];
        System.arraycopy(all, headerSize, body, 0, bodyLen);
        return body;
    }

    private static String stripExt(String path) {
        int dot = path.lastIndexOf('.');
        return dot < 0 ? path : path.substring(0, dot);
    }

    private static void deleteDirRecursive(File dir) {
        if (dir == null || !dir.exists()) return;
        File[] entries = dir.listFiles();
        if (entries != null) {
            for (File f : entries) {
                if (f.isDirectory()) deleteDirRecursive(f);
                else f.delete();
            }
        }
        dir.delete();
    }
}
