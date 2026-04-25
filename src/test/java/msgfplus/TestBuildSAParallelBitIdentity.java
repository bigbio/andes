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
 * Bit-identity test: the parallel sort path must produce byte-identical
 * .csarr/.cnlcp output to the single-thread path between the header and footer
 * (header id and footer mtime are non-deterministic between builds).
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

    /** Copies the FASTA fixture into a fresh temp dir so build artifacts don't pollute test resources. */
    private static File stageFastaIntoTempDir(String prefix) throws Exception {
        Path tempDir = Files.createTempDirectory(prefix);
        File source = new File(TestBuildSAParallelBitIdentity.class.getClassLoader().getResource(FIXTURE).toURI());
        File dest = new File(tempDir.toFile(), source.getName());
        Files.copy(source.toPath(), dest.toPath(), StandardCopyOption.REPLACE_EXISTING);
        return dest;
    }

    /**
     * Read the file with the 8-byte header (size + id) and 12-byte footer
     * (lastModified + formatId) trimmed off. Both are non-deterministic between
     * runs; the body in between is the actual sort output to compare.
     */
    private static byte[] readBodyBytes(File f) throws IOException {
        byte[] all = Files.readAllBytes(f.toPath());
        int headerSize = 8;
        int footerSize = 8 + 4;
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
