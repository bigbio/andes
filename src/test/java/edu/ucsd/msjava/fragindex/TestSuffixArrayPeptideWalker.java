package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msutil.Enzyme;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class TestSuffixArrayPeptideWalker {

    private Path tempDir;

    @Before
    public void setup() throws IOException {
        tempDir = Files.createTempDirectory("sapw-test-");
    }

    @After
    public void teardown() throws IOException {
        if (tempDir != null) {
            // Best-effort cleanup of generated .cseq / .canno / .fasta files.
            File[] files = tempDir.toFile().listFiles();
            if (files != null) {
                for (File f : files) {
                    // noinspection ResultOfMethodCallIgnored
                    f.delete();
                }
            }
            // noinspection ResultOfMethodCallIgnored
            tempDir.toFile().delete();
        }
    }

    private CompactFastaSequence writeFasta(String name, String body) throws IOException {
        File f = tempDir.resolve(name + ".fasta").toFile();
        try (PrintWriter pw = new PrintWriter(f)) {
            pw.print(body);
        }
        return new CompactFastaSequence(f.getAbsolutePath());
    }

    /** FASTA with two small proteins used by several tests. */
    private static final String TWO_PROTEIN_FASTA =
            ">sp|P00001|PROT1 Test protein 1\n" +
            "MAEKVLR\n" +
            ">sp|P00002|PROT2 Test protein 2\n" +
            "KPIR\n";

    @Test
    public void basicTrypticWalkYieldsExpectedPeptides() throws IOException {
        CompactFastaSequence seq = writeFasta("basic", TWO_PROTEIN_FASTA);
        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 3, 10, 0);

        List<String> peptides = walker.collect();
        // Expected strict-tryptic peptides (length >= 3, no missed cleavages):
        //   PROT1 (MAEKVLR): MAEK, VLR   (K-length-1 filtered; no missed cleavages allowed)
        //   PROT2 (KPIR):    PIR          (leading K is length-1)
        Assert.assertEquals(Arrays.asList("MAEK", "VLR", "PIR"), peptides);
    }

    @Test
    public void missedCleavagesEmitSpanningPeptides() throws IOException {
        CompactFastaSequence seq = writeFasta("mmc", TWO_PROTEIN_FASTA);
        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 3, 20, 1);

        List<String> peptides = walker.collect();
        // With 1 missed cleavage:
        //   PROT1: MAEK, MAEKVLR, VLR
        //   PROT2: KPIR, PIR           (K alone is length 1, skipped)
        Assert.assertEquals(Arrays.asList("MAEK", "MAEKVLR", "VLR", "KPIR", "PIR"), peptides);
    }

    @Test
    public void lengthFilterExcludesShortPeptides() throws IOException {
        CompactFastaSequence seq = writeFasta("lenfilter", TWO_PROTEIN_FASTA);
        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 5, 20, 1);

        List<String> peptides = walker.collect();
        // minLen=5 filters out MAEK(4), VLR(3), KPIR(4), PIR(3). Only MAEKVLR(7) survives.
        Assert.assertEquals(Collections.singletonList("MAEKVLR"), peptides);
    }

    @Test
    public void decoyProteinsAreSkipped() throws IOException {
        String fasta =
                ">sp|P00001|PROT1 target\n" +
                "MAEKVLR\n" +
                ">XXX_sp|P00001|PROT1 decoy\n" +
                "RLVKEAM\n";
        CompactFastaSequence seq = writeFasta("decoy", fasta);
        seq.setDecoyProteinPrefix("XXX_");

        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 3, 10, 0);
        List<String> peptides = walker.collect();
        // Only PROT1 peptides: MAEK, VLR. Decoy RLVK / EAM should be absent.
        Assert.assertEquals(Arrays.asList("MAEK", "VLR"), peptides);
    }

    @Test
    public void residueOutsideAlphabetRejectsPeptide() throws IOException {
        // '*' is the stop-codon symbol; it is NOT in the default 26-letter alphabet
        // and gets stored as INVALID_CHAR_CODE in the compact buffer.
        String fasta =
                ">sp|PROT|contains-stop\n" +
                "MAEK*LR\n";
        CompactFastaSequence seq = writeFasta("nonalpha", fasta);
        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 2, 10, 1);

        List<String> peptides = walker.collect();
        // Buffer: M A E K ? L R  (the '*' becomes '?' on decode and is not in the alphabet)
        // Cleavage sites after K (pos 5) only; '?' is not cleavable.
        // Candidate peptides: MAEK(1..5) -> clean, emit. ?LR(5..end) contains '?'
        //   -> rejected. MAEK?LR would also contain '?' -> rejected.
        Assert.assertTrue("MAEK should be emitted",   peptides.contains("MAEK"));
        for (String p : peptides) {
            Assert.assertFalse("peptide must not contain invalid char: " + p, p.indexOf('?') >= 0);
        }
    }

    @Test
    public void proteinShorterThanMinLenProducesNoOutput() throws IOException {
        String fasta =
                ">sp|TINY|short\n" +
                "MK\n";
        CompactFastaSequence seq = writeFasta("tiny", fasta);
        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 5, 20, 0);
        Assert.assertTrue(walker.collect().isEmpty());
    }

    @Test
    public void collectMatchesForEachPeptide() throws IOException {
        CompactFastaSequence seq = writeFasta("consistent", TWO_PROTEIN_FASTA);
        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 2, 20, 2);

        List<String> viaCollect = walker.collect();

        List<String> viaForEach = new ArrayList<>();
        walker.forEachPeptide(viaForEach::add);

        Assert.assertEquals(viaForEach, viaCollect);
    }

    @Test
    public void constructorValidatesArguments() throws IOException {
        CompactFastaSequence seq = writeFasta("args", TWO_PROTEIN_FASTA);
        try {
            new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 0, 10, 0);
            Assert.fail("expected IllegalArgumentException for minLen < 1");
        } catch (IllegalArgumentException expected) {
            // ok
        }
        try {
            new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 5, 3, 0);
            Assert.fail("expected IllegalArgumentException for maxLen < minLen");
        } catch (IllegalArgumentException expected) {
            // ok
        }
        try {
            new SuffixArrayPeptideWalker(seq, Enzyme.TRYPSIN, 3, 10, -1);
            Assert.fail("expected IllegalArgumentException for maxMissedCleavages < 0");
        } catch (IllegalArgumentException expected) {
            // ok
        }
    }

    @Test
    public void lysCEnzymeCleavesOnlyAfterK() throws IOException {
        // MAEKVLRKPIR with LysC: cleaves after each K only.
        // Sites in PROT1 (MAEKVLR): [start, after-K=after pos4, end].
        //   -> MAEK, VLR (no cleavage after R for LysC).
        // Sites in PROT2 (KPIR): [start, after-K, end]
        //   -> K (length 1, filtered), PIR.
        CompactFastaSequence seq = writeFasta("lysc", TWO_PROTEIN_FASTA);
        SuffixArrayPeptideWalker walker =
                new SuffixArrayPeptideWalker(seq, Enzyme.LysC, 3, 20, 0);
        List<String> peptides = walker.collect();
        Assert.assertEquals(Arrays.asList("MAEK", "VLR", "PIR"), peptides);
    }
}
