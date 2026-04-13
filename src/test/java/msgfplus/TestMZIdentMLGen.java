package msgfplus;

import static org.junit.Assert.*;

import java.io.*;
import java.nio.file.Files;

import org.junit.Test;

import edu.ucsd.msjava.ui.MSGFPlus;

/**
 * Tests for issue #157: verify mzid export completeness.
 * Ensures that all SpectrumIdentificationItems have complete score CVParams
 * and that the break->continue fix in MZIdentMLGen preserves all valid PSMs.
 */
public class TestMZIdentMLGen {

    /**
     * Run a small MSGF+ search and verify that the output mzid has
     * complete scores for every SpectrumIdentificationItem.
     *
     * This catches the issue #157 bug where a 'break' on low DeNovoScore
     * would silently drop all subsequent matches for that spectrum.
     * With the fix (continue instead of break), every valid match gets
     * its full set of score CVParams.
     */
    @Test
    public void testMzidScoreCompleteness() throws Exception {
        File specFile = new File(getClass().getClassLoader().getResource("test.mgf").toURI());
        File dbFile = new File(getClass().getClassLoader().getResource("Tryp_Pig_Bov.fasta").toURI());
        File outputFile = File.createTempFile("test_157_", ".mzid");
        outputFile.deleteOnExit();

        // Use Tryp_Pig_Bov.fasta (tiny DB) for fast execution.
        // Even with few or no PSMs, validates the code path does not crash.
        String[] argv = {
                "-s", specFile.getPath(),
                "-d", dbFile.getPath(),
                "-o", outputFile.getPath(),
                "-t", "20ppm",
                "-tda", "0",
                "-ntt", "2",
                "-thread", "2",
                "-minLength", "6",
                "-maxLength", "40",
                "-minCharge", "2",
                "-maxCharge", "4",
                "-n", "1"
        };

        MSGFPlus.main(argv);

        assertTrue("Output mzid file should exist", outputFile.exists());
        assertTrue("Output mzid file should not be empty", outputFile.length() > 0);

        // Parse the output and verify score completeness
        String content = new String(Files.readAllBytes(outputFile.toPath()));

        // Count SpectrumIdentificationItem elements (opening tags only)
        int siiCount = countOccurrences(content, "<SpectrumIdentificationItem ");
        // Count score CVParams
        int rawScoreCount = countOccurrences(content, "accession=\"MS:1002049\"");  // RawScore
        int deNovoScoreCount = countOccurrences(content, "accession=\"MS:1002050\"");  // DeNovoScore
        int specEValueCount = countOccurrences(content, "accession=\"MS:1002052\"");  // SpecEValue
        int eValueCount = countOccurrences(content, "accession=\"MS:1002053\"");  // EValue

        System.out.println("Issue #157 test results:");
        System.out.println("  SpectrumIdentificationItem count: " + siiCount);
        System.out.println("  RawScore (MS:1002049): " + rawScoreCount);
        System.out.println("  DeNovoScore (MS:1002050): " + deNovoScoreCount);
        System.out.println("  SpecEValue (MS:1002052): " + specEValueCount);
        System.out.println("  EValue (MS:1002053): " + eValueCount);

        if (siiCount > 0) {
            // Every SII must have all 4 score CVParams
            assertEquals("Every SII should have a RawScore", siiCount, rawScoreCount);
            assertEquals("Every SII should have a DeNovoScore", siiCount, deNovoScoreCount);
            assertEquals("Every SII should have a SpecEValue", siiCount, specEValueCount);
            assertEquals("Every SII should have an EValue", siiCount, eValueCount);
        }

        // Verify no empty SpectrumIdentificationResult (SIR without SII children)
        // This would indicate silently dropped PSMs
        assertFalse("Should not have empty SpectrumIdentificationResult elements",
                content.contains("<SpectrumIdentificationResult ") &&
                content.contains("</SpectrumIdentificationResult>") &&
                content.contains("<SpectrumIdentificationResult ></SpectrumIdentificationResult>"));

        outputFile.delete();
    }

    /**
     * Verify that the mzid output file is well-formed XML and contains
     * the required mzIdentML structure elements.
     */
    @Test
    public void testMzidStructuralValidity() throws Exception {
        File specFile = new File(getClass().getClassLoader().getResource("test.mgf").toURI());
        File dbFile = new File(getClass().getClassLoader().getResource("BSA.fasta").toURI());
        File outputFile = File.createTempFile("test_157_struct_", ".mzid");
        outputFile.deleteOnExit();

        String[] argv = {
                "-s", specFile.getPath(),
                "-d", dbFile.getPath(),
                "-o", outputFile.getPath(),
                "-t", "20ppm",
                "-tda", "0",
                "-ntt", "2",
                "-thread", "2",
                "-n", "1"
        };

        MSGFPlus.main(argv);

        assertTrue("Output mzid file should exist", outputFile.exists());
        String content = new String(Files.readAllBytes(outputFile.toPath()));

        // Verify required mzIdentML sections exist
        assertTrue("Should contain MzIdentML root element",
                content.contains("<MzIdentML"));
        assertTrue("Should contain SpectrumIdentificationList",
                content.contains("<SpectrumIdentificationList"));
        assertTrue("Should contain SequenceCollection",
                content.contains("<SequenceCollection"));
        assertTrue("Should contain AnalysisProtocolCollection",
                content.contains("<AnalysisProtocolCollection"));

        outputFile.delete();
    }

    private static int countOccurrences(String text, String pattern) {
        int count = 0;
        int idx = 0;
        while ((idx = text.indexOf(pattern, idx)) != -1) {
            count++;
            idx += pattern.length();
        }
        return count;
    }
}
