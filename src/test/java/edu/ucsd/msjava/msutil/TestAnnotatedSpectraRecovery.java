package edu.ucsd.msjava.msutil;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Integration tests for the restored {@link AnnotatedSpectra} TSV/spectrum
 * loader. Covers:
 *
 *  - Extension validation: .mzid is rejected with the explicit "not supported
 *    in this build" message (the mzid/ package was removed); .tsv is the
 *    only accepted format.
 *  - TSV header validation: required columns must be present; missing columns
 *    produce specific errors.
 *  - End-to-end TSV + spectrum-file integration using the existing
 *    src/test/resources/iprg-2013/F13.mgf fixture (real MS2 data, 146K lines).
 *    Exercises SpectraAccessor lookup against an MGF spec dir without
 *    requiring any specific instrument data; the recovered code is
 *    instrument-agnostic. Astral-format mzML files are too large to ship as
 *    test fixtures, but the integration boundary is identical: parse a TSV
 *    row, look up the named SpecID via SpectraAccessor, return a graceful
 *    error if not found.
 *
 * Does NOT exercise the actual scoring-parameter trainer
 * (ScoringParameterGeneratorWithErrors.generateParameters); that needs
 * thousands of confidently-annotated PSMs and is out of scope for a recovery
 * smoke test.
 */
public class TestAnnotatedSpectraRecovery {

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    private static final String F13_MGF_RESOURCE_PATH = "src/test/resources/iprg-2013/F13.mgf";

    private File writeTsv(String content) throws IOException {
        File f = tmp.newFile("results.tsv");
        Files.write(f.toPath(), content.getBytes(StandardCharsets.UTF_8));
        return f;
    }

    private AnnotatedSpectra newWith(File... resultFiles) {
        return new AnnotatedSpectra(resultFiles, tmp.getRoot(), AminoAcidSet.getStandardAminoAcidSet());
    }

    @Test
    public void rejectsMzidExtensionWithExplicitMessage() throws IOException {
        File mzid = tmp.newFile("results.mzid");
        Files.write(mzid.toPath(), "<placeholder/>".getBytes(StandardCharsets.UTF_8));

        AnnotatedSpectra parser = newWith(mzid);
        String err = parser.parseFile(mzid);

        assertNotNull("mzid input should be rejected", err);
        assertTrue("error should call out that mzid is unsupported in this build; got: " + err,
                err.contains("mzid input is not supported in this build"));
    }

    @Test
    public void rejectsUnknownExtension() throws IOException {
        File txt = tmp.newFile("results.txt");
        Files.write(txt.toPath(), "irrelevant".getBytes(StandardCharsets.UTF_8));

        AnnotatedSpectra parser = newWith(txt);
        String err = parser.parseFile(txt);

        assertNotNull(err);
        assertTrue("error should mention the unrecognized extension; got: " + err,
                err.startsWith("Unrecognized result file extension"));
    }

    @Test
    public void rejectsTsvWithoutHeaderHash() throws IOException {
        // Header line must start with "#" per MS-GF+ result format. A blank
        // first line (or a line that doesn't start with #) is rejected.
        File tsv = writeTsv("SpecFile\tSpecID\tPeptide\tCharge\tQValue\nfoo\t1\tK.PEPTIDE.K\t2\t0.001\n");

        AnnotatedSpectra parser = newWith(tsv);
        String err = parser.parseFile(tsv);

        assertEquals("Not a valid tsv result file", err);
    }

    @Test
    public void rejectsTsvMissingRequiredColumns() throws IOException {
        // Header has only #SpecFile and Charge; SpecID and Peptide are required.
        File tsv = writeTsv("#SpecFile\tCharge\nfoo\t2\n");

        AnnotatedSpectra parser = newWith(tsv);
        String err = parser.parseFile(tsv);

        // Note: the upstream error message reads "Not a valid mzid file" even
        // though the input is .tsv — that string predates the mzID excision.
        // Keeping it as-is to stay faithful to the upstream parser; the test
        // documents the actual behavior so future readers aren't surprised.
        assertNotNull(err);
        assertTrue("missing required columns should be flagged; got: " + err,
                err.contains("Not a valid"));
    }

    @Test
    public void rejectsTsvMissingFdrColumn() throws IOException {
        // Header has SpecID/SpecFile/Peptide but no QValue/EFDR/SpecQValue/FDR.
        File tsv = writeTsv("#SpecFile\tSpecID\tPeptide\tCharge\nfoo\t1\tK.PEPTIDE.K\t2\n");

        AnnotatedSpectra parser = newWith(tsv);
        String err = parser.parseFile(tsv);

        assertEquals("QValue is missing", err);
    }

    @Test
    public void acceptsAnyOfTheFdrColumnAliases() throws IOException {
        // The parser recognises FDR / EFDR / QValue / SpecQValue interchangeably
        // as the FDR column. Verify each alias produces a passing header check
        // (no error). Use parse() (not parseFile() directly) because parseFile
        // expects the container to have been initialised by parse() first.
        for (String fdrAlias : new String[]{"FDR", "EFDR", "QValue", "SpecQValue"}) {
            File tsv = tmp.newFile("results-" + fdrAlias + ".tsv");
            Files.write(tsv.toPath(),
                    ("#SpecFile\tSpecID\tPeptide\tCharge\t" + fdrAlias + "\n")
                            .getBytes(StandardCharsets.UTF_8));
            AnnotatedSpectra parser = newWith(tsv);
            String err = parser.parse();
            // No data rows -> no lookups -> parser returns null (success).
            assertNull("FDR column alias '" + fdrAlias + "' should be accepted; got: " + err, err);
        }
    }

    @Test
    public void filtersRowsAboveFdrThreshold() throws IOException {
        // Default fdrThreshold is 0.01. A row with QValue=0.99 is dropped;
        // no spectrum lookup is attempted, so the parser succeeds without
        // touching the spec dir.
        File tsv = writeTsv(
                "#SpecFile\tSpecID\tPeptide\tCharge\tQValue\n"
                        + "irrelevant.mgf\tindex=0\tK.PEPTIDE.K\t2\t0.99\n");

        AnnotatedSpectra parser = newWith(tsv);
        String err = parser.parse();

        assertNull("rows above FDR threshold should be dropped silently; got: " + err, err);
        assertNotNull(parser.getAnnotatedSpecContainer());
        assertTrue("no PSMs should have passed the filter",
                parser.getAnnotatedSpecContainer().isEmpty());
    }

    @Test
    public void looksUpSpectraInRealMgfFixtureAndReportsMissingId() throws IOException {
        // Drive the full TSV-row parse path against the real F13.mgf fixture.
        // The TSV references a SpecID that does not exist in F13.mgf, so the
        // SpectraAccessor returns null and the parser produces a clean
        // "is not available!" error rather than crashing.
        Path mgfSrc = new File(F13_MGF_RESOURCE_PATH).toPath();
        assertTrue("test fixture missing: " + F13_MGF_RESOURCE_PATH, Files.exists(mgfSrc));

        File specDir = tmp.newFolder("specs");
        Path mgfDest = new File(specDir, "F13.mgf").toPath();
        Files.copy(mgfSrc, mgfDest);

        File tsv = writeTsv(
                "#SpecFile\tSpecID\tPeptide\tCharge\tQValue\n"
                        + "F13.mgf\tindex=99999999\tK.PEPTIDE.K\t2\t0.001\n");

        AnnotatedSpectra parser = new AnnotatedSpectra(
                new File[]{tsv}, specDir, AminoAcidSet.getStandardAminoAcidSet());
        String err = parser.parseFile(tsv);

        assertNotNull("missing spec ID should yield an error", err);
        assertTrue("error should name the file and id; got: " + err,
                err.contains("F13.mgf:index=99999999") && err.contains("is not available"));
    }

    @Test
    public void parseTopLevelHandlesValidEmptyTsvCleanly() throws IOException {
        // The high-level parse() method aggregates per-file results; an empty
        // (header-only) TSV produces no PSMs and no error.
        File tsv = writeTsv("#SpecFile\tSpecID\tPeptide\tCharge\tQValue\n");

        AnnotatedSpectra parser = newWith(tsv);
        String err = parser.parse();

        assertNull("empty but well-formed TSV should not error; got: " + err, err);
        assertNotNull(parser.getAnnotatedSpecContainer());
        assertTrue(parser.getAnnotatedSpecContainer().isEmpty());
    }

    @Test
    public void parseAggregatesMultipleFiles() throws IOException {
        // Two empty-but-valid TSVs — both succeed, container is empty,
        // no exceptions propagate. Proves the multi-file iteration loop works.
        File tsv1 = writeTsv("#SpecFile\tSpecID\tPeptide\tCharge\tQValue\n");
        File tsv2 = tmp.newFile("results2.tsv");
        Files.write(tsv2.toPath(),
                "#SpecFile\tSpecID\tPeptide\tCharge\tQValue\n".getBytes(StandardCharsets.UTF_8));

        AnnotatedSpectra parser = newWith(tsv1, tsv2);
        String err = parser.parse();

        assertNull("aggregating empty TSVs should not error; got: " + err, err);
        assertTrue(parser.getAnnotatedSpecContainer().isEmpty());
    }

    @Test
    public void parseFailFastOnFirstErrorWhenDropErrorsDisabled() throws IOException {
        // Two TSVs, first is invalid (.mzid) → parse() returns the first
        // error and stops.
        File mzid = tmp.newFile("first.mzid");
        Files.write(mzid.toPath(), "<placeholder/>".getBytes(StandardCharsets.UTF_8));
        File tsv = writeTsv("#SpecFile\tSpecID\tPeptide\tCharge\tQValue\n");

        AnnotatedSpectra parser = new AnnotatedSpectra(
                new File[]{mzid, tsv}, tmp.getRoot(), AminoAcidSet.getStandardAminoAcidSet());
        String err = parser.parse();

        assertNotNull(err);
        assertTrue("should surface the mzid-rejection message; got: " + err,
                err.contains("mzid input is not supported in this build"));
    }
}
