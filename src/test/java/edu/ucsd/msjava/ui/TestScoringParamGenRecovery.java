package edu.ucsd.msjava.ui;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Integration tests for the restored {@link ScoringParamGen} CLI front door.
 *
 * Scope: the recovery surface only — argument parsing, validation of required
 * options, error paths for invalid enums (activation method, instrument,
 * enzyme, protocol), and the TSV/spectrum-dir existence checks. The actual
 * .param model generation is not exercised; that requires real annotated
 * training PSMs which are not part of this fixture set.
 */
public class TestScoringParamGenRecovery {

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    private File makeEmptyTsv() throws IOException {
        File f = tmp.newFile("results.tsv");
        Files.write(f.toPath(), new byte[0]);
        return f;
    }

    private File makeSpecDir() throws IOException {
        return tmp.newFolder("spectra");
    }

    @Test
    public void runWithoutArgsReportsMissingI() {
        assertEquals("missing -i (training result TSV files)",
                ScoringParamGen.run(new String[]{}));
    }

    @Test
    public void runWithSingleDanglingArgReportsInvalidParameter() {
        // The arg parser requires key/value pairs; a lone "-i" is malformed.
        assertEquals("Invalid parameter: -i",
                ScoringParamGen.run(new String[]{"-i"}));
    }

    @Test
    public void runWithBareTokenReportsInvalidParameter() {
        // First arg must start with "-". "junk" is a bare positional, rejected.
        assertEquals("Invalid parameter: junk",
                ScoringParamGen.run(new String[]{"junk", "value"}));
    }

    @Test
    public void runWithUnknownOptionRejected() throws IOException {
        File tsv = makeEmptyTsv();
        String err = ScoringParamGen.run(new String[]{"-i", tsv.getPath(), "-bogus", "x"});
        assertEquals("Unknown option: -bogus", err);
    }

    @Test
    public void runWithMissingInputFileRejected() {
        String err = ScoringParamGen.run(new String[]{"-i", "/nonexistent/path/to/results.tsv"});
        assertNotNull(err);
        assertTrue("error should mention the bad path; got: " + err,
                err.startsWith("Input file does not exist:"));
    }

    @Test
    public void runWithMissingSpecDirRejected() throws IOException {
        File tsv = makeEmptyTsv();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", "/nonexistent/directory"});
        assertNotNull(err);
        assertTrue("error should mention the missing dir; got: " + err,
                err.startsWith("Spectrum directory does not exist"));
    }

    @Test
    public void runWithInvalidActivationMethodRejected() throws IOException {
        File tsv = makeEmptyTsv();
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath(),
                "-m", "FAKE_ACTIVATION"});
        assertEquals("Unrecognized activation method: FAKE_ACTIVATION", err);
    }

    @Test
    public void runWithInvalidInstrumentRejected() throws IOException {
        File tsv = makeEmptyTsv();
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath(),
                "-m", "HCD",
                "-inst", "FAKE_INSTRUMENT"});
        assertEquals("Unrecognized instrument type: FAKE_INSTRUMENT", err);
    }

    @Test
    public void runWithInvalidEnzymeRejected() throws IOException {
        File tsv = makeEmptyTsv();
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath(),
                "-m", "HCD",
                "-inst", "QExactive",
                "-e", "FAKE_ENZYME"});
        assertEquals("Unrecognized enzyme: FAKE_ENZYME", err);
    }

    @Test
    public void runWithInvalidProtocolRejected() throws IOException {
        File tsv = makeEmptyTsv();
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath(),
                "-m", "HCD",
                "-inst", "QExactive",
                "-e", "Tryp",
                "-protocol", "FAKE_PROTOCOL"});
        assertEquals("Unrecognized protocol: FAKE_PROTOCOL", err);
    }

    @Test
    public void runWithNonNumericThreadRejected() throws IOException {
        File tsv = makeEmptyTsv();
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath(),
                "-m", "HCD",
                "-inst", "QExactive",
                "-e", "Tryp",
                "-thread", "abc"});
        assertEquals("-thread must be an integer", err);
    }

    @Test
    public void runWithMissingDReportedAfterIIsValid() throws IOException {
        File tsv = makeEmptyTsv();
        String err = ScoringParamGen.run(new String[]{"-i", tsv.getPath()});
        assertEquals("missing -d (spectrum directory)", err);
    }

    @Test
    public void runWithMissingMReportedAfterIDIsValid() throws IOException {
        File tsv = makeEmptyTsv();
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath()});
        assertEquals("missing -m (activation method)", err);
    }

    @Test
    public void runWithAllRequiredArgsButEmptyTsvFailsAtParse() throws IOException {
        // All CLI args validate; the empty TSV (no header line) fails inside
        // AnnotatedSpectra.parseFile with "Not a valid tsv result file".
        // This proves the CLI hands off to the parser correctly post-recovery.
        File tsv = makeEmptyTsv();
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath(),
                "-m", "HCD",
                "-inst", "QExactive",
                "-e", "Tryp"});
        assertNotNull("expected a parse error to bubble up; got null", err);
        assertTrue("error should come from TSV parsing; got: " + err,
                err.contains("Not a valid tsv result file")
                        || err.startsWith("Error while parsing"));
    }

    @Test
    public void runWithNoMatchingPsmsExitsCleanlyBeforeTraining() throws IOException {
        // Build a valid-shaped TSV whose only row has FDR > 0.01 (the trainer's
        // FDR floor). AnnotatedSpectra accepts the file, filters the row out,
        // and the CLI returns "No results to train on" *without* invoking the
        // ScoringParameterGeneratorWithErrors trainer. Proves the full
        // arg-parse + handoff + early-exit chain works end to end.
        File tsv = tmp.newFile("noresults.tsv");
        String content = String.join("\t", "#SpecFile", "SpecID", "Peptide", "Charge", "QValue") + "\n"
                + String.join("\t", "synth.mgf", "index=0", "K.PEPTIDE.K", "2", "0.99") + "\n";
        Files.write(tsv.toPath(), content.getBytes(StandardCharsets.UTF_8));
        File dir = makeSpecDir();
        String err = ScoringParamGen.run(new String[]{
                "-i", tsv.getPath(),
                "-d", dir.getPath(),
                "-m", "HCD",
                "-inst", "QExactive",
                "-e", "Tryp"});
        assertEquals("No results to train on. Exiting.", err);
    }

    @Test
    public void mainWithEmptyArgsDoesNotThrow() {
        // Smoke check: empty args triggers printUsageInfo and returns. Should
        // not throw or call System.exit.
        ScoringParamGen.main(new String[]{});
    }

    @Test
    public void mainWithHelpFlagDoesNotThrow() {
        ScoringParamGen.main(new String[]{"-h"});
        ScoringParamGen.main(new String[]{"--help"});
        ScoringParamGen.main(new String[]{"-help"});
    }
}
