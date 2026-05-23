package msgfplus;

import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msscorer.ScoringParameterGeneratorWithErrors;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpectraContainer;
import edu.ucsd.msjava.msutil.Spectrum;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * End-to-end smoke test for the recovered scorer trainer.
 *
 * The other recovery tests
 * ({@link edu.ucsd.msjava.ui.TestScoringParamGenRecovery},
 *  {@link edu.ucsd.msjava.msutil.TestAnnotatedSpectraRecovery})
 * cover argument parsing and TSV-error paths but stop short of actually
 * invoking
 * {@link ScoringParameterGeneratorWithErrors#generateParameters(SpectraContainer,
 *  SpecDataType, AminoAcidSet, File, boolean, boolean)}.
 *
 * This test closes that gap. It builds a synthetic in-memory
 * {@link SpectraContainer} of annotated spectra (no MGF/mzML on disk, no
 * MS-GF+ search), feeds it directly to the trainer, and verifies that:
 *
 *   1. {@code generateParameters} runs to completion without throwing.
 *   2. A {@code <dataType>.param} file is produced on disk and is non-empty.
 *   3. The companion {@code .param.txt} debug dump is also produced.
 *   4. The binary {@code .param} parses cleanly via {@link NewRankScorer},
 *      proving the on-disk format is round-trippable by the same code that
 *      MS-GF+ uses to load shipped scoring models at search time.
 *
 * Approach rationale (Approach A in the planning notes): synthetic spectra
 * are sufficient because the trainer's invariants are all numeric --- it
 * computes b/y ion frequencies, rank distributions, and precursor offset
 * frequencies from (peptide, peak list) pairs without verifying that the
 * peptide is the "correct" assignment. A grid of distinct tryptic peptides
 * with realistic b/y peak placement at varied intensities exercises every
 * code path in the trainer (partition, precursorOFF, ion-type selection,
 * rank distribution, smoothing, write) at low cost.
 *
 * The {@code MIN_NUM_SPECTRA_PER_PARTITION = 400} threshold inside the
 * trainer is the lower bound on input size, but the trainer also dedupes
 * to {@code numSpecsPerPeptide = 3} spectra per (peptide, charge) pair
 * before partitioning. So the unique-peptide count is the dominant knob:
 * 200 unique peptides * 3 reps = 600 spectra retained, comfortably above
 * the 360 floor. The {@code IonProbability.getIonProb} path further
 * requires {@code numObservedPeaks + numMissingPeaks > 1000} per ion type
 * per segment for that ion type to be selected; with peptides of length
 * 8-12 and full b+y ion ladders, 200 unique peptides produce well over
 * that threshold.
 *
 * We use {@link InstrumentType#LOW_RESOLUTION_LTQ} ("LowRes") because:
 *   - it sets {@code errorScalingFactor = 0} inside the trainer, skipping
 *     the (slow, optional) ion-error and noise-error distribution stages,
 *   - it sets {@code applyDeconvolution = false}, skipping a full container
 *     rebuild,
 *   - the resulting .param file still exercises every required write path
 *     (charge histogram, partition table, precursor OFF, fragment OFF, rank
 *     distributions, no error dist).
 *
 * Runtime target: under 60s on a developer laptop.
 */
public class TestScoringParamGenSmoke {

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    /**
     * Residues used as the variable "body" of synthetic tryptic peptides.
     * Excludes K and R (those are reserved for the C-terminus to keep the
     * generated sequences tryptic) and excludes I/L collisions / U/O
     * non-standard codes. Cysteine is also omitted to avoid the standard
     * AA set's lack of fixed carbamidomethylation distorting masses in a
     * way that the trainer's nominal-mass binning might handle awkwardly.
     */
    private static final char[] BODY_RESIDUES =
            {'A', 'D', 'E', 'F', 'G', 'H', 'L', 'M', 'N', 'P',
             'Q', 'S', 'T', 'V', 'W', 'Y'};

    /**
     * Generates a deterministic pool of {@code count} distinct synthetic
     * tryptic peptides. Each peptide has length 8-12, body residues drawn
     * from {@link #BODY_RESIDUES}, and ends in K or R.
     *
     * The trainer keeps at most {@code numSpecsPerPeptide} (= 3 for
     * fixtures with under 2000 unique peptides) {@code (peptide, charge)}
     * pairs in its training set, so the unique-peptide count is the real
     * upper bound on usable training spectra after dedup. To clear the
     * {@code MIN_NUM_SPECTRA_PER_PARTITION} = 400 partitioning floor we
     * need ~140+ unique peptides (140 * 3 = 420), which this generator
     * supplies easily.
     */
    private static List<String> generateSyntheticTrypticPeptides(int count, long seed) {
        Random rnd = new Random(seed);
        java.util.LinkedHashSet<String> peps = new java.util.LinkedHashSet<String>();
        // Defensive guard: we ask for distinct peptides, so make sure the
        // alphabet is large enough that collisions don't loop forever.
        // 16^8 ~= 4e9 distinct length-8 bodies; we won't realistically
        // exceed that.
        while (peps.size() < count) {
            int len = 8 + rnd.nextInt(5);  // 8..12 residues total
            StringBuilder sb = new StringBuilder(len);
            for (int i = 0; i < len - 1; i++) {
                sb.append(BODY_RESIDUES[rnd.nextInt(BODY_RESIDUES.length)]);
            }
            sb.append(rnd.nextBoolean() ? 'K' : 'R');
            peps.add(sb.toString());
        }
        return new ArrayList<String>(peps);
    }

    /**
     * Writes a synthetic annotated MS2 spectrum into the container.
     *
     * @param peptideStr     the unmodified peptide sequence (e.g. "AAGDPLQNK")
     * @param charge         precursor charge state (typically 2)
     * @param noisePeakCount number of random noise peaks to scatter alongside
     *                       the b/y ions; small (3-5) is enough to keep the
     *                       trainer's noise-distribution stage exercised
     * @param rnd            shared RNG for reproducible noise placement
     */
    private static Spectrum buildSyntheticSpectrum(String peptideStr,
                                                   int charge,
                                                   AminoAcidSet aaSet,
                                                   int noisePeakCount,
                                                   Random rnd) {
        Peptide pep = new Peptide(peptideStr, aaSet);
        float parentMass = pep.getParentMass();              // M = sum(residues) + H2O
        float proton = (float) Composition.ChargeCarrierMass();
        float precursorMz = (parentMass + charge * proton) / charge;

        Spectrum spec = new Spectrum();
        spec.setPrecursor(new Peak(precursorMz, 1.0f, charge));
        spec.setActivationMethod(ActivationMethod.HCD);
        spec.setAnnotation(pep);

        // b ion ladder: prefix mass + proton.
        float prefix = 0f;
        for (int i = 0; i < pep.size() - 1; i++) {
            prefix += pep.get(i).getMass();
            float bMz = prefix + proton;
            // Vary intensities so setRanksOfPeaks produces a non-trivial ranking.
            float intensity = 100f + (i * 17f) + rnd.nextFloat() * 50f;
            spec.add(new Peak(bMz, intensity, 1));
        }

        // y ion ladder: suffix mass + H2O + proton.
        float suffix = 0f;
        float h2oPlusProton = (float) Composition.H2O + proton;
        for (int i = 0; i < pep.size() - 1; i++) {
            suffix += pep.get(pep.size() - 1 - i).getMass();
            float yMz = suffix + h2oPlusProton;
            float intensity = 200f + (i * 13f) + rnd.nextFloat() * 50f;
            spec.add(new Peak(yMz, intensity, 1));
        }

        // A few low-intensity noise peaks to keep the noise distribution
        // (rankDistTable's IonType.NOISE entry, generated by generateRankDist)
        // populated with non-zero counts.
        for (int i = 0; i < noisePeakCount; i++) {
            float mz = 100f + rnd.nextFloat() * (precursorMz - 100f);
            float intensity = 5f + rnd.nextFloat() * 15f;
            spec.add(new Peak(mz, intensity, 1));
        }

        // Peaks must be in ascending m/z order for getPeakByMass binary search.
        spec.setRanksOfPeaks();
        java.util.Collections.sort(spec, new Peak.MassComparator());
        return spec;
    }

    /**
     * Build {@code uniquePeptides * spectraPerPeptide} synthetic annotated
     * spectra at the given charge state.
     *
     * The trainer dedupes to {@code numSpecsPerPeptide} per (peptide,
     * charge) pair internally (= 3 for fixtures under 2000 unique
     * peptides), so {@code spectraPerPeptide} should be 3 to ensure all
     * generated spectra are retained after dedup. {@code uniquePeptides}
     * therefore directly determines how many spectra survive into
     * partition().
     */
    private static SpectraContainer buildSyntheticContainer(int uniquePeptides,
                                                            int spectraPerPeptide,
                                                            int charge,
                                                            AminoAcidSet aaSet,
                                                            long seed) {
        // Deterministic RNG so the test is reproducible across runs / CI.
        Random rnd = new Random(seed);
        List<String> peptides = generateSyntheticTrypticPeptides(uniquePeptides, seed);
        SpectraContainer container = new SpectraContainer();
        for (String pep : peptides) {
            for (int rep = 0; rep < spectraPerPeptide; rep++) {
                container.add(buildSyntheticSpectrum(pep, charge, aaSet, 4, rnd));
            }
        }
        return container;
    }

    @Test
    public void generatesParamFileFromSyntheticAnnotatedSpectra() throws Exception {
        // The trainer's MIN_NUM_SPECTRA_PER_PARTITION (400) is the dominant
        // input-size constraint, and the trainer dedupes to 3 spectra per
        // (peptide, charge) pair before partitioning. 200 unique peptides
        // * 3 reps = 600 retained spectra, comfortably above the floor and
        // well under 60s end-to-end.
        final int uniquePeptides = 200;
        final int spectraPerPeptide = 3;
        final int charge = 2;
        final int totalSpectra = uniquePeptides * spectraPerPeptide;

        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        SpectraContainer container = buildSyntheticContainer(
                uniquePeptides, spectraPerPeptide, charge, aaSet, 42L);

        // Sanity: the container is well-populated and every spectrum is annotated.
        assertEquals(totalSpectra, container.size());
        for (Spectrum s : container) {
            assertNotNull("every synthetic spectrum must carry a Peptide annotation", s.getAnnotation());
            assertEquals(charge, s.getCharge());
            assertTrue("every synthetic spectrum must have b/y peaks", s.size() > 0);
        }

        // LowRes: errorScalingFactor=0, no deconvolution. Tryp + Standard
        // protocol = the simplest happy path through the trainer.
        SpecDataType dataType = new SpecDataType(
                ActivationMethod.HCD,
                InstrumentType.LOW_RESOLUTION_LTQ,
                Enzyme.TRYPSIN,
                Protocol.STANDARD);

        File outDir = tmp.newFolder("trainerOut");

        long t0 = System.currentTimeMillis();
        ScoringParameterGeneratorWithErrors.generateParameters(
                container,
                dataType,
                aaSet,
                outDir,
                /* isText  = */ false,
                /* verbose = */ false);
        long elapsedMs = System.currentTimeMillis() - t0;

        // Soft runtime guard: if we ever drift over ~45s on this 600-spectrum
        // input something has likely regressed badly. The constraint in the
        // test plan is <60s end-to-end.
        assertTrue("trainer took too long: " + elapsedMs + " ms (expected < 45000)",
                elapsedMs < 45_000);

        // The .param file is written as "<dataType>.param" inside outDir; the
        // dataType.toString() form for non-Standard protocols includes a
        // trailing "_<protocol>". Standard protocol omits the suffix per
        // SpecDataType.toString().
        String expectedName = dataType.toString() + ".param";
        File paramFile = new File(outDir, expectedName);
        assertTrue("expected .param file was not produced at " + paramFile.getAbsolutePath(),
                paramFile.exists() && paramFile.isFile());
        assertTrue(".param file is empty", paramFile.length() > 0);

        File paramTextFile = new File(outDir, expectedName + ".txt");
        assertTrue("companion .param.txt debug dump was not produced",
                paramTextFile.exists() && paramTextFile.length() > 0);

        // Round-trip: reload via NewRankScorer. This exercises the same code
        // path MS-GF+ uses at search time to read a shipped .param resource,
        // so a successful constructor here proves the trainer's binary
        // output is structurally valid (charge histogram, partitions,
        // precursor OFF, fragment OFF, rank distributions, error dist
        // section absent because errorScalingFactor=0).
        NewRankScorer reloaded = new NewRankScorer(paramFile.getPath());
        assertNotNull(reloaded);
        assertNotNull("reloaded scorer must expose its dataType", reloaded.getSpecDataType());
        assertEquals("activation method round-tripped",
                ActivationMethod.HCD, reloaded.getSpecDataType().getActivationMethod());
        assertEquals("instrument type round-tripped",
                InstrumentType.LOW_RESOLUTION_LTQ, reloaded.getSpecDataType().getInstrumentType());
        assertEquals("enzyme round-tripped",
                Enzyme.TRYPSIN, reloaded.getSpecDataType().getEnzyme());
        assertNotNull("reloaded scorer must expose a non-null partition set",
                reloaded.getParitionSet());
        assertTrue("reloaded scorer must have at least one partition",
                reloaded.getParitionSet().size() > 0);
    }

    @Test
    public void debugTextDumpReferencesActivationAndEnzyme() throws Exception {
        // Smaller second @Test that re-runs the trainer with a slightly
        // different shape and inspects the human-readable .param.txt to
        // confirm the trainer wrote sensible header lines. This is cheap
        // because the trainer dominates the runtime and we already paid
        // for it above; this test pays for it again with a deliberately
        // smaller input to keep the suite fast.
        final int uniquePeptides = 180;
        final int spectraPerPeptide = 3;
        final int charge = 2;
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        SpectraContainer container = buildSyntheticContainer(
                uniquePeptides, spectraPerPeptide, charge, aaSet, 1337L);

        SpecDataType dataType = new SpecDataType(
                ActivationMethod.CID,
                InstrumentType.LOW_RESOLUTION_LTQ,
                Enzyme.TRYPSIN,
                Protocol.STANDARD);
        File outDir = tmp.newFolder("trainerOut2");

        ScoringParameterGeneratorWithErrors.generateParameters(
                container, dataType, aaSet, outDir, false, false);

        File paramTextFile = new File(outDir, dataType.toString() + ".param.txt");
        assertTrue("companion .param.txt debug dump was not produced",
                paramTextFile.exists() && paramTextFile.length() > 0);

        String contents = new String(java.nio.file.Files.readAllBytes(paramTextFile.toPath()));
        assertTrue(".param.txt should mention the activation method; got first 200 chars: "
                        + contents.substring(0, Math.min(200, contents.length())),
                contents.contains("CID"));
        assertTrue(".param.txt should mention the enzyme",
                contents.contains("Tryp"));
        assertTrue(".param.txt should include a charge histogram section",
                contents.contains("ChargeHistogram"));
    }

    /**
     * Documents the trainer's effective lower bound on input size. After the
     * (peptide, charge) dedup ({@code numSpecsPerPeptide = 3} for fixtures
     * under 2000 unique peptides), the trainer needs ~360 surviving spectra
     * (= {@code MIN_NUM_SPECTRA_PER_PARTITION * 0.9}) at one charge to
     * produce any partition at all; below that, partition() skips the
     * charge entirely and {@code writeParameters} dies via assert.
     *
     * 130 unique peptides * 3 reps = 390 retained spectra -- just clears
     * the floor with margin for the rounding in the {@code 0.9f} factor.
     * The test makes the implicit contract explicit so future readers don't
     * accidentally shrink fixtures below it.
     */
    @Test
    public void trainerHandlesInputSizeJustAboveMinimumPartitionThreshold() throws Exception {
        final int uniquePeptides = 130;
        final int spectraPerPeptide = 3;
        final int charge = 2;
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        SpectraContainer container = buildSyntheticContainer(
                uniquePeptides, spectraPerPeptide, charge, aaSet, 7L);

        SpecDataType dataType = new SpecDataType(
                ActivationMethod.HCD,
                InstrumentType.LOW_RESOLUTION_LTQ,
                Enzyme.TRYPSIN,
                Protocol.STANDARD);
        File outDir = tmp.newFolder("trainerOut3");

        ScoringParameterGeneratorWithErrors.generateParameters(
                container, dataType, aaSet, outDir, false, false);

        File paramFile = new File(outDir, dataType.toString() + ".param");
        assertTrue("trainer should still emit a .param at the input-size boundary",
                paramFile.exists() && paramFile.length() > 0);

        NewRankScorer reloaded = new NewRankScorer(paramFile.getPath());
        assertNotNull(reloaded);
        assertNotNull(reloaded.getParitionSet());
        // With 380 spectra at one charge and numSegments=2, the partition
        // routine produces 2 partitions (one per segment).
        assertTrue("expected at least one partition; got " + reloaded.getParitionSet().size(),
                reloaded.getParitionSet().size() >= 1);
    }
}
