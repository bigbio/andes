package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.Spectrum;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Unit tests for {@link FragmentIndexCandidateGenerator}.
 *
 * <p>Builds a tiny 3-peptide {@link FragmentIndex} (PEPTIDER, SAMPLERK, ANOTHERK)
 * and synthesizes a mock {@link Spectrum} whose peaks sit at exact b/y m/z of
 * one target peptide. Asserts the generator returns the target peptide as the
 * rank-1 candidate and handles K-bounds + empty-spectrum edge cases.
 *
 * <p>Zero file I/O beyond the ionstat param resource that
 * {@link NewScorerFactory} loads from the classpath.
 */
public class TestFragmentIndexCandidateGenerator {

    private static final double SLAB_MIN_MASS = 100.0;
    private static final double SLAB_MAX_MASS = 4000.0;
    private static final double SLAB_WIDTH = 50.0;
    private static final double BOUNDARY_OVERLAP = 0.5;
    private static final double FRAG_BIN_WIDTH = 1.0005;   // 1 Da-ish, matches typical low-res

    private AminoAcidSet aaSet;
    private FragmentIndex index;
    private NewRankScorer rankScorer;
    private List<String> peptides;

    @Before
    public void setUp() {
        aaSet = AminoAcidSet.getStandardAminoAcidSet();
        SlabAssigner assigner = new SlabAssigner(SLAB_MIN_MASS, SLAB_MAX_MASS, SLAB_WIDTH, BOUNDARY_OVERLAP);
        FragmentIndexBuilder builder = new FragmentIndexBuilder(aaSet, assigner, FRAG_BIN_WIDTH);

        // Three hand-crafted peptides whose b/y fragment masses are all distinct.
        // Using K/R-terminal fragments (tryptic-ish) so the mass is realistic.
        // Note: avoid non-standard residues (O, U) — FragmentIndexBuilder silently
        // skips peptides containing characters not in the standard amino-acid set.
        peptides = Arrays.asList("PEPTIDER", "SAMPLERK", "VLLLLLKR");
        index = builder.build(peptides);

        rankScorer = NewScorerFactory.get(
                ActivationMethod.CID,
                InstrumentType.LOW_RESOLUTION_LTQ,
                Enzyme.TRYPSIN,
                Protocol.STANDARD);
    }

    private double precursorMass(String peptide) {
        double sum = Composition.H2O;
        for (int i = 0; i < peptide.length(); i++) {
            AminoAcid aa = aaSet.getAminoAcid(peptide.charAt(i));
            sum += aa.getAccurateMass();
        }
        return sum;
    }

    /**
     * Constructs a Spectrum whose peaks are placed at the exact theoretical
     * b/y m/z values of {@code peptide}, at charge 2. Calls setRanksOfPeaks()
     * so the generator's {@link Peak#getRank()} reads return the intensity rank.
     */
    private Spectrum buildSpectrumFor(String peptide) {
        double pm = precursorMass(peptide);
        int charge = 2;
        float mz = (float) ((pm + charge * Composition.PROTON) / charge);
        Spectrum spec = new Spectrum();
        spec.setPrecursor(new Peak(mz, 1000f, charge));
        TheoreticalFragmentGenerator fgen = new TheoreticalFragmentGenerator(aaSet);
        TheoreticalFragmentGenerator.Fragment[] frags = fgen.fragmentsFor(peptide);
        // Intensities: decrease slightly with index so all ranks are distinct.
        float base = 1000f;
        for (int i = 0; i < frags.length; i++) {
            TheoreticalFragmentGenerator.Fragment f = frags[i];
            spec.add(new Peak((float) f.mass(), base - i, 1));
        }
        spec.setRanksOfPeaks();
        return spec;
    }

    @Test
    public void indexContainsAllThreePeptides() {
        // Sanity: all three land in their slabs.
        Assert.assertEquals(3, index.totalPeptideEntries());
    }

    @Test
    public void topOneReturnsSingleCandidate() {
        Spectrum spec = buildSpectrumFor("PEPTIDER");
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 1);
        Assert.assertEquals("K=1 must return 1 hit", 1, hits.size());
    }

    @Test
    public void topOneIsThePeptideMatchingTheSpectrum() {
        Spectrum spec = buildSpectrumFor("PEPTIDER");
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 5);
        Assert.assertFalse("expected at least one hit", hits.isEmpty());
        CandidateHit top = hits.get(0);
        String seq = index.peptideTable(top.slabId()).sequence(top.localPeptideId());
        Assert.assertEquals("top-1 candidate must be PEPTIDER", "PEPTIDER", seq);
        Assert.assertTrue("top-1 score must be positive (log-ratio of ion vs. noise)",
                top.newRankSum() > 0f);
    }

    @Test
    public void scoresAreDescendingInOutputOrder() {
        Spectrum spec = buildSpectrumFor("PEPTIDER");
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 5);
        for (int i = 1; i < hits.size(); i++) {
            Assert.assertTrue(
                    "hits must be sorted by score descending: " + hits.get(i - 1).newRankSum()
                            + " vs " + hits.get(i).newRankSum(),
                    hits.get(i - 1).newRankSum() >= hits.get(i).newRankSum());
        }
    }

    @Test
    public void kBoundedAboveDoesNotThrowAndCapsAtCandidatePoolSize() {
        // K larger than the number of peptides in the matched slab must return
        // <= slab peptide count (capped by the Tier-1a fingerprint survivors).
        // The FP pre-filter is the correctness hero here: peptides whose
        // fingerprint doesn't overlap the spectrum's fingerprint are excluded
        // regardless of K.
        Spectrum spec = buildSpectrumFor("PEPTIDER");
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 100);
        Assert.assertFalse("must return at least one hit", hits.isEmpty());
        int slabId = hits.get(0).slabId();
        Assert.assertTrue(
                "K=100 must cap at slab peptide count",
                hits.size() <= index.slab(slabId).peptideCount());
    }

    @Test
    public void emptySpectrumReturnsEmptyList() {
        Spectrum spec = new Spectrum();
        double pm = precursorMass("PEPTIDER");
        int charge = 2;
        float mz = (float) ((pm + charge * Composition.PROTON) / charge);
        spec.setPrecursor(new Peak(mz, 1000f, charge));
        // no peaks added; setRanksOfPeaks is a no-op on empty
        spec.setRanksOfPeaks();
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 5);
        Assert.assertTrue("empty spectrum must yield no candidates", hits.isEmpty());
    }

    @Test
    public void zeroKReturnsEmptyList() {
        Spectrum spec = buildSpectrumFor("PEPTIDER");
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 0);
        Assert.assertTrue(hits.isEmpty());
    }

    @Test
    public void peptideMatchingSpectrumScoresStrictlyHigherThanAnyInSameSlab() {
        // Feed a PEPTIDER spectrum; any other peptide in the SAME slab that
        // survives the FP pre-filter must score strictly lower than PEPTIDER.
        Spectrum spec = buildSpectrumFor("PEPTIDER");
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 100);
        Assert.assertFalse(hits.isEmpty());
        CandidateHit top = hits.get(0);
        Assert.assertEquals("PEPTIDER",
                index.peptideTable(top.slabId()).sequence(top.localPeptideId()));
        for (int i = 1; i < hits.size(); i++) {
            CandidateHit other = hits.get(i);
            String seq = index.peptideTable(other.slabId()).sequence(other.localPeptideId());
            Assert.assertTrue(
                    "PEPTIDER must outscore " + seq
                            + " (top=" + top.newRankSum() + ", other=" + other.newRankSum() + ")",
                    top.newRankSum() > other.newRankSum());
        }
    }

    @Test
    public void candidateResolvesToSequenceAndMassViaPeptideTable() {
        Spectrum spec = buildSpectrumFor("PEPTIDER");
        FragmentIndexCandidateGenerator gen = new FragmentIndexCandidateGenerator(index);
        List<CandidateHit> hits = gen.topKForSpectrum(spec, 1);
        Assert.assertEquals(1, hits.size());
        CandidateHit h = hits.get(0);
        PeptideTable table = index.peptideTable(h.slabId());
        Assert.assertEquals("PEPTIDER", table.sequence(h.localPeptideId()));
        Assert.assertEquals(
                precursorMass("PEPTIDER"),
                table.precursorMass(h.localPeptideId()),
                1e-3);
    }
}
