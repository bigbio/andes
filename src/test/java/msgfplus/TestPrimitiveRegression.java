package msgfplus;

import edu.ucsd.msjava.msgf.GeneratingFunction;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.PrimitiveAminoAcidGraph;
import edu.ucsd.msjava.msgf.PrimitiveGeneratingFunction;
import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msgf.FlexAminoAcidGraph;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Modification;
import edu.ucsd.msjava.msutil.Peak;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

public class TestPrimitiveRegression {

    private static final class StubScoredSpectrum implements ScoredSpectrum<NominalMass> {
        @Override
        public int getNodeScore(NominalMass prm, NominalMass srm) {
            return 0;
        }

        @Override
        public float getNodeScore(NominalMass node, boolean isPrefix) {
            return 0;
        }

        @Override
        public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float edgeMass) {
            return 0;
        }

        @Override
        public boolean getMainIonDirection() {
            return true;
        }

        @Override
        public Peak getPrecursorPeak() {
            return new Peak(500.0f, 1.0f, 2);
        }

        @Override
        public ActivationMethod[] getActivationMethodArr() {
            return new ActivationMethod[]{ActivationMethod.CID};
        }

        @Override
        public int[] getScanNumArr() {
            return new int[]{1};
        }
    }

    @Test
    public void testPrimitiveGraphSupportsNegativeNominalMassStates() {
        Modification negativeTermMod = Modification.register("TestNegativeNTerm", -200.0);
        ArrayList<Modification.Instance> mods = new ArrayList<>();
        mods.add(new Modification.Instance(negativeTermMod, '*', Modification.Location.N_Term));
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods);

        StubScoredSpectrum scoredSpectrum = new StubScoredSpectrum();
        FlexAminoAcidGraph legacyGraph = new FlexAminoAcidGraph(aaSet, 100, Enzyme.TRYPSIN, scoredSpectrum, false, false);
        PrimitiveAminoAcidGraph primitiveGraph = new PrimitiveAminoAcidGraph(aaSet, 100, Enzyme.TRYPSIN, scoredSpectrum, false, false);

        boolean legacyHasNegativeNode = false;
        for (NominalMass node : legacyGraph.getIntermediateNodeList()) {
            if (node.getNominalMass() < 0) {
                legacyHasNegativeNode = true;
                break;
            }
        }

        boolean primitiveHasNegativeNode = Arrays.stream(primitiveGraph.getActiveNodes()).anyMatch(mass -> mass < 0);

        Assert.assertTrue("Legacy graph should include a negative nominal-mass state", legacyHasNegativeNode);
        Assert.assertTrue("Primitive graph should preserve negative nominal-mass states", primitiveHasNegativeNode);
    }

    @Test
    public void testPrimitiveGeneratingFunctionMatchesLegacyWithNegativeNominalMassStates() {
        Modification negativeTermMod = Modification.register("TestNegativeGFNTerm", -200.0);
        ArrayList<Modification.Instance> mods = new ArrayList<>();
        mods.add(new Modification.Instance(negativeTermMod, '*', Modification.Location.N_Term));
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods);

        StubScoredSpectrum scoredSpectrum = new StubScoredSpectrum();
        FlexAminoAcidGraph legacyGraph = new FlexAminoAcidGraph(aaSet, 100, Enzyme.TRYPSIN, scoredSpectrum, false, false);
        PrimitiveAminoAcidGraph primitiveGraph = new PrimitiveAminoAcidGraph(aaSet, 100, Enzyme.TRYPSIN, scoredSpectrum, false, false);

        GeneratingFunction<NominalMass> legacyGF = new GeneratingFunction<>(legacyGraph).doNotBacktrack();
        PrimitiveGeneratingFunction primitiveGF = new PrimitiveGeneratingFunction(primitiveGraph);

        Assert.assertTrue("Legacy GF should compute successfully", legacyGF.computeGeneratingFunction());
        Assert.assertTrue("Primitive GF should compute successfully", primitiveGF.computeGeneratingFunction());
        Assert.assertEquals("Primitive graph should keep the source node first for DP ordering", 0, primitiveGraph.getActiveNodes()[0]);
        Assert.assertEquals("Primitive GF min score should match the legacy GF", legacyGF.getMinScore(), primitiveGF.getMinScore());
        Assert.assertEquals("Primitive GF max score should match the legacy GF", legacyGF.getMaxScore(), primitiveGF.getMaxScore());
        Assert.assertEquals("Primitive GF spectral probability should match the legacy GF at score 0",
                legacyGF.getSpectralProbability(0),
                primitiveGF.getSpectralProbability(0),
                1.0e-12);
    }
}
