package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import org.junit.Assert;
import org.junit.Test;

public class TestTheoreticalFragmentGenerator {

    @Test
    public void peptideYieldsLminusOneBAndLminusOneYIons() {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        TheoreticalFragmentGenerator gen = new TheoreticalFragmentGenerator(aaSet);
        // "PEPTIDE" has length 7; expect 6 b-ions and 6 y-ions.
        TheoreticalFragmentGenerator.Fragment[] fragments = gen.fragmentsFor("PEPTIDE");
        int bCount = 0, yCount = 0;
        for (TheoreticalFragmentGenerator.Fragment f : fragments) {
            if (f.isB()) bCount++;
            else yCount++;
        }
        Assert.assertEquals(6, bCount);
        Assert.assertEquals(6, yCount);
    }

    @Test
    public void bIonMassesMatchHandComputed() {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        TheoreticalFragmentGenerator gen = new TheoreticalFragmentGenerator(aaSet);
        // b1(P) = 97.05276 + 1.00728 = 98.06004
        // b2(PE) = 97.05276 + 129.04259 + 1.00728 = 227.10263
        // b3(PEP) = 227.10263 + 97.05276 = 324.15539
        double[] expectedB = {98.06004, 227.10263, 324.15539};
        TheoreticalFragmentGenerator.Fragment[] fragments = gen.fragmentsFor("PEP");
        // PEP has length 3 → 2 b-ions (b1, b2). Test uses longer peptide.
        TheoreticalFragmentGenerator.Fragment[] longer = gen.fragmentsFor("PEPT");
        // PEPT has length 4 → 3 b-ions (b1, b2, b3).
        double[] actualB = new double[3];
        int idx = 0;
        for (TheoreticalFragmentGenerator.Fragment f : longer) {
            if (f.isB()) actualB[idx++] = f.mass();
        }
        Assert.assertEquals(3, idx);
        for (int i = 0; i < expectedB.length; i++) {
            Assert.assertEquals("b" + (i + 1), expectedB[i], actualB[i], 0.001);
        }
    }

    @Test
    public void yIonMassesMatchHandComputed() {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        TheoreticalFragmentGenerator gen = new TheoreticalFragmentGenerator(aaSet);
        // For PEPT:
        // y1(T) = 101.04768 + 18.01056 + 1.00728 = 120.06552
        // y2(PT) = 101.04768 + 97.05276 + 18.01056 + 1.00728 = 217.11828
        // y3(EPT) = 101.04768 + 97.05276 + 129.04259 + 18.01056 + 1.00728 = 346.16086
        double[] expectedY = {120.06552, 217.11828, 346.16086};
        TheoreticalFragmentGenerator.Fragment[] fragments = gen.fragmentsFor("PEPT");
        double[] actualY = new double[3];
        int idx = 0;
        for (TheoreticalFragmentGenerator.Fragment f : fragments) {
            if (!f.isB()) actualY[idx++] = f.mass();
        }
        Assert.assertEquals(3, idx);
        for (int i = 0; i < expectedY.length; i++) {
            Assert.assertEquals("y" + (i + 1), expectedY[i], actualY[i], 0.001);
        }
    }

    @Test
    public void positionRecordsTheFragmentIndex() {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        TheoreticalFragmentGenerator gen = new TheoreticalFragmentGenerator(aaSet);
        TheoreticalFragmentGenerator.Fragment[] fragments = gen.fragmentsFor("PEPTIDE");
        // every b-ion position must be between 1 and L-1 = 6
        // every y-ion position must be between 1 and L-1 = 6
        for (TheoreticalFragmentGenerator.Fragment f : fragments) {
            Assert.assertTrue(f.position() >= 1 && f.position() <= 6);
        }
    }
}
