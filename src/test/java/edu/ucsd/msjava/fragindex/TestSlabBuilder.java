package edu.ucsd.msjava.fragindex;

import org.junit.Assert;
import org.junit.Test;

public class TestSlabBuilder {

    @Test
    public void buildSlabWithTwoPeptidesAndQueryByBucket() {
        SlabBuilder b = new SlabBuilder(/*slabId=*/0, /*minMassDa=*/500.0, /*maxMassDa=*/550.0);

        // peptide 0: b-ion at bucket 10, y-ion at bucket 20
        int p0 = b.addPeptide(/*precursorMassDa=*/510.0);
        b.addFragment(p0, /*bucket=*/10, /*isB=*/true);
        b.addFragment(p0, /*bucket=*/20, /*isB=*/false);

        // peptide 1: b-ion at bucket 10 (shared with p0), y-ion at bucket 30
        int p1 = b.addPeptide(520.0);
        b.addFragment(p1, 10, true);
        b.addFragment(p1, 30, false);

        Slab slab = b.finish();
        Assert.assertEquals(2, slab.peptideCount());

        // bucket 10 should contain both peptides; bucket 20 only p0
        int[] bucket10 = slab.peptidesInBucket(10);
        Assert.assertArrayEquals(new int[]{p0, p1}, bucket10);

        int[] bucket20 = slab.peptidesInBucket(20);
        Assert.assertArrayEquals(new int[]{p0}, bucket20);

        // fingerprints should reflect the fragments we added
        Fingerprint128 fp0 = slab.fingerprint(p0);
        Assert.assertEquals(1, fp0.popcountB()); // 1 b-ion
        Assert.assertEquals(1, fp0.popcountY()); // 1 y-ion
    }

    @Test
    public void fingerprintReturnsSnapshot_mutationDoesNotLeakIntoSlab() {
        SlabBuilder b = new SlabBuilder(0, 500.0, 550.0);
        int p0 = b.addPeptide(510.0);
        b.addFragment(p0, 10, true);
        Slab slab = b.finish();

        Fingerprint128 fp = slab.fingerprint(p0);
        fp.setYIonBucket(99); // must not affect slab

        Fingerprint128 fresh = slab.fingerprint(p0);
        Assert.assertEquals("original popcountB preserved", 1, fresh.popcountB());
        Assert.assertEquals("slab unchanged by external mutation", 0, fresh.popcountY());
    }

    @Test(expected = IllegalStateException.class)
    public void finishTwiceThrows() {
        SlabBuilder b = new SlabBuilder(0, 500.0, 550.0);
        b.addPeptide(510.0);
        b.finish();
        b.finish();
    }

    @Test(expected = IllegalStateException.class)
    public void addPeptideAfterFinishThrows() {
        SlabBuilder b = new SlabBuilder(0, 500.0, 550.0);
        b.finish();
        b.addPeptide(520.0);
    }

    @Test(expected = IllegalStateException.class)
    public void addFragmentAfterFinishThrows() {
        SlabBuilder b = new SlabBuilder(0, 500.0, 550.0);
        int p0 = b.addPeptide(510.0);
        b.finish();
        b.addFragment(p0, 10, true);
    }
}
