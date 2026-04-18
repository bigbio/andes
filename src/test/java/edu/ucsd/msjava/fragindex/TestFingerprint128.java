package edu.ucsd.msjava.fragindex;

import org.junit.Assert;
import org.junit.Test;

public class TestFingerprint128 {

    @Test
    public void newFingerprintIsEmpty() {
        Fingerprint128 fp = new Fingerprint128();
        Assert.assertEquals(0, fp.popcountAnd(fp));
    }

    @Test
    public void bIonBitSetsLowHalf() {
        Fingerprint128 fp = new Fingerprint128();
        fp.setBIonBucket(3);
        Assert.assertEquals(1, fp.popcountB());
        Assert.assertEquals(0, fp.popcountY());
    }

    @Test
    public void yIonBitSetsHighHalf() {
        Fingerprint128 fp = new Fingerprint128();
        fp.setYIonBucket(3);
        Assert.assertEquals(0, fp.popcountB());
        Assert.assertEquals(1, fp.popcountY());
    }

    @Test
    public void intersectionCountsOnlySharedBits() {
        Fingerprint128 a = new Fingerprint128();
        a.setBIonBucket(1); a.setBIonBucket(2); a.setBIonBucket(3);
        a.setYIonBucket(1); a.setYIonBucket(2);

        Fingerprint128 b = new Fingerprint128();
        b.setBIonBucket(2); b.setBIonBucket(4);
        b.setYIonBucket(1); b.setYIonBucket(5);

        // shared b-ion buckets: {2}; shared y-ion buckets: {1}; total popcount = 2
        Assert.assertEquals(2, a.popcountAnd(b));
    }

    @Test
    public void largeBucketIndicesWrapModulo64() {
        Fingerprint128 fp = new Fingerprint128();
        fp.setBIonBucket(0);
        fp.setBIonBucket(64);  // should collide with bucket 0 mod 64
        Assert.assertEquals(1, fp.popcountB());
    }
}
