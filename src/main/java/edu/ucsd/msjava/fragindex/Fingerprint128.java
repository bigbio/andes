package edu.ucsd.msjava.fragindex;

/**
 * 128-bit fragment fingerprint split into b-ion / y-ion halves.
 *
 * <p>Each theoretical fragment hashes into one bit in the appropriate half via
 * {@code bucket_index % 64}. At search time a spectrum's fingerprint is ANDed
 * with each candidate peptide's fingerprint and popcounted — peptides whose
 * fragment set doesn't share enough bits with the spectrum are pruned before
 * any fragment-index lookup runs.
 *
 * <p>Threshold tuning (default {@code popcountAnd ≥ 8}) lives in the caller,
 * not this class.
 */
public final class Fingerprint128 {
    private long lo; // b-ion bits
    private long hi; // y-ion bits

    public Fingerprint128() {}

    public Fingerprint128(long lo, long hi) {
        this.lo = lo;
        this.hi = hi;
    }

    public void setBIonBucket(int bucketIndex) {
        lo |= 1L << (bucketIndex & 63);
    }

    public void setYIonBucket(int bucketIndex) {
        hi |= 1L << (bucketIndex & 63);
    }

    public int popcountB() { return Long.bitCount(lo); }
    public int popcountY() { return Long.bitCount(hi); }

    public int popcountAnd(Fingerprint128 other) {
        return Long.bitCount(lo & other.lo) + Long.bitCount(hi & other.hi);
    }

    public long loBits() { return lo; }
    public long hiBits() { return hi; }
}
