package edu.ucsd.msjava.fragindex;

/**
 * Immutable read-only view over one precursor-mass slab of the fragment index.
 *
 * <p>Returned by {@link SlabBuilder#finish()} once all peptides and fragments
 * are loaded. Immutable by construction: the fingerprint array is never
 * mutated after construction, and {@link #fingerprint(int)} returns a fresh
 * snapshot rather than the internal object. Safe for concurrent readers.
 */
public final class Slab {
    private final int slabId;
    private final double minMassDa;
    private final double maxMassDa;
    private final int peptideCount;
    private final Fingerprint128[] fingerprints;
    private final byte[][] bucketEncoded;   // bucket -> Elias-Fano-encoded peptide-id list

    Slab(int slabId, double minMassDa, double maxMassDa,
         Fingerprint128[] fingerprints, byte[][] bucketEncoded) {
        this.slabId = slabId;
        this.minMassDa = minMassDa;
        this.maxMassDa = maxMassDa;
        this.peptideCount = fingerprints.length;
        this.fingerprints = fingerprints;
        this.bucketEncoded = bucketEncoded;
    }

    public int slabId() { return slabId; }
    public double minMassDa() { return minMassDa; }
    public double maxMassDa() { return maxMassDa; }
    public int peptideCount() { return peptideCount; }

    /**
     * Returns the fingerprint bits for the given peptide as an immutable
     * 2-long snapshot. The returned Fingerprint128 is a fresh object built
     * from the peptide's lo/hi bit-words; mutating it has no effect on the
     * slab's internal state. Callers that only need bit-level AND+popcount
     * can use {@link #fingerprintLoBits(int)} / {@link #fingerprintHiBits(int)}
     * for zero-allocation access.
     */
    public Fingerprint128 fingerprint(int peptideId) {
        Fingerprint128 src = fingerprints[peptideId];
        return new Fingerprint128(src.loBits(), src.hiBits());
    }

    /** Zero-allocation read of the b-ion fingerprint word for a peptide. */
    public long fingerprintLoBits(int peptideId) {
        return fingerprints[peptideId].loBits();
    }

    /** Zero-allocation read of the y-ion fingerprint word for a peptide. */
    public long fingerprintHiBits(int peptideId) {
        return fingerprints[peptideId].hiBits();
    }

    public int[] peptidesInBucket(int bucket) {
        if (bucket < 0 || bucket >= bucketEncoded.length) return new int[0];
        byte[] enc = bucketEncoded[bucket];
        if (enc == null) return new int[0];
        return EliasFano.decode(enc);
    }

    public EliasFano.Cursor bucketCursor(int bucket) {
        if (bucket < 0 || bucket >= bucketEncoded.length || bucketEncoded[bucket] == null) {
            return EliasFano.open(EliasFano.encode(new int[0]));
        }
        return EliasFano.open(bucketEncoded[bucket]);
    }
}
