package edu.ucsd.msjava.fragindex;

/**
 * Elias-Fano compression for sorted (monotonically non-decreasing) int[] lists.
 * Used by the fragment index to store peptide-id lists per fragment-mass bucket
 * at ~0.5-1 byte per entry.
 *
 * Empty-list case handled here; monotonic-list encoding lands in a later task.
 */
public final class EliasFano {
    private EliasFano() {}

    public static byte[] encode(int[] values) {
        if (values.length == 0) return new byte[]{0, 0, 0, 0}; // length prefix only
        throw new UnsupportedOperationException("non-empty encode not implemented yet");
    }

    public static int[] decode(byte[] encoded) {
        int len = (encoded[0] & 0xff) | ((encoded[1] & 0xff) << 8)
                | ((encoded[2] & 0xff) << 16) | ((encoded[3] & 0xff) << 24);
        if (len == 0) return new int[0];
        throw new UnsupportedOperationException("non-empty decode not implemented yet");
    }
}
