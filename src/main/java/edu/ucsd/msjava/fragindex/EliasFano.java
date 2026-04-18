package edu.ucsd.msjava.fragindex;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Simple Elias-Fano-inspired codec for sorted non-decreasing int[] lists.
 *
 * Layout (little-endian):
 *   [4 bytes: length N]
 *   [4 bytes: max value U, or 0 if N==0]
 *   [for each value i: 4 bytes raw int]
 *
 * This first cut is correctness-only — plain int array encoding. A compact
 * Elias-Fano layout replaces this in Task 7 once the API shape is stable.
 */
public final class EliasFano {
    private EliasFano() {}

    public static byte[] encode(int[] values) {
        int n = values.length;
        ByteBuffer buf = ByteBuffer.allocate(8 + 4 * n).order(ByteOrder.LITTLE_ENDIAN);
        buf.putInt(n);
        buf.putInt(n == 0 ? 0 : values[n - 1]);
        for (int v : values) buf.putInt(v);
        return buf.array();
    }

    public static int[] decode(byte[] encoded) {
        ByteBuffer buf = ByteBuffer.wrap(encoded).order(ByteOrder.LITTLE_ENDIAN);
        int n = buf.getInt();
        buf.getInt(); // max value; unused in this naive layout
        int[] out = new int[n];
        for (int i = 0; i < n; i++) out[i] = buf.getInt();
        return out;
    }
}
