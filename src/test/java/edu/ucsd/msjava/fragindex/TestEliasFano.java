package edu.ucsd.msjava.fragindex;

import org.junit.Assert;
import org.junit.Test;

public class TestEliasFano {

    @Test
    public void emptyListRoundTrip() {
        byte[] encoded = EliasFano.encode(new int[0]);
        Assert.assertNotNull(encoded);
        int[] decoded = EliasFano.decode(encoded);
        Assert.assertEquals(0, decoded.length);
    }

    @Test
    public void singleValueRoundTrip() {
        int[] original = {42};
        int[] decoded = EliasFano.decode(EliasFano.encode(original));
        Assert.assertArrayEquals(original, decoded);
    }

    @Test
    public void monotonicListRoundTrip() {
        int[] original = {0, 1, 5, 12, 12, 18, 31, 47};
        int[] decoded = EliasFano.decode(EliasFano.encode(original));
        Assert.assertArrayEquals(original, decoded);
    }

    @Test
    public void largeRangeRoundTrip() {
        int[] original = new int[1000];
        for (int i = 0; i < original.length; i++) original[i] = i * 53;
        int[] decoded = EliasFano.decode(EliasFano.encode(original));
        Assert.assertArrayEquals(original, decoded);
    }
}
