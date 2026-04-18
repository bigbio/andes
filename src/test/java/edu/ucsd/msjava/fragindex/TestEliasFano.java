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

    @Test
    public void iteratorMatchesArray() {
        int[] original = {2, 3, 5, 7, 11, 13, 17, 19};
        byte[] encoded = EliasFano.encode(original);
        EliasFano.Cursor it = EliasFano.open(encoded);
        int i = 0;
        while (it.hasNext()) {
            Assert.assertEquals(original[i++], it.next());
        }
        Assert.assertEquals(original.length, i);
    }

    @Test
    public void iteratorOnEmpty() {
        EliasFano.Cursor it = EliasFano.open(EliasFano.encode(new int[0]));
        Assert.assertFalse(it.hasNext());
    }
}
