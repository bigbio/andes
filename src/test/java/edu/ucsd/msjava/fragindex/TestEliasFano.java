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
}
