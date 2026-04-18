package edu.ucsd.msjava.fragindex;

import org.junit.Assert;
import org.junit.Test;

public class TestDirectStore {

    @Test
    public void putAndOpenSlab() {
        DirectStore store = new DirectStore(/*slabCount=*/2);

        SlabBuilder b0 = new SlabBuilder(0, 500.0, 550.0);
        int pid = b0.addPeptide(510.0);
        b0.addFragment(pid, 10, true);
        store.putSlab(0, b0.finish());

        Slab read = store.openSlab(0);
        Assert.assertEquals(0, read.slabId());
        Assert.assertEquals(1, read.peptideCount());
        Assert.assertArrayEquals(new int[]{0}, read.peptidesInBucket(10));
    }

    @Test
    public void openUnsetSlabReturnsNull() {
        DirectStore store = new DirectStore(2);
        Assert.assertNull(store.openSlab(1));
    }
}
