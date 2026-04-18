package edu.ucsd.msjava.fragindex;

import org.junit.Assert;
import org.junit.Test;

public class TestSlabAssigner {

    // Slabs: [500,550), [550,600), [600,650), ... up to [3950,4000).
    // Overlap: 0.5 Da.
    private SlabAssigner newDefaultAssigner() {
        return new SlabAssigner(500.0, 4000.0, 50.0, 0.5);
    }

    @Test
    public void numSlabsMatchesMassRange() {
        SlabAssigner a = newDefaultAssigner();
        Assert.assertEquals(70, a.numSlabs());
        Assert.assertEquals(500.0, a.slabLowMass(0), 1e-9);
        Assert.assertEquals(550.0, a.slabHighMass(0), 1e-9);
        Assert.assertEquals(3950.0, a.slabLowMass(69), 1e-9);
        Assert.assertEquals(4000.0, a.slabHighMass(69), 1e-9);
    }

    @Test
    public void peptideInMiddleOfSlabPlacedInOneSlabOnly() {
        SlabAssigner a = newDefaultAssigner();
        // 525.0 is in slab 0 ([500, 550)), 25 Da away from both boundaries.
        int[] ids = a.slabsFor(525.0);
        Assert.assertArrayEquals(new int[]{0}, ids);
    }

    @Test
    public void peptideNearUpperBoundaryReplicatesIntoNextSlab() {
        SlabAssigner a = newDefaultAssigner();
        // 549.8 is in slab 0, 0.2 Da from upper boundary 550.0 → replicate into slab 1.
        int[] ids = a.slabsFor(549.8);
        Assert.assertArrayEquals(new int[]{0, 1}, ids);
    }

    @Test
    public void peptideNearLowerBoundaryReplicatesIntoPreviousSlab() {
        SlabAssigner a = newDefaultAssigner();
        // 550.3 is in slab 1, 0.3 Da from lower boundary 550.0 → replicate into slab 0.
        int[] ids = a.slabsFor(550.3);
        Assert.assertArrayEquals(new int[]{0, 1}, ids);
    }

    @Test
    public void peptideAtFirstSlabLowerBoundaryDoesNotReplicateBelow() {
        SlabAssigner a = newDefaultAssigner();
        // 500.2 is in slab 0, near lower boundary BUT no slab below to replicate into.
        int[] ids = a.slabsFor(500.2);
        Assert.assertArrayEquals(new int[]{0}, ids);
    }

    @Test
    public void peptideAtLastSlabUpperBoundaryDoesNotReplicateAbove() {
        SlabAssigner a = newDefaultAssigner();
        // 3999.8 is in slab 69, near upper boundary BUT no slab above to replicate into.
        int[] ids = a.slabsFor(3999.8);
        Assert.assertArrayEquals(new int[]{69}, ids);
    }

    @Test
    public void peptideBelowRangeReturnsEmpty() {
        SlabAssigner a = newDefaultAssigner();
        Assert.assertArrayEquals(new int[0], a.slabsFor(499.0));
    }

    @Test
    public void peptideAboveRangeReturnsEmpty() {
        SlabAssigner a = newDefaultAssigner();
        Assert.assertArrayEquals(new int[0], a.slabsFor(4001.0));
    }

    @Test
    public void peptideExactlyOnBoundaryGoesIntoUpperSlab() {
        SlabAssigner a = newDefaultAssigner();
        // 550.0 exactly: primary slab = 1 (floor-based). Distance 0 from lower → replicate into slab 0.
        int[] ids = a.slabsFor(550.0);
        Assert.assertArrayEquals(new int[]{0, 1}, ids);
    }
}
