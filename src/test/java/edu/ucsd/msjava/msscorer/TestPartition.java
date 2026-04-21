package edu.ucsd.msjava.msscorer;

import org.junit.Assert;
import org.junit.Test;

public class TestPartition {

    @Test
    public void equalPartitionsHaveEqualHashCode() {
        Partition a = new Partition(2, 1234.5f, 1);
        Partition b = new Partition(2, 1234.5f, 1);

        Assert.assertEquals(a, b);
        Assert.assertEquals(a.hashCode(), b.hashCode());
    }

    @Test
    public void hashCodeTracksMutableFields() {
        Partition p = new Partition(2, 1234.5f, 1);
        int initialHash = p.hashCode();

        p.setCharge(3);
        Assert.assertNotEquals(initialHash, p.hashCode());

        int hashAfterCharge = p.hashCode();
        p.setParentMass(1235.5f);
        Assert.assertNotEquals(hashAfterCharge, p.hashCode());

        int hashAfterMass = p.hashCode();
        p.setPosIndex(2);
        Assert.assertNotEquals(hashAfterMass, p.hashCode());
    }
}
