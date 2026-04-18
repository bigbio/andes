package edu.ucsd.msjava.fragindex;

import org.junit.Assert;
import org.junit.Test;

public class TestPeptideTable {

    @Test
    public void emptyTableHasZeroSize() {
        PeptideTable t = new PeptideTable();
        Assert.assertEquals(0, t.size());
    }

    @Test
    public void addPeptideAssignsSequentialIds() {
        PeptideTable t = new PeptideTable();
        int p0 = t.addPeptide("PEPTIDE", 799.36);
        int p1 = t.addPeptide("KEEPSAKE", 902.45);
        int p2 = t.addPeptide("SAMPLE", 594.27);
        Assert.assertEquals(0, p0);
        Assert.assertEquals(1, p1);
        Assert.assertEquals(2, p2);
        Assert.assertEquals(3, t.size());
    }

    @Test
    public void retrieveBackSequenceAndMass() {
        PeptideTable t = new PeptideTable();
        t.addPeptide("PEPTIDE", 799.36);
        Assert.assertEquals("PEPTIDE", t.sequence(0));
        Assert.assertEquals(799.36, t.precursorMass(0), 1e-6);
    }

    @Test
    public void multipleEntriesKeepTheirOwnFields() {
        PeptideTable t = new PeptideTable();
        t.addPeptide("A", 100.0);
        t.addPeptide("BC", 200.5);
        t.addPeptide("DEF", 350.75);
        Assert.assertEquals("A", t.sequence(0));
        Assert.assertEquals("BC", t.sequence(1));
        Assert.assertEquals("DEF", t.sequence(2));
        Assert.assertEquals(100.0, t.precursorMass(0), 1e-9);
        Assert.assertEquals(200.5, t.precursorMass(1), 1e-9);
        Assert.assertEquals(350.75, t.precursorMass(2), 1e-9);
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void outOfBoundsReadThrows() {
        PeptideTable t = new PeptideTable();
        t.addPeptide("A", 100.0);
        t.sequence(1);
    }
}
