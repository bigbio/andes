package msgfplus;

import static org.junit.Assert.*;

import edu.ucsd.msjava.params.IntRangeParameter;
import org.junit.Test;

/**
 * Tests for IntRangeParameter, which supports single values and ranges.
 * Part of issue #159: the -msLevel parameter uses IntRangeParameter.
 */
public class TestIntRangeParameter {

    private IntRangeParameter createInclusiveParam() {
        IntRangeParameter p = new IntRangeParameter("test", "Test", "desc");
        p.setMaxInclusive();
        return p;
    }

    @Test
    public void testSingleValue() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("2");
        assertNull("Single value should parse successfully", err);
        assertEquals(2, (int) p.getMin());
        assertEquals(2, (int) p.getMax());
    }

    @Test
    public void testRange() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("2,3");
        assertNull("Range should parse successfully", err);
        assertEquals(2, (int) p.getMin());
        assertEquals(3, (int) p.getMax());
    }

    @Test
    public void testWideRange() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("1,5");
        assertNull(err);
        assertEquals(1, (int) p.getMin());
        assertEquals(5, (int) p.getMax());
    }

    @Test
    public void testSameMinMax() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("3,3");
        assertNull(err);
        assertEquals(3, (int) p.getMin());
        assertEquals(3, (int) p.getMax());
    }

    @Test
    public void testSingleValueExclusiveMaxRejects() {
        // Default constructor has isMaxInclusive=false, so single value "2"
        // produces min=2,max=2 but effective maxNumber=1 < minNumber=2 -> invalid
        IntRangeParameter p = new IntRangeParameter("test", "Test", "desc");
        String err = p.parse("2");
        assertNotNull("Single value with exclusive max should fail", err);
    }

    @Test
    public void testInvalidReversedRange() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("5,2");
        assertNotNull("Reversed range should fail", err);
    }

    @Test
    public void testInvalidTooManyValues() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("1,2,3");
        assertNotNull("Three values should fail", err);
        assertEquals("illegal syntax", err);
    }

    @Test
    public void testInvalidNonNumeric() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("abc");
        assertNotNull("Non-numeric should fail", err);
        assertEquals("not a valid integer or integer range", err);
    }

    @Test
    public void testInvalidEmpty() {
        IntRangeParameter p = createInclusiveParam();
        String err = p.parse("");
        assertNotNull("Empty string should fail", err);
    }
}
