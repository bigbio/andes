package msgfplus;

import static org.junit.Assert.*;

import edu.ucsd.msjava.params.IntRangeParameter;
import edu.ucsd.msjava.params.ParamManager;
import org.junit.Test;

/**
 * Tests for the -msLevel parameter (issue #159).
 * Verifies that MS level filtering is properly wired through ParamManager.
 */
public class TestMSLevelFiltering {

    private ParamManager createParamManager() {
        ParamManager pm = new ParamManager("MS-GF+", "test", "2024.01.01", "test");
        pm.addMSGFPlusParams();
        return pm;
    }

    @Test
    public void testMSLevelParameterExists() {
        ParamManager pm = createParamManager();
        IntRangeParameter msLevel = pm.getMSLevelParameter();
        assertNotNull("MS_LEVEL parameter should exist", msLevel);
    }

    @Test
    public void testMSLevelDefaultIsMS2() {
        ParamManager pm = createParamManager();
        IntRangeParameter msLevel = pm.getMSLevelParameter();
        // Default should be MS2 only (2,2)
        assertEquals("Default min MS level should be 2", 2, (int) msLevel.getMin());
        assertEquals("Default max MS level should be 2", 2, (int) msLevel.getMax());
    }

    @Test
    public void testMSLevelParseSingleValue() {
        ParamManager pm = createParamManager();
        IntRangeParameter msLevel = pm.getMSLevelParameter();
        String err = msLevel.parse("2");
        assertNull("Parsing '2' should succeed", err);
        assertEquals(2, (int) msLevel.getMin());
        assertEquals(2, (int) msLevel.getMax());
    }

    @Test
    public void testMSLevelParseRange() {
        ParamManager pm = createParamManager();
        IntRangeParameter msLevel = pm.getMSLevelParameter();
        String err = msLevel.parse("2,3");
        assertNull("Parsing '2,3' should succeed", err);
        assertEquals(2, (int) msLevel.getMin());
        assertEquals(3, (int) msLevel.getMax());
    }

    @Test
    public void testMSLevelParseMS3Only() {
        ParamManager pm = createParamManager();
        IntRangeParameter msLevel = pm.getMSLevelParameter();
        String err = msLevel.parse("3");
        assertNull("Parsing '3' should succeed", err);
        assertEquals(3, (int) msLevel.getMin());
        assertEquals(3, (int) msLevel.getMax());
    }

    @Test
    public void testMSLevelParseWideRange() {
        ParamManager pm = createParamManager();
        IntRangeParameter msLevel = pm.getMSLevelParameter();
        String err = msLevel.parse("1,5");
        assertNull("Parsing '1,5' should succeed", err);
        assertEquals(1, (int) msLevel.getMin());
        assertEquals(5, (int) msLevel.getMax());
    }
}
