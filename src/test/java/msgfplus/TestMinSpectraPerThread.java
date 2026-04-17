package msgfplus;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.params.Parameter;
import org.junit.Assert;
import org.junit.Test;

public class TestMinSpectraPerThread {

    private static final String KEY =
            ParamManager.ParamNameEnum.MIN_SPECTRA_PER_THREAD.getKey();

    @Test
    public void defaultIs250() {
        ParamManager pm = new ParamManager("MS-GF+", "test", "test", "java -jar MSGFPlus.jar");
        pm.addMSGFPlusParams();
        Assert.assertEquals(250, pm.getMinSpectraPerThread());
    }

    @Test
    public void overrideAppliesThroughGetter() {
        ParamManager pm = new ParamManager("MS-GF+", "test", "test", "java -jar MSGFPlus.jar");
        pm.addMSGFPlusParams();
        Parameter param = pm.getParameter(KEY);
        Assert.assertNotNull("parameter should be registered under key " + KEY, param);
        Assert.assertNull("'50' should parse as a valid minSpectraPerThread", param.parse("50"));
        Assert.assertEquals(50, pm.getMinSpectraPerThread());
    }

    @Test
    public void rejectsZero() {
        ParamManager pm = new ParamManager("MS-GF+", "test", "test", "java -jar MSGFPlus.jar");
        pm.addMSGFPlusParams();
        Parameter param = pm.getParameter(KEY);
        Assert.assertNotNull(param);
        Assert.assertNotNull("'0' must be rejected (minValue is 1)", param.parse("0"));
    }

    @Test
    @SuppressWarnings("deprecation")
    public void msgfdbEntryPointAlsoRegistersTheParam() {
        ParamManager pm = new ParamManager("MSGFDB", "test", "test", "java -jar MSGFDB.jar");
        pm.addMSGFDBParams();
        Assert.assertEquals(250, pm.getMinSpectraPerThread());
    }
}
