package msgfplus;

import edu.ucsd.msjava.misc.MSGFLogger;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Method;

public class TestMSGFLogger {

    private ByteArrayOutputStream outBuf;
    private ByteArrayOutputStream errBuf;
    private PrintStream capturedOut;
    private PrintStream capturedErr;

    @Before
    public void captureStreams() throws Exception {
        outBuf = new ByteArrayOutputStream();
        errBuf = new ByteArrayOutputStream();
        capturedOut = new PrintStream(outBuf);
        capturedErr = new PrintStream(errBuf);
        // setStreams is package-private; reflect since the test lives in msgfplus, not misc.
        Method m = MSGFLogger.class.getDeclaredMethod("setStreams", PrintStream.class, PrintStream.class);
        m.setAccessible(true);
        m.invoke(null, capturedOut, capturedErr);
    }

    @After
    public void restoreStreams() throws Exception {
        Method m = MSGFLogger.class.getDeclaredMethod("setStreams", PrintStream.class, PrintStream.class);
        m.setAccessible(true);
        m.invoke(null, System.out, System.err);
        MSGFLogger.setVerbose(false);
    }

    @Test
    public void infoAlwaysPrintsToStdout() {
        MSGFLogger.setVerbose(false);
        MSGFLogger.info("hello");
        Assert.assertTrue(outBuf.toString().contains("hello"));
        Assert.assertEquals("", errBuf.toString());
    }

    @Test
    public void debugIsSuppressedWhenVerboseOff() {
        MSGFLogger.setVerbose(false);
        MSGFLogger.debug("internal chatter");
        Assert.assertEquals("", outBuf.toString());
    }

    @Test
    public void debugPrintsWhenVerboseOn() {
        MSGFLogger.setVerbose(true);
        MSGFLogger.debug("internal chatter");
        Assert.assertTrue(outBuf.toString().contains("internal chatter"));
    }

    @Test
    public void warnGoesToStderrWithPrefix() {
        MSGFLogger.warn("disk getting full");
        Assert.assertTrue(errBuf.toString().contains("[Warning] disk getting full"));
        Assert.assertEquals("", outBuf.toString());
    }

    @Test
    public void errorGoesToStderrWithPrefix() {
        MSGFLogger.error("crashed");
        Assert.assertTrue(errBuf.toString().contains("[Error] crashed"));
    }

    @Test
    public void formatArgumentsAreInterpolated() {
        MSGFLogger.info("hit %d / %d at %.1f%%", 3, 10, 30.0f);
        Assert.assertTrue(outBuf.toString().contains("hit 3 / 10 at 30.0%"));
    }

    @Test
    public void isVerboseReflectsFlag() {
        MSGFLogger.setVerbose(false);
        Assert.assertFalse(MSGFLogger.isVerbose());
        MSGFLogger.setVerbose(true);
        Assert.assertTrue(MSGFLogger.isVerbose());
    }
}
