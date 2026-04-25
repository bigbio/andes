package edu.ucsd.msjava.msdbsearch;

import org.junit.Assert;
import org.junit.Test;

import java.util.Collections;
import java.util.concurrent.atomic.AtomicInteger;

public class TestConcurrentMSGFPlus {

    @Test
    public void defersScoredSpectraMapConstructionUntilRun() {
        AtomicInteger buildCount = new AtomicInteger();
        ConcurrentMSGFPlus.RunMSGFPlus task = new ConcurrentMSGFPlus.RunMSGFPlus(
                () -> {
                    buildCount.incrementAndGet();
                    throw new IllegalStateException("sentinel");
                },
                null,
                null,
                Collections.<MSGFPlusMatch>emptyList(),
                1
        );

        Assert.assertEquals(0, buildCount.get());

        try {
            task.run();
            Assert.fail("Expected the ScoredSpectraMap supplier to run inside run().");
        } catch (IllegalStateException expected) {
            Assert.assertEquals("sentinel", expected.getMessage());
        }

        Assert.assertEquals(1, buildCount.get());
    }
}
