package edu.ucsd.msjava.msdbsearch;

import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
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
                1
        );

        Assert.assertEquals(0, buildCount.get());
        Assert.assertNotNull("Per-task result buffer must exist before run()", task.getResults());
        Assert.assertTrue("Per-task result buffer starts empty", task.getResults().isEmpty());

        try {
            task.run();
            Assert.fail("Expected the ScoredSpectraMap supplier to run inside run().");
        } catch (IllegalStateException expected) {
            Assert.assertEquals("sentinel", expected.getMessage());
        }

        Assert.assertEquals(1, buildCount.get());
    }

    @Test
    public void drainsTaskLocalResultsIntoCallerBuffer() {
        ConcurrentMSGFPlus.RunMSGFPlus task = new ConcurrentMSGFPlus.RunMSGFPlus(
                () -> null,
                null,
                null,
                1
        );

        task.getResults().add(null);
        task.getResults().add(null);

        List<MSGFPlusMatch> merged = new ArrayList<>();
        task.drainResultsTo(merged);

        Assert.assertEquals(2, merged.size());
        Assert.assertTrue("Drain should clear the task-local buffer", task.getResults().isEmpty());
    }
}
