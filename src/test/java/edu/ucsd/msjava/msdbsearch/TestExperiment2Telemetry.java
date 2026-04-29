package edu.ucsd.msjava.msdbsearch;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class TestExperiment2Telemetry {

    @Before
    public void resetCounters() {
        Experiment2Telemetry.reset();
    }

    @Test
    public void countsEvaluationsAndPrunes() {
        Experiment2Telemetry.recordEvaluation(false);
        Experiment2Telemetry.recordEvaluation(true);
        Experiment2Telemetry.recordEvaluation(true);
        Experiment2Telemetry.recordEvaluation(false);

        assertEquals(4L, Experiment2Telemetry.getPrefixesEvaluated());
        assertEquals(2L, Experiment2Telemetry.getPrefixesPruned());
        assertEquals(0.5, Experiment2Telemetry.pruneRatio(), 1e-9);
    }

    @Test
    public void pruneRatioIsZeroWhenNoEvaluations() {
        assertEquals(0.0, Experiment2Telemetry.pruneRatio(), 0.0);
    }

    @Test
    public void resetClearsCounters() {
        Experiment2Telemetry.recordEvaluation(true);
        Experiment2Telemetry.recordEvaluation(true);
        Experiment2Telemetry.reset();
        assertEquals(0L, Experiment2Telemetry.getPrefixesEvaluated());
        assertEquals(0L, Experiment2Telemetry.getPrefixesPruned());
    }

    @Test
    public void countersAreThreadSafe() throws InterruptedException {
        final int threads = 8;
        final int perThread = 10_000;
        Thread[] workers = new Thread[threads];
        for (int i = 0; i < threads; i++) {
            final boolean prune = (i % 2 == 0);
            workers[i] = new Thread(() -> {
                for (int j = 0; j < perThread; j++) {
                    Experiment2Telemetry.recordEvaluation(prune);
                }
            });
        }
        for (Thread w : workers) w.start();
        for (Thread w : workers) w.join();

        assertEquals((long) threads * perThread, Experiment2Telemetry.getPrefixesEvaluated());
        assertEquals((long) (threads / 2) * perThread, Experiment2Telemetry.getPrefixesPruned());
    }

    @Test
    public void enabledReflectsSystemPropertyAtClassLoad() {
        // ENABLED is captured at class-load time. With the property unset
        // (default in tests), enabled() must be false.
        assertEquals(
                Boolean.parseBoolean(System.getProperty(Experiment2Telemetry.SYSTEM_PROPERTY, "false")),
                Experiment2Telemetry.enabled());
        assertEquals("msgfplus.experiment2Telemetry", Experiment2Telemetry.SYSTEM_PROPERTY);
    }
}
