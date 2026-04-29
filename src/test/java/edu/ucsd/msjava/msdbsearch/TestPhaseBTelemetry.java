package edu.ucsd.msjava.msdbsearch;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TestPhaseBTelemetry {

    @Before
    public void resetCounters() {
        PhaseBTelemetry.reset();
    }

    @Test
    public void countsPairingCallsAndMatchedKeys() {
        PhaseBTelemetry.recordPairing(3);
        PhaseBTelemetry.recordPairing(5);
        PhaseBTelemetry.recordPairing(0);  // zero-matched calls still count

        assertEquals(3L, PhaseBTelemetry.getPairingCalls());
        assertEquals(8L, PhaseBTelemetry.getMatchedSpecKeys());
        assertEquals(8.0 / 3.0, PhaseBTelemetry.meanMatchedPerCall(), 1e-9);
    }

    @Test
    public void meanIsZeroWhenNoCallsRecorded() {
        assertEquals(0.0, PhaseBTelemetry.meanMatchedPerCall(), 0.0);
    }

    @Test
    public void resetClearsCounters() {
        PhaseBTelemetry.recordPairing(7);
        PhaseBTelemetry.reset();

        assertEquals(0L, PhaseBTelemetry.getPairingCalls());
        assertEquals(0L, PhaseBTelemetry.getMatchedSpecKeys());
        assertEquals(0.0, PhaseBTelemetry.meanMatchedPerCall(), 0.0);
    }

    @Test
    public void countersAreThreadSafe() throws InterruptedException {
        final int threads = 8;
        final int perThread = 10_000;
        Thread[] workers = new Thread[threads];
        for (int i = 0; i < threads; i++) {
            workers[i] = new Thread(() -> {
                for (int j = 0; j < perThread; j++) {
                    PhaseBTelemetry.recordPairing(2);
                }
            });
        }
        for (Thread w : workers) w.start();
        for (Thread w : workers) w.join();

        assertEquals((long) threads * perThread, PhaseBTelemetry.getPairingCalls());
        assertEquals((long) threads * perThread * 2, PhaseBTelemetry.getMatchedSpecKeys());
    }

    @Test
    public void enabledIsControlledBySystemProperty() {
        // The static ENABLED is captured at class-load time. We can't reliably
        // toggle it after the fact in a single JVM, but we can at least verify
        // the contract: when the property is unset (the test default), the
        // method returns false. This is the no-op invariant the recordPairing
        // call site relies on for OFF-mode bit-identical behaviour.
        assertEquals("PhaseBTelemetry should be disabled when -Dmsgfplus.phaseBTelemetry is unset",
                Boolean.parseBoolean(System.getProperty(PhaseBTelemetry.SYSTEM_PROPERTY, "false")),
                PhaseBTelemetry.enabled());
        // Sanity: the SYSTEM_PROPERTY constant is the documented name.
        assertEquals("msgfplus.phaseBTelemetry", PhaseBTelemetry.SYSTEM_PROPERTY);
        // Sanity: after enabling, recordPairing still works (purely additive).
        PhaseBTelemetry.recordPairing(1);
        assertTrue(PhaseBTelemetry.getPairingCalls() >= 1);
    }
}
