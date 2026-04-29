package edu.ucsd.msjava.msdbsearch;

import java.util.concurrent.atomic.LongAdder;

/**
 * Opt-in counter for Phase B (calibrated precursor-window tightening) verification.
 *
 * <p>Records two aggregate metrics across all worker tasks:
 * <ul>
 *   <li>{@code pairingCalls} — number of times {@code DBScanner} hit the
 *       {@code pepMassSpecKeyMap.subMap(leftThr, rightThr)} pairing site for
 *       a candidate peptide.</li>
 *   <li>{@code matchedSpecKeys} — total number of SpecKeys returned across
 *       those pairing calls. Mean per-call = matched / pairingCalls reflects
 *       the post-tightening pairing fan-out the plan asks us to verify.</li>
 * </ul>
 *
 * <p>Enable via {@code -Dmsgfplus.phaseBTelemetry=true}. Off by default; OFF
 * mode is bit-identical (the {@code if (enabled())} guard short-circuits to
 * a single load+branch). Intentionally not a CLI flag: this is a developer
 * diagnostic for the Phase B retrospective, not a user feature.
 *
 * <p>Designed to live one-instance-per-JVM since each {@code java -jar
 * MSGFPlus.jar} invocation is its own process. Tests should call
 * {@link #reset()} between cases.
 */
public final class PhaseBTelemetry {

    static final String SYSTEM_PROPERTY = "msgfplus.phaseBTelemetry";

    private static final boolean ENABLED =
            Boolean.parseBoolean(System.getProperty(SYSTEM_PROPERTY, "false"));

    private static final LongAdder pairingCalls = new LongAdder();
    private static final LongAdder matchedSpecKeys = new LongAdder();

    private PhaseBTelemetry() {}

    public static boolean enabled() {
        return ENABLED;
    }

    /** Records one pairing call and the size of its result set. */
    public static void recordPairing(int matched) {
        pairingCalls.increment();
        matchedSpecKeys.add(matched);
    }

    public static long getPairingCalls() {
        return pairingCalls.sum();
    }

    public static long getMatchedSpecKeys() {
        return matchedSpecKeys.sum();
    }

    /** Mean matched SpecKeys per pairing call, or 0.0 if no calls recorded. */
    public static double meanMatchedPerCall() {
        long calls = pairingCalls.sum();
        if (calls == 0) return 0.0;
        return (double) matchedSpecKeys.sum() / calls;
    }

    /** Tests should call this between cases since the counters are static. */
    public static void reset() {
        pairingCalls.reset();
        matchedSpecKeys.reset();
    }
}
