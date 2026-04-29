package edu.ucsd.msjava.msdbsearch;

import java.util.concurrent.atomic.LongAdder;

/**
 * Opt-in counters for Experiment 2 (exact prefix mass-interval pruning)
 * Checkpoint 1. Records two aggregates across all worker tasks:
 *
 * <ul>
 *   <li>{@code prefixesEvaluated} — number of times the SA walk reached the
 *       pruning hook (i.e. addResidue succeeded and we considered whether
 *       this branch can produce any spectrum match).</li>
 *   <li>{@code prefixesPruned} — number of those evaluations where the
 *       reachable final-mass interval did not intersect any spectrum
 *       window, so the branch could be safely killed.</li>
 * </ul>
 *
 * <p>Checkpoint 1 reports the would-be-prune ratio without actually
 * breaking out of the SA walk. The decision gate (per
 * {@code experiment-2-mass-interval-pruning.md} §6):
 * ratio ≥ 5 % → proceed to Checkpoint 2 with the actual {@code break};
 * ratio &lt; 1 % → bookkeeping cost likely exceeds savings; kill.
 *
 * <p>Enable via {@code -Dmsgfplus.experiment2Telemetry=true}. Off by
 * default; OFF-mode is bit-identical (single load+branch when disabled).
 * Mirrors the {@link PhaseBTelemetry} pattern.
 */
public final class Experiment2Telemetry {

    static final String SYSTEM_PROPERTY = "msgfplus.experiment2Telemetry";
    static final String PRUNING_PROPERTY = "msgfplus.experiment2Pruning";

    private static final boolean ENABLED =
            Boolean.parseBoolean(System.getProperty(SYSTEM_PROPERTY, "false"));
    /** Checkpoint 2: when true, the bound test in {@code DBScanner.dbSearch}
     *  actually breaks out of the residue-extension loop instead of just
     *  recording would-be prunes. Independent of {@link #ENABLED}; either or
     *  both can be set. Default: off (Checkpoint 1 telemetry only). */
    private static final boolean PRUNING_ENABLED =
            Boolean.parseBoolean(System.getProperty(PRUNING_PROPERTY, "false"));

    private static final LongAdder prefixesEvaluated = new LongAdder();
    private static final LongAdder prefixesPruned = new LongAdder();

    private Experiment2Telemetry() {}

    public static boolean enabled() {
        return ENABLED;
    }

    /** Returns true when {@code -Dmsgfplus.experiment2Pruning=true} —
     *  i.e. the bound test should break out of the SA walk on a hit. */
    public static boolean pruningEnabled() {
        return PRUNING_ENABLED;
    }

    /** True when the bound must be computed at all (either for telemetry
     *  or for actual pruning). Used to short-circuit OFF-mode cleanly. */
    public static boolean boundComputationActive() {
        return ENABLED || PRUNING_ENABLED;
    }

    public static void recordEvaluation(boolean wouldPrune) {
        prefixesEvaluated.increment();
        if (wouldPrune) prefixesPruned.increment();
    }

    public static long getPrefixesEvaluated() {
        return prefixesEvaluated.sum();
    }

    public static long getPrefixesPruned() {
        return prefixesPruned.sum();
    }

    public static double pruneRatio() {
        long evaluated = prefixesEvaluated.sum();
        if (evaluated == 0) return 0.0;
        return (double) prefixesPruned.sum() / evaluated;
    }

    public static void reset() {
        prefixesEvaluated.reset();
        prefixesPruned.reset();
    }
}
