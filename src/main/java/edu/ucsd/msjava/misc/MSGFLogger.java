package edu.ucsd.msjava.misc;

import java.io.PrintStream;

/**
 * Lightweight leveled logger for MS-GF+ console output.
 *
 * <p>The runtime verbose flag (from {@code -verbose 0/1}) gates {@link #debug}; all other
 * levels print unconditionally. Call {@link #setVerbose(boolean)} once at startup after
 * parsing CLI arguments; the default is {@code false} (compatible with today's behaviour).
 *
 * <p>Designed to replace ad-hoc {@code System.out.println} calls at the top-level entry
 * points without pulling in slf4j / log4j. Info/debug write to {@code stdout}; warn/error
 * write to {@code stderr}.
 */
public final class MSGFLogger {

    private static volatile boolean verbose = false;
    private static PrintStream out = System.out;
    private static PrintStream err = System.err;

    private MSGFLogger() {}

    public static void setVerbose(boolean v) {
        verbose = v;
    }

    public static boolean isVerbose() {
        return verbose;
    }

    /** Testing hook: swap the output streams. Package-private. */
    static void setStreams(PrintStream outStream, PrintStream errStream) {
        out = outStream;
        err = errStream;
    }

    /** Always printed; for top-level progress the user should see. */
    public static void info(String msg) {
        out.println(msg);
    }

    public static void info(String fmt, Object... args) {
        out.println(String.format(fmt, args));
    }

    /** Printed only when {@code -verbose 1}. Use for per-thread / per-task chatter. */
    public static void debug(String msg) {
        if (verbose) {
            out.println(msg);
        }
    }

    public static void debug(String fmt, Object... args) {
        if (verbose) {
            out.println(String.format(fmt, args));
        }
    }

    /** Always printed to stderr, prefixed with {@code [Warning]}. */
    public static void warn(String msg) {
        err.println("[Warning] " + msg);
    }

    public static void warn(String fmt, Object... args) {
        err.println("[Warning] " + String.format(fmt, args));
    }

    /** Always printed to stderr, prefixed with {@code [Error]}. */
    public static void error(String msg) {
        err.println("[Error] " + msg);
    }

    public static void error(String fmt, Object... args) {
        err.println("[Error] " + String.format(fmt, args));
    }
}
