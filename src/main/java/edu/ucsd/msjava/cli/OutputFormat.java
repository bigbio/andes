package edu.ucsd.msjava.cli;

/**
 * Search output format selected by {@code -outputFormat}. Picocli matches
 * incoming values case-insensitively (see
 * {@code @Command(caseInsensitiveEnumValuesAllowed = true)}).
 *
 * <p>Numeric forms ({@code 0} / {@code 1}) accepted by older releases are
 * intentionally not supported. Users on legacy invocations should switch
 * to the named values.
 */
public enum OutputFormat {
    /** Percolator {@code .pin} (default). */
    PIN,
    /** Tab-separated values, direct inspection / downstream tools. */
    TSV
}
