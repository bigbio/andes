package edu.ucsd.msjava.cli;

import picocli.CommandLine.ITypeConverter;
import picocli.CommandLine.TypeConversionException;

/**
 * Inclusive integer range parsed from CLI/config-file syntax
 * {@code "min,max"} or single value {@code "n"} (interpreted as
 * {@code n,n}). Used by {@code -ti}, {@code -msLevel}, {@code -index}.
 */
public record IntRange(int min, int max) {

    public IntRange {
        if (min > max) {
            throw new IllegalArgumentException("min (" + min + ") > max (" + max + ")");
        }
    }

    public static IntRange parse(String value) {
        String[] tok = value.split(",");
        try {
            if (tok.length == 1) {
                int v = Integer.parseInt(tok[0].trim());
                return new IntRange(v, v);
            }
            if (tok.length == 2) {
                return new IntRange(
                        Integer.parseInt(tok[0].trim()),
                        Integer.parseInt(tok[1].trim()));
            }
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("invalid range: " + value, e);
        }
        throw new IllegalArgumentException("invalid range syntax (expected 'min,max' or single int): " + value);
    }

    @Override public String toString() {
        return min == max ? Integer.toString(min) : min + "," + max;
    }

    /** picocli {@link ITypeConverter} that wraps {@link #parse(String)}. */
    public static final class Converter implements ITypeConverter<IntRange> {
        @Override public IntRange convert(String value) {
            try {
                return parse(value);
            } catch (IllegalArgumentException e) {
                throw new TypeConversionException(e.getMessage());
            }
        }
    }
}
