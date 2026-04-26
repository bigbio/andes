package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.msgf.Tolerance;
import picocli.CommandLine.ITypeConverter;
import picocli.CommandLine.TypeConversionException;

/**
 * Typed precursor mass tolerance: a left and a right
 * {@link Tolerance}. Supports symmetric form ({@code "20ppm"}) and
 * asymmetric form ({@code "0.5Da,2.5Da"}). Both sides must use the
 * same unit and be non-negative.
 */
public final class PrecursorTolerance {

    public final Tolerance left;
    public final Tolerance right;

    private PrecursorTolerance(Tolerance left, Tolerance right) {
        this.left = left;
        this.right = right;
    }

    public static PrecursorTolerance parse(String value) {
        String[] tok = value.split(",");
        Tolerance l, r;
        if (tok.length == 1) {
            l = r = Tolerance.parseToleranceStr(tok[0]);
        } else if (tok.length == 2) {
            l = Tolerance.parseToleranceStr(tok[0]);
            r = Tolerance.parseToleranceStr(tok[1]);
        } else {
            throw new IllegalArgumentException("invalid tolerance value: " + value);
        }
        if (l == null || r == null) {
            throw new IllegalArgumentException("invalid tolerance value: " + value);
        }
        if (l.isTolerancePPM() != r.isTolerancePPM()) {
            throw new IllegalArgumentException("left and right tolerance units must be the same");
        }
        if (l.getValue() < 0 || r.getValue() < 0) {
            throw new IllegalArgumentException("parent mass tolerance must not be negative");
        }
        return new PrecursorTolerance(l, r);
    }

    @Override public String toString() {
        return left.equals(right) ? left.toString() : left + "," + right;
    }

    /** picocli {@link ITypeConverter} that wraps {@link #parse(String)}. */
    public static final class Converter implements ITypeConverter<PrecursorTolerance> {
        @Override public PrecursorTolerance convert(String value) {
            try {
                return parse(value);
            } catch (IllegalArgumentException e) {
                throw new TypeConversionException(e.getMessage());
            }
        }
    }
}
