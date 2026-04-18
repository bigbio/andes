package edu.ucsd.msjava.msutil;

import java.io.File;

public class DBSearchIOFiles {
    private File specFile;
    private SpecFileFormat specFileFormat;
    private File outputFile;

    /**
     * Per-file precursor mass shift learned by two-pass calibration (P2-cal).
     * Expressed in ppm; defaults to 0.0 (no calibration).
     *
     * The learned shift is the median of (observed - theoretical) / theoretical * 1e6
     * across high-confidence pre-pass PSMs. It is applied later in
     * {@code ScoredSpectraMap} as {@code mass * (1 - shiftPpm * 1e-6)} to
     * remove a systematic positive bias.
     *
     * This field is written once on the orchestrator thread before any
     * {@code ScoredSpectraMap} is constructed for the file, and is read
     * (immutable) by worker threads thereafter. No synchronization needed.
     */
    private double precursorMassShiftPpm = 0.0;

    public DBSearchIOFiles(File specFile, SpecFileFormat specFileFormat, File outputFile) {
        this.specFile = specFile;
        this.specFileFormat = specFileFormat;
        this.outputFile = outputFile;
    }

    public File getSpecFile() {
        return specFile;
    }

    public SpecFileFormat getSpecFileFormat() {
        return specFileFormat;
    }

    public File getOutputFile() {
        return outputFile;
    }

    public double getPrecursorMassShiftPpm() {
        return precursorMassShiftPpm;
    }

    public void setPrecursorMassShiftPpm(double precursorMassShiftPpm) {
        this.precursorMassShiftPpm = precursorMassShiftPpm;
    }
}
