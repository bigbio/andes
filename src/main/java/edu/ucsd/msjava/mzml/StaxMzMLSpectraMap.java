package edu.ucsd.msjava.mzml;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;

import java.util.ArrayList;

/**
 * StAX-based implementation of SpectrumAccessorBySpecIndex for mzML files.
 * Drop-in replacement for MzMLSpectraMap (jmzml-based).
 */
public class StaxMzMLSpectraMap implements SpectrumAccessorBySpecIndex {
    private final StaxMzMLParser parser;
    private final int minMSLevel;
    private final int maxMSLevel;

    public StaxMzMLSpectraMap(StaxMzMLParser parser, int minMSLevel, int maxMSLevel) {
        this.parser = parser;
        this.minMSLevel = minMSLevel;
        this.maxMSLevel = maxMSLevel;
    }

    @Override
    public Spectrum getSpectrumBySpecIndex(int specIndex) {
        Spectrum spec = parser.getSpectrumBySpecIndex(specIndex);
        if (spec != null && (spec.getMSLevel() < minMSLevel || spec.getMSLevel() > maxMSLevel))
            return null;
        return spec;
    }

    @Override
    public Spectrum getSpectrumById(String specId) {
        Spectrum spec = parser.getSpectrumById(specId);
        if (spec != null && (spec.getMSLevel() < minMSLevel || spec.getMSLevel() > maxMSLevel))
            return null;
        return spec;
    }

    @Override
    public String getID(int specIndex) {
        return parser.getID(specIndex);
    }

    @Override
    public Float getPrecursorMz(int specIndex) {
        return parser.getPrecursorMz(specIndex);
    }

    @Override
    public String getTitle(int specIndex) {
        return null;
    }

    @Override
    public ArrayList<Integer> getSpecIndexList() {
        return parser.getSpecIndexList(minMSLevel, maxMSLevel);
    }
}
