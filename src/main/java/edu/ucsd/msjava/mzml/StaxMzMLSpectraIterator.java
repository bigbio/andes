package edu.ucsd.msjava.mzml;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.mgf.SpectrumParser;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * StAX-based mzML spectrum iterator with MS level filtering.
 * Drop-in replacement for MzMLSpectraIterator (jmzml-based).
 */
public class StaxMzMLSpectraIterator implements Iterator<Spectrum>, Iterable<Spectrum> {
    private final Iterator<Spectrum> delegate;
    private Spectrum currentSpectrum;
    private boolean hasNext;
    private long negativePolarityWarningCount = 0;

    public StaxMzMLSpectraIterator(StaxMzMLParser parser, int minMSLevel, int maxMSLevel) {
        this.delegate = parser.iterator(minMSLevel, maxMSLevel);
        this.currentSpectrum = delegate.hasNext() ? delegate.next() : null;
        this.hasNext = currentSpectrum != null;
    }

    @Override
    public boolean hasNext() {
        return hasNext;
    }

    @Override
    public Spectrum next() {
        if (!hasNext) throw new NoSuchElementException("No more spectra");

        Spectrum cur = currentSpectrum;
        currentSpectrum = delegate.hasNext() ? delegate.next() : null;
        if (currentSpectrum == null) hasNext = false;

        if (cur.getScanPolarity() == Spectrum.Polarity.NEGATIVE) {
            warnNegativePolarity(cur);
        }
        return cur;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("StaxMzMLSpectraIterator.remove() not implemented");
    }

    @Override
    public Iterator<Spectrum> iterator() {
        return this;
    }

    private void warnNegativePolarity(Spectrum spec) {
        negativePolarityWarningCount++;
        if (negativePolarityWarningCount > SpectrumParser.MAX_NEGATIVE_POLARITY_WARNINGS)
            return;

        if (negativePolarityWarningCount == 1) {
            System.out.println("Warning: negative polarity spectrum found; you likely need to use a negative charge carrier");
        }
        System.out.println("Negative polarity spectrum found, scan " + spec.getScanNum());

        if (negativePolarityWarningCount == SpectrumParser.MAX_NEGATIVE_POLARITY_WARNINGS) {
            System.out.println("Additional warnings regarding negative polarity will not be shown");
        }
    }
}
