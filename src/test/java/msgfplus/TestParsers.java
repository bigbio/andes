package msgfplus;

import java.io.File;
import java.net.URISyntaxException;
import java.util.Iterator;

import org.junit.Assert;
import org.junit.Test;

import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.mzml.StaxMzMLParser;

import javax.xml.stream.XMLStreamException;
import java.io.IOException;

public class TestParsers {

    @Test
    public void testReadingMgf() throws URISyntaxException {
        File mgfFile = new File(TestParsers.class.getClassLoader().getResource("test.mgf").toURI());
        SpectraAccessor specAcc = new SpectraAccessor(mgfFile);
        Iterator<Spectrum> itr = specAcc.getSpecItr();
        int numSpecs = 0;
        while(itr.hasNext()) {
            itr.next();
            numSpecs++;
        }
        Assert.assertTrue(numSpecs == 5760);
    }

    @Test
    public void testMzML() throws URISyntaxException, IOException, XMLStreamException {
        File mzMLFile = new File(TestParsers.class.getClassLoader().getResource("tiny.pwiz.mzML").toURI());
        StaxMzMLParser parser = new StaxMzMLParser(mzMLFile);
        Assert.assertTrue("Should have at least 1 spectrum", parser.getSpectrumCount() > 0);
    }

    @Test
    public void testMzMLSpectraAccessor() throws URISyntaxException {
        File mzMLFile = new File(TestParsers.class.getClassLoader().getResource("tiny.pwiz.mzML").toURI());
        SpectraAccessor specAcc = new SpectraAccessor(mzMLFile);
        Iterator<Spectrum> itr = specAcc.getSpecItr();
        int numSpecs = 0;
        while(itr.hasNext()) {
            itr.next();
            numSpecs++;
        }
        Assert.assertTrue("Should parse spectra from mzML", numSpecs > 0);
    }

}
