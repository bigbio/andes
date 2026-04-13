package msgfplus;

import edu.ucsd.msjava.mzml.StaxMzMLParser;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;

import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * Tests for the StAX-based mzML parser.
 * Uses tiny.pwiz.mzML which has 4 spectra:
 *   index 0 (scan=19): MS1, 15 peaks, RT=5.89 min
 *   index 1 (scan=20): MS2, 10 peaks, precursor m/z=445.34, charge=2, CID
 *   index 2 (scan=21): MS1, 0 peaks
 *   index 3 (scan=22): MS1, 15 peaks, RT=42.05 sec
 */
public class TestStaxMzMLParser {

    private File getMzMLFile() throws Exception {
        return new File(getClass().getClassLoader().getResource("tiny.pwiz.mzML").toURI());
    }

    @Test
    public void testSpectrumCount() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Assert.assertEquals("Should have 4 spectra", 4, parser.getSpectrumCount());
    }

    @Test
    public void testSpecIndexList() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        ArrayList<Integer> indices = parser.getSpecIndexList();
        Assert.assertEquals(4, indices.size());
        // 1-based indices
        Assert.assertEquals(Integer.valueOf(1), indices.get(0));
        Assert.assertEquals(Integer.valueOf(2), indices.get(1));
        Assert.assertEquals(Integer.valueOf(3), indices.get(2));
        Assert.assertEquals(Integer.valueOf(4), indices.get(3));
    }

    @Test
    public void testSpecIndexListByMSLevel() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        // MS1 only: scans 19, 21, 22 → indices 1, 3, 4
        ArrayList<Integer> ms1 = parser.getSpecIndexList(1, 1);
        Assert.assertEquals(3, ms1.size());

        // MS2 only: scan 20 → index 2
        ArrayList<Integer> ms2 = parser.getSpecIndexList(2, 2);
        Assert.assertEquals(1, ms2.size());
        Assert.assertEquals(Integer.valueOf(2), ms2.get(0));
    }

    @Test
    public void testIndexMetadata() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());

        // Check scan=19 (index 1)
        StaxMzMLParser.SpectrumIndex si1 = parser.getSpectrumIndex(1);
        Assert.assertNotNull(si1);
        Assert.assertEquals(1, si1.specIndex);
        Assert.assertEquals("scan=19", si1.id);
        Assert.assertEquals(19, si1.scanNum);
        Assert.assertEquals(1, si1.msLevel);
        Assert.assertEquals(15, si1.defaultArrayLength);

        // Check scan=20 (index 2) - MS2 with precursor
        StaxMzMLParser.SpectrumIndex si2 = parser.getSpectrumIndex(2);
        Assert.assertNotNull(si2);
        Assert.assertEquals(2, si2.specIndex);
        Assert.assertEquals("scan=20", si2.id);
        Assert.assertEquals(20, si2.scanNum);
        Assert.assertEquals(2, si2.msLevel);
        Assert.assertEquals(10, si2.defaultArrayLength);
    }

    @Test
    public void testGetID() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Assert.assertEquals("scan=19", parser.getID(1));
        Assert.assertEquals("scan=20", parser.getID(2));
        Assert.assertEquals("scan=21", parser.getID(3));
        Assert.assertEquals("scan=22", parser.getID(4));
        Assert.assertNull(parser.getID(99));
    }

    @Test
    public void testMS1SpectrumParsing() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Spectrum spec = parser.getSpectrumBySpecIndex(1);
        Assert.assertNotNull(spec);

        Assert.assertEquals("scan=19", spec.getID());
        Assert.assertEquals(1, spec.getSpecIndex());
        Assert.assertEquals(19, spec.getScanNum());
        Assert.assertEquals(1, spec.getMSLevel());
        Assert.assertTrue(spec.isCentroided());
        Assert.assertEquals(Spectrum.Polarity.POSITIVE, spec.getScanPolarity());

        // 15 peaks, 64-bit uncompressed
        Assert.assertEquals(15, spec.size());

        // First peak: m/z=0.0, intensity=15.0
        Assert.assertEquals(0.0f, spec.get(0).getMz(), 0.01f);
        Assert.assertEquals(15.0f, spec.get(0).getIntensity(), 0.01f);

        // RT = 5.89 minutes
        Assert.assertTrue(spec.getRt() > 5.8f && spec.getRt() < 6.0f);
    }

    @Test
    public void testMS2SpectrumParsing() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Spectrum spec = parser.getSpectrumBySpecIndex(2);
        Assert.assertNotNull(spec);

        Assert.assertEquals("scan=20", spec.getID());
        Assert.assertEquals(2, spec.getSpecIndex());
        Assert.assertEquals(2, spec.getMSLevel());

        // 10 peaks
        Assert.assertEquals(10, spec.size());

        // Precursor info
        Assert.assertNotNull(spec.getPrecursorPeak());
        Assert.assertEquals(2, spec.getCharge());
        Assert.assertTrue(spec.getPrecursorPeak().getIntensity() > 120000f);

        // Activation method: CID
        Assert.assertEquals(ActivationMethod.CID, spec.getActivationMethod());
    }

    @Test
    public void testEmptySpectrum() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Spectrum spec = parser.getSpectrumBySpecIndex(3);
        Assert.assertNotNull(spec);
        Assert.assertEquals("scan=21", spec.getID());
        Assert.assertEquals(1, spec.getMSLevel());
        // Empty spectrum should have 0 peaks
        Assert.assertEquals(0, spec.size());
    }

    @Test
    public void testRetentionTimeSeconds() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        // scan=22 has RT in seconds (UO:0000010)
        Spectrum spec = parser.getSpectrumBySpecIndex(4);
        Assert.assertNotNull(spec);
        Assert.assertEquals(42.05f, spec.getRt(), 0.01f);
        Assert.assertTrue(spec.getRtIsSeconds());
    }

    @Test
    public void testGetSpectrumById() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Spectrum spec = parser.getSpectrumById("scan=20");
        Assert.assertNotNull(spec);
        Assert.assertEquals(2, spec.getMSLevel());
        Assert.assertEquals(10, spec.size());
    }

    @Test
    public void testCacheHit() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        // First access
        Spectrum spec1 = parser.getSpectrumBySpecIndex(2);
        // Second access should hit cache (same object)
        Spectrum spec2 = parser.getSpectrumBySpecIndex(2);
        Assert.assertSame("Cache should return same object", spec1, spec2);
    }

    @Test
    public void testIteratorWithMSLevelFilter() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Iterator<Spectrum> itr = parser.iterator(2, 2);

        int count = 0;
        while (itr.hasNext()) {
            Spectrum spec = itr.next();
            Assert.assertEquals(2, spec.getMSLevel());
            count++;
        }
        Assert.assertEquals("Should have 1 MS2 spectrum", 1, count);
    }

    @Test
    public void testIteratorAllSpectra() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Iterator<Spectrum> itr = parser.iterator(1, Integer.MAX_VALUE);

        int count = 0;
        while (itr.hasNext()) {
            itr.next();
            count++;
        }
        Assert.assertEquals("Should have 4 spectra total", 4, count);
    }

    @Test
    public void testSpectraAccessorIntegration() throws Exception {
        File mzMLFile = getMzMLFile();
        SpectraAccessor specAcc = new SpectraAccessor(mzMLFile);
        // Default MS level range is 2,2
        Iterator<Spectrum> itr = specAcc.getSpecItr();

        int count = 0;
        while (itr.hasNext()) {
            Spectrum spec = itr.next();
            Assert.assertEquals(2, spec.getMSLevel());
            count++;
        }
        Assert.assertEquals("Should have 1 MS2 spectrum via SpectraAccessor", 1, count);
    }

    @Test
    public void testSpectraAccessorRandomAccess() throws Exception {
        File mzMLFile = getMzMLFile();
        SpectraAccessor specAcc = new SpectraAccessor(mzMLFile);
        specAcc.setMSLevelRange(1, Integer.MAX_VALUE);

        SpectrumAccessorBySpecIndex specMap = specAcc.getSpecMap();
        Assert.assertNotNull(specMap);

        // Get MS2 spectrum by index
        Spectrum spec = specMap.getSpectrumBySpecIndex(2);
        Assert.assertNotNull(spec);
        Assert.assertEquals(2, spec.getMSLevel());
        Assert.assertEquals(10, spec.size());

        // Get spectrum ID
        Assert.assertEquals("scan=20", specMap.getID(2));
    }

    @Test
    public void testPeakValuesAccuracy() throws Exception {
        StaxMzMLParser parser = new StaxMzMLParser(getMzMLFile());
        Spectrum spec = parser.getSpectrumBySpecIndex(2);
        Assert.assertNotNull(spec);

        // scan=20: 10 peaks, 64-bit float, no compression
        // Expected m/z values: 0, 2, 4, 6, 8, 10, 12, 14, 16, 18
        // Expected intensity values: 20, 18, 16, 14, 12, 10, 8, 6, 4, 2
        Assert.assertEquals(10, spec.size());

        // Peaks should be sorted by m/z
        for (int i = 0; i < spec.size() - 1; i++) {
            Assert.assertTrue("Peaks should be sorted by m/z",
                    spec.get(i).getMz() <= spec.get(i + 1).getMz());
        }

        // Check first peak
        Assert.assertEquals(0.0f, spec.get(0).getMz(), 0.01f);
        Assert.assertEquals(20.0f, spec.get(0).getIntensity(), 0.01f);

        // Check last peak
        Assert.assertEquals(18.0f, spec.get(9).getMz(), 0.01f);
        Assert.assertEquals(2.0f, spec.get(9).getIntensity(), 0.01f);
    }

    @Test
    public void testBinaryDataDecoding() {
        // Test the static decodeBinaryData method directly
        // 64-bit float, no compression, 3 values: 1.0, 2.0, 3.0
        // Base64 of little-endian 64-bit doubles
        String base64 = "AAAAAAAA8D8AAAAAAAAAQAAAAAAAAAhA";
        float[] values = StaxMzMLParser.decodeBinaryData(base64, 64, false, 3);
        Assert.assertEquals(3, values.length);
        Assert.assertEquals(1.0f, values[0], 0.001f);
        Assert.assertEquals(2.0f, values[1], 0.001f);
        Assert.assertEquals(3.0f, values[2], 0.001f);
    }

    @Test
    public void testScanNumberParsing() {
        Assert.assertEquals(19, StaxMzMLParser.parseScanNumber("scan=19"));
        Assert.assertEquals(20, StaxMzMLParser.parseScanNumber("controllerType=0 controllerNumber=1 scan=20"));
        Assert.assertEquals(-1, StaxMzMLParser.parseScanNumber("no_scan_here"));
        Assert.assertEquals(-1, StaxMzMLParser.parseScanNumber(null));
    }
}
