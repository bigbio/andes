package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.cli.MSGFPlusOptions;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;

public class SearchParamsTest {

    @Test
    public void parse() throws URISyntaxException {
        MSGFPlusOptions opts = new MSGFPlusOptions();

        URI url = SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI();
        opts.configFile = new File(url);

        url = SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI();
        opts.spectrumFile = new File(url);

        url = SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI();
        opts.databaseFile = new File(url);

        SearchParams params = new SearchParams();
        String err = params.parse(opts);
        Assert.assertNull("SearchParams.parse returned: " + err, err);

        Assert.assertEquals("HighRes", opts.effectiveInstrumentType().getName());
        Assert.assertEquals("20.0 ppm", params.getLeftPrecursorMassTolerance().toString());
        Assert.assertEquals("20.0 ppm", params.getRightPrecursorMassTolerance().toString());
    }
}
