package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URISyntaxException;

import org.junit.Ignore;
import org.junit.Test;
import picocli.CommandLine;

import edu.ucsd.msjava.cli.MSGFPlus;
import edu.ucsd.msjava.cli.MSGFPlusOptions;

public class TestPercolator {

    @Test
    @Ignore
    public void testAddFeatures() throws URISyntaxException {
        File specFile = new File(TestPercolator.class.getClassLoader().getResource("iprg-2013/F13.mgf").toURI());
        File dbFile = new File(TestPercolator.class.getClassLoader().getResource("iprg-2013/Homo_sapiens_non-redundant.GRCh37.68.pep.all_FPKM-cRAP.fasta").toURI());
        String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-addFeatures", "1", "-m", "3"};

        MSGFPlusOptions opts = new MSGFPlusOptions();
        new CommandLine(opts).parseArgs(argv);

        assertTrue(MSGFPlus.runMSGFPlus(opts) == null);
    }
}
