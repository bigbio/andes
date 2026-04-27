package msgfplus;

import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.Ignore;
import org.junit.Test;

import edu.ucsd.msjava.cli.MSGFPlusOptions;
import picocli.CommandLine;
import edu.ucsd.msjava.cli.MSGFPlus;

@Ignore
public class TestCollaboration {

    @Test
    @Ignore
    public void testSujunLiIndiana()
    {
        File dir = new File("C:\\cygwin\\home\\kims336\\Data\\Sujun");

        File specFile = new File(dir.getPath()+File.separator+"scan22564.mgf");
        File dbFile = new File(dir.getPath()+File.separator+"scan22564.fasta");
        File modFile = new File(dir.getPath()+File.separator+"Mods.txt");
        String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-t", "2.5Da", "-mod", modFile.getPath()
                }; 

        MSGFPlusOptions paramManager = new MSGFPlusOptions();
        
        String msg = null; MSGFPlusOptions.commandLine(paramManager).parseArgs(argv);
        if(msg != null)
            System.out.println(msg);
        assertTrue(msg == null);
        
        assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
    }
}
