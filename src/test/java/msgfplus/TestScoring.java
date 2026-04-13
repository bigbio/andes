package msgfplus;

import java.io.File;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.nio.file.Paths;
import edu.ucsd.msjava.msutil.*;
import org.junit.Ignore;
import org.junit.Test;

import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.ScoringParameterGeneratorWithErrors;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.mzml.StaxMzMLParser;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.ScoringParamGen;


public class TestScoring {

    @Test
    public void testReadingParamFile() throws URISyntaxException {
        String paramFile = new File(TestScoring.class.getClassLoader().getResource("HCD_HighRes_Tryp_TMT.param").toURI()).getAbsolutePath();
        NewRankScorer scorer = new NewRankScorer(paramFile);
    }
        
        @Test
        public void testWritingParamAsPlainText() throws URISyntaxException {
            String paramFile = new File(TestScoring.class.getClassLoader().getResource("HCD_QExactive_Tryp.param").toURI()).getAbsolutePath();
            NewRankScorer scorer = new NewRankScorer(paramFile);
            
            Path paramPath = Paths.get(paramFile);
            String fileName = paramPath.getFileName().toString();
            String path = paramPath.getParent().toString();
            
            String output = String.format("%1s\\%2s.txt", path, fileName);
            scorer.writeParametersPlainText(new File(output));
        }
    
    @Test
    @Ignore
    public void testScoringParamGen()
    {
        File resultPath = new File("C:\\cygwin\\home\\kims336\\Data\\Scoring");
        File specPath = new File("C:\\cygwin\\home\\kims336\\Data\\Scoring");

        String[] argv = {"-i", resultPath.getPath(), "-d", specPath.getPath(), "-m", "2", "-inst", "3", "-e", "0" 
                ,"-protocol", "5"
                };
        
        ParamManager paramManager = new ParamManager("ScoringParamGen", "Test", "Test",
                "java -Xmx2000M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen");
            
        StaxMzMLParser.turnOffLogs();
        paramManager.addScoringParamGenParams();
        paramManager.parseParams(argv);
        ScoringParamGen.runScoringParamGen(paramManager);
        System.out.println("Done");        
    }        
    
    @Test
    @Ignore
    public void testScoringParamGenFromMgf()
    {
        ActivationMethod actMethod = ActivationMethod.HCD;
        InstrumentType instType = InstrumentType.QEXACTIVE;
        Enzyme enzyme = Enzyme.TRYPSIN;
        Protocol protocol = Protocol.STANDARD;
        File specFile = new File("D:\\Research\\Data\\TrainingMSGFPlus\\AnnotatedSpectra\\HCD_QExactive_Tryp.mgf");
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        
        SpecDataType dataType = new SpecDataType(actMethod, instType, enzyme, protocol);
        System.out.println("Processing " + dataType.toString());
        ScoringParameterGeneratorWithErrors.generateParameters(
                specFile,
                dataType,
                aaSet, 
                new File("C:\\cygwin\\home\\kims336\\Data\\Scoring"),
                false, 
                false,
                false);
        
    }
    
}
