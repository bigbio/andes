package edu.ucsd.msjava.mzid;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.List;

/**
 * Regression tests locking in the fix for
 * <a href="https://github.com/MSGFPlus/msgfplus/issues/72">MSGFPlus/msgfplus#72</a>:
 * whatever the user passed via {@code -maxMissedCleavages} (or the default when
 * omitted) must be propagated to the mzIdentML output
 * {@code /SpectrumIdentificationProtocol/Enzymes/Enzyme/@missedCleavages} field.
 *
 * Prior to the #72 fix the attribute was left unset, so downstream consumers
 * that trusted the mzid (OpenMS / PeptideShaker / IDPicker) saw
 * "missedCleavages=0" instead of whatever the search actually used.
 */
public class AnalysisProtocolCollectionGenTest {

    private ParamManager buildParamManagerWith(int maxMissedCleavages) throws URISyntaxException {
        ParamManager manager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE,
                "java -Xmx3500M -jar MSGFPlus.jar");
        manager.addMSGFPlusParams();

        URI paramUri = SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI();
        manager.getParameter("conf").parse(new File(paramUri).getAbsolutePath());

        URI specUri = SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI();
        manager.getParameter("s").parse(new File(specUri).getAbsolutePath());

        URI dbUri = SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI();
        manager.getParameter("d").parse(new File(dbUri).getAbsolutePath());

        manager.getParameter("maxMissedCleavages").parse(String.valueOf(maxMissedCleavages));
        return manager;
    }

    private int missedCleavagesInGeneratedMzid(int maxMissedCleavages) throws URISyntaxException {
        ParamManager manager = buildParamManagerWith(maxMissedCleavages);
        SearchParams params = new SearchParams();
        params.parse(manager);

        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        AnalysisProtocolCollectionGen gen = new AnalysisProtocolCollectionGen(params, aaSet);

        List<Enzyme> enzymeList = gen.getSpectrumIdentificationProtocol().getEnzymes().getEnzyme();
        Assert.assertFalse("mzid enzyme list should not be empty", enzymeList.isEmpty());
        return enzymeList.get(0).getMissedCleavages();
    }

    @Test
    public void missedCleavagesIsPropagatedToMzid() throws URISyntaxException {
        Assert.assertEquals(
                "User-supplied -maxMissedCleavages must reach mzid (regression for MSGFPlus/msgfplus#72)",
                2, missedCleavagesInGeneratedMzid(2));
    }

    @Test
    public void missedCleavagesOfZeroIsPropagated() throws URISyntaxException {
        Assert.assertEquals(
                "-maxMissedCleavages 0 must be carried to mzid, not silently replaced",
                0, missedCleavagesInGeneratedMzid(0));
    }

    @Test
    public void missedCleavagesUnlimitedIsPropagated() throws URISyntaxException {
        // -1 means "no limit" per ParamManager.MAX_MISSED_CLEAVAGES default.
        Assert.assertEquals(
                "-maxMissedCleavages -1 (unlimited) must be carried to mzid",
                -1, missedCleavagesInGeneratedMzid(-1));
    }
}
