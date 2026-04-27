package edu.ucsd.msjava.cli;

import java.io.File;
import java.net.URISyntaxException;

/** Shared test helpers for the standard search fixture set
 *  ({@code MSGFDB_Param.txt} + {@code test.mgf} + {@code human-uniprot-contaminants.fasta}). */
public final class SearchTestFixtures {

    private SearchTestFixtures() {}

    /** Build an {@link MSGFPlusOptions} pointing at the bundled
     *  {@code MSGFDB_Param.txt} config, {@code test.mgf} spectra, and
     *  {@code human-uniprot-contaminants.fasta} database. */
    public static MSGFPlusOptions standardOpts() throws URISyntaxException {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        opts.configFile   = resource("MSGFDB_Param.txt");
        opts.spectrumFile = resource("test.mgf");
        opts.databaseFile = resource("human-uniprot-contaminants.fasta");
        return opts;
    }

    private static File resource(String name) throws URISyntaxException {
        return new File(SearchTestFixtures.class.getClassLoader().getResource(name).toURI());
    }
}
