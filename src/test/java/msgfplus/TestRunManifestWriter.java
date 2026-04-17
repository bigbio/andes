package msgfplus;

import edu.ucsd.msjava.misc.RunManifestWriter;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msdbsearch.SearchParamsTest;
import edu.ucsd.msjava.msutil.DBSearchIOFiles;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Map;

/**
 * Shape and contract tests for {@link RunManifestWriter#buildManifestMap}.
 *
 * These don't actually run a search — they construct a {@link SearchParams}
 * from the standard test fixtures and verify the manifest map contains the
 * expected keys and that the values match what was passed on the CLI.
 * End-to-end write-to-disk is exercised by the last test.
 */
public class TestRunManifestWriter {

    private SearchParams parsedSearchParams() throws URISyntaxException {
        ParamManager manager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE,
                "java -Xmx3500M -jar MSGFPlus.jar");
        manager.addMSGFPlusParams();

        URI paramUri = SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI();
        manager.getParameter("conf").parse(new File(paramUri).getAbsolutePath());

        URI specUri = SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI();
        manager.getParameter("s").parse(new File(specUri).getAbsolutePath());

        URI dbUri = SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI();
        manager.getParameter("d").parse(new File(dbUri).getAbsolutePath());

        manager.getParameter("maxMissedCleavages").parse("2");

        SearchParams params = new SearchParams();
        String err = params.parse(manager);
        Assert.assertNull("SearchParams.parse should succeed: " + err, err);
        return params;
    }

    private DBSearchIOFiles firstIo(SearchParams params) {
        return params.getDBSearchIOList().get(0);
    }

    @Test
    public void manifestMapHasRequiredIdentityFields() throws URISyntaxException {
        SearchParams params = parsedSearchParams();
        DBSearchIOFiles io = firstIo(params);

        Map<String, Object> m = RunManifestWriter.buildManifestMap(
                io, params, "Release (v-test)", new String[]{"-s", "x.mgf", "-d", "y.fasta"});

        Assert.assertEquals("Release (v-test)", m.get("msgfplus_version"));
        Assert.assertNotNull("run_timestamp_utc must be set", m.get("run_timestamp_utc"));
        Assert.assertEquals(System.getProperty("java.version"), m.get("java_version"));
        Assert.assertEquals(System.getProperty("os.name"), m.get("os_name"));
        Assert.assertNotNull("max_heap_mb must be set", m.get("max_heap_mb"));
        Assert.assertTrue("available_processors must be positive",
                ((Number) m.get("available_processors")).intValue() > 0);
    }

    @Test
    public void manifestMapEchoesSearchParams() throws URISyntaxException {
        SearchParams params = parsedSearchParams();
        DBSearchIOFiles io = firstIo(params);

        Map<String, Object> m = RunManifestWriter.buildManifestMap(
                io, params, MSGFPlus.VERSION, new String[0]);

        Assert.assertEquals(2, m.get("max_missed_cleavages"));
        Assert.assertEquals(params.getMinCharge(), m.get("min_charge"));
        Assert.assertEquals(params.getMaxCharge(), m.get("max_charge"));
        Assert.assertEquals(params.getMinPeptideLength(), m.get("min_peptide_length"));
        Assert.assertEquals(params.getMaxPeptideLength(), m.get("max_peptide_length"));
        Assert.assertEquals(params.getEnzyme().getName(), m.get("enzyme"));
        Assert.assertEquals(io.getSpecFile().getAbsolutePath(), m.get("spec_file"));
        Assert.assertEquals(io.getOutputFile().getAbsolutePath(), m.get("output_file"));
        Assert.assertEquals(params.getDatabaseFile().getAbsolutePath(), m.get("fasta_file"));
    }

    @Test
    public void manifestMapPreservesCliArgs() throws URISyntaxException {
        SearchParams params = parsedSearchParams();
        DBSearchIOFiles io = firstIo(params);
        String[] argv = {"-s", "demo.mgf", "-d", "demo.fasta", "-t", "10ppm", "-e", "1"};

        Map<String, Object> m = RunManifestWriter.buildManifestMap(
                io, params, MSGFPlus.VERSION, argv);

        Object cli = m.get("cli_args");
        Assert.assertTrue("cli_args should be iterable", cli instanceof Iterable);
        int i = 0;
        for (Object token : (Iterable<?>) cli) {
            Assert.assertEquals(argv[i++], token);
        }
        Assert.assertEquals(argv.length, i);
    }

    @Test
    public void nullArgvIsToleratedAsEmptyList() throws URISyntaxException {
        SearchParams params = parsedSearchParams();
        DBSearchIOFiles io = firstIo(params);

        Map<String, Object> m = RunManifestWriter.buildManifestMap(
                io, params, MSGFPlus.VERSION, null);

        Object cli = m.get("cli_args");
        Assert.assertTrue(cli instanceof Iterable);
        Assert.assertFalse("null argv should serialise as empty list",
                ((Iterable<?>) cli).iterator().hasNext());
    }

    @Test
    public void writeProducesValidJsonSidecar() throws Exception {
        SearchParams params = parsedSearchParams();
        DBSearchIOFiles io = firstIo(params);

        // Override the DBSearchIOFiles output path so we don't write next to the
        // real test resources. Easiest way: create a fresh DBSearchIOFiles that
        // points at a temp mzid path but reuses the spec file.
        File tmpDir = java.nio.file.Files.createTempDirectory("msgfplus-manifest-test").toFile();
        File tmpOut = new File(tmpDir, "sidecar.mzid");
        DBSearchIOFiles tmpIo = new DBSearchIOFiles(io.getSpecFile(), io.getSpecFileFormat(), tmpOut);

        try {
            RunManifestWriter.write(tmpIo, params, "Release (v-test)", new String[]{"-s", "x.mgf"});

            File manifest = new File(tmpOut.getPath() + ".manifest.json");
            Assert.assertTrue("Manifest sidecar should exist at " + manifest, manifest.exists());

            String content = new String(java.nio.file.Files.readAllBytes(manifest.toPath()),
                    java.nio.charset.StandardCharsets.UTF_8);
            Assert.assertTrue("Manifest should start with '{'", content.trim().startsWith("{"));
            Assert.assertTrue("Manifest should end with '}'", content.trim().endsWith("}"));
            Assert.assertTrue("Manifest should contain msgfplus_version key",
                    content.contains("\"msgfplus_version\""));
            Assert.assertTrue("Manifest should echo the supplied version",
                    content.contains("\"Release (v-test)\""));
        } finally {
            new File(tmpOut.getPath() + ".manifest.json").delete();
            tmpOut.delete();
            tmpDir.delete();
        }
    }
}
