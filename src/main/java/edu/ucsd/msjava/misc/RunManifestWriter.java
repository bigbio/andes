package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msutil.DBSearchIOFiles;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Writes a JSON run-manifest sidecar alongside each mzIdentML output.
 *
 * <p>The manifest captures the run context — MS-GF+ version, Java version and
 * heap, host OS, thread count, enzyme / instrument / activation / protocol,
 * precursor tolerance, isotope range, length / charge / mod bounds, FASTA
 * path and size, original CLI argv — so that downstream pipelines
 * (quantms, Galaxy-P, custom scripts) can reproduce or verify a search
 * without re-parsing logs.
 *
 * <p>Output path is {@code <outputMzid>.manifest.json}. The JSON is hand-rolled
 * with a stable key order; no new dependencies are pulled in.
 *
 * <p>Failures to write are logged as warnings via {@link MSGFLogger} and never
 * abort the search — the manifest is advisory metadata, not search output.
 */
public final class RunManifestWriter {

    private RunManifestWriter() {}

    /**
     * Write a manifest for the given IO pair. Caller is responsible for
     * invoking this after the mzid has been written successfully.
     *
     * @param io        spectrum/output pair from {@link SearchParams#getDBSearchIOList()}
     * @param params    parsed search parameters
     * @param version   MS-GF+ version string (e.g. {@code "v2024.07.27"})
     * @param argv      original CLI argv (used verbatim under {@code "cli_args"})
     */
    public static void write(DBSearchIOFiles io, SearchParams params, String version, String[] argv) {
        File outputFile = io.getOutputFile();
        File manifestFile = new File(outputFile.getPath() + ".manifest.json");
        try {
            Map<String, Object> m = buildManifestMap(io, params, version, argv);
            try (BufferedWriter w = Files.newBufferedWriter(manifestFile.toPath(), StandardCharsets.UTF_8)) {
                writeJson(w, m, 0);
                w.write("\n");
            }
            MSGFLogger.debug("Run manifest written to " + manifestFile.getPath());
        } catch (IOException | RuntimeException e) {
            MSGFLogger.warn("Could not write run manifest to %s: %s", manifestFile.getPath(), e.getMessage());
        }
    }

    /** Testing and inspection hook. Builds the manifest map without writing to disk. */
    public static Map<String, Object> buildManifestMap(DBSearchIOFiles io, SearchParams params, String version, String[] argv) {
        Map<String, Object> m = new LinkedHashMap<String, Object>();
        m.put("msgfplus_version", version);
        m.put("run_timestamp_utc", Instant.now().toString());

        m.put("java_version", System.getProperty("java.version"));
        m.put("java_vendor", System.getProperty("java.vendor"));
        m.put("os_name", System.getProperty("os.name"));
        m.put("os_version", System.getProperty("os.version"));
        m.put("os_arch", System.getProperty("os.arch"));

        Runtime rt = Runtime.getRuntime();
        m.put("max_heap_mb", rt.maxMemory() / (1024L * 1024L));
        m.put("available_processors", rt.availableProcessors());
        m.put("requested_threads", params.getNumThreads());
        m.put("num_tasks", params.getNumTasks());
        m.put("min_spectra_per_thread", params.getMinSpectraPerThread());

        File specFile = io.getSpecFile();
        m.put("spec_file", specFile.getAbsolutePath());
        m.put("spec_file_size_bytes", specFile.length());
        m.put("spec_file_format", io.getSpecFileFormat() == null ? null : io.getSpecFileFormat().toString());

        File fastaFile = params.getDatabaseFile();
        if (fastaFile != null) {
            m.put("fasta_file", fastaFile.getAbsolutePath());
            m.put("fasta_file_size_bytes", fastaFile.length());
        }

        File outputFile = io.getOutputFile();
        m.put("output_file", outputFile.getAbsolutePath());

        m.put("enzyme", params.getEnzyme() == null ? null : params.getEnzyme().getName());
        m.put("activation_method", params.getActivationMethod() == null ? null : params.getActivationMethod().getName());
        m.put("instrument", params.getInstType() == null ? null : params.getInstType().getName());
        m.put("protocol", params.getProtocol() == null ? null : params.getProtocol().getName());

        m.put("precursor_tol_left", params.getLeftPrecursorMassTolerance() == null ? null : params.getLeftPrecursorMassTolerance().toString());
        m.put("precursor_tol_right", params.getRightPrecursorMassTolerance() == null ? null : params.getRightPrecursorMassTolerance().toString());
        m.put("isotope_error_min", params.getMinIsotopeError());
        m.put("isotope_error_max", params.getMaxIsotopeError());

        m.put("num_tolerable_termini", params.getNumTolerableTermini());
        m.put("min_peptide_length", params.getMinPeptideLength());
        m.put("max_peptide_length", params.getMaxPeptideLength());
        m.put("min_charge", params.getMinCharge());
        m.put("max_charge", params.getMaxCharge());
        m.put("max_missed_cleavages", params.getMaxMissedCleavages());
        m.put("num_matches_per_spec", params.getNumMatchesPerSpec());
        m.put("min_ms_level", params.getMinMSLevel());
        m.put("max_ms_level", params.getMaxMSLevel());

        m.put("cli_args", argv == null ? new ArrayList<String>() : java.util.Arrays.asList(argv));
        return m;
    }

    // --- tiny hand-rolled JSON writer -----------------------------------
    // Keeps the jar dep-free. Supports String, Number, Boolean, null,
    // List/Iterable of the same, and Map<String, ?> via nested emit.

    private static void writeJson(BufferedWriter w, Object value, int indent) throws IOException {
        if (value == null) {
            w.write("null");
            return;
        }
        if (value instanceof Map) {
            @SuppressWarnings("unchecked")
            Map<String, Object> map = (Map<String, Object>) value;
            w.write("{");
            boolean first = true;
            for (Map.Entry<String, Object> e : map.entrySet()) {
                if (!first) w.write(",");
                first = false;
                w.write("\n");
                indent(w, indent + 1);
                w.write(jsonString(e.getKey()));
                w.write(": ");
                writeJson(w, e.getValue(), indent + 1);
            }
            if (!first) {
                w.write("\n");
                indent(w, indent);
            }
            w.write("}");
            return;
        }
        if (value instanceof Iterable) {
            w.write("[");
            boolean first = true;
            for (Object item : (Iterable<?>) value) {
                if (!first) w.write(", ");
                first = false;
                writeJson(w, item, indent + 1);
            }
            w.write("]");
            return;
        }
        if (value instanceof Number || value instanceof Boolean) {
            w.write(value.toString());
            return;
        }
        w.write(jsonString(value.toString()));
    }

    private static void indent(BufferedWriter w, int level) throws IOException {
        for (int i = 0; i < level; i++) w.write("  ");
    }

    private static String jsonString(String s) {
        StringBuilder sb = new StringBuilder(s.length() + 2);
        sb.append('"');
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '"':  sb.append("\\\""); break;
                case '\\': sb.append("\\\\"); break;
                case '\n': sb.append("\\n"); break;
                case '\r': sb.append("\\r"); break;
                case '\t': sb.append("\\t"); break;
                case '\b': sb.append("\\b"); break;
                case '\f': sb.append("\\f"); break;
                default:
                    if (c < 0x20) {
                        sb.append(String.format("\\u%04x", (int) c));
                    } else {
                        sb.append(c);
                    }
            }
        }
        sb.append('"');
        return sb.toString();
    }
}
