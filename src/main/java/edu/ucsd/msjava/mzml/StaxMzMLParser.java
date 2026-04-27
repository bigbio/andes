package edu.ucsd.msjava.mzml;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.CvParamInfo;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;

import org.slf4j.LoggerFactory;
import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * StAX-based mzML parser optimized for MS-GF+ usage patterns.
 *
 * Design:
 * - Single-pass index build: scans the file once to record byte offsets and
 *   lightweight metadata (MS level, precursor m/z) for every spectrum.
 * - Random access: seeks to the byte offset and parses only the requested spectrum.
 * - Full preload cache: on first random access, all spectra are parsed and cached
 *   in memory to avoid repeated XML parsing during the database search phase.
 * - Extracts only the 11 fields MSGF+ needs (no full JAXB object model).
 */
public class StaxMzMLParser {

    /** Indexed metadata for each spectrum, built during the index pass. */
    public static class SpectrumIndex {
        public final int specIndex;       // 1-based
        public final String id;
        public final int scanNum;
        public final int msLevel;
        public final float precursorMz;
        public final long byteOffset;     // byte offset of <spectrum> element
        public final int defaultArrayLength;

        SpectrumIndex(int specIndex, String id, int scanNum, int msLevel,
                      float precursorMz, long byteOffset, int defaultArrayLength) {
            this.specIndex = specIndex;
            this.id = id;
            this.scanNum = scanNum;
            this.msLevel = msLevel;
            this.precursorMz = precursorMz;
            this.byteOffset = byteOffset;
            this.defaultArrayLength = defaultArrayLength;
        }
    }

    private final File specFile;
    private final List<SpectrumIndex> indexList;              // ordered by specIndex
    private final Map<Integer, SpectrumIndex> indexBySpecIdx; // specIndex -> index entry
    private final Map<String, SpectrumIndex> indexById;       // id -> index entry

    // Referenceable param groups: group ID -> list of [accession, name, value, unitAccession, unitName]
    private final Map<String, List<String[]>> refParamGroups;

    /** MS-level filter: spectra outside this range are never decoded or cached. */
    private final int minMSLevel;
    private final int maxMSLevel;

    /** Synchronized cache of in-filter spectra. Returns defensive copies on read
     *  so pre-pass mutations cannot leak to the main pass. */
    private final Map<Integer, Spectrum> cache;
    private volatile boolean allLoaded = false;

    // Reusable XMLInputFactory (thread-safe for creation)
    private static final XMLInputFactory XML_INPUT_FACTORY;
    static {
        XML_INPUT_FACTORY = XMLInputFactory.newInstance();
        XML_INPUT_FACTORY.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, false);
        XML_INPUT_FACTORY.setProperty(XMLInputFactory.IS_VALIDATING, false);
        XML_INPUT_FACTORY.setProperty(XMLInputFactory.SUPPORT_DTD, false);
        XML_INPUT_FACTORY.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, false);
    }

    /**
     * Construct a parser for the given mzML file with no MS-level filter.
     * Prefer {@link #StaxMzMLParser(File, int, int)} so MS1 spectra can be
     * skipped during the binary-decode preload.
     */
    public StaxMzMLParser(File specFile) throws IOException, XMLStreamException {
        this(specFile, 1, Integer.MAX_VALUE);
    }

    /**
     * Construct a parser for the given mzML file, decoding/caching only spectra
     * with MS level inside {@code [minMSLevel, maxMSLevel]}. Immediately builds
     * the spectrum index (single sequential pass; no peak decode).
     */
    public StaxMzMLParser(File specFile, int minMSLevel, int maxMSLevel) throws IOException, XMLStreamException {
        this.specFile = specFile;
        this.minMSLevel = minMSLevel;
        this.maxMSLevel = maxMSLevel;
        this.indexList = new ArrayList<>();
        this.indexBySpecIdx = new HashMap<>();
        this.indexById = new HashMap<>();
        this.refParamGroups = new HashMap<>();
        this.cache = Collections.synchronizedMap(new HashMap<>());
        buildIndex();
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    public int getSpectrumCount() {
        return indexList.size();
    }

    public ArrayList<Integer> getSpecIndexList() {
        ArrayList<Integer> list = new ArrayList<>(indexList.size());
        for (SpectrumIndex si : indexList) list.add(si.specIndex);
        return list;
    }

    public ArrayList<Integer> getSpecIndexList(int minMSLevel, int maxMSLevel) {
        ArrayList<Integer> list = new ArrayList<>();
        for (SpectrumIndex si : indexList) {
            if (si.msLevel >= minMSLevel && si.msLevel <= maxMSLevel)
                list.add(si.specIndex);
        }
        return list;
    }

    public SpectrumIndex getSpectrumIndex(int specIndex) {
        return indexBySpecIdx.get(specIndex);
    }

    public String getID(int specIndex) {
        SpectrumIndex si = indexBySpecIdx.get(specIndex);
        return si != null ? si.id : null;
    }

    public Float getPrecursorMz(int specIndex) {
        SpectrumIndex si = indexBySpecIdx.get(specIndex);
        if (si == null) return null;
        return si.precursorMz > 0 ? si.precursorMz : null;
    }

    /**
     * Parse and return the full spectrum (with peaks) for the given 1-based index.
     * On first cache miss, performs a bulk preload of all in-filter spectra; every
     * subsequent call returns a defensive copy from the cache. Returns {@code null}
     * for unknown indices and for spectra outside the configured MS-level filter.
     */
    public Spectrum getSpectrumBySpecIndex(int specIndex) {
        SpectrumIndex si = indexBySpecIdx.get(specIndex);
        if (si == null) return null;
        if (si.msLevel < minMSLevel || si.msLevel > maxMSLevel) return null;

        if (!allLoaded && !cache.containsKey(specIndex)) {
            try {
                preloadAllSpectra();
            } catch (Exception e) {
                throw new RuntimeException("Failed to preload spectra while retrieving spectrum index " + specIndex, e);
            }
        }
        return cloneSpectrum(cache.get(specIndex));
    }

    /**
     * Walk the file once and cache every in-filter spectrum. Out-of-filter
     * spectra are skipped without binary decode — no {@link Spectrum} or
     * {@code Peak} objects allocated.
     */
    private synchronized void preloadAllSpectra() throws IOException, XMLStreamException {
        if (allLoaded) return;
        long startTime = System.currentTimeMillis();
        int loaded = 0, skipped = 0;
        try (InputStream is = new BufferedInputStream(new FileInputStream(specFile), 256 * 1024)) {
            XMLStreamReader reader = XML_INPUT_FACTORY.createXMLStreamReader(is);
            try {
                while (reader.hasNext()) {
                    int event = reader.next();
                    if (event == XMLStreamConstants.START_ELEMENT && "spectrum".equals(reader.getLocalName())) {
                        // Skip out-of-filter spectra without binary decode by consulting
                        // the pre-built index via the spectrum's id attribute.
                        String id = reader.getAttributeValue(null, "id");
                        SpectrumIndex si = id != null ? indexById.get(id) : null;
                        if (si != null && (si.msLevel < minMSLevel || si.msLevel > maxMSLevel)) {
                            skipElement(reader, "spectrum");
                            skipped++;
                            continue;
                        }
                        Spectrum spec = parseOneSpectrum(reader);
                        if (spec != null) {
                            int ms = spec.getMSLevel();
                            if (ms < minMSLevel || ms > maxMSLevel) {
                                // Index lookup missed (malformed mzML id mismatch); drop post-parse.
                                skipped++;
                                continue;
                            }
                            cache.put(spec.getSpecIndex(), spec);
                            loaded++;
                        }
                    }
                }
            } finally {
                reader.close();
            }
        } catch (XMLStreamException e) {
            throw annotate(e, "preload");
        }
        allLoaded = true;
        long elapsed = System.currentTimeMillis() - startTime;
        System.out.println("StAX mzML preload: " + loaded + " spectra loaded (" + skipped + " filtered out by MS level) in " + elapsed + " ms");
    }

    /**
     * Defensive copy of a cached {@link Spectrum}. Mirrors the field set
     * populated by {@link #parseOneSpectrum}; keep the two in lock-step.
     */
    private static Spectrum cloneSpectrum(Spectrum src) {
        if (src == null) return null;
        Spectrum dst = new Spectrum();
        dst.setID(src.getID());
        dst.setSpecIndex(src.getSpecIndex());
        if (src.getScanNum() > 0) dst.setScanNum(src.getScanNum());
        dst.setMsLevel(src.getMSLevel());
        dst.setIsCentroided(src.isCentroided());
        if (src.getScanPolarity() != null) dst.setScanPolarity(src.getScanPolarity());
        dst.setRt(src.getRt());
        dst.setRtIsSeconds(src.getRtIsSeconds());
        dst.setIsolationWindowTargetMz(src.getIsolationWindowTargetMz());
        if (src.getPrecursorPeak() != null) {
            Peak p = src.getPrecursorPeak();
            dst.setPrecursor(new Peak(p.getMz(), p.getIntensity(), p.getCharge()));
        }
        if (src.getActivationMethod() != null) dst.setActivationMethod(src.getActivationMethod());
        if (src.getAddlCvParams() != null) {
            for (CvParamInfo cv : src.getAddlCvParams()) dst.addAddlCvParam(cv);
        }
        for (Peak p : src) {
            dst.add(new Peak(p.getMz(), p.getIntensity(), p.getCharge()));
        }
        return dst;
    }

    /**
     * Rethrow an {@link XMLStreamException} with a context-rich message. If
     * the underlying error looks like a BOM or XML-prolog / encoding issue
     * (the most common cause of "ParseError in XML prolog" on Windows),
     * suggest the concrete fix.
     *
     * @param e      the original Stax exception; wrapped as cause
     * @param phase  short tag identifying the parse phase ("index", "preload")
     */
    private XMLStreamException annotate(XMLStreamException e, String phase) {
        String msg = e.getMessage() == null ? "" : e.getMessage();
        StringBuilder sb = new StringBuilder();
        sb.append("Could not parse mzML file '").append(specFile.getAbsolutePath()).append("' during ").append(phase).append(".");
        if (looksLikeBomOrPrologIssue(msg)) {
            sb.append(" This usually means the file has a byte-order mark (BOM) or an encoding mismatch in the XML prolog. Verify that the file starts with `<?xml version=\"1.0\" encoding=\"UTF-8\"?>` with no leading whitespace or BOM (on Linux/macOS: `head -c 3 \"")
                    .append(specFile.getName()).append("\" | xxd`; a BOM shows as `ef bb bf`). Re-converting the raw file with ThermoRawFileParser or MSConvert usually resolves it. See docs/troubleshooting.md for details.");
        }
        sb.append(" Underlying parser error: ").append(msg);
        // Note: XMLStreamException(msg, location, nested) stores the cause as a
        // "nested exception" but does NOT invoke Throwable.initCause, so
        // getCause() returns null. Call initCause() explicitly so standard
        // Java chaining (printStackTrace, causal frames) works.
        XMLStreamException wrapped = new XMLStreamException(sb.toString(), e.getLocation());
        wrapped.initCause(e);
        return wrapped;
    }

    private static boolean looksLikeBomOrPrologIssue(String msg) {
        if (msg == null) return false;
        String m = msg.toLowerCase(java.util.Locale.ROOT);
        return m.contains("prolog")
                || m.contains("bom")
                || m.contains("byte order mark")
                || m.contains("encoding")
                || m.contains("invalid character")
                || m.contains("content is not allowed");
    }

    public Spectrum getSpectrumById(String specId) {
        SpectrumIndex si = indexById.get(specId);
        if (si == null) return null;
        return getSpectrumBySpecIndex(si.specIndex);
    }

    /**
     * Returns an iterator that streams spectra sequentially (no random seeks).
     * More efficient than random access when all spectra are needed.
     * Applies MS level filtering.
     */
    public Iterator<Spectrum> iterator(int minMSLevel, int maxMSLevel) {
        return new StaxSequentialIterator(minMSLevel, maxMSLevel);
    }

    public List<SpectrumIndex> getIndexList() {
        return Collections.unmodifiableList(indexList);
    }

    /**
     * Detect the spectrum ID format CV param by scanning file header.
     * Returns a 2-element array [accession, name] or null if not found.
     */
    public String[] detectSpectrumIDFormat() {
        try (InputStream is = new BufferedInputStream(new FileInputStream(specFile), 64 * 1024)) {
            XMLStreamReader reader = XML_INPUT_FACTORY.createXMLStreamReader(is);
            try {
                while (reader.hasNext()) {
                    int event = reader.next();
                    if (event == XMLStreamConstants.START_ELEMENT) {
                        String eName = reader.getLocalName();
                        if ("spectrumList".equals(eName) || "run".equals(eName))
                            break; // past file description, stop
                        if ("cvParam".equals(eName)) {
                            String acc = reader.getAttributeValue(null, "accession");
                            if (acc != null && isSpectrumIDFormatAccession(acc)) {
                                String cvName = reader.getAttributeValue(null, "name");
                                return new String[]{acc, cvName != null ? cvName : "nativeID format"};
                            }
                        }
                    }
                }
            } finally {
                reader.close();
            }
        } catch (Exception e) {
            // fall through
        }
        return null;
    }

    // -----------------------------------------------------------------------
    // Index building (single sequential pass)
    // -----------------------------------------------------------------------

    private void buildIndex() throws IOException, XMLStreamException {
        try (CountingInputStream cis = new CountingInputStream(
                new BufferedInputStream(new FileInputStream(specFile), 256 * 1024))) {
            XMLStreamReader reader = XML_INPUT_FACTORY.createXMLStreamReader(cis);
            try {
                buildIndexFromReader(reader, cis);
            } finally {
                reader.close();
            }
        } catch (XMLStreamException e) {
            throw annotate(e, "index");
        }
    }

    private void buildIndexFromReader(XMLStreamReader reader, CountingInputStream cis)
            throws XMLStreamException {
        boolean inSpectrum = false;
        boolean inPrecursor = false;
        boolean inSelectedIon = false;
        boolean inScan = false;
        boolean inRefParamGroup = false;
        String curRefGroupId = null;

        // Current spectrum being indexed
        int curIndex = -1;
        String curId = null;
        int curScanNum = -1;
        int curMsLevel = 0;
        float curPrecursorMz = -1;
        long curByteOffset = 0;
        int curArrayLength = 0;

        while (reader.hasNext()) {
            int event = reader.next();

            if (event == XMLStreamConstants.START_ELEMENT) {
                String name = reader.getLocalName();

                if ("referenceableParamGroup".equals(name)) {
                    inRefParamGroup = true;
                    curRefGroupId = reader.getAttributeValue(null, "id");
                    if (curRefGroupId != null)
                        refParamGroups.put(curRefGroupId, new ArrayList<>());
                } else if (inRefParamGroup && "cvParam".equals(name)) {
                    if (curRefGroupId != null) {
                        String acc = reader.getAttributeValue(null, "accession");
                        String cvName = reader.getAttributeValue(null, "name");
                        String val = reader.getAttributeValue(null, "value");
                        String unitAcc = reader.getAttributeValue(null, "unitAccession");
                        String unitName = reader.getAttributeValue(null, "unitName");
                        refParamGroups.get(curRefGroupId).add(
                                new String[]{acc, cvName, val, unitAcc, unitName});
                    }
                } else if ("spectrum".equals(name)) {
                    inSpectrum = true;
                    inRefParamGroup = false;
                    curByteOffset = cis.getBytesRead();
                    curId = reader.getAttributeValue(null, "id");
                    String indexStr = reader.getAttributeValue(null, "index");
                    curIndex = indexStr != null ? Integer.parseInt(indexStr) + 1 : indexList.size() + 1;
                    String arrLen = reader.getAttributeValue(null, "defaultArrayLength");
                    curArrayLength = arrLen != null ? Integer.parseInt(arrLen) : 0;
                    curScanNum = parseScanNumber(curId);
                    curMsLevel = 0;
                    curPrecursorMz = -1;
                } else if (inSpectrum && "referenceableParamGroupRef".equals(name)) {
                    // Resolve referenced param group during indexing for MS level
                    String ref = reader.getAttributeValue(null, "ref");
                    if (ref != null) {
                        List<String[]> params = refParamGroups.get(ref);
                        if (params != null) {
                            for (String[] p : params) {
                                if ("MS:1000511".equals(p[0]) && p[2] != null)
                                    curMsLevel = Integer.parseInt(p[2]);
                            }
                        }
                    }
                } else if (inSpectrum && "cvParam".equals(name)) {
                    String acc = reader.getAttributeValue(null, "accession");
                    if (acc != null) {
                        if ("MS:1000511".equals(acc)) {
                            String val = reader.getAttributeValue(null, "value");
                            curMsLevel = val != null ? Integer.parseInt(val) : 0;
                        } else if (inSelectedIon && "MS:1000744".equals(acc)) {
                            String val = reader.getAttributeValue(null, "value");
                            if (val != null) curPrecursorMz = Float.parseFloat(val);
                        } else if (inScan && "MS:1000016".equals(acc)) {
                            // retention time - skip during indexing, parse during full parse
                        }
                    }
                } else if (inSpectrum && "precursorList".equals(name)) {
                    inPrecursor = true;
                } else if (inPrecursor && "selectedIon".equals(name)) {
                    inSelectedIon = true;
                } else if (inSpectrum && "scan".equals(name)) {
                    inScan = true;
                } else if ("binaryDataArrayList".equals(name)) {
                    // Skip binary data during index pass
                    skipElement(reader, "binaryDataArrayList");
                }
            } else if (event == XMLStreamConstants.END_ELEMENT) {
                String name = reader.getLocalName();
                if ("referenceableParamGroup".equals(name)) {
                    inRefParamGroup = false;
                    curRefGroupId = null;
                } else if ("spectrum".equals(name)) {
                    SpectrumIndex si = new SpectrumIndex(
                            curIndex, curId, curScanNum, curMsLevel,
                            curPrecursorMz, curByteOffset, curArrayLength);
                    indexList.add(si);
                    indexBySpecIdx.put(curIndex, si);
                    if (curId != null) indexById.put(curId, si);
                    inSpectrum = false;
                    inPrecursor = false;
                    inSelectedIon = false;
                    inScan = false;
                } else if ("selectedIon".equals(name)) {
                    inSelectedIon = false;
                } else if ("precursorList".equals(name)) {
                    inPrecursor = false;
                } else if ("scan".equals(name)) {
                    inScan = false;
                }
            }
        }
    }

    private void skipElement(XMLStreamReader reader, String elementName) throws XMLStreamException {
        int depth = 1;
        while (reader.hasNext() && depth > 0) {
            int event = reader.next();
            if (event == XMLStreamConstants.START_ELEMENT) depth++;
            else if (event == XMLStreamConstants.END_ELEMENT) depth--;
        }
    }

    // -----------------------------------------------------------------------
    // Full spectrum parsing (random access)
    // -----------------------------------------------------------------------

    /**
     * Parse a single &lt;spectrum&gt; element. Reader is positioned just after
     * the START_ELEMENT of &lt;spectrum&gt;.
     */
    Spectrum parseOneSpectrum(XMLStreamReader reader) throws XMLStreamException {
        Spectrum spec = new Spectrum();

        // Attributes from <spectrum> element
        String id = reader.getAttributeValue(null, "id");
        String indexStr = reader.getAttributeValue(null, "index");
        String arrLenStr = reader.getAttributeValue(null, "defaultArrayLength");

        spec.setID(id);
        int specIndex = indexStr != null ? Integer.parseInt(indexStr) + 1 : 0;
        spec.setSpecIndex(specIndex);

        int scanNum = parseScanNumber(id);
        if (scanNum > 0) spec.setScanNum(scanNum);

        int defaultArrayLength = arrLenStr != null ? Integer.parseInt(arrLenStr) : 0;

        // Parse content
        boolean inScan = false;
        boolean inPrecursor = false;
        boolean inSelectedIon = false;
        boolean inActivation = false;
        boolean inIsolationWindow = false;
        boolean inBinaryDataArray = false;
        boolean inBinary = false;

        int msLevel = 0;
        boolean isCentroided = false;
        Spectrum.Polarity polarity = Spectrum.Polarity.POSITIVE;
        float scanStartTime = -1;
        boolean scanStartTimeIsSeconds = true;
        float precursorMz = -1;
        int precursorCharge = 0;
        float precursorIntensity = 0;
        Float isolationWindowTargetMz = null;
        ActivationMethod activationMethod = null;
        boolean isETD = false;
        float thermoMonoMz = -1;

        // Binary data array state
        int binaryArrayCount = 0;
        int precision = 32;           // bits (32 or 64)
        boolean compressed = false;   // zlib
        boolean isMzArray = false;
        boolean isIntensityArray = false;
        StringBuilder binaryText = null;

        float[] mzValues = null;
        float[] intensityValues = null;

        int depth = 1; // inside <spectrum>

        while (reader.hasNext() && depth > 0) {
            int event = reader.next();

            if (event == XMLStreamConstants.START_ELEMENT) {
                depth++;
                String name = reader.getLocalName();

                if ("cvParam".equals(name)) {
                    String acc = reader.getAttributeValue(null, "accession");
                    String val = reader.getAttributeValue(null, "value");

                    if (acc == null) continue;

                    // Spectrum-level CV params
                    if (!inScan && !inPrecursor && !inBinaryDataArray) {
                        switch (acc) {
                            case "MS:1000511": msLevel = parseInt(val, 0); break;
                            case "MS:1000127": isCentroided = true; break;
                            case "MS:1000128": isCentroided = false; break;
                            case "MS:1000129": polarity = Spectrum.Polarity.NEGATIVE; break;
                            case "MS:1000130": polarity = Spectrum.Polarity.POSITIVE; break;
                        }
                    }
                    // Scan-level CV params
                    else if (inScan && !inPrecursor) {
                        if ("MS:1000016".equals(acc)) {
                            scanStartTime = parseFloat(val, -1);
                            String unitAcc = reader.getAttributeValue(null, "unitAccession");
                            if ("UO:0000031".equals(unitAcc)) scanStartTimeIsSeconds = false;
                            else if ("UO:0000010".equals(unitAcc)) scanStartTimeIsSeconds = true;
                        }
                        // Ion mobility params
                        else if ("MS:1001581".equals(acc) || "MS:1002476".equals(acc) || "MS:1002815".equals(acc)) {
                            String cvName = reader.getAttributeValue(null, "name");
                            String unitAcc = reader.getAttributeValue(null, "unitAccession");
                            String unitName = reader.getAttributeValue(null, "unitName");
                            CvParamInfo info = (unitAcc != null && !unitAcc.isEmpty())
                                    ? new CvParamInfo(acc, cvName, val, unitAcc, unitName)
                                    : new CvParamInfo(acc, cvName, val);
                            spec.addAddlCvParam(info);
                        }
                    }
                    // Isolation window CV params
                    else if (inIsolationWindow) {
                        if ("MS:1000827".equals(acc))
                            isolationWindowTargetMz = parseFloat(val, -1);
                    }
                    // Selected ion CV params
                    else if (inSelectedIon) {
                        switch (acc) {
                            case "MS:1000744": // selected ion m/z
                            case "MS:1000040": // m/z (generic, used in some older files)
                                if (precursorMz < 0.01f) precursorMz = parseFloat(val, -1);
                                break;
                            case "MS:1000041": precursorCharge = parseInt(val, 0); break;
                            case "MS:1000042": precursorIntensity = parseFloat(val, 0); break;
                        }
                    }
                    // Activation CV params
                    else if (inActivation) {
                        ActivationMethod am = ActivationMethod.getByCV(acc);
                        if (am != null) {
                            if (am == ActivationMethod.ETD) {
                                isETD = true;
                            } else if (activationMethod == null) {
                                activationMethod = am;
                            }
                        }
                    }
                    // Binary data array CV params
                    else if (inBinaryDataArray && !inBinary) {
                        switch (acc) {
                            case "MS:1000523": precision = 64; break; // 64-bit float
                            case "MS:1000521": precision = 32; break; // 32-bit float
                            case "MS:1000574": compressed = true; break; // zlib
                            case "MS:1000576": compressed = false; break; // no compression
                            case "MS:1000514": isMzArray = true; break;
                            case "MS:1000515": isIntensityArray = true; break;
                        }
                    }
                }
                else if ("referenceableParamGroupRef".equals(name)) {
                    // Resolve referenced param group - apply its CV params in current context
                    String ref = reader.getAttributeValue(null, "ref");
                    if (ref != null) {
                        List<String[]> params = refParamGroups.get(ref);
                        if (params != null) {
                            for (String[] p : params) {
                                String pAcc = p[0];
                                String pVal = p[2];
                                String pUnitAcc = p[3];
                                if (pAcc == null) continue;

                                // Apply in current context (spectrum-level or scan-level)
                                if (!inScan && !inPrecursor && !inBinaryDataArray) {
                                    switch (pAcc) {
                                        case "MS:1000511": msLevel = parseInt(pVal, 0); break;
                                        case "MS:1000127": isCentroided = true; break;
                                        case "MS:1000128": isCentroided = false; break;
                                        case "MS:1000129": polarity = Spectrum.Polarity.NEGATIVE; break;
                                        case "MS:1000130": polarity = Spectrum.Polarity.POSITIVE; break;
                                    }
                                } else if (inScan && !inPrecursor) {
                                    if ("MS:1000016".equals(pAcc)) {
                                        scanStartTime = parseFloat(pVal, -1);
                                        if ("UO:0000031".equals(pUnitAcc)) scanStartTimeIsSeconds = false;
                                        else if ("UO:0000010".equals(pUnitAcc)) scanStartTimeIsSeconds = true;
                                    }
                                }
                            }
                        }
                    }
                }
                else if ("userParam".equals(name)) {
                    if (inScan) {
                        String paramName = reader.getAttributeValue(null, "name");
                        if ("[Thermo Trailer Extra]Monoisotopic M/Z:".equals(paramName)) {
                            String val = reader.getAttributeValue(null, "value");
                            thermoMonoMz = parseFloat(val, -1);
                        }
                    }
                }
                else if ("scan".equals(name)) {
                    inScan = true;
                }
                else if ("precursor".equals(name)) {
                    inPrecursor = true;
                }
                else if ("isolationWindow".equals(name)) {
                    inIsolationWindow = true;
                }
                else if ("selectedIon".equals(name)) {
                    inSelectedIon = true;
                }
                else if ("activation".equals(name)) {
                    inActivation = true;
                }
                else if ("binaryDataArray".equals(name)) {
                    inBinaryDataArray = true;
                    binaryArrayCount++;
                    precision = 32;
                    compressed = false;
                    isMzArray = false;
                    isIntensityArray = false;
                }
                else if ("binary".equals(name)) {
                    inBinary = true;
                    binaryText = new StringBuilder();
                }
            }
            else if (event == XMLStreamConstants.CHARACTERS || event == XMLStreamConstants.CDATA) {
                if (inBinary && binaryText != null) {
                    binaryText.append(reader.getText());
                }
            }
            else if (event == XMLStreamConstants.END_ELEMENT) {
                depth--;
                String name = reader.getLocalName();

                if ("binary".equals(name)) {
                    if (binaryText != null && binaryText.length() > 0) {
                        float[] values = decodeBinaryData(binaryText.toString(), precision, compressed, defaultArrayLength);
                        if (isMzArray) mzValues = values;
                        else if (isIntensityArray) intensityValues = values;
                    }
                    inBinary = false;
                    binaryText = null;
                }
                else if ("binaryDataArray".equals(name)) {
                    inBinaryDataArray = false;
                }
                else if ("scan".equals(name)) {
                    inScan = false;
                }
                else if ("selectedIon".equals(name)) {
                    inSelectedIon = false;
                }
                else if ("isolationWindow".equals(name)) {
                    inIsolationWindow = false;
                }
                else if ("activation".equals(name)) {
                    inActivation = false;
                }
                else if ("precursor".equals(name)) {
                    inPrecursor = false;
                }
                // "spectrum" end is handled by depth == 0
            }
        }

        // Assemble the Spectrum object
        spec.setMsLevel(msLevel);
        spec.setIsCentroided(isCentroided);
        spec.setScanPolarity(polarity);
        spec.setRt(scanStartTime);
        spec.setRtIsSeconds(scanStartTimeIsSeconds);
        spec.setIsolationWindowTargetMz(isolationWindowTargetMz);

        // Precursor: prefer Thermo monoisotopic M/Z if available
        if (thermoMonoMz > 0.01f) precursorMz = thermoMonoMz;
        if (precursorMz > 0) {
            spec.setPrecursor(new Peak(precursorMz, precursorIntensity, precursorCharge));
        }

        // Activation method
        if (isETD) activationMethod = ActivationMethod.ETD;
        if (activationMethod != null) spec.setActivationMethod(activationMethod);

        // Peak list
        if (mzValues != null && intensityValues != null) {
            int len = Math.min(mzValues.length, intensityValues.length);
            if (mzValues.length != intensityValues.length) {
                System.err.println("Warning: different sizes for m/z (" + mzValues.length
                        + ") and intensity (" + intensityValues.length + ") arrays for spectrum " + id);
            }
            for (int i = 0; i < len; i++) {
                spec.add(new Peak(mzValues[i], intensityValues[i], 1));
            }
        }

        // Sort peaks by m/z
        Collections.sort(spec);
        spec.determineIsCentroided();

        return spec;
    }

    // -----------------------------------------------------------------------
    // Binary data decoding
    // -----------------------------------------------------------------------

    public static float[] decodeBinaryData(String base64Text, int precision, boolean compressed, int expectedCount) {
        // Strip whitespace from base64
        byte[] encoded = stripWhitespace(base64Text);

        // Base64 decode
        byte[] decoded = java.util.Base64.getDecoder().decode(encoded);

        // Decompress if zlib
        if (compressed) {
            decoded = zlibDecompress(decoded, expectedCount * (precision / 8));
            if (decoded == null) return new float[0];
        }

        ByteBuffer buffer = ByteBuffer.wrap(decoded).order(ByteOrder.LITTLE_ENDIAN);
        int count = precision == 64 ? decoded.length / 8 : decoded.length / 4;
        float[] values = new float[count];

        if (precision == 64) {
            for (int i = 0; i < count; i++)
                values[i] = (float) buffer.getDouble();
        } else {
            for (int i = 0; i < count; i++)
                values[i] = buffer.getFloat();
        }
        return values;
    }

    private static byte[] stripWhitespace(String text) {
        // Fast path: check if there's any whitespace
        boolean hasWhitespace = false;
        for (int i = 0; i < text.length(); i++) {
            char c = text.charAt(i);
            if (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
                hasWhitespace = true;
                break;
            }
        }
        if (!hasWhitespace) return text.getBytes(java.nio.charset.StandardCharsets.ISO_8859_1);

        byte[] result = new byte[text.length()];
        int pos = 0;
        for (int i = 0; i < text.length(); i++) {
            char c = text.charAt(i);
            if (c != ' ' && c != '\n' && c != '\r' && c != '\t')
                result[pos++] = (byte) c;
        }
        return java.util.Arrays.copyOf(result, pos);
    }

    private static byte[] zlibDecompress(byte[] data, int estimatedSize) {
        Inflater inflater = new Inflater();
        try {
            inflater.setInput(data);
            byte[] result = new byte[Math.max(estimatedSize, data.length * 2)];
            int offset = 0;
            while (!inflater.finished()) {
                int remaining = result.length - offset;
                if (remaining == 0) {
                    result = java.util.Arrays.copyOf(result, result.length * 2);
                    remaining = result.length - offset;
                }
                try {
                    int count = inflater.inflate(result, offset, remaining);
                    if (count == 0 && inflater.needsInput()) break;
                    offset += count;
                } catch (DataFormatException e) {
                    System.err.println("Error decompressing binary data: " + e.getMessage());
                    return null;
                }
            }
            return java.util.Arrays.copyOf(result, offset);
        } finally {
            inflater.end();
        }
    }

    // -----------------------------------------------------------------------
    // Utility
    // -----------------------------------------------------------------------

    public static int parseScanNumber(String id) {
        if (id == null) return -1;
        // Parse "scan=NNN" from the id string
        int idx = id.lastIndexOf("scan=");
        if (idx >= 0) {
            int start = idx + 5;
            int end = start;
            while (end < id.length() && Character.isDigit(id.charAt(end))) end++;
            if (end > start) {
                try { return Integer.parseInt(id.substring(start, end)); }
                catch (NumberFormatException e) { /* fall through */ }
            }
        }
        return -1;
    }

    private static int parseInt(String s, int defaultVal) {
        if (s == null) return defaultVal;
        try { return Integer.parseInt(s); }
        catch (NumberFormatException e) { return defaultVal; }
    }

    private static float parseFloat(String s, float defaultVal) {
        if (s == null) return defaultVal;
        try { return Float.parseFloat(s); }
        catch (NumberFormatException e) { return defaultVal; }
    }

    private static boolean isSpectrumIDFormatAccession(String acc) {
        if (!acc.startsWith("MS:")) return false;
        try {
            long num = Long.parseLong(acc.substring(3));
            return (num >= 1000768 && num <= 1000777)
                    || num == 1000823 || num == 1000824 || num == 1000929
                    || num == 1001508 || num == 1001526 || num == 1001528
                    || num == 1001531 || num == 1001532 || num == 1001559
                    || num == 1001562 || num == 1002818 || num == 1001480
                    || num == 1002303 || num == 1002532 || num == 1002898;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    // -----------------------------------------------------------------------
    // CountingInputStream — tracks bytes read for offset recording
    // -----------------------------------------------------------------------

    static class CountingInputStream extends InputStream {
        private final InputStream in;
        private long bytesRead = 0;

        CountingInputStream(InputStream in) { this.in = in; }
        long getBytesRead() { return bytesRead; }

        @Override public int read() throws IOException {
            int b = in.read();
            if (b >= 0) bytesRead++;
            return b;
        }
        @Override public int read(byte[] buf, int off, int len) throws IOException {
            int n = in.read(buf, off, len);
            if (n > 0) bytesRead += n;
            return n;
        }
        @Override public void close() throws IOException { in.close(); }
    }

    // -----------------------------------------------------------------------
    // Sequential iterator (efficient single-pass)
    // -----------------------------------------------------------------------

    private class StaxSequentialIterator implements Iterator<Spectrum> {
        private final int minMSLevel, maxMSLevel;
        private XMLStreamReader reader;
        private InputStream inputStream;
        private Spectrum nextSpectrum;
        private boolean done;

        StaxSequentialIterator(int minMSLevel, int maxMSLevel) {
            this.minMSLevel = minMSLevel;
            this.maxMSLevel = maxMSLevel;
            this.done = false;
            try {
                inputStream = new BufferedInputStream(new FileInputStream(specFile), 256 * 1024);
                reader = XML_INPUT_FACTORY.createXMLStreamReader(inputStream);
                advance();
            } catch (Exception e) {
                done = true;
                System.err.println("Error creating mzML iterator: " + e.getMessage());
            }
        }

        @Override
        public boolean hasNext() {
            return nextSpectrum != null;
        }

        @Override
        public Spectrum next() {
            if (nextSpectrum == null) throw new NoSuchElementException();
            Spectrum current = nextSpectrum;
            advance();
            return current;
        }

        private void advance() {
            nextSpectrum = null;
            if (done) return;

            try {
                while (reader.hasNext()) {
                    int event = reader.next();
                    if (event == XMLStreamConstants.START_ELEMENT && "spectrum".equals(reader.getLocalName())) {
                        Spectrum spec = parseOneSpectrum(reader);
                        if (spec != null) {
                            int ms = spec.getMSLevel();
                            if (ms >= minMSLevel && ms <= maxMSLevel) {
                                // Cache it for potential random access later
                                cache.put(spec.getSpecIndex(), spec);
                                nextSpectrum = spec;
                                return;
                            }
                        }
                    }
                }
                // End of file
                done = true;
                cleanup();
            } catch (XMLStreamException e) {
                done = true;
                cleanup();
                System.err.println("Error iterating mzML: " + e.getMessage());
            }
        }

        private void cleanup() {
            try {
                if (reader != null) reader.close();
                if (inputStream != null) inputStream.close();
            } catch (Exception e) { /* ignore */ }
        }
    }

    // -----------------------------------------------------------------------
    // Logging utility (replaces MzMLAdapter.turnOffLogs)
    // -----------------------------------------------------------------------

    private static boolean logOff = false;

    /**
     * Suppress all logback logging. Called at startup to silence noisy
     * library output.
     */
    public static void turnOffLogs() {
        if (!logOff) {
            LoggerContext context = (LoggerContext) LoggerFactory.getILoggerFactory();
            context.reset();
            Logger rootLogger = context.getLogger(Logger.ROOT_LOGGER_NAME);
            rootLogger.detachAndStopAllAppenders();
            logOff = true;
        }
    }
}
