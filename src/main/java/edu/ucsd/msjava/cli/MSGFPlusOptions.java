package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Typed command-line options for MS-GF+. Picocli reads {@code argv} into
 * the {@code @Option}-annotated fields below; {@link #applyConfigFile}
 * fills in any field the CLI did not set from a {@code -conf} file
 * (CLI takes precedence). {@link #validate} enforces required-input
 * and numeric/enum range invariants. Each {@code effectiveXxx()} accessor
 * returns the user-supplied value or the legacy default.
 *
 * Flag inventory: see {@code .claude/plans/parameter-modernization-flag-inventory.md}.
 */
@Command(
        name = "MS-GF+",
        mixinStandardHelpOptions = true,
        sortOptions = false,
        description = "MS-GF+: peptide identification by database search of mass spectra.")
public final class MSGFPlusOptions {

    // ---------- input (required at runtime, but may be provided via -conf) ----------

    @Option(names = "-s", paramLabel = "SpectrumFile",
            description = "Input spectrum file (*.mzML, *.mgf) or directory of spectra. "
                    + "Required, unless provided via -conf as SpectrumFile=...")
    public File spectrumFile;

    @Option(names = "-d", paramLabel = "DatabaseFile",
            description = "Database file (*.fasta, *.fa, *.faa). "
                    + "Required, unless provided via -conf as DatabaseFile=...")
    public File databaseFile;

    // ---------- optional config + output ----------

    @Option(names = "-conf", paramLabel = "ConfigFile",
            description = "Configuration file path; CLI flags override config file values")
    public File configFile;

    @Option(names = "-o", paramLabel = "OutputFile",
            description = "Output file (*.pin or *.tsv); Default: <SpectrumFileName>.pin")
    public File outputFile;

    @Option(names = "-decoy", paramLabel = "Prefix",
            description = "Decoy protein prefix; Default: XXX")
    public String decoyPrefix;

    // ---------- precursor mass tolerance ----------

    @Option(names = "-t", paramLabel = "Tolerance",
            converter = PrecursorTolerance.Converter.class,
            description = "Precursor mass tolerance, e.g. 20ppm or 0.5Da or 0.5Da,2.5Da; Default: 20ppm. " +
                    "Asymmetric form sets left tolerance (ObsMass < TheoMass) and right tolerance (ObsMass > TheoMass).")
    public PrecursorTolerance precursorTolerance;

    @Option(names = "-u", paramLabel = "Units", hidden = true,
            description = "Tolerance units (legacy): 0=Da, 1=ppm, 2=as written in -t (Default: 2)")
    public Integer precursorToleranceUnits;

    @Option(names = "-ti", paramLabel = "Range",
            converter = IntRange.Converter.class,
            description = "Isotope-error range, e.g. -1,2 (both inclusive); Default: 0,1")
    public IntRange isotopeErrorRange;

    // ---------- threading / parallelism ----------

    @Option(names = "-thread", paramLabel = "N",
            description = "Number of worker threads; Default: number of available cores")
    public Integer numThreads;

    @Option(names = "-tasks", paramLabel = "N",
            description = "Number of tasks: 0=auto, >0=fixed, <0=N*threads; Default: 0")
    public Integer numTasks;

    @Option(names = "-minSpectraPerThread", paramLabel = "N",
            description = "Minimum spectra per thread/task; Default: 250")
    public Integer minSpectraPerThread;

    @Option(names = "-verbose", paramLabel = "N",
            description = "Verbosity: 0=total progress only (Default), 1=per-thread")
    public Integer verbose;

    // ---------- target/decoy + scoring shape ----------

    @Option(names = "-tda", paramLabel = "N",
            description = "Target-decoy strategy: 0=off (Default), 1=concatenated decoy search")
    public Integer tdaStrategy;

    @Option(names = "-m", paramLabel = "ID",
            description = "Fragmentation method ID: 0=as written/CID (Default), 1=CID, 2=ETD, 3=HCD, 4=UVPD")
    public Integer fragMethodId;

    @Option(names = "-inst", paramLabel = "ID",
            description = "Instrument type ID; default depends on registry")
    public Integer instrumentTypeId;

    @Option(names = "-e", paramLabel = "ID",
            description = "Enzyme ID; default depends on registry")
    public Integer enzymeId;

    @Option(names = "-protocol", paramLabel = "ID",
            description = "Protocol ID; default depends on registry")
    public Integer protocolId;

    @Option(names = "-ntt", paramLabel = "N",
            description = "Number of tolerable termini (0..2); Default: 2 (fully tryptic)")
    public Integer numTolerableTermini;

    // ---------- modifications ----------

    @Option(names = "-mod", paramLabel = "ModFile",
            description = "Modification file (also accepts StaticMod=, DynamicMod=, CustomAA= entries via -conf)")
    public File modificationFile;

    // ---------- peptide / charge bounds ----------

    @Option(names = "-minLength", paramLabel = "N",
            description = "Minimum peptide length; Default: 6")
    public Integer minPeptideLength;

    @Option(names = "-maxLength", paramLabel = "N",
            description = "Maximum peptide length; Default: 40")
    public Integer maxPeptideLength;

    @Option(names = "-minCharge", paramLabel = "N",
            description = "Minimum precursor charge; Default: 2")
    public Integer minCharge;

    @Option(names = "-maxCharge", paramLabel = "N",
            description = "Maximum precursor charge; Default: 3")
    public Integer maxCharge;

    @Option(names = "-n", paramLabel = "N",
            description = "Number of matches reported per spectrum; Default: 1")
    public Integer numMatchesPerSpec;

    // ---------- output / features / calibration ----------

    @Option(names = "-addFeatures", paramLabel = "N",
            description = "Include extra features for Percolator: 0=basic (Default), 1=+features")
    public Integer addFeatures;

    @Option(names = "-outputFormat", paramLabel = "Format",
            description = "Output format: pin (Default) or tsv")
    public String outputFormat;

    @Option(names = "-precursorCal", paramLabel = "Mode",
            description = "Precursor calibration mode: auto (Default), on, off")
    public String precursorCalMode;

    @Option(names = "-ccm", paramLabel = "Mass",
            description = "Charge carrier mass; Default: 1.00727649 (proton)")
    public Double chargeCarrierMass;

    @Option(names = "-maxMissedCleavages", paramLabel = "N",
            description = "Max missed cleavages per peptide; -1 = unlimited (Default)")
    public Integer maxMissedCleavages;

    @Option(names = "-numMods", paramLabel = "N",
            description = "Max dynamic mods per peptide; Default: 3")
    public Integer maxNumMods;

    @Option(names = "-allowDenseCentroidedPeaks", paramLabel = "N",
            description = "Allow centroid scans with dense peaks: 0=skip (Default), 1=allow")
    public Integer allowDenseCentroidedPeaks;

    @Option(names = "-msLevel", paramLabel = "Range",
            converter = IntRange.Converter.class,
            description = "MS level or range, e.g. 2 or 2,3; Default: 2,2")
    public IntRange msLevel;

    // ---------- hidden flags ----------

    @Option(names = "-dd", paramLabel = "Dir", hidden = true,
            description = "Database index directory")
    public File dbIndexDir;

    @Option(names = "-index", paramLabel = "Range", hidden = true,
            converter = IntRange.Converter.class,
            description = "Spectrum index range, e.g. 1,1000 (both inclusive)")
    public IntRange specIndexRange;

    @Option(names = "-edgeScore", paramLabel = "N", hidden = true,
            description = "Edge scoring: 0=use (Default), 1=skip")
    public Integer edgeScore;

    @Option(names = "-minNumPeaks", paramLabel = "N", hidden = true,
            description = "Minimum number of peaks per spectrum")
    public Integer minNumPeaks;

    @Option(names = "-iso", paramLabel = "N", hidden = true,
            description = "Number of isoforms to consider per peptide")
    public Integer numIsoforms;

    @Option(names = "-ignoreMetCleavage", paramLabel = "N", hidden = true,
            description = "Ignore N-terminal Met cleavage: 0=consider (Default), 1=ignore")
    public Integer ignoreMetCleavage;

    @Option(names = "-minDeNovoScore", paramLabel = "N", hidden = true,
            description = "Minimum de novo score")
    public Integer minDeNovoScore;

    // ---------- config-file-only entries (populated by applyConfigFile) ----------

    /** {@code DynamicMod=...} entries from the config file (or {@code -mod} file). */
    public final List<String> dynamicMods = new ArrayList<>();
    /** {@code StaticMod=...} entries from the config file (or {@code -mod} file). */
    public final List<String> staticMods = new ArrayList<>();
    /** {@code CustomAA=...} entries from the config file (or {@code -mod} file). */
    public final List<String> customAAs = new ArrayList<>();

    /** Set when {@link #applyConfigFile(File)} encounters {@code MaxNumModsPerPeptide=}
     *  via the legacy alias path; allows the config-file value to feed the
     *  {@link #effectiveMaxNumMods()} default. */
    private Integer configMaxNumMods;

    // ---------- effective-value resolvers (CLI value, else config-file value, else default) ----------

    public int effectiveMinPeptideLength()        { return minPeptideLength        != null ? minPeptideLength        : 6; }
    public int effectiveMaxPeptideLength()        { return maxPeptideLength        != null ? maxPeptideLength        : 40; }
    public int effectiveMinCharge()               { return minCharge               != null ? minCharge               : 2; }
    public int effectiveMaxCharge()               { return maxCharge               != null ? maxCharge               : 3; }
    public int effectiveNumMatchesPerSpec()       { return numMatchesPerSpec       != null ? numMatchesPerSpec       : 1; }
    public int effectiveNumThreads()              { return numThreads              != null ? numThreads              : Runtime.getRuntime().availableProcessors(); }
    public int effectiveNumTasks()                { return numTasks                != null ? numTasks                : 0; }
    public int effectiveMinSpectraPerThread()     { return minSpectraPerThread     != null ? minSpectraPerThread     : 250; }
    public int effectiveVerbose()                 { return verbose                 != null ? verbose                 : 0; }
    public int effectiveTdaStrategy()             { return tdaStrategy             != null ? tdaStrategy             : 0; }
    public int effectiveAddFeatures()             { return addFeatures             != null ? addFeatures             : 0; }
    public int effectiveMaxMissedCleavages()      { return maxMissedCleavages      != null ? maxMissedCleavages      : -1; }
    public int effectiveMaxNumMods()              { return maxNumMods              != null ? maxNumMods              : (configMaxNumMods != null ? configMaxNumMods : 3); }
    public int effectiveAllowDenseCentroidedPeaks() { return allowDenseCentroidedPeaks != null ? allowDenseCentroidedPeaks : 0; }
    public int effectiveNumTolerableTermini()     { return numTolerableTermini     != null ? numTolerableTermini     : 2; }
    public int effectiveEdgeScore()               { return edgeScore               != null ? edgeScore               : 0; }
    public int effectiveIgnoreMetCleavage()       { return ignoreMetCleavage       != null ? ignoreMetCleavage       : 0; }
    public int effectiveMinNumPeaks()             { return minNumPeaks             != null ? minNumPeaks             : edu.ucsd.msjava.sequences.Constants.MIN_NUM_PEAKS_PER_SPECTRUM; }
    public int effectiveNumIsoforms()             { return numIsoforms             != null ? numIsoforms             : edu.ucsd.msjava.sequences.Constants.NUM_VARIANTS_PER_PEPTIDE; }
    public int effectiveMinDeNovoScore()          { return minDeNovoScore          != null ? minDeNovoScore          : edu.ucsd.msjava.sequences.Constants.MIN_DE_NOVO_SCORE; }
    public int effectiveToleranceUnits()          { return precursorToleranceUnits != null ? precursorToleranceUnits : 2; }
    public double effectiveChargeCarrierMass()    { return chargeCarrierMass       != null ? chargeCarrierMass       : 1.00727649; }

    public String effectiveDecoyPrefix()          { return decoyPrefix             != null ? decoyPrefix             : "XXX"; }
    public String effectivePrecursorCalRaw()      { return precursorCalMode        != null ? precursorCalMode        : "auto"; }

    /** 0 = pin (default), 1 = tsv. */
    public int effectiveOutputFormat() {
        if (outputFormat == null) return 0;
        String n = outputFormat.trim().toLowerCase();
        if (n.equals("tsv") || n.equals("1")) return 1;
        return 0;
    }

    public PrecursorTolerance effectivePrecursorTolerance() {
        return precursorTolerance != null ? precursorTolerance : PrecursorTolerance.parse("20ppm");
    }

    public IntRange effectiveIsotopeErrorRange() {
        return isotopeErrorRange != null ? isotopeErrorRange : new IntRange(0, 1);
    }

    public IntRange effectiveMSLevel() {
        return msLevel != null ? msLevel : new IntRange(2, 2);
    }

    public IntRange effectiveSpecIndexRange() {
        return specIndexRange != null ? specIndexRange : new IntRange(1, Integer.MAX_VALUE - 1);
    }

    /** Resolves {@code -m} index to {@link ActivationMethod}. MSGFPlus exposes
     *  0=ASWRITTEN, 1=CID, 2=ETD, 3=HCD, 4=UVPD. The registry also defines
     *  FUSION (merge-mode synthetic method) and PQD, but neither is exposed
     *  as a user-selectable index by MSGFPlus -- FUSION was hidden by the
     *  legacy {@code addFragMethodParam(..., doNotAddMergeMode=true)}, which
     *  shifted UVPD from registry slot 5 down to user-facing index 4. */
    public ActivationMethod effectiveActivationMethod() {
        int idx = fragMethodId != null ? fragMethodId : 0;
        switch (idx) {
            case 0: return ActivationMethod.ASWRITTEN;
            case 1: return ActivationMethod.CID;
            case 2: return ActivationMethod.ETD;
            case 3: return ActivationMethod.HCD;
            case 4: return ActivationMethod.UVPD;
            default: throw new IllegalArgumentException("invalid -m index: " + idx);
        }
    }

    public InstrumentType effectiveInstrumentType() {
        InstrumentType[] all = InstrumentType.getAllRegisteredInstrumentTypes();
        int idx = instrumentTypeId != null ? instrumentTypeId : 0;
        if (idx < 0 || idx >= all.length) throw new IllegalArgumentException("invalid -inst index: " + idx);
        return all[idx];
    }

    public Enzyme effectiveEnzyme() {
        Enzyme[] all = Enzyme.getAllRegisteredEnzymes();
        // TRYPSIN is registered at index 1 (UnspecificCleavage at 0). See Enzyme static init.
        int idx = enzymeId != null ? enzymeId : 1;
        if (idx < 0 || idx >= all.length) throw new IllegalArgumentException("invalid -e index: " + idx);
        return all[idx];
    }

    public Protocol effectiveProtocol() {
        Protocol[] all = Protocol.getAllRegisteredProtocols();
        int idx = protocolId != null ? protocolId : 0;
        if (idx < 0 || idx >= all.length) throw new IllegalArgumentException("invalid -protocol index: " + idx);
        return all[idx];
    }

    // ---------- config-file overlay ----------

    /**
     * Read {@code -conf} config file and populate any fields the CLI did not
     * already set. Recognizes legacy aliases (IsotopeError → IsotopeErrorRange,
     * etc.) and collects repeated {@code DynamicMod=}, {@code StaticMod=},
     * {@code CustomAA=} entries.
     *
     * @return null on success, error string otherwise.
     */
    public String applyConfigFile(File file) {
        unrecognizedConfigEntries = 0;
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;
            int lineNum = 0;
            while ((line = reader.readLine()) != null) {
                lineNum++;
                String trimmed = stripComment(line);
                if (trimmed.isEmpty()) continue;
                int eq = trimmed.indexOf('=');
                if (eq <= 0) continue;
                String rawKey = trimmed.substring(0, eq).trim();
                String value = trimmed.substring(eq + 1).trim();
                String key = canonicalConfigKey(rawKey);
                String err = applyConfigEntry(key, value, file.getName());
                if (err != null) {
                    return "Error parsing line " + lineNum + " of " + file.getName() + ": " + err;
                }
            }
        } catch (IOException e) {
            return "Error reading config file " + file.getPath() + ": " + e.getMessage();
        }
        if (unrecognizedConfigEntries > 0) {
            System.out.println("Valid parameters are described in the example parameter file at " +
                    "https://github.com/MSGFPlus/msgfplus/blob/master/docs/examples/MSGFPlus_Params.txt");
        }
        return null;
    }

    /** Counter incremented inside {@link #applyConfigEntry} whenever an unknown
     *  config-file key is seen; surfaced via the end-of-file URL hint and
     *  reset at the start of each {@link #applyConfigFile} call. */
    private int unrecognizedConfigEntries;

    private String applyConfigEntry(String key, String value, String fileName) {
        // Config-file matching is case-insensitive. canonicalConfigKey()
        // already returns lowercase canonical names, so the switch labels
        // are lowercase too. Repeated mod entries are matched first since
        // they accumulate rather than overwrite.
        switch (key) {
            case "dynamicmod":  if (!value.equalsIgnoreCase("none")) dynamicMods.add(value); return null;
            case "staticmod":   if (!value.equalsIgnoreCase("none")) staticMods.add(value); return null;
            case "customaa":    if (!value.equalsIgnoreCase("none")) customAAs.add(value); return null;
            default: break;
        }
        // Single-valued entries: only fill in if CLI did not set the field.
        try {
            switch (key) {
                case "spectrumfile":               if (spectrumFile == null)              spectrumFile = new File(value); return null;
                case "databasefile":               if (databaseFile == null)              databaseFile = new File(value); return null;
                case "outputfile":                 if (outputFile == null)                outputFile = new File(value); return null;
                case "modificationfilename":
                case "modificationfile":           if (modificationFile == null)          modificationFile = new File(value); return null;
                case "dbindexdir":                 if (dbIndexDir == null)                dbIndexDir = new File(value); return null;
                case "decoyprefix":                if (decoyPrefix == null)               decoyPrefix = value; return null;
                case "precursormasstolerance":     if (precursorTolerance == null)        precursorTolerance = PrecursorTolerance.parse(value); return null;
                case "precursormasstoleranceunits":if (precursorToleranceUnits == null)   precursorToleranceUnits = Integer.parseInt(value); return null;
                case "isotopeerrorrange":          if (isotopeErrorRange == null)         isotopeErrorRange = IntRange.parse(value); return null;
                case "fragmentationmethodid":      if (fragMethodId == null)              fragMethodId = Integer.parseInt(value); return null;
                case "instrumentid":               if (instrumentTypeId == null)          instrumentTypeId = Integer.parseInt(value); return null;
                case "enzymeid":                   if (enzymeId == null)                  enzymeId = Integer.parseInt(value); return null;
                case "protocolid":                 if (protocolId == null)                protocolId = Integer.parseInt(value); return null;
                case "ntt":                        if (numTolerableTermini == null)       numTolerableTermini = Integer.parseInt(value); return null;
                case "minpeplength":               if (minPeptideLength == null)          minPeptideLength = Integer.parseInt(value); return null;
                case "maxpeplength":               if (maxPeptideLength == null)          maxPeptideLength = Integer.parseInt(value); return null;
                case "mincharge":                  if (minCharge == null)                 minCharge = Integer.parseInt(value); return null;
                case "maxcharge":                  if (maxCharge == null)                 maxCharge = Integer.parseInt(value); return null;
                case "nummatchesperspec":          if (numMatchesPerSpec == null)         numMatchesPerSpec = Integer.parseInt(value); return null;
                case "numthreads":                 if (numThreads == null && !value.equalsIgnoreCase("all"))
                                                       numThreads = Integer.parseInt(value); return null;
                case "numtasks":                   if (numTasks == null)                  numTasks = Integer.parseInt(value); return null;
                case "minspectraperthread":        if (minSpectraPerThread == null)       minSpectraPerThread = Integer.parseInt(value); return null;
                case "verbose":                    if (verbose == null)                   verbose = Integer.parseInt(value); return null;
                case "tda":                        if (tdaStrategy == null)               tdaStrategy = Integer.parseInt(value); return null;
                case "addfeatures":                if (addFeatures == null)               addFeatures = Integer.parseInt(value); return null;
                case "outputformat":               if (outputFormat == null)              outputFormat = value; return null;
                case "precursorcal":               if (precursorCalMode == null)          precursorCalMode = value; return null;
                case "chargecarriermass":          if (chargeCarrierMass == null)         chargeCarrierMass = Double.parseDouble(value); return null;
                case "maxmissedcleavages":         if (maxMissedCleavages == null)        maxMissedCleavages = Integer.parseInt(value); return null;
                case "nummods":                    if (maxNumMods == null)                configMaxNumMods = Integer.parseInt(value); return null;
                case "allowdensecentroidedpeaks":  if (allowDenseCentroidedPeaks == null) allowDenseCentroidedPeaks = Integer.parseInt(value); return null;
                case "mslevel":                    if (msLevel == null)                   msLevel = IntRange.parse(value); return null;
                case "specindex":                  if (specIndexRange == null)            specIndexRange = IntRange.parse(value); return null;
                case "edgescore":                  if (edgeScore == null)                 edgeScore = Integer.parseInt(value); return null;
                case "minnumpeaksperspectrum":     if (minNumPeaks == null)               minNumPeaks = Integer.parseInt(value); return null;
                case "numisoforms":                if (numIsoforms == null)               numIsoforms = Integer.parseInt(value); return null;
                case "ignoremetcleavage":          if (ignoreMetCleavage == null)         ignoreMetCleavage = Integer.parseInt(value); return null;
                case "mindenovoscore":             if (minDeNovoScore == null)            minDeNovoScore = Integer.parseInt(value); return null;
                default:
                    if (!key.startsWith("enzymedef")) {
                        System.out.println("Warning, unrecognized parameter '" + key + "=" + value + "' in config file " + fileName);
                        unrecognizedConfigEntries++;
                    }
                    return null;
            }
        } catch (IllegalArgumentException e) {
            return "invalid value for '" + key + "': " + value + " (" + e.getMessage() + ")";
        }
    }

    private static String stripComment(String line) {
        int hash = line.indexOf('#');
        return (hash >= 0 ? line.substring(0, hash) : line).trim();
    }

    /** Normalize legacy / alternate config-file keys to canonical form.
     *  Returns lowercase so {@link #applyConfigEntry} can match
     *  case-insensitively (the legacy {@code ParamManager.parseConfigParamFile}
     *  matched names with {@code equalsIgnoreCase}). Mirrors the alias
     *  rewrites previously in {@code ParamNameEnum.getParamNameFromLine}. */
    private static String canonicalConfigKey(String key) {
        String norm = key.toLowerCase(java.util.Locale.ROOT);
        switch (norm) {
            case "isotopeerror":         return "isotopeerrorrange";
            case "targetdecoyanalysis":  return "tda";
            case "fragmentationmethod":  return "fragmentationmethodid";
            case "instrument":           return "instrumentid";
            case "enzyme":               return "enzymeid";
            case "protocol":             return "protocolid";
            case "numtolerabletermini":  return "ntt";
            case "minnumpeaks":          return "minnumpeaksperspectrum";
            case "maxnummods":           return "nummods";
            case "maxnummodsperpeptide": return "nummods";
            case "minlength":            return "minpeplength";
            case "minpeptidelength":     return "minpeplength";
            case "maxlength":            return "maxpeplength";
            case "maxpeptidelength":     return "maxpeplength";
            case "pmtolerance":          return "precursormasstolerance";
            case "parentmasstolerance":  return "precursormasstolerance";
            default:                     return norm;
        }
    }

    /** Validates required-input invariants and the numeric/enum range
     *  constraints the legacy {@code IntParameter.minValue}/{@code maxValue}
     *  and {@code EnumParameter} machinery used to enforce. Returns
     *  {@code null} on success or a user-facing error string otherwise.
     *
     *  <p>Required: {@code -s} and {@code -d} (either via CLI or {@code -conf}).
     *  Numeric flags must satisfy their original lower bounds; enum-shaped
     *  flags must fall in their defined index range. */
    public String validate() {
        if (spectrumFile == null) return "Spectrum file is not defined; use -s at the command line or SpectrumFile in a config file";
        if (databaseFile == null) return "Database file is not defined; use -d at the command line or DatabaseFile in a config file";

        String err;
        if ((err = checkMin("-thread",                    numThreads,                1))    != null) return err;
        if ((err = checkMin("-tasks",                     numTasks,                  -10))  != null) return err;
        if ((err = checkMin("-minSpectraPerThread",       minSpectraPerThread,       1))    != null) return err;
        if ((err = checkMin("-minLength",                 minPeptideLength,          1))    != null) return err;
        if ((err = checkMin("-maxLength",                 maxPeptideLength,          1))    != null) return err;
        if ((err = checkMin("-minCharge",                 minCharge,                 1))    != null) return err;
        if ((err = checkMin("-maxCharge",                 maxCharge,                 1))    != null) return err;
        if ((err = checkMin("-n",                         numMatchesPerSpec,         1))    != null) return err;
        if ((err = checkMin("-maxMissedCleavages",        maxMissedCleavages,        -1))   != null) return err;
        if ((err = checkMin("-numMods",                   maxNumMods,                0))    != null) return err;
        if ((err = checkMin("-minNumPeaks",               minNumPeaks,               0))    != null) return err;
        if ((err = checkMin("-iso",                       numIsoforms,               0))    != null) return err;
        if ((err = checkMin("-minDeNovoScore",            minDeNovoScore,            Integer.MIN_VALUE)) != null) return err;

        if ((err = checkRange("-ntt",                     numTolerableTermini,        0, 2)) != null) return err;
        if ((err = checkRange("-tda",                     tdaStrategy,                0, 1)) != null) return err;
        if ((err = checkRange("-verbose",                 verbose,                    0, 1)) != null) return err;
        if ((err = checkRange("-addFeatures",             addFeatures,                0, 1)) != null) return err;
        if ((err = checkRange("-allowDenseCentroidedPeaks", allowDenseCentroidedPeaks, 0, 1)) != null) return err;
        if ((err = checkRange("-edgeScore",               edgeScore,                  0, 1)) != null) return err;
        if ((err = checkRange("-ignoreMetCleavage",       ignoreMetCleavage,          0, 1)) != null) return err;
        if ((err = checkRange("-u",                       precursorToleranceUnits,    0, 2)) != null) return err;

        if (chargeCarrierMass != null && chargeCarrierMass <= 0.1) {
            return "Invalid value for parameter -ccm: " + chargeCarrierMass + " (must be > 0.1)";
        }

        if (fragMethodId != null && (fragMethodId < 0 || fragMethodId > 4)) {
            return "Invalid value for parameter -m: " + fragMethodId + " (valid: 0..4)";
        }
        int instMax = ActivationMethodAvailability.instCount() - 1;
        if (instrumentTypeId != null && (instrumentTypeId < 0 || instrumentTypeId > instMax)) {
            return "Invalid value for parameter -inst: " + instrumentTypeId + " (valid: 0.." + instMax + ")";
        }
        int enzMax = Enzyme.getAllRegisteredEnzymes().length - 1;
        if (enzymeId != null && (enzymeId < 0 || enzymeId > enzMax)) {
            return "Invalid value for parameter -e: " + enzymeId + " (valid: 0.." + enzMax + ")";
        }
        int protMax = Protocol.getAllRegisteredProtocols().length - 1;
        if (protocolId != null && (protocolId < 0 || protocolId > protMax)) {
            return "Invalid value for parameter -protocol: " + protocolId + " (valid: 0.." + protMax + ")";
        }
        return null;
    }

    private static String checkMin(String flag, Integer value, int min) {
        if (value == null) return null;
        if (value < min) return "Invalid value for parameter " + flag + ": " + value + " (must be >= " + min + ")";
        return null;
    }

    private static String checkRange(String flag, Integer value, int min, int max) {
        if (value == null) return null;
        if (value < min || value > max) return "Invalid value for parameter " + flag + ": " + value + " (valid: " + min + ".." + max + ")";
        return null;
    }

    /** Helper that hides the {@link InstrumentType#getAllRegisteredInstrumentTypes}
     *  call from {@code validate()} so the import block stays minimal. */
    private static final class ActivationMethodAvailability {
        static int instCount() { return InstrumentType.getAllRegisteredInstrumentTypes().length; }
    }

    /** Mutator used by {@code AminoAcidSet} when the parsed mod metadata
     *  changes the effective max-num-mods (the AA set is authoritative once
     *  loaded). Mirrors the legacy {@code ParamManager.setMaxNumMods}. */
    public void setMaxNumModsFromMetadata(int n) {
        this.maxNumMods = n;
    }
}
