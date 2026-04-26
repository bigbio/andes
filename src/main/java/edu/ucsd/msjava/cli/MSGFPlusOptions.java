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
 * Typed command-line options for MS-GF+. Replaces the imperative
 * {@code addParameter()} calls in {@code ParamManager.addMSGFPlusParams()}
 * with declarative picocli annotations.
 *
 * Phase 1 scope: every flag from {@link ParamNameEnum} that
 * {@code addMSGFPlusParams()} registers, parsed into typed fields.
 * Complex domain types (Tolerance, IntRange, dynamic enums) are
 * captured here as raw strings; the adapter at
 * {@code MSGFPlusOptionsAdapter} round-trips them through the existing
 * {@code params.Parameter#parse(String)} hierarchy to populate a
 * {@code ParamManager} that {@code SearchParams.parse(ParamManager)}
 * can consume unchanged. Phase 3 collapses that round-trip away.
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
            description = "Input spectrum file (*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl, *_dta.txt) or directory of spectra. "
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
            description = "Decoy protein prefix; Default: DECOY_")
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
            description = "Fragmentation method ID: 0=as written/CID (Default), 1=CID, 2=ETD, 3=HCD")
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
     *  0=ASWRITTEN, 1=CID, 2=ETD, 3=HCD (FUSION is excluded by
     *  {@code addFragMethodParam(..., doNotAddMergeMode=true)}). */
    public ActivationMethod effectiveActivationMethod() {
        int idx = fragMethodId != null ? fragMethodId : 0;
        switch (idx) {
            case 0: return ActivationMethod.ASWRITTEN;
            case 1: return ActivationMethod.CID;
            case 2: return ActivationMethod.ETD;
            case 3: return ActivationMethod.HCD;
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
        return null;
    }

    private String applyConfigEntry(String key, String value, String fileName) {
        // Repeated entries: collect into lists. "none" is treated as no entry.
        if (key.equalsIgnoreCase("DynamicMod")) {
            if (!value.equalsIgnoreCase("none")) dynamicMods.add(value);
            return null;
        }
        if (key.equalsIgnoreCase("StaticMod")) {
            if (!value.equalsIgnoreCase("none")) staticMods.add(value);
            return null;
        }
        if (key.equalsIgnoreCase("CustomAA")) {
            if (!value.equalsIgnoreCase("none")) customAAs.add(value);
            return null;
        }
        // Single-valued entries: only fill in if CLI did not set the field.
        try {
            switch (key) {
                case "SpectrumFile":           if (spectrumFile == null)             spectrumFile = new File(value); return null;
                case "DatabaseFile":           if (databaseFile == null)             databaseFile = new File(value); return null;
                case "OutputFile":             if (outputFile == null)               outputFile = new File(value); return null;
                case "ModificationFileName":
                case "ModificationFile":       if (modificationFile == null)         modificationFile = new File(value); return null;
                case "DBIndexDir":             if (dbIndexDir == null)               dbIndexDir = new File(value); return null;
                case "DecoyPrefix":            if (decoyPrefix == null)              decoyPrefix = value; return null;
                case "PrecursorMassTolerance": if (precursorTolerance == null)       precursorTolerance = PrecursorTolerance.parse(value); return null;
                case "PrecursorMassToleranceUnits":
                                               if (precursorToleranceUnits == null)  precursorToleranceUnits = Integer.parseInt(value); return null;
                case "IsotopeErrorRange":      if (isotopeErrorRange == null)        isotopeErrorRange = IntRange.parse(value); return null;
                case "FragmentationMethodID":  if (fragMethodId == null)             fragMethodId = Integer.parseInt(value); return null;
                case "InstrumentID":           if (instrumentTypeId == null)         instrumentTypeId = Integer.parseInt(value); return null;
                case "EnzymeID":               if (enzymeId == null)                 enzymeId = Integer.parseInt(value); return null;
                case "ProtocolID":             if (protocolId == null)               protocolId = Integer.parseInt(value); return null;
                case "NTT":                    if (numTolerableTermini == null)      numTolerableTermini = Integer.parseInt(value); return null;
                case "MinPepLength":           if (minPeptideLength == null)         minPeptideLength = Integer.parseInt(value); return null;
                case "MaxPepLength":           if (maxPeptideLength == null)         maxPeptideLength = Integer.parseInt(value); return null;
                case "MinCharge":              if (minCharge == null)                minCharge = Integer.parseInt(value); return null;
                case "MaxCharge":              if (maxCharge == null)                maxCharge = Integer.parseInt(value); return null;
                case "NumMatchesPerSpec":      if (numMatchesPerSpec == null)        numMatchesPerSpec = Integer.parseInt(value); return null;
                case "NumThreads":             if (numThreads == null)               { if (!value.equalsIgnoreCase("all")) numThreads = Integer.parseInt(value); } return null;
                case "NumTasks":               if (numTasks == null)                 numTasks = Integer.parseInt(value); return null;
                case "MinSpectraPerThread":    if (minSpectraPerThread == null)      minSpectraPerThread = Integer.parseInt(value); return null;
                case "Verbose":                if (verbose == null)                  verbose = Integer.parseInt(value); return null;
                case "TDA":                    if (tdaStrategy == null)              tdaStrategy = Integer.parseInt(value); return null;
                case "AddFeatures":            if (addFeatures == null)              addFeatures = Integer.parseInt(value); return null;
                case "OutputFormat":           if (outputFormat == null)             outputFormat = value; return null;
                case "PrecursorCal":           if (precursorCalMode == null)         precursorCalMode = value; return null;
                case "ChargeCarrierMass":      if (chargeCarrierMass == null)        chargeCarrierMass = Double.parseDouble(value); return null;
                case "MaxMissedCleavages":     if (maxMissedCleavages == null)       maxMissedCleavages = Integer.parseInt(value); return null;
                case "NumMods":                if (maxNumMods == null)               configMaxNumMods = Integer.parseInt(value); return null;
                case "AllowDenseCentroidedPeaks":
                                               if (allowDenseCentroidedPeaks == null) allowDenseCentroidedPeaks = Integer.parseInt(value); return null;
                case "MSLevel":                if (msLevel == null)                  msLevel = IntRange.parse(value); return null;
                case "SpecIndex":              if (specIndexRange == null)           specIndexRange = IntRange.parse(value); return null;
                case "EdgeScore":              if (edgeScore == null)                edgeScore = Integer.parseInt(value); return null;
                case "MinNumPeaksPerSpectrum": if (minNumPeaks == null)              minNumPeaks = Integer.parseInt(value); return null;
                case "NumIsoforms":            if (numIsoforms == null)              numIsoforms = Integer.parseInt(value); return null;
                case "IgnoreMetCleavage":      if (ignoreMetCleavage == null)        ignoreMetCleavage = Integer.parseInt(value); return null;
                case "MinDeNovoScore":         if (minDeNovoScore == null)           minDeNovoScore = Integer.parseInt(value); return null;
                default:
                    if (!key.toLowerCase().startsWith("enzymedef")) {
                        System.out.println("Warning, unrecognized parameter '" + key + "=" + value + "' in config file " + fileName);
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
     *  Mirrors the rewrites previously in {@code ParamNameEnum.getParamNameFromLine}. */
    private static String canonicalConfigKey(String key) {
        if (key.equalsIgnoreCase("IsotopeError"))          return "IsotopeErrorRange";
        if (key.equalsIgnoreCase("TargetDecoyAnalysis"))   return "TDA";
        if (key.equalsIgnoreCase("FragmentationMethod"))   return "FragmentationMethodID";
        if (key.equalsIgnoreCase("Instrument"))            return "InstrumentID";
        if (key.equalsIgnoreCase("Enzyme"))                return "EnzymeID";
        if (key.equalsIgnoreCase("Protocol"))              return "ProtocolID";
        if (key.equalsIgnoreCase("NumTolerableTermini"))   return "NTT";
        if (key.equalsIgnoreCase("MinNumPeaks"))           return "MinNumPeaksPerSpectrum";
        if (key.equalsIgnoreCase("MaxNumMods"))            return "NumMods";
        if (key.equalsIgnoreCase("MaxNumModsPerPeptide"))  return "NumMods";
        if (key.equalsIgnoreCase("minLength"))             return "MinPepLength";
        if (key.equalsIgnoreCase("MinPeptideLength"))      return "MinPepLength";
        if (key.equalsIgnoreCase("maxLength"))             return "MaxPepLength";
        if (key.equalsIgnoreCase("MaxPeptideLength"))      return "MaxPepLength";
        if (key.equalsIgnoreCase("PMTolerance"))           return "PrecursorMassTolerance";
        if (key.equalsIgnoreCase("ParentMassTolerance"))   return "PrecursorMassTolerance";
        return key;
    }

    /** Validates required-input invariants that the CLI alone can't enforce
     *  (since {@code -s}/{@code -d} may come from {@code -conf}). */
    public String validateRequired() {
        if (spectrumFile == null) return "Spectrum file is not defined; use -s at the command line or SpectrumFile in a config file";
        if (databaseFile == null) return "Database file is not defined; use -d at the command line or DatabaseFile in a config file";
        return null;
    }

    /** Mutator used by {@code AminoAcidSet} when the parsed mod metadata
     *  changes the effective max-num-mods (the AA set is authoritative once
     *  loaded). Mirrors the legacy {@code ParamManager.setMaxNumMods}. */
    public void setMaxNumModsFromMetadata(int n) {
        this.maxNumMods = n;
    }
}
