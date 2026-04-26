package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.params.ParamManager.ParamNameEnum;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;

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

    // ---------- required input ----------

    @Option(names = "-s", required = true, paramLabel = "SpectrumFile",
            description = "Input spectrum file (*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl, *_dta.txt) or directory of spectra")
    public File spectrumFile;

    @Option(names = "-d", required = true, paramLabel = "DatabaseFile",
            description = "Database file (*.fasta, *.fa, *.faa)")
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
            description = "Precursor mass tolerance, e.g. 20ppm or 0.5Da or 0.5Da,2.5Da; Default: 20ppm. " +
                    "Asymmetric form sets left tolerance (ObsMass < TheoMass) and right tolerance (ObsMass > TheoMass).")
    public String precursorTolerance;

    @Option(names = "-u", paramLabel = "Units", hidden = true,
            description = "Tolerance units (legacy): 0=Da, 1=ppm, 2=as written in -t (Default: 2)")
    public Integer precursorToleranceUnits;

    @Option(names = "-ti", paramLabel = "Range",
            description = "Isotope-error range, e.g. -1,2 (both inclusive); Default: 0,1")
    public String isotopeErrorRange;

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
            description = "MS level or range, e.g. 2 or 2,3; Default: 2,2")
    public String msLevel;

    // ---------- hidden flags ----------

    @Option(names = "-dd", paramLabel = "Dir", hidden = true,
            description = "Database index directory")
    public File dbIndexDir;

    @Option(names = "-index", paramLabel = "Range", hidden = true,
            description = "Spectrum index range, e.g. 1,1000 (both inclusive)")
    public String specIndexRange;

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
}
