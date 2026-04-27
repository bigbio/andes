package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.cli.IntRange;
import edu.ucsd.msjava.cli.MSGFPlusOptions;
import edu.ucsd.msjava.cli.OutputFormat;
import edu.ucsd.msjava.cli.PrecursorTolerance;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static edu.ucsd.msjava.msutil.Composition.POTASSIUM_CHARGE_CARRIER_MASS;
import static edu.ucsd.msjava.msutil.Composition.PROTON;
import static edu.ucsd.msjava.msutil.Composition.SODIUM_CHARGE_CARRIER_MASS;

public class SearchParams {

    /**
     * Two-pass precursor mass calibration (P2-cal) mode.
     *
     * <ul>
     *     <li>{@link #AUTO} (default) — run the pre-pass, apply the learned shift
     *         only if at least 200 high-confidence PSMs are collected; otherwise
     *         fall through with a 0 ppm shift.</li>
     *     <li>{@link #ON} — run the pre-pass and always apply the learned shift,
     *         even when fewer than 200 confident PSMs are collected.</li>
     *     <li>{@link #OFF} — skip calibration entirely. The code path MUST be
     *         bit-identical to a baseline build without the flag.</li>
     * </ul>
     */
    public enum PrecursorCalMode {
        AUTO,
        ON,
        OFF
    }

    private List<DBSearchIOFiles> dbSearchIOList;
    private File databaseFile;
    private String decoyProteinPrefix;
    private Tolerance leftPrecursorMassTolerance;
    private Tolerance rightPrecursorMassTolerance;
    private int minIsotopeError;
    private int maxIsotopeError;
    private Enzyme enzyme;
    private int numTolerableTermini;
    private ActivationMethod activationMethod;
    private InstrumentType instType;
    private Protocol protocol;
    private AminoAcidSet aaSet;
    private int numMatchesPerSpec;
    private int startSpecIndex;
    private int endSpecIndex;
    private boolean useTDA;
    private boolean ignoreMetCleavage;
    private int minPeptideLength;
    private int maxPeptideLength;
    private int maxNumVariantsPerPeptide;
    private int minCharge;
    private int maxCharge;
    private int numThreads;
    private int numTasks;
    private int minSpectraPerThread;
    private boolean verbose;
    private boolean doNotUseEdgeScore;
    private File dbIndexDir;
    private boolean outputAdditionalFeatures;
    private int minNumPeaksPerSpectrum;
    private int minDeNovoScore;
    private double chargeCarrierMass;
    private int maxMissedCleavages;
    private int maxNumMods;
    private boolean allowDenseCentroidedPeaks;
    private int minMSLevel;
    private int maxMSLevel;
    private OutputFormat outputFormat;
    private PrecursorCalMode precursorCalMode = PrecursorCalMode.AUTO;

    public SearchParams() {
    }

    /**
     * Returns the configured precursor mass calibration mode; defaults
     * to {@link PrecursorCalMode#AUTO}.
     */
    public PrecursorCalMode getPrecursorCalMode() {
        return precursorCalMode;
    }

    // Used by MS-GF+
    public List<DBSearchIOFiles> getDBSearchIOList() {
        return dbSearchIOList;
    }

    // Used by MS-GF+
    public File getDatabaseFile() {
        return databaseFile;
    }

    // Used by MS-GF+
    public String getDecoyProteinPrefix() {
        return decoyProteinPrefix;
    }

    // Used by MS-GF+
    public Tolerance getLeftPrecursorMassTolerance() {
        return leftPrecursorMassTolerance;
    }

    // Used by MS-GF+
    public Tolerance getRightPrecursorMassTolerance() {
        return rightPrecursorMassTolerance;
    }

    // Used by MS-GF+
    public int getMinIsotopeError() {
        return minIsotopeError;
    }

    // Used by MS-GF+
    public int getMaxIsotopeError() {
        return maxIsotopeError;
    }

    // Used by MS-GF+
    public Enzyme getEnzyme() {
        return enzyme;
    }

    public int getNumTolerableTermini() {
        return numTolerableTermini;
    }

    // Used by MS-GF+
    public ActivationMethod getActivationMethod() {
        return activationMethod;
    }

    // Used by MS-GF+
    public InstrumentType getInstType() {
        return instType;
    }

    // Used by MS-GF+
    public Protocol getProtocol() {
        return protocol;
    }

    // Used by MS-GF+
    public AminoAcidSet getAASet() {
        return aaSet;
    }

    // Used by MS-GF+
    public int getNumMatchesPerSpec() {
        return numMatchesPerSpec;
    }

    // Used by MS-GF+
    public int getStartSpecIndex() {
        return startSpecIndex;
    }

    // Used by MS-GF+
    public int getEndSpecIndex() {
        return endSpecIndex;
    }

    // Used by MS-GF+
    public boolean useTDA() {
        return useTDA;
    }

    // Used by MS-GF+
    public boolean ignoreMetCleavage() {
        return ignoreMetCleavage;
    }

    // Used by MS-GF+
    public int getMinPeptideLength() {
        return minPeptideLength;
    }

    // Used by MS-GF+
    public int getMaxPeptideLength() {
        return maxPeptideLength;
    }

    // Used by MS-GF+
    public int getMaxNumVariantsPerPeptide() {
        return maxNumVariantsPerPeptide;
    }

    // Used by MS-GF+
    public int getMinCharge() {
        return minCharge;
    }

    // Used by MS-GF+
    public int getMaxCharge() {
        return maxCharge;
    }

    // Used by MS-GF+
    public int getNumThreads() {
        return numThreads;
    }

    public int getNumTasks() {
        return numTasks;
    }

    public int getMinSpectraPerThread() {
        return minSpectraPerThread;
    }

    public boolean getVerbose() {
        return verbose;
    }

    // Used by MS-GF+
    public boolean doNotUseEdgeScore() {
        return doNotUseEdgeScore;
    }

    // Used by MS-GF+
    public File getDBIndexDir() {
        return dbIndexDir;
    }

    public boolean outputAdditionalFeatures() {
        return outputAdditionalFeatures;
    }

    // Used by MS-GF+
    public int getMinNumPeaksPerSpectrum() {
        return minNumPeaksPerSpectrum;
    }

    // Used by MS-GF+
    public int getMinDeNovoScore() {
        return minDeNovoScore;
    }

    public double getChargeCarrierMass() {
        return chargeCarrierMass;
    }

    // Used by MS-GF+
    public int getMaxMissedCleavages() {
        return maxMissedCleavages;
    }

    // Used by MS-GF+
    public boolean getAllowDenseCentroidedPeaks() {
        return allowDenseCentroidedPeaks;
    }

    // Used by MS-GF+
    public int getMinMSLevel() {
        return minMSLevel;
    }

    // Used by MS-GF+
    public int getMaxMSLevel() {
        return maxMSLevel;
    }

    public boolean writeTsv() {
        return outputFormat == OutputFormat.TSV;
    }

    /**
     * Look for # in dataLine
     * If present, remove that character and any comment after it
     *
     * @param dataLine
     * @return dataLine without the comment
     */
    public static String getConfigLineWithoutComment(String dataLine) {
        return MSGFPlusOptions.stripComment(dataLine);
    }

    /**
     * Build a SearchParams from the typed CLI/config-file model. Reads {@code -conf}
     * (when set) via {@link MSGFPlusOptions#applyConfigFile(File)} so any unset CLI
     * fields are filled from the config file before the rest of the build runs.
     *
     * @return null on success; user-facing error string otherwise.
     */
    public String parse(MSGFPlusOptions opts) {
        // Apply config-file overlay first: fills in any opts.* fields the CLI did
        // not set, plus collects DynamicMod/StaticMod/CustomAA into opts.*Mods lists.
        if (opts.configFile != null) {
            String err = opts.applyConfigFile(opts.configFile);
            if (err != null) return err;
        }

        // Required-input + numeric/enum range check now that CLI +
        // config-file have both run. Catches things like -m 99 with a
        // user-facing error instead of the IllegalArgumentException
        // the resolver would otherwise raise during search setup.
        String requiredErr = opts.validate();
        if (requiredErr != null) return requiredErr;

        chargeCarrierMass = opts.chargeCarrierMass != null ? opts.chargeCarrierMass : 1.00727649;
        Composition.setChargeCarrierMass(chargeCarrierMass);

        // Read outputFormat up-front so the default-output-file extension logic
        // below sees the user-supplied value, not the field's zero initializer.
        outputFormat = opts.effectiveOutputFormat();

        File specPath = opts.spectrumFile;
        if (!specPath.exists()) {
            return "Spectrum file not found: " + specPath.getPath();
        }

        dbSearchIOList = new ArrayList<>();
        String defaultExt = outputFormat == OutputFormat.TSV ? ".tsv" : ".pin";

        if (!specPath.isDirectory()) {
            SpecFileFormat specFormat = SpecFileFormat.getSpecFileFormat(specPath.getName());
            if (!isSupportedSpectrumFormat(specFormat)) {
                return "Spectrum file extension does not match a supported format (*.mzML, *.mgf): " + specPath.getName();
            }
            File outputFile = opts.outputFile;
            if (outputFile == null) {
                String outputFilePath = specPath.getPath().substring(0, specPath.getPath().lastIndexOf('.')) + defaultExt;
                outputFile = new File(outputFilePath);
            }
            dbSearchIOList.add(new DBSearchIOFiles(specPath, specFormat, outputFile));
        } else {
            for (File f : specPath.listFiles()) {
                SpecFileFormat specFormat = SpecFileFormat.getSpecFileFormat(f.getName());
                if (isSupportedSpectrumFormat(specFormat)) {
                    String outputFileName = f.getName().substring(0, f.getName().lastIndexOf('.')) + defaultExt;
                    File outputFile = new File(outputFileName);
                    dbSearchIOList.add(new DBSearchIOFiles(f, specFormat, outputFile));
                }
            }
        }

        databaseFile = opts.databaseFile;
        decoyProteinPrefix = opts.decoyPrefix != null ? opts.decoyPrefix : "XXX";

        PrecursorTolerance tol = opts.precursorTolerance != null ? opts.precursorTolerance : PrecursorTolerance.parse("20ppm");
        leftPrecursorMassTolerance = tol.left;
        rightPrecursorMassTolerance = tol.right;

        int toleranceUnit = opts.precursorToleranceUnits != null ? opts.precursorToleranceUnits : 2;
        if (toleranceUnit != 2) {
            boolean isTolerancePPM = toleranceUnit != 0;
            leftPrecursorMassTolerance = new Tolerance(leftPrecursorMassTolerance.getValue(), isTolerancePPM);
            rightPrecursorMassTolerance = new Tolerance(rightPrecursorMassTolerance.getValue(), isTolerancePPM);
        }

        IntRange isotope = opts.isotopeErrorRange != null ? opts.isotopeErrorRange : new IntRange(0, 1);
        this.minIsotopeError = isotope.min;
        this.maxIsotopeError = isotope.max;

        if (rightPrecursorMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f ||
                leftPrecursorMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f) {
            minIsotopeError = maxIsotopeError = 0;
        }

        enzyme = opts.effectiveEnzyme();
        numTolerableTermini = opts.numTolerableTermini != null ? opts.numTolerableTermini : 2;
        activationMethod = opts.effectiveActivationMethod();
        instType = opts.effectiveInstrumentType();
        if (activationMethod == ActivationMethod.HCD
                && instType != InstrumentType.HIGH_RESOLUTION_LTQ
                && instType != InstrumentType.QEXACTIVE) {
            instType = InstrumentType.QEXACTIVE; // default to Q-Exactive for HCD
        }
        protocol = opts.effectiveProtocol();

        aaSet = null;
        File modFile = opts.modificationFile;
        boolean hasConfigMods = !opts.dynamicMods.isEmpty()
                || !opts.staticMods.isEmpty()
                || !opts.customAAs.isEmpty();

        if (modFile == null && !hasConfigMods) {
            aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        } else {
            if (modFile != null) {
                String modFileName = modFile.getName();
                String ext = modFileName.substring(modFileName.lastIndexOf('.') + 1);
                if (ext.equalsIgnoreCase("xml")) {
                    aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
                } else {
                    aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath(), opts);
                }
            } else {
                List<String> mods = new ArrayList<>(opts.staticMods.size() + opts.dynamicMods.size());
                mods.addAll(opts.staticMods);
                mods.addAll(opts.dynamicMods);
                aaSet = AminoAcidSet.getAminoAcidSetFromModEntries(
                        opts.configFile != null ? opts.configFile.getName() : "config",
                        opts.customAAs, mods, opts);
            }

            if (protocol == Protocol.AUTOMATIC) {
                if (aaSet.containsITRAQ()) {
                    protocol = aaSet.containsPhosphorylation() ? Protocol.ITRAQPHOSPHO : Protocol.ITRAQ;
                } else if (aaSet.containsTMT()) {
                    protocol = Protocol.TMT;
                } else {
                    protocol = aaSet.containsPhosphorylation() ? Protocol.PHOSPHORYLATION : Protocol.STANDARD;
                }
            }
        }

        numMatchesPerSpec = opts.numMatchesPerSpec != null ? opts.numMatchesPerSpec : 1;

        IntRange specIdx = opts.specIndexRange != null ? opts.specIndexRange : new IntRange(1, Integer.MAX_VALUE - 1);
        startSpecIndex = specIdx.min;
        endSpecIndex = specIdx.max;

        useTDA = opts.effectiveTdaStrategy() == 1;
        ignoreMetCleavage = (opts.ignoreMetCleavage != null ? opts.ignoreMetCleavage : 0) == 1;
        outputAdditionalFeatures = (opts.addFeatures != null ? opts.addFeatures : 0) == 1;

        minPeptideLength = opts.effectiveMinPeptideLength();
        maxPeptideLength = opts.effectiveMaxPeptideLength();
        maxNumVariantsPerPeptide = opts.numIsoforms != null ? opts.numIsoforms : edu.ucsd.msjava.sequences.Constants.NUM_VARIANTS_PER_PEPTIDE;

        if (minPeptideLength > maxPeptideLength) {
            return "MinPepLength must not be larger than MaxPepLength";
        }

        minCharge = opts.effectiveMinCharge();
        maxCharge = opts.effectiveMaxCharge();
        if (minCharge > maxCharge) {
            return "MinCharge must not be larger than MaxCharge";
        }

        numThreads = opts.numThreads != null ? opts.numThreads : Runtime.getRuntime().availableProcessors();
        numTasks = opts.numTasks != null ? opts.numTasks : 0;
        minSpectraPerThread = opts.effectiveMinSpectraPerThread();
        verbose = opts.effectiveVerbose() == 1;
        doNotUseEdgeScore = (opts.edgeScore != null ? opts.edgeScore : 0) == 1;

        dbIndexDir = opts.dbIndexDir;
        minNumPeaksPerSpectrum = opts.minNumPeaks != null ? opts.minNumPeaks : edu.ucsd.msjava.sequences.Constants.MIN_NUM_PEAKS_PER_SPECTRUM;
        minDeNovoScore = opts.minDeNovoScore != null ? opts.minDeNovoScore : edu.ucsd.msjava.sequences.Constants.MIN_DE_NOVO_SCORE;

        maxMissedCleavages = opts.maxMissedCleavages != null ? opts.maxMissedCleavages : -1;
        if (maxMissedCleavages > -1 && enzyme.getName().equals("UnspecificCleavage")) {
            return "Cannot specify a MaxMissedCleavages when using unspecific cleavage enzyme";
        } else if (maxMissedCleavages > -1 && enzyme.getName().equals("NoCleavage")) {
            return "Cannot specify a MaxMissedCleavages when using no cleavage enzyme";
        }

        allowDenseCentroidedPeaks = (opts.allowDenseCentroidedPeaks != null ? opts.allowDenseCentroidedPeaks : 0) == 1;
        precursorCalMode = opts.precursorCalMode != null ? opts.precursorCalMode : PrecursorCalMode.AUTO;

        IntRange ms = opts.msLevel != null ? opts.msLevel : new IntRange(2, 2);
        minMSLevel = ms.min;
        maxMSLevel = ms.max;

        maxNumMods = opts.effectiveMaxNumMods();
        int maxNumModsCompare = aaSet.getMaxNumberOfVariableModificationsPerPeptide();
        if (maxNumMods != maxNumModsCompare) {
            System.err.println("Error, code bug: MaxNumModsPerPeptide tracked by MSGFPlusOptions ("
                    + maxNumMods + ") does not match value tracked by AminoAcidSet ("
                    + maxNumModsCompare + ")");
            System.exit(-1);
        }

        Modification.setModIdentifiers();
        return null;
    }

    /** Spectrum-format whitelist: only mzML and MGF are supported. */
    private static boolean isSupportedSpectrumFormat(SpecFileFormat fmt) {
        return fmt == SpecFileFormat.MZML
                || fmt == SpecFileFormat.MGF;
    }


    @Override
    public String toString() {
        StringBuffer buf = new StringBuffer();

//		buf.append("Spectrum File(s):\n");
//		for(DBSearchIOFiles ioFile : this.dbSearchIOList)
//		{
//			buf.append("\t"+ioFile.getSpecFile().getAbsolutePath()+"\n");
//		}
//		buf.append("Database File: " + this.databaseFile.getAbsolutePath() + "\n");

        buf.append("\tPrecursorMassTolerance: ");
        if (leftPrecursorMassTolerance.equals(rightPrecursorMassTolerance)) {
            buf.append(leftPrecursorMassTolerance);
        } else {
            buf.append("[" + leftPrecursorMassTolerance + "," + rightPrecursorMassTolerance + "]");
        }
        buf.append("\n");

        buf.append("\tIsotopeError: " + this.minIsotopeError + "," + this.maxIsotopeError + "\n");
        buf.append("\tTargetDecoyAnalysis: " + this.useTDA + "\n");
        buf.append("\tFragmentationMethod: " + this.activationMethod + "\n");
        buf.append("\tInstrument: " + (instType == null ? "null" : this.instType.getNameAndDescription()) + "\n");
        buf.append("\tEnzyme: " + (enzyme == null ? "null" : this.enzyme.getName()) + "\n");

        String customEnzymeFile = Enzyme.getCustomEnzymeFilePath();
        if (customEnzymeFile != null && !customEnzymeFile.isEmpty()) {
            buf.append("\tEnzyme file: " + customEnzymeFile + "\n");
        }

        ArrayList<String> customEnzymeMessages = Enzyme.getCustomEnzymeMessages();
        for (String message : customEnzymeMessages) {
            buf.append("\tEnzyme info: " + message + "\n");
        }

        buf.append("\tProtocol: " + (protocol == null ? "null" : this.protocol.getName()) + "\n");
        buf.append("\tNumTolerableTermini: " + this.numTolerableTermini + "\n");
        buf.append("\tIgnoreMetCleavage: " + this.ignoreMetCleavage + "\n");
        buf.append("\tMinPepLength: " + this.minPeptideLength + "\n");
        buf.append("\tMaxPepLength: " + this.maxPeptideLength + "\n");
        buf.append("\tMinCharge: " + this.minCharge + "\n");
        buf.append("\tMaxCharge: " + this.maxCharge + "\n");
        buf.append("\tNumMatchesPerSpec: " + this.numMatchesPerSpec + "\n");
        buf.append("\tMaxMissedCleavages: " + this.maxMissedCleavages + "\n");
        buf.append("\tMaxNumModsPerPeptide: " + this.maxNumMods + "\n");
        buf.append("\tChargeCarrierMass: " + this.chargeCarrierMass);

        if (Math.abs(this.chargeCarrierMass - PROTON) < 0.005) {
            buf.append(" (proton)\n");
        } else if (Math.abs(this.chargeCarrierMass - POTASSIUM_CHARGE_CARRIER_MASS) < 0.005) {
            buf.append(" (potassium)\n");
        } else if (Math.abs(this.chargeCarrierMass - SODIUM_CHARGE_CARRIER_MASS) < 0.005) {
            buf.append(" (sodium)\n");
        } else {
            buf.append(" (custom)\n");
        }

        buf.append("\tMSLevel: " + this.minMSLevel + "," + this.maxMSLevel + "\n");
        buf.append("\tMinNumPeaksPerSpectrum: " + this.minNumPeaksPerSpectrum + "\n");
        buf.append("\tNumIsoforms: " + this.maxNumVariantsPerPeptide + "\n");

        ArrayList<String> modificationsInUse = aaSet.getModificationsInUse();

        if (modificationsInUse.size() == 0) {
            buf.append("No static or dynamic post translational modifications are defined.\n");
        } else {
            buf.append("Post translational modifications in use:\n");
            for (String modInfo : modificationsInUse)
                buf.append("\t" + modInfo + "\n");
        }

        return buf.toString();
    }
}
