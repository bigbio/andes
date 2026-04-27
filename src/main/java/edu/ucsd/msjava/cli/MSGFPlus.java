package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.fdr.ComputeFDR;
import edu.ucsd.msjava.misc.MSGFLogger;
import edu.ucsd.msjava.misc.RunManifestWriter;
import edu.ucsd.msjava.misc.ThreadPoolExecutorWithExceptions;
import edu.ucsd.msjava.msdbsearch.*;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.output.DirectPinWriter;
import edu.ucsd.msjava.output.DirectTSVWriter;
import edu.ucsd.msjava.mzml.StaxMzMLParser;
import edu.ucsd.msjava.sequences.Constants;
import picocli.CommandLine;
import picocli.CommandLine.ParameterException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;


public class MSGFPlus {
    public static final String VERSION = "Release (v2026.03.25)";
    public static final String RELEASE_DATE = "25 March 2026";

    public static final String DECOY_DB_EXTENSION = ".revCat.fasta";
    public static final String DEFAULT_DECOY_PROTEIN_PREFIX = "XXX";

    // Set this to true when debugging
    private static final boolean DISABLE_THREADING = false;

    /** Default numTasks-per-thread multiplier when {@code -tasks} is not
     *  passed. Users can override at the CLI via {@code -tasks -N}. */
    private static final int DEFAULT_TASKS_PER_THREAD = 3;
    private static final String USE_FORK_JOIN_PROPERTY = "msgfplus.useForkJoin";

    // Snapshot of the original CLI argv, captured in main() so that
    // RunManifestWriter can record it alongside the mzid without
    // threading argv through runMSGFPlus's many call sites.
    private static volatile String[] argvSnapshot = new String[0];

    public static void main(String argv[]) {
        long startTime = System.currentTimeMillis();
        argvSnapshot = argv == null ? new String[0] : argv.clone();

        MSGFPlusOptions opts = new MSGFPlusOptions();
        CommandLine cl = MSGFPlusOptions.commandLine(opts);

        if (argv.length == 0) {
            printToolInfo();
            cl.usage(System.out);
            return;
        }

        StaxMzMLParser.turnOffLogs();

        try {
            cl.parseArgs(argv);
        } catch (ParameterException e) {
            MSGFLogger.error(e.getMessage());
            System.out.println();
            cl.usage(System.out);
            System.exit(-1);
        }

        if (cl.isUsageHelpRequested()) {
            cl.usage(System.out);
            return;
        }
        if (cl.isVersionHelpRequested()) {
            System.out.println(VERSION);
            return;
        }

        // Propagate verbose flag to the shared logger before any downstream code logs.
        MSGFLogger.setVerbose(opts.effectiveVerbose() == 1);

        printToolInfo();
        printJVMInfo();

        String errorMessage = null;
        try {
            errorMessage = runMSGFPlus(opts);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }

        if (errorMessage != null) {
            MSGFLogger.error(errorMessage);
            System.out.println();
            System.exit(-1);
        } else
            MSGFLogger.info("MS-GF+ complete (total elapsed time: %.2f sec)", (System.currentTimeMillis() - startTime) / (float) 1000);
    }

    private static void printToolInfo() {
        System.out.println("MS-GF+ " + VERSION + " (" + RELEASE_DATE + ")");
    }

    private static void printJVMInfo() {
        System.out.println("Java " + System.getProperty("java.version") + " (" + System.getProperty("java.vendor") + ")");
        System.out.println(System.getProperty("os.name") + " (" + System.getProperty("os.arch") + ", version " + System.getProperty("os.version") + ")");
    }

    public static String runMSGFPlus(MSGFPlusOptions opts) {
        SearchParams params = new SearchParams();
        String errorMessage = params.parse(opts);

        if (errorMessage != null) {
            return errorMessage;
        }

        List<DBSearchIOFiles> ioList = params.getDBSearchIOList();
        boolean multiFiles = false;
        if (ioList.size() >= 2) {
            MSGFLogger.info("Processing " + ioList.size() + " spectra");
            for (DBSearchIOFiles ioFiles : ioList) {
                MSGFLogger.debug("\t" + ioFiles.getSpecFile().getName());
            }
            multiFiles = true;
        }

        int ioIndex = -1;
        for (DBSearchIOFiles ioFiles : ioList) {
            ++ioIndex;
            File specFile = ioFiles.getSpecFile();
            SpecFileFormat specFormat = ioFiles.getSpecFileFormat();
            File outputFile = ioFiles.getOutputFile();

            if (multiFiles) {
                if (!outputFile.exists()) {
                    MSGFLogger.info("\nProcessing " + specFile.getPath());
                    MSGFLogger.debug("Writing results to " + outputFile.getPath());
                    String errMsg = runMSGFPlus(ioIndex, specFormat, outputFile, params);
                    if (errMsg != null) {
                        return errMsg;
                    }
                    RunManifestWriter.write(ioFiles, params, VERSION, argvSnapshot);
                } else {
                    MSGFLogger.info("\nIgnoring " + specFile.getPath());
                    MSGFLogger.debug("Output file " + outputFile.getPath() + " exists.");
                }
            } else {
                String errMsg = runMSGFPlus(ioIndex, specFormat, outputFile, params);
                if (errMsg != null) {
                    return errMsg;
                }
                RunManifestWriter.write(ioFiles, params, VERSION, argvSnapshot);
            }
        }

        return null;
    }

    private static String runMSGFPlus(int ioIndex, SpecFileFormat specFormat, File outputFile, SearchParams params) {
        long startTime = System.currentTimeMillis();

        // Verify that the output directory exists and can be written to
        File outputDirectory = outputFile.getParentFile();
        if (outputDirectory != null) {
            if (!outputDirectory.exists()) {
                System.out.println("Creating directory " + outputDirectory.getPath());
                boolean success = outputDirectory.mkdirs();
                if (!success) {
                    return "Unable to create the missing directory: " + outputDirectory.getPath();
                }
            } else if (!outputDirectory.isDirectory()) {
                return "Invalid output file path (file path instead of directory path?): " + outputDirectory.getPath();
            }

            // An easy way to test for write access is outputDirectory.canWrite()
            // However, on Windows this is not always accurate
            // Thus, create a temporary file then delete it
            try {
                File testFile = File.createTempFile("MSGFPlus", ".tmp", outputDirectory);
                testFile.delete();
            } catch (java.io.IOException e) {
                return "Cannot create files in the output directory: " + e.getMessage();
            } catch (SecurityException e) {
                return "Cannot create files in the output directory; permission denied for: " + outputDirectory.getPath();
            }
        }

        // DB file
        File databaseFile = params.getDatabaseFile();

        if (databaseFile == null) {
            return "Database file is not defined; use -d at the command line or DatabaseFile in a config file";
        }

        if (!databaseFile.exists()) {
            return "Database file not found: " + databaseFile.getPath();
        }

        // Precursor mass tolerance
        Tolerance leftPrecursorMassTolerance = params.getLeftPrecursorMassTolerance();
        Tolerance rightPrecursorMassTolerance = params.getRightPrecursorMassTolerance();

        int minIsotopeError = params.getMinIsotopeError();    // inclusive
        int maxIsotopeError = params.getMaxIsotopeError();    // inclusive

        Enzyme enzyme = params.getEnzyme();

        ActivationMethod activationMethod = params.getActivationMethod();
        InstrumentType instType = params.getInstType();
        Protocol protocol = params.getProtocol();

        AminoAcidSet aaSet = params.getAASet();

        int startSpecIndex = params.getStartSpecIndex();
        int endSpecIndex = params.getEndSpecIndex();

        boolean useTDA = params.useTDA();

        int minCharge = params.getMinCharge();
        int maxCharge = params.getMaxCharge();

        int numThreads = params.getNumThreads();
        boolean doNotUseEdgeScore = params.doNotUseEdgeScore();
        boolean allowDenseCentroidedPeaks = params.getAllowDenseCentroidedPeaks();

        int minNumPeaksPerSpectrum = params.getMinNumPeaksPerSpectrum();
        if (minNumPeaksPerSpectrum == -1)    // not specified
        {
            if (instType == InstrumentType.TOF)
                minNumPeaksPerSpectrum = Constants.MIN_NUM_PEAKS_PER_SPECTRUM_TOF;
            else
                minNumPeaksPerSpectrum = Constants.MIN_NUM_PEAKS_PER_SPECTRUM;
        }

        String decoyProteinPrefix = params.getDecoyProteinPrefix();

        System.out.println("Loading database files...");

        File dbIndexDir = params.getDBIndexDir();
        if (dbIndexDir != null) {

            File newDBFile = new File(Paths.get(dbIndexDir.getPath(), databaseFile.getName()).toString());
            if (!useTDA) {
                if (!newDBFile.exists()) {
                    System.out.println("Creating " + newDBFile.getPath() + ".");
                    ReverseDB.copyDB(databaseFile.getPath(), newDBFile.getPath());
                }
            }
            databaseFile = newDBFile;
        }

        if (useTDA) {
            String dbFileName = databaseFile.getName();
            String concatDBFileName = dbFileName.substring(0, dbFileName.lastIndexOf('.')) + DECOY_DB_EXTENSION;

            String concatDBFilePath = Paths.get(databaseFile.getAbsoluteFile().getParent(), concatDBFileName).toString();
            File concatTargetDecoyDBFile = new File(concatDBFilePath);

            if (!concatTargetDecoyDBFile.exists()) {
                System.out.println("Creating " + concatTargetDecoyDBFile.getPath() + ".");
                if (ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true, decoyProteinPrefix) == false) {
                    return "Cannot create a decoy database file!";
                }
            }
            databaseFile = concatTargetDecoyDBFile;
        }

        DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
        aaSet.registerEnzyme(enzyme);

        CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getPath());
        fastaSequence.setDecoyProteinPrefix(decoyProteinPrefix);

        if (useTDA) {
            float ratioUniqueProteins = fastaSequence.getRatioUniqueProteins();
            if (ratioUniqueProteins < 0.5f) {
                fastaSequence.printTooManyDuplicateSequencesMessage(databaseFile.getName(), "MS-GF+");
                System.exit(-1);
            }

            float fractionDecoyProteins = fastaSequence.getFractionDecoyProteins();
            if (fractionDecoyProteins < 0.4f || fractionDecoyProteins > 0.6f) {
                MSGFLogger.error("Error while reading: " + databaseFile.getName() + " (fraction of decoy proteins: " + fractionDecoyProteins + ")");
                MSGFLogger.error("Delete " + databaseFile.getName() + " and run MS-GF+ again.");
                MSGFLogger.error("Decoy protein names should start with " + fastaSequence.getDecoyProteinPrefix());
                System.exit(-1);
            }
        }

        CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, params.getMaxPeptideLength());
        System.out.print("Loading database finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - startTime) / 1000);

        System.out.println("Reading spectra...");

        File specFile = params.getDBSearchIOList().get(ioIndex).getSpecFile();

        // Show a message of the form "Opening mzML file QC_Mam_19_01_PNNL_10_06Jan21_Arwen_WBEH-20-12-01.mzML"
        System.out.printf("Opening %s %s\n", specFormat.getPSIName(), specFile.getName());

        SpectraAccessor specAcc = new SpectraAccessor(specFile, specFormat);
        int minMSLevel = params.getMinMSLevel();
        int maxMSLevel = params.getMaxMSLevel();
        specAcc.setMSLevelRange(minMSLevel, maxMSLevel);

        if (specAcc.getSpecMap() == null || specAcc.getSpecItr() == null)
            return "Error while parsing spectrum file: " + specFile.getPath();

        ArrayList<SpecKey> specKeyList = SpecKey.getSpecKeyList(specAcc,
                startSpecIndex, endSpecIndex, minCharge, maxCharge, activationMethod, minNumPeaksPerSpectrum, allowDenseCentroidedPeaks,
                minMSLevel, maxMSLevel);

        int specSize = specKeyList.size();
        if (specSize == 0)
            return specFile.getPath() + " does not have any valid spectra";

        System.out.print("Reading spectra finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - startTime) / 1000);

        if (numThreads <= 0)
            numThreads = 1;

        // Minimum spectra/task(or thread) floor for efficiency; going smaller slows down processing.
        // Configurable via -minSpectraPerThread for users on many-core hosts with small inputs (see #52).
        int spectraPerTaskMinimum = params.getMinSpectraPerThread();
        int maxThreads = Math.max(1, Math.round((float) specSize / spectraPerTaskMinimum));
        if (maxThreads < numThreads) {
            if (maxThreads == 1) {
                System.out.println("Note: under " + spectraPerTaskMinimum + " spectra; using 1 thread instead of " + numThreads);
            } else {
                System.out.println("Note: " + spectraPerTaskMinimum + " spectra per thread minimum; using " + maxThreads + " threads instead of " + numThreads);
            }

            numThreads = maxThreads;
        }

        System.out.println("Using " + numThreads + (numThreads == 1 ? " thread." : " threads."));

        // Print out parameters
        System.out.println("Search Parameters:");
        System.out.println(params.toString());

        SpecDataType specDataType = new SpecDataType(activationMethod, instType, enzyme, protocol);

        // Achievement B — two-pass precursor mass calibration (P2-cal).
        // Runs a sampled pre-pass over the current file's SpecKeys to learn
        // a per-file ppm shift, then stores it on DBSearchIOFiles so every
        // task-local ScoredSpectraMap picks it up. OFF mode is a strict
        // no-op: we skip the pre-pass entirely and never call the setter,
        // so DBSearchIOFiles.precursorMassShiftPpm stays at its 0.0 default
        // and ScoredSpectraMap.applyShift() takes its exact-zero fast path.
        DBSearchIOFiles currentIoFiles = params.getDBSearchIOList().get(ioIndex);
        if (params.getPrecursorCalMode() != SearchParams.PrecursorCalMode.OFF) {
            long calStart = System.currentTimeMillis();
            MassCalibrator calibrator = new MassCalibrator(
                    specAcc,
                    sa,
                    aaSet,
                    params,
                    specKeyList,
                    leftPrecursorMassTolerance,
                    rightPrecursorMassTolerance,
                    minIsotopeError,
                    maxIsotopeError,
                    specDataType);
            double shiftPpm = calibrator.learnPrecursorShiftPpm(ioIndex);
            boolean applyLearnedShift = shiftPpm != 0.0
                    || params.getPrecursorCalMode() == SearchParams.PrecursorCalMode.ON;
            if (applyLearnedShift) {
                currentIoFiles.setPrecursorMassShiftPpm(shiftPpm);
                System.out.printf("Precursor mass shift learned: %.3f ppm (elapsed: %.2f sec)%n",
                        shiftPpm, (System.currentTimeMillis() - calStart) / 1000.0);
            } else {
                System.out.printf("Precursor mass calibration skipped (insufficient confident PSMs; elapsed: %.2f sec)%n",
                        (System.currentTimeMillis() - calStart) / 1000.0);
            }
        }
        double precursorMassShiftPpm = currentIoFiles.getPrecursorMassShiftPpm();

        List<MSGFPlusMatch> resultList;

        int toIndexGlobal = specSize;
        while (toIndexGlobal < specSize) {
            SpecKey lastSpecKey = specKeyList.get(toIndexGlobal - 1);
            SpecKey nextSpecKey = specKeyList.get(toIndexGlobal);

            if (lastSpecKey.getSpecIndex() == nextSpecKey.getSpecIndex())
                toIndexGlobal++;
            else
                break;
        }

        System.out.println("Spectrum 0-" + (toIndexGlobal - 1) + " (total: " + specSize + ")");

        boolean useForkJoin = Boolean.getBoolean(USE_FORK_JOIN_PROPERTY);

        ThreadPoolExecutorWithExceptions executor =
                useForkJoin ? null : ThreadPoolExecutorWithExceptions.newFixedThreadPool(numThreads);
        if (executor != null) executor.setTaskName("Search");
        ForkJoinPool fjp = useForkJoin ? new ForkJoinPool(numThreads) : null;
        List<Future<?>> fjpFutures = useForkJoin ? new ArrayList<>() : null;

        int numTasks = Math.min(numThreads * DEFAULT_TASKS_PER_THREAD, Math.round((float) specSize / spectraPerTaskMinimum));
        if (numThreads <= 1) {
            numTasks = 1;
        }

        if (params.getNumTasks() != 0) {
            numTasks = params.getNumTasks();
            if (numTasks < 0) {
                numTasks = numThreads * (numTasks * -1);
            }
            if (numTasks < numThreads) {
                System.out.println("Changing specified tasks from " + numTasks + " to " + numThreads + " to provide the minimum of one task per thread.");
                numTasks = numThreads;
            }
        }
        if (numTasks > 1) {
            System.out.println("Splitting work into " + numTasks + " tasks.");
        } else {
            System.out.println("Searching using a single task.");
        }

        // Partition specKeyList
        int size = toIndexGlobal;
        int residue = size % numTasks;

        int[] startIndex = new int[numTasks];
        int[] endIndex = new int[numTasks];

        int subListSize = size / numTasks;
        for (int i = 0; i < numTasks; i++) {
            startIndex[i] = i > 0 ? endIndex[i - 1] : 0;
            endIndex[i] = startIndex[i] + subListSize + (i < residue ? 1 : 0);

            subListSize = size / numTasks;
            while (endIndex[i] < specKeyList.size()) {
                SpecKey lastSpecKey = specKeyList.get(endIndex[i] - 1);
                SpecKey nextSpecKey = specKeyList.get(endIndex[i]);

                if (lastSpecKey.getSpecIndex() == nextSpecKey.getSpecIndex()) {
                    ++endIndex[i];
                    --subListSize;
                } else
                    break;
            }
        }

        List<ConcurrentMSGFPlus.RunMSGFPlus> submittedTasks = new ArrayList<>(numTasks);

        try {
            for (int i = 0; i < numTasks; i++) {
                final int taskStartIndex = startIndex[i];
                final int taskEndIndex = endIndex[i];
                final boolean storeRankScorer = params.outputAdditionalFeatures();
                final int taskNum = i + 1;

                // Defer ScoredSpectraMap construction to the worker so the
                // per-task spectrum heap isn't queued up front.
                ConcurrentMSGFPlus.RunMSGFPlus msgfplusExecutor = new ConcurrentMSGFPlus.RunMSGFPlus(
                        () -> {
                            ScoredSpectraMap specScanner = new ScoredSpectraMap(
                                    specAcc,
                                    specKeyList.subList(taskStartIndex, taskEndIndex),
                                    leftPrecursorMassTolerance,
                                    rightPrecursorMassTolerance,
                                    minIsotopeError,
                                    maxIsotopeError,
                                    specDataType,
                                    storeRankScorer,
                                    false,
                                    precursorMassShiftPpm
                            );
                            if (doNotUseEdgeScore)
                                specScanner.turnOffEdgeScoring();
                            return specScanner;
                        },
                        sa,
                        params,
                        taskNum
                );

                submittedTasks.add(msgfplusExecutor);

                if (DISABLE_THREADING) {
                    msgfplusExecutor.run();
                } else if (useForkJoin) {
                    fjpFutures.add(fjp.submit(msgfplusExecutor));
                } else {
                    executor.execute(msgfplusExecutor);
                }

            }

            if (useForkJoin) {
                fjp.shutdown();
                try {
                    fjp.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, e.getMessage(), e);
                }
                for (Future<?> f : fjpFutures) {
                    try { f.get(); }
                    catch (java.util.concurrent.ExecutionException ex) {
                        Throwable cause = ex.getCause();
                        Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, cause.getMessage(), cause);
                        fjp.shutdownNow();
                        return "Search failed: " + cause.getMessage();
                    }
                    catch (InterruptedException ex) { Thread.currentThread().interrupt(); }
                }
            } else {
                executor.outputProgressReport();
                executor.shutdown();
                try {
                    executor.awaitTerminationWithExceptions(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                } catch (InterruptedException e) {
                    if (!executor.HasThrownData()) {
                        e.printStackTrace();
                        Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, e.getMessage(), e);
                    }
                }
                executor.outputProgressReport();
            }

            // awaitTermination above establishes happens-before on every
            // task's writes (JLS §17.4.5), so the per-task ArrayLists can
            // be drained single-threaded with no synchronization.
            int totalResults = 0;
            for (ConcurrentMSGFPlus.RunMSGFPlus t : submittedTasks) {
                totalResults += t.getResultCount();
            }
            resultList = new ArrayList<>(totalResults);
            for (ConcurrentMSGFPlus.RunMSGFPlus t : submittedTasks) {
                t.drainResultsTo(resultList);
            }

            if (numTasks > 1) {
                printTaskWallSummary(submittedTasks);
            }
            submittedTasks.clear();

        } catch (OutOfMemoryError ex) {
            ex.printStackTrace();
            Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, null, ex);
            shutdownPoolNow(executor, fjp);
            int taskMult = numTasks / numThreads;
            return "Task terminated; results incomplete. Please run again with a greater amount of memory, using \"-Xmx4G\", for example.\n" +
                    "\tYou can also use less memory by increasing the number of tasks used for the search, at the cost of more time.\n" +
                    "\tTry doubling the number used for this search with \"-tasks -" + (taskMult * 2) + "\" or \"-tasks " + (numTasks * 2) + "\".";
        } catch (Exception ex) {
            ex.printStackTrace();
            Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, null, ex);
            shutdownPoolNow(executor, fjp);
            return "Task terminated; results incomplete. Please run again.";
        } catch (Throwable ex) {
            ex.printStackTrace();
            Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, null, ex);
            shutdownPoolNow(executor, fjp);
            return "Task terminated; results incomplete. Please run again.";
        }

        long qValueStartTime = System.currentTimeMillis();

        if (params.useTDA()) {
            // Compute Q-values
            System.out.println("Computing q-values...");
            ComputeFDR.addQValues(resultList, sa, false, decoyProteinPrefix);
            System.out.print("Computing q-values finished ");
            System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - qValueStartTime) / 1000);
        }

        // Sort by spectral E-values then write to disk

        long saveResultsStartTime = System.currentTimeMillis();

        System.out.println("Writing results...");
        Collections.sort(resultList);

        if (params.writeTsv()) {
            DirectTSVWriter tsvWriter = new DirectTSVWriter(params, aaSet, sa, specAcc, ioIndex);
            try {
                tsvWriter.writeResults(resultList, outputFile);
            } catch (IOException e) {
                return "Error writing TSV output: " + e.getMessage();
            }
            System.out.println("TSV file: " + outputFile.getPath());
        }

        if (!params.writeTsv()) {
            DirectPinWriter pinWriter = new DirectPinWriter(params, aaSet, sa, specAcc, ioIndex);
            try {
                pinWriter.writeResults(resultList, outputFile);
            } catch (IOException e) {
                return "Error writing pin output: " + e.getMessage();
            }
            System.out.println("PIN file: " + outputFile.getPath());
        }

        System.out.print("Writing results finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - saveResultsStartTime) / 1000);
        return null;
    }

    private static void shutdownPoolNow(ThreadPoolExecutorWithExceptions executor, ForkJoinPool fjp) {
        if (executor != null) executor.shutdownNow();
        else if (fjp != null) fjp.shutdownNow();
    }

    /**
     * One-line wall-time summary across completed tasks. tail_gap (max -
     * median) is the load-balance signal; high values point at uneven
     * SpecKey distribution and motivate raising the {@code -tasks -N} multiplier.
     */
    private static void printTaskWallSummary(List<ConcurrentMSGFPlus.RunMSGFPlus> tasks) {
        List<Long> walls = new ArrayList<>(tasks.size());
        for (ConcurrentMSGFPlus.RunMSGFPlus t : tasks) {
            ConcurrentMSGFPlus.TaskWallStats s = t.getWallStats();
            if (s != null) walls.add(s.totalMs());
        }
        if (walls.isEmpty()) return;
        Collections.sort(walls);
        long min = walls.get(0);
        long max = walls.get(walls.size() - 1);
        long median = walls.get(walls.size() / 2);
        long p95 = walls.get(Math.min(walls.size() - 1, (int) Math.ceil(walls.size() * 0.95) - 1));
        long sum = 0L;
        for (long w : walls) sum += w;
        System.out.format(
                "Task wall summary (n=%d): min=%.1fs median=%.1fs p95=%.1fs max=%.1fs total=%.1fs tail_gap=%.1fs (%.0f%% of median)%n",
                walls.size(), min / 1000.0, median / 1000.0, p95 / 1000.0, max / 1000.0,
                sum / 1000.0, (max - median) / 1000.0,
                median > 0 ? 100.0 * (max - median) / median : 0.0);
    }
}
