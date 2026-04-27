package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.misc.ProgressData;
import edu.ucsd.msjava.misc.ProgressReporter;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

public class ConcurrentMSGFPlus {
    private static final PrintStream NULL_PRINT_STREAM = new PrintStream(OutputStream.nullOutputStream());

    /** Per-task wall stats in milliseconds. {@code null} if the task didn't
     *  complete (interrupted). */
    public record TaskWallStats(int taskNum, long preprocessMs, long dbSearchMs,
                                long computeEvalueMs, long totalMs) {}

    public static class RunMSGFPlus implements Runnable, ProgressReporter {
        private final Supplier<ScoredSpectraMap> specScannerSupplier;
        private final CompactSuffixArray sa;
        SearchParams params;
        private final List<MSGFPlusMatch> resultList;
        private final int taskNum;
        private ProgressData progress;
        private ScoredSpectraMap specScanner;
        private DBScanner scanner;
        // Written once at end of run(); read by the main thread only after
        // executor.awaitTermination, which establishes happens-before.
        private TaskWallStats wallStats;

        public List<MSGFPlusMatch> getResults() {
            return resultList;
        }

        public int getResultCount() {
            return resultList.size();
        }

        public void drainResultsTo(List<MSGFPlusMatch> destination) {
            destination.addAll(resultList);
            resultList.clear();
        }

        public TaskWallStats getWallStats() {
            return wallStats;
        }

        @Override
        public void setProgressData(ProgressData data) {
            progress = data;
        }

        @Override
        public ProgressData getProgressData() {
            return progress;
        }

        public RunMSGFPlus(
                Supplier<ScoredSpectraMap> specScannerSupplier,
                CompactSuffixArray sa,
                SearchParams params,
                int taskNum
        ) {
            this.resultList = new ArrayList<>();
            this.specScannerSupplier = specScannerSupplier;
            this.sa = sa;
            this.params = params;
            this.taskNum = taskNum;
            progress = null;
        }

        @Override
        public void run() {
            long taskStartNs = System.nanoTime();
            long preprocessMs = 0, dbSearchMs = 0, computeEvalueMs = 0;
            if (progress == null) {
                progress = new ProgressData();
            }

            if (specScanner == null) {
                specScanner = specScannerSupplier.get();
                scanner = new DBScanner(
                        specScanner,
                        sa,
                        params.getEnzyme(),
                        params.getAASet(),
                        params.getNumMatchesPerSpec(),
                        params.getMinPeptideLength(),
                        params.getMaxPeptideLength(),
                        params.getMaxNumVariantsPerPeptide(),
                        params.getMinDeNovoScore(),
                        params.ignoreMetCleavage(),
                        params.getMaxMissedCleavages()
                );
            }

            PrintStream output;
            if (params.getVerbose()) {
                output = System.out;
            } else {
                output = NULL_PRINT_STREAM;
            }

            progress.stepRange(5.0);
            String threadName = Thread.currentThread().getName();
            output.println(threadName + ": Starting task " + taskNum);

            specScanner.setProgressObj(new ProgressData(progress));

            // Pre-process spectra
            long startTimePreprocess = System.currentTimeMillis();
            if (Thread.currentThread().isInterrupted()) {
                return;
            }

            if (specScanner.getPepMassSpecKeyMap().size() == 0)
                specScanner.makePepMassSpecKeyMap();

            output.println(threadName + ": Preprocessing spectra...");
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            specScanner.preProcessSpectra();
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            preprocessMs = System.currentTimeMillis() - startTimePreprocess;
            output.print(threadName + ": Preprocessing spectra finished ");
            output.format("(elapsed time: %.2f sec)\n", preprocessMs / 1000.0f);

            specScanner.getProgressObj().setParentProgressObj(null);
            progress.report(5.0);
            progress.stepRange(80.0);
            scanner.setProgressObj(new ProgressData(progress));

            long startTimeDbSearch = System.currentTimeMillis();

            // DB search
            output.println(threadName + ": Database search...");
            scanner.setThreadName(threadName);
            scanner.setPrintStream(output);

            int ntt = params.getNumTolerableTermini();
            if (params.getEnzyme() == null)
                ntt = 0;
            int nnet = 2 - ntt;
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            scanner.dbSearch(nnet);
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            dbSearchMs = System.currentTimeMillis() - startTimeDbSearch;
            output.print(threadName + ": Database search finished ");
            output.format("(elapsed time: %.2f sec)\n", dbSearchMs / 1000.0f);

            progress.stepRange(95.0);

            long startTimeComputeEvalue = System.currentTimeMillis();
            output.println(threadName + ": Computing spectral E-values...");
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            scanner.computeSpecEValue(false);
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            computeEvalueMs = System.currentTimeMillis() - startTimeComputeEvalue;
            output.print(threadName + ": Computing spectral E-values finished ");
            output.format("(elapsed time: %.2f sec)\n", computeEvalueMs / 1000.0f);

            scanner.getProgressObj().setParentProgressObj(null);
            progress.stepRange(100);

            if (Thread.currentThread().isInterrupted()) {
                return;
            }

            scanner.generateSpecIndexDBMatchMap();

            progress.report(30.0);

            if (params.outputAdditionalFeatures())
                scanner.addAdditionalFeatures();

            progress.report(60.0);

            scanner.addResultsToList(resultList);

            progress.report(100.0);
//			gen.addSpectrumIdentificationResults(scanner.getSpecIndexDBMatchMap());
            long totalMs = (System.nanoTime() - taskStartNs) / 1_000_000L;
            wallStats = new TaskWallStats(taskNum, preprocessMs, dbSearchMs, computeEvalueMs, totalMs);
            scanner = null;
            specScanner = null;
            output.println(threadName + ": Task " + taskNum + " completed.");
        }
    }
}
