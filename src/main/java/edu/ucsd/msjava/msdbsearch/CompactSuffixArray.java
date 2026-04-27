package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.sequences.Constants;
import edu.ucsd.msjava.suffixarray.ByteSequence;
import edu.ucsd.msjava.suffixarray.SuffixFactory;
import it.unimi.dsi.fastutil.ints.IntArrays;

import java.io.*;
import java.nio.file.Files;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * SuffixArray class for fast exact matching.
 *
 * @author Sangtae Kim
 */
public class CompactSuffixArray {

    public static final int COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID = 8294;

    /***** CONSTANTS *****/
    /**
     * Default extension of a suffix array file.
     */
    protected static final String EXTENSION_INDICES = ".csarr";

    /**
     * Default extension of a neighboring longest common prefix file
     */
    protected static final String EXTENSION_NLCPS = ".cnlcp";

    /**
     * Size of the bucket for the suffix array creation
     */
    protected static final int BUCKET_SIZE = 5;

    /**
     * Size of an int primitive type in bytes
     */
    protected static final int INT_BYTE_SIZE = Integer.SIZE / Byte.SIZE;

    /***** MEMBERS *****/
    /**
     * Tracks indices of the sorted suffixes
     */
    private final File indexFile;

    /**
     * Tracks precomputed LCPs (longest common prefixes) of neighboring suffixes
     */
    private final File nlcpFile;

    /**
     * Sequence representing all the suffixes
     */
    private CompactFastaSequence sequence;

    /**
     * Class that generates suffixes from the given adapter
     */
    private SuffixFactory factory;

    /**
     * Number of suffixes in this suffix array
     */
    private int size;

    /**
     * Maximum peptide length
     */
    private int maxPeptideLength;

    /**
     * number of distinct peptides
     */
    private int[] numDistinctPeptides;


    /**
     * Constructor that attempts to read the suffix array from the provided file.
     *
     * @param sequence the sequence object.
     */
    public CompactSuffixArray(CompactFastaSequence sequence) {
        // infer the suffix array file from the sequence.
        this.sequence = sequence;
        this.size = (int) sequence.getSize();
        this.factory = new SuffixFactory(sequence);
        indexFile = new File(sequence.getBaseFilepath() + EXTENSION_INDICES);
        nlcpFile = new File(sequence.getBaseFilepath() + EXTENSION_NLCPS);

        // create the file if it doesn't exist or the metadata differs
        if (!indexFile.exists() || !nlcpFile.exists() || !isCompactSuffixArrayValid(sequence.getLastModified())) {
            createSuffixArrayFiles(sequence, indexFile, nlcpFile);
        }

        // check the ids of indexFile and nlcpFile
        int id = checkID();

        // check that the files are consistent
        if (id != sequence.getId()) {
            System.err.println("Suffix array files are not consistent: " + indexFile + ", " + nlcpFile + " (" + id + "!=" + sequence.getId() + ")");
            System.err.println("Please recreate the suffix array file by deleting the .canno, .cseq, and .csarr files.");
            System.exit(-1);
        }
    }

    /**
     * Constructor that attempts to read the suffix array from the provided file.
     *
     * @param sequence the sequence object.
     */
    public CompactSuffixArray(CompactFastaSequence sequence, int maxPeptideLength) {
        this(sequence);
        this.maxPeptideLength = maxPeptideLength;
        computeNumDistinctPeptides();
    }

    public File getIndexFile() {
        return this.indexFile;
    }

    public File getNeighboringLcpFile() {
        return this.nlcpFile;
    }

    public CompactFastaSequence getSequence() {
        return sequence;
    }

    public int getSize() {
        return size;
    }

    public int getNumDistinctPeptides(int length) {
        // no boundary check
        return numDistinctPeptides[length];
    }

    public String getAnnotation(long index) {
        return sequence.getAnnotation(index);
    }

    private boolean isCompactSuffixArrayValid(long lastModified) {
        File[] files = {indexFile, nlcpFile};

        for (File f : files) {
            try {
                RandomAccessFile raf = new RandomAccessFile(f, "r");
                raf.seek(raf.length() - Integer.SIZE / 8 - Long.SIZE / 8);
                long lastModifiedRecorded = raf.readLong();
                int id = raf.readInt();
                raf.close();

                if (!NearlyEqualFileTimes(lastModifiedRecorded, lastModified)) {
                    Date suffixArrayModificationTime = new Date(lastModifiedRecorded);
                    Date fastaFileModificationTime = new Date(lastModified);
                    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss", Locale.US);

                    System.out.println("Re-creating suffix array files since the cached LastModified time is not within 2 seconds " +
                            "of the LastModified time of the sequence file:\n" +
                            " Time cached in " + f.getName() + " is " + lastModifiedRecorded +
                            " (" + dateFormat.format(suffixArrayModificationTime) + ")" +
                            " while the sequence file has " + dateFormat.format(fastaFileModificationTime));
                    return false;
                }

                if (id != COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID) {
                    System.out.println("Re-creating suffix array files since " + f.getName() +
                            " has file format ID " + id + " instead of " + COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID);
                    return false;
                }

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return true;
    }

    // TODO: this method has a bug (according to Sangtae in 2011)
    // The only evident bug is no checks for reading past the end of a file
    private void computeNumDistinctPeptides() {
        boolean[] isValidResidue = new boolean[128];
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        for (AminoAcid aa : aaSet)
            isValidResidue[aa.getResidue()] = true;

        // This array keeps track of the number of possible peptides of each length
        numDistinctPeptides = new int[maxPeptideLength + 2];
        try {
            File indexFile = getIndexFile();
            System.out.printf("Counting number of distinct peptides in %s using %s\n", indexFile.getName(), nlcpFile.getName());

            DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(indexFile)));
            indices.skip(CompactSuffixArray.INT_BYTE_SIZE * 2);    // skip size and id

            DataInputStream neighboringLcps = new DataInputStream(new BufferedInputStream(new FileInputStream(nlcpFile)));
            int size = neighboringLcps.readInt();
            neighboringLcps.readInt();    // skip id

            long lastStatusTime = System.currentTimeMillis();

            for (int i = 0; i < size; i++) {
                // print progress
                if (i % 100000 == 0 && System.currentTimeMillis() - lastStatusTime > 2000) {
                    lastStatusTime = System.currentTimeMillis();
                    System.out.printf("Counting distinct peptides: %.2f%% complete.\n", i * 100.0 / size);
                }

                int index = indices.readInt();
                byte lcp = neighboringLcps.readByte();
                int idx = sequence.getCharAt(index);
                if (isValidResidue[idx] == false)
                    continue;

                for (int l = lcp + 1; l < numDistinctPeptides.length; l++) {
                    numDistinctPeptides[l]++;
                }
            }
            neighboringLcps.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    /**
     * Helper method that initializes the suffixArray object from the file.
     * Initializes indices, leftMiddleLcps, middleRightLcps and neighboringLcps.
     *
     * @return returns the id of this file for consistency check.
     */
    private int checkID() {
        //		System.out.println("SAForMSGFDB Reading " + suffixFile);
        try {
            DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(indexFile)));
            // read the first integer which encodes for the size of the file
            int sizeIndexFile = indices.readInt();
            // the second integer is the id
            int idIndexFile = indices.readInt();

            DataInputStream neighboringLcps = new DataInputStream(new BufferedInputStream(new FileInputStream(nlcpFile)));
            int sizeNLcp = neighboringLcps.readInt();
            int idNLcp = neighboringLcps.readInt();

            indices.close();
            neighboringLcps.close();

            if (sizeIndexFile == sizeNLcp && idIndexFile == idNLcp)
                return idIndexFile;
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        return 0;
    }

    /** Sysprop overriding the number of threads used during the sort+LCP phase. */
    static final String SA_BUILD_THREADS_PROPERTY = "msgfplus.buildsa.threads";

    /** Cap on default thread count: higher values give diminishing returns and thrash IO. */
    private static final int MAX_DEFAULT_SA_BUILD_THREADS = 8;

    /**
     * Build the suffix-array index files. Two-phase radix-then-sort: each suffix
     * is hashed by its first {@link #BUCKET_SIZE} residues into a bucket, then
     * sorted lexicographically from offset {@code BUCKET_SIZE} onward. The
     * sort+LCP phase is parallelised across contiguous bucket-id ranges; the
     * write step is single-threaded to preserve on-disk ordering.
     */
    private void createSuffixArrayFiles(CompactFastaSequence sequence, File indexFile, File nlcpFile) {
        System.out.println("Creating the suffix array indexed file... Size: " + sequence.getSize());

        // the size of the alphabet to make the hashes
        int hashBase = sequence.getAlphabetSize();
        System.out.println("AlphabetSize: " + sequence.getAlphabetSize());
        if (hashBase > 30) {
            System.err.println("Suffix array construction failure: alphabet size is too large: " + sequence.getAlphabetSize());
            System.exit(-1);
        }

        // this number is to efficiently calculate the next hash
        int denominator = 1;
        for (int i = 0; i < BUCKET_SIZE - 1; i++)
            denominator *= hashBase;

        // the number of buckets  required to encode for all hashes
        int numBuckets = denominator * hashBase;

        // initial value of the hash
        int currentHash = 0;
        for (int i = 0; i < BUCKET_SIZE - 1; i++) {
            currentHash = currentHash * hashBase + sequence.getByteAt(i);
        }

        // the main array that stores the sorted buckets of suffixes
        Bucket[] bucketSuffixes = new Bucket[numBuckets];

        long lastStatusTime = System.currentTimeMillis();
        int numResiduesInSequence = (int) sequence.getSize();

        // main loop for putting suffixes into the buckets
        for (int i = BUCKET_SIZE - 1, j = 0; j < numResiduesInSequence; i++, j++) {
            // print progress
            if (j % 100000 == 0 && System.currentTimeMillis() - lastStatusTime > 2000) {
                lastStatusTime =  System.currentTimeMillis();
                System.out.printf("Suffix creation: %.2f%% complete.\n", j * 100.0 / numResiduesInSequence);
            }

            // quick wait to derive the next hash, since we are reading the sequence in order
            byte b = Constants.TERMINATOR;
            if (i < numResiduesInSequence)
                b = sequence.getByteAt(i);

            currentHash = (currentHash % denominator) * hashBase + b;

            // first bucket at this position
            if (bucketSuffixes[currentHash] == null) bucketSuffixes[currentHash] = new Bucket();

            // insert suffix
            bucketSuffixes[currentHash].add(j);
        }

        try {
            DataOutputStream indexOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(indexFile)));
            DataOutputStream nlcpOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(nlcpFile)));
            indexOut.writeInt(numResiduesInSequence);
            indexOut.writeInt(sequence.getId());
            nlcpOut.writeInt(numResiduesInSequence);
            nlcpOut.writeInt(sequence.getId());

            System.out.println("Sorting suffixes... Size: " + bucketSuffixes.length);
            sortAndWriteBuckets(sequence, bucketSuffixes, indexFile, indexOut, nlcpOut);

            long lastModified = sequence.getLastModified();
            indexOut.writeLong(lastModified);
            indexOut.writeInt(CompactSuffixArray.COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID);
            indexOut.flush();
            indexOut.close();

            nlcpOut.writeLong(lastModified);
            nlcpOut.writeInt(CompactSuffixArray.COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID);
            nlcpOut.flush();
            nlcpOut.close();

            // Do not compute Llcps and Rlcps
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return;
    }

    /**
     * Sort + LCP compute phase. Parallelises across contiguous bucket-id
     * ranges; each worker streams its sorted indices + intra-range LCPs into
     * per-range temp files. The merge step fixes up the cross-range boundary
     * LCP byte and streams the temp files into the final output sequentially
     * (writing single-threaded preserves on-disk ordering). Temp files are
     * deleted in the {@code finally} block, with {@link File#deleteOnExit} as
     * a fallback for hard crashes.
     */
    private static void sortAndWriteBuckets(CompactFastaSequence sequence,
                                            Bucket[] bucketSuffixes,
                                            File indexFile,
                                            DataOutputStream indexOut,
                                            DataOutputStream nlcpOut) throws IOException {
        int numThreads = resolveSortThreads();
        int[][] ranges = partitionBucketIds(bucketSuffixes, numThreads);

        if (ranges.length == 1) {
            writeBucketsDirect(sequence, bucketSuffixes, ranges[0][0], ranges[0][1], indexOut, nlcpOut);
            return;
        }

        File parentDir = indexFile.getAbsoluteFile().getParentFile();
        if (parentDir == null) parentDir = new File(".");
        String tempBasename = indexFile.getName() + ".buildsa-tmp." + ProcessHandle.current().pid() + "." + System.nanoTime();

        List<RangeMetadata> rangeMetadatas = new ArrayList<>(ranges.length);
        try {
            ExecutorService pool = Executors.newFixedThreadPool(ranges.length, r -> {
                Thread t = new Thread(r, "buildsa-sort");
                t.setDaemon(true);
                return t;
            });
            try {
                List<Future<RangeMetadata>> futures = new ArrayList<>(ranges.length);
                for (int idx = 0; idx < ranges.length; idx++) {
                    final int from = ranges[idx][0];
                    final int to = ranges[idx][1];
                    final File tempIndices = new File(parentDir, tempBasename + ".indices." + idx);
                    final File tempLcps = new File(parentDir, tempBasename + ".lcps." + idx);
                    tempIndices.deleteOnExit();
                    tempLcps.deleteOnExit();
                    futures.add(pool.submit(() -> processBucketRangeToTempFiles(
                            sequence, bucketSuffixes, from, to, tempIndices, tempLcps)));
                }
                for (Future<RangeMetadata> f : futures) {
                    rangeMetadatas.add(f.get());
                }
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new IOException("Interrupted while building suffix array", e);
            } catch (ExecutionException e) {
                Throwable cause = e.getCause();
                if (cause instanceof RuntimeException) throw (RuntimeException) cause;
                if (cause instanceof IOException) throw (IOException) cause;
                throw new IOException("Suffix array sort worker failed", cause != null ? cause : e);
            } finally {
                pool.shutdown();
            }

            int prevRangeLastBucketFirst = -1;
            for (RangeMetadata md : rangeMetadatas) {
                if (md.numEntries == 0) continue;
                mergeRangeIntoOutput(sequence, md, prevRangeLastBucketFirst, indexOut, nlcpOut);
                prevRangeLastBucketFirst = md.lastBucketFirstSuffix;
            }
        } finally {
            for (RangeMetadata md : rangeMetadatas) {
                deleteQuietly(md.tempIndicesFile);
                deleteQuietly(md.tempLcpsFile);
            }
            // Sweep debris from workers that died before returning a RangeMetadata.
            File[] orphans = parentDir.listFiles((dir, name) -> name.startsWith(tempBasename));
            if (orphans != null) {
                for (File f : orphans) deleteQuietly(f);
            }
        }
    }

    private static void deleteQuietly(File f) {
        if (f == null) return;
        try { Files.deleteIfExists(f.toPath()); } catch (IOException ignored) { }
    }

    /**
     * Stream one range's temp files into the final output. The first LCP byte
     * is rewritten against {@code prevRangeLastBucketFirst} to bridge the
     * cross-range boundary; for the globally-first range
     * {@code prevRangeLastBucketFirst} is -1 and the placeholder 0 written by
     * the worker passes through.
     */
    private static void mergeRangeIntoOutput(CompactFastaSequence sequence,
                                             RangeMetadata md,
                                             int prevRangeLastBucketFirst,
                                             DataOutputStream indexOut,
                                             DataOutputStream nlcpOut) throws IOException {
        try (DataInputStream idxIn = new DataInputStream(new BufferedInputStream(new FileInputStream(md.tempIndicesFile)));
             DataInputStream lcpIn = new DataInputStream(new BufferedInputStream(new FileInputStream(md.tempLcpsFile)))) {
            int firstIndex = idxIn.readInt();
            byte firstLcp = lcpIn.readByte();
            if (prevRangeLastBucketFirst >= 0) {
                firstLcp = computeLcpByte(sequence, firstIndex, prevRangeLastBucketFirst, 0);
            }
            indexOut.writeInt(firstIndex);
            nlcpOut.writeByte(firstLcp);

            for (int i = 1; i < md.numEntries; i++) {
                indexOut.writeInt(idxIn.readInt());
                nlcpOut.writeByte(lcpIn.readByte());
            }
        }
    }

    private static int resolveSortThreads() {
        String configured = System.getProperty(SA_BUILD_THREADS_PROPERTY);
        if (configured != null) {
            try {
                int n = Integer.parseInt(configured.trim());
                if (n > 0) return n;
            } catch (NumberFormatException ignored) { }
        }
        int procs = Runtime.getRuntime().availableProcessors();
        return Math.max(1, Math.min(procs, MAX_DEFAULT_SA_BUILD_THREADS));
    }

    /**
     * Split bucket ids into contiguous ranges balanced by total suffix count
     * (so each worker has roughly equal sort+LCP work, not equal bucket count).
     */
    private static int[][] partitionBucketIds(Bucket[] buckets, int numThreads) {
        if (numThreads <= 1 || buckets.length == 0) {
            return new int[][]{{0, buckets.length}};
        }
        long totalSuffixes = 0L;
        for (Bucket b : buckets) {
            if (b != null) totalSuffixes += b.size;
        }
        if (totalSuffixes == 0L) {
            return new int[][]{{0, buckets.length}};
        }
        long perThread = (totalSuffixes + numThreads - 1) / numThreads;

        int[][] ranges = new int[numThreads][];
        int rangeStart = 0;
        int rangeIdx = 0;
        long running = 0L;
        for (int i = 0; i < buckets.length; i++) {
            Bucket b = buckets[i];
            if (b != null) running += b.size;
            if (running >= perThread && rangeIdx < numThreads - 1) {
                ranges[rangeIdx++] = new int[]{rangeStart, i + 1};
                rangeStart = i + 1;
                running = 0L;
            }
        }
        ranges[rangeIdx++] = new int[]{rangeStart, buckets.length};
        if (rangeIdx != numThreads) {
            int[][] trimmed = new int[rangeIdx][];
            System.arraycopy(ranges, 0, trimmed, 0, rangeIdx);
            ranges = trimmed;
        }
        return ranges;
    }

    /**
     * Sort each bucket in the range, compute intra-range LCPs, and stream the
     * output into per-worker temp files. The first LCP byte is a placeholder
     * (0) — the merge step rewrites it against the previous range's last
     * bucket. Each bucket's storage is released as soon as it is sorted, so
     * peak heap is bounded by the largest in-flight bucket per thread.
     */
    private static RangeMetadata processBucketRangeToTempFiles(CompactFastaSequence sequence,
                                                               Bucket[] buckets,
                                                               int from,
                                                               int to,
                                                               File tempIndicesFile,
                                                               File tempLcpsFile) throws IOException {
        long count = 0L;
        for (int i = from; i < to; i++) {
            if (buckets[i] != null) count += buckets[i].size;
        }
        if (count == 0L) {
            return new RangeMetadata(null, null, 0, -1);
        }
        if (count > Integer.MAX_VALUE) {
            throw new IllegalStateException("Suffix array bucket range exceeds Integer.MAX_VALUE entries");
        }

        int lastBucketFirstSuffix = -1;
        int prevIntraBucketLast = -1;
        int prevBucketFirst = -1;
        int numEntries = 0;
        boolean firstBucketSeen = false;

        try (DataOutputStream idxOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(tempIndicesFile)));
             DataOutputStream lcpOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(tempLcpsFile)))) {
            for (int bucketId = from; bucketId < to; bucketId++) {
                Bucket bucket = buckets[bucketId];
                if (bucket == null) continue;

                int[] sorted = bucket.trimmedArray();
                buckets[bucketId] = null;
                IntArrays.quickSort(sorted, (a, b) -> compareSuffixesFrom(sequence, a, b, BUCKET_SIZE));

                int first = sorted[0];
                idxOut.writeInt(first);
                byte lcp = firstBucketSeen ? computeLcpByte(sequence, first, prevBucketFirst, 0) : 0;
                lcpOut.writeByte(lcp);
                numEntries++;
                firstBucketSeen = true;
                prevIntraBucketLast = first;

                for (int j = 1; j < sorted.length; j++) {
                    int thisIndex = sorted[j];
                    idxOut.writeInt(thisIndex);
                    lcpOut.writeByte(computeLcpByte(sequence, thisIndex, prevIntraBucketLast, BUCKET_SIZE));
                    numEntries++;
                    prevIntraBucketLast = thisIndex;
                }

                prevBucketFirst = first;
                lastBucketFirstSuffix = first;
            }
        }

        return new RangeMetadata(tempIndicesFile, tempLcpsFile, numEntries, lastBucketFirstSuffix);
    }

    /**
     * Single-thread direct-write path: sort each bucket, compute LCPs, and
     * write to disk in one pass. Used when {@link #SA_BUILD_THREADS_PROPERTY}
     * resolves to 1.
     */
    private static void writeBucketsDirect(CompactFastaSequence sequence,
                                           Bucket[] buckets,
                                           int from,
                                           int to,
                                           DataOutputStream indexOut,
                                           DataOutputStream nlcpOut) throws IOException {
        int prevBucketFirstIndex = -1;
        long lastStatusTime = System.currentTimeMillis();
        for (int i = from; i < to; i++) {
            if (i % 100000 == 0 && System.currentTimeMillis() - lastStatusTime > 2000) {
                lastStatusTime = System.currentTimeMillis();
                System.out.printf("Sorting: %.2f%% complete.%n", (i - from) * 100.0 / (to - from));
            }

            Bucket bucket = buckets[i];
            if (bucket == null) continue;

            int[] sorted = bucket.trimmedArray();
            buckets[i] = null;
            IntArrays.quickSort(sorted, (a, b) -> compareSuffixesFrom(sequence, a, b, BUCKET_SIZE));

            int first = sorted[0];
            byte lcp = 0;
            if (prevBucketFirstIndex >= 0) {
                lcp = computeLcpByte(sequence, first, prevBucketFirstIndex, 0);
            }
            indexOut.writeInt(first);
            nlcpOut.writeByte(lcp);
            int prev = first;

            for (int j = 1; j < sorted.length; j++) {
                int thisIndex = sorted[j];
                indexOut.writeInt(thisIndex);
                lcp = computeLcpByte(sequence, thisIndex, prev, BUCKET_SIZE);
                nlcpOut.writeByte(lcp);
                prev = thisIndex;
            }
            prevBucketFirstIndex = first;
        }
    }

    /** Per-worker sort+LCP output handle. Indices/LCPs live on disk; this carries
     *  the small metadata the merge step needs. Empty ranges return {@code null}
     *  file paths. */
    static final class RangeMetadata {
        final File tempIndicesFile;
        final File tempLcpsFile;
        final int numEntries;
        final int lastBucketFirstSuffix;

        RangeMetadata(File tempIndicesFile, File tempLcpsFile, int numEntries, int lastBucketFirstSuffix) {
            this.tempIndicesFile = tempIndicesFile;
            this.tempLcpsFile = tempLcpsFile;
            this.numEntries = numEntries;
            this.lastBucketFirstSuffix = lastBucketFirstSuffix;
        }
    }

    /** Growable {@code int[]} bucket of suffix indices. Shared between the
     *  bucketing phase (sequential {@link #add}) and the per-range worker
     *  threads (concurrent {@link #trimmedArray} — safe because bucketing
     *  completes before any worker starts). */
    private static final class Bucket {
        private int[] items;
        private int size;

        Bucket() {
            this.items = new int[10];
            this.size = 0;
        }

        void add(int item) {
            if (this.size >= items.length) {
                this.items = Arrays.copyOf(this.items, this.size * 2);
            }
            this.items[this.size++] = item;
        }

        /** Return a fresh int[] of exactly {@code size} entries. The bucket's
         *  internal storage can then be dropped. */
        int[] trimmedArray() {
            return (this.size == this.items.length) ? this.items : Arrays.copyOf(this.items, this.size);
        }
    }

    /**
     * Compare two suffixes of {@code sequence} starting at the given offset.
     * Sign semantics match {@link Comparable#compareTo} and {@link ByteSequence#compareTo};
     * magnitude is not preserved.
     */
    private static int compareSuffixesFrom(CompactFastaSequence sequence, int idxA, int idxB, int startOffset) {
        if (idxA == idxB) return 0;
        long seqSize = sequence.getSize();
        long remainA = seqSize - idxA;
        long remainB = seqSize - idxB;
        long limitLong = Math.min(remainA, remainB);
        int limit = limitLong > ByteSequence.MAX_COMPARISON_LENGTH
                ? ByteSequence.MAX_COMPARISON_LENGTH
                : (int) limitLong;
        for (int offset = startOffset; offset < limit; offset++) {
            byte a = sequence.getByteAt(idxA + offset);
            byte b = sequence.getByteAt(idxB + offset);
            if (a != b) return Byte.compare(a, b); // signed compare, matches ByteSequence.compareTo
        }
        // Shorter suffix sorts first (matches ByteSequence.compareTo semantics).
        return Long.compare(remainA, remainB);
    }

    /** LCP of two suffixes starting from {@code startOffset}, capped at {@link Byte#MAX_VALUE}. */
    private static byte computeLcpByte(CompactFastaSequence sequence, int idxA, int idxB, int startOffset) {
        long seqSize = sequence.getSize();
        long remainA = seqSize - idxA;
        long remainB = seqSize - idxB;
        long limitLong = Math.min(remainA, remainB);
        int limit = limitLong > Byte.MAX_VALUE ? Byte.MAX_VALUE : (int) limitLong;
        int offset = startOffset;
        for (; offset < limit; offset++) {
            byte a = sequence.getByteAt(idxA + offset);
            byte b = sequence.getByteAt(idxB + offset);
            if (a != b) return (byte) offset;
        }
        return (byte) offset;
    }

    @Override
    public String toString() {
        return "Size of the suffix array: " + this.size + "\n";
    }

    public void measureNominalMassError(AminoAcidSet aaSet) throws Exception {
        //		  ArrayList<Pair<Float,Integer>> pepList = new ArrayList<Pair<Float,Integer>>();
        double[] aaMass = new double[128];
        int[] nominalAAMass = new int[128];
        for (int i = 0; i < aaMass.length; i++) {
            aaMass[i] = -1;
            nominalAAMass[i] = -1;
        }

        for (AminoAcid aa : aaSet) {
            aaMass[aa.getResidue()] = aa.getAccurateMass();
            nominalAAMass[aa.getResidue()] = aa.getNominalMass();
        }
        double[] prm = new double[maxPeptideLength];
        int[] nominalPRM = new int[maxPeptideLength];
        int i = Integer.MAX_VALUE - 1000;
        int[] numPeptides = new int[maxPeptideLength];
        int[][] numPepWithError = new int[maxPeptideLength][11];

        DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(getIndexFile())));
        indices.skip(CompactSuffixArray.INT_BYTE_SIZE * 2);    // skip size and id

        DataInputStream nlcps = new DataInputStream(new BufferedInputStream(new FileInputStream(getNeighboringLcpFile())));
        nlcps.skip(CompactSuffixArray.INT_BYTE_SIZE * 2);

        int size = this.getSize();
        int index = -1;
        for (int bufferIndex = 0; bufferIndex < size; bufferIndex++) {
            index = indices.readInt();
            int lcp = nlcps.readByte();

            int idx = sequence.getCharAt(index);
            if (aaMass[idx] <= 0)
                continue;

            if (lcp > i)
                continue;
            for (i = lcp; i < maxPeptideLength; i++) {
                char residue = sequence.getCharAt(index + i);
                double m = aaMass[residue];
                if (m <= 0) {
                    break;
                }
                if (i != 0) {
                    prm[i] = prm[i - 1] + m;
                    nominalPRM[i] = nominalPRM[i - 1] + nominalAAMass[residue];
                } else {
                    prm[i] = m;
                    nominalPRM[i] = nominalAAMass[residue];
                }
                if (i + 1 <= maxPeptideLength) {
                    numPeptides[i]++;
                    int error = (int) Math.round(prm[i] * 0.9995) - nominalPRM[i];
                    error += 5;
                    numPepWithError[i][error]++;
//					System.out.println(index+"\t"+(float)prm[i]+"\t"+sequence.getSubsequence(index, index+i+1));
                }
            }
        }

        long total = 0;
        long totalErr = 0;
        System.out.println("Length\tNumDistinctPeptides\tNumPeptides\tNumPeptidesWithErrors");
        for (i = 0; i < maxPeptideLength; i++) {
            System.out.print((i + 1) + "\t" + this.numDistinctPeptides[i + 1] + "\t" + numPeptides[i]);
            total += numPeptides[i];
            for (int j = 0; j < 11; j++) {
                if (numPepWithError[i][j] > 0) {
                    System.out.print("\t" + (j - 5) + ":" + numPepWithError[i][j]);
                    if (j != 5)
                        totalErr += numPepWithError[i][j];
                }
            }
            System.out.println("\t" + total + "\t" + totalErr + "\t" + (totalErr / (double) total));
        }
        System.out.println("Total #Peptides\t" + total);
        System.out.println("Total #Peptides with nominalMass errors\t" + totalErr + "\t" + totalErr / (double) total);

        indices.close();
        nlcps.close();
    }

    /**
     * Compares two timestamps (typically the lastModified value for a file)
     * If they agree within 2 seconds, returns True, otherwise false
     * @param time1 First file time (milliseconds since 1/1/1970)
     * @param time2 Second file time (milliseconds since 1/1/1970)
     * @return True if the times agree within 2 seconds
     */
    public static boolean NearlyEqualFileTimes(long time1, long time2)
    {
        double timeDiffSeconds = (time1 - time2) / 1000.0;
        if (Math.abs(timeDiffSeconds) <= 2.05)
        {
            return true;
        }

        return false;
    }

}
