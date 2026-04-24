package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.sequences.Constants;
import edu.ucsd.msjava.suffixarray.ByteSequence;
import edu.ucsd.msjava.suffixarray.SuffixFactory;
import it.unimi.dsi.fastutil.ints.IntArrays;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Locale;

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

        //System.out.println("LastModified times in the existing csarr and cnlcp files " +
        //        "match the LastModified time of the sequence file (" + lastModified + ")");

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

    /**
     * Helper method that creates the suffixFile.
     *
     * <p>Algorithm: two-phase radix-then-sort. Each suffix is hashed by its
     * first {@link #BUCKET_SIZE} residues into one of {@code alphabetSize^BUCKET_SIZE}
     * buckets; within a bucket, suffixes are sorted lexicographically from offset
     * {@code BUCKET_SIZE} onward (the bucket prefix is shared by construction).
     *
     * <p>Scaling notes:
     *
     * <ul>
     *   <li>Buckets store raw {@code int[]} indices and are sorted in-place via
     *       {@link IntArrays#quickSort} with a comparator that reads from the
     *       sequence directly. Prior revisions materialised a
     *       {@code SuffixFactory.Suffix[]} per bucket before sorting; that
     *       allocated ~32 bytes of Java object overhead per suffix, which on a
     *       100 MB+ FASTA added multiple GB of transient heap churn and tripped
     *       OOM on an 8 GB JVM. The int-based path stores 4 bytes per suffix
     *       and reuses the sequence directly during compare/LCP.</li>
     *   <li>LCP between adjacent suffixes is computed by an int-based helper
     *       on the sequence — no {@code ByteSequence} wrappers allocated in the
     *       sort/write loop.</li>
     *   <li>Output is written sequentially: one bucket finalised + freed before
     *       the next begins, so peak memory during sort is dominated by the
     *       largest individual bucket rather than the full set.</li>
     * </ul>
     *
     * @param sequence  the Adapter object that represents the database (text).
     * @param indexFile newly created index file.
     * @param nlcpFile  newly created nlcp file.
     */
    private void createSuffixArrayFiles(CompactFastaSequence sequence, File indexFile, File nlcpFile) {
        System.out.println("Creating the suffix array indexed file... Size: " + sequence.getSize());

        // helper local class: a growable int[] buffer for suffix indices.
        class Bucket {

            private int[] items;
            private int size;

            public Bucket() {
                this.items = new int[10];
                this.size = 0;
            }

            public void add(int item) {
                if (this.size >= items.length) {
                    this.items = Arrays.copyOf(this.items, this.size * 2);
                }
                this.items[this.size++] = item;
            }

            /**
             * Trim + sort the bucket's indices lexicographically by suffix, starting
             * from offset {@link #BUCKET_SIZE} since all members of a bucket share
             * the first {@code BUCKET_SIZE} residues by hash construction. Returns
             * a new {@code int[]} of exactly {@link #size} entries; the internal
             * storage can be released after this call.
             */
            public int[] sortedIndexArray() {
                int[] out = (this.size == this.items.length) ? this.items : Arrays.copyOf(this.items, this.size);
                IntArrays.quickSort(out, (a, b) -> compareSuffixesFrom(sequence, a, b, BUCKET_SIZE));
                return out;
            }
        }

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
            // Track the last index emitted across bucket boundaries so we can
            // compute the LCP between the last suffix of the previous bucket
            // and the first suffix of the current one. -1 = no predecessor.
            int prevBucketFirstIndex = -1;

            System.out.println("Sorting suffixes... Size: " + bucketSuffixes.length);
            lastStatusTime = System.currentTimeMillis();

            for (int i = 0; i < bucketSuffixes.length; i++) {
                // print progress
                if (i % 100000 == 0 && System.currentTimeMillis() - lastStatusTime > 2000) {
                    lastStatusTime = System.currentTimeMillis();
                    System.out.printf("Sorting: %.2f%% complete.\n", i * 100.0 / bucketSuffixes.length);
                }

                Bucket bucket = bucketSuffixes[i];
                if (bucket == null) continue;

                int[] sortedIndices = bucket.sortedIndexArray();
                bucketSuffixes[i] = null; // release the Bucket wrapper ASAP

                int first = sortedIndices[0];
                byte lcp = 0;
                if (prevBucketFirstIndex >= 0) {
                    // Inter-bucket LCP starts at offset 0 (buckets only share prefix
                    // with their OWN members, not neighbours).
                    lcp = computeLcpByte(sequence, first, prevBucketFirstIndex, 0);
                }
                indexOut.writeInt(first);
                nlcpOut.writeByte(lcp);
                int prev = first;

                for (int j = 1; j < sortedIndices.length; j++) {
                    int thisIndex = sortedIndices[j];
                    indexOut.writeInt(thisIndex);
                    // Intra-bucket: shared prefix of BUCKET_SIZE bytes is guaranteed.
                    lcp = computeLcpByte(sequence, thisIndex, prev, BUCKET_SIZE);
                    nlcpOut.writeByte(lcp);
                    prev = thisIndex;
                }
                prevBucketFirstIndex = first;
            }

            bucketSuffixes = null;

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
     * Compare two suffixes of {@code sequence} starting at the given offset. Sign
     * semantics match {@link Comparable#compareTo}; magnitude is not preserved
     * (we only need the sort order, not the LCP — that's computed separately).
     * Comparison length is capped at {@link ByteSequence#MAX_COMPARISON_LENGTH}
     * to match the legacy {@code Suffix.compareTo} behaviour.
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

    /**
     * Compute the LCP (longest common prefix length) of two suffixes starting
     * from {@code startOffset}. Capped at {@link Byte#MAX_VALUE} so the result
     * fits in a single byte for on-disk storage.
     */
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
        String retVal = "Size of the suffix array: " + this.size + "\n";
//		int rank = 0;
//		while(indices.hasRemaining()) {
//			int index = indices.get();
//			int lcp = this.neighboringLcps.get(rank);
//			retVal += rank + "\t" + index + "\t" + lcp + "\t" + sequence.toString(factory.makeSuffix(index).getSequence()) + "\n";
//			rank++;
//		}
//		indices.rewind();        // reset marks after iteration
//		neighboringLcps.rewind();
        return retVal;
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
