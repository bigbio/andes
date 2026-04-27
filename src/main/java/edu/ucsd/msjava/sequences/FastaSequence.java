package edu.ucsd.msjava.sequences;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.Map.Entry;

/** Sequence implementation backed by a FASTA file. */
public class FastaSequence implements Sequence {

    //this is the file in which the sequence was generated
    private String baseFilepath;

    // used for writing the encoded binary sequence.
    private final String seqExtension;

    // maps the terminator character position of this sequence to its annotation
    private TreeMap<Integer, String> annotations;

    // maps the header strings of the fasta entries to the position of the terminators
    private TreeMap<String, Integer> header2ends;

    // the contents of the sequence concatenated into a long string
    private ByteBuffer sequence;

    // the original serialized fasta file
    private ByteBuffer original;

    // the number of characters in the buffer
    private int size;

    // the alphabet map
    private HashMap<Character, Byte> alpha2byte;

    // the reverse translation map
    private HashMap<Byte, Character> byte2alpha;

    // the string representation of the alphabet
    private String alphabetString;

    // the identifier for this sequence
    private int id;


    // initialize alphabet from a colon-separated string
    private void initializeAlphabet(String s) {
        String[] tokens = s.split(":");
        this.alpha2byte = new HashMap<Character, Byte>();
        this.byte2alpha = new HashMap<Byte, Character>();
        this.byte2alpha.put(Constants.TERMINATOR, Constants.TERMINATOR_CHAR);
        for (byte i = 0, value = 1; i < tokens.length; i++, value++) {
            for (int j = 0; j < tokens[i].length(); j++) {
                alpha2byte.put(tokens[i].charAt(j), value);
            }
            byte2alpha.put(value, tokens[i].charAt(0));
        }
    }

    private void createObjectFromRawFile(String filepath) {

        // a rough estimate of the space required to hold everything
        int bufferSize = (int) new File(filepath).length();
        ByteBuffer sequence = ByteBuffer.allocate(bufferSize);
        StringBuffer original = new StringBuffer();
        HashMap<Integer, String> annotations = new HashMap<Integer, String>();
        HashMap<Character, Byte> alpha2byte = new HashMap<Character, Byte>();
        String alphabet = "";
        byte alphabetSize = 1;
        int size = 0;
        int id = UUID.randomUUID().hashCode();

        // read the fasta file
        try {
            BufferedReader in = new BufferedReader(new FileReader(filepath));

            Integer offset = 0;
            String annotation = null;
            String s;              //
            while ((s = in.readLine()) != null) {

                // this is a regular fasta line
                if (!s.startsWith(">")) {
                    for (int index = 0; index < s.length(); index++) {
                        Byte encoded = alpha2byte.get(s.charAt(index));
                        if (encoded != null) {
                            sequence.put(encoded);
                        } else {
                            sequence.put(alphabetSize);
                            alpha2byte.put(s.charAt(index), alphabetSize++);
                            alphabet += ":" + s.charAt(index);
                        }
                        original.append(s.charAt(index));
                    }
                    offset += s.length();
                }

                // annotation line
                else {
                    sequence.put(Constants.TERMINATOR);
                    original.append('_');
                    // the offset always points to the terminator of this sequence
                    if (annotation != null) annotations.put(offset, annotation);

                    // remember for the next annotation
                    offset++;
                    annotation = s.substring(1);
                }
            }
            sequence.put(Constants.TERMINATOR);
            original.append('_');
            offset++;
            // the offset always points to the terminator of this sequence
            annotations.put(offset, annotation);
            size = offset;
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        writeMetaInfo(annotations, alphabet.substring(1), size, id);
        writeSequence(original, sequence, size, id);
    }

    private void createObjectFromRawFile(String filepath, String alphabet) {

        // estimate the length of the buffer
        int bufferSize = (int) new File(filepath).length();
        ByteBuffer sequence = ByteBuffer.allocate(bufferSize);
        StringBuffer original = new StringBuffer();
        HashMap<Integer, String> annotations = new HashMap<Integer, String>();
        int size = 0;
        int id = UUID.randomUUID().hashCode();

        // initialization
        initializeAlphabet(alphabet);

        // read the fasta file
        try {
            BufferedReader in = new BufferedReader(new FileReader(filepath));

            Integer offset = 0;
            String annotation = null;
            String s;              //
            while ((s = in.readLine()) != null) {

                // this is a regular fasta line (not annotation)
                if (!s.startsWith(">")) {
                    for (int index = 0; index < s.length(); index++) {
                        Byte encoded = this.alpha2byte.get(s.charAt(index));
                        if (encoded != null) {
                            sequence.put(encoded);
                        } else {
                            sequence.put(Constants.TERMINATOR);
                        }
                        original.append(s.charAt(index));
                    }
                    offset += s.length();
                }

                // annotation line
                else {

                    // terminate the last sequence
                    sequence.put(Constants.TERMINATOR);
                    original.append(Constants.TERMINATOR_CHAR);

                    // the offset always points to the terminator of this sequence
                    if (annotation != null) annotations.put(offset, annotation);

                    // remember for the next annotation
                    offset++;
                    annotation = s.substring(1);
                }
            }

            // process the last sequence
            sequence.put(Constants.TERMINATOR);
            original.append(Constants.TERMINATOR_CHAR);
            offset++;
            // the offset always points to the terminator of this sequence
            annotations.put(offset, annotation);
            size = offset;
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        writeMetaInfo(annotations, alphabet, size, id);
        writeSequence(original, sequence, size, id);
    }

    private void writeMetaInfo(HashMap<Integer, String> annotations, String alphabet, int size, int id) {
        String filepath = this.baseFilepath + this.seqExtension + "anno";
        try {
            PrintWriter out = new PrintWriter(filepath);
            out.println(size);
            out.println(id);
            out.println(alphabet);
            Set<Integer> keys = annotations.keySet();
            for (Integer key : keys) {
                out.println(key + ":" + annotations.get(key));
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private int readMetaInfo() {
        String filepath = this.baseFilepath + this.seqExtension + "anno";
        try {
            BufferedReader in = new BufferedReader(new FileReader(filepath));
            this.size = Integer.parseInt(in.readLine());
            int id = Integer.parseInt(in.readLine());
            this.alphabetString = in.readLine().trim();
            this.annotations = new TreeMap<Integer, String>();
            for (String line = in.readLine(); line != null; line = in.readLine()) {
                String[] tokens = line.split(":", 2);
                this.annotations.put(Integer.parseInt(tokens[0]), tokens[1]);
            }
            in.close();
            return id;
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return 0;
    }

    private void writeSequence(StringBuffer original, ByteBuffer sequence, int size, int id) {
        String filepath = this.baseFilepath + this.seqExtension;
        try {
            DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filepath)));
            out.writeInt(size);
            out.writeInt(id);
            for (int i = 0; i < size; i++) {
                out.writeByte(sequence.get(i));
            }
            out.write(original.toString().getBytes());
            out.flush();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private int readSequence() {
        String filepath = this.baseFilepath + this.seqExtension;
        try {
            DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(filepath)));
            int size = in.readInt();
            int id = in.readInt();
            byte[] sequenceArr = new byte[size];
            in.read(sequenceArr);
            sequence = ByteBuffer.wrap(sequenceArr).asReadOnlyBuffer();
            byte[] originalArr = new byte[size];
            in.read(originalArr);
            original = ByteBuffer.wrap(originalArr).asReadOnlyBuffer();
            in.close();
            return id;
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return 0;
    }


    public FastaSequence(String filepath) {
        this(filepath, null);
    }

    /** Letters not in the alphabet are encoded as TERMINATOR. */
    public FastaSequence(String filepath, String alphabet) {
        this(filepath, alphabet, Constants.FILE_EXTENSION);
    }

    /** Letters not in the alphabet are encoded as TERMINATOR. */
    public FastaSequence(String filepath, String alphabet, String seqExtension) {

        this.seqExtension = seqExtension;

        String[] tokens = filepath.split("\\.");
        String extension = tokens[tokens.length - 1];
        String basepath = filepath.substring(0, filepath.length() - extension.length() - 1);

        this.baseFilepath = basepath;
        if (!extension.equalsIgnoreCase("fasta") && !extension.equalsIgnoreCase("fa")) {
            System.err.println("Input error: not a fasta file");
            System.exit(-1);
        }

        String metaFile = basepath + this.seqExtension + Constants.ANNO_FILE_SUFFIX;
        String sequenceFile = basepath + seqExtension;
        if (!new File(metaFile).exists() || !new File(sequenceFile).exists()) {
            if (alphabet != null) createObjectFromRawFile(filepath, alphabet);
            else createObjectFromRawFile(filepath);

        }

        int metaId = readMetaInfo();
        int seqId = readSequence();

        if (metaId == seqId) {
            initializeAlphabet(this.alphabetString);
            //initializeAlphabet(alphabet);
            this.id = metaId;
        } else {
            System.err.println("The files " + metaFile + " and " + sequenceFile + " have different ids.");
            System.err.println("The problem can be solved by recreating the files");
            System.exit(-1);
        }

        // populate the header2ends map
        this.header2ends = new TreeMap<String, Integer>();
        for (int position : this.annotations.keySet()) {
            this.header2ends.put(this.annotations.get(position), position);
        }
    }


    public Set<Byte> getAlphabetAsBytes() {
        return this.byte2alpha.keySet();
    }

    public Collection<Character> getAlphabet() {
        ArrayList<Character> results = new ArrayList<Character>();
        for (char c : this.byte2alpha.values())
            if (c != '_') results.add(c);
        return results;
    }

    public boolean isTerminator(long position) {
        return getByteAt(position) == Constants.TERMINATOR;
    }

    public char toChar(byte b) {
        if (byte2alpha.containsKey(b)) return byte2alpha.get(b);
        return '?';
    }

    public int getAlphabetSize() {
        return this.byte2alpha.size();
    }

    public long getSize() {
        return this.size;
    }

    public byte getByteAt(long position) {
        // forget boundary check for faster access
        if (position >= this.size) return Constants.TERMINATOR;
        return this.sequence.get((int) position);
    }

    public String getSubsequence(long start, long end) {
        if (start >= end || end > this.size) return null;
        char[] seq = new char[(int) (end - start)];
        for (long i = start; i < end; i++) {
            seq[(int) (i - start)] = (char) this.original.get((int) i);
        }
        return new String(seq);
    }

    public char getCharAt(long position) {
        return (char) this.original.get((int) position);
    }

    public String toString(byte[] sequence) {
        String retVal = "";
        for (byte item : sequence) {
            Character c = byte2alpha.get(item);
            if (c != null) retVal += c;
            else retVal += '?';
        }
        return retVal;
    }

    public byte toByte(char c) {
        return alpha2byte.get(c);
    }

    public byte[] getBytes(int start, int end) {
        byte[] result = new byte[end - start];
        for (int i = start; i < end; i++) {
            result[i - start] = getByteAt(i);
        }
        return result;
    }

    public boolean isInAlphabet(char c) {
        return alpha2byte.containsKey(c);
    }

    public boolean isValid(long position) {
        if (isTerminator(position)) return false;
        return isInAlphabet(getCharAt(position));
    }

    public int getId() {
        return this.id;
    }

    public String getAnnotation(long position) {
        Entry<Integer, String> entry = annotations.higherEntry((int) position);
        if (entry != null)
            return entry.getValue();
        else
            return null;
    }

    public long getStartPosition(long position) {
        Integer startPos = annotations.floorKey((int) position);
        if (startPos == null) {
            return 0;
        }
        return startPos;
    }

    public String getMatchingEntry(long position) {
        Integer start = annotations.floorKey((int) position);     // always "_" at start
        Integer end = annotations.higherKey((int) position);       // exclusive
        if (start == null) start = 0;
        if (end == null) end = (int) this.getSize();
        while (!isValid(end - 1)) end--;     // ensure that the last character is valid (exclusive)
        return this.getSubsequence(start + 1, end);
    }

    public String getMatchingEntry(String name) {
        String key = this.header2ends.ceilingKey(name);
        if (key == null || !key.startsWith(name)) return null;
        int position = this.header2ends.get(key) - 1;
        Integer start = annotations.floorKey(position);   // always "_" at start
        Integer end = annotations.higherKey(position);    // exclusive
        if (start == null) start = 0;
        if (end == null) end = (int) this.getSize();
        while (!isValid(end - 1)) end--;     // ensure that the last character is valid (exclusive)
        return this.getSubsequence(start + 1, end);
    }

    public void setBaseFilepath(String baseFilepath) {
        this.baseFilepath = baseFilepath;
    }

    public String getBaseFilepath() {
        return this.baseFilepath;
    }

    public void set(long start, char c) {
        this.sequence.put((int) start, this.alpha2byte.get(c));
        this.original.put((int) start, (byte) c);
    }

    /** Must be called before set() — read-only ByteBuffers do not support put(). */
    public void makeModifiable() {
        ByteBuffer sequenceCopy = ByteBuffer.allocateDirect(this.size);
        ByteBuffer originalCopy = ByteBuffer.allocateDirect(this.size);
        sequenceCopy.put(this.sequence);
        originalCopy.put(this.original);
        this.sequence = sequenceCopy;
        this.original = originalCopy;
    }

    public List<String> getAnnotations() {
        return new ArrayList<String>(annotations.values());
    }
}
