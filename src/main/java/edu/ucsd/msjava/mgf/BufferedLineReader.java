package edu.ucsd.msjava.mgf;

import java.io.*;

/**
 * Buffered line reader. Wraps the file in {@link UnicodeBOMInputStream}
 * and consumes the BOM via {@code skipBOM()} so the first line returned by
 * {@link #readLine()} never contains the BOM glyph -- this matters for
 * config / mod / FASTA files saved by Windows editors that prepend a UTF-8
 * BOM.
 */
public class BufferedLineReader extends BufferedReader implements LineReader {

    public BufferedLineReader(String fileName) throws IOException {
        super(new InputStreamReader(new UnicodeBOMInputStream(new FileInputStream(fileName)).skipBOM()));
    }

    @Override
    public String readLine() {
        try {
            return super.readLine();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
}
