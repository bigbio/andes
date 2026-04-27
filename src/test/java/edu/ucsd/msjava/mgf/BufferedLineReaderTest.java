package edu.ucsd.msjava.mgf;

import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Regression test for the BOM-strip fix on {@link BufferedLineReader}: the
 * constructor must invoke {@link UnicodeBOMInputStream#skipBOM()} so the
 * leading byte-order-mark bytes are consumed before the first
 * {@link BufferedLineReader#readLine()} call. Caught by the Copilot review on
 * PR #25.
 */
public class BufferedLineReaderTest {

    private static final byte[] UTF8_BOM = new byte[] {(byte) 0xEF, (byte) 0xBB, (byte) 0xBF};

    @Test
    public void firstLineDoesNotContainUtf8Bom() throws IOException {
        Path tmp = Files.createTempFile("msgfplus-bom-", ".txt");
        try {
            byte[] payload = ("ParentMassTolerance=20ppm\n").getBytes(StandardCharsets.UTF_8);
            byte[] withBom = new byte[UTF8_BOM.length + payload.length];
            System.arraycopy(UTF8_BOM, 0, withBom, 0, UTF8_BOM.length);
            System.arraycopy(payload, 0, withBom, UTF8_BOM.length, payload.length);
            Files.write(tmp, withBom);

            try (BufferedLineReader reader = new BufferedLineReader(tmp.toString())) {
                String first = reader.readLine();
                Assert.assertEquals("BOM bytes must not appear in line 1", "ParentMassTolerance=20ppm", first);
                Assert.assertNull("only one line in fixture", reader.readLine());
            }
        } finally {
            Files.deleteIfExists(tmp);
        }
    }

    @Test
    public void firstLineUnchangedWhenNoBomPresent() throws IOException {
        Path tmp = Files.createTempFile("msgfplus-no-bom-", ".txt");
        try {
            Files.writeString(tmp, "Header\nbody\n");
            try (BufferedLineReader reader = new BufferedLineReader(tmp.toString())) {
                Assert.assertEquals("Header", reader.readLine());
                Assert.assertEquals("body", reader.readLine());
                Assert.assertNull(reader.readLine());
            }
        } finally {
            Files.deleteIfExists(tmp);
        }
    }
}
