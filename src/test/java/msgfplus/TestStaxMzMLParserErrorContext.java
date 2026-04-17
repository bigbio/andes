package msgfplus;

import edu.ucsd.msjava.mzml.StaxMzMLParser;
import org.junit.Assert;
import org.junit.Test;

import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Covers Q8: when the mzML has a byte-order mark (BOM) or a malformed XML
 * prolog, the constructor's {@link XMLStreamException} is re-thrown with an
 * actionable message instead of Stax's terse "ParseError in XML prolog".
 */
public class TestStaxMzMLParserErrorContext {

    private File writeBytesToTempMzml(byte[] bytes) throws IOException {
        Path tmp = Files.createTempFile("msgfplus-stax-context-", ".mzML");
        Files.write(tmp, bytes);
        tmp.toFile().deleteOnExit();
        return tmp.toFile();
    }

    @Test
    public void bomPrefixedMzmlGivesActionableMessage() throws Exception {
        // UTF-8 BOM (EF BB BF) followed by a plausible-looking mzML prolog.
        byte[] bom = new byte[]{(byte) 0xEF, (byte) 0xBB, (byte) 0xBF};
        byte[] prolog = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><mzML/>".getBytes(StandardCharsets.UTF_8);
        byte[] content = new byte[bom.length + prolog.length];
        System.arraycopy(bom, 0, content, 0, bom.length);
        System.arraycopy(prolog, 0, content, bom.length, prolog.length);

        File mzml = writeBytesToTempMzml(content);

        try {
            new StaxMzMLParser(mzml);
            // Note: some Stax implementations tolerate a UTF-8 BOM. If this one
            // does, the test becomes a no-op — we can't force the parser to
            // fail, so just return.
        } catch (XMLStreamException e) {
            String msg = e.getMessage();
            Assert.assertNotNull("Wrapped XMLStreamException should carry a message", msg);
            Assert.assertTrue("Message should include the full file path for context",
                    msg.contains(mzml.getAbsolutePath()));
            Assert.assertTrue("Message should mention the BOM / prolog / encoding hint",
                    msg.contains("byte-order mark") || msg.contains("BOM")
                            || msg.contains("XML prolog") || msg.contains("encoding"));
            Assert.assertTrue("Message should point at Troubleshooting.md",
                    msg.contains("Troubleshooting.md"));
        }
    }

    @Test
    public void garbledPrologAlwaysProducesAnnotatedMessage() throws Exception {
        // Definitely-malformed XML (just random text, no prolog at all).
        // Every Stax impl rejects this.
        byte[] garbage = "this is not xml at all".getBytes(StandardCharsets.UTF_8);
        File mzml = writeBytesToTempMzml(garbage);

        try {
            new StaxMzMLParser(mzml);
            Assert.fail("Parsing random bytes as mzML should not succeed");
        } catch (XMLStreamException e) {
            String msg = e.getMessage();
            Assert.assertNotNull(msg);
            Assert.assertTrue("Message should include the index phase tag",
                    msg.contains("during index"));
            Assert.assertTrue("Message should include the file path",
                    msg.contains(mzml.getAbsolutePath()));
            Assert.assertTrue("Original parser error should be preserved in the message",
                    msg.contains("Underlying parser error"));
            Assert.assertSame("Original exception should be the cause",
                    e.getCause().getClass(), XMLStreamException.class);
        }
    }
}
