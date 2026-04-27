package msgfplus;

import edu.ucsd.msjava.cli.MSGFPlusOptions;
import org.junit.Assert;
import org.junit.Test;
import picocli.CommandLine;

public class TestMinSpectraPerThread {

    @Test
    public void defaultIs250() {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        Assert.assertEquals(250, opts.effectiveMinSpectraPerThread());
    }

    @Test
    public void overrideAppliesThroughGetter() {
        MSGFPlusOptions opts = new MSGFPlusOptions();
        MSGFPlusOptions.commandLine(opts).parseArgs("-minSpectraPerThread", "50");
        Assert.assertEquals(50, opts.effectiveMinSpectraPerThread());
    }

    @Test
    public void parsesZero() {
        // Picocli has no min-value enforcement on Integer fields by default,
        // so '0' is parseable here. Range checks moved to SearchParams.parse
        // (which would reject zero earlier in the search-engine flow if needed).
        MSGFPlusOptions opts = new MSGFPlusOptions();
        MSGFPlusOptions.commandLine(opts).parseArgs("-minSpectraPerThread", "0");
        Assert.assertEquals(0, opts.effectiveMinSpectraPerThread());
    }
}
