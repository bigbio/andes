package edu.ucsd.msjava.ui;

import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msscorer.ScoringParameterGeneratorWithErrors;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.AnnotatedSpectra;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Trainer entry point. The MS-GF+ {@code params/ParamManager} CLI framework
 * was removed from this fork; this restored entry point uses a simple
 * {@code args[]} parser instead. Inputs that depended on the deleted
 * {@code mzid/} package are not accepted - supply pre-converted .tsv files.
 */
public class ScoringParamGen {

    public static final int VERSION = 8831;
    public static final String DATE = "02/04/2013";

    public static void main(String argv[]) {
        if (argv.length == 0 || hasHelpFlag(argv)) {
            printUsageInfo(null);
            return;
        }

        String errMsg = run(argv);
        if (errMsg != null) {
            System.err.println("[Error] " + errMsg);
            System.out.println();
            printUsageInfo(null);
        }
    }

    private static boolean hasHelpFlag(String[] argv) {
        for (String a : argv) {
            if (a.equalsIgnoreCase("-h") || a.equalsIgnoreCase("--help") || a.equalsIgnoreCase("-help"))
                return true;
        }
        return false;
    }

    static String run(String[] argv) {
        // -i <results.tsv>[,results2.tsv,...]   (one or more comma-separated training result TSVs)
        // -d <specDir>                          (directory holding the spectrum files referenced by the TSVs)
        // -m <activationMethodName>             (e.g. CID, ETD, HCD, UVPD)
        // -inst <instrumentTypeName>            (e.g. LowRes, HighRes)
        // -e <enzymeName>                       (e.g. Tryp)
        // -protocol <protocolName>              (e.g. NoProtocol, optional, default automatic)
        // -thread <numThreads>                  (default 1)
        // -dropErrors 0|1                       (default 0)
        // -mgf 0|1                              (default 0; if 1, also writes <dataType>.mgf)
        File[] resultFiles = null;
        File specDir = null;
        ActivationMethod activationMethod = null;
        InstrumentType instType = null;
        Enzyme enzyme = null;
        Protocol protocol = null;
        int numThreads = 1;
        boolean dropErrors = false;
        boolean createMgf = false;

        for (int i = 0; i < argv.length; i += 2) {
            if (!argv[i].startsWith("-") || i + 1 >= argv.length) {
                return "Invalid parameter: " + argv[i];
            }
            String key = argv[i];
            String val = argv[i + 1];
            if (key.equalsIgnoreCase("-i")) {
                String[] paths = val.split(",");
                List<File> files = new ArrayList<File>(paths.length);
                for (String p : paths) {
                    File f = new File(p);
                    if (!f.exists())
                        return "Input file does not exist: " + p;
                    files.add(f);
                }
                resultFiles = files.toArray(new File[0]);
            } else if (key.equalsIgnoreCase("-d")) {
                specDir = new File(val);
                if (!specDir.exists() || !specDir.isDirectory())
                    return "Spectrum directory does not exist or is not a directory: " + val;
            } else if (key.equalsIgnoreCase("-m")) {
                activationMethod = ActivationMethod.get(val);
                if (activationMethod == null)
                    return "Unrecognized activation method: " + val;
            } else if (key.equalsIgnoreCase("-inst")) {
                instType = InstrumentType.get(val);
                if (instType == null)
                    return "Unrecognized instrument type: " + val;
            } else if (key.equalsIgnoreCase("-e")) {
                enzyme = Enzyme.getEnzymeByName(val);
                if (enzyme == null)
                    return "Unrecognized enzyme: " + val;
            } else if (key.equalsIgnoreCase("-protocol")) {
                protocol = Protocol.get(val);
                if (protocol == null)
                    return "Unrecognized protocol: " + val;
            } else if (key.equalsIgnoreCase("-thread")) {
                try {
                    numThreads = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    return "-thread must be an integer";
                }
            } else if (key.equalsIgnoreCase("-dropErrors")) {
                dropErrors = "1".equals(val);
            } else if (key.equalsIgnoreCase("-mgf")) {
                createMgf = "1".equals(val);
            } else {
                return "Unknown option: " + key;
            }
        }

        if (resultFiles == null) return "missing -i (training result TSV files)";
        if (specDir == null) return "missing -d (spectrum directory)";
        if (activationMethod == null) return "missing -m (activation method)";
        if (instType == null) return "missing -inst (instrument type)";
        if (enzyme == null) return "missing -e (enzyme)";
        if (protocol == null) protocol = Protocol.AUTOMATIC;

        long time = System.currentTimeMillis();
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();

        AnnotatedSpectra annotatedSpec = new AnnotatedSpectra(resultFiles, specDir, aaSet);
        System.out.println("Reading training PSMs...");
        String parseErr = annotatedSpec.parse(numThreads, dropErrors);
        if (parseErr != null) {
            if (dropErrors) {
                System.out.println("Datasets with errors (dropped): " + parseErr);
            } else {
                return parseErr;
            }
        }
        if (annotatedSpec.getAnnotatedSpecContainer() == null
                || annotatedSpec.getAnnotatedSpecContainer().isEmpty()) {
            return "No results to train on. Exiting.";
        }
        System.out.println("Done.");

        SpecDataType dataType = new SpecDataType(activationMethod, instType, enzyme, protocol);

        if (createMgf) {
            String mgfFileName = dataType.toString() + ".mgf";
            File mgfFile = new File(mgfFileName);
            System.out.println("Creating " + mgfFile.getPath());
            try {
                PrintStream mgfOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(mgfFile)));
                annotatedSpec.writeToMgf(mgfOut);
                mgfOut.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        ScoringParameterGeneratorWithErrors.generateParameters(
                annotatedSpec.getAnnotatedSpecContainer(),
                dataType,
                aaSet,
                new File("."),
                false,
                true);

        System.out.format("ScoringParamGen complete (total elapsed time: %.2f sec)\n",
                (System.currentTimeMillis() - time) / (float) 1000);
        return null;
    }

    public static void printUsageInfo(String message) {
        if (message != null) {
            System.err.println(message);
        }
        System.out.println("ScoringParamGen v" + VERSION + " (" + DATE + ")");
        System.out.println("Trains MS-GF+ scoring parameters (.param) from annotated PSM training data.");
        System.out.println();
        System.out.println("Usage: java -Xmx2000M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen [options]");
        System.out.println();
        System.out.println("Required:");
        System.out.println("  -i  <tsv1[,tsv2,...]>  Training result TSV files (mzID input not supported in this build)");
        System.out.println("  -d  <specDir>          Directory holding the spectrum files referenced by the TSVs");
        System.out.println("  -m  <activation>       Activation method (e.g. CID, ETD, HCD, UVPD)");
        System.out.println("  -inst <instrument>     Instrument type (e.g. LowRes, HighRes, QExactive)");
        System.out.println("  -e  <enzyme>           Enzyme name (e.g. Tryp, Chymotryp, LysC, AspN)");
        System.out.println();
        System.out.println("Optional:");
        System.out.println("  -protocol <name>       Protocol (default: NoProtocol/automatic)");
        System.out.println("  -thread <int>          Worker threads for parsing PSMs (default: 1)");
        System.out.println("  -dropErrors 0|1        Drop datasets with errors instead of failing (default: 0)");
        System.out.println("  -mgf 0|1               Also emit aggregated <dataType>.mgf (default: 0)");
        System.out.println();
        System.out.println("Notes:");
        System.out.println("  * The output .param file is written to the current working directory.");
        System.out.println("  * mzID input was supported by upstream MS-GF+ but the mzid/ package was removed");
        System.out.println("    from this fork; pre-convert .mzid to .tsv if needed.");
    }
}
