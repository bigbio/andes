package edu.ucsd.msjava.msutil;

import java.util.ArrayList;


public class SpecFileFormat extends FileFormat {
    private final String psiAccession;
    private final String psiName;

    private SpecFileFormat(String suffix, String psiAccession, String psiName) {
        super(suffix);
        this.psiAccession = psiAccession;
        this.psiName = psiName;
    }

    public String getPSIAccession() {
        return psiAccession;
    }

    public String getPSIName() {
        return psiName;
    }

    public static final SpecFileFormat MGF;
    public static final SpecFileFormat MZML;

    public static SpecFileFormat getSpecFileFormat(String specFileName) {
        String lowerCaseFileName = specFileName.toLowerCase();
        for (SpecFileFormat f : specFileFormatList) {
            for (String suffix : f.getSuffixes()) {
                if (lowerCaseFileName.endsWith(suffix.toLowerCase()))
                    return f;
            }
        }
        return null;
    }

    private static ArrayList<SpecFileFormat> specFileFormatList;

    static {
        MGF = new SpecFileFormat(".mgf", "MS:1001062", "Mascot MGF file");
        MZML = new SpecFileFormat(".mzML", "MS:1000584", "mzML file");

        specFileFormatList = new ArrayList<SpecFileFormat>();
        specFileFormatList.add(MGF);
        specFileFormatList.add(MZML);
    }
}
