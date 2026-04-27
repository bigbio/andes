package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.mzml.StaxMzMLParser;
import edu.ucsd.msjava.mzml.StaxMzMLSpectraIterator;
import edu.ucsd.msjava.mzml.StaxMzMLSpectraMap;
import edu.ucsd.msjava.mgf.MgfSpectrumParser;
import edu.ucsd.msjava.mgf.SpectrumParser;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class SpectraAccessor {
    private final File specFile;
    private final SpecFileFormat specFormat;

    private SpectrumParser spectrumParser;

    private StaxMzMLParser staxParser = null;

    private int minMSLevel = 2;
    private int maxMSLevel = 2;

    SpectrumAccessorBySpecIndex specMap = null;
    Iterator<Spectrum> specItr = null;

    /**
     * Constructor that accepts a file
     * Determines the file format based on the file extension
     *
     * @param specFile
     */
    public SpectraAccessor(File specFile) {
        this(specFile, SpecFileFormat.getSpecFileFormat(specFile.getName()));
    }

    /**
     * Constructor that accepts a file and a file format
     *
     * @param specFile
     * @param specFormat
     */
    public SpectraAccessor(File specFile, SpecFileFormat specFormat) {
        if (specFormat == null) {
            throw new IllegalArgumentException("Unsupported spectrum file format: " + specFile.getName());
        }
        this.specFile = specFile;
        this.specFormat = specFormat;
        this.spectrumParser = null;
    }

    /**
     * Set the MS level range for spectrum filtering (both inclusive).
     *
     * @param minMSLevel minimum MS level to consider (inclusive).
     * @param maxMSLevel maximum MS level to consider (inclusive).
     */
    public void setMSLevelRange(int minMSLevel, int maxMSLevel) {
        this.minMSLevel = minMSLevel;
        this.maxMSLevel = maxMSLevel;
    }

    public SpectrumAccessorBySpecIndex getSpecMap() {
        if (specMap == null) {
            if (specFormat == SpecFileFormat.MZML) {
                if (staxParser == null) {
                    try {
                        staxParser = new StaxMzMLParser(specFile, minMSLevel, maxMSLevel);
                    } catch (Exception e) {
                        throw new RuntimeException("Failed to parse mzML file: " + specFile.getAbsolutePath(), e);
                    }
                }
                specMap = new StaxMzMLSpectraMap(staxParser, minMSLevel, maxMSLevel);
            } else if (specFormat == SpecFileFormat.MGF) {
                SpectrumParser parser = new MgfSpectrumParser();
                spectrumParser = parser;
                specMap = new SpectraMap(specFile.getPath(), parser);
            } else {
                return null;
            }
        }

        if (specMap == null) {
            System.out.println("No spectra were found");
            System.out.println("File: " + specFile.getAbsolutePath());
            System.out.println("Format: " + specFormat.getPSIName());
        }
        return specMap;
    }

    public Iterator<Spectrum> getSpecItr() {
        if (specItr == null) {
            if (specFormat == SpecFileFormat.MZML) {
                if (staxParser == null) {
                    try {
                        staxParser = new StaxMzMLParser(specFile, minMSLevel, maxMSLevel);
                    } catch (Exception e) {
                        throw new RuntimeException("Failed to parse mzML file: " + specFile.getAbsolutePath(), e);
                    }
                }
                specItr = new StaxMzMLSpectraIterator(staxParser, minMSLevel, maxMSLevel);
            } else if (specFormat == SpecFileFormat.MGF) {
                SpectrumParser parser = new MgfSpectrumParser();
                spectrumParser = parser;
                try {
                    specItr = new SpectraIterator(specFile.getPath(), parser);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else {
                return null;
            }
        }

        return specItr;
    }

    public Spectrum getSpectrumBySpecIndex(int specIndex) {
        return getSpecMap().getSpectrumBySpecIndex(specIndex);
    }

    public Spectrum getSpectrumById(String specId) {
        return getSpecMap().getSpectrumById(specId);
    }

    /**
     * Get the current spectrum parser, or null if no parser
     * @return
     */
    public SpectrumParser getSpectrumParser() {
        return spectrumParser;
    }

    public String getID(int specIndex) {
        return getSpecMap().getID(specIndex);
    }

    public float getPrecursorMz(int specIndex) {
        return getSpecMap().getPrecursorMz(specIndex);
    }

    public String getTitle(int specIndex) {
        return getSpecMap().getTitle(specIndex);
    }

    public CvParamInfo getSpectrumIDFormatCvParam() {
        CvParamInfo cvParam = null;
        if (specFormat == SpecFileFormat.MGF)
            cvParam = new CvParamInfo("MS:1000774", "multiple peak list nativeID format", null);
        else if (specFormat == SpecFileFormat.MZML) {
            if (staxParser == null) {
                try {
                    staxParser = new StaxMzMLParser(specFile);
                } catch (Exception e) {
                    throw new RuntimeException("Failed to parse mzML file: " + specFile.getAbsolutePath(), e);
                }
            }
            String[] idFormat = staxParser.detectSpectrumIDFormat();
            if (idFormat != null) {
                cvParam = new CvParamInfo(idFormat[0], idFormat[1], null);
            } else {
                throw new IllegalStateException("Unsupported mzML format: " + specFile.getAbsolutePath()
                        + " does not contain a child term of MS:1000767 (native spectrum identifier format)");
            }
        }

        return cvParam;
    }

}
