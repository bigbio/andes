package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.mzid.Constants;
import edu.ucsd.msjava.mzml.StaxMzMLParser;
import edu.ucsd.msjava.mzml.StaxMzMLSpectraIterator;
import edu.ucsd.msjava.mzml.StaxMzMLSpectraMap;
import edu.ucsd.msjava.parser.*;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;

import java.io.File;
import java.io.FileNotFoundException;
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
                        staxParser = new StaxMzMLParser(specFile);
                    } catch (Exception e) {
                        throw new RuntimeException("Failed to parse mzML file: " + specFile.getAbsolutePath(), e);
                    }
                }
                specMap = new StaxMzMLSpectraMap(staxParser, minMSLevel, maxMSLevel);
            } else if (specFormat == SpecFileFormat.DTA_TXT)
                specMap = new PNNLSpectraMap(specFile.getPath());
            else {
                SpectrumParser parser = null;
                if (specFormat == SpecFileFormat.MGF)
                    parser = new MgfSpectrumParser();
                else if (specFormat == SpecFileFormat.MS2)
                    parser = new MS2SpectrumParser();
                else if (specFormat == SpecFileFormat.PKL)
                    parser = new PklSpectrumParser();
                else
                    return null;

                spectrumParser = parser;
                specMap = new SpectraMap(specFile.getPath(), parser);
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
                        staxParser = new StaxMzMLParser(specFile);
                    } catch (Exception e) {
                        throw new RuntimeException("Failed to parse mzML file: " + specFile.getAbsolutePath(), e);
                    }
                }
                specItr = new StaxMzMLSpectraIterator(staxParser, minMSLevel, maxMSLevel);
            } else if (specFormat == SpecFileFormat.DTA_TXT)
                try {
                    specItr = new PNNLSpectraIterator(specFile.getPath());
                } catch (IOException e) {
                    e.printStackTrace();
                }
            else {
                SpectrumParser parser = null;
                if (specFormat == SpecFileFormat.MGF)
                    parser = new MgfSpectrumParser();
                else if (specFormat == SpecFileFormat.MS2)
                    parser = new MS2SpectrumParser();
                else if (specFormat == SpecFileFormat.PKL)
                    parser = new PklSpectrumParser();
                else
                    return null;

                spectrumParser = parser;
                try {
                    specItr = new SpectraIterator(specFile.getPath(), parser);
                } catch (IOException e) {
                    e.printStackTrace();
                }
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

    public CvParam getSpectrumIDFormatCvParam() {
        CvParam cvParam = null;
        if (specFormat == SpecFileFormat.DTA_TXT
                || specFormat == SpecFileFormat.MGF
                || specFormat == SpecFileFormat.PKL
                || specFormat == SpecFileFormat.MS2
        )
            cvParam = Constants.makeCvParam("MS:1000774", "multiple peak list nativeID format");
        else if (specFormat == SpecFileFormat.MZDATA)
            cvParam = Constants.makeCvParam("MS:1000777", "spectrum identifier nativeID format");
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
                cvParam = Constants.makeCvParam(idFormat[0], idFormat[1]);
            } else {
                throw new IllegalStateException("Unsupported mzML format: " + specFile.getAbsolutePath()
                        + " does not contain a child term of MS:1000767 (native spectrum identifier format)");
            }
        }

        return cvParam;
    }

}
