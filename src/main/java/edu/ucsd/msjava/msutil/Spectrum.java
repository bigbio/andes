package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.msgf.Tolerance;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * Representation of a mass spectrum object.
 *
 * @author Sangtae Kim
 */
public class Spectrum extends ArrayList<Peak> implements Comparable<Spectrum> {

    public enum Polarity {
        POSITIVE,
        NEGATIVE
    }

    //this is recommended for Serializable objects
    static final private long serialVersionUID = 1L;

    // required members
    private Peak precursor = null;

    // optional members
    private String id;    // unique identifier of the spectrum
    private int startScanNum = -1;
    private int endScanNum = -1;
    private int specIndex = -1;    //
    private String title = null;

    private Peptide annotation = null;
    private ArrayList<String> seqList = null;    // SEQ fields of mgf spectrum
    private float rt = -1;                    // retention time
    private boolean rtIsSeconds = true;      // retention time units - false, is minutes, true, is seconds
    private ActivationMethod activationMethod = null;    // fragmentation method
    private int msLevel = 2;    // ms level
    private Polarity scanPolarity = Polarity.POSITIVE;

    private Boolean isCentroided = true;
    private Boolean externalSetIsCentroided = false;
    private Boolean isCentroidedWithDensePeaks = false;

    private boolean isHighPrecision = false;

    private ArrayList<CvParamInfo> addlCvParams;

    private Float isolationWindowTargetMz = null;

    public Spectrum() {
    }

    public Spectrum(Peak precursorPeak) {
        this.precursor = precursorPeak;
    }

    public Spectrum(float precursorMz, int charge, float precursorIntensity) {
        this.precursor = new Peak(precursorMz, precursorIntensity, charge);
    }

    public String getID() {
        return id;
    }

    public Peptide getAnnotation() {
        return annotation;
    }

    public String getAnnotationStr() {
        if (annotation != null) return annotation.toString();
        return null;
    }

    public ArrayList<String> getSeqList() {
        return seqList;
    }

    public int getCharge() {
        return precursor.getCharge();
    }

    public int getEndScanNum() {
        return endScanNum;
    }

    @Deprecated
    public float getParentMass() {
        return getPrecursorMass();
    }

    /**
     * Gets the monoisotopic (de-charged) precursor mass of this spectrum.
     *
     * @return the mass in Daltons.
     */
    public float getPrecursorMass() {
        return precursor.getMass();
    }

    /**
     * Gets the peptide mass of this spectrum: parentMass-mass(H2O)
     *
     * @return the peptide mass in Daltons.
     */
    public float getPeptideMass() {
        Float peptideMass = precursor.getMass();
        if (peptideMass > 0)
            return peptideMass - (float)Composition.H2O;
        else
            return 0;
    }

    public Peak getPrecursorPeak() {
        return precursor;
    }

    public int getScanNum() {
        return getStartScanNum();
    }

    public int getSpecIndex() {
        return specIndex;
    }

    public int getStartScanNum() {
        return startScanNum;
    }

    public String getTitle() {
        return title;
    }

    public float getRt() {
        return this.rt;
    }

    /** Returns true if retention time is in seconds, false if in minutes. */
    public boolean getRtIsSeconds() {
        return this.rtIsSeconds;
    }

    public ActivationMethod getActivationMethod() {
        return this.activationMethod;
    }

    public Polarity getScanPolarity() {
        return this.scanPolarity;
    }

    public boolean isCentroided() {
        return this.isCentroided;
    }

    /**
     * Whether this spectrum is centroided according to the reader, but failed determineIfCentroided() because peaks are too dense.
     *
     * @return false unless the reader called setIsCentroided(true) and determineIfCentroided() failed
     */
    public boolean isCentroidedWithDensePeaks() {
        return this.isCentroidedWithDensePeaks;
    }

    public boolean isHighPrecision() {
        return this.isHighPrecision;
    }

    public int getMSLevel() {
        return this.msLevel;
    }

    /** Returns additional cvParams to output under the mzIdentML SpectrumIdentificationResult. */
    public ArrayList<CvParamInfo> getAddlCvParams() {
        return this.addlCvParams;
    }

    public void setID(String id) {
        this.id = id;
    }

    public void setAnnotation(Peptide annotation) {
        this.annotation = annotation;
    }

    public void addSEQ(String seq) {
        if (seqList == null)
            seqList = new ArrayList<String>();
        this.seqList.add(seq);
    }

    public void setPrecursor(Peak precursor) {
        this.precursor = precursor;
    }

    public void setStartScanNum(int startScanNum) {
        this.startScanNum = startScanNum;
    }

    public void setEndScanNum(int endScanNum) {
        this.endScanNum = endScanNum;
    }

    public void setScanNum(int scanNum) {
        this.startScanNum = scanNum;
    }

    public void setSpecIndex(int specIndex) {
        this.specIndex = specIndex;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    /** @param rt retention time; see {@link #setRtIsSeconds} for units. */
    public void setRt(float rt) {
        this.rt = rt;
    }

    /** Sets retention time units: true = seconds, false = minutes. */
    public void setRtIsSeconds(boolean isSeconds) {
        this.rtIsSeconds = isSeconds;
    }

    public void setActivationMethod(ActivationMethod fragMethod) {
        this.activationMethod = fragMethod;
    }

    public void setMsLevel(int msLevel) {
        this.msLevel = msLevel;
    }

    public void setScanPolarity(Polarity scanPolarity) {
        this.scanPolarity = scanPolarity;
    }

    public void setIsCentroided(boolean isCentroided) {
        this.isCentroided = isCentroided;
        // track that isCentroided was set from external reader (mzML/mzXML)
        this.externalSetIsCentroided = true;
    }

    public void setIsHighPrecision(boolean isHighPrecision) {
        this.isHighPrecision = isHighPrecision;
    }

    public void setIsolationWindowTargetMz(Float isolationWindowTargetMz) {
        this.isolationWindowTargetMz = isolationWindowTargetMz;
    }

    public Float getIsolationWindowTargetMz() {
        return isolationWindowTargetMz;
    }

    public void determineIsCentroided() {
        boolean centroidedCheckPass = true;

        if (this.size() > 0) {
            ArrayList<Float> diff = new ArrayList<Float>();
            float prevMz = this.get(0).getMz();
            for (int i = 1; i < this.size(); i++) {
                if (this.get(i).getIntensity() == 0)
                    continue;
                float curMz = this.get(i).getMz();
                diff.add((curMz - prevMz) / curMz * 1e6f);
                prevMz = curMz;
            }
            Collections.sort(diff);
            if (diff.size() > 0 && diff.get(diff.size() / 2) < 50) {
                // Check failed - the median PPM distance between peaks is less than 50 PPM
                centroidedCheckPass = false;
            }
        }
        
        if (centroidedCheckPass) {
            this.isCentroided = true;
        } else {
            if (this.isCentroided && this.externalSetIsCentroided) {
                // set a flag to notify the user
                this.isCentroidedWithDensePeaks = true;
            }

            this.isCentroided = false;
        }
    }

    public void setChargeIfSinglyCharged() {
        if (precursor == null || precursor.getCharge() != 0)
            return;
        float tic = 0;
        float ticBelowPrecursor = 0;
        float precursorMz = this.precursor.getMz();
        for (Peak p : this) {
            tic += p.getIntensity();
            if (p.getMz() < precursorMz)
                ticBelowPrecursor += p.getIntensity();
        }

        if (ticBelowPrecursor / tic > 0.9f)
            precursor.setCharge(1);
    }
    
    /**
     * Add an additional cvParam to output as a cvParam under the mzIdentML SpectrumIdentificationResult
     * @param cvParam
     */
    public void addAddlCvParam(CvParamInfo cvParam) {
        if (addlCvParams == null){
            addlCvParams = new ArrayList<CvParamInfo>();
        }

        addlCvParams.add(cvParam);
    }

    @Override
    public String toString() {
        return "Spectrum - mz: " + getPrecursorPeak().getMz() + ", peaks: " + size();
    }

    public Spectrum getCloneWithoutPeakList() {
        Spectrum newSpec = new Spectrum();
        newSpec.precursor = this.precursor.clone();
        newSpec.startScanNum = this.startScanNum;
        newSpec.endScanNum = this.endScanNum;
        newSpec.title = this.title;
        newSpec.seqList = this.seqList;
        newSpec.annotation = this.annotation;
        newSpec.seqList = this.seqList;
        return newSpec;
    }


    public Spectrum getDeconvolutedSpectrum(float toleranceBetweenIsotopes) {
        int charge = this.getCharge();
        if (charge == 0)
            return null;

        Spectrum deconvSpec = this.getCloneWithoutPeakList();
        boolean[] ignore = new boolean[this.size()];
        for (int i = 0; i < this.size(); i++) {
            if (ignore[i])
                continue;
            Peak p = this.get(i);
            float pMz = p.getMz();
            for (int ionCharge = 2; ionCharge < charge && ionCharge < 4; ionCharge++) {
                boolean isDeconvoluted = false;
                for (int j = i + 1; j < this.size(); j++) {
                    Peak p2 = this.get(j);
                    float diff = p2.getMz() - pMz - (float) Composition.ISOTOPE / ionCharge;
                    if (diff > -toleranceBetweenIsotopes && diff < toleranceBetweenIsotopes) {
                        ignore[j] = true;
                        p.setMz(ionCharge * p.getMz() - (ionCharge - 1) * (float) Composition.ChargeCarrierMass());
                        isDeconvoluted = true;
                        float p2Mz = p2.getMz();
                        for (int k = j + 1; k < this.size(); k++) {
                            Peak p3 = this.get(k);
                            float diff2 = p3.getMz() - p2Mz - (float) (Composition.C14 - Composition.C13) / ionCharge;
                            if (diff2 > -toleranceBetweenIsotopes && diff2 < toleranceBetweenIsotopes) {
                                ignore[k] = true;
                                p3.setMz(ionCharge * p3.getMz() - (ionCharge - 1) * (float) Composition.ChargeCarrierMass());
                                deconvSpec.add(p3);
                                break;
                            } else if (diff2 > toleranceBetweenIsotopes)
                                break;
                        }
                        p2.setMz(ionCharge * p2.getMz() - (ionCharge - 1) * (float) Composition.ChargeCarrierMass());
                        deconvSpec.add(p2);
                        break;
                    } else if (diff > toleranceBetweenIsotopes)
                        break;
                }
                if (isDeconvoluted)
                    break;
            }
            deconvSpec.add(p);
        }
        Collections.sort(deconvSpec, new Peak.MassComparator());
        return deconvSpec;
    }

    public void addPeak(Peak peak) {
        this.add(peak);
    }


    public void correctParentMass() {
        if (this.annotation == null || this.getCharge() <= 0)
            return;
        else
            this.precursor.setMz((annotation.getParentMass() + precursor.getCharge() * (float) Composition.ChargeCarrierMass()) / precursor.getCharge());
    }

    public void correctParentMass(float parentMass) {
        this.precursor.setMz((parentMass + precursor.getCharge() * (float) Composition.ChargeCarrierMass()) / precursor.getCharge());
    }

    public void correctParentMass(Peptide pep) {
        if (this.getCharge() <= 0)
            return;
        else
            this.precursor.setMz((pep.getParentMass() + precursor.getCharge() * (float) Composition.ChargeCarrierMass()) / precursor.getCharge());
    }

    public void setCharge(int charge) {
        this.precursor.setCharge(charge);
    }

    public void setPrecursorCharge(int charge) {
        this.precursor.setCharge(charge);
    }

    /**
     * Returns a list of peaks that match the target mass within the tolerance
     * value. The absolute distance between the target mass and a returned peak
     * is less or equal that the tolerance value. The current implementation
     * cycles through all peaks per call.
     *
     * @param mass      target mass.
     * @param tolerance tolerance.
     * @return an ArrayList object of the matching peaks. The array will be empty
     * if there are no peaks within tolerance.
     */
    public ArrayList<Peak> getPeakListByMass(float mass, Tolerance tolerance) {
        float toleranceDa = tolerance.getToleranceAsDa(mass, getCharge());
        return getPeakListByMassRange(mass - toleranceDa, mass + toleranceDa);
    }

    public ArrayList<Peak> getPeakListByMz(float mz, Tolerance tolerance) {
        float toleranceDa = tolerance.getToleranceAsDa(mz);
        return getPeakListByMassRange(mz - toleranceDa, mz + toleranceDa);
    }

    /**
     * Returns the most intense peak that is within tolerance of the target mass.
     * The current implementation takes linear time.
     *
     * @param mass      target mass.
     * @param tolerance tolerance.
     * @return a Peak object if there is match or null otherwise.
     */
    public Peak getPeakByMass(float mass, Tolerance tolerance) {
        ArrayList<Peak> matchList = getPeakListByMass(mass, tolerance);
        if (matchList == null || matchList.size() == 0)
            return null;
        else
            return Collections.max(matchList, new IntensityComparator());
    }

    /**
     * Returns a list of peaks that match the target mass within the specified range.
     * Assuming spectrum is sorted by mass!!!
     *
     * @param minMass minimum mass.
     * @param maxMass maximum mass.
     * @return an ArrayList object of the matching peaks. The array will be empty
     * if there are no peaks within tolerance.
     */
    public ArrayList<Peak> getPeakListByMassRange(float minMass, float maxMass) {
        ArrayList<Peak> matchList = new ArrayList<Peak>();
        int start = Collections.binarySearch(this, new Peak(minMass, 0, 0));
        if (start < 0)
            start = -start - 1;
        for (int i = start; i < this.size(); i++) {
            Peak p = this.get(i);
            if (p.getMz() > maxMass)
                break;
            else
                matchList.add(p);
        }
        return matchList;
    }

    /** Ranks peaks by intensity descending; rank 1 = highest intensity. */
    public void setRanksOfPeaks() {
        ArrayList<Peak> intensitySorted = new ArrayList<Peak>(this);
        Collections.sort(intensitySorted, Collections.reverseOrder(new IntensityComparator()));
        for (int i = 0; i < intensitySorted.size(); i++) {
            intensitySorted.get(i).setRank(i + 1);
        }
    }

    /**
     * Sets intensities of the charge two parent ion and its water loss to 0
     *
     */
    @Deprecated
    public void filterPrecursorPeaks(Tolerance tolerance) {
        filterPrecursorPeaks(tolerance, 0, 0);
    }

    /**
     * Filter (charge-reduced) precursor peaks with the specified offset
     */
    public void filterPrecursorPeaks(Tolerance tolerance, int reducedCharge, float offset) {
        int c = this.getCharge() - reducedCharge;
        float mass = (this.getPrecursorMass() + c * (float) Composition.ChargeCarrierMass()) / c + offset;
        for (Peak p : getPeakListByMass(mass, tolerance))
            p.setIntensity(0);
    }

    public void filterPrecursorPeaksAroundPM() {
        for (int i = 0; i < this.size(); i++) {
            float m = get(i).getMass();
            int nominalMass = Math.round(m * Constants.INTEGER_MASS_SCALER);
            if (nominalMass < 38)
                this.get(i).setIntensity(0);
        }

        // Remove all peaks with masses >= M+H - 38
        int nominalPM = Math.round((getPrecursorMass() - (float) Composition.H2O) * Constants.INTEGER_MASS_SCALER);
        for (int i = this.size() - 1; i >= 0; i--) {
            float m = get(i).getMass();
            int nominalMass = Math.round(m * Constants.INTEGER_MASS_SCALER);
            if (nominalPM - nominalMass >= 38)
                break;
            this.get(i).setIntensity(0);
        }

    }


    public int compareTo(Spectrum s) {
        if (getPrecursorMass() > s.getPrecursorMass())
            return 1;
        else if (getPrecursorMass() < s.getPrecursorMass())
            return -1;
        return 0;
    }

    /**
     * Output this spectrum to the input PrintStream as the mgf format.
     * It needs to be changed later.
     *
     * @param out PrintStream object that the mgf spectrum will be written.
     */
    public void outputMgf(PrintStream out) {
        outputMgf(out, true);
    }

    /**
     * Output this spectrum to the input PrintStream as the mgf format.
     * It needs to be changed later.
     *
     * @param out                   PrintStream object that the mgf spectrum will be written.
     * @param writeActivationMethod don't write ACTIVATION field if false
     */
    public void outputMgf(PrintStream out, boolean writeActivationMethod) {
        out.println("BEGIN IONS");
        if (this.title != null)
            out.println("TITLE=" + getTitle());
        else {
            out.println("TITLE=" + id);
        }
        if (this.annotation != null)
            out.println("SEQ=" + getAnnotationStr());
        if (this.getActivationMethod() != null && writeActivationMethod)
            out.println("ACTIVATION=" + this.getActivationMethod().getName());
        float precursorMz = precursor.getMz();
        out.println("PEPMASS=" + precursorMz);
        if (startScanNum > 0)
            out.println("SCANS=" + startScanNum);
        int charge = getCharge();
        out.println("CHARGE=" + charge + (charge > 0 ? "+" : ""));
        for (Peak p : this)
            if (p.getIntensity() > 0)
                out.println(p.getMz() + "\t" + p.getIntensity());
        out.println("END IONS");
    }

    /**
     * Output this spectrum to the input PrintStream as the dta format.
     * It needs to be changed later.
     *
     * @param fileName dta file name.
     */
    public void outputDta(String fileName) {
        PrintStream out = null;
        try {
            out = new PrintStream(new BufferedOutputStream(new FileOutputStream(fileName)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        out.println(this.getPrecursorMass() + Composition.ChargeCarrierMass() + "\t" + this.getPrecursorPeak().getCharge());
        for (Peak p : this)
            out.println(p.getMz() + "\t" + p.getIntensity());
        out.close();
    }

    /**
     * Convert this spectrum into a dta string representation.
     *
     * @return the dta representation.
     */
    public String toDta() {
        StringBuffer sb = new StringBuffer();
        sb.append(this.getPrecursorMass() + Composition.ChargeCarrierMass() + "\t" + this.getPrecursorPeak().getCharge() + "\n");
        for (Peak p : this) sb.append(p.getMz() + "\t" + p.getIntensity() + "\n");
        return sb.toString();
    }

    class IntensityComparator implements Comparator<Peak> {

        public int compare(Peak o1, Peak o2) {
            if (o1.getIntensity() > o2.getIntensity()) return 1;
            if (o2.getIntensity() > o1.getIntensity()) return -1;
            if (o1.getMz() > o2.getMz()) return 1;
            if (o2.getMz() > o1.getMz()) return -1;
            return 0;
        }

        public boolean equals(Peak o1, Peak o2) {
            return compare(o1, o2) == 0;
        }

    }

    public static SpecFileFormat getSpectrumFileFormat(String specFileName) {
        SpecFileFormat specFormat = null;

        int posDot = specFileName.lastIndexOf('.');
        if (posDot >= 0) {
            String extension = specFileName.substring(posDot);
            if (extension.equalsIgnoreCase(".mzML"))
                specFormat = SpecFileFormat.MZML;
            else if (extension.equalsIgnoreCase(".mgf"))
                specFormat = SpecFileFormat.MGF;
        }

        return specFormat;
    }
}