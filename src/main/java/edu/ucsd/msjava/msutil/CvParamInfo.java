package edu.ucsd.msjava.msutil;

/**
 * Lightweight controlled-vocabulary parameter metadata used by parsers and
 * runtime metadata plumbing without depending on mzIdentML model classes.
 *
 * @author Bryson Gibbons
 */
public class CvParamInfo {
    private final String accession;
    private final String name;
    private final String value;
    private final String unitAccession;
    private final String unitName;
    private final Boolean hasUnit;

    public CvParamInfo(String accession, String name, String value) {
        this.accession = accession;
        this.name = name;
        this.value = value;
        this.unitAccession = null;
        this.unitName = null;
        this.hasUnit = false;
    }

    public CvParamInfo(String accession, String name, String value, String unitAccession, String unitName) {
        this.accession = accession;
        this.name = name;
        this.value = value;
        this.hasUnit = true;
        this.unitAccession = unitAccession;
        this.unitName = unitName;
    }

    public String getAccession() {
        return this.accession;
    }

    public String getName() {
        return this.name;
    }

    public String getValue() {
        return this.value;
    }

    public Boolean getHasUnit() {
        return this.hasUnit;
    }

    public String getUnitAccession() {
        return this.unitAccession;
    }

    public String getUnitName() {
        return this.unitName;
    }
}
