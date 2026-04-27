package edu.ucsd.msjava.msutil;

/**
 * Lightweight controlled-vocabulary parameter metadata used by parsers and
 * runtime metadata plumbing without depending on mzIdentML model classes.
 *
 * @author Bryson Gibbons
 */
public record CvParamInfo(String accession, String name, String value,
                          String unitAccession, String unitName) {

    public CvParamInfo(String accession, String name, String value) {
        this(accession, name, value, null, null);
    }

    public boolean hasUnit() {
        return unitAccession != null;
    }

    public String getAccession()     { return accession; }
    public String getName()          { return name; }
    public String getValue()         { return value; }
    public Boolean getHasUnit()      { return hasUnit(); }
    public String getUnitAccession() { return unitAccession; }
    public String getUnitName()      { return unitName; }
}
