package edu.ucsd.msjava.msutil;

public record Annotation(AminoAcid prevAA, Peptide peptide, AminoAcid nextAA) {

    public Annotation(String annotationStr, AminoAcidSet aaSet) {
        this(
                aaSet.getAminoAcid(annotationStr.charAt(0)),
                aaSet.getPeptide(annotationStr.substring(annotationStr.indexOf('.') + 1, annotationStr.lastIndexOf('.'))),
                aaSet.getAminoAcid(annotationStr.charAt(annotationStr.length() - 1))
        );
    }

    public boolean isProteinNTerm() { return prevAA == null; }
    public boolean isProteinCTerm() { return nextAA == null; }

    public AminoAcid getPrevAA() { return prevAA; }
    public Peptide   getPeptide() { return peptide; }
    public AminoAcid getNextAA() { return nextAA; }

    @Override public String toString() {
        if (peptide == null) return null;
        StringBuilder output = new StringBuilder();
        if (prevAA != null) output.append(prevAA.getResidueStr());
        output.append('.').append(peptide).append('.');
        if (nextAA != null) output.append(nextAA.getResidueStr());
        return output.toString();
    }
}
