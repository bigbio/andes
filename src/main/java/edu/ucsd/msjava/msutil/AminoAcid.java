package edu.ucsd.msjava.msutil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;


/**
 * @author Sangtae Kim
 */
public class AminoAcid extends Matter {

    // this is recommended for Serializable objects
    static final private long serialVersionUID = 1L;

    private double mass;
    private int nominalMass;
    private char residue;    // 1 letter code for standard amino acid
    private String name;
    private float probability = 0.05f;
    private Composition composition;

    protected AminoAcid(char residue, String name, Composition composition) {
        this.mass = composition.getAccurateMass();
        this.nominalMass = composition.getNominalMass();
        this.residue = residue;
        this.name = name;
        this.composition = composition;
    }

    protected AminoAcid(char residue, String name, double mass) {
        this.mass = mass;
        this.nominalMass = Math.round(Constants.INTEGER_MASS_SCALER * (float) mass);
        this.residue = residue;
        this.name = name;
    }

    public AminoAcid setProbability(float probability) {
        this.probability = probability;
        return this;
    }

    public String toString() {
        return String.valueOf(residue) + ": " + String.format("%.2f", mass);
    }

    /** Returns false; overridden by {@code ModifiedAminoAcid}. */
    public boolean isModified() {
        return false;
    }

    /** Returns 0; overridden by {@code ModifiedAminoAcid}. */
    public int getNumVariableMods() {
        return 0;
    }

    /** Returns false; overridden by {@code ModifiedAminoAcid}. */
    public boolean hasTerminalVariableMod() {
        return false;
    }

    /** Returns false; overridden by {@code ModifiedAminoAcid}. */
    public boolean hasResidueSpecificVariableMod() {
        return false;
    }

    @Override
    public float getMass() {
        return (float) mass;
    }

    @Override
    public double getAccurateMass() {
        return mass;
    }

    @Override
    public int getNominalMass() {
        return nominalMass;
    }

    public float getProbability() {
        return probability;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof AminoAcid))
            return false;
        AminoAcid aa = (AminoAcid) obj;
        return this == aa;
    }

    public String getResidueStr() {
        return String.valueOf(residue);
    }

    public char getResidue() {
        return residue;
    }

    /** Returns the unmodified residue letter; overridden by ModifiedAminoAcid. */
    public char getUnmodResidue() {
        return residue;
    }

    public String getName() {
        return name;
    }

    public Composition getComposition() {
        return composition;
    }

    public static AminoAcid getStandardAminoAcid(char residue) {
        return residueMap.get(residue);
    }

    public static AminoAcid[] getStandardAminoAcids() {
        return standardAATable;
    }

    public AminoAcid getAAWithFixedModification(Modification mod) {
        String name = mod.getName() + " " + this.getName();
        AminoAcid modAA;
        if (mod.getComposition() == null)
            modAA = getCustomAminoAcid(residue, name, mass + mod.getAccurateMass());
        else
            modAA = getAminoAcid(residue, name, composition.getAddition(mod.getComposition()));
        return modAA;
    }

    public static AminoAcid getCustomAminoAcid(char residue, String name, double mass) {
        AminoAcid standardAA = AminoAcid.getStandardAminoAcid(residue);
        if (standardAA != null && Math.abs(mass - standardAA.getMass()) < 0.001f)
            return standardAA;
        else
            return new AminoAcid(residue, name, mass);
    }

    public static AminoAcid getCustomAminoAcid(char residue, float mass) {
        return new AminoAcid(residue, "Custom amino acid", mass);
    }

    public static AminoAcid getAminoAcid(char residue, String name, Composition composition) {
        AminoAcid standardAA = AminoAcid.getStandardAminoAcid(residue);
        if (standardAA != null && composition.getAccurateMass() == standardAA.getAccurateMass())
            return standardAA;
        else
            return new AminoAcid(residue, name, composition);
    }

    @Override
    public int hashCode() {
        return (int) residue;
    }

    private static Hashtable<Character, AminoAcid> residueMap;
    // Standard amino acids sorted by increasing nominal mass
    private static final AminoAcid[] standardAATable =
            {
                    //                                                   C  H  N  O  S
                    new AminoAcid('G', "Glycine",        new Composition(2, 3, 1, 1, 0)),   // 57.0215
                    new AminoAcid('A', "Alanine",        new Composition(3, 5, 1, 1, 0)),   // 71.0371
                    new AminoAcid('S', "Serine",         new Composition(3, 5, 1, 2, 0)),   // 87.032
                    new AminoAcid('P', "Proline",        new Composition(5, 7, 1, 1, 0)),   // 97.0528
                    new AminoAcid('V', "Valine",         new Composition(5, 9, 1, 1, 0)),   // 99.0684
                    new AminoAcid('T', "Threonine",      new Composition(4, 7, 1, 2, 0)),   // 101.0477
                    new AminoAcid('C', "Cystine",        new Composition(3, 5, 1, 1, 1)),   // 103.0092
                    // new AminoAcid('O', "Hydroxyproline", new Composition(5, 7, 1, 2, 0)),   // 113.0477; note that O could be Hydroxyproline, Ornithine, or Pyrrolysine
                    new AminoAcid('L', "Leucine",        new Composition(6, 11, 1, 1, 0)),  // 113.0841
                    new AminoAcid('I', "Isoleucine",     new Composition(6, 11, 1, 1, 0)),  // 113.0841
                    new AminoAcid('N', "Asparagine",     new Composition(4, 6, 2, 2, 0)),   // 114.0429
                    new AminoAcid('D', "Aspartate",      new Composition(4, 5, 1, 3, 0)),   // 115.0269
                    new AminoAcid('Q', "Glutamine",      new Composition(5, 8, 2, 2, 0)),   // 128.0586
                    new AminoAcid('K', "Lysine",         new Composition(6, 12, 2, 1, 0)),  // 128.095
                    new AminoAcid('E', "Glutamate",      new Composition(5, 7, 1, 3, 0)),   // 129.0426
                    new AminoAcid('M', "Methionine",     new Composition(5, 9, 1, 1, 1)),   // 131.0405
                    new AminoAcid('H', "Histidine",      new Composition(6, 7, 3, 1, 0)),   // 137.0589
                    new AminoAcid('F', "Phenylalanine",  new Composition(9, 9, 1, 1, 0)),   // 147.0684
                    // new AminoAcid('U',  "Selenocysteine", 150.0379),                                    // 150.9536
                    new AminoAcid('R', "Arginine",       new Composition(6, 12, 4, 1, 0)),  // 156.1011
                    new AminoAcid('Y', "Tyrosine",       new Composition(9, 9, 1, 2, 0)),   // 163.0633
                    new AminoAcid('W', "Tryptophan",     new Composition(11, 10, 2, 1, 0)), // 186.0793
            };

    static {
        residueMap = new Hashtable<Character, AminoAcid>();
        for (AminoAcid aa : standardAATable)
            residueMap.put(aa.getResidue(), aa);
    }

    public static ArrayList<AminoAcid> getAminoAcids(int mass) {
        if (mass2aa.containsKey(mass)) return mass2aa.get(mass);
        return new ArrayList<AminoAcid>();
    }

    public static boolean isStdAminoAcid(char c) {
        return residueMap.containsKey(c);
    }

    private static HashMap<Integer, ArrayList<AminoAcid>> mass2aa;

    static {
        mass2aa = new HashMap<Integer, ArrayList<AminoAcid>>();
        for (AminoAcid aa : getStandardAminoAcids()) {
            if (!mass2aa.containsKey(aa.getNominalMass())) {
                mass2aa.put(aa.getNominalMass(), new ArrayList<AminoAcid>());
            }
            mass2aa.get(aa.getNominalMass()).add(aa);
        }
    }
}

