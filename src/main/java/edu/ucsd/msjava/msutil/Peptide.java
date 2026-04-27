package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.msgf.IntMassFactory;
import edu.ucsd.msjava.msgf.IntMassFactory.IntMass;
import edu.ucsd.msjava.msgf.MassListComparator;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.Modification.Location;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class Peptide extends Sequence<AminoAcid> implements Comparable<Peptide> {

    //this is recommended for Serializable objects
    static final private long serialVersionUID = 1L;
    // maximum length of a peptide
    static final int MAX_LENGTH = 30;

    // fields
    private boolean isModified; // Indicates the peptide has a modified amino acid

    static final boolean FAIL_WHEN_PEPTIDE_IS_MODIFIED = false; // Fail loudly

    // true if this peptide contains invalid amino acid
    private boolean isInvalid = false;

    /** Parses a sequence string, supporting N-term mods (e.g. +42ACDEFGR) and inline mods (e.g. QSV+2.12QLK). Not fully implemented for all edge cases. */
    public Peptide(String sequence, AminoAcidSet aaSet) {
        isModified = false;
        int seqLen = sequence.length();
        int index = 0;

        float nTermModMass = 0;

        // sequence has an N-term fixed mod
        while (index < seqLen) {
            char c = sequence.charAt(index);
            if (c == '-' || c == '+')    // sequence has an N-term mod (e.g. +42ACDEFGR)
            {
                int startIndex = index;
                while (++index < seqLen) {
                    c = sequence.charAt(index);
                    if (!Character.isDigit(c) && c != '.')
                        break;
                }
                nTermModMass += Float.parseFloat(sequence.substring(startIndex, index));
            } else
                break;
        }

        boolean isNTerm = true;
        for (; index < seqLen; index++) {
            char c = sequence.charAt(index);
            assert (Character.isLetter(c)) : "Error in string at index " + index;
            float mod = 0f;
            if (index + 1 < seqLen) { // Check for modification (e.g. +17, -12.5)
                char sign = sequence.charAt(index + 1);
                while (sign == '-' || sign == '+') { // Modification found
                    assert (index + 2 < seqLen) : "Missing value after \"" + sign + "\"";
                    assert (c >= 'A' && c <= 'Z' || c >= 'a' && c <= 'z') : "Error in string at index " + index + 2;
                    int startModIdx = index + 2;
                    int endModIdx = startModIdx + 1;
                    // Extends substring to find modification value
                    while (endModIdx < seqLen &&
                            (sequence.charAt(endModIdx) == '.' ||
                                    sequence.charAt(endModIdx) >= '0' && sequence.charAt(endModIdx) <= '9')) {
                        endModIdx++; // A+76
                    }
                    float modMass = Float.parseFloat(sequence.substring(startModIdx, endModIdx));
                    if (sign == '-') modMass *= -1f;
                    mod += modMass;
                    index = endModIdx - 1;
                    if (endModIdx < sequence.length())
                        sign = sequence.charAt(endModIdx);
                    else
                        break;
                }
                if (index + 4 < seqLen && sign == 'p' && sequence.charAt(index + 2) == 'h')    // phos
                {
                    assert (sequence.charAt(index + 3) == 'o');
                    assert (sequence.charAt(index + 4) == 's');
                    mod = 79.966331f;
                    index += 4;
                } else if (index + 4 < seqLen && sign >= 'a' && sign <= 'z' && (Character.toUpperCase(sign) == c) && (sequence.charAt(index + 2) == '-'))    // mutation or phosphorylation
                {
                    assert (sequence.charAt(index + 3) == '>');
                    char mutatedResidue = sequence.charAt(index + 4);
                    assert (mutatedResidue >= 'a' && mutatedResidue <= 'z');
                    c = Character.toUpperCase(mutatedResidue);
                    index += 4;
                }
            }

            AminoAcid aa;
            if (isNTerm) {
                aa = aaSet.getAminoAcid(Location.N_Term, c);
                isNTerm = false;
            } else
                aa = aaSet.getAminoAcid(c);

            // TODO: how to deal C-term fixed mods
            if (!Character.isUpperCase(c) || aa == null)    // not a valid amino acid
            {
                this.isInvalid = true;
                return;
            }
            if (this.size() == 0)
                mod += nTermModMass;

            if (mod == 0f) this.add(aa);
            else { // modified
                isModified = true; // Now peptide is modified
                float mass = aa.getMass() + mod;
                AminoAcid modAA = VolatileAminoAcid.getVolatileAminoAcid(mass);
                this.add(modAA);
            }
        }
    }

    public Peptide(String sequence) {
        this(sequence, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
    }

    public Peptide(ArrayList<AminoAcid> aaArray) {
        for (AminoAcid aa : aaArray) {
            assert (aa != null) : "Null amino acid";
            this.add(aa);
        }
    }

    public Peptide(List<AminoAcid> aaArray) {
        for (AminoAcid aa : aaArray) {
            assert (aa != null) : "Null amino acid";
            this.add(aa);
        }
    }

    public Peptide(AminoAcid[] aaArray) {
        for (AminoAcid aa : aaArray) this.add(aa);
    }


    public Peptide subPeptide(int fromIndex, int toIndex) {
        return (Peptide) super.subSequence(fromIndex, toIndex);
    }

    public Peptide setModified() {
        isModified = true;
        return this;
    }

    public Peptide setModified(boolean isModified) {
        this.isModified = isModified;
        return this;
    }

    /** Returns boolean array indexed by nominal mass; true at each prefix-mass position. */
    public boolean[] getBooleanPeptide() {
        boolean[] boolPeptide = new boolean[this.getNominalMass() + 1];
        int mass = 0;
        for (AminoAcid aa : this) {
            mass += aa.getNominalMass();
            boolPeptide[mass] = true;
        }
        return boolPeptide;
    }


    public boolean isGappedPeptideTrue(ArrayList<Integer> gp) {
        boolean[] boolPeptide = getBooleanPeptide();
        boolean isTrue = true;
        for (int m : gp)
            if (boolPeptide[m] == false)
                isTrue = boolPeptide[m];
        return isTrue;
    }

    public boolean isInvalid() {
        return this.isInvalid;
    }

    public boolean isCTermModified() {
        return get(this.size() - 1).isModified();
    }


    public boolean hasTrypticCTerm() {
        AminoAcid cTerm = this.get(this.size() - 1);
        return !isCTermModified() &&
                (cTerm == AminoAcid.getStandardAminoAcid('K') || cTerm == AminoAcid.getStandardAminoAcid('R'));
    }

    public boolean hasCleavageSite(Enzyme enzyme) {
        AminoAcid target;
        if (enzyme.isCTerm())
            target = this.get(this.size() - 1);
        else
            target = this.get(0);
        return enzyme.isCleavable(target);
    }

    public AminoAcid get(int i) {
        if (i <= -1) // N-terminal
            return null;
        else if (i >= this.size()) // C-terminal
            return null;
        return super.get(i);
    }


    public int compareTo(Peptide other) {
        // funky ordering
        int minSize = java.lang.Math.min(this.size(), other.size());

        for (int i = 0; i < minSize; i++) {
            int r = get(i).compareTo(other.get(i));
            if (r != 0) {
                return r;
            }
        }

        int r = size() - other.size();
        if (r > 0) {
            return 1;
        } else if (r < 0) {
            return -1;
        }
        return 0;
    }

    public boolean equalsIgnoreIL(Peptide pep) {
        if (this.size() != pep.size())
            return false;
        for (int i = 0; i < this.size(); i++) {
            Composition c1 = this.get(i).getComposition();
            Composition c2 = pep.get(i).getComposition();
            if (!c1.equals(c2))
                return false;
        }
        return true;
    }

    public String toString() {
        StringBuffer output = new StringBuffer();
        for (AminoAcid aa : this) {
            output.append(aa.getResidueStr());
        }
        return output.toString();
    }

    public Sequence<Composition> toCumulativeCompositionSequence(boolean isPrefix, Composition offset) {
        Sequence<Composition> seq = new Sequence<Composition>();
        Composition c = offset;
        for (int i = 0; i < this.size(); i++) {
            if (isPrefix) {
                c = c.getAddition(this.get(i).getComposition());
                seq.add(c);
            } else {
                c = c.getAddition(this.get(this.size() - 1 - i).getComposition());
                seq.add(c);
            }
        }
        return seq;
    }

    public Sequence<Composition> toCompositionSequence() {
        Sequence<Composition> seq = new Sequence<Composition>();
        for (AminoAcid aa : this)
            seq.add(aa.getComposition());
        return seq;
    }

    public Sequence<Composition> toReverseCompositionSequence() {
        Sequence<Composition> seq = new Sequence<Composition>();
        for (int i = this.size() - 1; i >= 0; i--)
            seq.add(this.get(i).getComposition());
        return seq;
    }

    public Sequence<IntMass> toPrefixIntMassSequence(IntMassFactory factory) {
        Sequence<IntMass> seq = new Sequence<IntMass>();
        for (int i = 0; i < this.size(); i++)
            seq.add(factory.getInstance(this.get(i).getMass()));
        return seq;
    }

    public Sequence<IntMass> toCumulativeIntMassSequence(boolean isPrefix, IntMassFactory factory) {
        Sequence<IntMass> seq = new Sequence<IntMass>();
        float mass = 0;
        for (int i = 0; i < this.size(); i++) {
            if (isPrefix) {
                mass += this.get(i).getMass();
                seq.add(factory.getInstance(mass));
            } else {
                mass += this.get(this.size() - 1 - i).getMass();
                seq.add(factory.getInstance(mass));
            }
        }
        return seq;
    }

    public Sequence<IntMass> toSuffixIntMassSequence(IntMassFactory factory) {
        Sequence<IntMass> seq = new Sequence<IntMass>();
        for (int i = this.size() - 1; i >= 0; i--)
            seq.add(factory.getInstance(this.get(i).getMass()));
        return seq;
    }

    /** Sum of residue masses plus H2O (neutral monoisotopic peptide mass). */
    public float getParentMass() {
        return getMass() + (float) Composition.H2O;
    }

    public int getNumSymmetricPeaks(Tolerance tolerance) {
        ArrayList<Composition> bIons = toCumulativeCompositionSequence(true, new Composition(0, 1, 0, 0, 0));
        ArrayList<Composition> yIons = toCumulativeCompositionSequence(false, new Composition(0, 3, 0, 1, 0));
        MassListComparator<Composition> comparator = new MassListComparator<Composition>(bIons, yIons);

        return comparator.getMatchedList(tolerance).length;
    }

    /** Uses nominal masses. */
    public int getNumSymmetricPeaks() {
        int numSymmPeaks = 0;
        HashSet<Integer> bIons = new HashSet<Integer>();
        int bMass = 1;
        for (int i = 0; i < this.size(); i++) {
            bMass += this.get(i).getNominalMass();
            bIons.add(bMass);
        }
        int yMass = 19;
        for (int i = this.size() - 1; i >= 0; i--) {
            yMass += this.get(i).getNominalMass();
            if (bIons.contains(yMass))
                numSymmPeaks++;
        }
        return numSymmPeaks;
    }

    public int getNominalMass() {
        int sum = 0;
        for (AminoAcid aa : this) {
            sum += aa.getNominalMass();
        }
        return sum;
    }

    public int getIntMassIndex(IntMassFactory factory) {
        int sum = 0;
        for (AminoAcid aa : this) {
            sum += factory.getMassIndex(aa.getMass());
        }
        return sum;
    }

    public Composition getComposition() {
        Composition c = new Composition(0);
        for (AminoAcid aa : this)
            c.add(aa.getComposition());
        return c;
    }

    public float getProbability() {
        float prob = 1;
        for (int i = 0; i < this.size(); i++) {
            AminoAcid aa = this.get(i);
            prob *= aa.getProbability();
        }
        return prob;
    }


    public float getNumber() {
        float number = 1;
        AminoAcid aaL = AminoAcid.getStandardAminoAcid('L');
        AminoAcid aaI = AminoAcid.getStandardAminoAcid('I');
        AminoAcid aaQ = AminoAcid.getStandardAminoAcid('Q');
        AminoAcid aaK = AminoAcid.getStandardAminoAcid('K');
        for (int i = 0; i < this.size(); i++) {
            AminoAcid aa = this.get(i);
            if (aa == aaL || aa == aaI || aa == aaQ || aa == aaK)
                number *= 2;
        }
        return number;
    }


    public Peptide slice(int from, int to) {
        from = java.lang.Math.max(0, from);
        to = java.lang.Math.min(this.size(), to);

        ArrayList<AminoAcid> aaList = new ArrayList<AminoAcid>();
        for (int i = from; i < to; i++)
            aaList.add(this.get(i));
        if (aaList.size() > 0) {
            return new Peptide(aaList);
        }
        return null;
    }


    public static Peptide getSequence(String seq) {
        ArrayList<AminoAcid> aaList = new ArrayList<AminoAcid>();
        int seqLen = seq.length();
        for (int i = 0; i < seqLen; i++) {
            aaList.add(AminoAcid.getStandardAminoAcid(seq.charAt(i)));
        }
        return new Peptide(aaList);
    }


    public boolean isCorrect(ArrayList<Integer> masses) {
        int cumMass = 0;
        int massIndex = 0;
        int targetMass = masses.get(massIndex++);
        for (AminoAcid aa : this) {
            cumMass += aa.getNominalMass();
            if (cumMass < targetMass) {
                continue;  // move to the next mass
            }

            if (cumMass == targetMass) {
                // we got a match
                if (massIndex < masses.size())
                    targetMass += masses.get(massIndex++);
                else
                    // we matched everything
                    return true;
            } else {
                // no match
                return false;
            }
        }

        return massIndex == masses.size();
    }


    public static boolean isCorrect(String sequence, ArrayList<Integer> masses, AminoAcidSet aaSet) {
        int cumMass = 0;
        int massIndex = 0;
        int targetMass = masses.get(massIndex++);
        for (int i = 0; i < sequence.length(); i++) {
            cumMass += aaSet.getAminoAcid(sequence.charAt(i)).getNominalMass();
            if (cumMass < targetMass) {
                continue;  // move to the next mass
            }

            if (cumMass == targetMass) {
                // we got a match
                if (massIndex < masses.size())
                    targetMass += masses.get(massIndex++);
                else
                    // we matched everything
                    return true;
            } else {
                // no match
                return false;
            }
        }

        return massIndex == masses.size();
    }


    public static boolean isCorrect(String sequence, ArrayList<Integer> masses) {
        return isCorrect(sequence, masses, AminoAcidSet.getStandardAminoAcidSet());
    }


    public float[] getPRMMasses(boolean isPrefix, float offset) {
        if (isModified) // TODO handle modified peptide
            return null;
        float[] masses = new float[this.size() - 1];
        float mass = offset;

        for (int i = 0; i < this.size() - 1; i++) {
            if (isPrefix)
                mass += this.get(i).getMass();
            else
                mass += this.get(this.size() - 1 - i).getMass();
            masses[i] = mass;
        }
        return masses;
    }

    public boolean isModified() {
        return isModified;
    }


    public static float getMassFromString(String peptide) {
        float cumMass = 0;
        for (int i = peptide.length(), j = 0; i > 0; i--, j++) {
            cumMass += AminoAcid.getStandardAminoAcid(peptide.charAt(j)).getMass();

        }
        return cumMass;
    }


}
