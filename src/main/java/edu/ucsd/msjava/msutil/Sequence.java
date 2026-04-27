package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.msgf.MassListComparator;
import edu.ucsd.msjava.msgf.Tolerance;

import java.util.ArrayList;
import java.util.HashSet;


/**
 * Superclass for a list of masses. Peptide, GappedPeptide, Tag should extend
 * this class.
 *
 * @author jung
 */
public class Sequence<T extends Matter> extends ArrayList<T> {

    //this is recommended for Serializable objects
    static final private long serialVersionUID = 1L;


    public float getMass() {
        return getMass(0, this.size());
    }

    public double getAccurateMass() {
        return getMass(0, this.size());
    }

    /** Sum of masses in [from, to), clamped to [0, size). */
    public float getMass(int from, int to) {
        from = java.lang.Math.max(from, 0);
        to = java.lang.Math.min(to, this.size());
        float sum = 0.f;
        for (int i = from; i < to; i++)
            sum += this.get(i).getMass();
        return sum;
    }

    public double getAccurateMass(int from, int to) {
        from = java.lang.Math.max(from, 0);
        to = java.lang.Math.min(to, this.size());
        double sum = 0;
        for (int i = from; i < to; i++)
            sum += this.get(i).getAccurateMass();
        return sum;
    }

    public Sequence<T> subSequence(int fromIndex, int toIndex) {
        return (Sequence<T>) super.subList(fromIndex, toIndex);
    }

    public String toString() {
        StringBuffer output = new StringBuffer();
        for (T matter : this) {
            output.append(matter.toString() + " ");
        }
        return output.toString();
    }

    public static <T extends Matter> Sequence<T> getIntersection(Sequence<T> seq1, Sequence<T> seq2) {
        Sequence<T> union = new Sequence<T>();
        HashSet<T> set = new HashSet<T>();
        for (T m : seq1)
            set.add(m);
        for (T m : seq2)
            if (set.contains(m))
                union.add(m);
        return union;
    }

    public boolean isMatchedTo(Peptide peptide, Tolerance tolerance, boolean isPrefix) {
        ArrayList<Mass> pepMassList = new ArrayList<Mass>();
        float mass = 0;
        for (int i = 0; i < peptide.size(); i++) {
            if (isPrefix)
                mass += peptide.get(i).getMass();
            else
                mass += peptide.get(peptide.size() - 1 - i).getMass();
            pepMassList.add(new Mass(mass));
        }
        ArrayList<Mass> massList = new ArrayList<Mass>();
        for (int i = 0; i < this.size(); i++)
            massList.add(new Mass(this.get(i).getMass()));
        MassListComparator<Mass> comparator = new MassListComparator<Mass>(pepMassList, massList);
        int matchSize = comparator.getMatchedList(tolerance).length;
        return (matchSize == this.size());
    }

    public boolean isMatchedToNominalMasses(Peptide peptide, boolean isPrefix) {
        HashSet<Integer> massList = new HashSet<Integer>();
        int mass = 0;
        for (int i = 0; i < peptide.size(); i++) {
            if (isPrefix)
                mass += peptide.get(i).getNominalMass();
            else
                mass += peptide.get(peptide.size() - 1 - i).getNominalMass();
            massList.add(mass);
        }
        for (Matter m : this) {
            if (!massList.contains(m.getNominalMass()))
                return false;
        }
        return true;
    }

    public float[] toMassArray() {
        float[] massArr = new float[this.size()];
        int index = 0;
        for (T m : this)
            massArr[index++] = m.getMass();
        return massArr;
    }
}
