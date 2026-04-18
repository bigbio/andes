package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;

/**
 * Computes theoretical singly-charged b-ion and y-ion m/z values for an
 * unmodified peptide sequence.
 *
 * <p>A peptide of length L yields L-1 b-ions (b1..b_{L-1}) and L-1 y-ions
 * (y1..y_{L-1}). Higher fragment charges and variable-mod aware fragments
 * are future work.
 *
 * <p>Mass formulas (singly-charged):
 * <ul>
 *   <li>b-ion at position i = sum(residues 1..i) + PROTON</li>
 *   <li>y-ion at position i = sum(residues from C-term, length i) + H2O + PROTON</li>
 * </ul>
 */
public final class TheoreticalFragmentGenerator {

    /** Immutable value object representing a single theoretical fragment ion. */
    public static final class Fragment {
        private final double mass;
        private final boolean isB;
        private final int position;  // 1..L-1

        public Fragment(double mass, boolean isB, int position) {
            this.mass = mass;
            this.isB = isB;
            this.position = position;
        }

        /** Monoisotopic m/z of the singly-charged ion. */
        public double mass() { return mass; }

        /** Returns {@code true} for a b-ion, {@code false} for a y-ion. */
        public boolean isB() { return isB; }

        /**
         * Fragment position index (1-based).
         * For b-ions: number of residues from N-terminus included.
         * For y-ions: number of residues from C-terminus included.
         */
        public int position() { return position; }
    }

    private final AminoAcidSet aaSet;

    public TheoreticalFragmentGenerator(AminoAcidSet aaSet) {
        this.aaSet = aaSet;
    }

    /**
     * Computes all b- and y-ions for the given unmodified peptide sequence.
     *
     * <p>Returns an empty array if the peptide is shorter than 2 residues or
     * contains an unrecognised amino acid character.
     *
     * <p>The returned array contains b-ions first (positions 1..L-1), then
     * y-ions (positions 1..L-1), where position 1 is the shortest fragment
     * from each terminus.
     *
     * @param peptide single-letter-code peptide string, upper-case
     * @return array of 2*(L-1) Fragment objects, or empty on error
     */
    public Fragment[] fragmentsFor(String peptide) {
        int L = peptide.length();
        if (L < 2) return new Fragment[0];

        double[] residueMasses = new double[L];
        for (int i = 0; i < L; i++) {
            AminoAcid aa = aaSet.getAminoAcid(peptide.charAt(i));
            if (aa == null) return new Fragment[0];
            residueMasses[i] = aa.getAccurateMass();
        }

        final double proton = Composition.PROTON;
        final double h2o    = Composition.H2O;

        Fragment[] out = new Fragment[2 * (L - 1)];
        int idx = 0;

        // b-ions: b_i covers residues 0..(i-1), position = i
        double bSum = 0;
        for (int i = 0; i < L - 1; i++) {
            bSum += residueMasses[i];
            out[idx++] = new Fragment(bSum + proton, true, i + 1);
        }

        // y-ions: y_i covers the last i residues (C-terminal), position = i
        // Iterate from C-terminus inward: residue[L-1], [L-2], ...
        double ySum = h2o;
        for (int i = L - 1; i >= 1; i--) {
            ySum += residueMasses[i];
            out[idx++] = new Fragment(ySum + proton, false, L - i);
        }

        return out;
    }
}
