package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;

import java.util.List;

/**
 * Orchestrates fragment-index construction from a list of unmodified peptide
 * sequences. Phase 2 scope: unmodified sequences only; variable-mod variant
 * enumeration lives in a later phase via a {@code VariantEnumerator} stage
 * upstream of this builder.
 *
 * <p>Build flow:
 * <pre>
 *   for each peptide p:
 *     precursorMass = sum(residueMasses) + H2O
 *     slabIds = slabAssigner.slabsFor(precursorMass)
 *     for each slabId in slabIds:
 *       localPid = slabBuilder[slabId].addPeptide(precursorMass)
 *       peptideTable[slabId].addPeptide(p, precursorMass)
 *       for each theoretical fragment f of p:
 *         slabBuilder[slabId].addFragment(localPid, bucket(f), f.isB)
 *   finalize every slabBuilder into a Slab and put it in the DirectStore
 * </pre>
 *
 * <p>Peptides whose precursor mass falls outside the SlabAssigner's range
 * are silently skipped.
 */
public final class FragmentIndexBuilder {

    private final AminoAcidSet aaSet;
    private final SlabAssigner slabAssigner;
    private final double fragmentBinWidthDa;
    private final TheoreticalFragmentGenerator fragmentGenerator;

    public FragmentIndexBuilder(AminoAcidSet aaSet,
                                 SlabAssigner slabAssigner,
                                 double fragmentBinWidthDa) {
        if (fragmentBinWidthDa <= 0) {
            throw new IllegalArgumentException("fragmentBinWidthDa must be > 0");
        }
        this.aaSet = aaSet;
        this.slabAssigner = slabAssigner;
        this.fragmentBinWidthDa = fragmentBinWidthDa;
        this.fragmentGenerator = new TheoreticalFragmentGenerator(aaSet);
    }

    public FragmentIndex build(List<String> peptides) {
        int numSlabs = slabAssigner.numSlabs();
        SlabBuilder[] slabBuilders = new SlabBuilder[numSlabs];
        PeptideTable[] peptideTables = new PeptideTable[numSlabs];
        for (int s = 0; s < numSlabs; s++) {
            slabBuilders[s] = new SlabBuilder(s, slabAssigner.slabLowMass(s), slabAssigner.slabHighMass(s));
            peptideTables[s] = new PeptideTable();
        }

        for (String peptide : peptides) {
            double precursorMass = computePrecursorMass(peptide);
            if (Double.isNaN(precursorMass)) continue;

            int[] slabIds = slabAssigner.slabsFor(precursorMass);
            if (slabIds.length == 0) continue;

            TheoreticalFragmentGenerator.Fragment[] fragments = fragmentGenerator.fragmentsFor(peptide);
            for (int slabId : slabIds) {
                int localPid = slabBuilders[slabId].addPeptide(precursorMass);
                peptideTables[slabId].addPeptide(peptide, precursorMass);
                for (TheoreticalFragmentGenerator.Fragment f : fragments) {
                    int bucket = (int) (f.mass() / fragmentBinWidthDa);
                    slabBuilders[slabId].addFragment(localPid, bucket, f.isB());
                }
            }
        }

        DirectStore store = new DirectStore(numSlabs);
        for (int s = 0; s < numSlabs; s++) {
            store.putSlab(s, slabBuilders[s].finish());
        }
        return new FragmentIndex(store, peptideTables, slabAssigner, fragmentBinWidthDa);
    }

    private double computePrecursorMass(String peptide) {
        double sum = Composition.H2O;
        for (int i = 0; i < peptide.length(); i++) {
            AminoAcid aa = aaSet.getAminoAcid(peptide.charAt(i));
            if (aa == null) return Double.NaN;
            sum += aa.getAccurateMass();
        }
        return sum;
    }
}
