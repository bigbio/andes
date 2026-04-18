package edu.ucsd.msjava.fragindex;

/**
 * Top-level fragment index: a store of precursor-mass slabs plus the
 * per-slab peptide metadata needed to resolve search hits back to
 * sequences + masses.
 *
 * <p>Constructed by {@link FragmentIndexBuilder}; held by the search
 * engine for the duration of a run. Immutable after construction.
 */
public final class FragmentIndex {

    private final FragmentIndexStore store;
    private final PeptideTable[] peptideTables;
    private final SlabAssigner slabAssigner;
    private final double fragmentBinWidthDa;

    FragmentIndex(FragmentIndexStore store,
                  PeptideTable[] peptideTables,
                  SlabAssigner slabAssigner,
                  double fragmentBinWidthDa) {
        this.store = store;
        this.peptideTables = peptideTables;
        this.slabAssigner = slabAssigner;
        this.fragmentBinWidthDa = fragmentBinWidthDa;
    }

    public FragmentIndexStore store() { return store; }
    public SlabAssigner slabAssigner() { return slabAssigner; }
    public double fragmentBinWidthDa() { return fragmentBinWidthDa; }
    public int numSlabs() { return slabAssigner.numSlabs(); }
    public PeptideTable peptideTable(int slabId) { return peptideTables[slabId]; }
    public Slab slab(int slabId) { return store.openSlab(slabId); }

    /** Convenience: total peptide count across all slabs (sum of table sizes). */
    public int totalPeptideEntries() {
        int total = 0;
        for (PeptideTable t : peptideTables) total += t.size();
        return total;
    }
}
