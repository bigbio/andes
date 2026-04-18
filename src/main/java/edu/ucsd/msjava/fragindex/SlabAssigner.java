package edu.ucsd.msjava.fragindex;

/**
 * Routes peptides to precursor-mass slabs. A slab covers a range
 * {@code [lo, hi)} of width {@code slabWidthDa}. Peptides within
 * {@code boundaryOverlapDa} of a slab boundary are replicated into the
 * adjacent slab so search-time queries always land in exactly one slab.
 *
 * <p>Stateless and thread-safe: every method is a pure function of the
 * input mass and the constructor parameters.
 */
public final class SlabAssigner {

    private final double minMassDa;
    private final double maxMassDa;
    private final double slabWidthDa;
    private final double boundaryOverlapDa;
    private final int numSlabs;

    public SlabAssigner(double minMassDa, double maxMassDa, double slabWidthDa, double boundaryOverlapDa) {
        if (minMassDa >= maxMassDa) throw new IllegalArgumentException("minMassDa must be < maxMassDa");
        if (slabWidthDa <= 0) throw new IllegalArgumentException("slabWidthDa must be > 0");
        if (boundaryOverlapDa < 0) throw new IllegalArgumentException("boundaryOverlapDa must be >= 0");
        this.minMassDa = minMassDa;
        this.maxMassDa = maxMassDa;
        this.slabWidthDa = slabWidthDa;
        this.boundaryOverlapDa = boundaryOverlapDa;
        this.numSlabs = (int) Math.ceil((maxMassDa - minMassDa) / slabWidthDa);
    }

    public int numSlabs() { return numSlabs; }

    public double slabLowMass(int slabId) {
        return minMassDa + slabId * slabWidthDa;
    }

    public double slabHighMass(int slabId) {
        return minMassDa + (slabId + 1) * slabWidthDa;
    }

    public int[] slabsFor(double peptideMassDa) {
        if (peptideMassDa < minMassDa || peptideMassDa >= maxMassDa) return new int[0];
        int primary = (int) Math.floor((peptideMassDa - minMassDa) / slabWidthDa);
        if (primary < 0) primary = 0;
        if (primary >= numSlabs) primary = numSlabs - 1;

        boolean replicatePrev = primary > 0
                && (peptideMassDa - slabLowMass(primary)) < boundaryOverlapDa;
        boolean replicateNext = primary < numSlabs - 1
                && (slabHighMass(primary) - peptideMassDa) < boundaryOverlapDa;

        if (replicatePrev && replicateNext) return new int[]{primary - 1, primary, primary + 1};
        if (replicatePrev) return new int[]{primary - 1, primary};
        if (replicateNext) return new int[]{primary, primary + 1};
        return new int[]{primary};
    }
}
