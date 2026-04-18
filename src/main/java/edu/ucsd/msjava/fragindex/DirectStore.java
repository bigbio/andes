package edu.ucsd.msjava.fragindex;

/**
 * In-memory fragment-index store. Slabs live as {@link Slab} objects on the
 * Java heap by default; Phase 4 will add an off-heap {@link java.nio.DirectByteBuffer}
 * backing so the working set stays out of {@code -Xmx}.
 *
 * <p>For Phase 1 this is the only available store. It's used by unit tests
 * and by in-memory search runs against small FASTA files.
 */
public final class DirectStore implements FragmentIndexStore {
    private final Slab[] slabs;

    public DirectStore(int slabCount) {
        this.slabs = new Slab[slabCount];
    }

    @Override
    public int slabCount() { return slabs.length; }

    @Override
    public void putSlab(int slabId, Slab slab) {
        slabs[slabId] = slab;
    }

    @Override
    public Slab openSlab(int slabId) {
        return slabs[slabId];
    }
}
