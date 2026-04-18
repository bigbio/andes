package edu.ucsd.msjava.fragindex;

/**
 * Storage backend for fragment-index slabs. Implementations differ in whether
 * slabs live off-heap in memory ({@link DirectStore}) or on disk via mmap
 * (to be added in Phase 4).
 *
 * <p>All methods must be thread-safe for concurrent readers after
 * {@link #putSlab(int, Slab)} has been called for each slab during build.
 */
public interface FragmentIndexStore {
    /** Total number of slabs this store holds (set at construction). */
    int slabCount();

    /** Install a slab at the given id. Called once per slab during build. */
    void putSlab(int slabId, Slab slab);

    /** Return the slab at the given id, or null if none has been put yet. */
    Slab openSlab(int slabId);

    /** Optional hint that the caller has finished with the slab. Default: no-op. */
    default void closeSlab(Slab slab) {}
}
