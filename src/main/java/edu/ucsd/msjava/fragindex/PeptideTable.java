package edu.ucsd.msjava.fragindex;

import java.util.ArrayList;
import java.util.List;

/**
 * Per-slab peptide metadata store.
 *
 * <p>Stores the peptide sequence and precursor monoisotopic mass for each
 * peptide assigned to one slab. Peptide ids are slab-local and assigned
 * sequentially starting at 0, matching the ids that {@link SlabBuilder}
 * uses. The sequence is the unmodified peptide string; variable-mod
 * configuration will be added as an additional column in a later phase.
 *
 * <p>Mutable during build; treated as read-only after the slab is finalized.
 */
public final class PeptideTable {

    private final List<String> sequences = new ArrayList<>();
    private final List<Double> precursorMasses = new ArrayList<>();

    public int addPeptide(String sequence, double precursorMassDa) {
        int id = sequences.size();
        sequences.add(sequence);
        precursorMasses.add(precursorMassDa);
        return id;
    }

    public int size() { return sequences.size(); }

    public String sequence(int peptideId) {
        return sequences.get(peptideId);
    }

    public double precursorMass(int peptideId) {
        return precursorMasses.get(peptideId);
    }
}
