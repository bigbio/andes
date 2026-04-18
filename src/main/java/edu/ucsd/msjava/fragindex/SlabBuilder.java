package edu.ucsd.msjava.fragindex;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

/**
 * Writable buffer for assembling a single slab of the fragment index.
 *
 * <p>Used during index build. Caller flow:
 * <pre>
 *   SlabBuilder b = new SlabBuilder(slabId, minMassDa, maxMassDa);
 *   int pid = b.addPeptide(peptideMassDa);
 *   b.addFragment(pid, fragmentBucket, isB);
 *   ...
 *   Slab slab = b.finish();
 * </pre>
 *
 * <p>Not thread-safe. Each builder is owned by a single build thread.
 */
public final class SlabBuilder {
    private final int slabId;
    private final double minMassDa;
    private final double maxMassDa;
    private final List<Fingerprint128> fingerprints = new ArrayList<>();
    private final TreeMap<Integer, List<Integer>> bucketToPeptides = new TreeMap<>();
    private boolean finished;

    public SlabBuilder(int slabId, double minMassDa, double maxMassDa) {
        this.slabId = slabId;
        this.minMassDa = minMassDa;
        this.maxMassDa = maxMassDa;
    }

    public int addPeptide(double precursorMassDa) {
        requireNotFinished();
        int pid = fingerprints.size();
        fingerprints.add(new Fingerprint128());
        return pid;
    }

    public void addFragment(int peptideId, int bucket, boolean isB) {
        requireNotFinished();
        Fingerprint128 fp = fingerprints.get(peptideId);
        if (isB) fp.setBIonBucket(bucket);
        else fp.setYIonBucket(bucket);
        bucketToPeptides.computeIfAbsent(bucket, k -> new ArrayList<>()).add(peptideId);
    }

    public Slab finish() {
        requireNotFinished();
        finished = true;
        int maxBucket = bucketToPeptides.isEmpty() ? 0 : bucketToPeptides.lastKey();
        byte[][] bucketEncoded = new byte[maxBucket + 1][];
        for (var entry : bucketToPeptides.entrySet()) {
            List<Integer> pids = entry.getValue();
            Collections.sort(pids);
            int[] arr = new int[pids.size()];
            for (int i = 0; i < arr.length; i++) arr[i] = pids.get(i);
            bucketEncoded[entry.getKey()] = EliasFano.encode(arr);
        }
        Fingerprint128[] fpArr = fingerprints.toArray(new Fingerprint128[0]);
        return new Slab(slabId, minMassDa, maxMassDa, fpArr, bucketEncoded);
    }

    private void requireNotFinished() {
        if (finished) {
            throw new IllegalStateException("SlabBuilder is single-use; finish() already called");
        }
    }
}
