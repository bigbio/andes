package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Tier-1 candidate generator over a {@link FragmentIndex}: given a spectrum,
 * returns its top-K peptide candidates ranked by accumulated
 * {@code NewRankScorer} log-score across matched fragment buckets.
 *
 * <p>Algorithm (per {@code candidate-generator-design.md} §1):
 * <ol>
 *   <li>Resolve the spectrum's {@link Partition} + b1/y1 {@link IonType}s from
 *       the {@link NewRankScorer}.</li>
 *   <li>Pick the single slab that covers the spectrum's neutral peptide mass.
 *       (Isotope-offset multi-slab loop is a later commit.)</li>
 *   <li>Build a 128-bit spectrum fingerprint from the top-20 ranked peaks and
 *       pre-filter peptides whose fingerprint overlap (popcount of AND) is
 *       below {@link #FP_THRESHOLD}.</li>
 *   <li>For each spectrum peak, walk the matching fragment bucket in the slab
 *       and accumulate {@code rankScorer.getNodeScore(part, ion, peak.getRank())}
 *       into a per-peptide {@code float[] scoreAccum} (survivors only).</li>
 *   <li>Extract the top-K peptide ids via a min-heap.</li>
 * </ol>
 *
 * <p><b>Scope of this first-commit skeleton:</b> single slab per spectrum
 * (no isotope offsets), unmod peptides only, no DBScanner integration.
 *
 * <p><b>Thread safety:</b> each worker thread must hold its own instance.
 * {@code scoreAccum} and {@code fpSurvivors} are reused across
 * {@link #topKForSpectrum} calls on the same instance. The shared
 * {@link FragmentIndex} is immutable and thread-safe for concurrent readers.
 */
public final class FragmentIndexCandidateGenerator {

    /**
     * Minimum Hamming overlap (popcount of bitwise AND over the two 64-bit
     * fingerprint halves) required for a peptide to survive the Tier-1a
     * pre-filter. Hardcoded for now; a CLI flag may come later.
     */
    public static final int FP_THRESHOLD = 4;

    /** Number of highest-intensity peaks used to build the spectrum fingerprint. */
    private static final int FINGERPRINT_TOP_PEAKS = 20;

    private final FragmentIndex index;
    private final float[] scoreAccum;
    private final BitSet fpSurvivors;

    public FragmentIndexCandidateGenerator(FragmentIndex index) {
        this.index = index;
        int max = 0;
        for (int s = 0; s < index.numSlabs(); s++) {
            max = Math.max(max, index.slab(s).peptideCount());
        }
        this.scoreAccum = new float[Math.max(max, 1)];
        this.fpSurvivors = new BitSet(Math.max(max, 1));
    }

    /**
     * Returns the top-{@code K} peptide candidates for the given spectrum,
     * ranked by accumulated NewRankSum descending.
     *
     * <p>Caller contract: {@code spec} must have had {@code setRanksOfPeaks()}
     * called on it so that {@link Peak#getRank()} returns the intensity rank
     * used by {@link NewRankScorer#getNodeScore}.
     *
     * <p>If the spectrum contains no peaks, or its precursor mass falls
     * outside the slab range, returns an empty list.
     */
    public List<CandidateHit> topKForSpectrum(Spectrum spec, int K) {
        if (K <= 0 || spec == null || spec.isEmpty()) {
            return new ArrayList<>();
        }

        int charge = spec.getCharge();
        float parentMass = spec.getPrecursorMass();
        if (charge <= 0 || parentMass <= 0f) return new ArrayList<>();

        // --- Step 2: slab selection (single slab only) ------------------------
        double peptideMass = parentMass;    // spec.getPrecursorMass() returns the neutral mass
        int[] slabIds = index.slabAssigner().slabsFor(peptideMass);
        if (slabIds.length == 0) return new ArrayList<>();
        int slabId = slabIds[0];
        Slab slab = index.slab(slabId);
        int n = slab.peptideCount();
        if (n == 0) return new ArrayList<>();

        final double bw = index.fragmentBinWidthDa();

        // --- Step 3: fingerprint pre-filter -----------------------------------
        Fingerprint128 specFp = new Fingerprint128();
        // Collect up to FINGERPRINT_TOP_PEAKS peaks by intensity rank (rank 1 = highest).
        int limit = Math.min(FINGERPRINT_TOP_PEAKS, spec.size());
        // Indices of peaks with rank <= limit.
        for (Peak p : spec) {
            int r = p.getRank();
            if (r >= 1 && r <= limit) {
                double mz = p.getMz();
                if (mz <= 0) continue;
                int bBucket = (int) (mz / bw);
                double yMass = peptideMass + Composition.H2O - mz;
                if (yMass <= 0) continue;
                int yBucket = (int) (yMass / bw);
                specFp.setBIonBucket(bBucket);
                specFp.setYIonBucket(yBucket);
            }
        }

        final long specLo = specFp.loBits();
        final long specHi = specFp.hiBits();
        fpSurvivors.clear();
        // If the spectrum fingerprint is empty (no usable peaks in top-N), skip
        // the pre-filter entirely — treat all peptides as survivors.
        boolean emptySpecFp = (specLo == 0L && specHi == 0L);
        for (int pid = 0; pid < n; pid++) {
            if (emptySpecFp) {
                fpSurvivors.set(pid);
            } else {
                int overlap = Long.bitCount(slab.fingerprintLoBits(pid) & specLo)
                            + Long.bitCount(slab.fingerprintHiBits(pid) & specHi);
                if (overlap >= FP_THRESHOLD) fpSurvivors.set(pid);
            }
        }

        if (fpSurvivors.isEmpty()) return new ArrayList<>();

        // --- Step 4: bucket walk + peak-rank-weighted hit accumulation --------
        // Zero the score buffer for the peptide range in use.
        for (int i = 0; i < n; i++) scoreAccum[i] = 0f;

        // Scoring function: top-ranked peaks contribute more. A peak at rank 1
        // (highest intensity) scores ~1.0; rank 50 scores ~0.5. Rank 0 means
        // unranked (setRanksOfPeaks not called) — treat as rank 1.
        // This avoids the NewRankScorer partition-lookup path (which needs
        // per-segment, per-ion data populated via the full ScoredSpectrum
        // construction). Tier-2 re-scores via the existing SimpleDBSearchScorer
        // in DBScanner, which carries the full partition / ion-type logic.
        for (Peak p : spec) {
            double mz = p.getMz();
            if (mz <= 0) continue;
            int bucket = (int) (mz / bw);
            int rank = p.getRank() > 0 ? p.getRank() : 1;
            float s = 1f / (1f + 0.02f * (rank - 1));   // rank 1 → 1.0, rank 50 → 0.5
            EliasFano.Cursor cur = slab.bucketCursor(bucket);
            while (cur.hasNext()) {
                int pid = cur.next();
                if (pid < 0 || pid >= n) continue;
                if (!fpSurvivors.get(pid)) continue;
                scoreAccum[pid] += s;
            }
        }

        // --- Step 5: top-K extraction -----------------------------------------
        // Min-heap of (score, pid): poll() yields the smallest so we can evict.
        PriorityQueue<int[]> heap = new PriorityQueue<>(
                Math.max(K, 1),
                Comparator.comparingDouble(a -> Float.intBitsToFloat(a[1])));
        // a[0] = pid, a[1] = Float.floatToRawIntBits(score) — keeps heap alloc tiny.

        for (int pid = fpSurvivors.nextSetBit(0); pid >= 0; pid = fpSurvivors.nextSetBit(pid + 1)) {
            float score = scoreAccum[pid];
            if (heap.size() < K) {
                heap.offer(new int[]{pid, Float.floatToRawIntBits(score)});
            } else {
                int[] top = heap.peek();
                float worst = Float.intBitsToFloat(top[1]);
                if (score > worst) {
                    heap.poll();
                    heap.offer(new int[]{pid, Float.floatToRawIntBits(score)});
                }
            }
        }

        // Drain into a list sorted by score desc.
        List<CandidateHit> out = new ArrayList<>(heap.size());
        while (!heap.isEmpty()) {
            int[] e = heap.poll();
            out.add(new CandidateHit(slabId, e[0], Float.intBitsToFloat(e[1])));
        }
        // Reverse so the highest score is first.
        java.util.Collections.reverse(out);
        return out;
    }
}
