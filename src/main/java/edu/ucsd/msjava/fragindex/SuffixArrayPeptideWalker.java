package edu.ucsd.msjava.fragindex;

import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msutil.Enzyme;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Walks a {@link CompactFastaSequence} and emits every tryptic (enzyme-specific)
 * peptide string whose length is in {@code [minLen, maxLen]} and that belongs
 * to a target (non-decoy) protein.
 *
 * <p>Protein boundaries in the flat sequence buffer count as cleavage sites on
 * both ends. Between a pair of cleavage sites the walker considers every
 * candidate peptide that spans at most {@code maxMissedCleavages} internal
 * sites. Peptides containing residues outside the sequence's configured
 * alphabet (e.g. terminators, invalid codes, or stop-codon symbols such as
 * {@code '*'}) are silently skipped.
 *
 * <p>Decoy-prefixed proteins are skipped entirely at walk time; the resulting
 * peptide stream is target-only. (Search-time scoring still uses the full
 * target+decoy database; decoys are handled downstream via the suffix array.)
 */
public final class SuffixArrayPeptideWalker {

    private final CompactFastaSequence sequence;
    private final Enzyme enzyme;
    private final int minLen;
    private final int maxLen;
    private final int maxMissedCleavages;

    public SuffixArrayPeptideWalker(CompactFastaSequence sequence,
                                     Enzyme enzyme,
                                     int minLen,
                                     int maxLen,
                                     int maxMissedCleavages) {
        if (sequence == null) throw new IllegalArgumentException("sequence must not be null");
        if (enzyme == null) throw new IllegalArgumentException("enzyme must not be null");
        if (minLen < 1) throw new IllegalArgumentException("minLen must be >= 1, got " + minLen);
        if (maxLen < minLen) throw new IllegalArgumentException("maxLen (" + maxLen + ") must be >= minLen (" + minLen + ")");
        if (maxMissedCleavages < 0) throw new IllegalArgumentException("maxMissedCleavages must be >= 0, got " + maxMissedCleavages);

        this.sequence = sequence;
        this.enzyme = enzyme;
        this.minLen = minLen;
        this.maxLen = maxLen;
        this.maxMissedCleavages = maxMissedCleavages;
    }

    /**
     * Iterate every qualifying peptide in the sequence. The callback is invoked
     * synchronously once per peptide in walk order.
     */
    public void forEachPeptide(Consumer<String> consumer) {
        if (consumer == null) throw new IllegalArgumentException("consumer must not be null");

        final int size = (int) sequence.getSize();
        if (size <= 0) return;

        final String decoyPrefix = sequence.getDecoyProteinPrefix();
        final boolean hasDecoyPrefix = decoyPrefix != null && !decoyPrefix.isEmpty();
        final boolean enzymeIsNTerm = enzyme.isNTerm();
        final boolean enzymeIsCTerm = enzyme.isCTerm();

        int proteinStart = -1;
        List<Integer> sites = null;
        boolean skipProtein = false;

        for (int pos = 0; pos < size; pos++) {
            boolean isTerm = sequence.isTerminator(pos);

            if (proteinStart < 0) {
                // Looking for the start of a protein.
                if (isTerm) {
                    continue;
                }
                // First residue of a new protein.
                proteinStart = pos;
                skipProtein = false;
                sites = null;

                if (hasDecoyPrefix) {
                    String annotation = sequence.getAnnotation(proteinStart);
                    if (annotation != null) {
                        String accession = firstToken(annotation);
                        if (accession.startsWith(decoyPrefix)) {
                            skipProtein = true;
                        }
                    }
                }

                if (!skipProtein) {
                    sites = new ArrayList<>();
                    sites.add(proteinStart); // protein N-terminus
                }
            }

            if (isTerm) {
                // End of current protein. Emit and reset.
                if (!skipProtein && sites != null) {
                    int last = sites.get(sites.size() - 1);
                    if (last != pos) sites.add(pos);
                    emitPeptides(sites, consumer);
                }
                proteinStart = -1;
                sites = null;
                skipProtein = false;
                continue;
            }

            // Residue byte inside a protein.
            if (skipProtein) continue;

            char residue = sequence.getCharAt(pos);

            if (enzymeIsCTerm && enzyme.isCleavable(residue)) {
                // Cleavage AFTER residue; site = pos + 1 (exclusive end of a peptide).
                int site = pos + 1;
                int lastSite = sites.get(sites.size() - 1);
                if (site != lastSite) {
                    sites.add(site);
                }
            } else if (enzymeIsNTerm && enzyme.isCleavable(residue)) {
                // Cleavage BEFORE residue; site = pos (inclusive start of a peptide).
                if (pos != proteinStart) {
                    int lastSite = sites.get(sites.size() - 1);
                    if (pos != lastSite) {
                        sites.add(pos);
                    }
                }
            }
        }

        // Defensive tail handling: if the buffer did not end with a terminator,
        // still close out the last protein.
        if (proteinStart >= 0 && !skipProtein && sites != null) {
            int last = sites.get(sites.size() - 1);
            if (last != size) sites.add(size);
            emitPeptides(sites, consumer);
        }
    }

    /**
     * Same as {@link #forEachPeptide} but materialises every peptide into a
     * list. Do not call on production-size FASTAs; intended for tests and
     * small inputs.
     */
    public List<String> collect() {
        List<String> out = new ArrayList<>();
        forEachPeptide(out::add);
        return out;
    }

    // --- internals -----------------------------------------------------------

    private void emitPeptides(List<Integer> sites, Consumer<String> consumer) {
        final int n = sites.size();
        for (int i = 0; i < n - 1; i++) {
            int start = sites.get(i);
            int maxJ = Math.min(n - 1, i + 1 + maxMissedCleavages);
            for (int j = i + 1; j <= maxJ; j++) {
                int end = sites.get(j);
                int len = end - start;
                if (len > maxLen) break; // further j only grows
                if (len < minLen) continue;
                if (!allInAlphabet(start, end)) continue;

                String peptide = sequence.getSubsequence(start, end);
                if (peptide != null) {
                    consumer.accept(peptide);
                }
            }
        }
    }

    private boolean allInAlphabet(int start, int end) {
        for (int p = start; p < end; p++) {
            char c = sequence.getCharAt(p);
            if (!sequence.isInAlphabet(c)) {
                return false;
            }
        }
        return true;
    }

    private static String firstToken(String annotation) {
        for (int i = 0; i < annotation.length(); i++) {
            char c = annotation.charAt(i);
            if (Character.isWhitespace(c)) {
                return annotation.substring(0, i);
            }
        }
        return annotation;
    }
}
