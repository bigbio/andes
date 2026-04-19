package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.fragindex.FragmentIndex;
import edu.ucsd.msjava.fragindex.FragmentIndexBuilder;
import edu.ucsd.msjava.fragindex.SlabAssigner;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import org.junit.Assert;
import org.junit.Test;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.Arrays;

/**
 * Pins the fragment-index dispatch predicate inside {@link DBScanner}.
 *
 * <p>The correctness gate for Phase V2 is: when {@code fragmentIndexMode == OFF}
 * (the default) <em>or</em> when {@code fragmentIndex == null}, the
 * {@code dbSearch} method runs the classic SA-walk body byte-identical to every
 * pre-Phase-3 build. That invariant is encoded in the constructor:
 * {@link DBScanner#candidateGenerator} is non-null if and only if both
 * conditions hold (index present AND mode != OFF).
 *
 * <p>We avoid exercising a full search here — that is the job of the next
 * commit's 3-arm benchmark. This test is scoped to the dispatch predicate
 * only so it stays in unit-test latency and never allocates a suffix array
 * or spectrum preprocessor.
 */
public class TestDBScannerFragmentIndexDispatch {

    /**
     * Builds a tiny, valid {@link FragmentIndex} with one peptide so we can
     * pass something non-null into the constructor. The index itself is
     * never scanned by this test — only the dispatch predicate is examined.
     */
    private static FragmentIndex buildTinyIndex() {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        SlabAssigner assigner = new SlabAssigner(100.0, 4000.0, 50.0, 0.5);
        FragmentIndexBuilder builder = new FragmentIndexBuilder(aaSet, assigner, 1.0005);
        return builder.build(Arrays.asList("PEPTIDER"));
    }

    /**
     * Reflectively resolve {@link DBScanner#candidateGenerator}. Its nullness
     * is the single source of truth for "did the constructor decide this run
     * will hit the fragment-index path?" — so we pin against it here rather
     * than against the method that consumes it.
     */
    private static Object candidateGeneratorField(DBScanner scanner) throws Exception {
        Field f = DBScanner.class.getDeclaredField("candidateGenerator");
        f.setAccessible(true);
        return f.get(scanner);
    }

    /**
     * Build a minimally-viable DBScanner via the Phase-3 constructor using
     * only reflection-friendly nulls for inputs that the dispatch predicate
     * does not read. The legacy constructor delegates through the same path,
     * so exercising the Phase-3 ctor covers both.
     */
    private static DBScanner newScanner(FragmentIndex idx,
                                        SearchParams.FragmentIndexMode mode) throws Exception {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        Enzyme enzyme = Enzyme.TRYPSIN;

        // A ScoredSpectraMap and a CompactSuffixArray are required fields on
        // DBScanner. Neither is exercised by the dispatch predicate — only
        // the ctor itself reads sa.getSize(). We reflectively build a
        // minimal CompactSuffixArray stand-in via an Unsafe-free path: grab
        // the declared ctor and allocate with nulls + a small placeholder
        // size. If that proves fragile across JDK vendors, the fallback is
        // to build a real SA from a 3-byte fasta, which is not worth the
        // test-fixture complexity for pinning a constructor predicate.
        //
        // The simpler, version-safe approach used here: build a
        // CompactSuffixArray stand-in by bypassing the normal ctor via the
        // JDK-standard "allocate" pattern using reflection on a private
        // constructor with a minimal size field. Failing cleanly means the
        // JDK hardened field access since; we surface that via the test.
        CompactSuffixArray sa = allocateStubSA();

        // ScoredSpectraMap is held as a final field but is only dereferenced
        // inside dbSearch body — the dispatch predicate does not touch it,
        // so null is safe for the constructor + field inspection below.
        // However, DBScanner stores specScanner as `final` — a null would
        // be stored and the dispatch check would still work.
        ScoredSpectraMap specScanner = null;

        Constructor<DBScanner> ctor = DBScanner.class.getDeclaredConstructor(
                ScoredSpectraMap.class,
                CompactSuffixArray.class,
                Enzyme.class,
                AminoAcidSet.class,
                int.class, int.class, int.class, int.class, int.class,
                boolean.class, int.class,
                FragmentIndex.class,
                SearchParams.FragmentIndexMode.class);
        return ctor.newInstance(
                specScanner, sa, enzyme, aaSet,
                1, 6, 40, 128, 0,
                false, -1,
                idx, mode);
    }

    /**
     * Produce a {@link CompactSuffixArray} that satisfies only the two reads
     * the DBScanner ctor performs (getSize + storing the reference). The
     * suffix-array machinery itself is never invoked by the dispatch test.
     */
    private static CompactSuffixArray allocateStubSA() throws Exception {
        // Walk every declared ctor; pick the one with a size-only-like
        // signature or allocate via sun.misc.Unsafe only if needed. Simpler:
        // try to find a no-arg ctor; if none, skip the stub and accept that
        // the legacy ctor branch test is not runnable.
        for (Constructor<?> c : CompactSuffixArray.class.getDeclaredConstructors()) {
            if (c.getParameterCount() == 0) {
                c.setAccessible(true);
                return (CompactSuffixArray) c.newInstance();
            }
        }
        // Fallback: use JDK's Unsafe to allocate a zero-initialized instance
        // (bypassing the ctor entirely). Java 17 still exposes
        // sun.misc.Unsafe via its module; the test tolerates failure here by
        // throwing — this marks the constructor-only tests unrunnable if
        // Unsafe is removed in a future JDK, at which point the stub needs
        // to be replaced with a real tiny SA.
        Class<?> unsafeCls = Class.forName("sun.misc.Unsafe");
        Field theUnsafe = unsafeCls.getDeclaredField("theUnsafe");
        theUnsafe.setAccessible(true);
        Object unsafe = theUnsafe.get(null);
        java.lang.reflect.Method allocateInstance =
                unsafeCls.getMethod("allocateInstance", Class.class);
        return (CompactSuffixArray) allocateInstance.invoke(unsafe, CompactSuffixArray.class);
    }

    @Test
    public void offModeWithNullIndexDoesNotAllocateGenerator() throws Exception {
        DBScanner scanner = newScanner(null, SearchParams.FragmentIndexMode.OFF);
        Assert.assertNull(
                "fragmentIndex=null + mode=OFF must leave candidateGenerator null",
                candidateGeneratorField(scanner));
    }

    /**
     * Hardest-to-break-later invariant: even if an index IS available but the
     * user explicitly opted out via {@code -useFragmentIndex off}, the scanner
     * must behave exactly as though no index had been passed. A regression
     * here would silently turn OFF-mode into ON-mode and break bit-identity.
     */
    @Test
    public void offModeWithNonNullIndexStillGoesThroughSaPath() throws Exception {
        FragmentIndex idx = buildTinyIndex();
        DBScanner scanner = newScanner(idx, SearchParams.FragmentIndexMode.OFF);
        Assert.assertNull(
                "fragmentIndex present but mode=OFF must still disable the fragment-index path",
                candidateGeneratorField(scanner));
    }

    @Test
    public void onModeWithNullIndexStillGoesThroughSaPath() throws Exception {
        // Defensive: if the pipeline forgets to build the index (e.g. a
        // regression in MSGFPlus.runMSGFPlus's build block) but the CLI still
        // says ON, we must NOT allocate a generator pointing at null.
        DBScanner scanner = newScanner(null, SearchParams.FragmentIndexMode.ON);
        Assert.assertNull(
                "mode=ON but fragmentIndex=null must still disable the fragment-index path",
                candidateGeneratorField(scanner));
    }

    @Test
    public void onModeWithIndexAllocatesGenerator() throws Exception {
        FragmentIndex idx = buildTinyIndex();
        DBScanner scanner = newScanner(idx, SearchParams.FragmentIndexMode.ON);
        Assert.assertNotNull(
                "fragmentIndex present + mode=ON must enable the fragment-index dispatch path",
                candidateGeneratorField(scanner));
    }

    @Test
    public void compareModeWithIndexAllocatesGenerator() throws Exception {
        FragmentIndex idx = buildTinyIndex();
        DBScanner scanner = newScanner(idx, SearchParams.FragmentIndexMode.COMPARE);
        Assert.assertNotNull(
                "mode=COMPARE with index must also allocate the generator",
                candidateGeneratorField(scanner));
    }

    @Test
    public void legacyConstructorDefaultsToOffPath() throws Exception {
        // The 11-arg legacy constructor is used by MSGFDB and must default
        // to the SA-walk path — it has no knowledge of fragment-index wiring.
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        CompactSuffixArray sa = allocateStubSA();
        Constructor<DBScanner> ctor = DBScanner.class.getDeclaredConstructor(
                ScoredSpectraMap.class,
                CompactSuffixArray.class,
                Enzyme.class,
                AminoAcidSet.class,
                int.class, int.class, int.class, int.class, int.class,
                boolean.class, int.class);
        DBScanner scanner = ctor.newInstance(
                null, sa, Enzyme.TRYPSIN, aaSet,
                1, 6, 40, 128, 0,
                false, -1);
        Assert.assertNull(
                "legacy 11-arg ctor must default to the OFF path (candidateGenerator null)",
                candidateGeneratorField(scanner));
    }
}
