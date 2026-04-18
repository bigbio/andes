package edu.ucsd.msjava.mzid;

import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msdbsearch.DatabaseMatch;
import edu.ucsd.msjava.msdbsearch.MSGFPlusMatch;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Modification;
import edu.ucsd.msjava.msutil.ModifiedAminoAcid;
import edu.ucsd.msjava.msutil.Pair;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;

/**
 * Writes MS-GF+ search results in Percolator {@code .pin} format, bypassing
 * the external {@code msgf2pin} converter. Emitted file is directly usable
 * by Percolator (<a href="https://github.com/percolator/percolator">percolator</a>)
 * and downstream MS²Rescore / Mokapot pipelines.
 *
 * <p>Column layout (tab-separated, header on first line):
 * <pre>
 *   SpecId  Label  ScanNr  ExpMass  CalcMass
 *   RawScore  DeNovoScore  lnSpecEValue  lnEValue  IsotopeError
 *   PepLen  dM  absdM
 *   Charge2 … ChargeK          (one-hot over params.getMinCharge()..params.getMaxCharge())
 *   NumMatchedMainIons  ExplainedIonCurrentRatio  NTermIonCurrentRatio
 *   CTermIonCurrentRatio  MS2IonCurrent  IsolationWindowEfficiency
 *   MeanErrorTop7  StdevErrorTop7  MeanRelErrorTop7  StdevRelErrorTop7
 *   lnDeltaSpecEValue  matchedIonRatio
 *   Peptide  Proteins
 * </pre>
 *
 * <p>{@code Label} is {@code 1} when at least one protein match is not a decoy,
 * {@code -1} when every match for the PSM is a decoy. PSMs with no real protein
 * are written with Label = -1 so Percolator can use them for the null
 * distribution.
 *
 * <p>The per-match additional-feature columns (rows 8-17 above) are zero-filled
 * when {@code -addFeatures 1} is not supplied — so the column count is stable
 * across runs. Downstream config files that reference the feature column index
 * therefore work regardless of whether the upstream search used {@code -addFeatures 1}.
 */
public class DirectPinWriter {

    private final SearchParams params;
    private final AminoAcidSet aaSet;
    private final CompactSuffixArray sa;
    private final SpectraAccessor specAcc;
    private final String decoyProteinPrefix;
    private final Map<String, List<Double>> fixedModMasses;

    /** Feature names sourced from {@code Match.getAdditionalFeatureList()}, in stable order. */
    private static final String[] PIN_FEATURES = {
            "NumMatchedMainIons",
            "ExplainedIonCurrentRatio", "NTermIonCurrentRatio", "CTermIonCurrentRatio",
            "MS2IonCurrent", "IsolationWindowEfficiency",
            "MeanErrorTop7", "StdevErrorTop7", "MeanRelErrorTop7", "StdevRelErrorTop7"
    };

    /**
     * Extra PSM-level features computed here (not sourced from the match list):
     *  - lnDeltaSpecEValue: log(rank1 SpecEValue / rank2 SpecEValue) for rank-1 PSMs; 0 otherwise.
     *  - matchedIonRatio:   NumMatchedMainIons / PepLen.
     */
    private static final String[] PIN_EXTRA_FEATURES = {
            "lnDeltaSpecEValue", "matchedIonRatio"
    };

    public DirectPinWriter(SearchParams params, AminoAcidSet aaSet,
                           CompactSuffixArray sa, SpectraAccessor specAcc, int ioIndex) {
        this.params = params;
        this.aaSet = aaSet;
        this.sa = sa;
        this.specAcc = specAcc;
        this.decoyProteinPrefix = params.getDecoyProteinPrefix();
        this.fixedModMasses = buildFixedModMap(aaSet);
        // ioIndex accepted for API symmetry with DirectTSVWriter; not
        // currently referenced but reserved for per-file logging later.
    }

    public void writeResults(List<MSGFPlusMatch> resultList, File outputFile) throws IOException {
        int minCharge = params.getMinCharge();
        int maxCharge = params.getMaxCharge();

        try (PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile), 256 * 1024))) {
            writeHeader(out, minCharge, maxCharge);

            for (MSGFPlusMatch mpMatch : resultList) {
                int specIndex = mpMatch.getSpecIndex();
                List<DatabaseMatch> matchList = mpMatch.getMatchList();
                if (matchList == null || matchList.isEmpty()) continue;

                Spectrum spec = specAcc.getSpecMap().getSpectrumBySpecIndex(specIndex);
                if (spec == null) continue;

                String specID = spec.getID();
                int scanNum = spec.getScanNum();
                float precursorMz = spec.getPrecursorPeak().getMz();

                double rank2SpecEValue = findRank2SpecEValue(matchList, params.getMinDeNovoScore());

                int rank = 0;
                double prevSpecEValue = Double.NaN;
                for (int i = matchList.size() - 1; i >= 0; --i) {
                    DatabaseMatch match = matchList.get(i);
                    if (match.getDeNovoScore() < params.getMinDeNovoScore()) continue;

                    if (match.getSpecEValue() != prevSpecEValue) ++rank;
                    prevSpecEValue = match.getSpecEValue();

                    writeRow(out, specID, scanNum, rank, precursorMz, match, minCharge, maxCharge,
                            rank2SpecEValue);
                }
            }
        }
    }

    private void writeHeader(PrintStream out, int minCharge, int maxCharge) {
        StringBuilder h = new StringBuilder(256);
        h.append("SpecId\tLabel\tScanNr\tExpMass\tCalcMass")
                .append("\tRawScore\tDeNovoScore\tlnSpecEValue\tlnEValue\tIsotopeError")
                .append("\tPepLen\tdM\tabsdM");
        for (int c = minCharge; c <= maxCharge; c++) {
            h.append("\tCharge").append(c);
        }
        for (String f : PIN_FEATURES) h.append('\t').append(f);
        for (String f : PIN_EXTRA_FEATURES) h.append('\t').append(f);
        h.append("\tPeptide\tProteins");
        out.println(h);
    }

    private void writeRow(PrintStream out, String specID, int scanNum, int rank,
                          float precursorMz, DatabaseMatch match, int minCharge, int maxCharge,
                          double rank2SpecEValue) {
        int length = match.getLength();
        int charge = match.getCharge();
        float peptideMass = match.getPeptideMass();
        float theoMz = (peptideMass + (float) Composition.H2O) / charge + (float) Composition.ChargeCarrierMass();

        double specEValue = match.getSpecEValue();
        int numPeptides = sa.getNumDistinctPeptides(params.getEnzyme() == null ? length - 2 : length - 1);
        double eValue = specEValue * numPeptides;

        float expMass = precursorMz * charge;
        float theoMass = theoMz * charge;
        int isotopeError = Math.round((expMass - theoMass) / (float) Composition.ISOTOPE);
        double adjustedExpMz = precursorMz - Composition.ISOTOPE * isotopeError / charge;
        double dM = adjustedExpMz - theoMz;

        String peptideSeq = formatPeptideWithMods(match.getPepSeq());
        ProteinFormatResult proteins = formatProteinsForPin(match, length);

        // Drop all-decoy matches? Percolator prefers to see them with Label=-1.
        int label = proteins.allDecoy ? -1 : 1;

        String psmId = specID + "_" + scanNum + "_" + rank;
        Map<String, String> features = collectFeatures(match);

        StringBuilder row = new StringBuilder(512);
        row.append(psmId)
                .append('\t').append(label)
                .append('\t').append(scanNum)
                .append('\t').append(formatDouble(expMass))
                .append('\t').append(formatDouble(theoMass))
                .append('\t').append(match.getScore())
                .append('\t').append(match.getDeNovoScore())
                .append('\t').append(formatDouble(specEValue > 0 ? Math.log(specEValue) : -Double.MAX_VALUE))
                .append('\t').append(formatDouble(eValue > 0 ? Math.log(eValue) : -Double.MAX_VALUE))
                .append('\t').append(isotopeError)
                .append('\t').append(length)
                .append('\t').append(formatDouble(dM))
                .append('\t').append(formatDouble(Math.abs(dM)));
        for (int c = minCharge; c <= maxCharge; c++) {
            row.append('\t').append(c == charge ? 1 : 0);
        }
        for (String f : PIN_FEATURES) {
            String v = features.get(f);
            row.append('\t').append(v == null || v.isEmpty() ? "0" : v);
        }
        double lnDeltaSpecEValue = computeLnDeltaSpecEValue(rank, specEValue, rank2SpecEValue);
        double matchedIonRatio = computeMatchedIonRatio(features.get("NumMatchedMainIons"), length);
        row.append('\t').append(formatDouble(lnDeltaSpecEValue))
                .append('\t').append(formatDouble(matchedIonRatio));
        // Peptide in Percolator "flanking.PEPTIDE.flanking" format.
        row.append('\t').append(proteins.pre).append('.').append(peptideSeq).append('.').append(proteins.post);
        for (String acc : proteins.accessions) row.append('\t').append(acc);
        out.println(row);
    }

    private static String formatDouble(double v) {
        if (Double.isNaN(v) || Double.isInfinite(v)) return "0";
        // Percolator is fine with plain scientific or fixed notation.
        return String.format(Locale.ROOT, "%.6g", v);
    }

    private static Map<String, String> collectFeatures(DatabaseMatch match) {
        Map<String, String> m = new HashMap<>();
        List<Pair<String, String>> featureList = match.getAdditionalFeatureList();
        if (featureList != null) {
            for (Pair<String, String> p : featureList) m.put(p.getFirst(), p.getSecond());
        }
        return m;
    }

    /**
     * Scans the match list (ordered worst-to-best like {@code writeResults}) and returns the
     * SpecEValue of the rank-2 PSM: the first distinct SpecEValue encountered after the
     * rank-1 value, skipping duplicates (ties share a rank) and matches below
     * {@code minDeNovoScore}. Returns {@link Double#NaN} if no rank-2 exists.
     */
    public static double findRank2SpecEValue(List<DatabaseMatch> matchList, int minDeNovoScore) {
        double rank1 = Double.NaN;
        for (int i = matchList.size() - 1; i >= 0; --i) {
            DatabaseMatch m = matchList.get(i);
            if (m.getDeNovoScore() < minDeNovoScore) continue;
            double se = m.getSpecEValue();
            if (Double.isNaN(rank1)) {
                rank1 = se;
            } else if (se != rank1) {
                return se;
            }
        }
        return Double.NaN;
    }

    /**
     * {@code log(rank1 SpecEValue / rank2 SpecEValue)} for rank-1 PSMs; {@code 0} otherwise
     * or when either SpecEValue is non-positive / NaN. Larger (more negative) values mean
     * the top hit is more separated from the next best, which Percolator / MS²Rescore /
     * Mokapot can exploit for rescoring.
     */
    public static double computeLnDeltaSpecEValue(int rank, double rank1SpecEValue, double rank2SpecEValue) {
        if (rank != 1) return 0.0;
        if (Double.isNaN(rank1SpecEValue) || Double.isNaN(rank2SpecEValue)) return 0.0;
        if (rank1SpecEValue <= 0 || rank2SpecEValue <= 0) return 0.0;
        return Math.log(rank1SpecEValue / rank2SpecEValue);
    }

    /** {@code NumMatchedMainIons / PepLen}: peptide-length-normalized ion-match density. */
    public static double computeMatchedIonRatio(String numMatchedMainIons, int pepLen) {
        if (pepLen <= 0) return 0.0;
        if (numMatchedMainIons == null || numMatchedMainIons.isEmpty()) return 0.0;
        try {
            double n = Double.parseDouble(numMatchedMainIons);
            return n / pepLen;
        } catch (NumberFormatException e) {
            return 0.0;
        }
    }

    // -----------------------------------------------------------------------
    // Protein flanking + decoy resolution (Percolator-specific)
    // -----------------------------------------------------------------------

    /** Flanking residues + accession list resolved from the suffix array. */
    private static final class ProteinFormatResult {
        char pre = '-';
        char post = '-';
        boolean allDecoy = true;
        List<String> accessions = new ArrayList<>();
    }

    private ProteinFormatResult formatProteinsForPin(DatabaseMatch match, int length) {
        ProteinFormatResult res = new ProteinFormatResult();
        SortedSet<Integer> indices = match.getIndices();
        CompactFastaSequence seq = sa.getSequence();
        HashSet<String> seen = new HashSet<>();

        boolean firstRealCaptured = false;
        for (int index : indices) {
            boolean isNTermMetCleaved = false;
            if (seq.getByteAt(index) == 0 && seq.getCharAt(index + 1) == 'M') {
                Peptide peptide = aaSet.getPeptide(match.getPepSeq());
                StringBuilder pepUnmod = new StringBuilder();
                for (AminoAcid aa : peptide) pepUnmod.append(aa.getUnmodResidue());
                String pepSeqStr = pepUnmod.toString();
                isNTermMetCleaved = match.isNTermMetCleaved() || pepSeqStr.charAt(0) != 'M';
                if (!isNTermMetCleaved) {
                    String matchSequence = seq.getSubsequence(index + 2, index + 3 + pepSeqStr.length());
                    isNTermMetCleaved = matchSequence.startsWith(pepSeqStr);
                }
            }

            char pre = seq.getCharAt(index);
            if (pre == '_') pre = isNTermMetCleaved ? 'M' : '-';
            char post = isNTermMetCleaved ? seq.getCharAt(index + length) : seq.getCharAt(index + length - 1);
            if (post == '_') post = '-';

            int protStart = (int) seq.getStartPosition(index);
            String annotation = seq.getAnnotation(protStart);
            String accession = annotation.split("\\s+")[0];

            boolean isDecoy = accession.startsWith(decoyProteinPrefix);
            if (!isDecoy) res.allDecoy = false;

            if (!seen.add(accession)) continue;
            res.accessions.add(accession);

            // Capture pre/post from the first non-decoy occurrence; fall back to the
            // first entry if every match is a decoy.
            if (!firstRealCaptured && !isDecoy) {
                res.pre = pre;
                res.post = post;
                firstRealCaptured = true;
            } else if (!firstRealCaptured && res.accessions.size() == 1) {
                res.pre = pre;
                res.post = post;
            }
        }
        return res;
    }

    // -----------------------------------------------------------------------
    // Peptide formatting — duplicated from DirectTSVWriter. Both should move
    // to a shared PeptideFormatter in a follow-up.
    // -----------------------------------------------------------------------

    private static Map<String, List<Double>> buildFixedModMap(AminoAcidSet aaSet) {
        Map<String, List<Double>> m = new HashMap<>();
        for (Modification.Instance mod : aaSet.getModifications()) {
            if (mod.isFixedModification()) {
                String key = modKey(mod.getResidue(), mod.getLocation());
                List<Double> list = m.get(key);
                if (list == null) { list = new ArrayList<>(); m.put(key, list); }
                list.add(mod.getModification().getAccurateMass());
            }
        }
        return m;
    }

    private static String modKey(char residue, Modification.Location location) {
        switch (location) {
            case N_Term:
            case Protein_N_Term:
                return "[" + residue;
            case C_Term:
            case Protein_C_Term:
                return residue + "]";
            default:
                return String.valueOf(residue);
        }
    }

    private String formatPeptideWithMods(String pepSeq) {
        Peptide peptide = aaSet.getPeptide(pepSeq);
        StringBuilder unmodSeq = new StringBuilder();
        String[] modArr = new String[peptide.size() + 2];

        int location = 1;
        for (AminoAcid aa : peptide) {
            unmodSeq.append(aa.getUnmodResidue());
            if (aa.isModified()) {
                ModifiedAminoAcid modAA = (ModifiedAminoAcid) aa;
                int modLoc = resolveModLocation(modAA, location, peptide.size());
                appendMassStr(modArr, modLoc, modAA.getModification().getAccurateMass());
                while (modAA.getTargetAA().isModified()) {
                    modAA = (ModifiedAminoAcid) modAA.getTargetAA();
                    int stackLoc = resolveModLocation(modAA, location, peptide.size());
                    appendMassStr(modArr, stackLoc, modAA.getModification().getAccurateMass());
                }
            }
            List<Double> fixedResMods = fixedModMasses.get(String.valueOf(aa.getUnmodResidue()));
            if (fixedResMods != null) {
                for (double mass : fixedResMods) appendMassStr(modArr, location, mass);
            }
            if (location == 1) appendTerminalFixedMods(modArr, 0, aa.getUnmodResidue(), "[");
            if (location == peptide.size()) appendTerminalFixedMods(modArr, peptide.size() + 1, aa.getUnmodResidue(), "]");
            location++;
        }

        StringBuilder buf = new StringBuilder();
        if (modArr[0] != null) buf.append(modArr[0]);
        for (int i = 0; i < unmodSeq.length(); i++) {
            buf.append(unmodSeq.charAt(i));
            if (modArr[i + 1] != null) buf.append(modArr[i + 1]);
        }
        if (modArr[modArr.length - 1] != null) buf.append(modArr[modArr.length - 1]);
        return buf.toString();
    }

    private static int resolveModLocation(ModifiedAminoAcid modAA, int location, int pepLen) {
        if (location == 1 && modAA.isNTermVariableMod()) return 0;
        if (location == pepLen && modAA.isCTermVariableMod()) return pepLen + 1;
        return location;
    }

    private static void appendMassStr(String[] modArr, int loc, double mass) {
        String str = mass >= 0 ? "+" + String.format(Locale.ROOT, "%.3f", mass)
                               : String.format(Locale.ROOT, "%.3f", mass);
        modArr[loc] = (modArr[loc] == null) ? str : modArr[loc] + str;
    }

    private void appendTerminalFixedMods(String[] modArr, int loc, char residue, String bracket) {
        String keyRes = bracket.equals("[") ? "[" + residue : residue + "]";
        List<Double> mods1 = fixedModMasses.get(keyRes);
        if (mods1 != null) for (double m : mods1) appendMassStr(modArr, loc, m);
        String keyAny = bracket.equals("[") ? "[*" : "*]";
        List<Double> mods2 = fixedModMasses.get(keyAny);
        if (mods2 != null) for (double m : mods2) appendMassStr(modArr, loc, m);
    }
}
