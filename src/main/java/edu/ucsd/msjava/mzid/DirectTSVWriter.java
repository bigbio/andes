package edu.ucsd.msjava.mzid;

import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msdbsearch.DatabaseMatch;
import edu.ucsd.msjava.msdbsearch.MSGFPlusMatch;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;

import java.io.*;
import java.util.*;

/**
 * Writes MS-GF+ search results directly to TSV format from in-memory objects,
 * bypassing mzIdentML serialization. Output is column-compatible with MzIDToTsv
 * so that OpenMS MSGFPlusAdapter can consume it without changes.
 */
public class DirectTSVWriter {

    private final SearchParams params;
    private final AminoAcidSet aaSet;
    private final CompactSuffixArray sa;
    private final SpectraAccessor specAcc;
    private final int ioIndex;
    private final boolean isPrecursorTolerancePPM;
    private final String decoyProteinPrefix;
    private final boolean isMgf;

    // Fixed mod map: residue key -> list of modification masses
    // Keys: "C" for residue-specific, "[C" or "[*" for N-term, "C]" or "*]" for C-term
    private final Map<String, List<Double>> fixedModMasses;

    public DirectTSVWriter(SearchParams params, AminoAcidSet aaSet,
                           CompactSuffixArray sa, SpectraAccessor specAcc, int ioIndex) {
        this.params = params;
        this.aaSet = aaSet;
        this.sa = sa;
        this.specAcc = specAcc;
        this.ioIndex = ioIndex;
        this.isPrecursorTolerancePPM = params.getRightPrecursorMassTolerance().isTolerancePPM();
        this.decoyProteinPrefix = params.getDecoyProteinPrefix();

        SpecFileFormat fmt = params.getDBSearchIOList().get(ioIndex).getSpecFileFormat();
        this.isMgf = (fmt == SpecFileFormat.MGF);

        // Build fixed modification mass map from AminoAcidSet
        this.fixedModMasses = new HashMap<>();
        for (Modification.Instance mod : aaSet.getModifications()) {
            if (mod.isFixedModification()) {
                String key = getModKey(mod.getResidue(), mod.getLocation());
                fixedModMasses.computeIfAbsent(key, k -> new ArrayList<>())
                        .add(mod.getModification().getAccurateMass());
            }
        }
    }

    private static String getModKey(char residue, Modification.Location location) {
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

    public void writeResults(List<MSGFPlusMatch> resultList, File outputFile) throws IOException {
        String specFileName = params.getDBSearchIOList().get(ioIndex).getSpecFile().getName();
        boolean showQValue = params.useTDA();

        try (PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile), 256 * 1024))) {
            // Header
            out.println("#SpecFile"
                    + "\tSpecID"
                    + "\tScanNum"
                    + (isMgf ? "\tTitle" : "")
                    + "\tFragMethod"
                    + "\tPrecursor"
                    + "\tIsotopeError"
                    + "\tPrecursorError(" + (isPrecursorTolerancePPM ? "ppm" : "Da") + ")"
                    + "\tCharge"
                    + "\tPeptide"
                    + "\tProtein"
                    + "\tDeNovoScore"
                    + "\tMSGFScore"
                    + "\tSpecEValue"
                    + "\tEValue"
                    + (showQValue ? "\tQValue\tPepQValue" : ""));

            for (MSGFPlusMatch mpMatch : resultList) {
                int specIndex = mpMatch.getSpecIndex();
                List<DatabaseMatch> matchList = mpMatch.getMatchList();
                if (matchList == null || matchList.isEmpty())
                    continue;

                Spectrum spec = specAcc.getSpecMap().getSpectrumBySpecIndex(specIndex);
                if (spec == null) continue;

                String specID = spec.getID();
                int scanNum = spec.getScanNum();
                float precursorMz = spec.getPrecursorPeak().getMz();
                String title = isMgf ? spec.getTitle() : null;

                int rank = 0;
                double prevSpecEValue = Double.NaN;
                for (int i = matchList.size() - 1; i >= 0; --i) {
                    DatabaseMatch match = matchList.get(i);

                    if (match.getDeNovoScore() < params.getMinDeNovoScore())
                        continue;

                    int length = match.getLength();
                    int charge = match.getCharge();
                    float peptideMass = match.getPeptideMass();
                    float theoMz = (peptideMass + (float) Composition.H2O) / charge + (float) Composition.ChargeCarrierMass();

                    int score = match.getScore();
                    double specEValue = match.getSpecEValue();
                    int numPeptides = sa.getNumDistinctPeptides(params.getEnzyme() == null ? length - 2 : length - 1);
                    double eValue = specEValue * numPeptides;

                    if (prevSpecEValue != specEValue) ++rank;
                    prevSpecEValue = specEValue;

                    String specEValueStr;
                    if (specEValue < Float.MIN_NORMAL)
                        specEValueStr = String.valueOf(specEValue);
                    else
                        specEValueStr = String.valueOf((float) specEValue);

                    String eValueStr;
                    if (specEValue < Float.MIN_NORMAL)
                        eValueStr = String.valueOf(eValue);
                    else
                        eValueStr = String.valueOf((float) eValue);

                    // Isotope error
                    float expMass = precursorMz * charge;
                    float theoMass = theoMz * charge;
                    int isotopeError = Math.round((expMass - theoMass) / (float) Composition.ISOTOPE);

                    // Precursor error
                    double adjustedExpMz = precursorMz - Composition.ISOTOPE * isotopeError / charge;
                    double precursorError = adjustedExpMz - theoMz;
                    if (isPrecursorTolerancePPM)
                        precursorError = precursorError / theoMz * 1e6;

                    // Fragmentation method
                    ActivationMethod[] actMethodArr = match.getActivationMethodArr();
                    String fragMethod = "";
                    if (actMethodArr != null) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(actMethodArr[0]);
                        for (int j = 1; j < actMethodArr.length; j++)
                            sb.append("/").append(actMethodArr[j]);
                        fragMethod = sb.toString();
                    }

                    // Peptide sequence with modifications
                    String peptideSeq = formatPeptideWithMods(match.getPepSeq());

                    // Protein accessions with pre/post
                    String proteinStr = formatProteins(match, length);

                    if (proteinStr.isEmpty()) continue; // all decoy, skip

                    out.print(specFileName
                            + "\t" + specID
                            + "\t" + scanNum
                            + (isMgf ? "\t" + (title != null ? title : "N/A") : "")
                            + "\t" + fragMethod
                            + "\t" + precursorMz
                            + "\t" + isotopeError
                            + "\t" + (float) precursorError
                            + "\t" + charge
                            + "\t" + peptideSeq
                            + "\t" + proteinStr
                            + "\t" + match.getDeNovoScore()
                            + "\t" + score
                            + "\t" + specEValueStr
                            + "\t" + eValueStr
                    );
                    if (showQValue) {
                        Float psmQValue = match.getPSMQValue();
                        Float pepQValue = match.getPepQValue();
                        out.print("\t" + (psmQValue != null ? psmQValue : "")
                                + "\t" + (pepQValue != null ? pepQValue : ""));
                    }
                    out.println();
                }
            }
        }
    }

    /**
     * Format peptide sequence with inline modification masses.
     * Matches the format produced by MzIDParser.getPeptideSeq():
     * e.g. "NLANPTSVILASIQM+15.995LEYLGMADK"
     */
    private String formatPeptideWithMods(String pepSeq) {
        edu.ucsd.msjava.msutil.Peptide peptide = aaSet.getPeptide(pepSeq);
        StringBuilder unmodSeq = new StringBuilder();
        // modArr indexed by location: 0 = N-term, 1..len = residues, len+1 = C-term
        String[] modArr = new String[peptide.size() + 2];

        int location = 1;
        for (AminoAcid aa : peptide) {
            unmodSeq.append(aa.getUnmodResidue());

            if (aa.isModified()) {
                ModifiedAminoAcid modAA = (ModifiedAminoAcid) aa;

                // Determine location for the mod
                int modLoc;
                if (location == 1 && modAA.isNTermVariableMod()) {
                    modLoc = 0; // N-term
                } else if (location == peptide.size() && modAA.isCTermVariableMod()) {
                    modLoc = peptide.size() + 1; // C-term
                } else {
                    modLoc = location;
                }

                double mass = modAA.getModification().getAccurateMass();
                String massStr = mass >= 0 ? "+" + String.format("%.3f", mass) : String.format("%.3f", mass);
                modArr[modLoc] = (modArr[modLoc] == null) ? massStr : modArr[modLoc] + massStr;

                // Handle stacked modifications
                while (modAA.getTargetAA().isModified()) {
                    modAA = (ModifiedAminoAcid) modAA.getTargetAA();
                    int stackModLoc;
                    if (location == 1 && modAA.isNTermVariableMod()) {
                        stackModLoc = 0;
                    } else if (location == peptide.size() && modAA.isCTermVariableMod()) {
                        stackModLoc = peptide.size() + 1;
                    } else {
                        stackModLoc = location;
                    }
                    double stackMass = modAA.getModification().getAccurateMass();
                    String stackMassStr = stackMass >= 0 ? "+" + String.format("%.3f", stackMass) : String.format("%.3f", stackMass);
                    modArr[stackModLoc] = (modArr[stackModLoc] == null) ? stackMassStr : modArr[stackModLoc] + stackMassStr;
                }
            }

            // Fixed modifications (residue-specific)
            List<Double> fixedResideMods = fixedModMasses.get(String.valueOf(aa.getUnmodResidue()));
            if (fixedResideMods != null) {
                for (double mass : fixedResideMods) {
                    String massStr = mass >= 0 ? "+" + String.format("%.3f", mass) : String.format("%.3f", mass);
                    modArr[location] = (modArr[location] == null) ? massStr : modArr[location] + massStr;
                }
            }

            // Fixed terminal modifications
            if (location == 1) {
                addFixedTerminalMods(modArr, 0, aa.getUnmodResidue(), "[");
            }
            if (location == peptide.size()) {
                addFixedTerminalMods(modArr, peptide.size() + 1, aa.getUnmodResidue(), "]");
            }

            location++;
        }

        // Build the modified peptide string
        StringBuilder buf = new StringBuilder();
        if (modArr[0] != null) buf.append(modArr[0]);
        for (int i = 0; i < unmodSeq.length(); i++) {
            buf.append(unmodSeq.charAt(i));
            if (modArr[i + 1] != null) buf.append(modArr[i + 1]);
        }
        if (modArr[modArr.length - 1] != null) buf.append(modArr[modArr.length - 1]);

        return buf.toString();
    }

    private void addFixedTerminalMods(String[] modArr, int loc, char residue, String bracket) {
        // Residue-specific terminal mod (e.g., "[C" for N-term on C)
        String key1 = bracket.equals("[") ? "[" + residue : residue + "]";
        List<Double> mods1 = fixedModMasses.get(key1);
        if (mods1 != null) {
            for (double mass : mods1) {
                String massStr = mass >= 0 ? "+" + String.format("%.3f", mass) : String.format("%.3f", mass);
                modArr[loc] = (modArr[loc] == null) ? massStr : modArr[loc] + massStr;
            }
        }
        // Wildcard terminal mod (e.g., "[*" for N-term on any residue)
        String key2 = bracket.equals("[") ? "[*" : "*]";
        List<Double> mods2 = fixedModMasses.get(key2);
        if (mods2 != null) {
            for (double mass : mods2) {
                String massStr = mass >= 0 ? "+" + String.format("%.3f", mass) : String.format("%.3f", mass);
                modArr[loc] = (modArr[loc] == null) ? massStr : modArr[loc] + massStr;
            }
        }
    }

    /**
     * Format protein accessions in merged mode:
     * "accession1(pre=X,post=Y);accession2(pre=X,post=Y)"
     * Mirrors MzIDParser merged-mode protein formatting.
     */
    private String formatProteins(DatabaseMatch match, int length) {
        SortedSet<Integer> indices = match.getIndices();
        CompactFastaSequence seq = sa.getSequence();
        StringBuilder proteinBuf = new StringBuilder();
        HashSet<String> proteinSet = new HashSet<>();
        boolean isAllDecoy = true;

        for (int index : indices) {
            boolean isNTermMetCleaved = false;

            // Check for N-terminal Met cleavage (same logic as MZIdentMLGen)
            if (seq.getByteAt(index) == 0 && seq.getCharAt(index + 1) == 'M') {
                edu.ucsd.msjava.msutil.Peptide peptide = aaSet.getPeptide(match.getPepSeq());
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
            if (pre == '_') {
                pre = isNTermMetCleaved ? 'M' : '-';
            }
            char post;
            if (isNTermMetCleaved)
                post = seq.getCharAt(index + length);
            else
                post = seq.getCharAt(index + length - 1);
            if (post == '_') post = '-';

            int protStartIndex = (int) seq.getStartPosition(index);
            String annotation = seq.getAnnotation(protStartIndex);
            String accession = annotation.split("\\s+")[0];

            boolean isDecoy = accession.startsWith(decoyProteinPrefix);
            if (!isDecoy) isAllDecoy = false;

            String key = pre + accession + post;
            if (proteinSet.add(key)) {
                if (proteinBuf.length() != 0) proteinBuf.append(";");
                proteinBuf.append(accession).append("(pre=").append(pre).append(",post=").append(post).append(")");
            }
        }

        if (isAllDecoy) return "";
        return proteinBuf.toString();
    }
}
