# MS-GFDB

[MS-GF+ Documentation home](README.md)

MS-GFDB is an old application that is no longer under development. It was supserseded by [MS-GF+](MSGFPlus.md).  
MS-GF+ has all the functionalities provided by MS-GFDB, plus numerous improvements.

### Differences between MS-GF+ and MS-GFDB

- **Input**
  - MS-GF+ supports mzML in addition to mzXML, mgf, ms2, pkl and \_dta.txt
  - "-t PrecursorMassTolerance" is optional with MS-GF+ (default 20ppm)
  - "-c13 0/1/2" was changed to "-ti IsotopeErrorRange" in MS-GF+
    - IsotopeErrorRange: MinIsotopeError,MaxIsotopeError (both are inclusive)
    - -c13 x == -ti 0,x
  - "-nnet" was changed to "-ntt" in MS-GF+
    - -nnet 0 == -ntt 2, -nnet 1 == -ntt 1, -nnet 2 == -ntt 0
  - Modification file format change
    - In MS-GF+, the name of the modification should match the PSI-MS name (accessible from [http://www.unimod.org](http://www.unimod.org/))
    - CompositionStr can take Br, Cl, Fe, Se in addition to C, H, N, O, S, and P
    - The sequence of the atoms can be arbitrary.
      - Previously C2H2O was valid but OH2C2 was invalid
      - With MS-GF+, both are valid
  - "-uniformAAProb 0/1" was deleted in MS-GF+
  - "-addFeatures 0/1" was added to MS-GF+; "-addFeatures 1" will output the following extra features for each PSM (will be useful to downstream tools like Percolator or IDPicker):
    - MS2IonCurrent: Summed intensity of all product ions
    - ExplainedIonCurrentRatio: Summed intensity of all matched product ions (e.g. b, b-H2O, y, etc.) divided by MS2IonCurrent
    - NTermIonCurrentRatio: Summed intensity of all matched prefix ions (e.g. b, b-H2O, etc.) divided by MS2IonCurrent
    - CTermIonCurrentRatio: Summed intensity of all matched suffix ions (e.g. y, y-H2O, etc.) divided by MS2IonCurrent
  - "-showQValue 0/1" was added to MS-GF+
  - "-showDecoy 0/1" was added to MS-GF+
- **Output**
  - Output format for MS-GF+ is the HUPO PSI mzIdentML version 1.1 (\*.mzid); see <http://www.psidev.info/mzidentml> for details.
    - Decoy protein prefix is "XXX" in MS-GF+ (vs. "REV" in MS-GFDB)
  - MS-GF+ provides a converter from mzIdentML to tsv (the resulting tsv file will be similar to the MS-GFDB output file).
    - The converter is included in the MSGFPlus.jar file
    - It can be run by "java -Xmx2000M edu.ucsd.msjava.ui.MzIDToTsv"
    - A faster converter that supports larger result files is the Mzid-To-Tsv-Converter, [available on GitHub](https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases). This is a C# application that works under Windows or on Linux with mono.
  - Difference between the MS-GFDB output and the MS-GF+ TSV output
    - MS-GF+ includes SpecID (native spectrum ID) instead of SpecIndex
    - MS-GF+ reports IsotopeError
    - When a peptide matches to multiple proteins, all protein accessions will be reported by MS-GF+
    - SpecProb was renamed to SpecEValue in MS-GF+
    - MS-GF+ reports EValue (database-level E-value) instead of PValue (database-level P-value)
    - FDR and PepFDR were renamed to QValue and PepQValue in MS-GF+

# MS-GFDB


```text
Usage: java -Xmx2000M -jar MSGFDB.jar
    -s SpectrumFile (*.mzXML, *.mzML, *.mgf, *.ms2, *.pkl or *_dta.txt)
    -d DatabaseFile (*.fasta or .fa)
    -t ParentMassTolerance (e.g. 2.5Da, 30ppm, or 0.5Da,2.5Da)
       Use comma to set asymmetric values. E.g. "-t 0.5Da,2.5Da" will set 0.5Da to the left (expMass<theoMass) and 2.5Da to the right (expMass>theoMass).
    [-o outputFileName] (Default: stdout)
    [-thread NumOfThreads] (Number of concurrent threads to be executed, Default: Number of available cores)
    [-tda 0/1] (0: don't search decoy database (default), 1: search decoy database to compute FDR)
    [-m FragmentationMethodID] (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD, 4: Merge spectra from the same precursor)
    [-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: High-res LTQ (Default for HCD), 2: TOF)
    [-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N, 8: aLP, 9: Endogenous peptides)
    [-c13 0/1/2] (Number of allowed C13, Default: 1)
    [-nnet 0/1/2] (Number of allowed non-enzymatic termini, Default: 1)
    [-mod ModificationFileName] (Modification file, Default: standard amino acids with fixed C+57)
    [-minLength MinPepLength] (Minimum peptide length to consider, Default: 6)
    [-maxLength MaxPepLength] (Maximum peptide length to consider, Default: 40)
    [-minCharge MinPrecursorCharge] (Minimum precursor charge to consider if not specified in the spectrum file, Default: 2)
    [-maxCharge MaxPrecursorCharge] (Maximum precursor charge to consider if not specified in the spectrum file, Default: 3)
    [-n NumMatchesPerSpec] (Number of matches per spectrum to be reported, Default: 1)
    [-uniformAAProb 0/1] (0: use amino acid probabilities computed from the input database (default), 1: use probability 0.05 for all amino acids)
```


### Parameters:

- **-s SpectrumFile** (\*.mzXML, \*.mzML, \*.mgf, \*.ms2, \*.pkl or \*\_dta.txt) - Required
  - Spectrum file name. Currently, MS-GFDB supports the following file formats: mzXML, mzML, mgf, ms2, pkl and \_dta.txt.
- **-d DatabaseFile** (\*.fasta or \*.fa) - Required
  - Path to the protein database file. If the database file does not have auxiliary index files (\*.canno, \*.cnlcp, \*.csarr, and \*.cseq), MS-GFDB will create them.
  - When "-tda 1" option is used, the database must contain only target protein sequences.

If multiple MS-GFDB processes access the same database file, it is strongly recommended to index the database prior to the database search by running BuildSA (see below).

- **-t ParentMassTolerance** - Required
  - Parent mass tolerance in Da. or ppm. There must be no space between the number and the unit. E.g. 2.5Da, 30ppm
  - To set asymmetric tolerances, use comma to separate left (experimental mass \< theoretical mass) or right (experimental mass \> theoretical mass) tolerances. E.g. 0.5Da,2.5Da
- **-o OutputFile** (Default: stdout)
  - Filename where the output will be written.
  - The output will be printed to standard out by default.
- **-thread NumOfThreads** (Number of concurrent threads to be executed, Default: Number of available cores)
  - Number of concurrent threads to be executed together.
  - Default value is the number of available logical cores (e.g. 8 for quad-core processor with hyper-threading support).
- **-tda 0/1** (0: don't search decoy database (default), 1: search decoy database to compute FDR)
  - Indicates whether to search the decoy database or not.
  - If 0, the decoy database is not searched and FDRs are theoretically derived from P-values (EFDR).
  - If 1, FDRs are computed based on the target-decoy approach (i.e. reversed database is appended to the target database and MS-GFDB searches the combined database)
    - FDR(t) = \#(DecoyPSMs with score equal or above t) / \#(TargetPSMs with score equal or above t).
    - PSM: Peptide-Spectrum Match
    - -log(SpecProb) is used as the score to compute FDR.

If -tda 1 is specified, MS-GFDB automatically creates a combined target/reversed database file (DBFileName.revConcat.fasta). Thus, when specifying "-d" parameter, DatabaseFile must contain only target proteins.

- **-m FragmentationMethodID** (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD, 4: Merge spectra from the same precursor)
  - Fragmentation method identifier (used to determine the scoring model).
  - If the identifier is 0 and fragmentation method is written in the spectrum file (e.g. activationMethod field in mzXML files), MS-GFDB will recognize the fragmentation method and use a relevant scoring model.
  - If the identifier is 0 and there is no fragmentation method information in the spectrum (e.g. mgf files), CID model will be used by default.
  - If the identifier is non-zero and the spectrum has fragmentation method information, only the spectra that match with the identifier will be processed.
  - If the identifier is non-zero and the spectrum has no fragmentation method information, MS-GFDB will process all spectra assuming the specified fragmentation method.
  - If the identifier is 4, MS/MS spectra from the same precursor ion (e.g. CID/ETD pairs, CID/HCD/ETD triplets) will be merged and the "merged" spectrum will be used for searching instead of individual spectra. See Kim et al., MCP 2010 for details.
- **-inst InstrumentID** (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: TOF , 2: High-res LTQ (Default for HCD))
  - Identifier of the instrument to generate MS/MS spectra (used to determine the scoring model).
  - For "hybrid" spectra with high-precision MS1 and low-precision MS2, use 0.
  - For usual low-precision instruments (e.g. Thermo LTQ), use 0.
  - For TOF instruments, use 1.
  - If MS/MS fragment ion peaks are of high-precision (e.g. tolerance = 10ppm), use 2.
- **-e EnzymeID** (Default: 1)
  - Enzyme identifier. Trypsin (1) will be used by default.
  - 0: No enzyme, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: Endogenous peptides
- **-c13 0/1/2** (Number of allowed isotope errors, Default: 1)
  - Instruments often choose 2nd or 3rd isotope peak instead of mono-isotope peak from MS1 spectrum.
  - If this value is non-zero, expPeptideMass-1.00335 (i.e. mass(13C)-mass(12C)) and expPeptideMass-2.00671 (i.e. 2\*(mass(C13)-mass(C12)) (only if -c13 2) will be considered along with expPeptideMass.
  - If accurate precursor ion mass is available (e.g. LTQ-Orbitrap), it is better to set a narrow parent mass tolerance and non-zero -c13 value (e.g. -t 30ppm -c13 1) than to set a wide tolerance (e.g. -t 0.5Da,2.5Da).
  - If the parent mass tolerance is equal to or larger than 0.5Da or 500ppm, this parameter will be ignored.
- **-nnet 0/1/2** (Number of allowed non-enzymatic termini, Default: 1)
  - This parameter is used to determine the enzyme cleavage rule.
  - Specifies the maximum number of peptide termini that are not cleaved by the enzyme.
  - For example, for trypsin, K.ACDEFGHR.C, G.ACDEFGHR.C, K.ACDEFGHI.C and G.ACDEFGHR.C have 0, 1, 1 and 2 non-enzymatic termini, accordingly.
  - By default, -nnet 1 will be used. Using -nnet 0 (or 2) will make the search significantly faster (slower).
- **-mod ModificationFile** (Default: standard amino acids with fixed C+57)\]
  - Modification file name. ModificationFile contains the modifications to be considered in the search.
  - If -mod option is not specified, standard amino acids with fixed Carbamidomethylation C will be used.
  - See an [example MS-GFDB modification file](MSGFDB_ModFile.md).
- **-minLength MinPepLength** (Default: 6)
  - Minimum length of the peptide to be considered.
- **-maxLength MaxPepLength** (Default: 40)
  - Maximum length of the peptide to be considered.
- **-minCharge MinPrecursorCharge** (Default: 2)
  - Minimum precursor charge to consider. This parameter is used only for spectra with no charge.
- **-maxCharge MinPrecursorCharge** (Default: 3)
  - Maximum precursor charge to consider. This parameter is used only for spectra with no charge.
- **-n NumMatchesPerSpec** (Default: 1)
  - Number of peptide matches per spectrum to report.
  - Expected false discovery rates (EFDRs) will be reported only when this value is 1.
- **-uniformAAProb** 0/1 (Default: 0)
  - If 0, compute amino acid frequencies from the input database and use them as amino acid probabilities.
  - If 1, use uniform amino acid probability (preferable when the database size is small).

### MS-GFDB output

MS-GFDB outputs a tab-delimited file with the following columns: \#SpecFile, Scan#, FragMethod, Precursor, PMError, Charge, Peptide, Protein, DeNovoScore, MSGFScore, SpecProb, P-value, EFDR.

- **SpecFile**: spectrum file name
- **SpecIndex**: spectrum index (1-based) in the file. The first spectrum has index 1, the second has index 2, and so on. For mzXML files this value is same as the scan number.
- **Scan#**: scan number of the spectrum. If the scan number is not available, the value will be -1.
- **FragMethod**: fragmentation method used to generate the spectrum (e.g. CID, ETD, etc.). When spectra from the same precursor are merged, fragmentation methods of merged spectra will be shown as a form "FragMethod1/FragMethod2/..." (e.g. CID/ETD, CID/HCD/ETD).
- **Precursor**: precursor mass in m/z or ppm
- **Charge**: precursor ion charge
- **Peptide**: peptide sequence with neighboring amino acids
- **Protein**: protein name
- **DeNovoScore**: the score of the optimal scoring peptide (not necessary in the database)
- **MSGFScore**: MS-GF raw score of the peptide-spectrum match (MSGFScore \<= DeNovoScore)
- **SpecProb**: spectral probability (spectrum level p-value) of the peptide-spectrum match
- **P-value**: database level p-value (probability that a random PSM have an equal or better score against a random database of the same size)
- **EFDR** or **FDR**: false discovery rate
  - If "-tda 1" is specified, FDRs are estimated using the target-decoy approach using the spectral probability (SpecProb) as the score (the lower, the better).
  - Otherwise, FDRs are estimated using P-values without searching the decoy database (EFDR). See Gupta et al., JASMS 2011 for details.
  - MS-GFDB reports EFDR only when it is configured to report 1 peptide match per spectrum (i.e. -n 1).
  - EFDR accurately estimates FDR when the parent mass tolerance is equal or larger than 0.5.
  - EFDR conservatively estimates FDR when the parent mass tolerance is small.
    - E.g. When parent mass tolerance is 30ppm, at EFDR 1% threshold, one identifies approximately 7% less peptide-spectrum matches (PSMs) compared to the case when the target-decoy approach is used to estimate the FDR.
- **PepFDR**
  - Peptide-level FDR estimated using the target-decoy approach.
  - Reported only if "-tda 1" is specified.
  - If multiple spectra are matched to the same peptide, only the best scoring PSM (lowest SpecProb) is retained. After that, PepFDR is calculated as \#DecoyPSMs\>s / \#TargetPSMs\>s among the retained PSMs. This approximates the FDR of the set of unique peptides. In the MS-GFDB output, the same PepFDR value is given to all PSMs sharing the peptide. So, even a low-quality PSM may get a low PepFDR value (if it has a high-quality "sibling" PSM sharing the peptide). Note that this should not be used to count the number of identified PSMs.

### MS-GFDB output example


| \#SpecFile | SpecIndex | Scan# | FragMethod | Precursor | PMError(ppm) | Charge | Peptide | Protein | DeNovoScore | MSGFScore | SpecProb | P-value | FDR | PepFDR |
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| 090121_NM_Trypsin_20.mzXML | 2838 | 2838 | CID | 964.7707 | 1.5199227 | 3 | K.TIQNSSVSPTSSSSSSSSTGETQTQSSSR.L | IPI:IPI00002349.2\|SWISS-PROT:Q7Z417\|TREMBL:A1L3A7\|ENSEMBL:ENSP00000225388\|REFSEQ:NP_065823\|H-INV:HIT000001036\|VEGA:OTTHUMP00000181037 Tax_Id=9606 Gene_Symbol=NUFIP2 Nuclear fragile X mental retardation-interacting protein 2 | 190 | 181 | 9.380133E-30 | 2.9333857E-22 | 0.0 | 0.0 |
| 090121_NM_Trypsin_20.mzXML | 3671 | 3671 | ETD | 1113.4758 | 0.6583758 | 2 | R.VGPADDGPAPSGEEEGEGGGEAGGK.E | IPI:IPI00016725.2\|SWISS-PROT:Q9UKN8\|TREMBL:B3KNH2;Q05CN7\|ENSEMBL:ENSP00000361219\|REFSEQ:NP_036336\|H-INV:HIT000071196\|VEGA:OTTHUMP00000022434 Tax_Id=9606 Gene_Symbol=GTF3C4 General transcription factor 3C polypeptide 4 | 162 | 158 | 1.9912463E-28 | 6.0892146E-21 | 0.0 | 0.0 |
| 090121_NM_Trypsin_20.mzXML | 3031 | 3031 | ETD | 651.64874 | 1.7510794 | 3 | K.GAAAAAAASGAAGGGGGGAGAGAPGGGR.L | IPI:IPI00644073.1\|VEGA:OTTHUMP00000038687 Tax_Id=9606 Gene_Symbol=INTS3 18 kDa protein | 214 | 202 | 6.7318633E-28 | 2.093763E-20 | 0.0 | 0.0 |
| 090121_NM_Trypsin_20.mzXML | 19088 | 19088 | CID | 1199.0916 | 10.392676 | 2 | K.VNFSPPGDTNSLFPGTWYLER.V | IPI:IPI00945760.1\|TREMBL:B7Z784;B7Z7M8;B7Z8R3\|REFSEQ:NP_001159579 Tax_Id=9606 Gene_Symbol=HMGCS2 hydroxymethylglutaryl-CoA synthase, mitochondrial isoform 2 precursor | 243 | 243 | 2.9611275E-27 | 8.838129E-20 | 0.0 | 0.0 |
| 090121_NM_Trypsin_20.mzXML | 3030 | 3030 | CID/ETD | 651.64874 | 1.7510794 | 3 | K.GAAAAAAASGAAGGGGGGAGAGAPGGGR.L | IPI:IPI00644073.1\|VEGA:OTTHUMP00000038687 Tax_Id=9606 Gene_Symbol=INTS3 18 kDa protein | 389 | 389 | 7.508096E-33 | 2.335189E-25 | 0.0 | 0.0 |


# BuildSA

Index a protein database for fast searching.


```text
Usage: java -cp MSGFDB.jar msdbsearch.BuildSA
    -d DatabaseFile (*.fasta or *.fa)
    [-tda 0/1/2] (0: target only, 1: target-decoy database only, 2: both)
```


**Parameters:**

- **-d DbPath**
  - Name of a protein database (\*.fasta or \*.fa)
  - Database file must ends with ".fasta" or ".fa".
- **-tda 0/1/2**
  - If 0, only "DatabaseFile" will be indexed.
  - If 1, a new database file (\*.revConcat.fasta) will be generated by appending reversed proteins. This forward-reverse database will be indexed.
  - If 2, both the original database and the forward-reverse database file will be indexed.

BuildSA creates a suffix array of the protein database. For an input database file DBFileName.fasta, BuildSA will generate 4 auxiliary files (DbFileName.canno, DBFileName.cnlcp, DBFileName.csarr, DBFileName.cseq).It needs to be executed only once per each database file.
