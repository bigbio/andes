# MS-GF+

[MS-GF+ Documentation home](README.md) · [ChangeLog](Changelog.md "MS-GF+ ChangeLog")


```text
Usage: java -Xmx3500M -jar MSGFPlus.jar

[-conf ConfigurationFile] (Configuration file path; options specified at the command line will override settings in the config file)
   An example parameter file is at https://github.com/MSGFPlus/msgfplus/blob/master/docs/examples/MSGFPlus_Params.txt
   Additional parameter files are at https://github.com/MSGFPlus/msgfplus/tree/master/docs/ParameterFiles

[-s SpectrumFile] (*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt)
   Spectra should be centroided (see below for MSConvert example). Profile spectra will be ignored.

[-d DatabaseFile] (*.fasta or *.fa or *.faa)

[-decoy DecoyPrefix] (Prefix for decoy protein names; Default: XXX)

[-o OutputFile (*.mzid)] (Default: [SpectrumFileName].mzid)

[-t PrecursorMassTolerance] (e.g. 2.5Da, 20ppm or 0.5Da,2.5Da; Default: 20ppm)
   Use a comma to define asymmetric values. 
   E.g. "-t 0.5Da,2.5Da" will set 0.5Da to the left (ObservedPepMass < TheoreticalPepMass) 
                              and 2.5Da to the right (ObservedPepMass > TheoreticalPepMass)

[-ti IsotopeErrorRange] (Range of allowed isotope peak errors; Default: 0,1)
   Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation.
   The combination of -t and -ti determines the precursor mass tolerance.
   E.g. "-t 20ppm -ti -1,2" tests abs(ObservedPepMass - TheoreticalPepMass - n * 1.00335Da) < 20ppm for n = -1, 0, 1, 2.

[-thread NumThreads] (Number of concurrent threads to be executed; Default: Number of available cores)
   This is best set to the number of physical cores in a single NUMA node.
   Generally a single NUMA node is 1 physical processor.
   The default will try to use hyperthreading cores, which can increase the amount of time this process will take.
   This is because the part of Scoring param generation that is multithreaded is also I/O intensive.

[-tasks NumTasks] (Override the number of tasks to use on the threads; Default: internally calculated based on inputs)
   More tasks than threads will reduce the memory requirements of the search, but will be slower (how much depends on the inputs).
   1 <= tasks <= numThreads: will create one task per thread, which is the original behavior.
   tasks = 0: use default calculation - minimum of: (threads*3) and (numSpectra/250).
   tasks < 0: multiply number of threads by abs(tasks) to determine number of tasks (i.e., -2 means "2 * numThreads" tasks).
   One task per thread will use the most memory, but will usually finish the fastest.
   2-3 tasks per thread will use comparably less memory, but may cause the search to take 1.5 to 2 times as long.

[-verbose 0/1] (Console output message verbosity; Default: 0)
   0: Report total progress only
   1: Report total and per-thread progress/status

[-tda 0/1] (Target decoy strategy; Default: 0)
   0: Don't use a decoy database
   1: Search with a decoy database (forward + reverse proteins)

[-m FragmentationMethodID] (Fragmentation Method; Default: 0)
   0: As written in the spectrum or CID if no info
   1: CID
   2: ETD
   3: HCD
   4: UVPD

[-inst InstrumentID] (Instrument ID; Default: 0)
   0: Low-res LCQ/LTQ
   1: Orbitrap/FTICR/Lumos
   2: TOF
   3: Q-Exactive

[-e EnzymeID] (Enzyme ID; Default: 1)
   0: Unspecific cleavage
   1: Trypsin
   2: Chymotrypsin
   3: Lys-C
   4: Lys-N
   5: glutamyl endopeptidase
   6: Arg-C
   7: Asp-N
   8: alphaLP
   9: no cleavage

[-protocol ProtocolID] (Protocol ID; Default: 0)
   0: Automatic
   1: Phosphorylation
   2: iTRAQ
   3: iTRAQPhospho
   4: TMT
   5: Standard

[-ntt 0/1/2] (Number of Tolerable Termini; Default: 2)
   When EnzymeID is 1 (trypsin),
     2: Only search for fully-tryptic peptides
     1: Search for semi-tryptic and fully-tryptic peptides
     0: Non-tryptic search

[-mod ModificationFileName] (Modification file; Default: standard amino acids with fixed C+57; only if -mod is not specified)

[-minLength MinPepLength] (Minimum peptide length to consider; Default: 6)

[-maxLength MaxPepLength] (Maximum peptide length to consider; Default: 40)

[-minCharge MinCharge] (Minimum precursor charge to consider if charges are not specified in the spectrum file; Default: 2)

[-maxCharge MaxCharge] (Maximum precursor charge to consider if charges are not specified in the spectrum file; Default: 3)

[-n NumMatchesPerSpec] (Number of matches per spectrum to be reported; Default: 1)

[-addFeatures 0/1] (Include additional features in the output; enable this to post-process results with Percolator; Default: 0)
   0: Output basic scores only
   1: Output additional features

[-ccm ChargeCarrierMass] (Mass of charge carrier; Default: mass of a proton, 1.00727649)

[-ignoreMetCleavage 0/1] (N-terminal methionine cleavage behavior; Default: 0)

[-maxMissedCleavages Count] (Exclude peptides with more than this number of missed cleavages from the search; Default: -1, meaning no limit)

[-minNumPeaks Count] (Minimum number of ions a spectrum must have to be examined; Default: 10)

[-iso NumIsoforms] (Number of isoforms to consider per peptide; Default: 128)

[-numMods Count] (Maximum number of dynamic (variable) modifications per peptide; Default: 3)

[-allowDenseCentroidedPeaks 0/1] (Default: 0 (disabled); 1: (for mzML/mzXML input only) allows inclusion of spectra with high-density centroid data in the search)
   MS-GF+ checks the distance between consecutive peaks in the spectrum, and if the median distance is less than 50 ppm, they are considered profile spectra regardless of the value provided in mzML and mzXML files.
   This parameter allows overriding this check when the mzML/mzXML file says the spectrum is centroided.
      
```


#### Examples:

Example command (using a parameter file):

`java -Xmx3500M -jar MSGFPlus.jar -s Dataset.mzML -d ProteinList.fasta -conf MSGFPlus_PartTryp_MetOx_20ppmParTol.txt`

Example command (high-precision spectra, using arguments):

`java -Xmx3500M -jar MSGFPlus.jar -s Dataset.mzML -d IPI_human_3.79.fasta -inst 1 -t 20ppm -ti -1,2 -ntt 2 -tda 1 -o PSMs.mzid`

Example command (low-precision spectra):

`java -Xmx3500M -jar MSGFPlus.jar -s Dataset.mzML -d IPI_human_3.79.fasta -inst 0 -t 0.5Da,2.5Da -ntt 2 -tda 1 -o PSMs.mzid`

### Parameters:

- **-s SpectrumFile** (.mzML\*, \*.mzXML, \*.mgf, \*.ms2, \*.pkl or \*\_dta.txt) - Required

  - Spectrum file name. Currently, MS-GF+ supports the following file formats: mzML, mzXML, mzML, mgf, ms2, pkl and \_dta.txt.
  - We recommend to use mzML, whenever possible.
  - For Thermo .raw files, obtain a centroided .mzML using MSConvert, which is part of [ProteoWizard](http://proteowizard.sourceforge.net/).

  `MSConvert.exe --mzML --32 --filter "peakPicking true 1-" DatasetName.raw`

- **-d DatabaseFile** (\*.fasta or \*.fa or \*.faa) - Required

  - Path to the protein database file. If the database file does not have auxiliary index files (\*.canno, \*.cnlcp, \*.csarr, and \*.cseq), MS-GF+ will create them.
  - When "-tda 1" option is used, the database specified here must contain only target protein sequences.

  If multiple MS-GF+ processes access the same database file, it is strongly recommended to index the database prior to the database search by [running BuildSA](BuildSA.md).

- **-conf ConfigurationFile**
  - Path to the configuration file (aka parameter file) that defines settings for MS-GF+
  - Options specified at the command line will override settings in the config file
  - Example parameter file: [MSGFPlus_Params.txt](https://github.com/MSGFPlus/msgfplus/blob/master/docs/examples/MSGFPlus_Params.txt)
  - See also these additional [example parameter files](https://github.com/MSGFPlus/msgfplus/tree/master/docs/ParameterFiles)

- **-decoy DecoyPrefix**
  - Text to prepend to protein names when including decoy (reverse sequence) proteins in the .revCat.fasta file and related index files
  - Defaults to XXX (though an underscore is also added, giving names like `XXX_Contaminant_TRYP_BOVIN`)
  - Use `-decoy REV` to get names like `REV_Contaminant_TRYP_BOVIN`

<!-- -->

- **-o OutputFile** (\*.mzid)
  - Filename where the output (mzIdentML 1.1 format) will be written.
  - File extension must be "mzid" (case sensitive).
  - By default, the output file name will be "\[SpectrumFileName\].mzid".
  - E.g. for the input spectrum file "test.mzML", the output will be written to "test.mzid" if this parameter is not specified.

- **-t PrecursorMassTolerance** (Default: 20ppm)
  - Precursor mass tolerance in Da. or ppm. There must be no space between the number and the unit. E.g. `2.5Da` or `20ppm`
  - To set asymmetric tolerances, use a comma to separate left (observed mass \< theoretical mass) and  
    right (observed mass \> theoretical mass) tolerances.  
    E.g. `0.5Da,2.5Da`
  - It is recommended to use a tight tolerance rather than a loose tolerance (e.g. for Orbitrap data, 10ppm or 20ppm usually identifies more spectra than 50ppm).

- **-ti IsotopeErrorRange** (Default: 0,1)
  - Takes into account of the error introduced by choosing non-monoisotopic peak for fragmentation.
  - If the precursor mass tolerance is equal to or larger than 0.5Da or 500ppm, this parameter will be ignored.
  - The combination of -t and -ti determines the precursor mass tolerance.
  - E.g. `-t 20ppm -ti -1,2` tests abs(ObservedPepMass - TheoreticalPepMass - n \* 1.00335Da) \< 20ppm for n = -1, 0, 1, 2

- **-thread NumOfThreads** (Default: Number of available cores)
  - Number of concurrent threads to be executed together.
  - Default value is the number of available logical cores (e.g. 8 for quad-core processor with hyper-threading support).

- **-tasks NumTasks** (Default: internally calculated based on inputs)
  - Use this to manually set the number of tasks to create for the search.
  - More tasks than threads will reduce the memory requirements of the search, but will be slower (how much depends on the inputs).
  - If the spectrum file is particularly large, a larger number of tasks will decrease the possibility of out of memory errors.
  - If the fasta file being searched is larger than 10MB, more tasks will cause a noticeably longer search time.
  - 1 \<= tasks \<= numThreads: will create one task per thread, which is the original behavior.
  - tasks = 0: use default calculation - minimum of :(threads\*3), and (numSpectra/250).
  - tasks \< 0: multiply number of threads by abs(tasks) to determine number of tasks (i.e., -2 =\> "2 \* numThreads" tasks).
  - One task per thread will use the most memory, but will usually finish the fastest.
  - 2-3 tasks per thread will use comparably less memory, but may cause the search to take 1.5 to 2 times as long with a 23MB fasta file.

- **-verbose 0/1** (Default: 0)
  - Changes the verbosity of the output
  - If 0, only the overall progress is reported, creating the minimal useful output to console.
  - If 1, you see all of the output of 0, but with additional console output from each thread and task.
  - 1 will produce console output that matches the console output of older versions.

- **-tda 0/1** (Default: 0)

  - Indicates whether to search normal (forward only) protein sequences, or a decoy file where the reversed protein sequences are appended to the normal protein sequences
    - 0: Search the protein sequences as listed in the FASTA file (the target database)
    - 1: Search a target-decoy database, allowing for the computation of QValues (FDR)
  - QValue is defined as the minimum false discovery rate (FDR) at which the test may be called significant (ReversePeptideCount / ForwardPeptideCount)
    - QValue(t) = (Number of DecoyPSMs with score equal or above t) ÷ (Number of TargetPSMs with score equal or above t)
    - PSM: Peptide-Spectrum Match
    - -log(SpecProb) is used as the score to compute QValue.

  If `-tda 1` is specified, MS-GF+ automatically creates a combined target/reversed database file (DBFileName.revConcat.fasta).  
  Thus, when specifying "-d" parameter, DatabaseFile must contain only target proteins.

<!-- -->

- **-m FragmentationMethodID** (Default: 0)
  - Fragmentation method identifier (used to determine the scoring model).
    - 0: As written in the spectrum or CID if no info (default)
    - 1: CID
    - 2: ETD
    - 3: HCD
    - 4: UVPD
  - If the identifier is 0 and fragmentation method is written in the spectrum file (e.g. mzML files), MS-GF+ will recognize the fragmentation method and use a relevant scoring model.
  - If the identifier is 0 and there is no fragmentation method information in the spectrum (e.g. mgf files), CID model will be used by default.
  - If the identifier is non-zero and the spectrum has fragmentation method information, only the spectra that match with the identifier will be processed.
  - If the identifier is non-zero and the spectrum has no fragmentation method information, MS-GF+ will process all spectra assuming the specified fragmentation method.
- **-inst InstrumentID**
  - Identifier of the instrument used to generate MS/MS spectra (this parameter defines the the scoring model).
    - 0: Low-res LCQ/LTQ (Default for CID and ETD)
    - 1: Orbitrap/FTICR/Lumos (Default for HCD)
    - 2: TOF
    - 3: Q-Exactive
  - For "hybrid" spectra with high-precision MS1 and low-precision MS2, use 0.
  - For usual low-precision instruments (e.g. Thermo LTQ), use 0.
  - If MS/MS fragment ion peaks are of high-precision (e.g. tolerance = 10ppm), use 2.
  - For TOF instruments, use 2.
  - For Q-Exactive HCD spectra, use 3.
  - For other HCD spectra, use 1.
- **-e EnzymeID** (Default: 1)
  - Enzyme identifier.
    - 0: unspecific cleavage (cleave after any residue)
    - 1: Trypsin (default)
    - 2: Chymotrypsin
    - 3: Lys-C
    - 4: Lys-N
    - 5: glutamyl endopeptidase (Glu-C)
    - 6: Arg-C
    - 7: Asp-N
    - 8: alphaLP
    - 9: no cleavage
  - Use 9 for peptidomics studies
  - Create file params\enzymes.txt (or params/enzymes.txt on Linux) below the working directory to define custom enzymes or override the cleavage residues for built-in enzymes
  - For more info, see [enzymes.txt](examples/enzymes.txt)
- **-protocol ProtocolID** (Default: 0)
  - Protocol identifier. Protocols are used to enable scoring parameters for enriched and/or labeled samples.
    - 0: Automatic (Default)
      - This will set the protocol based on the names of the modifications in Mods.txt
      - It looks for names (case insensitive) that start with "itraq", "phospho", and "tmt"
    - 1: Phosphorylation: for phosphopeptide enriched samples
    - 2: iTRAQ: for iTRAQ-labeled samples
    - 3: iTRAQPhospho: for phosphopeptide enriched and iTRAQ-labeled samples
    - 4: TMT: for TMT-labeled samples
    - 5: Standard: for samples not in the above categories (no protocol)
- **-ntt 0/1/2** (Default: 2)
  - Number of tolerable termini (aka tryptic termini)
  - This parameter is used to apply the enzyme cleavage specificity rule when searching the database.
  - Specifies the minimum number of termini matching the enzyme specificity rule.
    - For example, for trypsin, K.ACDEFGHR.C (NTT=2), G.ACDEFGHR.C (NTT=1), K.ACDEFGHI.C (NTT=1) and G.ACDEFGHR.C (NTT=0).
    - `-ntt 2` will search for fully tryptic peptides only.
  - By default, `-ntt 2` is used.
  - Using `-ntt 1` or `-ntt 0` can make the search significantly slower.

<!-- -->

- **-mod ModificationFile** (Default: standard amino acids with fixed C+57, though only if `-mod` is not specified)
  - Modification file name. ModificationFile contains the modifications to be considered in the search.
  - If `-mod` is not specified, standard amino acids with fixed Carbamidomethylation C will be used.
  - See an [example MS-GF+ modification file](examples/Mods.txt).
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
- **-addFeatures 0/1** (Default: 0)
  - If 0, only basic scores are reported.
  - If 1, the following features are reported
    - MS2IonCurrent: Summed intensity of all product ions
    - ExplainedIonCurrentRatio: Summed intensity of all matched product ions (e.g. b, b-H2O, y, etc.) divided by MS2IonCurrent
    - NTermIonCurrentRatio: Summed intensity of all matched prefix ions (e.g. b, b-H2O, etc.) divided by MS2IonCurrent
    - CTermIonCurrentRatio: Summed intensity of all matched suffix ions (e.g. y, y-H2O, etc.) divided by MS2IonCurrent
- **-ccm ChargeCarrierMass** (Default: 1.00727649)
  - Override the default charge carrier mass
- **-ignoreMetCleavage 0/1** (Default: 0)
  - 0: consider cleavage of methionine from the protein's N-terminus, even when NTT=2
  - 1: disable N-terminal methionine cleavage
- **-maxMissedCleavages Count** (Default: -1, meaning no limit)
  - Exclude peptides with more than this number of missed cleavages
- **-numMods Count** (Default: 3)
  - Maximum number of dynamic (variable) modifications per peptide
  - If this value is large and multiple dynamic modifications are defined, the search will be slow (depending on FASTA file size)

### MS-GF+ output

MS-GF+ outputs results as an mzIdentML (version 1.1) file. See <http://www.psidev.info/mzidentml/> for details on the mzIdentML format. For every PSM, MS-GF+ reports the following scores:

- **MS-GF:RawScore**: MS-GF+ raw score of the peptide-spectrum match
- **MS-GF:DeNovoScore:** the score of the optimal scoring peptide for the spectrum (not necessary in the database) (MS-GF:RawScore \<= MS-GF:DeNovoScore)
- **MS-GF:SpecEValue**: spectral E-value (spectrum level E-value) of the peptide-spectrum match - the lower the better
- **MS-GF:EValue**: database level E-value (expected number of peptides in a random database having equal or better scores than the PSM score) - the lower the better
- **MS-GF:QValue**
  - PSM-level Q-value estimated using the target-decoy approach.
  - MS-GF:QValue is computed solely based on MS-GF:SpecEValue.
- **MS-GF:PepQValue**
  - Peptide-level Q-value estimated using the target-decoy approach.
  - Reported only if "-tda 1" is specified.
  - If multiple spectra are matched to the same peptide, only the best scoring PSM (lowest SpecProb) is retained.  
    After that, MS-GF:PepQValue is calculated as \#DecoyPSMs\>s / \#TargetPSMs\>s among the retained PSMs.  
    This approximates the Q-value of the set of unique peptides.
  - In the MS-GF+ output, the same PepQValue value is given to all PSMs sharing the peptide.
    - Thus, even a low-quality PSM may get a low PepQValue (if it has a high-quality "sibling" PSM sharing the peptide).
    - Note that this should not be used to count the number of identified PSMs.

### MS-GF+ output example

Shown below is a sample of the MS-GF+ output in table form, as extracted from a simple MzIdentML file: [test.mzid](examples/test.mzid)

There are two options for converting an MS-GF+ output file (**.mzid**) into a tab-separated file (**.tsv**).

1.  The MzIDToTsv utility built into MSGFPlus.jar (see the [MzIDToTsv](MzidToTsv.md) page)
    - Easy to access (though syntax is a bit tricky)
    - Can be slow for large .mzid files
2.  The Mzid-To-Tsv-Converter standalone application, [available on GitHub](https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases)
    - Fast conversion
    - Handles large .mzid files
    - Runs natively on Windows, but requires mono to use on Linux


| \#SpecFile | SpecID | ScanNum | FragMethod | Precursor | IsotopeError | PrecursorError(ppm) | Charge | Peptide | Protein | DeNovoScore | MSGFScore | SpecEValue | EValue | QValue | PepQValue |
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| test.mgf | index=0 | 26559 | CID | 1285.3457 | 1 | -5.049801 | 3 | K.IGAYLFVDMAHVAGLIAAGVYPNPVPHAHVVTSTTHK.T | test | 299 | 244 | 1.4807088E-31 | 3.2871733E-29 | 0.0 | 0.0 |
| test.mgf | index=0 | 26559 | CID | 1285.3457 | 1 | -5.049801 | 3 | K.IGAYLFVDMAHVAGLIAAGVYPNPVPHAHVVTSTTHK.T | test_isoform | 299 | 244 | 1.4807088E-31 | 3.2871733E-29 | 0.0 | 0.0 |
| test.mgf | index=1 | -1 | CID | 870.11743 | 0 | 0.14029178 | 3 | K.NLANPTSVILASIQM+15.995LEYLGMADK.A | test2 | 156 | 136 | 2.2559852E-22 | 4.4217308E-20 | 0.0 | 0.0 |


(Text file of this table: [test_Unrolled.tsv](examples/test_Unrolled.tsv))
