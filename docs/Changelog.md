# MS-GF+ ChangeLog

[MS-GF+ Documentation home](README.md)

**vNEXT — Unreleased (breaking change)**

- **BREAKING:** Removed mzIdentML (`.mzid`) output format. MS-GF+ now feeds directly into Percolator via the `.pin` format, which is the expected downstream pipeline. Users who previously relied on `.mzid` output should:
  - Switch to `-outputFormat pin` (the new default) and process the `.pin` with Percolator.
  - Or use `-outputFormat tsv` for a simple tabular summary.
  - Or use an older release (v2026.03.25 or earlier) if legacy `.mzid` output is strictly required.
- `-outputFormat` enum values changed: `0=pin` (new default), `1=tsv`. The previous `0=mzid`, `2=both`, `3=pin` layout is retired. Integer `2` and `3` are now rejected.
- Added precursor mass calibration: `-precursorCal auto|on|off` (default `auto`). Merged via PR #22.
- Added experimental fragment-index-based Tier-1 search: `-useFragmentIndex off|on|compare` (default `off`). Off path is bit-identical to the SA-walk baseline; on path is experimental and not yet recommended for production use.
- Added Percolator `.pin` output with OpenMS-parity features (`enzN`, `enzC`, `enzInt`, `mass`, `lnDeltaSpecEValue`, `matchedIonRatio`) and lowercase renames (`peplen`, `charge2/3/4`, `dm`, `absdm`, `isotope_error`) for OpenMS `PercolatorAdapter` interoperability. Merged via PR #22.

**v2026.03.25**

- Add EnzymeID 10: TrypPlusC (cleave after K, R, or C)

**v2023.01.12**

- Add parameter and output messages for working with particularly dense centroided data (read from mzML or mzXML)

**v2022.04.18**

- Check for negative masses when generating the list of candidate peptides

**v2022.01.07**

- Add support for a time range for RTINSECONDS in .mgf files

**v2021.09.06**

- Show warnings when scores are out of the expected range, but continue searching

**v2021.03.22**

- When displaying parameters, show the value for IgnoreMetCleavage
- Update online documentation, including example parameter files

**v2021.01.15**

- Allow integer, float, and double values in a parameter file to end with an exclamation mark (which is removed)

**v2021.01.07**

- Show the spectrum file name and format when opening the input file
- Update .mgf parsing logic to support negative charge states
- Warn of negative polarity spectra in .mzML and .mgf files

**v2020.08.05**

- Add .faa as a supported extension for FASTA files

**v2020.07.02**

- Update publication links and error message

**v2020.06.22**

- Show progress messages every 2 seconds when indexing the FASTA file

**v2020.06.16**

- Fix bug parsing a custom charge carrier mass from the command line (-ccm) or from a parameter file (e.g. ChargeCarrierMass=1.00727649)

**v2020.03.14**

- Update the embedded Unimod 'brick' definitions to fix the mzid generation failure after updating Unimod.obo

**v2020.03.12**

- Update the embedded Unimod.obo file to support recent Unimod accessions

**v2020.01.15**

- Resolve a new issue with duplicate peptide ids in the mzid output (problem with 2 modifications rounding to the same integer mass)
- Don't exit with exit code '0' when there is a parameter parsing error
- Add additional example files to the documentation

**v2019.08.26**

- Improve support for searching mgf files, particularly those containing spectra with no charge in the title

**v2019.07.03**

- `-s` and `-d` are now optional at the command line, allowing the spectrum file and database file to be defined in the configuration file

**v2019.06.28**

- Auto-create the output directory if missing
- Verify that the user has write access to the output directory

**v2019.06.20**

- Allow the title line of spectra in a .mgf file to have DatasetName.ScanStart.ScanEnd.Charge followed by a space and extra info
- Output the Java version and OS information to the console
- Update jmzReader and jmzidentml dependencies

**v2019.04.18**

- Fix to issue \#62, which prevents dropping of the highest m/z data point from most MS/MS spectra

**v2019.02.28**

- Fix a math error that introduced new erroneous output.

**v2019.02.27**

- Fix a bug that would output a single peptide evidence twice, under certain circumstances
- Prevent inaccurate information in peptide evidences due to differences in how protein N-Term Methionine cleavage is handled and tracked
  - This prevents an issue where a peptide that did have N-Term Methionine cleavage would say it started at residue 1, instead of residue 2

**v2019.02.20**

When reading settings from a configuration file, store custom amino acid definitions first

- Allows for dynamic modifications to be applied to custom amino acids

**v2019.02.14**

Add support for reading settings from a configuration file

- Use the `-conf` switch
- Example parameter file: [MSGFPlus_Params.txt](examples/MSGFPlus_Params.txt)

**v2019.02.07**

Improved support for customizing enzyme definitions

- Create file params\enzymes.txt (or params/enzymes.txt on Linux) below the working directory to define custom enzymes or override the cleavage residues for built-in enzymes
- For more info, see [enzymes.txt](examples/enzymes.txt)

**v2019.02.05**

Add option to use a custom prefix for decoy proteins (default is XXX\_)

- Both MSGFPlus.jar and BuildSA now support the `-decoy` switch

**v2019.02.01**

- Add validation of user-defined dynamic and static modifiations to check for duplicates (same name, different mass)
- Add validation to check for the same dynamic or static modification being defined for the same residue (or N- or C- terminus)
- Validate user-defined modification masses against default values for commonly used modifications
  - Warn if the mass is more than 0.01 Da larger or smaller than the default
  - Common modifications tracked via array [defaultModList](https://github.com/MSGFPlus/msgfplus/blob/8cba976300a651bab9104b73bc79030912b51185/src/main/java/edu/ucsd/msjava/msutil/Modification.java#L146) in file [Modification.java](https://github.com/MSGFPlus/msgfplus/blob/master/src/main/java/edu/ucsd/msjava/msutil/Modification.java)

**v2019.01.22**

- Fix bug reading spectra from .mgf files that start with a byte order mark (BOM)

**v2019.01.18**

- Skip spectra with a 0 m/z precursor ion

**v2019.01.11**

- Add support for reading Unicode files that start with a byte order mark (BOM)

**v2018.09.12**

- In BuildSA, validate that the output directory exists and is writable
- Update BuildSA syntax to include -o and mention that -d can be a directory with several FASTA files (addresses [Issue \#47](https://github.com/MSGFPlus/msgfplus/issues/47))

**v2018.07.17**

- Change how rank is set in mzid results to properly comply with the mzIdentML specification (noticed thanks to adder in [Issue \#40](https://github.com/MSGFPlus/msgfplus/issues/40))
- Get the scan retention time (and other cvParams) from the mzML ScanList for non-Thermo instruments (removed a check for a userParam that would prevent accessing the cvParams)
- Add support for newer nativeID formats (in part thanks to [Matt Chambers](https://github.com/chambm)' [pull request](https://github.com/MSGFPlus/msgfplus/pull/41))
- Read some additional CV Params from mzML and output them to mzIdentML (fixes [Issue \#32](https://github.com/MSGFPlus/msgfplus/issues/32))

**v2018.06.28**

- Add option to specify the maximum number of missed cleavages (thanks to Sean at [CWRU-CPB](http://proteomics.case.edu/))
- Make command line argument names case-insensitive
- Add additional checks to prevent the output of duplicate peptide evidences (fixes [Issue \#24](https://github.com/MSGFPlus/msgfplus/issues/24))

**v2018.04.09**

- When the EnzymeID is 9 (NoCleavage) do not cleavage after any residue (thanks to Sean at [CWRU-CPB](http://proteomics.case.edu/))
  - Previously EnzymeID 0 and EnzymeID 9 were identical.
  - Now EnzymeID 0 cleaves after every residue while EnzymeID 9 won't cleave after any residue (useful for peptidomics)
- Allow the fasta file to have extension .faa (in addition to .fasta and .fa)

**v2018.01.30**

- MzIdentML creation: Don't output an empty ModificationParams element when there are no modifications (because not including it complies with the XML Schema, while including it with no child nodes does not)
- Java 9 compatibility: add dependencies on com.sun.xml.bind:jaxb-impl and jaxb:activation to support use with Java 9 (with some warnings). Without this change, users must supply the Java VM argument "--add-modules=java.xml.bind,java.activation"
- Fix a type in the usage output of ComputeFDR (thanks to GitHub user Jong-hun-Park)

**v2017.08.23**

- MzIdentML creation: Change how the peptide ids are created, to further reduce possibility of duplicate peptide ids

**v2017.07.21**

- Performance improvements when reading mzML and mzIdentML files

**v2017.05.18**

- Add scoring parameters for UVPD TMT6Plex
- ScoringParamGen: Write and keep the converted TSV files to reduce time to re-run training on the same data
- ScoringParamGen: Write both the binary and text versions of the training results
- ScoringParamGen: Run multithreaded
- Improve the error message displayed when there are too many duplicate proteins

**v2017.01.27**

- Reduce the number of exceptions seen when a thread exits with an exception - There were multiple exceptions being displayed that were the result of killing the other threads.
- Change how the overall progress reporting and exception handling is managed, to minimize the thread overhead burden.
- Allow the user to specify a number of tasks on the command line using "-tasks n", and modify the default calculation. The previous default had a hard-coded limit of 64, and multiplied the number of threads by 10; new default removes the hard-coded limit, and multiplies by only 3.
- Change the output behavior - by default, the per-thread status messages will no longer be output. Per-thread status messages are a large amount of output, and aren't usually necessary with good progress reporting. The old output can be re-enabled using the command line argument "-verbose 1"
- Add the elapsed time to the overall progress report
- Fix the overall progress reporting - get the active task progress again with the overall task progress.

**v2017.01.10**

- Allow CustomAminoAcid formulas to not specify a number if it is 1.
- Updating UVPD scoring parameters

**v2016.12.12**

- Properly output the unitCvRef when outputting the scan start time. Since v2016.08.31

**v2016.12.08**

- Internal PeptideId change for mzid files  
  Last version had a possibility of duplicate peptide ID strings, when the same mod was possible on the first residue and the N-terminal, or on the last residue and the C-terminal. N-terminal mods are now added prior to the first residue, prefixed with '\[', and C-terminal mods are added after the last residue, prefixed by '\]'

**v2016.12.02**

- Fix some oddities with the mzid file peptide id strings

**v2016.11.29**

- Clean up some debugging messages, and limit the number of times in a single search

**v2016.10.26**

- Minimum spectra per thread reduced from 1000 to 250

**v2016.10.24**

- Return a non-zero exit code if an error occurs
- Better exception handling for multi-threaded searches

**v2016.10.14**

- Fix: handling of -m (FragmentMethodID) when processing .mgf files

**v2016.10.10**

- Fix: mzid output - Peptide IDs were being mishandled and led to creating incorrect PeptideEvidence references in SpectrumIdentificationItems. Also, added more information to the Peptide IDs to decrease the possibility of ID collision on edge cases.

**v2016.09.22**

- New: Add the ability to set the charge carrier mass to something besides the mass of a proton.

**v2016.08.31**

- Fix: output the scan start time units to the mzid file when the input is mzML; previously only the value (without the units specified) was being output, which was ambiguous and did not comply with the CV specification
- New: When input is a \_dta.txt file, with supporting \_ScanType.txt file, a 4th column listing scan start time (in minutes) is supported, and will be output to the mzid if present.

**v2016.07.26**

- Add UVPD as a dissociation method
- Add new ions: a. and x.
- Fix off-by-one error in IonType code that ignores c-NH3

**v2016.06.29**

- Clean up the mzid output a little, reducing file size (v2016.06.15 introduced a change that resulted in larger mzid files)
- Reduced the amount of unnecessary data in the .jar file, cutting size by ~1/2

**v2016.06.15**

- Fix the mzid output when the modification is unknown to unimod. (cvRef now correctly references PSI-MS, and the value will be the name provided in the Mods.txt file)
- Update the unimod.obo to date 2016:02:01, and change code to parse it properly after updating
- Bump some dependencies to a newer version, and change code to properly work with those changes.

**v2016.05.25**

- Output the residue letter and mass of any custom amino acids to the "MassTable" portion of the mzid.
- Update the URIs for the ontologies to their current locations (psi-ms.obo and unit.obo).

**v2016.02.12**

- Ensure that ETciD and EThcD are handled as ETD when using a \_dta.txt file with a \_ScanType.txt file

**v2016.01.29**

- Added ability to enter custom amino acids using the Mods.txt file.

**v2016.01.21**

- Changed the versioning system - SVN revisions don't work for Git repositories. Now commit date is used for the version.
- Improve the overall progress reporting.
- Cause entire search to fail if a thread exits with a failure.
- Retention Time/Scan Start Time is now output to the mzid file for searches performed on mzML files, as long as the data is available.

**2014-07-16 v10089**

- Fixed a bug that crashes when C-term mod mass is below -57Da.
- Added a test that checks whether the output path is valid before processing the data.

**2014-06-30 v10072**

- Optimization for multithreaded performance.
- Several minor bug fixes.

**2014-02-10 v9949**

- New scoring parameters are added for HCD/Q-Exactive/Trypsin/TMT. Parameters for HCD/HighRes/Trypsin/TMT have also been changed. As a result, for HCD spectra of TMT peptides, the number of identifications has been significantly increased.
- MS-GF+ now automatically recognize "TMT6plex" in the modification file and change the protocol to TMT.
- For mzML files converted from Thermo raw data using msconvert, MS-GF+ now reads monoisotopic precursor m/z first from a user parameter (\[Thermo Trailer Extra\]Monoisotopic M/Z:) instead of the CvParam "MS:1000744".
- The maximum number of variable modifications (per peptide) is written to the mzIdentML output file.
- Previously MS-GF+ ignored spectra having less than **20 peaks** for non-TOF spectra. Now it ignores spectra having less than **10 peaks**.

**2013-08-28 v9881**

- Change in database indexing format. The index file keeps non-standard amino acids (characters other than 20 standard residue characters). Previously, non-standard amino acids were converted into ‘?’ while indexing. It caused a problem when converting mzid into pepXML using idconvert (ProteoWizard).
- Fixed a bug that creates duplicate PeptideEvidence items when more than one matches are reported for a spectrum (e.g. -n 10),
- Fixed a bug that crashes MS-GF+ when a modification makes the amino acid mass 0 or smaller (e.g. M-131).

**2013-08-28 v9733**

- Added parameters for CID-LowRes-NoCleavage.
- Misc bug fixes.

**2013-04-03 v9517**

- Previously separate SpectrumIdentificationItems were created for the same peptide if "pre" is different (e.g. R.SIPDSMNYGDEEENK and K.SIPDSMNYGDEEENK). Now they show up as the same SpectrumIdentificationItem. If the score is different due to different NTT (e.g. G.SIPDSMNYGDEEENK), a separate SpectrumIdentificationItem is created.
- Bug fix: falling into an infinite loop while reading an mzML file if it contains \<binaryDataArray\> of encodedLength="0"

**2013-04-03 v9501**

- Previously, spectra are ignored in the search if the number of peaks is less than 20 for all types spectra. Now, for TOF spectra (i.e. -inst 2), this number has been changed to 3.

**2013-04-02 v9494**

- The following features are added in the SpectrumIdentificationItem when "-addFeatures 1"
  - MeanErrorAll: Mean of mass errors of explained fragment ion peaks
  - StdevErrorAll: Standard deviation of mass errors of explained fragment ion peaks
  - MeanErrorTop7: Mean of mass errors of 7 fragment ion peaks with highest intensities
  - StdevErrorTop7: Standard deviation of mass errors of 7 fragment ion peaks with highest intensities

**2013-03-25 v9436**

- Added TMT scoring model for HCD/HighRes (-m 3 -inst 1)
- Bug fix: duplicate PeptideEvidence id in mzid output
- Bug fix: crashing while using user-trained scoring parameters

**2013-03-05 v9324**

- Added Q-Exactive unlabeled phosphorylation parameters
- Bug fix: recognizing ETD+SA spectra in mzML as CID

**2013-02-27 v9324**

- Minor bug fix: mistakenly assigning 3 mods to a N-term amino acid.

**2013-02-15 v9312**

- Added scoring parameter sets for Q-Exactive iTRAQ (-inst 3 -protocol 2) and iTRAQ phosphopeptide enriched (-inst 3 -protocol 3) samples.

**2013-02-15 v9284**

- When "-protocol" parameter is missing, MS-GF+ automatically selects an appropriate protocol depending on the modification file.
- In the mzid output, Enzyme has attributes "missedCleavages".

**2013-02-14 v9249**

- Scoring parameters for Q-Exactive (-inst 3) have been added.

**2013-02-04 v9244**

- A bug (crash with an exception) in edu.ucsd.msjava.ui.ScoringParamGen has been fixed.

**2013-01-03 v9176**

- "-ti" parameter accepts only two comma separated integers, i.e., "-ti 1" will be rejected.
- MzIDToTsv adds "Title" column to the tsv file if the input spectrum is the mgf format.

**2012-12-19 v9107**

- "No enzyme" (-e 1) was renamed to "unspecific cleavage".
- "No enzyme (peptidomics)" (-e 9) was renamed to "no cleavage".
- \<Enzyme\>...\</Enzyme\> in the mzid output now shows right CVs for -e 1 and -e 9.
- The following bug has been fixed
  - Precursor errors are negated when MzIDToTsv was run.
  - C-term specific modifications do not show up in result files.

**2012-12-10 v9014**

- The following bug has been fixed
  - Ignoring some spectra when -n value is larger than 1

**2012-11-30 v9012**

- The following bugs have been fixed
  - Profile spectra should be ignored but not.
  - Peptides shorter than MinLength are reported.
  - Redundant proteins are reported in the same line in tsv files
  - Some SpectrumIdentificationResult elements have no SpectrumIdentificationItem element.
  - Different PeptideEvident elements have the same ID.

**2012-11-09 v8884**

- Bug fix: crashing while reading mzML files converted from wiff

**2012-11-09 v8873**

- Fixed a bug of ignoring N-term peptide when N-term Met cleaved peptide exists.
- Add scan numbers for NIST mgf files

**2012-10-30 v8806**

- Fixed a bug in edu.ucsd.msjava.misc.MS2ToMgf

**2012-10-29 v8792**

- Fixed bugs to ignore N-term fixed mods in the output
- Added edu.ucsd.msjava.misc.MS2ToMgf

**2012-10-11 v8719**

- Fixed bugs reporting NaN as additional features when spectrum has charge 0

**2012-10-04 v8605**

- Updates to handle multiple charge states in the ms2 format

**2012-10-03 v8597**

- Fixed the bug to erroneously report precursorMz in converted tsv files

**2012-09-26 v8540**

- Update MS2 parser to support <http://noble.gs.washington.edu/proj/crux/ms2-format.html>
- Old MS2 format (<http://fields.scripps.edu/sequest/SQTFormat.html>) is no longer supported

**2012-09-20 v8490**

- Fix the following bugs:
  - Same PeptideEvidence occurs more than once
  - SearchModification/SpecificityRules is always empty
  - Unknown modification cvParam does not have cvRef
  - SpectrumIdentificationResult has no SpectrumIdentificationItem
- MzIDToTsv -o TSVFile is now an optional parameter.

**2012-09-18 v8477**

- The extension of target/decoy concatenated database file has changed from .revConcat.fasta to .revCat.fasta.
- \*.revConcat.fasta: decoy proteins have prefix "REV\_" (used by MS-GFDB)
- \*.revCat.fasta: decoy proteins have prefix "XXX\_" (used by MS-GF+)

**2012-09-18 v8472**

- Bug fix: N-term or C-term residue-specific modifications (e.g. pyro-glu from Q) will have locations 1 (N-term) or length (C-term).
- Fix a bug to write SpectrumIdentificationResult twice for some spectra without precursor charges

**2012-09-17 v8449**

- Fix a bug in parsing Agilent mzML files

**2012-09-13 v8442**

- Fix the bug ignoring PSMs with IsotopeError=0, when MinIsotopeError was negative
- Change the '-ntt' parameter to set the minimum number of termini following the enzyme specificity rule
- Change the default value of '-ntt' from 1 to 2
- Fix the bug to set "fixedMod=true" for variable modifications
- Fix the bug not showing fixed modifications in the peptide list
- Fix the bug to show the location of N-term modification as 1 instead of 0
- Fix the bug to calculate massToChargeRatio as neutral masses instead of charged masses
- Replace the cvParam showing the dissociation method with a userParam named AssumedDissociationMethod in the SpectrumIdentificationItem

**2012-08-30 v8299**

- Provide MzID to TSV converter (edu.ucsd.msjava.ui.MzIDToTsv)
- Ignore profile spectra in the analysis
- Support JRE 1.6

**2012-08-29 v8297**

- Fix minor bugs
- Add the scan number CV in the output

**2012-08-27 v8283**

- Fix minor bugs
