# MzIDToTsv

[MS-GF+ Documentation home](README.md)

There are two options for converting an MS-GF+ output file (**.mzid**) into a tab-separated file (**.tsv**)

1.  The MzIDToTsv utility built into MSGFPlus.jar
    - Easy to access (though syntax is a bit tricky; see below)
    - Can be slow for large .mzid files

2.  The Mzid-To-Tsv-Converter standalone application, [available on GitHub](https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases)

    - Fast conversion
    - Handles large .mzid files
    - Runs natively on Windows, but requires mono to use on Linux
    - Example command:

    `MzidToTsvConverter.exe -mzid:SearchResults.mzid -unroll -showDecoy`

## MzIDToTsv Utility

Converts MS-GF+ output (**.mzid**) into the tsv format (**.tsv**)


```text
Usage: java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv
    -i MzIDFile (MS-GF+ output file (*.mzid))
    [-o TSVFile] (TSV output file (*.tsv) (Default: MzIDFileName.tsv))
    [-showQValue 0/1] (0: do not show Q-values, 1: show Q-values (Default))
    [-showDecoy 0/1] (0: do not show decoy PSMs (Default), 1: show decoy PSMs)
    [-unroll 0/1] (0: merge shared peptides (Default), 1: unroll shared peptides)
```


**Parameters:**

- **-i MzIDFile**
  - Path to the MS-GF+ result file (\*.mzid)
- **-o TSVFile**
  - Path to the tsv output file (\*.tsv)
  - If not specified, for input MyFile.mzid, the output will be MyFile.tsv.
- **-showQValue 0/1**
  - If 0, QValue and PepQValue are not be reported.
  - If 1, QValue and PepQValue are reported (Default).
- **-showDecoy 0/1**
  - If 0, decoy PSMs will not be reported (Default).
  - If 1, decoy PSMs will be reported.
- **-unroll 0/1**
  - This parameter controls the output format for shared peptides (peptides matched to multiple proteins).
  - When `-unroll 0` (Default), a PSM matched to a shared peptide will be printed as a single line.
    - Peptide column does not show neighboring amino acids (e.g. QVHPDTGISSK).
    - Protein column shows all proteins in a single line.
    - Example: MyProtein(pre=K,post=T);MyProteinIsoform(pre=K,post=T)
    - [Download example file](examples/test.tsv)
  - When `-unroll 1`, a PSM matched to a shared peptide will be printed in multiple lines.
    - Peptide column shows neighboring amino acids (e.g. K.QVHPDTGISSK.A).
    - Different peptide-protein matches are printed in different lines.
    - [Download example file](examples/test_Unrolled.tsv)
