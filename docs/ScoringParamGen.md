# ScoringParamGen

[MS-GF+ Documentation home](README.md)

### Creates a scoring parameter file (UserParam) to be used for MS-GF+.

**Installation**

- ScoringParamGen is included in MSGFPlus.jar.


```text
Usage: java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen
        -i ResultPath (MSGFDBResultFile (*.mzid) or MSGFDBResultDir)
        -d SpecDir (Path to directory containing spectrum files)
        -m FragmentMethodID (0: CID, 1: ETD, 2: HCD)
        -inst InstrumentID (0: Low-res LCQ/LTQ, 1: High-res LTQ, Orbitrap, Lumos, 2: TOF, 3: Q-Exactive)
        -e EnzymeID (0: Trypsin, 1: Chymotrypsin, 2: Lys-C, 3: Lys-N, 4: glutamyl endopeptidase, 5: Arg-C, 6: Asp-N, 7: alphaLP, 8: no cleavage)
        [-protocol ProtocolID] (0: NoProtocol (Default), 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho, 4: TMT)
```


### Parameters:

- **-i ResultPath** - Required
  - Path to the folder containing MS-GF+ result files (\*.mzid or \*.tsv).
- **-d SpecDir** - Required
  - Path to the folder containing spectrum files used to generate MS-GF+ results.
- **-m FragmentationMethodID** - Required
  - Fragmentation method ID (0: CID, 1: ETD, 2: HCD)
- **-inst InstrumentID** - Required
  - ID specifying the instrument to measure product ions, i.e. MS/MS
  - 0: Low-res LCQ/LTQ, 1: High-res LTQ, 2: TOF, 3: Q-Exactive
- **-e EnzymeID** - Required
  - EnzymeID
  - 0: Trypsin, 1: Chymotrypsin, 2: Lys-C, 3: Lys-N, 4: glutamyl endopeptidase, 5: Arg-C, 6: Asp-N, 7: alphaLP, 8: no cleavage
- **-protocol ProtocolID** - optional
  - Protocol ID
  - 0: NoProtocol (Default), 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho, 4: TMT

**Output**

- A scoring parameter file containing scoring parameters (\*.param).
- The name of the scoring parameter file represents the spectrum type (FragmentationMethod, InstrumentType, Enzyme, Protocol) used to generate the spectrum.
  - FragmentationMethod_InstrumentType_Enzyme.param
  - FragmentationMethod_InstrumentType_Enzyme_Protocol.param
  - E.g. ETD_LowRes_GluC.param, HCD_HighRes_Tryp_Phosphorylation.param

**How to use a custom scoring parameter file?**

- Place scoring parameter files in the "params" directory.
- MS-GF+ uses the (user-defined) scoring parameter files to analyze spectra of the specified type.

**How to make a parameter file for a new fragmentation method or a new enzyme?**

- Define new fragmentation methods in [params/activationMethods.txt](examples/activationMethods.txt "activationMethods.txt") or new enzymes in [params/enzymes.txt](examples/enzymes.txt "enzymes.txt"). The params directory must exist below the working directory.
- Run ScoringParamGen with no parameter: "-cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen". It will now print out the ID of the new enzyme or the new fragmentation method.
- Run ScoringParamGen with new FragmentationMethodID or new EnzymeID

**What is "protocol" and how to define it?**

- Protocol is a user-defined tag associated with a scoring parameter file.
- For example, one may want to create a scoring parameter file for phosphorylation-enriched sample because the fragmentation propensities of phosphopeptides are different from other peptides even if the exact same condition (CID_LowRes_Tryp) is used to generate them. In such a case, one can define a protocol "Phosphorylation" and create CID_LowRes_Tryp_Phosphorylation.param.
- Define new protocols in [params/protocols.txt](examples/protocols.txt "protocols.txt").
- New protocols will show up in the usage information.

### Example

- Purpose: creating a scoring parameter file for HCD, HighRes, CNBr
- Add "CNBr,M,C,CNBr" in the params/enzyme.txt file.
- Run MS-GF+ with no parameter and check the id of CNBr (assume the id is 10).
- Run MS-GF+ by specifying -m 3 (HCD), -inst 2 (HighRes), -e 10 (CNBr). Suppose that the spectrum file name is test.mzML and the MS-GF+ result file is "test.mzid".
- Create a folder "results" and move test.mzid into it.
- Run ScoringParamGen:


```text
java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen -i results -d . -m 2 -inst 1 -e 10
```


- Check if "HCD_HighRes_CNBr.param" is created. Create a directory named "params" and move "HCD_HighRes_CNBr.param" into "params".
- From now on, running MS-GF+ with "-m 3 -inst 2 -e 10" will use HCD_HighRes_CNBr.param.
