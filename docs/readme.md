# MS-GF+ documentation (Markdown)

Static HTML under `docs/` was replaced with these Markdown pages so they read well on GitHub and in editors.

### Contacts

- PNNL Proteomics: proteomics@pnnl.gov
- Sangtae Kim: sangtae.kim (at) gmail.com

### Summary

- MS-GF+ is an MS/MS database search tool that is sensitive (it identifies more peptides than other database search tools and as many peptides as spectral library search tools) and universal (works well for diverse types of spectra, different configurations of MS instruments and different experimental protocols).
- Input: HUPO PSI standard mzML and MGF only (mzXML, MS2, PKL, and `_dta.txt` are not supported in this fork).
- Output: Percolator `.pin` (default, for rescoring) or TSV. mzIdentML (`.mzid`) output has been removed — MS-GF+ now feeds downstream Percolator pipelines directly via `.pin`. See [Changelog](changelog.md) for migration notes.

### Usage and help

- [MS-GF+ usage](msgfplus.md) — [Change log](changelog.md)
- [Output formats — pin / tsv column reference](output.md)
- [MS-GF+ parameter files](https://github.com/MSGFPlus/msgfplus/tree/master/docs/parameterfiles) (config examples on the upstream repo)
- [Suffix array builder (BuildSA)](buildsa.md)
- [Isobaric labelling: TMT / TMTpro / iTRAQ recipes](isobariclabeling.md)
- [Troubleshooting & common errors](troubleshooting.md)

### Publications

**MS-GF+ makes progress towards a universal database search tool for proteomics.**  
Sangtae Kim and Pavel A Pevzner. *Nat Commun.* 2014 Oct 31; 5:5277. doi: 10.1038/ncomms6277.  
[PubMed ID 25358478](https://pubmed.ncbi.nlm.nih.gov/25358478/)

**Spectral probabilities and generating functions of tandem mass spectra: a strike against decoy databases.**  
Sangtae Kim, Nitin Gupta, and Pavel A Pevzner. *J Proteome Res.* 2008 Aug; 7(8):3354-63. doi: 10.1021/pr8001244.  
[PubMed ID 18597511](https://pubmed.ncbi.nlm.nih.gov/18597511/)
