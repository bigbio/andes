# Mass spectrometry correctness

## Constants and masses

- Use the constants in `crates/model/src/mass.rs` (proton, water, element and residue
  masses); never retype a mass literal (see `parameters-and-constants`). `crates/model/tests/chemistry_constants.rs` guards
  them.
- Distinguish monoisotopic and average mass, neutral mass and m/z, `[M+H]+` and `M`.
  Name variables accordingly (`neutral_mass`, `mz`, `mh`).
- `m/z = (M + z * proton) / z`. A charge must be at least 1 before dividing.

## Tolerances

- ppm is relative to the measured value: `da = mz * ppm * 1e-6`. Never mix a ppm window
  with a Da window in one comparison.
- Fragment tolerance in training is rebuilt from nominal masses (TRAIN.md); do not switch
  it to ppm there.
- Isotope errors (-1..+2 13C peaks) shift the precursor by multiples of 1.003355 Da / z.

## Modifications

- Resolve names and masses through Unimod (`resources/unimod.obo`, `unimod_aliases.tsv`);
  check site and position together (N-term Q for pyro-Glu, protein N-term for acetyl).
- Fixed vs variable, label channels (TMT, SILAC, dimethyl) and their limits are counted
  separately from other variable mods.
- Report modifications with Unimod accessions in outputs.

## Peptides and decoys

- Digestion rules: cleavage sites, missed cleavages, non-specific ends, protein N-terminal
  Met removal. Handle `X`, `B`, `Z`, `U`, `O` residues explicitly.
- Decoys must have the same length, composition and mass distribution as targets; a decoy
  change needs an entrapment check (`fdr-validation`).
- I/L are isobaric: peptide-to-protein mapping and inference must treat them as equal.

## Spectra

- Activation (HCD, CID, ETD, EThcD) decides ion series (b/y vs c/z•). Instrument
  resolution decides centroiding and tolerances.
- Missing or zero charge, unsorted peaks and profile-mode data are real inputs; handle
  them, do not assume clean files.
- timsTOF: ion mobility per precursor; FAIMS: compensation voltage per scan. Keep them
  attached to the spectrum through every stage.

## Test

Hand-compute at least one expected mass, m/z or fragment ladder in the test, from a source
you can cite, rather than asserting the code's own output.
