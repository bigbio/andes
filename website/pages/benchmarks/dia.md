# DIA benchmark: ProteoBench Astral

andes `--dia` was run on the dataset and evaluation that OpenDIA uses to compare itself with
DIA-NN, so the three tools are counted the same way.

| | |
|---|---|
| Dataset | ProteoBench `LFQ_Astral_DIA_15min_50ng_Human`, runs 01 and 02: Orbitrap Astral, 300 isolation windows of 2 Th over 380–980 m/z, 1.65 s cycle, about 293,000 MS2 scans per run |
| Database | OpenDIA's PAIRLIB FASTA: 100,634 one-peptide entries, human targets and an equal number of `_p_target` entrapment peptides |
| andes | `--dia --rescore --predict rt,ms2 --enzyme nocleavage --charge 1..5`, Carbamidomethyl (C) fixed, Oxidation (M) variable, Percolator 3.7.1, 2026-11-17 |
| DIA-NN | 2.2 with its own decoys on the PAIRLIB targets; OpenDIA arm E from OpenDIA's results |
| Counted as | precursors (peptidoform, charge) of the PAIRLIB library, each scored by its best PSM |

## Precursors at a true FDP of 1%

![Precursors at an entrapment true FDP of 1%](../figures/dia_astral.svg)

Each tool's precursors are ranked by its own score and counted up to an entrapment-estimated
false discovery proportion of 1%, FDP = E (1 + 1/r) / N with r the entrapment-to-target ratio of
the library. This compares ranking quality at a matched true error, independent of how each
tool calibrates its q-values. The andes bars are the mean of five Percolator seeds (42, 1, 2, 3,
7); the error bars span the seeds.

| | raw01 | raw02 | at 2% true FDP | nominal q ≤ 0.01 | its true FDP |
|---|---:|---:|---:|---:|---:|
| andes `--dia` | 1,777 | 1,797 | 1,886 / 1,926 | 1,585 / 1,719 | 0.3% / 0.7% |
| DIA-NN 2.2 | 2,067 | 2,301 | 2,522 / 2,543 | 2,361 / 2,322 | 1.4% / 1.0% |
| OpenDIA arm E | 1,694 | 1,628 | 1,908 / 1,903 | 1,667 / 1,457 | 1.0% / 0.7% |

andes reaches 86% / 78% of DIA-NN and 105% / 110% of OpenDIA arm E. Its nominal 1% is
conservative on both runs (0.1–0.4% and 0.5–0.9% over the seeds).

## Time

Wall time per run, 14 threads: andes 5.7–6.8 min, of which the search, PIN and DIA columns take
about 60–65 s and the two Percolator passes 4.4–5.4 min; DIA-NN 8.8 min (raw02) and 13.5 min
(raw01). The two tools ran on different machines (andes on a cluster node, DIA-NN in Docker on a
workstation), so the times are indicative only.

## Where the gap to DIA-NN is

93% of DIA-NN's raw01 precursors have a row in the andes PIN, and for 84% of them the best row is
inside DIA-NN's chromatographic peak. The ones andes misses are weak: their spectra are faint and
the precursor spans only 1–3 cycles, too short for chromatogram correlation to separate them from
noise. andes reaches 25% of DIA-NN's lowest-intensity quartile and 80% of its highest. A
peptide-centric peak-group pass tested on this set did not beat the spectrum search (284 against
285 precursors on a 3-minute slice). Every variant tested, with its numbers, is in the
[DIA benchmark notes](https://github.com/bigbio/andes/blob/main/docs/benchmarks/dia/README.md).
