# Open search on a mouse-liver fraction

`--open-search` accepts any database peptide within a mass-offset window of the observed
precursor and reports the offset as the modification. This page measures two things on one
real HCD file: **how many modified peptides the open search recovers** against a narrow search
that declares the modification, and **whether the nominal 1% q holds** once the window is wide.

The localization accuracy — the paper's own benchmark design, with a substituted database and
unmodified spectra — is measured separately in [`substitution/`](substitution/README.md):
90.0% of called sites are the substituted residue, and the shifted-ion shortlist recovers 20%
more truth peptides than the regular vote alone.

## Setup

| | |
|---|---|
| Spectra | `MouseLiver-Z-T-1.raw` (PXD005553), 45,905 MS2 scans, read natively |
| Database | `mouse_entrap.fasta` (mouse reviewed + peptide-level shuffled entrapment), T/E = 0.965 |
| andes | 0.4.0, branch `feat/open-search` (`7e6d71d4`, `d303e538`, `532064ab`, and the P4 changes in this commit) |
| Percolator | 3.7.1, `--seed 42 -Y --only-psms=false` |
| Host | 112 threads; searches 35 s (narrow) and 53 s (open, re-scores included), index build included |

Two searches of the same spectra and database:

- **truth** — the defaults: Carbamidomethyl-C fixed, Oxidation-M and protein-N-term-acetyl
  variable. Its oxidised-M PSMs are the ground truth.
- **open** — fixed mods only (`--mods` with Carbamidomethyl-C) plus `--open-search` with the
  default `-100..300` Da window and `--open-search-top-k 100`. Oxidation is then an unknown
  offset.

Both runs are scored by Percolator at `q ≤ 0.01`; the true FDP is the entrapment estimate with
the database's T/E factor (`(entrapment hits / accepted) × (1 + T/E)`).

## Results

| metric | truth (declared mods) | open (mass offset) |
|---|---:|---:|
| spectra with a PSM (pre-FDR) | 33,744 | 44,241 |
| PSMs at `q ≤ 0.01` | 2,652 | 4,260 |
| true FDP at `q ≤ 0.01` (entrapment) | 0.96% | 0.83% |
| of which, `\|offset\| ≤ 0.5 Da` | — | 1,714 rows, **0.00%** FDP |
| of which, `\|offset\| > 0.5 Da` | — | 2,546 rows, **1.39%** FDP |
| rows with a localised offset | — | 427 (10.0%) |

The open run finds more PSMs than the narrow one because its candidate set is one form per
peptide instead of every variable-mod combination; the truth run pays for oxidation and
N-terminal acetylation enumeration.

**Modification-aware scoring and the shifted-ion shortlist.** Every shortlisted candidate with a
non-zero offset is re-scored as the peptide with that mass on its best-localising residue, trying
the isotope offsets `0..=2`; the better of the plain and shifted scores is kept and all per-PSM
features are then computed on the shifted form. The inferred offset carries the precursor's mass
error, so it is first re-centred on the shifted fragment ions (a coarse window per residue, then
the residues are separated on the tight fragment score). The shortlist itself votes on the
complement (shifted) ions as well as the regular ones. This is the localization-aware open
search of Yu et al., *Nat. Commun.* 2020. On this file the pass ran **2,987,845 re-scores and
improved the node score 1,070,077 times (36%)**.

**Error control.** The zero-offset rows read none, but the offset-bearing rows now read **1.39%
true FDP at a nominal 1%** — above nominal. Measured on the same file with the regular shortlist
only, they read 0.72%: the shifted votes admit more offset candidates and lift the recovery of
the truth set (32 → 40 of 141, below), so this file shows the sensitivity/error trade-off of
admitting shifted evidence. 18 entrapment hits among 2,546 offset rows is about 1.5σ above the
nominal expectation, so read it as a direction, not a calibration. Open-search FDR is still not
*conditioned* on the offset by the engine — treat this as a measured example, not a guarantee.

**Recovery of a declared modification.** Of the 141 oxidised-M PSMs the narrow search finds at
1%, the open search reports **40 (28%)** with the same scan and bare peptide and an offset
within 0.02 Da of +15.9949; the median absolute offset error is **0.0019 Da**. 21 of the 40
carry a localised site and **20 are the truth's M**.

The form-based re-score alone (without the shifted votes) did not raise this number — 32 of 141
— and the reason is visible in the data: on many of these PSMs the **unmodified** form scores
higher on the fragments than the modified form (74 vs 69 node score on the scan inspected in the
run log), so the re-score correctly declines to replace it. The narrow search assigns the
modification from the precursor mass; the fragments alone do not prefer it. The shifted-ion
shortlist is what recovers the extra eight: it reaches modified peptides whose own ladder is
shifted, which the regular vote was dropping before scoring. A modification whose site has weak
b/y coverage is still not recovered by any form-based re-score.

**Isotopes.** The offset is measured from the monoisotopic precursor (`--precursor-mono` is
active in this mode; it shifted 447 of 45,905 precursors on this file). Without MS1 a recorded
M+1 peak is indistinguishable from a ~1.003 Da modification, so when the fragments give no
evidence for an isotope step the search reports offset 0 and `isotope_error`, rather than
biasing the reported mass down by one isotope.

**Where the offsets are.** Bin of the modification mass over every accepted open-search PSM
(1 Da bins; the count is not deduplicated by peptide):

| offset (Da) | rows | likely origin |
|---:|---:|---|
| 0 | 1,714 | unmodified |
| +1 | 399 | deamidation (+0.984) and noise |
| +17 / −17 | 125 / 79 | ammonia (+17.027) / loss |
| +53 | 87 | — |
| +16 | 63 | oxidation (+15.995) |
| +80 | 60 | phosphorylation (+79.966) |

## Reproduce

```bash
# release binary with native .raw support
cargo build --release -p andes --bin andes --features thermo
export DOTNET_ROOT=~/.dotnet

./run_open_bench.sh /path/to/andes-bench /path/to/MouseLiver-Z-T-1.raw
# T/E is a property of the database build; pass TE=<ratio> to skip the
# (minutes-long) rescan of every tryptic peptide of the FASTA.
```

`run_open_bench.sh` runs both searches, both Percolator runs and `analyze.py`, and writes
`report.json` with the counts above. The scripts take no hardcoded host paths.

## Limits

- One file, one Percolator seed; the oxidised-M ground truth is 141 PSMs.
- The sample is glycopeptide-enriched, so the peptide yield is low for both runs.
- The offset is not localised for most rows; localisation is a b/y heuristic (see `DOCS.md`).
- A declared modification is recovered only partly (23% here). The re-score does not raise it
  because the modified form is not fragment-supported on many of these PSMs; a precursor-δ
  prior or MS1-calibrated mono masses would be needed to prefer it.
