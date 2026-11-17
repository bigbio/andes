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
| andes | 0.4.0, branch `feat/open-search` (pull request #397) |
| Percolator | 3.7.1, `--seed 42 -Y --only-psms=false` |
| Host | 112 threads; searches 39 s (narrow) and 84 s (open, re-scores included), index build included |

Two searches of the same spectra and database:

- **truth** — the defaults: Carbamidomethyl-C fixed, Oxidation-M and protein-N-term-acetyl
  variable. Its oxidised-M PSMs are the ground truth.
- **open** — fixed mods only (`--fixed-mods "Carbamidomethyl (C)" --variable-mods ""`) plus
  `--open-search` with the default `-100..300` Da window and `--open-search-top-k 100`.
  Oxidation is then an unknown offset.

Both runs are scored by Percolator at `q ≤ 0.01`; the true FDP is the entrapment estimate with
the database's T/E factor (`(entrapment hits / accepted) × (1 + T/E)`).

## Results

| metric | truth (declared mods) | open (mass offset) |
|---|---:|---:|
| spectra with a PSM (pre-FDR) | 33,745 | 44,415 |
| PSMs at `q ≤ 0.01` | 2,660 | 5,213 |
| true FDP at `q ≤ 0.01` (entrapment) | 0.66% | 0.79% |
| of which, `\|offset\| ≤ 0.5 Da` | — | 1,683 rows, **0.00%** FDP |
| of which, `\|offset\| > 0.5 Da` | — | 3,530 rows, **1.17%** FDP |
| rows with a localised offset | — | 569 (10.9%) |
| true FDP, zero-offset rows / offset rows | — | 0.00% / **1.17%** |

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
search of Yu et al., *Nat. Commun.* 2020. On this file the pass ran **12,931,409 re-scores and
improved the node score 4,168,971 times (32%)**.

**Error control.** The zero-offset rows read none, but the offset-bearing rows read **1.17%
true FDP at a nominal 1%** — above nominal. A development run with the regular shortlist only
read 0.72% on this file, so the shifted votes admit more offset candidates and lift the recovery
of the truth set (below); this file shows the sensitivity/error trade-off of admitting shifted
evidence. 21 entrapment hits among 3,530 offset rows is about 1σ above the nominal expectation,
so read it as a direction, not a calibration. Open-search FDR is still not
*conditioned* on the offset by the engine — treat this as a measured example, not a guarantee.

**Recovery of a declared modification.** Of the 138 oxidised-M PSMs the narrow search finds at
1%, the open search reports **41 (30%)** with the same scan and bare peptide and an offset
within 0.02 Da of +15.9949; the median absolute offset error is **0.0019 Da**. 22 of the 41
carry a localised site and **21 are the truth's M**.

The form-based re-score alone (without the shifted votes, on the pre-review code) did not raise
this number — 32 of 141 — and the reason is visible in the data: on many of these PSMs the **unmodified** form scores
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
| 0 | 1,683 | unmodified |
| +1 | 347 | deamidation (+0.984) and noise |
| +17 / −17 | 129 / 87 | ammonia (+17.027) / loss |
| +53 | 92 | — |
| +257 | 77 | — |
| −1 | 72 | deamidation (+0.984) / noise |
| +16 | 64 | oxidation (+15.995) |

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

- One file, one Percolator seed; the oxidised-M ground truth is 138 PSMs.
- The sample is glycopeptide-enriched, so the peptide yield is low for both runs.
- The offset is not localised for most rows; localisation is a b/y heuristic (see `DOCS.md`).
- A declared modification is recovered only partly (30% here). The re-score alone does not raise it
  because the modified form is not fragment-supported on many of these PSMs; a precursor-δ
  prior or MS1-calibrated mono masses would be needed to prefer it.
