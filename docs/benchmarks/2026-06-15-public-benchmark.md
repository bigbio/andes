# Public benchmark — andes vs the open-source field (2026-06-15)

Head-to-head of **andes** (top-1 and `--chimeric`), **Java MS-GF+ v20240326**,
**Sage 0.14.7**, **Comet 2025.01** (via OpenMS), and **ProSE** (OpenMS), across
three datasets spanning high-res HCD and low-res ion-trap CID. Every engine is
re-scored through **one uniform Percolator** (3.7.1, `--seed 42 -Y`), on the same
8-thread VM. FDR honesty is verified independently with a 1:1 entrapment search.

## Results — PSMs / peptides @ 1% FDR

| Engine | Astral (high-res HCD) | TMT a05058 (low-res CID) | UPS1 (low-res CID) |
|---|---:|---:|---:|
| **andes** (`--chimeric`) | **69,968 / 28,442** | **12,043 / 10,654** | **17,879 / 4,706** |
| **andes** (top-1) | **36,782 / 23,898** | **11,957 / 10,636** | 17,143 / 4,451 |
| Java MS-GF+ v20240326 | 26,542 / 17,954 | 11,555 / 9,863 | 17,305 / 4,421 |
| Sage 0.14.7 | 32,091 / 21,398 | 11,232 / 9,481 | 15,653 / 4,202 |
| Comet 2025.01 | 31,435 / 20,608 | 10,876 / 9,290 | 15,809 / 4,219 |
| ProSE (OpenMS) | 30,590 / 20,590 | 7,659 / 7,066 | 8,901 / 2,960 |

- **Astral (high-res HCD):** andes top-1 alone beats every competitor on both PSMs
  and peptides; `--chimeric` (co-isolated second peptides) is +118% over the next
  engine. Reproduced from native `.raw` and converted mzML alike.
- **TMT a05058 (low-res ion-trap CID, TMT-labeled):** andes top-1 leads on PSMs
  **and** peptides. ProSE underperforms — it caps fragment tolerance at 0.1 Da and
  is designed for high-res fragmentation.
- **UPS1 / PXD001819 (low-res CID LFQ):** andes top-1 is within 1% of Java on PSMs
  and ahead on peptides; `--chimeric` takes the PSM lead.

Wall time: andes finishes each run in ~1–4 min vs Java MS-GF+'s 9 min – 2.5 h
(≈10–40×), on par with the C++/Rust engines (Comet/Sage 1–4 min). andes is the
only engine here that reads Thermo `.raw` and Bruker timsTOF `.d` natively.

## The 1% FDR is real — entrapment validation

Target-decoy q-values are self-consistent *by construction*, so they can't prove
FDR honesty alone. We checked it independently with a **1:1 entrapment** database
on Astral (real HYE proteins + an equal set of `ENT_` foreign sequences; andes
adds `XXX_` reversed decoys as usual). A target PSM mapping only to `ENT_`
proteins is a known false positive, so **true FDP = 2 × ENT-hits / total**,
*independent* of the decoy machinery.

| q ≤ | top-1 true FDP | chimeric true FDP |
|---|---:|---:|
| 0.5% | 0.54% | 0.52% |
| **1.0%** | **1.06%** | **1.14%** |
| 2.0% | 2.27% | 2.34% |
| 5.0% | 5.20% | 5.55% |

True FDP tracks the nominal q across the whole range; at the 1% line the actual
error is ~1.1% (very slightly optimistic, a benign Percolator characteristic —
not an FDR violation). Crucially the **chimeric near-doubling holds at 1.14%**,
so the extra co-isolated IDs are genuine. The non-tryptic LysC and GluC+Trypsin
runs hold ≈1% true FDP too.

## Datasets & parameters

| | Astral | TMT a05058 | UPS1 / PXD001819 |
|---|---|---|---|
| Run | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1` | `a05058` (PXD007683) | `UPS1_5000amol_R1` |
| Instrument | Orbitrap Astral, high-res HCD | ion-trap CID-MS2, TMT, low-res | LTQ-Orbitrap, ion-trap CID, low-res |
| Database | ProteoBench HYE (human+yeast+E.coli) | human+yeast reviewed | yeast + UPS1 |
| Precursor tol | 10 ppm | 20 ppm | 20 ppm |
| Fragment tol | 20 ppm (high-res) | 0.4 Da (low-res; ProSE 0.1 Da) | 0.4 Da (low-res; ProSE 0.1 Da) |
| Fixed mods | Cam-C | Cam-C, TMT6plex (K + N-term) | Cam-C |
| Variable mods | Ox-M, Acetyl (protein N-term) | Ox-M | Ox-M |

Common: trypsin, fully specific, ≤2 missed cleavages; length 6–40 (7–40 Astral);
charge 2–4; isotope error −1…+2; FDR 1% via Percolator 3.7.1 (`--seed 42 -Y`).

## Methodology notes

- **Uniform Percolator.** Each engine emits a Percolator PIN; all PINs go through
  the same `quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2` (`--seed 42 -Y`),
  counts at q ≤ 0.01. PSMs and unique peptides derived from the same target list.
- **PIN building.** andes / Sage / Comet write PIN directly. Java MS-GF+ →
  `MzIDToTsv` + `build_pins.py java` (its concatenated-TDA mzid crashes
  `msgf2pin`). ProSE → OpenMS idXML → `build_pins.py prose`.
- **Java MS-GF+ is deterministic.** Comet, Sage, and ProSE each reproduced their
  prior numbers to the digit; the Astral Java count is reused from a prior run
  (its `msgf2pin` step crashes regardless of input, and the count is
  pin-builder-independent).
- **Protein counts are intentionally omitted** from the comparison: a fair
  protein-level comparison needs uniform parsimony grouping, since raw
  `proteinIds` lists differ by engine output format (e.g. ungrouped accession
  lists inflate the unique count well past the peptide count).
- **Precursor calibration** is off (the andes default).

## Reproduce

VM scripts: `scripts/bench_astral_competitors.sh`,
`scripts/bench_tmt_ups_competitors.sh` (andes + Java + Sage + Comet + ProSE →
uniform Percolator), `scripts/astral_entrapment_experiment.sh` +
`scripts/entrap_fdp.py` (the entrapment-FDP validation). Per-engine configs in
[`configs/`](configs/).
