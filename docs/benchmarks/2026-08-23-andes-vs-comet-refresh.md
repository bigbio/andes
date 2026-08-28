# andes vs Comet — rebuilt benchmark, Astral + TMT (2026-08-23)

Rebuilt from scratch: the previous benchmark data had been purged from scratch storage on
both hosts, so every raw file, database and result below was re-acquired and re-run.

**Scope: the head-to-head covers Astral and TMT only.** UPS1 is listed in the dataset table
below because it was re-acquired alongside the others, but it was NOT re-run against Comet in
this refresh, so it appears in no results table. Do not read the two-row identification and
speed tables as a three-dataset comparison.

## 1. Datasets

| benchmark | accession | file | database |
|---|---|---|---|
| Astral (high-res LFQ) | PXD070049 | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` (2.58 GB) | ProteoBench mixed-species HYE (31,889 seqs) |
| TMT (low-res ion-trap CID) | PXD007683, run `a05058` | `a05058.raw` (0.54 GB) | UniProt human + yeast reviewed (26,483 seqs) |
| UPS1 (low-res) — *re-acquired, not re-run here* | PXD001819 | `UPS1_5000amol_R1.mzML` | `yeast_entrap.fasta` (entrapment) |

The HYE FASTA is served from `proteobench.cubimed.rub.de/fasta/`, not GitHub.

> **Caveat.** The TMT database is a 2026 reconstruction of the depositors' 2018
> `PXD007683_UP000005640_UP000002311_reviewed.fasta`, not that file. Cross-engine comparison
> is valid (all engines share it); absolute TMT counts are not comparable to earlier tables.

## 2. Method

All engines: same mzML, same database, uniform Percolator 3.7.1 `--seed 42 -Y
--only-psms=false`, 1% FDR. Percolator mode reported by the log: **Separate target and decoy**.
Comet runs from `ghcr.io/bigbio/openms-tools-thirdparty` with the harmonized parameters of the
earlier benchmark (10 ppm precursor / 0.02 Da fragment bins on Astral; 20 ppm / 0.4 Da on the
low-res TMT run).

Protein counts are **distinct non-decoy protein IDs**, not the uniform-parsimony grouping used
in the 2026-06 table. Comet's Astral proteins therefore read 4,951 here vs 4,203 there; this is
a counting convention, not a change in behaviour.

**Harness faithfulness.** Comet's Astral PSMs reproduce the published figure exactly (31,435;
peptides 20,608 vs 20,607), and andes TMT peptides land at 10,690 vs a published 10,691.

## 3. Identifications @ 1% FDR

| dataset | engine | PSMs | peptides | proteins |
|---|---|---|---|---|
| Astral | **andes** | **38,437** | **24,981** | **5,641** |
| Astral | Comet | 31,435 | 20,608 | 4,951 |
| TMT a05058 | **andes** | **12,316** | **10,977** | **5,519** |
| TMT a05058 | Comet | 10,504 | 9,008 | 5,031 |

andes leads by **+22.3% PSMs / +21.2% peptides** on Astral and **+17.2% / +21.9%** on TMT.

## 4. Speed — and a correction to earlier claims

| dataset | andes (before) | andes (after opt) | Comet |
|---|---|---|---|
| Astral | 899 s | **555 s** | 217 s |
| TMT | 193 s | **143 s** default / 117 s with `--gbdt-max-trees 100` | 80 s |

**Comet was 4.1x faster on Astral and 2.4x on TMT** as originally measured (899/217 and
193/80 from the table above). This contradicts
the claim in earlier benchmark notes that andes is the fastest engine; that figure was for
native `.raw` input, whereas these runs are mzML. The numbers above are as measured.

Two optimisations closed most of the gap (see `docs/` history and the commits):

1. **GBDT de-duplication** — the fragment-intensity ensemble was being walked twice per
   candidate for identical inputs. Removing the second walk is **byte-identical output** and
   gave −38.3% (Astral) / −23.8% (TMT) in the dedicated A/B. (The headline table's own
   figures, 193 s → 143 s, imply −25.9%; the A/B and the headline runs are separate
   executions, and the gap is run-to-run variance, not a different optimisation.)
2. **`--gbdt-max-trees`** (default off) — truncating to 100 of 300 trees. **Superseded
   measurement (2026-08-27):** both regimes have since been swept at 5 seeds each against the
   current flag, which truncates both shipped ensembles:

   | regime | K=0 | K=100 | ΔPSMs | effect/SE | wall |
   |---|---|---|---|---|---|
   | Astral (high-res) | 38,436.8 | 38,444.4 | +7.6 (+0.020%) | +0.25 | 574 s → 338 s (−41%) |
   | TMT (low-res) | 11,935.2 | 11,900.6 | −34.6 (−0.290%) | −1.14 | 167 s → 99 s (−41%) |

   The speedup replicates exactly (−41% both) but the identification effect **does not**:
   Astral is flat, TMT's point estimate is negative. The TMT difference is below this
   design's resolution (~0.71% at n=5), but "below the detection floor" is not "zero", and
   the sign differs from Astral. **This is why the flag stays off by default** — buying 41%
   at a possible ~0.3% low-res identification loss, for every user, is the wrong default when
   they can opt in. The earlier "no measurable identification cost" claim rested on TMT alone
   and on a narrower version of the flag; it is superseded by the table above. The flag should be
   faster than −18.2% and its identification cost at that setting is UNMEASURED. It ships off.

   On high-res the flag is also riskier than on low-res: `--score auto` resolves to `strong`
   there, and `StrongScore` consumes the frag-intensity GBDT, so truncation changes the
   *ranking* and not merely PIN feature values. The Astral K-sweep confirms this — row counts
   shifted (1,213,883 → 1,213,884 at K=200) where TMT's were identical at every K.

Gap after optimisation **at shipped defaults** (de-duplication only; `--gbdt-max-trees` is off):
**2.6x (Astral, 555/217)** and **1.8x (TMT, 143/80)**, with the identification lead intact.
Enabling `--gbdt-max-trees 100` takes TMT to 1.5x (117/80), but that is not the default and is
not what a user gets out of the box.

## 5. Reproducing

Scripts used live on the benchmark VM under `/srv/data/andes-bench/bench2/`:
`fetch.sh` (download + ThermoRawFileParser conversion), `run_andes.sh`, `run_comet.sh`,
`score.sh` (uniform Percolator + counting), `seeds.sh` / `ksweep.sh` (multi-seed sweeps).
