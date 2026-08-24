# andes vs Comet — rebuilt 3-dataset benchmark (2026-08-23)

Rebuilt from scratch: the previous benchmark data had been purged from scratch storage on
both hosts, so every raw file, database and result below was re-acquired and re-run.

## 1. Datasets

| benchmark | accession | file | database |
|---|---|---|---|
| Astral (high-res LFQ) | PXD070049 | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` (2.58 GB) | ProteoBench mixed-species HYE (31,889 seqs) |
| TMT (low-res ion-trap CID) | PXD007683, run `a05058` | `a05058.raw` (0.54 GB) | UniProt human + yeast reviewed (26,483 seqs) |
| UPS1 (low-res) | PXD001819 | `UPS1_5000amol_R1.mzML` | `yeast_entrap.fasta` (entrapment) |

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
| TMT | 193 s | **117 s** (K=100) | 80 s |

**Comet was 3.9x faster on Astral and 2.0x on TMT** as originally measured. This contradicts
the claim in earlier benchmark notes that andes is the fastest engine; that figure was for
native `.raw` input, whereas these runs are mzML. The numbers above are as measured.

Two optimisations closed most of the gap (see `docs/` history and the commits):

1. **GBDT de-duplication** — the fragment-intensity ensemble was being walked twice per
   candidate for identical inputs. Removing the second walk is **byte-identical output** and
   gave −38.3% (Astral) / −23.8% (TMT).
2. **`--gbdt-max-trees`** (default off) — truncating the ensemble to 100 of 300 trees gave a
   further −18.2% on TMT with **no measurable identification cost** (5 seeds; largest deviation
   7.8 PSMs against a control spread of 25.7). Measured on TMT only; not yet validated on
   high-res, which is why it ships off.

Gap after optimisation: **2.6x (Astral)**, **1.5x (TMT)**, with the identification lead intact.

## 5. Reproducing

Scripts used live on the benchmark VM under `/srv/data/andes-bench/bench2/`:
`fetch.sh` (download + ThermoRawFileParser conversion), `run_andes.sh`, `run_comet.sh`,
`score.sh` (uniform Percolator + counting), `seeds.sh` / `ksweep.sh` (multi-seed sweeps).
