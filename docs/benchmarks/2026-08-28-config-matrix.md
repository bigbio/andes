# Benchmark refresh — configuration matrix (2026-08-28)

Every number below was produced on one Linux benchmark host from a single build of
`main` at `e81df523`, on the same day. Comet figures are the 2026-08-23 rebuild, which
reproduced Comet's own published Astral count exactly (31,435 PSMs; 20,608 peptides
against a published 20,607) and is unchanged since.

FDR is Percolator 3.7.1 (`--seed 42 -Y`) at q ≤ 0.01 for the peptide rows. Where an
entrapment database exists, the **true** false-discovery proportion is reported next to
the count, because a q-value is a claim and an entrapment FDP is a measurement.

## 1. Peptide identification

| dataset | config | wall | PSMs @ q≤0.01 | vs default |
|---|---|---|---|---|
| Astral (high-res LFQ, PXD070049) | default | 450 s | 38,402 | — |
| | `--chimeric` | 478 s | **64,958** | +69.2% |
| | `--refine` | 627 s | 46,410 | +20.9% |
| | *Comet 2025.01* | *217 s* | *31,435* | *−18.1%* |
| TMT (low-res CID, PXD007683 `a05058`) | default | 124 s | 12,278 | — |
| | `--chimeric` | 82 s | 12,477 | +1.6% |
| | `--refine` | 123 s | 12,278 | skipped, see §3 |
| | *Comet 2025.01* | *80 s* | *10,504* | *−14.5%* |
| UPS1 (low-res LFQ, PXD001819) | default | 89 s | 15,938 | — |
| | `--chimeric` | — † | 17,294 | +8.5% |
| | `--refine` | 82 s | 15,938 | skipped, see §3 |

† The UPS1 chimeric run shared the host with another job, so its wall time is not
comparable and is omitted rather than reported misleadingly. Identification counts are
unaffected by concurrency.

**andes leads Comet by 22.1% (Astral) and 16.9% (TMT) on PSMs at matched FDR, and is
2.07x and 1.55x slower respectively.** Both halves are the result.

## 2. Is the `--chimeric` gain real?

`--chimeric` emits a structurally different PIN, so its q-values are computed over a
different population — worth checking before quoting the gain as identifications:

| | rows per scan | decoy fraction |
|---|---|---|
| Astral default | 9.98 | 47.3% |
| Astral `--chimeric` | 2.49 | 34.2% |

A 4x smaller candidate pool with a lower decoy fraction can inflate the count at a fixed
q-value without finding anything new. UPS1 ships an entrapment database (6,733 yeast
real + 4,531 `ENTRAP_` E. coli), so the question is decidable:

| UPS1 arm | PSMs @1% | entrapment hits | corrected FDP |
|---|---|---|---|
| default | 15,938 | 165 | **2.57%** |
| `--chimeric` | 17,294 | 170 | **2.44%** |

**The gain survives.** +1,356 PSMs at slightly *lower* true error, so the additional
identifications are real rather than an artifact of the changed decoy population.

Two honest qualifications:

- This validates the mechanism on **UPS1 only**. The Astral +69.2% is measured against a
  mixed-species database with no entrapment component, so it is **not** entrapment-validated
  and should not be quoted as though it were.
- Both UPS1 arms sit at **~2.5% true FDP against a nominal 1%**. That gap is present in the
  default configuration too — it is a property of this dataset and the rescoring, not
  something `--chimeric` introduced — but it means the q≤0.01 counts in §1 are optimistic
  by roughly 2.5x in absolute terms. They remain valid for comparing configurations and
  engines under one methodology, which is what the table is for.

## 3. `--refine` is high-res only, by design

On TMT and UPS1, `--refine` produced output identical to the default. That is correct
behaviour, not a silent failure: the run prints

```
WARN: refine is high-res-only and the data is low-res; skipping refinement.
```

and skips the second pass. On Astral it runs (14,813 confident-protein anchors over
101,664 unidentified spectra; the refinement phase took 178.7 s of the 627 s total).

The Astral `--refine` gain (+20.9%) is **not** entrapment-validated: the entrapment metric
is blind to a peptide-anchored second pass, which is why refine ships as a capability
rather than a headline number.

## 4. Intact N-glycopeptide search

Human plasma (PXD030622), three fractions searched separately and **pooled** before
scoring, then re-scored under **five Percolator seeds** — mandatory in this regime,
where `q_min = 1/T_top` makes the count at 1% a step function.

| | mean glycoPSMs @1% | sd | range | mean entrapment FDP |
|---|---|---|---|---|
| current `main` | 244.8 | 20.4 | 225–273 | 2.59% |
| prior run (2026-08-27) | 223.8 | 56.4 | 128–267 | 0.81% |

The +21.0 PSM difference is **effect/SE +0.78** — inside noise. At five seeds per arm this
design resolves only ~117 PSMs (~58% relative), so runs of this size cannot be separated
and should not be read as progress or regression.

The FDP is likewise seed-unstable (single seeds returned 8.69%, 4.26%, and three 0.00%),
so neither the 2.59% nor the earlier 0.81% is a reliable point estimate.

`--glyco-taxon auto` correctly excluded NeuGc on all three fractions, agreeing on both
evidence arms (oxonium ratio and FASTA taxon).

## Reproducing

Scripts used: `bench2/fullbench.sh` (the matrix) and `bench2/chim_entrap.sh` (the
entrapment check) on the benchmark host. They are not yet committed to the repository,
which is a real gap for anyone wanting to reproduce this from a checkout.
