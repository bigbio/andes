# Soft fragment matching

> **Status: on by default, parameter-free.** There is no flag and no knob — the
> softening width is the model's own match tolerance.

## What it does

Standard fragment matching is a hard cutoff: a theoretical ion's m/z is matched
to the highest-intensity observed peak within the fragment tolerance, and every
peak inside that window contributes its full rank log-likelihood ratio (LLR)
regardless of *where* in the window it sits. At low resolution (a wide tolerance,
e.g. 0.5 Da ion-trap CID) this rewards a chance/noise match near the tolerance
edge as much as a true match at the window centre — the source of low-resolution
score brittleness.

Soft fragment matching removes that cliff. A matched peak's contribution is
weighted by a Gaussian of its mass error and blended toward the *missing-ion*
score, so a centred peak keeps full credit while an off-centre (likely-noise)
peak is smoothly discounted. The score becomes a smooth function of mass error
instead of a step function of the tolerance.

## Formula

For a matched ion, let:

- `matched` = the hard per-ion score (rank LLR + any per-peak GBDT term),
- `missing` = the model's missing-ion ("absent") score for that ion,
- `Δm = peak_mz − theo_mz` = the matched peak's mass error,
- `σ = tol_da` = the Gaussian width = the model's own effective match tolerance
  at that m/z.

Then the per-ion score is:

```
w = exp(−½ (Δm / σ)²)
score = w · matched + (1 − w) · missing
```

A peak at the window centre (`Δm = 0`) keeps full credit (`w = 1`); a peak at the
tolerance edge (`Δm = σ`) keeps `w = exp(−½) ≈ 0.61`; peaks farther out approach
`missing`.

## Why parameter-free (no `σ` knob)

`σ` is the model's own match tolerance (`mme`), not a tuned coefficient. This:

- **scales per regime automatically** — a low-res model carries a wide tolerance
  (wide σ, meaningful softening); a high-res model deconvolves to a tight window
  (Δm ≪ σ ⇒ `w ≈ 1`), so soft matching is ~inert on high-res and **needs no
  per-model configuration**;
- **adds no parameter** — there is nothing for a user to set; the one tolerance
  the model already has drives the softening.

A *tuned* width (e.g. σ = tolerance / 2) extracts a larger low-res gain but trades
it off against other regimes (it regresses low-res TMT); the parameter-free
σ = tolerance is the setting that is **net-positive across regimes with no knob**.

## Validation (1% true entrapment-FDP, shipped models)

| dataset | hard matching | soft (σ = tolerance) |
|---|---:|---:|
| UPS1 (low-res LFQ) | 14,948 | **15,061 (+0.8%)** |
| TMT a05058 (low-res CID) | 11,125 | **11,163 (+0.3%)** |
| Astral (high-res HCD) | 36,700 | **36,873 (+0.5%)** |

Net-positive on **all three** regimes — no regression anywhere, no parameter.
