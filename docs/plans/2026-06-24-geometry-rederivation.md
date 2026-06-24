# Geometry re-derivation — making andes's partition/segment geometry its own

Status: design → implementation (branch `chore/independence-100`).
Goal: remove the last MS-GF+ structural carryover so the independence claim can
truthfully say "shares no code, constants, **or geometry**." Result-changing →
every model change is benchmark-gated at 1% true entrapment-FDP before adoption.

## 1. The problem (what is still seed-inherited)

`Estimator::estimate(counts, template)` (`crates/model-train/src/estimate.rs:117`)
learns only the *table values* (`rank_dist_table`, error/existence tables,
`charge_hist`). It copies **all geometry from the seed `template`**:

| id | quantity | source today | file |
|----|----------|--------------|------|
| G1 | `num_segments` (=2) | `template.num_segments` | estimate.rs:180; also in the LLR denominator `rank_scorer.rs:92` |
| G2 | segment formula `floor(mz/parent_mass·N)` | hardcoded | param_model.rs:154-160 |
| G3 | `max_rank` (=150) | `template.max_rank` | estimate.rs:156,185 |
| G4 | partition keys + per-charge `parent_mass` tier boundaries | `template.frag_off_table.keys()` | estimate.rs:168,221-237 |
| G5 | charge span | data-derived (✓), template fallback | estimate.rs:165 |
| G6 | `frag_off_table` ion-type membership per partition | `template.frag_off_table.clone()` | estimate.rs:184 |

Prior art (sourced): the **top-150 peak cap** is field-standard (Sage
`max_peaks`=150, MSFragger `use_topN_peaks`=150, *verbatim*); coarse m/z
segmentation is generic. The **distinctive** MS-GF+ element is the
*precursor-mass-quantile × charge × segment partitioned trained tables*. So the
real independence work is to **derive the partition geometry from andes's own
corpus** instead of copying the seed's — even if optimization re-lands on the
same numbers (legitimately convergent), it must be *derived*, not *inherited*.

## 2. Architecture — derive a geometry-only template, reuse the learner

Key realization: `estimate()` already takes the geometry from a `template`
`Param` and learns everything else. So we do **not** rewrite the learner. We add
a function that *constructs the template from the corpus* and pass it in place of
the seed:

```
derive_geometry(corpus_stats, GeometryConfig) -> Param   // geometry only, empty tables
        │
        ├─ charge span + per-charge parent_mass tier boundaries  (G4/G5, from data)
        ├─ num_segments, max_rank                                (G1/G3, from config — sweepable)
        ├─ partition skeleton = {charge} × {mass-tier} × {0..num_segments-1}
        └─ frag_off_table: per partition, candidate b/y ions (chemistry: charges
           1..=precursor_charge, offsets = PROTON / H2O+PROTON) kept when observed
           match-frequency > ion_freq_threshold                 (G6, from data)
                  │
                  ▼
         accumulate(derived_template) → estimate(counts, derived_template)   // UNCHANGED
```

`Partition { charge: i32, parent_mass: f32, seg_num: i32 }` where `parent_mass`
is the tier **lower bound** (floor lookup in `find_partition`,
param_model.rs:97). `FragmentOffsetFrequency { ion_type, frequency }`. Ion
offsets (PROTON, H2O+PROTON) are chemistry, already CODATA-sourced — not IP.

## 3. The derivation (concrete)

`GeometryConfig { num_segments: i32, max_rank: i32, n_mass_tiers: usize, ion_freq_threshold: f32 }`
(defaults that reproduce seed-like geometry: `2, 150, ~4, 0.15`).

**Step A — charge span + mass tiers (G4/G5).** From the corpus labeled-PSM
`(charge, parent_mass)` set:
- charge span = observed min/max charge (already what `charge_range` does).
- for each charge, parent_mass tier lower-bounds = **equal-occupancy quantiles**
  of that charge's precursor-mass distribution into `n_mass_tiers` bins. (This is
  the "data-derived boundaries" the first audit wrongly assumed already existed.)

**Step B — partition skeleton.** Cartesian product
`{charge} × {tier lower-bounds for that charge} × {0..num_segments-1}`, sorted by
the `Partition` lex order (the loader/`find_partition` invariant).

**Step C — frag_off_table ion membership (G6).** Candidate ions per partition =
b/y (`Prefix`/`Suffix`, loss_class 0) at fragment charges `1..=max(1, charge-1)`
(chemistry rule, not IP). Run a **light count pass** over the corpus bucketing
each matched fragment into its partition; keep an ion type for a partition when
its observed match-frequency `> ion_freq_threshold`; store the frequency. Always
include a `Noise` entry (RankScorer requires it per populated partition).

**Step D — assemble** a geometry-only `Param` (derived G1–G6 + empty learned
tables, `rebuild_cache()`), to hand to `accumulate`/`estimate`.

Steps A+C need one light pass over labeled PSMs (no full search). Once the
template exists, the existing `accumulate → estimate` flow is untouched.

### Pinned integration details (verified in-tree)

- **Corpus `parent_mass`** = `peptide.mass()` (f64) — the same `neutral_mass =
  (precursor_mz − PROTON)·charge` the partitioner uses (`scored_spectrum.rs:347`);
  `LabeledMatch { peptide, charge }` (`labeled.rs:73`) supplies both. Step A
  iterates the label set → `(charge, peptide.mass())`.
- **Ion set (G6)** is chemistry, not seed data: per partition of precursor charge
  `C`, emit `Prefix`(b) and `Suffix`(y) at fragment charges `1..=max(1, C−1)`,
  `loss_class 0`, offsets `model::mass::PROTON` (b) and `H2O+PROTON` (y), plus a
  `Noise` entry (RankScorer requires one per populated partition). Frequency:
  from a count pass, or uniform-then-learned (the rank/existence tables carry the
  real distribution either way). Keep the threshold prune (>0.15) as a later
  refinement.
- **Non-geometry metadata** (`data_type`, `mme` tolerance, `apply_deconvolution`,
  `version`, `precursor_off_map`) is *not* geometry/IP — `derive_geometry` takes
  it from a base `Param` (the seed, or a minimal config). Only G1/G3/G4/G6 are
  derived. (A fully config-driven metadata path is a later step.)
- **Train wiring**: `train` (`andes.rs:2697`) builds its template via
  `derive_geometry(&labels, &seed_param, geo_cfg)` and passes it to
  `RankScorer::new` / `accumulate` / `estimate` unchanged. **No CLI flag** — see §4.

## 4. Integration — derived geometry is the DEFAULT, no new CLI flag

Geometry being own-derived is the *goal*, not an option, so it ships as the
default `train` behavior — not an opt-in flag (a flag nobody should ever turn
off, and pure bloat on the just-minimized CLI). Concretely:

- `train` always builds its template via `derive_geometry(&labels, &seed_param,
  geo_cfg)`; the seed supplies only **non-geometry** metadata (tolerance /
  activation / deconv / version) until seedless training removes even that. There
  is **no `--derive-geometry` flag**.
- The seed-vs-derived A/B needed to *validate* the change is driven by the
  **benchmark harness / a throwaway internal switch** (an undocumented `ANDES_*`
  env), deleted once derived geometry is the validated default. It never becomes
  CLI surface.
- `GeometryConfig` (`num_segments`, `max_rank`, `n_mass_tiers`) is a **sweep-harness**
  concern, not CLI flags; once the §5 sweep picks the optimum it collapses to
  fixed defaults in code.
- Flipping the default is **result-changing for every retrain**, so it is gated:
  derivation becomes the default only after the §5 sweep proves derived ≥ seed at
  1% entrapment-FDP on all 3 regimes. Until then it lives behind the internal
  switch.
- Bonus: making geometry derivable also closes the "seed structurally required"
  gap, enabling fully seedless training later.

## 5. Sweep (after the code lands) — benchmark-gated, safest-first

Prototype on 3 regimes spanning the axes: `hcd_astral_tryp` (high-res),
a low-res TMT (**regression-risk canary**), UPS1 (`hcd_qexactive_tryp`,
FDP-honest). Most variants need only `accumulate+estimate` (minutes), not a
re-search. Order:
0. anchor (reproduce current numbers)
1. `max_rank ∈ {100,125,150,200}` (cheapest; expect flat ≈150, Sage-corroborated)
2. mass-tier count + equal-occupancy boundaries `∈ {3,4,5,6}` (the independence win; gate hard on the low-res canary for sparsity regression)
3. `num_segments ∈ {1,2,3,4}` (riskiest — it's in the LLR normalizer; per-regime allowed, expect 2 low-res / maybe 3 Astral)
4. segment-formula variant (optional)

Gate: wins-or-ties high-res **and** no low-res regression, at entrapment-FDP ≤1%.
Expected: values mostly re-land on the MS-GF+ numbers — but now **derived from
andes's corpus**, satisfying the independence gate (and possibly a real Astral
win at `num_segments=3`).

## 6. Testing (TDD)

- unit: equal-occupancy quantiles split a known mass distribution into
  balanced-occupancy tiers; partition skeleton has `|charges|·n_tiers·num_segments`
  keys, sorted, lex-correct; frequency threshold includes/excludes the right ions.
- integration: `derive_geometry → accumulate → estimate` yields a `Param` that
  `RankScorer::new` accepts (every populated partition has a Noise entry; all
  rank vectors strictly positive) and that scores a fixture spectrum.
- parity sanity: deriving with seed-matching config produces geometry of the same
  shape/cardinality as the seed (not byte-identical — the boundaries are now
  data-derived).
- adoption gated on the §5 benchmark, never on target-decoy alone.

## 7. Independence outcome

Lands the last structural carryover: after this, geometry is computed from
andes's own data, satisfying the Phase-4/E3 gate ("no MS-GF+ table values **and**
no MS-GF+ partition geometry"). Only then strengthen `NOTICE` to add "constants
or geometry" to the shared-nothing claim.
