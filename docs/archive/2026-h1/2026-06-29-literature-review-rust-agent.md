# Rust / Code-Aware Agent Brief — DDA Scoring Implementation

**Parent:** [`2026-06-29-literature-review-brief.md`](2026-06-29-literature-review-brief.md)  
**Web-research companion:** [`2026-06-29-literature-review-web-research-agent.md`](2026-06-29-literature-review-web-research-agent.md)  
**Repo root:** `msgf-rust/` (andes engine)

**Gate:** UPS1 low-res @ **1% paired entrapment-FDP** (+ no Astral/TMT regression).

---

## Current architecture (what exists)

### Score decomposition

```
RawScore (rank) = Σ_splits round(node_score(prefix, suffix)) + cleavage_credit + optional_loss
StrongScore     = fuse(intensity_signal, chance_match_surprise, mass_competition, entropy, listwise_gap)
TailorScore     = RawScore / Q99_denominator   (additive PIN only)
RawScoreCal     = z-score(StrongScore) per spectrum
```

**Primary ranking:** `RankScore` / `RawScore` (integer rank LLR path) — see `ScoreMode` in `andes.rs`.

### Code map (integration hooks)

| Concern | File | Lines (approx) | Notes |
|---------|------|----------------|-------|
| Integer RawScore sum | `crates/scoring/src/scoring/psm_score.rs` | 234–368 | `score_psm()` — per-split `round()` |
| Float rank (additive) | same | 370+ | `score_psm_float()` — no per-split round |
| Soft fragment LLR | `crates/scoring/src/scoring/scored_spectrum.rs` | 909–987 | `directional_node_score_inner` — σ from model mme |
| Tailor denominator | `crates/search/src/psm.rs` | 236–280 | `TAILOR_QUANTILE=0.01`, `tailor_denominator()` |
| Tailor histogram fold | `crates/search/src/match_engine.rs` | ~565, ~915, ~1033 | During candidate scoring |
| Tailor feature emit | `match_engine.rs` | 1073–1077 | `features.tailor_score = psm.score / tailor_denom` |
| Chance-match surprise | `match_engine.rs` | 1759–1775 | `ρ·Δ` local density null |
| Strong calibration | `crates/scoring/src/scoring/strong_score.rs` | 399–438 | LOO z-score vs `OnlineStats` null |
| PsmFeatures struct | `crates/search/src/psm.rs` | 12–234 | All PIN columns |
| PIN writer | `crates/output/src/pin.rs` | 135–443 | Column order contract |
| Percolator rescore | `crates/andes/src/rescore.rs` | — | GBDT fallback; production uses `--rescore` Percolator |
| Co-isolation pass | `crates/search/src/coisolation.rs` | 230–353 | Secondary Tailor/strong cal |

### Deliberately removed (do not reintroduce without FTO)

- **Generating function / SpecEValue** — US 8639447; stripped from PIN (`pin.rs` ~315–318 comments).

---

## Literature → code mapping

| Literature method | Engine status | Action |
|-------------------|---------------|--------|
| Tailor (Sulimov 2020) | **Shipped** (`TailorScore`) | A/B: does Percolator use it on low-res? |
| OMSSA Poisson E-value | Partial (`ChanceMatchSurprise`) | Compare formulas; consider `-log10` PIN column |
| Andromeda local q | Not shipped | Prototype windowed top-q on **rank landscape** |
| Exact XCorr DP (Howbert 2014) | Not applicable to integer LLR | See novel ideas below |
| MS-GF+ SpecProb | **Off-limits** | — |
| Keich MC 10K decoys (JPR 2014) | Not shipped | **Highest-value new feature** |
| Empirical null (RS³ Gate 0) | Design only | `docs/.../rs3-spectral-significance-design.md` |
| Percolator | Shipped | Keep PSM-level only |

---

## Implementation sequence (ordered by risk × expected PSM gain)

### Phase 0 — Measure (no new math)

1. **PIN feature importance** on UPS1 low-res: which of `TailorScore`, `ChanceMatchSurprise`, `RawScoreCal`, `DeltaRankScore`, `ListwiseScoreGap` does Percolator weight?
2. **Ablation:** run with/without `TailorScore` column at fixed entrapment gate.
3. **Confirm hypothesis:** MS-GF+ gap is calibration not rank model — compare score distributions per-spectrum variance.

Files: run benchmark harness; inspect Percolator weights from `--rescore` output.

### Phase 1 — Cheap calibration (1–3 days)

#### 1a. `NegLog10Tailor` / `LogTailorScore` PIN column

Percolator may prefer log-scale calibrated features (matches `lnSpecEValue` real estate removed).

- **Hook:** `pin.rs` after `TailorScore` definition (~378–382)
- **Math:** `-log10(1 / TailorScore)` or `-log10(empirical_rank)` — keep monotonic with TailorScore

#### 1b. Calibrate **StrongScore** path harder

`RawScoreCal` exists; ensure `strong_score_calibrated` null pool uses **all scored candidates** not just retained queue (check `pin_score` OnlineStats in `match_engine.rs` scoring loop).

- **Hook:** `strong_score.rs` `STRONG_CAL_MIN_CANDIDATES` (mirror Tailor `TAILOR_MIN_CANDIDATES` in `psm.rs`)

#### 1c. `DeltaTailor` — calibrated deltaCn analog

```
DeltaTailor = -log10(p1) - (-log10(p2))  using Tailor-scaled scores for rank-1 vs rank-2
```

- **Hook:** `match_engine.rs` where `delta_raw_score` captured (~1070); mirror in `PsmFeatures`

### Phase 2 — Empirical per-spectrum null (patent-free SpecProb substitute)

**Reference:** Keich & Noble JPR 2014 (MC decoys); RS³ Gate 0 in design doc.

#### Module: `crates/scoring/src/scoring/empirical_null.rs` (new)

```rust
// Per spectrum S, after scoring N candidates:
// scores[] = all candidate RawScores (or StrongScores)
// p_emp(T) = #{s >= T} / N   (include target candidates as null proxy — Tailor uses same pool)
// emit: NegLog10EmpiricalP, EmpiricalStdScore
```

**Integration:**
1. During `run_chunk_inner` candidate loop (`match_engine.rs` ~915), accumulate `Vec<f32>` scores per spectrum (already have Tailor histogram — extend).
2. In `fill_post_topn` (~1050), compute `p_emp` for rank-1 and rank-2.
3. Add to `PsmFeatures` + `pin.rs` column list (~930).

**Cost:** O(N) per spectrum — already paying for Tailor histogram.

**Validation:** Gate 0 — correlate `NegLog10EmpiricalP` with brute-force random-peptide MC on 20 UPS1 spectra.

### Phase 3 — Novel: score-landscape DP (patent-distinct from GF)

**Key insight [workspace]:** `g(m)` = node kernel evaluated on mass axis is **spectrum-determined once**, independent of peptide identity.

For peptide of mass M with splits at masses m₁…m_{k-1}:
```
Score = Σ g(m_i) + cleavage_credit
```

**Null:** random cleavage sites → partial sums of AA random walk conditioned on total M.

#### 3a. Fast approximate p-value (no saddlepoint)

Use **empirical null over decoy peptides** but score each decoy via **cached g(m)** table (O(length) per decoy, not full rescore):

- Build `g[m_bin]` once: loop mass bins, call `score_landscape()` using existing node kernel
- For each null peptide: sum `g[m_bin(cleavage)]` — vectorized

**Hook:** new fn adjacent to `score_psm` in `psm_score.rs`; call from `match_engine.rs` post-topn.

#### 3b. Andromeda-style local q on landscape

Within each 100 Th window, only peaks with rank ≤ q contribute to Σ g(m). Optimize q per spectrum by maximizing calibrated sum. **Low-res suited** (Cox 2011 local density).

**Hook:** `scored_spectrum.rs` — optional preprocessing pass before node scoring.

### Phase 4 — Cleavage-credit calibration

RawScore adds cleavage credit separate from node sum. MS-GF+ includes enzyme cleavage in graph; your split may decouple.

- **Audit:** `grep cleavage` in `psm_score.rs` / `match_engine.rs`
- **Idea:** apply separate Tailor denominator to cleavage-only component vs node-only component → two PIN features

### Phase 5 — Own-model rescoring (no external weights)

| Feature | Source | License |
|---------|--------|---------|
| Rank LLR tables | Own training | OK |
| Rich-ion GBDT | `rich_ion_llr()` in `strong_score.rs` | OK |
| MS2PIP-style | Train XGBoost on ProteomeTools CID trap | Apache code |

**Do not** ship Prosit/InstaNovo weights in Apache release.

---

## Original ideas (not in literature — validate by experiment)

### Idea 1: **Calibration stack** (multiplicative in log space)

```
log_sig = log(RawScore) - log(Q99) + α·ChanceMatchSurprise - β·CandidateRankEntropy
```

Percolator learns α,β but starting features are physically motivated. **Hypothesis:** entropy corrects dense-spectrum null inflation where Tailor alone fails (few candidates spread across many bins).

### Idea 2: **Listwise null**

When `top_n > 1`, use retained queue scores as **conditional null** for rank-1 (already partially in `listwise_score_gap`). Extend to:

```
p_listwise = P(s_1 > s_2 + margin | both from same candidate pool)
```

### Idea 3: **Peak-rank shuffle null** (Monte Carlo, spectrum-only)

Shuffle intensity ranks within ±0.5 Da windows, rescore top peptide → empirical null without peptide sampling. Tests whether signal is in **rank pattern** vs **mass accidents**.

### Idea 4: **Integer vs float A/B as ranking** (not just PIN)

Run `ScoreMode` variant ranking by `rank_score_float` on low-res only. Literature (Frank) says ranks matter; engine keeps integer for Java parity. **Low-res-only flag** `--score float-rank` may close gap without touching high-res.

### Idea 5: **σ-adaptive Tailor**

Use Q99 from candidates with `|Δm| < 2·σ` only (soft-match core) vs all candidates — reduces contamination from accidental wide matches in dense spectra.

### Idea 6: **Percolator feature: `SpecDensity`**

`peaks_per_da * candidate_count` — encodes heterogeneity MS-GF+ SpecProb corrects analytically; let GBDT learn correction if empirical null too expensive.

---

## Tests to add

| Test | Location | Assert |
|------|----------|--------|
| `empirical_null_uniform` | `empirical_null.rs` | p=0.5 for median of uniform scores |
| `tailor_monotonic` | `psm.rs` existing | extend for `DeltaTailor` |
| `landscape_cache` | `psm_score.rs` | g(m) independent of peptide sequence |
| `gate0_vs_bruteforce` | `tests/` integration | \|log10 p_emp - log10 p_MC\| < 0.2 |

---

## Files to touch (expected diff)

```
crates/scoring/src/scoring/empirical_null.rs   [NEW]
crates/scoring/src/scoring/mod.rs              [export]
crates/search/src/psm.rs                       [PsmFeatures fields]
crates/search/src/match_engine.rs              [accumulate + compute]
crates/output/src/pin.rs                       [columns]
crates/output/src/qpx.rs                       [mirror columns]
docs/superpowers/experiment-protocol.md        [UPS1 low-res gate]
```

---

## Do NOT implement (FTO / license)

- `gf/` generating function DP (`docs/.../phase6-generating-function-plan.md` — superseded)
- MSFragger hyperscore tail copy for commercial use
- InstaNovo weights
- Paragon feature probabilities

---

## Benchmark commands (fill from experiment-protocol)

```bash
# Low-res UPS1 entrapment gate (paired r=1)
# astral-speed/rust or msgf-rust — use project's prov.sh / benchmark script
# Compare: baseline vs +NegLog10EmpiricalP vs +DeltaTailor
# Report: PSMs @ 1% paired entrapment-FDP
```

---

## Cross-reference: RS³ campaign ruling

The analytic saddlepoint RS³ path was **superseded** by decoy-calibrated empirical null — see [`2026-06-29-rs3-spectral-significance-design.md`](2026-06-29-rs3-spectral-significance-design.md) and campaign plan. Implement **empirical null first**; saddlepoint only if Gate 0 passes against brute-force MC.
