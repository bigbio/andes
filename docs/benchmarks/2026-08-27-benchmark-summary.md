# Benchmark summary — August 2026

Where andes stands on public data, measured against established engines. Comparison
engines are anonymised as **Software 1/2/3**; the datasets, FASTAs and settings are all
named, so every figure here is reproducible.

Two rules govern everything below:

1. **Yield is never reported alone.** Every glyco figure is paired with an entrapment
   false-discovery proportion (FDP), so "more IDs" cannot hide "more wrong IDs".
2. **A single Percolator replicate is not a measurement.** See
   [Measurement methodology](#measurement-methodology) — this bit us repeatedly and is
   the most transferable thing in this document.

---

## 1. Peptide identification (closed search)

| dataset | andes PSMs / peptides | Software 3 PSMs / peptides | andes delta | andes time / Software 3 time |
|---|---|---|---|---|
| Astral, high-res LFQ (PXD070049) | **38,437 / 24,981** | 31,435 / 20,608 | **+22.3% / +21.2%** | 899 s / **217 s** |
| TMT, low-res ion-trap CID (PXD007683, run `a05058`) | **12,316 / 10,977** | 10,504 / 9,008 | **+17.2% / +21.9%** | 193 s / **80 s** |

**andes identifies 17–22% more, and takes 2–4x longer.** Both halves are the result.

Harness validation: Software 3 reproduced its own published Astral figure exactly
(31,435 PSMs; 20,608 peptides against a published 20,607), and andes reproduced its
published TMT peptide count (10,977 against 10,691 on the depositors' database). A harness
that cannot reproduce known numbers cannot be trusted to produce new ones.

Caveat: the TMT protein database is a 2026 reconstruction rather than the depositors'
2018 file, so absolute TMT counts are comparable *between engines in this table* but not
against published 2018 values.

---

## 2. Intact N-glycopeptide identification

Human plasma (PXD030622), entrapment-controlled. **This is where andes loses**, and the
gap is architectural rather than a tuning deficit.

### 2a. Against Software 1 — like-for-like at matched 1% FDR

| metric | andes | Software 1 |
|---|---|---|
| glycoPSMs @1% FDR | **254** | **587** |
| unique glyco peptides | **64** | **91** |
| relative speed | ~16x slower | — |

**andes recovers 43% of Software 1's glycoPSMs and 70% of its glyco peptides.**

A counting note that matters: Software 1's raw accepted count is 3,126, but 2,539 of
those are *unmodified* peptides. Only 587 carry a glycan. Comparing 3,126 against andes's
glyco-only 254 would overstate the gap by roughly 5x. Glyco comparisons must filter on
delta mass (> 800 Da) before any ratio is computed.

### 2b. Against Software 2 — the depositors' own published results

Different run and different definitions from 2a, so these are **not** interchangeable
with the numbers above.

| metric | Software 2 | andes | andes / Software 2 |
|---|---|---|---|
| glycoPSMs @1% | 539 | 247 | 46% |
| unique glycopeptides | 185 | 152 | 82% |
| bare peptide sequences | 78 | 58 | 74% |
| peptide-sequence overlap | — | 48 / 78 | **62% recovered** |

The shape of the deficit is informative: andes recovers a **minority of the PSMs but a
majority of the peptides**. It is finding much of the same biology with less depth per
glycopeptide, rather than missing the analytes outright.

Corroborating that: **74% of andes's glycoPSMs carry a glycan mass Software 2 also
observed**, and the median glycan mass agrees (2,205 Da in both). Glycan *mass*
determination is largely correct; glycan *composition* assignment was the weaker link
(see below).

---

## 3. What moved the glyco numbers, and what did not

Roughly a dozen levers were measured this campaign. **Almost all were null or negative.**
Recording that is the point — it is why the remaining gap is described as architectural.

**Negative — generation-side expansion.** Three independent attempts to enlarge the
candidate pool (a wider glycan mass window, two-axis Y-ion retention, an isobar resolver)
each moved yield *down*. A larger surviving pool gives decoys more places to fit, which
tightens the accepted threshold. Generation is not the bottleneck.

**Negative — over-gating.** Combining gates (core-Y, minimum matched ions, selector
reweighting) lost on both axes at once: fewer PSMs *and* roughly 6x higher entrapment FDP.
Tightening acceptance did not trade yield for calibration; it forfeited both.

**Not demonstrable — selector reweighting.** Individual weight changes produced apparent
double-digit gains in single runs that did not survive replication (Section 4).

**Positive, on non-yield evidence — species-appropriate glycan lists.** The default
composition list included N-glycolylneuraminic acid (NeuGc). Humans lack a functional
*CMAH* gene and do not synthesise it, and NeuGc is exactly degenerate in precursor mass
with a NeuAc/Hex/Fuc combination — so including it injects isobaric decoy-like
compositions that human data can never populate. Excluding it on human samples raised
agreement with Software 2's independently-assigned glycan masses from **74% to 91%**.

That agreement figure, not a yield delta, is the evidence. The corresponding yield change
(+11.6%, five seeds) sits inside measurement noise. The change is retained because the
mechanism is externally verifiable — every mainstream human-labelled glycan list also
ships zero NeuGc — not because it produced more IDs.

---

## 4. Measurement methodology

The most reusable result here is about *measuring*, not about andes.

### Target-decoy q-values are a step function at low counts

Percolator's minimum achievable q-value is `1 / T_top`, so identifications at a 1% FDR
threshold move in discrete jumps. At glycopeptide counts (a few hundred targets) this
makes the reported yield strongly dependent on the random seed. **Fractions must be
pooled before scoring, and even pooled results must be replicated across seeds.**

### The detection floor

Re-scoring identical pooled inputs under five Percolator seeds:

- within-arm standard deviation ≈ **66 PSMs**
- within-arm range routinely spans **130–230 PSMs**

Required replicates per arm for 80% power at α = 0.05:

| true effect | seeds per arm |
|---|---|
| 23 PSMs | ~129 |
| 50 PSMs | ~27 |
| 100 PSMs | ~7 |
| 150 PSMs | ~3 |

**With five seeds per arm the smallest detectable effect is ~117 PSMs — about a 58%
relative change.** Everything smaller is unresolvable at any sane compute budget.

### Consequences

- Several single-replicate results in the +25% to +50% range were **retracted** once
  replicated. They were draws from a distribution, not effects.
- "Not demonstrable at this power" is the correct verdict for most of them — *not*
  "refuted". This design cannot distinguish a genuine +10–20% improvement from zero.
  Absence of evidence is not evidence of absence.
- For effects below ~50%, use an instrument without the step function: agreement with an
  orthogonal engine's assignments, entrapment FDP at matched yield, or fixed-score
  threshold counts. The 74% → 91% agreement figure in Section 3 is an example.
- An entrapment FDP of 0.00% is **not** a success signal. At these counts it usually
  means the arm is over-conservative, and it is statistically indistinguishable from a
  baseline reporting 0.5%.

---

## 5. Summary

| area | standing |
|---|---|
| Peptide IDs, high-res and low-res | **Ahead**, +17–22% |
| Speed, peptide search | **Behind**, 2–4x |
| Intact N-glycopeptide IDs | **Behind**, 43–46% of established tools |
| Speed, glyco search | **Behind**, ~16x |
| Glycan mass determination | Largely correct (74% external agreement; 91% with a species-appropriate list) |
| Glycan composition assignment | Was the weak link; improved by species-aware lists |

The glyco gap is not a tuning problem. Generation is adequate — the candidates are
largely present — and gate and weight adjustments have been measured to exhaustion
without closing it. The remaining deficit is in **separation**: the ranking function's
most heavily weighted terms vary little between the peptides competing for a single
spectrum, so they cannot discriminate among them. Closing the gap requires restructuring
that ranking, not reweighting it.
