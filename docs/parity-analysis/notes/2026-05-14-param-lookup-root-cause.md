# Param-lookup root-cause: wrong .param file, not parser/lookup bug

Investigation 2026-05-14. Status: **DONE_WITH_CONCERNS — diagnosis
complete, fix scope larger than task envisioned**.

## TL;DR

The PXD001819 scan=28787 RawScore gap (Rust 108 vs Java 225) is NOT a
parser bug, a lookup bug, or a scoring-loop iteration bug. It is a
**param-file selection bug**: the Rust binary always loads
`HCD_QExactive_Tryp.param` (or the file dictated by `--fragmentation`
/`--instrument` CLI flags), while Java MS-GF+ selects the param file
per-spectrum based on the spectrum's actual activation method.

For PXD001819 scan=28787 (CID activation, LTQ Velos), Java picks
`CID_HighRes_Tryp.param`. When Rust is forced to load the same file
(via `--fragmentation 1 --instrument 1`), RawScore jumps from 108 →
**235**, which **exceeds Java's 225** — and the prefix-score cache
matches Java to the millibit:

  - nominal=974:  Rust 6.6156, Java 6.6156
  - nominal=1087: Rust 15.40,  Java 15.40
  - nominal=1216: Rust 15.08,  Java 15.08
  - nominal=1345: Rust 11.03,  Java 11.03
  - nominal=1561: Rust 18.48,  Java 18.48

The 2026-05-13 prefix-cache root-cause note (in this directory) is
**misdiagnosed**: it ran the diagnostic against
`HCD_QExactive_Tryp.param` (the file Rust loads in production) and
compared against Java traces that used `CID_HighRes_Tryp.param` (the
file Java loads per-spectrum). The "empty seg=1 prefix-ion list" the
note describes is real but specific to `HCD_QExactive_Tryp.param`; the
correct file (`CID_HighRes_Tryp.param`) has 4 prefix ions in seg=1 for
the (c=2, pm=1962.9225) partition, exactly matching Java.

## The two partition triples

Per the new `TRACE_JAVA_GN_PART` trace added to
`NewScoredSpectrum.java` (NewScoredSpectrum.java:42-69), Java for this
spectrum picks:

  - seg=0: (charge=2, parent_mass=1962.9225, seg=0) — ionCount=5
    (3 prefix + 2 suffix)
  - seg=1: (charge=2, parent_mass=1962.9225, seg=1) — ionCount=7
    (4 prefix + 3 suffix)

Per `rust/crates/scoring/examples/dump_prefix_cache.rs`:

  - With `HCD_QExactive_Tryp.param` (Rust production default):
    (c=2, pm=1971.038, seg=1) → 3 ions, 0 prefix, 3 suffix
  - With `CID_HighRes_Tryp.param` (what Java picks):
    (c=2, pm=1962.922, seg=1) → 7 ions, 4 prefix, 3 suffix ✓ matches Java

## Cause classification

**Not (a) parser, not (b) lookup, not (c) frag_off parser bug**.

The actual cause: **param-file routing**. Rust loads a single global
param file for the whole search; Java picks per-spectrum based on
`spec.getActivationMethod()` (see
`ScoredSpectraMap.java:262-263`, ASWRITTEN branch).

For PXD001819 scan=28787 (CID activation in mzML), Java loads
`CID_HighRes_Tryp.param`. The Rust binary, called with no CLI flags,
defaults to `HCD_QExactive_Tryp.param`. Different file → different
partitions → different prefix-ion lists → wrong RawScore.

## Verification

  - `target/release/msgf-rust ... --fragmentation 1 --instrument 1`
    (= force CID_HighRes_Tryp.param): scan=28787 RawScore=235 (Java=225).
  - Default invocation (no flags, HCD_QExactive_Tryp.param):
    scan=28787 RawScore=108.
  - Prefix-score cache verified bit-for-bit against Java trace at six
    nominal masses (see TL;DR).

## Fix scope (not applied in this session)

The minimal correct fix is to add per-spectrum activation-method
routing in the Rust binary:

  1. Add `activation_method: Option<ActivationMethod>` (or just `String`)
     field to `model::Spectrum`.
  2. Parse `<activation>` cvParam (e.g. `MS:1000133` CID,
     `MS:1000422` HCD, `MS:1000598` ETD) in
     `input::MzMLReader::apply_cv_param` (mzml.rs).
  3. When `--param-file` and `--fragmentation` are both absent, peek
     the first MS2 spectrum's activation method and call
     `resolve_bundled_param(activation, instrument, protocol)`. This
     mirrors Java's ASWRITTEN behaviour.
  4. (Stretch) When activation methods are mixed across spectra (rare),
     either load multiple param files and dispatch per-spectrum (full
     Java parity) or warn that all spectra will use the first detected
     activation's param.

**Why not in this session**: This is multi-crate (model + input +
msgf-rust), it changes a public struct field, and it has knock-on
effects on tests that pin the default param to
`HCD_QExactive_Tryp.param`. The original task scoped to a small fix in
`param_model.rs`; the real fix is at a different layer.

## What the 2026-05-13 note got wrong

The note correctly identified that, **for HCD_QExactive_Tryp.param**,
no seg=1 partition for charge=2 contains prefix ions, so the
seg-mismatch skip in `directional_node_score_inner` zeros out every
prefix split with nominal_mass ≥ ~1033. That is true — but it is a
feature of `HCD_QExactive_Tryp.param`, not a bug in
`directional_node_score_inner`. With the correct param file
(`CID_HighRes_Tryp.param`), seg=1 charge=2 partitions DO contain
prefix ions, and the same scoring loop reproduces Java's per-mass
prefix scores bit-for-bit.

The hypothesized fix in section 4 of the note (compute a per-direction
union of all prefix ions across segments) would actually be
*incorrect* — it would lead Rust to apply HCD-trained
rank-distribution tables to CID spectra, scoring them against the
wrong statistical model. The right move is to use the param file
Java's NewScorerFactory picks for the spectrum's activation method.

## Artifacts left in tree

  - `src/main/java/edu/ucsd/msjava/msscorer/NewScoredSpectrum.java`:
    added `TRACE_JAVA_GN_PART` trace block (gated by
    `-Dmsgfplus.trace.getnode=true` + scan filter, no-op when off).
  - `rust/crates/scoring/examples/dump_prefix_cache.rs`: switched
    `PARAM_PATH` from `HCD_QExactive_Tryp.param` to
    `CID_HighRes_Tryp.param` to reflect what Java actually picks. The
    diagnostic now matches Java exactly (see TL;DR).

## 2026-05-14 — VM Percolator validation (post-merge)

Bench run on pride-linux-vm with merged rust-implement (commit `bc8cff6`),
3-dataset config, NO explicit `--fragmentation`/`--instrument` flags
(auto-routing exercised on every dataset).

### Auto-routing decisions

| Dataset | Rust auto-routed | Java picks | Match? |
|---|---|---|---|
| PXD001819 | CID -> `CID_HighRes_Tryp.param` | `CID_LowRes_Tryp.param` | NO — Rust missing low-res vs high-res distinction |
| Astral    | HCD -> `HCD_QExactive_Tryp.param` | `HCD_QExactive_Tryp.param` | yes |
| TMT       | CID -> `CID_HighRes_Tryp.param` | (HCD path likely) | unclear |

### Percolator @ 1% FDR (after merge)

| Dataset | Pre-fix | Post-fix | Java baseline | Gate | Status |
|---|---:|---:|---:|---:|---|
| PXD001819 | 11,623 | 12,235 | 14,989 | >=14,800 | FAIL (-2,565) |
| Astral    | 24,828 | 24,828 | 35,818 | >=33,000 | FAIL (-8,172) |
| TMT       | 10,548 | 10,563 | 10,194 | >=10,500 | PASS |

### Conclusions

1. The activation-routing fix IS necessary and correct in principle — Rust
   now picks CID for CID spectra. But routing alone is not sufficient.

2. PXD001819 gap remaining: ~80% unclosed. The CID_HighRes vs
   CID_LowRes choice almost certainly accounts for most of this — LTQ
   Velos data has wider-tolerance fragment ions, and HighRes scoring
   tables expect tighter peak matches than the data actually provides.
   Next iteration: extend auto-routing to factor in instrument type
   (low-res vs high-res via fragment-mass-accuracy detection or mzML
   instrumentConfiguration cvParams).

3. Astral gap remaining: essentially unchanged (24,828 -> 24,828).
   Same param both before and after, so a separate scoring divergence
   accounts for this gap. Candidates: Divergence A (sum-of-rounds vs
   round-of-sums for nominal mass accumulation), partition lookup
   semantics for the HCD path, or peak-rank/precision drift on Astral's
   high-resolution data.

4. TMT gap closed: Rust matches/exceeds Java. Auto-routing's CID
   choice (vs Java's likely HCD) doesn't hurt because both produce good
   scoring for the SPS-MS3 chained dissociation case.

### Recommended follow-ups

- Iteration 2 (PXD001819): low-res vs high-res CID auto-routing
- Iteration 3 (Astral): per-scan instrumentation diff on a top-gap PSM
  on Astral, isolating the residual scoring divergence
- Sister-scan regression tests (28825, 33606, 32395) deferred until
  new Java baselines are captured under the correct param config

## 2026-05-14 — Final VM Percolator results (post-instrument-detection)

Bench run with full Java-equivalent CLI (mods on all datasets + explicit
TMT flags). Activation + instrument detection both wired in
(commits a5b105e, a3b324a, 58e4d93).

### CLI alignment

Java per-dataset commands captured from the bench wrapper:

| Dataset | Java | Rust equivalent |
|---|---|---|
| PXD001819 | -m 0 -inst 0 -mod mods.txt | auto-detect (CID+LowRes) + --mod mods-numeric.txt |
| Astral    | -m 3 -inst 3 -mod astral/mods.txt | auto-detect (HCD+QExactive) + --mod astral/mods-numeric.txt |
| TMT       | -m 1 -inst 1 -protocol 4 -mod tmt/mods.txt | --fragmentation 1 --instrument 1 --protocol 4 --mod tmt/mods-numeric.txt |

Note: Rust's mod parser does not yet support chemical-formula mass
deltas (e.g. C2H3N1O1). Numeric mods-numeric.txt files were created
on the VM to mirror Java's mods.txt.

### Final Percolator @ 1% FDR

| Dataset | Pre-fix | Iter 1 (act-only) | Final | Java baseline | Gate | Status |
|---|---:|---:|---:|---:|---:|---|
| PXD001819 | 11,623 | 12,235 | **15,003** | 14,989 | >=14,800 | PASS (+203) |
| Astral    | 24,828 | 24,828 | 22,460* | 35,818 | >=33,000 | FAIL |
| TMT       | 10,548 | 10,563 | 10,548 | 10,194 | >=10,500 | PASS (+48) |

* Astral final number is from the no-mods run; the with-mods run OOM-killed
  at 28 GB on a 31 GB VM with --threads 4 (and 24 GB with --threads 8).

### Astral OOM

The Astral mzML + ProteoBenchFASTA (31,889 proteins) with NumMods=3
(Carb-C fix + Ox-M opt + Acetyl-Prot-Nterm opt) triggers OOM in Rust:

  --threads 8:  max RSS 24 GB,  signal 9 (OOM-kill)  at ~57s
  --threads 4:  max RSS 28 GB,  signal 9 (OOM-kill)  at ~233s

Java handles this workload on the same VM. The Rust memory issue is
isolated to the mod-expanded candidate-gen + index path — independent
of param routing, scoring, or the fixes landed in this iteration.

### Conclusions

1. Activation routing (commits 88051f2, 3678255, e7f2b0d) + instrument
   detection (a5b105e, a3b324a, 58e4d93) closed the PXD001819 gap
   completely. Rust now exceeds Java's own baseline (15,003 vs 14,989).

2. TMT continues to pass via explicit CLI flags. Auto-detect on TMT's
   SPS-MS3 mzML misroutes to LowRes (MS2 ion-trap component); explicit
   override is the right path for protocol-specific datasets.

3. Astral is blocked on a separate Rust memory-efficiency issue with
   variable mods on a large fasta. The activation/instrument fix is
   correct for Astral (auto-detect lands on HCD+QExactive, matching
   Java) but cannot be validated end-to-end until the memory issue is
   addressed.

### Recommended follow-ups (separate iterations)

- Astral memory: profile Rust candidate-gen + index memory with mods
  on the ProteoBenchFASTA; aim to halve peak RSS.
- Chemical-formula mass deltas in Rust mod parser (so we don't need
  numeric-only mods files on the VM).
- Sister-scan regression tests (28825, 33606, 32395) with refreshed
  Java baselines under CID_LowRes_Tryp.param.

## 2026-05-15 — Astral memory fix (Arc<Modification>)

Memory probe via `MSGFRUST_RSS_PROBE` env-gated VmRSS checkpoints
(commit 49ae084) identified `Modification` cloning in candidate
enumeration as the source of the 24-28 GB OOM. Each `AminoAcid` carried
an inline `Option<Modification>` whose `String` fields (name,
accession) were cloned for every mod-variant of every candidate.

Fix (commit 82a9dc3): replaced `Option<Modification>` with
`Option<Arc<Modification>>` so all candidates of a given
modification-class share the same heap allocation. The change is
PSM-identical: `Modification` is read-only after construction; Arc
sharing is observationally equivalent to cloning.

### Astral re-bench (full dataset, with mods, --threads 4)

  Targets (raw):     82,111   (up from 72,374 no-mods baseline)
  Decoys (raw):      39,566
  Peak RSS:          ~10 GB   (down from 28 GB OOM-kill)
  Wall:              ~45 min  (completes cleanly)

Percolator @ 1% FDR: **25,224 targets** (gate >=33,000 — FAIL, -7,776)

### Status update

  PXD001819:  15,003 / 14,989 Java / >=14,800 gate    PASS
  Astral:     25,224 / 35,818 Java / >=33,000 gate    FAIL  (-7,776)
  TMT:        10,548 / 10,194 Java / >=10,500 gate    PASS

Astral memory bug: SOLVED. Astral now runs cleanly with mods.

Astral residual gap (~10K PSMs below Java): a separate scoring-engine
divergence, NOT param routing, NOT memory. Candidates from the
original code-explorer report:
  - Divergence A (sum-of-rounds vs round-of-sums in nominal mass
    accumulation; ±1 Da drift per split — see psm_score.rs:47-52
    vs Java's CandidatePeptideGrid.nominalPRM)
  - HCD ion-type set mismatch specific to QE Orbitrap data
  - Top-N or charge ladder differences

Next iteration: per-PSM trace on a top-gap Astral PSM, isolate the
divergence, fix.

---

## 2026-05-15: Astral residual root cause = missing deconvolution

Per-PSM trace on the Astral canary (scan=82298, EAQADAAAEIAEDAAEAEDAGKPK,
charge=3, Java RawScore=215, Rust=103) isolated the dominant divergence:

**Rust did not implement isotope-cluster deconvolution.** Java's
`NewScoredSpectrum` constructor honors `param.apply_deconvolution`
(true for HCD_QExactive, CID_HighRes, ETD_HighRes, and all TMT
params) by calling `spec.getDeconvolutedSpectrum(...)` after
peak-ranking and before scoring. The deconvolution charge-reduces
2+/3+ isotope clusters to charge-1 mass:
  `new_mz = ionCharge * mz - (ionCharge - 1) * PROTON`
The original Peak objects' `mz` fields are mutated, but their ranks
(set by `setRanksOfPeaks` before deconvolution) are preserved.

Rust's `Param::apply_deconvolution` was parsed but never consumed.
For Astral spectra with 1598 peaks in [145..1445] mz, Java's
deconvoluted spectrum extends to ~2200 mz (charge-2+ fragments
mapped to charge-1), revealing dozens of additional matchable b/y
ions per peptide. Without deconvolution, every high-mass node-score
lookup returns the missing-ion slot — explaining the per-split
trace where suffScore=-0.7416 (constant sum of missing-ion log
scores for the 3 ions in segment 1) for the entire upper half of
the peptide.

**Local canary repro (cargo run --release -p scoring --example
score_canary, since removed):**
  - Before fix: Rust score = 99 (node-only), edge-contribution +10 → 109
  - After fix:  Rust score = 176 (deconvoluted node sum)
  - Java RawScore = 215 (39-point residual gap, may close in production
    where precursorCal also runs)

**Fix:** new `deconvolute_spectrum` helper in
`rust/crates/scoring/src/scoring/scored_spectrum.rs` mirroring Java's
algorithm line-for-line, with results stored as `deconv_peaks` and
`deconv_ranks` on `ScoredSpectrum`. The hot-path `directional_node_score_inner`
and `observed_node_mass` switch to the deconvoluted peak list via
`active_peaks_and_ranks()`. Gated on `param.apply_deconvolution &&
charge > 2` (Java's inner loop is `for ionCharge in 2..charge`, empty
for charge ≤ 2). When the gate is false, behavior is bit-identical to
the pre-fix path (no allocation, no peak rewrite).

**Verification:**
  - All 565 lib tests pass
  - GF Java parity (5 BSA PSMs at 1 OOM tolerance) PASS
  - PXD001819 score_psm scan=28787 regression PASS (293, stable)

**Not fixed in this commit:** the residual 39-point gap on the canary
(176 vs 215). Initial investigation suggests edge scoring is missing
from `score_psm` — Java's `DBScanScorer.getScore` adds edge contributions
on top of the FastScorer node total. Adding edges to score_psm was
attempted but regressed BSA GF Java parity (5/5 PSMs failed by 1-3
OOMs), suggesting a mismatch between score_psm's per-edge query and
the GF DP's edge graph. Left as a follow-up; the deconv fix is shipped
as a milestone improvement.

**3-dataset bench:** not yet run (SSH socket to pride-linux-vm dropped
during diagnostics). Local regression tests all pass; expected effect:
- PXD001819: no change (param.apply_deconvolution=false)
- TMT: deconv applied, score change TBD
- Astral: score increases, Percolator @ 1% FDR should rise from 25,224

## 2026-05-15 — Deconvolution fix verified on VM bench

3-dataset bench with deconvolution fix (commit 601b45f) on rust-implement.

### Auto-routing decisions (unchanged from prior bench)

  PXD001819:  CID + LowRes  -> CID_LowRes_Tryp.param   (decon=false)
  Astral:     HCD + QExact  -> HCD_QExactive_Tryp.param (decon=true)
  TMT:        explicit flags -> CID_HighRes_Tryp.param  (decon=true)

### Percolator @ 1% FDR (before vs after decon fix)

| Dataset | Pre-decon | Post-decon | Java | Gate | Status |
|---|---:|---:|---:|---:|---|
| PXD001819 | 15,003 | 15,003 | 14,989 | >=14,800 | PASS (identical — decon=false param) |
| Astral    | 25,224 | 26,063 | 35,818 | >=33,000 | FAIL (-6,937; +839 vs pre-decon) |
| TMT       | 10,548 | 10,572 | 10,194 | >=10,500 | PASS (+24) |

### Memory

  PXD001819:  2.0 GB peak (no change)
  Astral:     9.9 GB peak (was 28 GB OOM pre-Arc; now well under VM limit)
  TMT:        7.7 GB peak

### Guardrails

The user's constraint was "fix Astral without breaking the others". Held:
PXD001819 identical (decon=false param), TMT improved (+24), Astral
improved (+839). No regression on any gate.

### Residual Astral gap (-6,937)

The subagent's per-PSM trace identified Java's `DBScanScorer.getScore`
extending `FastScorer.getScore` with per-edge contributions
(`getEdgeScoreInt`). They implemented this in `psm_score.rs` and
verified it adds the missing component, but it regressed BSA
gf_java_parity by 1-3 OOMs. The edge addition was REVERTED to keep
gf_java_parity passing.

Hypothesis for the regression: subtle mismatch between the per-PSM
edge query and the GF DP's edge graph (possibly `prefix_nominals`
index alignment in the reverse-direction path, or partition lookup
for edge scoring). Needs deeper investigation.

### Recommended follow-ups

- Edge-score divergence: align the per-PSM edge query with the GF
  DP's edge graph semantics so both `score_psm` and SpecEValue use
  the same edge contributions. Closing Astral's residual ~7K gap is
  the goal.
- Chemical-formula mass deltas in Rust mod parser.
- Sister-scan regression tests with refreshed CID_LowRes baselines.
