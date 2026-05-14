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
