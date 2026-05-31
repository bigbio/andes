# TMT −5.5% gap: diagnosed — REAL, but NOT a gross GF-DP bug (2026-05-31)

After the chimeric two-pass cascade was optimized and gate-validated on all 3 datasets,
TMT remained the lone merge-gate blocker (Rust 9,628 vs Java 10,194 @1% FDR, −5.5%).
The standing hypothesis (memory: "Lever-2a") was that the gap is a divergence in the
GF SpecEValue DP itself. This note **refutes the gross-GF-bug framing** and localizes
the gap precisely, using the cascade binary (`feat/chimeric-dda-plus` HEAD `3d941a02`).

Config: `tmt-data/a05058.mzML`, `PXD007683_..._reviewed.fasta`, CID / high-res / TMT,
`-t 20ppm -ti -1,2`, Java `-precursorCal off`, Rust `--precursor-cal auto --chimeric`.
Same-machine A/B + entrapment (FDRBench shuffled-target 1:1) on both engines.

## 1. The gap is REAL (entrapment-confirmed, not a coincidental-target artifact)

| TMT entrapment ruler | Rust | Java |
|---|---:|---:|
| real PSMs (entrapment-validated @1%) | 8,436 | 9,224 |
| entrapment FDP (combined) | 0.80% | 1.08% |

Both the reversed-decoy ruler (9,628 vs 10,194) and the entrapment ruler (8,436 vs
9,224 real) agree: Java genuinely identifies ~788 more real TMT PSMs. Java's FDP is
slightly looser (1.08% vs 0.80%) but its real count still leads. Unlike the chimeric
"reversal," this gap survives entrapment scrutiny.

## 2. The GF SpecEValue is in GLOBAL PARITY (no gross bug)

Over all **15,595** scans where Rust and Java pick the **same** top peptide:

| metric (median) | Rust | Java | Δ (Rust−Java) |
|---|---:|---:|---:|
| RawScore | 58 | 49 | +9 |
| DeNovoScore | 79 | 73 | +6 |
| **lnSpecEValue** | **−20.46** | **−20.39** | **+0.05** |
| RawScore − DeNovoScore (headroom) | −20 | −21 | +1 |

RawScore/DeNovoScore are scale-shifted (+9/+6) but that is cosmetic: the FDR-relevant
**lnSpecEValue is effectively identical** (Δ +0.05), and the score headroom matches.
The GF DP is fundamentally correct on TMT. **There is no gross GF-shape bug to fix.**

## 3. Where the 1,238 Java-wins / Rust-misses scans actually go

Of scans Java IDs @1% FDR that Rust misses (Java\Rust = 1,238; Rust\Java = 611;
net Java advantage 627 ≈ the −566..−788 gap):

| bucket | count | share | nature |
|---|---:|---:|---|
| SAME peptide, Rust scored it but it failed FDR | 800 | 65% | boundary SpecEValue: Rust lnSpecEValue +0.54 worse (median), enough to flip pass→fail under Percolator |
| LABEL FLIP to a different target | 256 | 21% | Rust ranks a different target peptide top-1 |
| LABEL FLIP to a decoy | 182 | 15% | Rust ranks a **decoy** top-1 (no target row) |
| (enumeration miss) | 0 | 0% | Rust scores every scan — no missing peptides |

So the entire gap is **scoring/ranking quality**, distributed across (a) boundary
SpecEValue cases at the FDR margin and (b) top-1 ranking divergences (438 flips total).
No additive enumeration fix exists.

## 4. Implication for the fix

This is the hard, Percolator-sensitive territory the n=8 audit (`MEMORY.md`) repeatedly
shows per-feature Java-parity fixes REGRESS. The only scoring-fix class that has ever
gained PSMs is **top-1-changing fixes that RESTORE Java's ranking** (iter29, iter33).
The 438 label-flips are the candidate lever: find why Rust ranks a different peptide
(or a decoy) top-1 on these TMT scans, and whether a ranking change there restores
Java's choice without perturbing the emitted distribution (the T/D-ratio canary).
This requires per-candidate score tracing (Rust #1 vs Java #1 on the same scan) —
`msgf-trace`-style instrumentation, multi-day, Rule-2 regression risk.

The boundary-SpecEValue 800 (Δ +0.54) are NOT separately fixable without moving the
global SpecEValue distribution (which is already in parity) — touching it risks the
15,595 same-peptide scans that currently agree.

## Status

TMT remains the merge-gate blocker. The cascade cannot help (CID narrow isolation →
Pass-2 = 2.66s, ~no co-isolation). The fix is a per-candidate ranking investigation on
the 438 label-flips, not a GF-DP rewrite. Astral (+55%) and PXD001819 (+21.4%) are
fully won on both axes (entrapment-clean); see `MEMORY.md` index.

Analysis scripts (local on VM `/srv/data/msgf-bench/`): `tmt_pin_diff.py`,
`tmt_pin_diff2.py`, `tmt_norow.py`, `java_tmt_entrap.sh`.

---

## Addendum: traced the 438 label-flips to the root (2026-05-31, user chose "trace flips")

Ran Rust narrow `--top-n 10` on TMT (cal-auto), joined Java's @1% peptide against
Rust's full ranking on the 438 flip scans (233 to-target, 205 to-decoy):

- **Java's peptide is in Rust's top-10 (by RawScore) for only 20/438 (5%).** For 418
  (95%) it is **not in Rust's top-10 at all** — Rust scores it well below contention.

Ruled out, one by one:
- **Mass window / calibration** — NOT it. Cal tightened TMT 20ppm→8.135ppm (robust
  sigma 2.545ppm), but `--precursor-cal off` (full 20ppm) gave FEWER PSMs (9,238 vs
  cal-auto 9,628), not more. Widening admits more noise than it recovers; tightening
  is net-helpful. The window is not excluding Java's peptides on net.
- **Modifications** — NOT it. `mods.txt` (Java) and `mods-numeric.txt` (Rust) are
  numerically identical: Carbamidomethyl-C 57.021464, TMT6plex-K & N-term 229.162932,
  Ox-M 15.994915, NumMods=1. Same masses → same enumerated peptide space.
- **Top-1 selection criterion** — NOT it. `PsmMatch::Ord` already ranks by
  (spec_e_value asc, rank_score desc). The issue isn't SpecE-vs-RawScore selection:
  Java's peptide isn't even in Rust's RawScore top-10, so a larger SpecE-eval pool
  (K>1) wouldn't include it.
- **Aggregate T/D discrimination** — NOT worse. Top-1 RawScore target−decoy
  separation: Rust 14 vs Java 13; lnSpecE separation Rust 2.71 vs Java 2.17; both 24%
  competitive decoys; T/D 2.16 vs 2.14. Rust's aggregate discrimination is equal-or-
  better.

**Root:** Rust computes a LOWER RawScore (node score) than Java for Java's winning
peptide on these specific CID/TMT spectra — the peptide falls out of Rust's
RawScore-ranked contention entirely. This is a **per-peptide CID node-scoring
divergence** (RankScorer ion-match / per-rank log-prob behavior on CID spectra), the
deepest scoring layer. Both engines load the same `CID_HighRes_Tryp.param`, so the
divergence is in how Rust APPLIES it (ion-type list, peak-rank assignment, or per-rank
probability tables) on CID fragmentation — the same three hypotheses as the 2026-05-20
Astral score_psm divergence, now CID-specific.

**Why no clean fix is visible:** aggregate Rust discrimination is already competitive
(equal/better T/D separation), the GF SpecEValue is in global parity, and there is no
window/mods/selection bug. The gap is emergent from per-peptide CID node-score
differences in the FDR tail — exactly the territory the n=8 audit shows resists
piecewise fixes. Closing it requires per-ion CID score tracing (dump Rust's matched
ions + per-rank scores for a flip peptide, diff against Java instrumentation) —
multi-day, Rule-2 regression risk, low probability of a clean Percolator-FDR win.

**Recommendation:** bank the cascade wins (Astral +55%, PXD +21.4%, entrapment-clean,
faster-or-tied) as an opt-in `--chimeric` PR; treat TMT as a documented deep-CID-
scoring limitation. The per-ion CID trace is a separate research project if pursued.
