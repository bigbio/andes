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
