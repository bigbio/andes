# iter24: pass Astral mods.txt to Rust — +384 PSMs @ 1% FDR

_2026-05-21. The Astral bench harness was running Java with `mods.txt` (Carbamidomethyl-C + Oxidation-M + **Acetyl-Prot-N-term**) but Rust with no `--mod` flag (defaults only: Cam-C + M-ox, no Acetyl). Adding the acetyl mod to Rust closes another chunk of the gap._

## The finding

While auditing label-flip scans on iter22b-vs-Java, every sampled flip showed Java picking peptides with `+42.011` N-terminal acetylation while Rust picked unmodified alternatives. Counted across 16,956 java_target_rust_decoy flips: **1,427 (8.4%)** had acetylated peptides in Java's top-1.

Root cause: the Astral bench harness (`benchmark/run_astral_3arm.sh`) configures Java with `MODS="$DATA_DIR/mods.txt"`:
```
NumMods=3
C2H3N1O1,C,fix,any,Carbamidomethyl
O1,M,opt,any,Oxidation
C2H2O1,*,opt,Prot-N-term,Acetyl
```
But the Rust bench command never passed `--mod`, so Rust ran with built-in defaults (Cam-C + M-ox only). Acetylation was entirely absent from the candidate enumeration.

## The fix

Rust's `--mod` flag accepts numeric mass deltas (Java composition strings like `C2H2O1` are not yet supported). Created `/tmp/astral_mods_rust.txt`:
```
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
15.99491,M,opt,any,Oxidation
42.01057,*,opt,Prot-N-term,Acetyl
```

iter24 bench used `--mod /srv/data/msgf-bench/astral_mods_rust.txt`.

## Results

| Metric | iter22b | iter24 (+acetyl) | Δ |
|---|---:|---:|---:|
| Targets | 92,900 | 93,717 | +817 |
| Decoys | 56,488 | 55,812 | -676 |
| T/D ratio | 1.647 | 1.679 | better |
| **1% FDR** | 31,006 | **31,390** | **+384** |
| 5% FDR | 34,465 | 34,903 | +438 |
| Acetyl PSMs in Rust | 0 | 3,281 | new |
| Protein-N-term peptides | 0 (`_.` flanking absent) | 4,313 | new |
| Label flips (Java target, Rust decoy) | 16,956 | 16,437 | -519 |

Rust now enumerates **85%** as many acetyl-modified PSMs as Java (3,281 vs 3,868) and **89%** as many protein-N-terminal peptides (4,313 vs 4,819).

Gap to Java's 35,818: **13.5% → 12.4%**.

## Why the PSM gain is "only" +384 (8.4% × 16K flips ≈ 1,400 expected)

Per the n=9 audit, adding candidates dilutes Percolator's signal because:
1. Many recovered acetyl PSMs don't cross 1% FDR under Percolator's learned weights
2. Some recovered acetyls displaced previously-top-1 (non-acetyl) PSMs that DID cross 1% FDR
3. NumMods=3 (vs Rust's default 2) adds even more variants → more decoy candidates to compete

Net: 16,956 → 16,437 jtrd flips (-519), but only +384 PSMs at 1% FDR. The other 135 acetyl-recovered PSMs cancel out via decoy displacement.

## Bonus: notation difference (not a bug)

My initial audit reported "0% Rust protein-N-terminal peptides" based on substring search for `-.` flanking (Java's format). Rust uses `_.` flanking for protein-N-terminus — same enumeration semantics, different printed character. Once accounted for, iter22b had 4,313 protein-N-terminal peptides (after iter24 fix), close to Java's 4,819.

## Where to ship

The fix is at the **bench harness** level, not Rust source code. The `astral_mods_rust.txt` file needs to be checked in (or generated from Java's `mods.txt` via a converter) and the Astral bench script updated to pass `--mod`. iter24's Rust binary is identical to iter22b (no code changes).

This means:
- Production-Astral users who copy Java's `mods.txt` semantics will gain ~1% FDR by passing `--mod`
- The Rust mod parser already supports numeric mass deltas; only the bench script was missing the flag

## Remaining gap analysis

Java: 35,818. Rust iter24: 31,390. **Gap: 4,428 PSMs (12.4%).**

Per the n=9 audit + iter17/18/23 evidence:
- Most remaining flips are SCORING divergences (RawScore median -2, lnSpecEValue -0.72)
- Closing them requires score_psm changes that have regressed Percolator (edge-scoring iter17/18)
- Per-feature parity also regresses Percolator (iter23)

Remaining levers:
1. **Native SpecEValue-only FDR** (skip Percolator) — would let score_psm fixes actually help
2. **More mod variants** (carbamidomethyl-M, deamidation, etc. if Java's full search uses them) — diminishing returns
3. **Audit the 14,059 rtjd flips** — cases where Rust picks target and Java picks decoy. Could be where Rust's looser candidate set adversely affects Percolator's calibration

## Commit

No Rust code change. iter22b binary used (`/srv/data/msgf-bench/track-iter22b-build/`). Update needed:
- `benchmark/run_astral_3arm.sh` (or equivalent local harness): add `--mod path/to/astral_mods_rust.txt` to the Rust invocation line.
