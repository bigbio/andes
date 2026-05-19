# Astral PIN diff: Rust (iter11 R-2 baseline) vs Java — empirical findings

_2026-05-19. First run of `benchmark/parity/analyze_rust_java_pin_diff.py`
on the Astral LFQ benchmark PINs: Java `astral-java.pin` (136,271 rows,
121,654 scans) vs Rust `astral-rust-iter11.pin` (149,497 rows, 121,665
scans, R-2 baseline reverted to b1d45bb era + R-2 retention)._

The report localizes the ~30% Astral 1% FDR gap into specific PSM
buckets and specific feature columns, replacing source-reading guesses.

## Top-1-per-scan label disagreement: 25% of scans

| Bucket | Count | % |
|---|---:|---:|
| both target, same peptide (agreement) | 45,014 | 37.0% |
| both decoy (agreement, FDR-neutral) | 26,299 | 21.6% |
| both target, different peptide (ranking flip) | 19,372 | 15.9% |
| Java target, Rust decoy (label flip) | 16,990 | 14.0% |
| Rust target, Java decoy (label flip) | 13,963 | 11.5% |
| Java-only target (Rust missed) | 12 | 0.01% |
| Rust-only target (Java missed) | 16 | 0.01% |
| one-sided decoys | 15 | 0.01% |
| **total scans** | **121,681** | 100% |

**Headlines:**

- Only 37% of scans agree on the top-1 target peptide. 16% have ranking
  flips (same scan, different peptide is "best" for each engine).
- **25.5% of scans have a label flip** (top-1 target on one side, top-1
  decoy on the other). This is the dominant FDR divergence — half of
  it (14%) is Java rescuing PSMs Rust labels decoy; the other half
  (11.5%) is Rust calling targets that Java calls decoys.
- The "one-sided" miss buckets (Java-only or Rust-only) are tiny
  (<0.05%). Both engines see essentially the same scan set — divergence
  is in ranking and scoring, not in candidate enumeration.

## Per-feature divergence on the agreement bucket (49,538 PSMs)

Restricting to scans where both engines picked the same top-1 target
strips ranking-flip + retention noise. The per-feature deltas below
measure **feature-formation divergence specifically** — for the same
PSM, what does each engine emit?

| Feature | median Δ | p5 → p95 | %frac \|relΔ\|>1% | Tier |
|---|---:|---|---:|---|
| **enzC** | **-1** | -1 → -1 | 100% | 🔥 zero-stubbed in Rust |
| **enzN** | **-1** | -1 → -1 | 100% | 🔥 zero-stubbed in Rust |
| **enzInt** | 0 | -2 → 0 | 100% | 🔥 zero-stubbed in Rust |
| **MS2IonCurrent** | +78,570 | +16k → +296k | 98.6% | ⚠️ Rust higher by 1.6× |
| **StdevRelErrorTop7** | +41 ppm | -0.8 → +380 | 99.4% | ⚠️ scattering, units? |
| **MeanRelErrorTop7** | +1.3 ppm | -118 → +117 | 99.9% | ⚠️ huge spread |
| **lnSpecEValue** | +0.74 | -7 → +30 | 97.8% | ⚠️ Rust ~2× larger spec_e_value typical |
| **lnEValue** | -3.88 | -11.6 → +20.9 | 99.7% | ⚠️ bimodal (HIGH-2 narrowed this) |
| **RawScore** | -2 | -21 → +12 | 96.4% | ⚠️ Rust scores 2 pts lower median |
| **DeNovoScore** | -9 | -29 → +1212 | 98.8% | ⚠️ Rust max_score lower median, tail wider |
| **MeanErrorTop7** | -4.7 | -10.6 → -3.1 | 100% | ❓ likely units mismatch |
| **NumMatchedMainIons** | +3 | -1 → +9 | 93.9% | ⚠️ Rust counts more ions |
| **StdevErrorTop7** | -2.0 | -6.0 → -0.7 | 100% | ❓ likely units mismatch |
| **longest_b** | +2 | -1 → +6 | 84.2% | ⚠️ Rust longer b-runs |
| **longest_y** | 0 | -3 → +5 | 61.8% | ⚠️ moderate divergence |
| **matchedIonRatio** | +0.46 | +0.03 → +0.92 | 99.3% | downstream of NumMatched+pepLen |
| **longest_y_pct** | -0.04 | -0.3 → +0.46 | 99.8% | downstream of longest_y |
| **lnDeltaSpecEValue** | 0 | 0 → 0 | — | match for rank-1 only PSMs |
| **ExplainedIonCurrentRatio** | -0.023 | -0.12 → +0.04 | 99.5% | ⚠️ Rust slightly lower |
| **CTermIonCurrentRatio** | -0.019 | -0.09 → +0.03 | 99.6% | ⚠️ Rust slightly lower |
| **NTermIonCurrentRatio** | -0.003 | -0.04 → +0.01 | 99.0% | ⚠️ Rust slightly lower |
| **dm / absdm** | ~0 | ~±5e-5 | 60.9% | ✅ basically agreeing |
| **isotope_error, peplen, charge2/3/4** | 0 | 0 → 0 | 0% | ✅ perfect match |
| **IsolationWindowEfficiency** | 0 | 0 → 0 | — | both stub to 0 |

## Diagnosis

### Showstopper (single biggest fix)

**Rust emits 0 for enzN, enzC, and enzInt on every PSM.** Java emits
1/1/(real count). These are three of Percolator's most informative
discriminator features (peptide N-cleavage, C-cleavage, internal
cleavage consistency). The earlier `docs(parity-analysis)` note
explicitly marked them "zero-stubbed; Would require `Enzyme::is_cleavage_site`
wiring; deferred." The diff data shows this isn't a minor deferral —
it's actively stripping discriminator signal Percolator was designed
to use.

Implementation cost: small (mirror Java's `DirectPinWriter.isEnzymaticBoundary`
+ `countInternalEnzymatic`, then thread pre/post flanking residues into
the PIN write site). Expected impact: large per the diff data, since
these are the only three constant-0 features. Was deferred per
piecewise-risk concerns; the diff harness eliminates that risk because
this is an additive feature change, not a modification of an existing
feature distribution.

### Major scoring divergence

`RawScore` median Δ = -2 (Rust scores 2 points lower than Java on
average for the SAME peptide on the SAME spectrum). `lnSpecEValue`
median Δ = +0.74 (Rust's spec_e_value is ~2× Java's for the same PSM).
`DeNovoScore` median Δ = -9 (Rust's GF max_score is 9 points lower).

These three covary: a lower RawScore feeds a lower DeNovoScore and a
larger spec_e_value (PSM looks less significant). The root cause is
upstream of the PIN formatting — somewhere in `score_psm` or the GF
DP. Not a single-line fix; this is the structural scoring gap the
divergence audit catalogued.

### Likely units/formula bugs

`MeanErrorTop7` and `StdevErrorTop7` show 100% relative-deviation rate
with median Δ = -4.7 and -2.0 respectively. The relative-error variants
(`MeanRelErrorTop7`, `StdevRelErrorTop7`) are even more divergent —
StdevRelErrorTop7 median +41 ppm, p95 +380. This pattern (different
scales, not just different values) suggests Rust and Java use different
units or different filter sets for these top-7 statistics. Identifying
which side is correct is a small per-feature investigation.

### Fragment-ion enumeration divergence

`NumMatchedMainIons` median Δ = +3 (Rust counts 3 more matched ions
than Java for the same peptide). `longest_b` +2, `longest_y` +0.
`matchedIonRatio` and `longest_y_pct` are downstream of these counts.

This suggests Rust's fragment-ion prediction or matching tolerance is
more permissive than Java's. Could be the deconvolution step (already
landed at 601b45ff) overlapping Java's enumeration in a non-bit-exact
way, or a different ion-type set (b/y only vs b/y/c/z), or a different
tolerance.

### Intensity-ratio drift

`ExplainedIonCurrentRatio` median -0.023, `CTermIonCurrentRatio` -0.019,
`NTermIonCurrentRatio` -0.003. Rust's matched-intensity fractions are
all SMALLER than Java's. Combined with `MS2IonCurrent` being LARGER on
Rust by 1.6×, the picture is: **Rust's denominator (total MS2 intensity)
is larger than Java's**, which compresses every ratio. This likely
reflects different MS2 peak selection — Rust uses raw peaks, Java may
filter or rank-truncate.

## Recommended next implementation target

**Implement C-4 (enzN/enzC/enzInt) as the next single change.** It is
the only divergence in the table that is structurally additive — Rust
literally does not compute these features. Implementing them is a
strict gain in information for Percolator (no risk of disrupting an
existing feature distribution, unlike C-5b which rescaled an existing
feature). The earlier deferral was tied to a piecewise-fix-risk concern
that doesn't apply here: this is an ADD, not a CHANGE.

After C-4, the next data-driven candidate is the units bug on
MeanErrorTop7 / StdevErrorTop7 — quick to investigate, likely a small
fix. After that, the structural scoring divergence (RawScore / GF
max_score / spec_e_value) is the biggest remaining gap but the
hardest fix.

## Files

- `benchmark/parity/analyze_rust_java_pin_diff.py` — the harness
- `benchmark/parity/README.md` — usage notes
- Per-PSM CSV: `/tmp/parity-diff/iter11-vs-java/per_psm_diff.csv` (8.7 MB,
  not committed)
- This findings doc: `docs/parity-analysis/notes/2026-05-19-pin-diff-findings.md`

## Reproducibility

```bash
scp root@pride-linux-vm.ebi.ac.uk:/srv/data/msgf-bench/bench-3ds-results/astral-java.pin /tmp/parity-diff/astral-java.pin
scp root@pride-linux-vm.ebi.ac.uk:/srv/data/msgf-bench/bench-iter11-results/astral-rust-iter11.pin /tmp/parity-diff/astral-rust-iter11.pin
python3 benchmark/parity/analyze_rust_java_pin_diff.py \
    --java /tmp/parity-diff/astral-java.pin \
    --rust /tmp/parity-diff/astral-rust-iter11.pin \
    --out-dir /tmp/parity-diff/iter11-vs-java
```
