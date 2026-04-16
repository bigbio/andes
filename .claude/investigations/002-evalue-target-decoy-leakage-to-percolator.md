# Investigation 002: E-value Leaks Target/Decoy Information to Percolator

**Status:** OPEN
**Date reported:** 2026-04-15
**Severity:** HIGH — affects FDR estimation for all downstream rescoring tools
**Source:** EuBIC-MS Symposium 04/2026, Copenhagen — Henry Emanuel Weber, Ruhr-Universität Bochum (Jun.-Prof. Julien Urchueguía group)
**Slide screenshot:** `assets/Screenshot_2026-04-15_at_13.23.09-*.png`

## What Was Observed

When MS-GF+ results are passed to rescoring tools (Percolator, MS2Rescore, Oktoberfest),
the target and decoy score distributions become **completely separated** — 100% separation.
This does NOT happen with Comet results on the same data.

The presenter found that **removing the E-value (MS:1002053) from the MS-GF+ features
fixed the problem**, confirming that the E-value is the source of information leakage.

Key observations from the slide:
- **Comet + TDA/Percolator/MS2Rescore/Oktoberfest:** Normal overlapping distributions
- **MS-GF+ + TDA:** Normal overlapping distributions (E-value not used as feature)
- **MS-GF+ + Percolator/MS2Rescore/Oktoberfest:** Perfect separation (E-value used as feature)

## The Mechanism

### How MS-GF+ computes the E-value

The E-value is computed as:

```
E-value = SpecEValue × numDistinctPeptides
```

See `MZIdentMLGen.java:347`:
```java
double eValue = specEValue * numPeptides;
```

Where:
- **SpecEValue** (`MS:1002052`) = spectral-level E-value from the generating function
  (computed per spectrum, independent of target/decoy status)
- **numDistinctPeptides** = count of distinct peptide sequences of the matched length
  in the **entire** concatenated target-decoy database
  (from `CompactSuffixArray.getNumDistinctPeptides()`)

### Why it leaks

The `numDistinctPeptides` multiplier is derived from the suffix array built over the
**concatenated target+decoy database** (`-tda 1` mode). The count includes both target
and decoy peptides.

However, the critical issue is that `numDistinctPeptides` is looked up by **peptide
length** (see `CompactSuffixArray.java:138-140`):

```java
public int getNumDistinctPeptides(int length) {
    return numDistinctPeptides[length];
}
```

This is the same multiplier for targets and decoys of the same length, so the
E-value itself doesn't directly encode target/decoy status. The leakage likely
comes from a subtler mechanism:

**Hypothesis 1: Database-size asymmetry**
When `-tda 1` is used, MS-GF+ generates reversed decoys internally. The number
of distinct peptides at each length may differ slightly between the target and
decoy halves. Since the E-value uses the combined count, it implicitly encodes
information about the database composition. Percolator, being a machine learning
model, can learn to exploit even tiny systematic differences.

**Hypothesis 2: Score distribution coupling**
The generating function that produces SpecEValue is computed using score
distributions that are calibrated on the full database. If the score distribution
shape differs systematically between target and decoy hits (which it does — true
matches exist only for targets), the SpecEValue already carries some target/decoy
signal that gets amplified by the numPeptides multiplier.

**Hypothesis 3: Q-value propagation**
The Q-value (`MS:1002054`) is explicitly computed from TDA and directly encodes
target/decoy ranking. If Q-value is also passed to Percolator alongside E-value,
the combined features create a perfect classifier. However, the presenter
specifically identified E-value (not Q-value) as the problematic score.

**Hypothesis 4: E-value scale differences**
SpecEValue is a per-spectrum probability; E-value is SpecEValue × database_size.
Since all peptides (target and decoy) use the same `numDistinctPeptides[length]`,
the E-value is a monotonic transform of SpecEValue for peptides of the same
length. But across different lengths, the scaling differs, and Percolator could
learn length-dependent patterns that correlate with target/decoy status.

## Relevant Code

### E-value computation

- `MZIdentMLGen.java:345-347` — `eValue = specEValue * numPeptides`
- `DirectTSVWriter.java:138-141` — same computation for TSV output
- `DBScanner.java:853-854` — same computation for MSGFDB output
- `MSGFDBResultGenerator.java:92-104` — `getPValue()` and `getEValue()` static methods

### numDistinctPeptides lookup

- `CompactSuffixArray.java:138-140` — `getNumDistinctPeptides(length)`
- `CompactSuffixArray.java:196-228` — counting logic over suffix array
- `SuffixArrayForMSGFDB.java:43-46` — wrapper

### Scores written to mzIdentML

- `MS:1002049` — RawScore (integer, safe)
- `MS:1002050` — DeNovoScore (integer, safe)
- `MS:1002052` — SpecEValue (spectral E-value, probably safe)
- `MS:1002053` — EValue (database E-value, **LEAKS**)
- `MS:1002054` — QValue (from TDA, **inherently encodes T/D**)

## Impact

- **All rescoring workflows are affected:** Any tool that uses MS-GF+ E-value as a
  feature (Percolator, MS2Rescore, Oktoberfest) will produce artificially inflated
  identification rates
- **Published results may be affected:** Studies using MS-GF+ → Percolator pipelines
  may report overly optimistic PSM counts
- **FDR estimates are unreliable:** The 100% target/decoy separation means FDR
  cannot be meaningfully estimated

## Which Scores Leak?

### Safe scores (no target/decoy information)
| CV Accession | Name        | Why safe |
|-------------|-------------|----------|
| MS:1002049  | RawScore    | Integer score from generating function, per-spectrum |
| MS:1002050  | DeNovoScore | Integer de novo score, per-spectrum |
| MS:1002052  | SpecEValue  | Spectral E-value from generating function, per-spectrum. No TDA dependency. |

### Unsafe scores (leak target/decoy information)
| CV Accession | Name       | Why it leaks |
|-------------|------------|--------------|
| MS:1002053  | EValue     | `SpecEValue × numDistinctPeptides` — database-size multiplier may introduce asymmetry. Confirmed as the leak source by the presenter. |
| MS:1002054  | QValue     | **Directly computed from TDA** via `TargetDecoyAnalysis.getPSMQValue()` — it IS the target/decoy separation. Passing this to Percolator is giving it the answer key. |
| MS:1002055  | PepQValue  | Same as QValue but at peptide level. Also directly from TDA. |

### Q-value is categorically worse than E-value

The Q-value (`MS:1002054`) is computed by `TargetDecoyAnalysis.getFDRMap()` which:
1. Separates PSMs into target and decoy lists (by protein prefix, e.g. `XXX_`)
2. Sorts both by score
3. Walks down the ranked list computing `FDR = decoyCount / targetCount`
4. Converts FDRs to Q-values (monotonic minimum)

This is a **direct encoding** of target vs decoy status. If Percolator receives
QValue as a feature, it can trivially reconstruct whether a PSM is target or
decoy — far more directly than the E-value leakage. The EValue leakage is subtle
(the presenter had to investigate to find it); QValue leakage is by definition.

In practice, most rescoring tools (Percolator, MS2Rescore) likely skip QValue
because it's already an FDR estimate. But EValue looks like a "normal" search
engine score and gets picked up as a feature — which is why the EValue leak
is the one that actually manifests.

## Proposed Fix: Only Output SpecEValue (Omit EValue and QValue)

Since the downstream workflow is always `MS-GF+ → Percolator/rescoring tool → FDR`,
MS-GF+ does not need to output its own EValue or QValue. The rescoring tool will
compute its own FDR.

### What to change
1. **Stop writing EValue (MS:1002053) to mzIdentML** — or make it optional via CLI flag
2. **Stop writing QValue (MS:1002054) and PepQValue (MS:1002055)** — same treatment
3. **Keep SpecEValue (MS:1002052)** — this is the per-spectrum score, safe for rescoring
4. **Keep RawScore (MS:1002049) and DeNovoScore (MS:1002050)** — integer scores, safe

### Where to change
- `MZIdentMLGen.java:346-421` — mzIdentML output (remove/gate EValue, QValue, PepQValue CV params)
- `DirectTSVWriter.java:140-208` — TSV output (same)
- `DBScanner.java:853` — MSGFDB TSV output (same)
- `MSGFPlus.java` / `MSGFDB.java` — add CLI flag (e.g. `--no-evalue` or `--percolator-safe`)

### Impact on MSGFPlusAdapter (OpenMS)
The OpenMS `MSGFPlusAdapter` extracts scores from MS-GF+ mzIdentML output. If we
stop outputting EValue by default, the adapter needs to be updated to use SpecEValue
instead. This should be coordinated with the OpenMS team, or we add a CLI flag
so existing workflows keep working.

### Backward compatibility
- Add a flag like `-rescoring 1` that omits EValue/QValue from output
- Default behavior unchanged (EValue/QValue still written) for backward compat
- Document clearly that `-rescoring 1` should be used when piping to Percolator

## Next Steps

- [ ] Reproduce the issue: run MS-GF+ on a benchmark dataset, feed to Percolator,
      plot target/decoy distributions with and without E-value
- [ ] Contact Henry Emanuel Weber / Julien Urchueguía group for their test dataset
      and exact Percolator configuration
- [ ] Analyze whether SpecEValue alone also leaks (likely not, but should verify)
- [ ] Check if the leakage magnitude depends on database size (small DB = more leakage?)
- [ ] Review what scores MS2Rescore/Percolator extract from MS-GF+ mzIdentML by default
- [ ] Implement `-rescoring 1` CLI flag to omit EValue/QValue/PepQValue from output
- [ ] Coordinate with OpenMS team on MSGFPlusAdapter changes (use SpecEValue instead of EValue)
- [ ] Add skill documentation (DONE — see `.claude/skills/score-output-safety.md`)

## References

- Slide: "Target and decoy distributions" — EuBIC-MS Symposium 04/2026, Copenhagen
- Presenter: Henry Emanuel Weber, Medical Bioinformatics, Ruhr-Universität Bochum
- Group: Jun.-Prof. Julien Urchueguía
- Talk: "Leveling the playing field" (slide 9)
