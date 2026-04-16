# Skill: MS-GF+ Score Output Safety for Rescoring Workflows

## Context

MS-GF+ outputs several scores in mzIdentML and TSV formats. When results are
passed to ML-based rescoring tools (Percolator, MS2Rescore, Oktoberfest), some
scores **leak target/decoy information**, making FDR estimation unreliable.

This was identified at EuBIC-MS Symposium 04/2026 by Henry Emanuel Weber
(Ruhr-Universität Bochum). See investigation 002 for full details.

## Score Safety Classification

### SAFE — can be passed to rescoring tools
- **RawScore** (`MS:1002049`) — integer score from generating function
- **DeNovoScore** (`MS:1002050`) — integer de novo score
- **SpecEValue** (`MS:1002052`) — spectral-level E-value from generating function

### UNSAFE — must NOT be passed to rescoring tools
- **EValue** (`MS:1002053`) — `SpecEValue × numDistinctPeptides`. The database-size
  multiplier introduces target/decoy asymmetry that Percolator can exploit for
  100% separation of target and decoy distributions.
- **QValue** (`MS:1002054`) — computed directly from TDA (target/decoy counting).
  This is literally the target/decoy separation encoded as a number.
- **PepQValue** (`MS:1002055`) — same as QValue but at peptide level.

## When Modifying Score Output Code

### Files that write scores
1. `MZIdentMLGen.java` — mzIdentML output (lines ~345-421)
2. `DirectTSVWriter.java` — TSV output (lines ~138-208)
3. `DBScanner.java` — MSGFDB TSV output (lines ~850-915)
4. `MSGFDBResultGenerator.java` — result generation (lines ~92-104)

### Rules
- Never add EValue, QValue, or PepQValue as features for ML-based rescoring
- When adding a `-rescoring` or `--percolator-safe` mode, omit MS:1002053/54/55
- SpecEValue (MS:1002052) is always safe — it's per-spectrum, no TDA dependency
- RawScore and DeNovoScore are always safe — integer scores, no database info

### E-value computation (for reference)
```java
// MZIdentMLGen.java:346-347
int numPeptides = sa.getNumDistinctPeptides(enzyme == null ? length - 2 : length - 1);
double eValue = specEValue * numPeptides;
```

The `numDistinctPeptides` comes from `CompactSuffixArray`, which counts over the
full concatenated target+decoy database suffix array.

### Q-value computation (for reference)
```java
// ComputeFDR.java:272-276
float psmQValue = tda.getPSMQValue((float) m.getSpecEValue());
Float pepQValue = tda.getPepQValue(m.getPepSeq());
m.setPSMQValue(psmQValue);
m.setPepQValue(pepQValue);
```

`TargetDecoyAnalysis` separates PSMs by protein prefix (target vs decoy),
sorts by score, and computes FDR = decoyCount / targetCount. This directly
encodes target/decoy status.
