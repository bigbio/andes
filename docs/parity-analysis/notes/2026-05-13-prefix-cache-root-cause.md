# prefix_score_cache[1087] == 0.0 root cause (PXD001819 scan=28787)

Localized 2026-05-14 via diagnostic example
`rust/crates/scoring/examples/dump_prefix_cache.rs`. Diagnostic took ~3 s
on a release build. No production source change required (read-only:
relies on `cached_prefix_score`, `nearest_peak_rank`, `partition_ion_logs`,
and `param.segment_num` — all already `pub`).

## 1. Diagnostic output (verbatim, key excerpts)

```
== Param ==
num_segments     = 2
max_rank         = 150
mme              = Da(0.5)

== Per-segment partitions for THIS spectrum ==
seg=0 partition=(c=2, pm=1971.038, seg=0) total_ions=5 prefix=3 suffix=2
    Suffix(c=1,off=19.01839)
    Prefix(c=1,off=1.00783)        ← b+ (PROTON)
    Suffix(c=1,off=20.02175)
    Prefix(c=1,off=-17.00274)      ← b-NH3
    Prefix(c=1,off=-26.98709)      ← b-???
seg=1 partition=(c=2, pm=1971.038, seg=1) total_ions=3 prefix=0 suffix=3
    Suffix(c=1,off=19.01839)
    Suffix(c=1,off=20.02175)
    Suffix(c=1,off=21.02163)

== nominal_mass = 974.0 (is_prefix=true) ==
  seg=0 ion=Prefix(c=1,off=1.00783)  theo_mz=975.4980  seg(theo)=0 (IN)  matched rank=27  score=2.6511
  seg=0 ion=Prefix(c=1,off=-17.00274) theo_mz=957.4874 seg(theo)=0 (IN)  matched rank=125 score=0.7350
  seg=0 ion=Prefix(c=1,off=-26.98709) theo_mz=947.5031 seg(theo)=0 (IN)  matched rank=74  score=1.0874
  -> replicated_total=4.4735  cached_prefix_score(974)=Some(4.4734993)

== nominal_mass = 1087.0 (is_prefix=true) ==
  seg=0 ion=Prefix(c=1,off=1.00783)  theo_mz=1088.5549 seg(theo)=1 (OUT) SKIP(seg mismatch)
  seg=0 ion=Prefix(c=1,off=-17.00274) theo_mz=1070.5443 seg(theo)=1 (OUT) SKIP(seg mismatch)
  seg=0 ion=Prefix(c=1,off=-26.98709) theo_mz=1060.5600 seg(theo)=1 (OUT) SKIP(seg mismatch)
  (seg=1 iteration emits NO prefix lines — seg=1 partition has zero prefix ions)
  -> replicated_total=0.0000  cached_prefix_score(1087)=Some(0.0)
```

Identical zero-output for nominal in {1216, 1345, 1561, 1920} — every
prefix ion's theo_mz lands in seg=1 yet the seg=1 partition has no
prefix ions in its `frag_off_table` entry.

A wider scan confirmed: for charge=2, **every seg=1 partition in this
param file holds 0 prefix ions** (pfx=0, sfx=5 across ~30 charge=2
seg=1 rows of `param.partitions`). Prefix ions for charge=2 only have
rank-distribution tables in seg=0.

## 2. Root cause (one sentence)

`directional_node_score_inner` at
`rust/crates/scoring/src/scoring/scored_spectrum.rs:555-589` iterates
prefix ions only inside the segment whose partition's `frag_off_table`
contains them (seg=0 for charge=2), and then `segment_num(theo_mz,
parent_mass) != seg` skip at line 574 throws away every prefix ion
whose theoretical m/z lies in the upper half of the precursor range —
the seg=1 partition has no prefix ions to recover them, so the
contribution silently bottoms to 0.

## 3. Why the 297→108 gap

For this peptide (parent_mass≈2067) the b1+ theo_mz transitions from
seg=0 to seg=1 at roughly `nominal ≈ (parent_mass/2) * INTEGER_MASS_SCALER
≈ 1033`. So *every* split with prefix_nominal ≥ ~1033 contributes 0 to
the prefix score in Rust, while Java produces sizeable per-split
contributions (Java reported 20.22, 17.13, 13.92, 17.31, 22.98, -5.44,
16.38, 13.61 at nominal 1087, 1216, 1345, 1460, 1561, 1658, 1757, 1920
respectively). Summing those splits closes most of the 297 - 108 = 189
RawScore gap. The split at 974 (still in seg=0) gives 3.66 in Rust vs
~11.30 in Java; the residual 974-side gap is a separate (smaller)
issue, not the dominant cause.

Confirms **hypothesis #2** verbatim. Hypotheses #1 (partition seg=0
DOES have prefix ions), #3 (`logs.len() ≤ max_rank_idx` overflow), and
#4 (some other edge condition) all falsified — the matched seg=0
ions at nominal=974 score with idx in [26, 73, 124] and logs are
length 151 (= max_rank+1), so the bounds path is fine.

## 4. Recommended fix shape

Make the `seg` loop body score *every* prefix ion in its **correct**
segment, not only ions present in the current outer segment's
partition. Two viable shapes:

a. **Pre-compute a per-direction union of all prefix ions seen across
   any segment** (essentially `param.ion_types_for_segment` unioned
   across segments, deduped) and iterate that single list with the
   existing `segment_num(theo_mz) == ?` clause directing each ion to
   the partition whose seg matches `segment_num(theo_mz, parent_mass)`.
   This mirrors what Java appears to be doing — the b+ ion at nominal
   1087 must be scored against the seg=1 partition's rank-dist (or
   missing-ion slot) even though the seg=1 partition's ion list does
   not advertise it.

b. **Per (charge, parent_mass), fold seg=0's prefix-ion frag_off_table
   entries into the seg=1 partition** at param-load time so the
   existing inner loop stays unchanged. Cheaper to implement but
   semantically conflates two partitions' rank-dist tables — likely
   wrong.

Option (a) is correct; option (b) is a perf-shortcut variant. Either
way, side effects to consider: (i) `gf_java_parity` will move — the
prefix-score cache will gain non-zero entries that may change
`SpecEValue` for many spectra (intended), (ii) the
`missing_ion_score` path inside seg=1 will fire for prefix ions where
no peak is found, which can yield small *negative* contributions
matching Java's -5.44 pattern, (iii) the perf cost is a tiny constant —
the inner loop iterates one extra direction's ions per outer
segment.

## 5. Confidence

**High.** The replicated math inside the diagnostic matched
`cached_prefix_score` exactly (4.4735 vs 4.4734993 at nominal=974; 0.0
at the other six). The mechanism is mechanical: empty seg=1 prefix-ion
list + seg-mismatch skip = forced 0. Confirmed across 6 distinct
nominal masses spanning the seg=0→seg=1 boundary and beyond.
