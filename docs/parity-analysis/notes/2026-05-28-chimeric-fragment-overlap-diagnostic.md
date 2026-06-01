# Chimeric fragment-overlap diagnostic + BSA preview

**Date:** 2026-05-28
**Branch:** `feat/chimeric-dda-plus`

## Why

The chimeric Phase-1/2 post-mortem concluded the missing primitive is
**fragment-level shared-fragment competition (Phase 3)** — the premise being that
co-emitted peptides on a chimeric scan inflate FDR because the spurious runner-up
*claims the same MS2 peaks* as the real top peptide ("fragment theft"). That claim
was **reasoned, never measured.** This diagnostic measures it.

## Tool

Env-gated (`MSGF_CHIMERIC_OVERLAP=1`), `--chimeric`-only diagnostic in
`match_engine` (`run_chunk_inner`): for each scan that emits ≥2 distinct peptides,
it computes each of the top-2 peptides' matched charge-1 b/y **peak set** (reusing
the exact `compute_psm_features` matching via the `matched_peak_keys` helper) and
prints `CHIM_OVERLAP spec_idx=.. nA=.. nB=.. shared=.. jacc=.. fracmin=..`
(jaccard and shared/min). Zero production impact (double-gated on `--chimeric` +
env; the helper is never called otherwise).

Aggregate: `grep CHIM_OVERLAP <stderr> | awk ...` over scans with `nA,nB>0`.

## BSA preview (validation run)

`MSGF_CHIMERIC_OVERLAP=1 msgf-rust --spectrum test.mgf --database BSA.fasta
--chimeric --top-n 5`:

- 309 scans where both peptides matched ≥1 peak.
- **mean jaccard 0.08, mean fraction-of-smaller 0.12; only 13% of scans ≥0.5 overlap.**

**Low overlap** — tentatively *challenges* the fragment-theft hypothesis: the
spurious runner-ups appear to match their **own** (coincidental) peaks, not steal
the top peptide's. If that holds on real data, the chimeric inflation is "random
peptides finding enough coincidental matches in peak-rich spectra," NOT
shared-fragment borrowing — in which case **Phase 3 (greedy shared-fragment
removal) would NOT fix it** (removing the top peptide's peaks barely touches the
runner-up).

**Caveat (important):** BSA is a single-protein fixture with **no real
co-isolation** and tiny matched-peak sets — unrepresentative. This validates the
*tool*, not the hypothesis. The decisive measurement is Astral chimeric (dense
spectra, real wide-window co-isolation).

## Next (needs VM)

Run on Astral chimeric and aggregate the `CHIM_OVERLAP` distribution:
```
MSGF_CHIMERIC_OVERLAP=1 <astral chimeric run> 2> astral-overlap.log
grep CHIM_OVERLAP astral-overlap.log | awk '{for(i=1;i<=NF;i++){split($i,a,"=");v[a[1]]=a[2]}
  if(v["nA"]>0&&v["nB"]>0){n++;sj+=v["jacc"];sf+=v["fracmin"];if(v["fracmin"]>=0.5)hi++}}
  END{printf "scans=%d mean_jacc=%.3f mean_fracmin=%.3f hi50=%.0f%%\n",n,sj/n,sf/n,100*hi/n}'
```
- **High overlap (fracmin ≳ 0.5 common)** → theft confirmed → Phase 3 is the validated fix.
- **Low overlap (BSA pattern holds)** → theft refuted → the inflation is multiple-testing
  on coincidental matches → Phase 3 won't help; chimeric needs per-scan/peptide-level
  FDR (or stays shelved).

## Astral result (2026-05-29, decisive run)

`MSGF_CHIMERIC_OVERLAP=1` Astral chimeric run (chimeric-build rebuilt from commit
`59421180`'s `match_engine.rs`; HCD/QExactive, cal=auto, --chimeric). 123,385
`CHIM_OVERLAP` lines emitted; **121,423 scans** with both top-2 peptides matching ≥1 peak.

**Aggregate: mean jaccard 0.173, mean fraction-of-smaller 0.367, fracmin≥0.5 = 38%,
fracmin≥0.25 = 61%, jacc≥0.5 = 9%.**

fracmin decile histogram (n=121,423) — **bimodal**:

| bin | count | % |
|---|---|---|
| [0.0-0.1) | 34,462 | 28.4% |
| [0.1-0.2) |  6,668 |  5.5% |
| [0.2-0.3) | 15,859 | 13.1% |
| [0.3-0.4) | 13,327 | 11.0% |
| [0.4-0.5) |  4,438 |  3.7% |
| [0.5-0.6) | 17,612 | 14.5% |
| [0.6-0.7) |  8,628 |  7.1% |
| [0.7-0.8) |  3,507 |  2.9% |
| [0.8-0.9) |  3,049 |  2.5% |
| [0.9-1.0) | 13,873 | 11.4% |

**Verdict: the BSA low-overlap pattern does NOT hold on real co-isolated data.**
Astral mean fracmin 0.367 is ~3× BSA's 0.12; fracmin≥0.5 is 38% vs BSA's 13%. The
distribution is bimodal — a coincidental-match mode (~28% at [0.0-0.1)) **plus** a
substantial high-overlap mode peaking near-total at [0.9-1.0) (11.4%). BSA was
unrepresentative exactly as the caveat warned.

**Fragment theft is real for a substantial fraction (~38% with ≥50% overlap, 61%
with ≥25%) of co-emitting chimeric scans on Astral** → **Phase 3 (greedy
shared-fragment competition) is validated** as a relevant mechanism. It is NOT
universal — the ~28% coincidental population means Phase 3 alone won't fully fix
chimeric FDR inflation; per-scan/peptide-level FDR is still needed for the
coincidental-match tail. But the core premise that motivated Phase 3 holds.
</content>
