# Gap-corpus pipeline — design (2026-06-08)

Build an **MS-GF+-free** training corpus for the slugs where andes underperforms,
train own scoring/intensity models, and merge them once they pass an
entrapment-validated yield gate. This is the data lever for Phase-3 independence
(Apache relicense), after every estimator-side lever was refuted (see
[`2026-06-04-phase3-own-models.md`](2026-06-04-phase3-own-models.md) and the
session memory: noise_pseudo / temperature / isotonic / dedup all refuted; the
~3% own-vs-curated TMT gap is a corpus-quality problem, not an estimator knob).

Related: [`2026-06-08-algorithm-comparison-decision.md`](../../../internal-docs/2026-06-08-algorithm-comparison-decision.md).

## 1. Goal & success criteria

- **Target the gap slugs only** (where we are weak), not the full 39-model store:
  1. **TMT-CID** (`cid_lowres_tryp_tmt`) — the headline gap (a05058 / PXD016999).
  2. **low-res CID LFQ** (`cid_lowres_tryp`).
  3. **timsTOF** (TOF instrument class, e.g. `cid_tof_tryp`).
- **Success** = a trained, MS-GF+-free model for a slug whose yield ≥ curated on
  that slug's benchmark **with entrapment-FDP ≤ curated** (TDC alone is
  insufficient — proven this session). Diverse-pool model vs curated is the A/B.

## 2. Key decisions (from brainstorming)

| Decision | Choice |
|---|---|
| Scope | Gap slugs first (CID/TMT/low-res, then timsTOF) — "cover where we're not doing well" |
| Diversity | **Diverse pool** across instruments/labs per slug; entrapment A/B decides if it helps vs blurs |
| Label QC | **Strict & clean**: q≤0.01 + score-margin (unambiguous rank-1) + ≥6 matched b/y ions + ≥50% explained intensity + precursor-mass sanity + per-(peptide,charge) dedup |
| Identification engine | **MSFragger** for TMT/CID (mzML); **Sage** for timsTOF `.d` (native, Apache-2.0) |
| Provenance | Engines produce *labels*, not model parameters → trained tables embed neither tool → Apache-independence preserved |

## 3. Hard operational constraint: disk-bounded streaming

**The bench VM has only a few GB free.** The pipeline MUST process **one file at a
time** and never accumulate large artifacts:

```
per dataset (≈3 files, breadth > depth):
  for each file:
    download ONE raw/.d  → /tmp
    identify (MSFragger mzML  |  Sage .d)
    strict-clean QC select
    convert to flat  (small: KB–few MB)
    append a trace row
    DELETE the raw/.d + mzML + pepXML/Sage outputs immediately   <-- before next file
  # only the per-file flat parquet(s) survive (small)
```

- Only the **flats** persist (and the final merged `models.parquet`). Raw `.raw`/
  `.d`/`.mzML`/`.pepXML` are deleted right after the flat is written.
- `df` guard: abort a file if free disk < (file size × safety factor); never two
  raws on disk at once.
- ~3 files/dataset balances per-instrument signal against download/disk cost; the
  diverse pool gets breadth from *many datasets*, not many files per dataset.

## 4. Components

1. **Corpus config** (`gap_corpus.tsv` / `.jsonl`) — declarative list of
   `slug, dataset_accession, file_urls (≤3), engine, search_params, mods`. Curated
   by extending the MSnet catalog with a hand-picked PRIDE gap-list (diverse
   CID-TMT / low-res CID / timsTOF datasets). One row drives one dataset's run.
2. **Corpus driver** (`build_gap_corpus.sh`) — extends `phase_v_cid_train_models.sh`:
   reads the config, runs the disk-bounded per-file loop, accumulates flats per
   slug, writes the trace. Idempotent / resumable (skips datasets whose flat exists).
3. **`sage_to_flat`** (new) — Sage results TSV + `.d` spectra → flat, reusing
   **andes's native `.d` reader** so timsTOF peaks are ranked the same way at train
   and search time. (mzML+pepXML→flat already exists and is fixed:
   `mzml_pepxml_to_flat.py`.)
4. **Trainer** — existing `train-intensity` + `train-from-msnet`, peak filter on
   (auto for isobaric); pool all of a slug's flats.
5. **Validator** — existing entrapment-FDP harness on the slug's held-out run.
6. **Merger** — write the trained slug's rows into a fresh `models.parquet` with
   **zero MS-GF+-derived rows for that slug**; keep curated rows for untouched slugs
   until they too are replaced.

## 5. Tracing (reproducibility)

A `corpus_trace.tsv` row per processed file:
`slug, dataset, file, url, bytes, engine, raw_psms, qc_kept, flat_rows, q_thresh,
deleted_ok, timestamp`. Plus a per-slug summary (`total flats, pooled PSMs,
unique peptides, instruments covered`). The trace is the audit trail proving each
flat's provenance and that raws were deleted — and lets us reason about coverage
without re-downloading.

## 6. Per-engine identification params

- **MSFragger (TMT/CID, mzML)**: low-res CID (`fragment_mass_tolerance` 0.4–0.6 Da,
  `theoretical_fragment_ions` low-res), trypsin 2MC, len 6–40, z2–4, iso −1..2;
  TMT slug adds fixed TMT6plex (K, n-term). Output pepXML → `mzml_pepxml_to_flat.py`.
- **Sage (timsTOF, `.d`)**: native `.d`, low-res/TOF fragment tol, trypsin, same
  enzyme/length/charge; `--write-pin`/TSV → `sage_to_flat`. Apache-2.0, clean.

## 7. Build order

1. Build `sage_to_flat` + the disk-bounded corpus driver + config schema; unit-test
   on **one** small file end-to-end (download→search→flat→delete), verify the trace
   and that disk returns to baseline.
2. **Phase 1 (TMT-CID):** curate ~8–15 diverse CID-TMT PRIDE datasets (×3 files),
   run, pool, train, entrapment-A/B vs curated on a05058 + PXD016999 held-out.
3. **Phase 2 (low-res CID LFQ):** same, validate vs curated (e.g. UPS1/PXD001819 +
   a held-out).
4. **Phase 3 (timsTOF):** Sage path, validate on a timsTOF holdout.
5. Merge each slug that passes the gate; update the NOTICE/independence status as
   slugs become MS-GF+-free.

## 8. Risks / open

- **Diversity may blur** rather than help (earlier same-dataset attempts plateaued
  ~3% short). The entrapment A/B in Phase 1 is the early go/no-go — if a diverse
  pool still trails curated, fall back to instrument-homogeneous models or escalate
  to ProteomeTools/MassIVE-KB-scale data.
- **timsTOF flat fidelity**: `sage_to_flat` must rank `.d` peaks identically to
  andes search; validate by re-scoring a Sage-labeled `.d` file with andes and
  checking PSM agreement before trusting the flats.
- Strict QC may starve rare partitions → rely on the existing backoff; watch
  per-partition counts in the trace.
