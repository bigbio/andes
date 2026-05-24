# msgf-rust bug review (2026-05-23)

Branch: `review/bug-hunt` (from `master` @ 18360a3d)

Systematic review of the Rust MS-GF+ port: static analysis of critical paths,
full `cargo test --release --workspace`, and targeted code reading.

## Fixed in this branch

| ID | Severity | Location | Issue | Fix |
|---|---|---|---|---|
| B1 | **Critical** | `msgf-rust.rs` `send_chunks` | Bench cap (`--max-spectra N`) truncated the final partial chunk to zero when `total == N` (e.g. N=100 with chunk size 5000 → empty output). | Removed erroneous tail `truncate` block; loop already stops at cap. |
| B2 | **High** | `msgf-rust.rs` param routing | Activation auto-detect was gated on `instrument == low-res`, so `--fragmentation auto --instrument QExactive` on mzML skipped peek and resolved to CID params for HCD data. | Gate auto-route on `fragmentation == auto` + mzML extension only. |
| B3 | **High** | `msgf-rust.rs` TSV write | `write_tsv(..., is_mgf=true)` always emitted MGF layout (extra `Title` column) even for mzML inputs. | Pass `!is_mzml`. |
| B4 | **High** | `match_engine.rs` GF | SpecE GF graph used `start_offset == 0` for protein N-term instead of `cand.is_protein_n_term`, breaking Met-cleaved N-termini at offset 1. | Use `cand.is_protein_n_term` / `is_protein_c_term`. |
| B5 | **Medium** | `tsv.rs` | `IsotopeError` column hardcoded to 0 while PIN writes `psm.isotope_offset`. | Thread isotope offset from PSM. |
| B6 | **Medium** | `msgf-rust.rs` CLI | Inverted `--charge-min/--charge-max` or isotope ranges produced empty ranges with no error. | Validate at startup and return clear error. |
| B7 | **High** | `match_engine.rs` dedup | Dedup used bare sequence + pin score; merged mod variants incorrectly. | Mod-aware pepSeq key + `rank_score`. |
| B8 | **Medium** | `match_engine.rs` dedup | HashMap survivor order was nondeterministic. | `BTreeMap` + best-`rank_score` survivor rule. |

## Open — not fixed (documented for follow-up)

| ID | Severity | Location | Issue |
|---|---|---|---|
| B9 | **Low** | `sa_walk.rs` | Test-only SA walk helper does not enforce `max_missed_cleavages`; production search uses `candidate_gen::enumerate_candidates`, which does. |
| B10 | **High** | `mzml.rs` `Iterator::next` | First per-spectrum parse error sets `done=true` and aborts the entire file; remaining spectra are silently skipped. |
| B11 | **Low** | `sa_walk.rs` Met pass | Dedupes Met-cleaved peptides on residue bytes only, collapsing distinct C-terminal contexts. |

## Known test failures (pre-existing, CI-skipped)

These fail on `master` without the 7 CI skip flags; tracked as parity/min_peaks regressions:

- `match_engine_smoke::known_peptide_appears_in_top_n`
- `match_engine_smoke::charge_missing_spectrum_uses_per_charge_scored_spec`
- `match_engine_smoke::spectrum_without_charge_tries_charge_range`
- Maven fixture loads, thread-determinism test (see `.github/workflows/ci.yml`)

## Verification

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts
```

## Performance (dedup pass)

- PepSeq dedup keys use integer mod units + `Arc` cache per candidate (avoids repeated string formatting).
- Per-charge `TopNQueue` map uses `FxHashMap<u8, _>` (typically 1–3 charges per spectrum).

## Documentation review (2026-05-24)

Fixes applied on this branch:

| Issue | Location | Fix |
|---|---|---|
| PIN column count said "28" | `README.md` | Corrected to 36 (default charge 2–3) + EdgeScore note |
| Auto-detect described "first spectrum" only | `README.md` | First 64 MS2 histogram; `--instrument` does not gate peek |
| Auto-detect required `--instrument low-res` | `DOCS.md` §4 | Matches code: only `--fragmentation auto` + mzML |
| TSV `IsotopeError` documented as always 0 | `DOCS.md` §3b | Updated after B5 fix |
| Broken `known-divergences.md` links | `README.md`, `DOCS.md` §8d | Legacy file removed in iter39; point to §8d / tests |
| Inverted charge/isotope ranges undocumented | `DOCS.md` §1 | Startup validation documented |

**Still stale (not fixed here):**

- `benchmark/ci/README.md` — references Java Maven workflow; no Rust benchmark workflow in `.github/workflows/` yet.
- `.claude/CLAUDE.md` — Java-tree context; accurate on `java-legacy` branch only.
