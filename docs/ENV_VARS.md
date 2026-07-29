# `ANDES_*` environment variables

Every environment variable the engine reads, generated from the source and enforced
by `crates/search/tests/env_registry.rs` — a variable that is not listed here fails
the build.

This file exists because the codebase repeatedly shipped *corrections* that were
written, validated, and then left disabled behind an undocumented variable. A user
could not discover them, and neither could we: the Y-ladder size-bias correction sat
off by default long enough to be re-tuned around, and a c/z charge guard was written
and never called at all.

**Gating form matters.** `presence` means the variable is read with `var_os`, so
setting it to `0` *enables* it — a footgun. `value` means it requires an explicit
`"1"`. New variables should use `value`.

| Variable | Form | Defined at | Purpose |
| --- | --- | --- | --- |
| `ANDES_CHIMERIC_OVERLAP` | value | `crates/search/src/match_engine.rs:526` | _(undocumented — needs a description)_ |
| `ANDES_DENSE_NOISE` | value | `crates/model-train/src/accumulate.rs:110` | ANDES_DENSE_NOISE=<n> = dense random-position noise sampling (Kim et al., Nat Commun 5:5277, 2014 — sharper missing-slot-dominated noise; see dense_no |
| `ANDES_DENSITY_RAW` | presence | `crates/scoring/src/scoring/scored_spectrum.rs:769` | _(undocumented — needs a description)_ |
| `ANDES_ETHCD_AS_ETD` | presence | `crates/input/src/mzml.rs:172` | _(undocumented — needs a description)_ |
| `ANDES_GEO_MAX_FRAG_CHARGE` | ? | `crates/andes/src/bin/andes.rs:1221` | _(undocumented — needs a description)_ |
| `ANDES_GEO_MAX_RANK` | ? | `crates/andes/src/bin/andes.rs:1218` | _(undocumented — needs a description)_ |
| `ANDES_GEO_MAX_TIERS` | ? | `crates/andes/src/bin/andes.rs:1220` | _(undocumented — needs a description)_ |
| `ANDES_GEO_OCCUPANCY` | ? | `crates/andes/src/bin/andes.rs:1219` | _(undocumented — needs a description)_ |
| `ANDES_GEO_SEGMENTS` | ? | `crates/andes/src/bin/andes.rs:1217` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_CHARGE_PM1` | presence | `crates/search/src/glyco_search.rs:908` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_CZ_FIX_OFF` | presence | `crates/scoring/src/scoring/psm_score.rs:428` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_CZ_GATE_OFF` | presence | `crates/search/src/glyco_search.rs:923` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_CZ_INTENSITY` | presence | `crates/scoring/src/scoring/psm_score.rs:604` | ANDES_GLYCO_CZ_INTENSITY set, weight each matched ion by its observed base-peak-normalised intensity, making this a true explained-INTENSITY ratio. Un |
| `ANDES_GLYCO_CZ_MULTISITE` | presence | `crates/search/src/glyco_search.rs:938` | by `cz_score_best_site` too — the max-over-sites is applied symmetrically, so the per-candidate null is matched). Enabling it is gated on a decoy-cont |
| `ANDES_GLYCO_CZ_REMNANT` | value | `crates/scoring/src/scoring/fragment_ions.rs:266` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_CZ_SITE_LEGACY` | value | `crates/search/src/glyco_search.rs:365` | Requires an explicit "1" (matching `ladder_norm_enabled`): with a bare `is_some()`, `ANDES_GLYCO_CZ_SITE_LEGACY=0` would ENABLE the legacy resolver, w |
| `ANDES_GLYCO_CZ_STRUCT` | value | `crates/search/src/glyco_search.rs:1853` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_CZ_ZMAX` | value | `crates/scoring/src/scoring/psm_score.rs:445` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_ENUM_FALLBACK` | value | `crates/output/src/glyco_pin.rs:485` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_ETD_DBFALLBACK_OFF` | presence | `crates/search/src/glyco_search.rs:902` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_FULL_GLYCANS` | presence | `crates/andes/src/bin/andes.rs:2754` | hunt (2026-07-16) showed the default ~612 list MISSES the mouse-brain glycome at high charge (z5 69%/z6 38% coverage) — a generation ceiling. ANDES_GL |
| `ANDES_GLYCO_ISOBAR_REP` | value | `crates/andes-glyco/src/hybrid.rs:445` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_ISO_NEG` | presence | `crates/andes/src/bin/andes.rs:1934` | negative shift), while the iso=-1 arm emitted 28.5% of all candidate rows for 0.29% of the correct answers at a ~53:47 target:decoy ratio - pure FDR d |
| `ANDES_GLYCO_ISO_WIDE` | presence | `crates/andes/src/bin/andes.rs:1945` | monoisotopic peak mis-picked several 13C low, so the true neutral mass falls outside the default -1..=2 sweep. Widen the upper bound for glyco so that |
| `ANDES_GLYCO_LADDER_NORM` | value | `crates/andes-glyco/src/backbone.rs:651` | **Correction, now DEFAULT ON.** Divides the glycan-Y ladder by its rung count. The raw sum's expectation grows with glycan size, so it rewarded oversized glycans. Set to `0` to restore the biased estimator for A/B. Measured: mouse 664 -> 683 IDs @1%. |
| `ANDES_GLYCO_PAIR_RANK_ETD` | value | `crates/search/src/glyco_search.rs:1900` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_PAIR_RANK_GLYCAN` | value | `crates/search/src/glyco_search.rs:1927` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_PAIR_Y_ON_GEN` | presence | `crates/search/src/glyco_search.rs:1040` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_SCANS` | presence | `crates/search/src/glyco_search.rs:588` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_SEQUON_BOUNDARY` | value | `crates/search/src/glyco_search.rs:653` | _(undocumented — needs a description)_ |
| `ANDES_GLYCO_Y_HICHARGE` | presence | `crates/andes-glyco/src/backbone.rs:447` | _(undocumented — needs a description)_ |
| `ANDES_PEAK_PER_WINDOW` | value | `crates/scoring/src/scoring/scored_spectrum.rs:227` | _(undocumented — needs a description)_ |
| `ANDES_PEAK_WINDOW` | value | `crates/scoring/src/scoring/scored_spectrum.rs:224` | _(undocumented — needs a description)_ |
| `ANDES_PRECOFF_NOCLAMP` | presence | `crates/scoring/src/scoring/scored_spectrum.rs:269` | _(undocumented — needs a description)_ |
| `ANDES_RSS_PROBE` | presence | `crates/andes/src/bin/andes.rs:1127` | _(undocumented — needs a description)_ |
| `ANDES_SEED_GEOMETRY` | ? | `crates/andes/tests/train_from_msnet.rs:288` | _(undocumented — needs a description)_ |
| `ANDES_TEST_D` | presence | `crates/input/tests/timstof_d_loads.rs:18` | _(undocumented — needs a description)_ |
| `ANDES_TEST_PERCOLATOR_BIN` | presence | `crates/output/tests/percolator_integration.rs:18` | _(undocumented — needs a description)_ |
| `ANDES_TEST_RAW` | value | `crates/input/tests/thermo_raw.rs:16` | _(undocumented — needs a description)_ |
| `ANDES_TIGHT_HIGHRES` | value | `crates/scoring/src/scoring/scored_spectrum.rs:1538` | _(undocumented — needs a description)_ |
| `ANDES_TRAIN_BENCH` | value | `crates/model-train/tests/yield_nonregression.rs:107` | ── Skip guard ──────────────────────────────────────────────────────────── |
| `ANDES_V1_OUT` | value | `crates/model-train/tests/partition_parity.rs:208` | _(undocumented — needs a description)_ |
| `ANDES_V1_STORE` | value | `crates/model-train/tests/partition_parity.rs:179` | _(undocumented — needs a description)_ |
