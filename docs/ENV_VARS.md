# `ANDES_*` environment variables

Enforced by `crates/search/tests/env_registry.rs`: a variable absent from this file
fails the build, and an entry the code no longer reads fails it too.

## Policy

**Algorithm behaviour is configured by CLI flags and `--config`, not by the environment.**
Environment variables are invisible in `--help`, untyped, unvalidated, and easy to tune
around and then forget — this codebase re-tuned selector weights around a bias whose
correction was sitting in the tree behind an unset variable.

So a switch that changes search results must be either:

* **deleted**, once its A/B is settled — the winning behaviour becomes the only one; or
* **a typed CLI flag**, discoverable in `--help` and settable from a config file.

What remains here is limited to test fixtures, model-training tools, and diagnostics.
`presence` gating (`var_os`) is a footgun — `VAR=0` *enables* — and new variables must
use `value` form.

| Variable | Kind | Form | Defined at |
| --- | --- | --- | --- |
| `ANDES_CHIMERIC_OVERLAP` | ALGORITHM | value | `crates/search/src/match_engine.rs:526` |
| `ANDES_DENSITY_RAW` | ALGORITHM | presence | `crates/scoring/src/scoring/scored_spectrum.rs:780` |
| `ANDES_ETHCD_AS_ETD` | ALGORITHM | presence | `crates/input/src/mzml.rs:172` |
| `ANDES_PEAK_PER_WINDOW` | ALGORITHM | value | `crates/scoring/src/scoring/scored_spectrum.rs:227` |
| `ANDES_PEAK_WINDOW` | ALGORITHM | value | `crates/scoring/src/scoring/scored_spectrum.rs:224` |
| `ANDES_PRECOFF_NOCLAMP` | ALGORITHM | presence | `crates/scoring/src/scoring/scored_spectrum.rs:269` |
| `ANDES_TIGHT_HIGHRES` | ALGORITHM | value | `crates/scoring/src/scoring/scored_spectrum.rs:1549` |
| `ANDES_GLYCO_SCANS` | diagnostic | presence | `crates/search/src/glyco_search.rs:592` |
| `ANDES_RSS_PROBE` | diagnostic | presence | `crates/andes/src/bin/andes.rs:1167` |
| `ANDES_DENSE_NOISE` | model training | value | `crates/model-train/src/accumulate.rs:110` |
| `ANDES_GEO_MAX_FRAG_CHARGE` | model training | ? | `crates/andes/src/bin/andes.rs:1261` |
| `ANDES_GEO_MAX_RANK` | model training | ? | `crates/andes/src/bin/andes.rs:1258` |
| `ANDES_GEO_MAX_TIERS` | model training | ? | `crates/andes/src/bin/andes.rs:1260` |
| `ANDES_GEO_OCCUPANCY` | model training | ? | `crates/andes/src/bin/andes.rs:1259` |
| `ANDES_GEO_SEGMENTS` | model training | ? | `crates/andes/src/bin/andes.rs:1257` |
| `ANDES_SEED_GEOMETRY` | model training | ? | `crates/andes/tests/train_from_msnet.rs:288` |
| `ANDES_V1_OUT` | model training | value | `crates/model-train/tests/partition_parity.rs:208` |
| `ANDES_V1_STORE` | model training | value | `crates/model-train/tests/partition_parity.rs:179` |
| `ANDES_TEST_D` | test-harness | presence | `crates/input/tests/timstof_d_loads.rs:18` |
| `ANDES_TEST_PERCOLATOR_BIN` | test-harness | presence | `crates/output/tests/percolator_integration.rs:18` |
| `ANDES_TEST_RAW` | test-harness | value | `crates/input/tests/thermo_raw.rs:16` |
| `ANDES_TRAIN_BENCH` | test-harness | value | `crates/model-train/tests/yield_nonregression.rs:107` |
