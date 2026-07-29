# `ANDES_*` environment variables

Enforced by `crates/search/tests/env_registry.rs`: a variable absent from this file fails
the build, and an entry the code no longer reads fails it too.

## Policy

**The engine reads no environment variables.** Everything that affects a search — or a
model-training run — is a CLI flag with a documented default, so a result can be
reproduced from the command line that produced it.

Environment variables were removed because they are invisible in `--help`, untyped,
unvalidated, and easy to tune around and then forget. This codebase re-tuned selector
weights around a Y-ladder size bias whose correction was sitting in the tree behind an
unset variable, and shipped a c/z charge guard that was never wired to a call site.

A switch that changes results is therefore either **deleted** once its A/B is settled —
the winning behaviour becomes the only behaviour — or **a typed CLI flag**. Values that
reach hot inner functions are installed once at startup from validated CLI input
(`ScoringSettings`, `CzSettings`, `init_y_max_charge`, `init_ethcd_as_etd`,
`init_dense_noise`) rather than read from the environment on each call.

What remains below is limited to **test-harness variables**, which select optional
fixtures at `cargo test` time and are never read by the shipped binary.

| Variable | Kind | Defined at |
| --- | --- | --- |
| `ANDES_SEED_GEOMETRY` | test-harness | `crates/andes/tests/train_from_msnet.rs:288` |
| `ANDES_TEST_D` | test-harness | `crates/input/tests/timstof_d_loads.rs:18` |
| `ANDES_TEST_PERCOLATOR_BIN` | test-harness | `crates/output/tests/percolator_integration.rs:18` |
| `ANDES_TEST_RAW` | test-harness | `crates/input/tests/thermo_raw.rs:16` |
| `ANDES_TRAIN_BENCH` | test-harness | `crates/model-train/tests/yield_nonregression.rs:107` |
| `ANDES_V1_OUT` | test-harness | `crates/model-train/tests/partition_parity.rs:208` |
| `ANDES_V1_STORE` | test-harness | `crates/model-train/tests/partition_parity.rs:179` |
