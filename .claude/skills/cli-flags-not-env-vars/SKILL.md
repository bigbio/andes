---
name: cli-flags-not-env-vars
description: Policy that no environment variable changes an andes result - every setting is a typed CLI flag with a documented default. Use when adding a setting, switch, tuning knob or A/B toggle, when code reads std::env, or when crates/search/tests/env_registry.rs fails.
---

Read and follow docs/dev/cli-flags-not-env-vars.md before starting.
