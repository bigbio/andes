---
name: parameters-and-constants
description: Decide whether a number in andes is a user parameter (CLI flag), a physical constant, an algorithm constant or an engineering constant, and how each kind is defined, named, scoped and changed. Use when adding or changing a numeric literal, a const, a threshold, a tolerance, a default value, or a mass.
---

# Parameters and constants

Every number in non-test code is one of five kinds. Find its kind first; the kind decides
where it lives and what it takes to change it.

## The decision

Ask in order:

1. **Would a user with another instrument, sample or experiment need a different value?**
   It is a **parameter**: a typed CLI flag with a documented default (see
   `cli-flags-not-env-vars`). Tolerances, enzyme, missed cleavages, modifications, charge
   range, FDR thresholds, peptide length, top-N. Never a `const`.
2. **Is it a fact of nature or of a standard?** It is a **physical constant**: masses,
   isotope spacings, the proton. One definition, in `crates/model/src/mass.rs`.
3. **Is it a fixed fact of a file format or external spec?** It is a **format constant**:
   Unimod accessions, mzML CV terms, PIN column names. Defined once next to the reader or
   writer that needs it.
4. **Does it change results, but should be the same for every user?** It is an
   **algorithm constant**: a scoring window, a sigma multiple, a minimum PSM count for
   calibration, a GBDT depth. A named `const`.
5. **Does it only affect speed or memory?** It is an **engineering constant**: chunk sizes,
   buffer capacities, batch counts. A named `const` that must not change any result.

If you cannot tell 4 from 1, it is 1 until a benchmark shows one value is right for every
dataset (then the flag can be removed, per `cli-flags-not-env-vars`).

## Rules for every constant

- **No bare literals** in non-test code beyond `0`, `1`, `2`, `-1` and `0.5` used as
  arithmetic. Everything else gets a named `const` or comes from a parameter.
- **The unit is in the name** whenever there is one: `_PPM`, `_DA`, `_MZ`, `_SECS`, `_MIN`,
  `_BYTES`, `_SIGMAS`. `FEATURE_WINDOW_DA`, not `FEATURE_WINDOW`.
- **Narrowest scope**: private by default, `pub(crate)` if the crate shares it, `pub` only
  when another crate must agree on the same value (and then it lives in the lower crate).
- **Defined once.** Two crates that must agree import one definition; they never each
  type the value.
- **A `///` doc comment** states what it bounds and its unit. Why that value was chosen
  goes in the commit message (CLAUDE.md), with the measurement.
- **No `static mut`**, no mutable globals. A value that comes from the CLI and reaches hot
  code is installed once at startup into a `OnceLock` (`init_scoring_settings`,
  `init_y_max_charge`), then read-only.

## Kind-specific rules

**Physical constants**
- Only in `crates/model/src/mass.rs`; everyone else imports them
  (`model::mass::PROTON`, `model::mass::ISOTOPE`). Do not write `1.00727649`, `1.00335483`,
  `18.010565` or `17.026549` anywhere else, and do not redefine them as local `const`s.
- Each has its source (CODATA, IUPAC atomic masses, Unimod) in its doc comment and a check
  in `crates/model/tests/chemistry_constants.rs`.
- Derived values are computed from the base ones (`H2O = 2 * H + O`), not typed.
- Replacing a retyped literal with the shared constant must keep the PIN bit-identical;
  if the digits differ, that difference is a bug to report, not to paper over.

**Algorithm constants**
- Changing one is a result change: it needs the `fdr-validation` evidence (entrapment FDP,
  several datasets, seeds) in the commit message.
- If its best value differs by instrument or experiment class, move it into the model
  store or make it a parameter, rather than a `match` on instrument in code.
- Group the constants of one method together at the top of its module.

**Engineering constants**
- Must not change results: the PIN is bit-identical for any value. A test or the
  out-of-core goldens should cover two different values where practical.
- Size from the memory budget (`memlimit.rs`) or the thread count when the right value
  depends on the machine; a constant is only the fallback or the cap.

**Parameters**
- The default lives in one place (the clap `default_value` or the params struct's
  `Default`), never repeated as a literal elsewhere.
- Validate range and unit at parse time (a negative ppm or a zero charge is rejected with
  a message naming the flag).
- The resolved value is printed in the `run` banner (see `cli-reports`) so a log shows what
  was used.

## Review checklist

- [ ] No new bare numeric literals outside tests.
- [ ] No mass or isotope value typed outside `model::mass`.
- [ ] New `const` has a unit suffix, minimal visibility and a `///` doc.
- [ ] Result-changing constant change carries `fdr-validation` evidence.
- [ ] Engineering constant change keeps the PIN bit-identical.
- [ ] A value users need to vary is a flag, documented in DOCS.md and `--help`.
