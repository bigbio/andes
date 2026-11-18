---
name: parameters-and-constants
description: Infer values from the data instead of adding parameters, and decide whether a remaining number is a user parameter (CLI flag), a physical, format, algorithm or engineering constant, and how each is defined, named, scoped and changed. Use when proposing an algorithm, adding a flag, or adding or changing a numeric literal, a const, a threshold, a tolerance, a default value, or a mass.
---

# Parameters and constants

## First rule: infer from the data, avoid parameters

andes is data-driven. Every parameter is a choice the user can get wrong, a value
benchmarks must cover, and a line in `--help`. Keep the number of parameters as small as
possible.

When proposing or changing an algorithm, first ask **what the data can tell us**, and
estimate the value from the run itself before defining any parameter or constant. andes
already does this:

- fragment resolution read from the peaks (isotope-pair enrichment) when the file has no
  analyzer metadata;
- precursor mass shift and tolerance learned from confident PSMs (`--precursor-cal`);
- the scoring model selected from the detected instrument and experiment class;
- the memory budget taken from the cgroup limit or `MemAvailable`.

An inferred value must:

- come from a robust estimator (median, MAD, quantiles, mixture fit), not from one PSM;
- need a minimum amount of evidence, and fall back to a safe, documented default when the
  data is too thin (that minimum and the fallback are algorithm constants, below);
- be deterministic for the same input, whatever the thread count;
- be reported in the log: what was inferred, from how much data, and when the fallback was
  used (see `clean-logs` and the `run` banner in `cli-reports`).

A proposal for a new method says, for each number it needs: inferred from the data (how),
a constant (which kind, below), or a parameter (why it cannot be inferred). "Make it a
flag" is the last option, not the first.

When a flag is unavoidable, prefer `auto` as its default, with the inferred value, and let
the flag override it. A flag whose A/B is settled is removed (`cli-flags-not-env-vars`).

## Kinds of numbers

Every number in non-test code is one of five kinds. Find its kind first; the kind decides
where it lives and what it takes to change it.

## The decision

Ask in order:

1. **Would a user with another instrument, sample or experiment need a different value?**
   Then first try to **infer it from the data** (above). Only if it cannot be inferred
   (it encodes the user's intent, not a property of the data) is it a **parameter**: a typed CLI flag with a documented default (see
   `cli-flags-not-env-vars`). Tolerances, enzyme, missed cleavages, modifications, charge
   range, FDR thresholds. Never a `const`.
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

Enzyme, modifications and the FDR threshold encode the user's intent and are parameters.
Tolerances, mass shifts, resolution, charge states and score distributions are properties
of the data and should be inferred. If you cannot tell 4 from 1, try inference first, then
an algorithm constant that a benchmark shows is right for every dataset.

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

- [ ] Every value that is a property of the data is inferred from it, with a minimum
      evidence count, a fallback and a log line.
- [ ] No new flag that could have been inferred; a needed flag defaults to `auto`.
- [ ] No new bare numeric literals outside tests.
- [ ] No mass or isotope value typed outside `model::mass`.
- [ ] New `const` has a unit suffix, minimal visibility and a `///` doc.
- [ ] Result-changing constant change carries `fdr-validation` evidence.
- [ ] Engineering constant change keeps the PIN bit-identical.
- [ ] A value users need to vary is a flag, documented in DOCS.md and `--help`.
