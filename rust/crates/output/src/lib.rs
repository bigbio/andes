//! Placeholder for the `output` crate. Implementation lands in a later phase.

/// Returns the crate name. Used only to keep `cargo build` from emitting an
/// empty-crate warning.
pub fn placeholder() -> &'static str {
    "output"
}
