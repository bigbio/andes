//! Placeholder for the `core` crate. Implementation lands in a later phase.

/// Returns the crate name. Used only to keep `cargo build` from emitting an
/// empty-crate warning during M0.
pub fn placeholder() -> &'static str {
    "core"
}
