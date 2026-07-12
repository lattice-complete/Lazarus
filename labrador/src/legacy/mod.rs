//! The original one-round LaBRADOR prototype (pre-rewrite), kept for reference.
//! Built on the legacy PolynomialRing types; superseded by the modules in the
//! crate root, which implement the full recursive protocol on algebra::ring.

pub mod core;
pub mod prover;
pub mod setup;
pub mod statement;
pub mod transcript;
pub mod verifier;
