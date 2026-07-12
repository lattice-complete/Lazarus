//! # algebra — the Lazarus foundation layer
//!
//! Portable pure-Rust port of the arithmetic core of the lazer / labrador
//! C libraries (lattice-based zero-knowledge proofs):
//!
//! * [`zq`] — the base field Z_q, q = 2^32 - 99 (LaBRADOR modulus)
//! * [`ring`] — fixed-degree polynomial rings R_q = Z_q[X]/(X^D + 1),
//!   vectors and matrices over them ([`ring::Rq64`] is the LaBRADOR ring)
//! * [`ntt`] — CRT/NTT multiplication engine for non-NTT-friendly q
//! * [`rng`] — deterministic, domain-separated XOFs (SHAKE128 / AES256-CTR)
//! * [`sampling`] — uniform, centered-binomial, discrete-Gaussian and
//!   rejection samplers
//!
//! The legacy types (`PolynomialRing`, `RqMatrix`, loose helper functions)
//! are kept intact for the existing labrador crate; new code should prefer
//! the `ring` module. `PolynomialRing` multiplication is transparently
//! accelerated through the NTT engine when both operands have full degree.

// ---- legacy modules (kept API-stable) ----
pub mod polynomial_ring;
pub mod rq_matrix;
pub mod utils;
pub mod zq;

// ---- modern foundation layer ----
pub mod ntt;
pub mod ring;
pub mod rng;
pub mod sampling;

pub use polynomial_ring::PolynomialRing;
pub use rq_matrix::RqMatrix;
pub use utils::*;
pub use zq::Zq;

pub use ring::{Poly, PolyMat, PolyVec, Rq64};
