//! Samplers over Z and R_q (ports of lazer's `urandom.c`, `brandom.c`,
//! `grandom.c`, `rejection.c`).
//!
//! Every sampler is deterministic given a 32-byte seed and a 64-bit domain
//! separator, mirroring lazer's `(seed, dom)` convention so that provers
//! and verifiers can re-derive identical randomness.

pub mod binomial;
pub mod challenge;
pub mod gaussian;
pub mod rejection;
pub mod uniform;

pub use binomial::{cbd_vec, poly_cbd};
pub use challenge::{challenge_poly, challenge_polys, opnorm, quarternary_poly, quarternary_polys};
pub use gaussian::{gaussian_vec, poly_gaussian, GaussianSampler};
pub use rejection::{rej_bimodal, rej_standard};
pub use uniform::{poly_uniform, polyvec_uniform, u64_uniform_mod, zq_uniform};
