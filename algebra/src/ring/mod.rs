//! Ring layer: fixed-degree polynomial arithmetic over R_q = Z_q[X]/(X^D+1).
//!
//! Modern replacements for the legacy `PolynomialRing`/`RqMatrix` pair,
//! ported from lazer's `poly.c` / `polyvec.c` / `polymat.c`.

pub mod poly;
pub mod polymat;
pub mod polyvec;

pub use poly::{Poly, Rq64};
pub use polymat::PolyMat;
pub use polyvec::PolyVec;
