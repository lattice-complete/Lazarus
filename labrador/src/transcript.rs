//! Fiat-Shamir transcript: a 32-byte rolling hash state with typed,
//! domain-tagged challenge derivations.
//!
//! Mirrors the reference's hash chaining (h <- SHAKE128(h || message))
//! structurally, but is NOT byte-compatible with the C implementation:
//! the state is 32 bytes (C: 16), messages are absorbed as canonical
//! little-endian Zq words, and every derivation carries a domain tag.

use crate::Rq;
use algebra::ring::PolyVec;
use algebra::rng::Shake128Rng;
use algebra::sampling::challenge::{challenge_polys, quarternary_polys};
use algebra::sampling::uniform::zq_uniform;
use algebra::zq::Zq;
use sha3::digest::{ExtendableOutput, Update, XofReader};

pub type Hash = [u8; 32];

/// h <- SHAKE128(tag || h || parts), 32 bytes out.
pub fn absorb(h: &mut Hash, tag: u8, parts: &[&[u8]]) {
    let mut hasher = sha3::Shake128::default();
    hasher.update(&[tag]);
    hasher.update(h.as_slice());
    for p in parts {
        hasher.update(p);
    }
    let mut reader = hasher.finalize_xof();
    reader.read(h);
}

pub fn absorb_polyvec(h: &mut Hash, tag: u8, v: &PolyVec<64>) {
    absorb(h, tag, &[&v.to_bytes()]);
}

pub fn absorb_poly(h: &mut Hash, tag: u8, p: &Rq) {
    absorb(h, tag, &[&p.to_bytes()]);
}

/// Derivation XOF seeded by the current state.
pub fn xof(h: &Hash, dom: u64) -> Shake128Rng {
    Shake128Rng::new(h, dom)
}

pub fn zq_scalars(h: &Hash, dom: u64, count: usize) -> Vec<Zq> {
    let mut x = xof(h, dom);
    (0..count).map(|_| zq_uniform(&mut x)).collect()
}

pub fn quarternary(h: &Hash, dom: u64, count: usize) -> Vec<Rq> {
    let mut x = xof(h, dom);
    quarternary_polys(&mut x, count)
}

pub fn challenges(h: &Hash, dom: u64, count: usize) -> Vec<Rq> {
    let mut x = xof(h, dom);
    challenge_polys(&mut x, count)
}

/// Domain tags (absorb) and XOF domains, kept distinct per use.
pub mod tags {
    pub const U1: u8 = 1;
    pub const JLPROJ: u8 = 2;
    pub const JLCOLLAPS: u8 = 3;
    pub const SPARSE: u8 = 4;
    pub const LIFT: u8 = 5;
    pub const AGG: u8 = 6;
    pub const U2: u8 = 7;
    pub const CHAL: u8 = 8;
    pub const GREYHOUND_XY: u8 = 9;
    pub const TAIL_H: u8 = 10;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chaining_is_deterministic_and_tag_separated() {
        let mut h1 = [0u8; 32];
        let mut h2 = [0u8; 32];
        absorb(&mut h1, 1, &[b"data"]);
        absorb(&mut h2, 1, &[b"data"]);
        assert_eq!(h1, h2);
        let mut h3 = [0u8; 32];
        absorb(&mut h3, 2, &[b"data"]);
        assert_ne!(h1, h3);
        // rolling: absorbing again changes the state
        let before = h1;
        absorb(&mut h1, 1, &[b"more"]);
        assert_ne!(before, h1);
    }

    #[test]
    fn derivations_reproducible() {
        let h = [7u8; 32];
        assert_eq!(zq_scalars(&h, 3, 5), zq_scalars(&h, 3, 5));
        assert_ne!(zq_scalars(&h, 3, 5), zq_scalars(&h, 4, 5));
        assert_eq!(challenges(&h, 0, 2), challenges(&h, 0, 2));
        assert_eq!(quarternary(&h, 1, 3), quarternary(&h, 1, 3));
    }
}
