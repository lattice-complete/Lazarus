//! Commitment key: lazily expanded uniform matrices over R_q, derived
//! from a seed (C: the global `comkey` expanded from a zero seed).
//!
//! Unlike the C implementation — which slices overlapping windows out of
//! one long uniform vector — every matrix here lives in its own domain
//! (A: inner commitments; B: outer commitment 1 keys for t-hat; C: outer
//! commitment 1 keys for g-hat; D: outer commitment 2 keys for h-hat;
//! W: greyhound's key for w-hat), removing key-reuse correlations while
//! staying deterministic in the seed.

use crate::Rq;
use algebra::ring::PolyVec;
use algebra::sampling::uniform::poly_uniform;
use std::collections::HashMap;
use std::sync::Mutex;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum KeyDomain {
    A = 0,
    B = 1,
    C = 2,
    D = 3,
    W = 4,
}

pub struct ComKey {
    pub seed: [u8; 32],
    cache: Mutex<HashMap<(u8, u32, u32), Rq>>,
}

impl ComKey {
    pub fn new(seed: [u8; 32]) -> Self {
        ComKey {
            seed,
            cache: Mutex::new(HashMap::new()),
        }
    }

    /// Uniform key polynomial at (domain, row, col).
    pub fn entry(&self, dom: KeyDomain, row: usize, col: usize) -> Rq {
        let key = (dom as u8, row as u32, col as u32);
        if let Some(p) = self.cache.lock().unwrap().get(&key) {
            return p.clone();
        }
        assert!(row < 1 << 20 && col < 1 << 20);
        let dom64 = ((dom as u64) << 40) | ((row as u64) << 20) | col as u64;
        let p: Rq = poly_uniform(&self.seed, dom64);
        self.cache.lock().unwrap().insert(key, p.clone());
        p
    }

    /// rows x len matrix-vector product: out[k] = sum_l K[k, col_off+l] s[l].
    pub fn mul(&self, dom: KeyDomain, rows: usize, col_off: usize, s: &[Rq]) -> PolyVec<64> {
        let mut out = PolyVec::zero(rows);
        for k in 0..rows {
            let mut acc = Rq::zero();
            for (l, sl) in s.iter().enumerate() {
                if sl.is_zero() {
                    continue;
                }
                acc = &acc + &(&self.entry(dom, k, col_off + l) * sl);
            }
            out.polys[k] = acc;
        }
        out
    }

    /// Challenge-collapsed key columns:
    /// out[l] = sum_k alpha[k] * K[k, col_off+l]  (C: polxvec_collaps_add_extension).
    pub fn collapse(&self, dom: KeyDomain, alphas: &[Rq], col_off: usize, len: usize) -> Vec<Rq> {
        (0..len)
            .map(|l| {
                let mut acc = Rq::zero();
                for (k, a) in alphas.iter().enumerate() {
                    acc = &acc + &(a * &self.entry(dom, k, col_off + l));
                }
                acc
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::zq::Zq;

    #[test]
    fn deterministic_and_domain_separated() {
        let k1 = ComKey::new([5u8; 32]);
        let k2 = ComKey::new([5u8; 32]);
        assert_eq!(k1.entry(KeyDomain::A, 1, 2), k2.entry(KeyDomain::A, 1, 2));
        assert_ne!(k1.entry(KeyDomain::A, 1, 2), k1.entry(KeyDomain::B, 1, 2));
        assert_ne!(k1.entry(KeyDomain::A, 1, 2), k1.entry(KeyDomain::A, 2, 1));
    }

    /// The identity every aggregation step relies on:
    /// <collapse(alpha, K), s> = <alpha, K s>.
    #[test]
    fn collapse_transpose_identity() {
        let key = ComKey::new([6u8; 32]);
        let mut c = [0i64; 64];
        c[0] = 3;
        c[5] = -2;
        let s: Vec<Rq> = (0..4)
            .map(|i| {
                let mut ci = c;
                ci[1] = i as i64;
                Rq::from_centered(&ci)
            })
            .collect();
        let alphas: Vec<Rq> = (0..3).map(|i| &Rq::monomial(i) * Zq::new(i + 1)).collect();

        let ks = key.mul(KeyDomain::B, 3, 7, &s);
        let lhs: Rq = alphas
            .iter()
            .zip(ks.polys.iter())
            .fold(Rq::zero(), |acc, (a, t)| &acc + &(a * t));

        let cols = key.collapse(KeyDomain::B, &alphas, 7, 4);
        let rhs: Rq = cols
            .iter()
            .zip(s.iter())
            .fold(Rq::zero(), |acc, (c, s)| &acc + &(c * s));

        assert_eq!(lhs, rhs);
    }
}
