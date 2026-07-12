//! Number-theoretic transform engine.
//!
//! Design ported from the lazer C library (`src/moduli.h`, `src/ntt.c`):
//! since the proof-system modulus q (here q = 2^32 - 99, which has 2-adicity
//! of q-1 equal to 2) is *not* NTT-friendly, multiplication in
//! R = Z_q[X]/(X^d + 1) is performed by
//!
//!   1. lifting both operands to centered integer polynomials,
//!   2. computing the exact negacyclic product over Z via NTTs modulo a few
//!      NTT-friendly primes (residue number system),
//!   3. recombining with the explicit CRT and reducing mod q.
//!
//! lazer uses a battery of ~50-bit primes tuned for AVX-512; this portable
//! scalar port uses two ~64/62-bit primes instead, which is the natural
//! choice for u64/u128 arithmetic. The CRT modulus P = p0 * p1 > 2^113
//! comfortably bounds the largest possible product coefficient
//! d * (q/2)^2 < 2^74 for every ring degree d <= 4096 used by lazer.
//!
//! The transforms are the standard in-place negacyclic NTT with
//! bit-reversed twiddle tables (Longa–Naehrig style): Cooley–Tukey
//! butterflies forward, Gentleman–Sande backward, so no explicit
//! bit-reversal pass is needed.

mod prime;

pub use prime::NttPrime;

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

/// The two NTT primes backing the CRT engine.
///
/// * `p0` = 2^64 - 2^32 + 1 (the "Goldilocks" prime), 2-adicity 32.
/// * `p1` is found by a deterministic downward search from 2^62 over
///   candidates of the form m * 2^20 + 1, giving 2-adicity >= 20.
///
/// Both are (re-)verified with a deterministic Miller–Rabin test at first
/// use, and their roots of unity are derived and order-checked at runtime,
/// so no unverifiable magic constants enter the arithmetic.
pub fn primes() -> &'static [NttPrime; 2] {
    static PRIMES: OnceLock<[NttPrime; 2]> = OnceLock::new();
    PRIMES.get_or_init(|| {
        let p0 = NttPrime::from_verified(0xFFFF_FFFF_0000_0001, 32);
        let p1 = NttPrime::search_below_2_62();
        assert_ne!(p0.p, p1.p);
        [p0, p1]
    })
}

/// Precomputed twiddle tables for one (degree, prime) pair.
pub struct NttTable {
    pub d: usize,
    pub p: u64,
    /// psi^brv(i) for i in 0..d, psi a primitive 2d-th root of unity mod p.
    psi_rev: Vec<u64>,
    /// psi^{-brv(i)} for i in 0..d.
    ipsi_rev: Vec<u64>,
    /// d^{-1} mod p, folded into the last inverse layer.
    d_inv: u64,
}

fn bitrev(x: usize, bits: u32) -> usize {
    x.reverse_bits() >> (usize::BITS - bits)
}

impl NttTable {
    fn new(d: usize, prime: &NttPrime) -> Self {
        assert!(d.is_power_of_two() && d >= 2, "NTT degree must be 2^k >= 2");
        let log2d = d.trailing_zeros();
        let psi = prime.root_of_order(2 * d as u64);
        let ipsi = prime.pow(psi, prime.p - 2);
        let mut pw: u64 = 1;
        let mut ipw: u64 = 1;
        let mut pows = Vec::with_capacity(d);
        let mut ipows = Vec::with_capacity(d);
        for _ in 0..d {
            pows.push(pw);
            ipows.push(ipw);
            pw = prime.mul(pw, psi);
            ipw = prime.mul(ipw, ipsi);
        }
        let psi_rev: Vec<u64> = (0..d).map(|i| pows[bitrev(i, log2d)]).collect();
        let ipsi_rev: Vec<u64> = (0..d).map(|i| ipows[bitrev(i, log2d)]).collect();
        let d_inv = prime.pow(d as u64, prime.p - 2);
        NttTable {
            d,
            p: prime.p,
            psi_rev,
            ipsi_rev,
            d_inv,
        }
    }

    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 {
        let (s, of) = a.overflowing_add(b);
        // p < 2^64, inputs < p, so a+b < 2p; subtract p on overflow or excess.
        if of || s >= self.p {
            s.wrapping_sub(self.p)
        } else {
            s
        }
    }

    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a.wrapping_sub(b).wrapping_add(self.p)
        }
    }

    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }

    /// In-place forward negacyclic NTT (Cooley–Tukey, natural -> bit-reversed).
    pub fn forward(&self, a: &mut [u64]) {
        debug_assert_eq!(a.len(), self.d);
        let mut t = self.d;
        let mut m = 1;
        while m < self.d {
            t >>= 1;
            for i in 0..m {
                let j1 = 2 * i * t;
                let s = self.psi_rev[m + i];
                for j in j1..j1 + t {
                    let u = a[j];
                    let v = self.mul(a[j + t], s);
                    a[j] = self.add(u, v);
                    a[j + t] = self.sub(u, v);
                }
            }
            m <<= 1;
        }
    }

    /// In-place inverse negacyclic NTT (Gentleman–Sande, bit-reversed -> natural).
    pub fn inverse(&self, a: &mut [u64]) {
        debug_assert_eq!(a.len(), self.d);
        let mut t = 1;
        let mut m = self.d;
        while m > 1 {
            let h = m >> 1;
            let mut j1 = 0;
            for i in 0..h {
                let s = self.ipsi_rev[h + i];
                for j in j1..j1 + t {
                    let u = a[j];
                    let v = a[j + t];
                    a[j] = self.add(u, v);
                    a[j + t] = self.mul(self.sub(u, v), s);
                }
                j1 += 2 * t;
            }
            t <<= 1;
            m = h;
        }
        for x in a.iter_mut() {
            *x = self.mul(*x, self.d_inv);
        }
    }
}

type NttCache = Mutex<HashMap<(usize, u64), Arc<NttTable>>>;

/// Global table cache keyed by (degree, prime).
fn table(d: usize, prime: &'static NttPrime) -> Arc<NttTable> {
    static CACHE: OnceLock<NttCache> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    let mut guard = cache.lock().expect("ntt table cache poisoned");
    guard
        .entry((d, prime.p))
        .or_insert_with(|| Arc::new(NttTable::new(d, prime)))
        .clone()
}

#[inline]
fn to_residues(a: &[i64], p: u64) -> Vec<u64> {
    a.iter()
        .map(|&x| {
            let r = x % p as i64;
            if r < 0 {
                (r + p as i64) as u64
            } else {
                r as u64
            }
        })
        .collect()
}

/// Exact negacyclic product over the integers.
///
/// Given centered coefficient vectors `a`, `b` of length d (a power of two),
/// returns the coefficients of a*b mod X^d + 1 as exact centered integers
/// (no reduction mod q — the caller chooses the ring).
///
/// Preconditions (checked): |a_i|, |b_i| <= 2^31 and d <= 2^12, which keeps
/// every product coefficient below d * 2^62 <= 2^74 << P/2.
pub fn negacyclic_mul_centered(a: &[i64], b: &[i64]) -> Vec<i128> {
    let d = a.len();
    assert_eq!(d, b.len());
    assert!(d.is_power_of_two() && d >= 2);
    assert!(d <= 1 << 12, "degree too large for the 2-prime CRT bound");
    assert!(
        a.iter()
            .chain(b.iter())
            .all(|&x| x.unsigned_abs() <= 1 << 31),
        "inputs must be centered representatives (|x| <= 2^31)"
    );

    let [pr0, pr1] = primes();
    let mut residues: [Vec<u64>; 2] = [Vec::new(), Vec::new()];
    for (k, pr) in [pr0, pr1].into_iter().enumerate() {
        let tb = table(d, pr);
        let mut fa = to_residues(a, pr.p);
        let mut fb = to_residues(b, pr.p);
        tb.forward(&mut fa);
        tb.forward(&mut fb);
        for (x, y) in fa.iter_mut().zip(fb.iter()) {
            *x = tb.mul(*x, *y);
        }
        tb.inverse(&mut fa);
        residues[k] = fa;
    }

    // Explicit CRT (Garner): x = r0 + p0 * ((r1 - r0) * p0^{-1} mod p1),
    // then center mod P = p0*p1.
    let p0 = pr0.p as u128;
    let p1 = pr1.p as u128;
    let p0_inv_p1 = {
        // p0^{-1} mod p1 via Fermat (p1 prime).
        let base = pr0.p % pr1.p;
        pr1.pow(base, pr1.p - 2) as u128
    };
    let big_p = p0 * p1;
    let half_p = big_p >> 1;

    residues[0]
        .iter()
        .zip(residues[1].iter())
        .map(|(&r0, &r1)| {
            let r0u = r0 as u128;
            let diff = (r1 as u128 + p1 - (r0 as u128 % p1)) % p1;
            let t = diff * p0_inv_p1 % p1;
            let x = r0u + p0 * t; // in [0, P)
            if x > half_p {
                x as i128 - big_p as i128
            } else {
                x as i128
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference negacyclic schoolbook product over i128.
    fn schoolbook(a: &[i64], b: &[i64]) -> Vec<i128> {
        let d = a.len();
        let mut out = vec![0i128; d];
        for i in 0..d {
            for j in 0..d {
                let prod = a[i] as i128 * b[j] as i128;
                if i + j < d {
                    out[i + j] += prod;
                } else {
                    out[i + j - d] -= prod;
                }
            }
        }
        out
    }

    #[test]
    fn primes_are_valid() {
        let [p0, p1] = primes();
        assert_eq!(p0.p, 0xFFFF_FFFF_0000_0001);
        assert!(p1.p > 1 << 61);
        // 2-adicity must support degree 4096 (needs 2^13th roots).
        assert!(p0.two_adicity >= 13 && p1.two_adicity >= 13);
    }

    #[test]
    fn ntt_roundtrip_identity() {
        let [pr0, _] = primes();
        let tb = table(64, pr0);
        let orig: Vec<u64> = (0..64u64).map(|i| i * i + 7).collect();
        let mut a = orig.clone();
        tb.forward(&mut a);
        assert_ne!(a, orig);
        tb.inverse(&mut a);
        assert_eq!(a, orig);
    }

    #[test]
    fn ntt_matches_schoolbook_small() {
        // (1 + 2X) * (3 + 4X) mod X^2 + 1 = 3 + 10X + 8X^2 = -5 + 10X
        let got = negacyclic_mul_centered(&[1, 2], &[3, 4]);
        assert_eq!(got, vec![-5, 10]);
    }

    #[test]
    fn ntt_matches_schoolbook_random_d64() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(0xA1B);
        for _ in 0..20 {
            let a: Vec<i64> = (0..64)
                .map(|_| rng.gen_range(-(1 << 31)..=1 << 31))
                .collect();
            let b: Vec<i64> = (0..64)
                .map(|_| rng.gen_range(-(1 << 31)..=1 << 31))
                .collect();
            assert_eq!(negacyclic_mul_centered(&a, &b), schoolbook(&a, &b));
        }
    }

    #[test]
    fn ntt_matches_schoolbook_other_degrees() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        for d in [2usize, 128, 512] {
            let a: Vec<i64> = (0..d)
                .map(|_| rng.gen_range(-(1 << 31)..=1 << 31))
                .collect();
            let b: Vec<i64> = (0..d)
                .map(|_| rng.gen_range(-(1 << 31)..=1 << 31))
                .collect();
            assert_eq!(negacyclic_mul_centered(&a, &b), schoolbook(&a, &b), "d={d}");
        }
    }
}
